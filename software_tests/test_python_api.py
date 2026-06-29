"""Python wrapper API tests.

The C++ binary tests (test_correctness, test_errors, test_schema) cover the
mapper. These tests cover the *wrapper*: that the public API around the binary
parses outputs into DataFrames, supports both single and batch calls, persists
outputs when asked, and round-trips through read_results_dir.

The plotter is only smoke-tested (file exists + non-empty) because matplotlib
output is not the wrapper's contract.

Requires the wrapper to be installed:
    pip install -e python/
"""

from __future__ import annotations

from pathlib import Path

import pytest

# Skip the whole file cleanly if the wrapper hasn't been installed —
# the C++-only Phase 1 tests can still run without it.
fastCDS = pytest.importorskip("fastCDS")


def test_build_index_python_mirror(binary, synthetic_gtfs, tmp_path):
    """`fc.build_index(...)` is the Python mirror of `fastCDS index`: it
    returns the .idx path and the result is a usable index."""
    py_idx = tmp_path / "py_built.idx"
    built = fastCDS.build_index(synthetic_gtfs["with_tags"], out=py_idx,
                                  binary=str(binary))
    assert Path(built) == py_idx and py_idx.exists()
    mapper = fastCDS.Mapper(index=str(py_idx), binary=str(binary))
    assert mapper.map("ENSP1", aa_start=1, aa_end=10, domain_id="BI").n_mapped == 1


def test_mapper_map_single(binary, with_tags_index, tmp_path):
    """`Mapper.map(...)` runs one query and returns a MappingResult with parsed
    summary + isoform DataFrames carrying the expected columns."""
    mapper = fastCDS.Mapper(index=str(with_tags_index), binary=str(binary))
    result = mapper.map("ENSP1", aa_start=1, aa_end=10, domain_id="API_ONE")

    assert result.n_mapped == 1
    assert result.n_unmapped == 0
    assert result.summary.iloc[0]["gene_name"] == "TEST1"
    assert result.summary.iloc[0]["status"] == "ok"
    # The isoform table should contain rows for the whole transcript.
    assert len(result.isoform) > 0
    assert "plot_group" in result.isoform.columns
    # BED12 has exactly one row for a single domain query.
    assert len(result.bed12) == 1


def test_mapper_map_batch_and_slice(binary, with_tags_index, tmp_path):
    """`Mapper.map_batch(...)` does a single binary invocation for many queries
    and supports per-query slicing via by_input_id."""
    mapper = fastCDS.Mapper(index=str(with_tags_index), binary=str(binary))
    queries = [
        {"protein_id": "ENSP1", "aa_start": 1, "aa_end": 10, "domain_id": "B1"},
        {"protein_id": "ENSP2", "aa_start": 1, "aa_end": 17, "domain_id": "B2"},
        {"protein_id": "ENSP1"},  # no aa range → structure_only
    ]
    batch = mapper.map_batch(queries)

    assert batch.n_total == 3
    assert batch.n_mapped == 3
    assert batch.n_unmapped == 0

    # Slice to a single query.
    just_b2 = batch.by_input_id("B2")
    assert len(just_b2.summary) == 1
    assert just_b2.summary.iloc[0]["protein_id"] == "ENSP2"
    # The structure_only row should have NA aa_start / aa_end.
    struct_row = batch.summary[batch.summary["status"] == "structure_only"]
    assert len(struct_row) == 1


def test_mapper_map_unmapped_query(binary, with_tags_index, tmp_path):
    """An ENSP that doesn't exist is reported in `result.unmapped`, and
    `n_unmapped` reflects it."""
    mapper = fastCDS.Mapper(index=str(with_tags_index), binary=str(binary))
    result = mapper.map("ENSP_NOTFOUND", aa_start=1, aa_end=10, domain_id="MISS")
    assert result.n_mapped == 0
    assert result.n_unmapped == 1
    assert len(result.unmapped) == 1
    assert result.unmapped.iloc[0]["reason"] == "protein_not_in_index"


def test_keep_outputs_then_read_results_dir(binary, with_tags_index, tmp_path):
    """`keep_outputs=...` persists the binary's TSVs/BEDs on disk; `read_results_dir`
    reads them back into a MappingResult identical to the in-memory one."""
    mapper = fastCDS.Mapper(index=str(with_tags_index), binary=str(binary))
    persist = tmp_path / "persist"
    fresh = mapper.map("ENSP1", aa_start=1, aa_end=10, domain_id="P",
                       keep_outputs=str(persist))
    # All standard files on disk.
    for name in ("domain_mapping_summary.tsv", "isoform_structure.tsv",
                 "domain_cds_segments.tsv", "run_metadata.json"):
        assert (persist / name).exists(), f"keep_outputs missed {name}"

    re_read = fastCDS.read_results_dir(str(persist))
    assert re_read.n_mapped == fresh.n_mapped
    # Same number of CDS-segment rows in both.
    assert len(re_read.cds_segments) == len(fresh.cds_segments)
    # Metadata round-trips.
    assert re_read.metadata.get("tool") == "fastCDS"


def test_plot_smoke(binary, with_tags_index, tmp_path):
    """`fc.plot(MappingResult, out=...)` and `fc.plot(DataFrame, ...)` both
    produce a non-empty PDF. Visual content is the CLI plotter's responsibility,
    not the wrapper's — we only check the file is written."""
    mapper = fastCDS.Mapper(index=str(with_tags_index), binary=str(binary))
    result = mapper.map("ENSP1", aa_start=1, aa_end=10, domain_id="PLOT_ME")
    pdf_a = tmp_path / "from_result.pdf"
    fastCDS.plot(result, out=str(pdf_a))
    assert pdf_a.exists() and pdf_a.stat().st_size > 0

    pdf_b = tmp_path / "from_df.pdf"
    fastCDS.plot(result.isoform, input_id="PLOT_ME", out=str(pdf_b))
    assert pdf_b.exists() and pdf_b.stat().st_size > 0
