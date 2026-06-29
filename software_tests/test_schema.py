"""Schema / coordinate-convention sanity checks. No golden files — these
verify documented invariants directly."""

from __future__ import annotations

import csv
import json

from conftest import run_mapping

EXPECTED_FILES_OUTPUT_ALL = [
    # TSVs (4)
    "domain_mapping_summary.tsv",
    "isoform_structure.tsv",
    "domain_cds_segments.tsv",
    "domain_introns.tsv",
    # BEDs (3 companion + 1 BED12)
    "domain_cds_segments.bed",
    "domain_introns.bed",
    "domain_span_with_introns.bed",
    "domain_blocks.bed12",
    # Metadata
    "run_metadata.json",
]


def test_output_files_exist(binary, with_tags_index, tmp_path):
    """`--output all` writes every documented file (4 TSV + 3 BED + 1 BED12 + 1 JSON)."""
    run_mapping(binary, with_tags_index,
                "ENSP1\t1\t10\tSCHEMA\n",
                tmp_path)
    missing = [name for name in EXPECTED_FILES_OUTPUT_ALL
               if not (tmp_path / name).exists()]
    assert not missing, f"missing files under --output all: {missing}"

    meta = json.loads((tmp_path / "run_metadata.json").read_text())
    assert meta["tool"] == "fastCDS"
    assert meta["output_kind"] == "all"
    assert meta["query_counts"]["total"] == 1
    assert meta["query_counts"]["mapped"] == 1


def test_coordinate_conventions(binary, with_tags_index, tmp_path):
    """BED is 0-based half-open; TSV is 1-based inclusive. For the same CDS slice
    both representations must agree: bed.start == tsv.start - 1, bed.end == tsv.end,
    and length is the same in both."""
    # ENSP1 aa 1..10 is exactly CDS_1 on chrA: genomic 120..149 (TSV), length 30.
    run_mapping(binary, with_tags_index,
                "ENSP1\t1\t10\tCOORDS\n",
                tmp_path)

    # TSV side — read the single coding_overlap row.
    with open(tmp_path / "domain_cds_segments.tsv") as f:
        cds_rows = [r for r in csv.DictReader(f, delimiter="\t")
                    if r["overlaps_domain"] == "coding_overlap"]
    assert len(cds_rows) == 1
    tsv_start = int(cds_rows[0]["feature_genomic_start"])
    tsv_end = int(cds_rows[0]["feature_genomic_end"])
    tsv_length = tsv_end - tsv_start + 1

    # BED side — the companion bed should have exactly one row for this query.
    bed_rows = [line.rstrip("\n").split("\t")
                for line in (tmp_path / "domain_cds_segments.bed").read_text().splitlines()
                if line and not line.startswith("#")]
    assert len(bed_rows) == 1
    bed_start = int(bed_rows[0][1])  # 0-based inclusive
    bed_end = int(bed_rows[0][2])    # 0-based exclusive
    bed_length = bed_end - bed_start

    assert bed_start == tsv_start - 1, "BED start should be TSV start − 1"
    assert bed_end == tsv_end, "BED end (exclusive) should equal TSV end (inclusive)"
    assert bed_length == tsv_length, "length must be identical in both conventions"
