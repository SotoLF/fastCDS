"""Error-handling tests — verify the documented `status` values and that the
unmapped_domains.tsv carries the right `reason`."""

from __future__ import annotations

import csv

from conftest import run_mapping


def _read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def test_protein_not_in_index(binary, with_tags_index, tmp_path):
    run_mapping(binary, with_tags_index,
                "ENSP99\t1\t10\tNOT_FOUND\n",
                tmp_path)
    summary = {r["input_id"]: r for r in _read_tsv(tmp_path / "domain_mapping_summary.tsv")}
    unmapped = {r["input_id"]: r for r in _read_tsv(tmp_path / "unmapped_domains.tsv")}

    assert summary["NOT_FOUND"]["status"] == "protein_not_in_index"
    assert unmapped["NOT_FOUND"]["reason"] == "protein_not_in_index"
    # No mapping coordinates for an unfound protein.
    assert summary["NOT_FOUND"]["domain_genomic_start"] == "NA"


def test_domain_beyond_protein_length(binary, with_tags_index, tmp_path):
    # ENSP1 is 44 aa long; aa 100..110 is past the end.
    run_mapping(binary, with_tags_index,
                "ENSP1\t100\t110\tBEYOND\n",
                tmp_path)
    summary = {r["input_id"]: r for r in _read_tsv(tmp_path / "domain_mapping_summary.tsv")}
    unmapped = {r["input_id"]: r for r in _read_tsv(tmp_path / "unmapped_domains.tsv")}

    assert summary["BEYOND"]["status"] == "domain_beyond_protein_length"
    assert unmapped["BEYOND"]["reason"] == "domain_beyond_protein_length"
    # Coordinate fields on the summary row are NA for unmapped queries.
    assert summary["BEYOND"]["domain_genomic_start"] == "NA"
    assert summary["BEYOND"]["n_coding_segments"] == "0"


def test_no_CDS_for_protein(binary, with_tags_index, tmp_path):
    # ENST6 is the non-coding transcript in the synthetic GTF (no CDS records).
    run_mapping(binary, with_tags_index,
                "ENST6\t1\t10\tNON_CODING\n",
                tmp_path)
    summary = {r["input_id"]: r for r in _read_tsv(tmp_path / "domain_mapping_summary.tsv")}
    unmapped = {r["input_id"]: r for r in _read_tsv(tmp_path / "unmapped_domains.tsv")}

    assert summary["NON_CODING"]["status"] == "no_CDS_for_protein"
    assert unmapped["NON_CODING"]["reason"] == "no_CDS_for_protein"
    # Coordinate / type fields are NA when the lookup fails before strand resolution.
    assert summary["NON_CODING"]["domain_genomic_start"] == "NA"
