"""Correctness tests — golden-file diffs of the mapper's TSV/BED outputs, plus
a couple of value-level checks (selenoprotein CDS mismatch, ENST≡ENSP) that read
from the shared `out_all` run."""

from __future__ import annotations

from conftest import assert_matches_golden, run_mapping, summary_by_id

GOLDEN_FILES = [
    "domain_mapping_summary.tsv",
    "isoform_structure.tsv",
    "domain_cds_segments.tsv",
    "domain_cds_segments.bed",
    "domain_introns.tsv",
    "domain_introns.bed",
    "domain_span_with_introns.bed",
    "domain_blocks.bed",
]


def test_plus_strand_multi_exon_domain(
    binary, with_tags_index, tmp_path, update_goldens,
):
    """ENSP1 (plus strand) aa 1..27 covers all of CDS_1 (aa 1..10) and CDS_2 (aa 11..27).
    Exercises the multi-exon-coding-segments path on a + strand transcript."""
    run_mapping(binary, with_tags_index,
                "ENSP1\t1\t27\tMULTI_PLUS\n",
                tmp_path)
    assert_matches_golden(tmp_path, "test_plus_strand_multi_exon_domain",
                          GOLDEN_FILES, update=update_goldens)


def test_minus_strand_multi_exon_domain(
    binary, with_tags_index, tmp_path, update_goldens,
):
    """ENSP2 (minus strand) aa 1..34 covers CDS_1 (aa 1..17, genomic 800..850) and
    CDS_2 (aa 18..34, genomic 600..650). Verifies translation-order vs genomic-order
    bookkeeping on the negative strand."""
    run_mapping(binary, with_tags_index,
                "ENSP2\t1\t34\tMULTI_MINUS\n",
                tmp_path)
    assert_matches_golden(tmp_path, "test_minus_strand_multi_exon_domain",
                          GOLDEN_FILES, update=update_goldens)


def test_codon_split_across_exons(
    binary, with_tags_index, tmp_path, update_goldens,
):
    """1+2 and 2+1 codon-split cases on the *same* domain aa (aa 2..2)."""
    bed = (
        "ENSP4\t2\t2\tSPLIT_1PLUS2\n"
        "ENSP5\t2\t2\tSPLIT_2PLUS1\n"
    )
    run_mapping(binary, with_tags_index, bed, tmp_path)
    assert_matches_golden(tmp_path, "test_codon_split_across_exons",
                          GOLDEN_FILES, update=update_goldens)


def test_no_domain_mode(
    binary, with_tags_index, tmp_path, update_goldens,
):
    """Protein-only input (no aa range) emits the whole-transcript structure with
    overlap columns = NA and empty companion BEDs."""
    run_mapping(binary, with_tags_index,
                "ENSP1\n",
                tmp_path)
    # bed12 is empty in no-domain mode (allowed); included so the test catches
    # a regression that would suddenly emit blocks there.
    assert_matches_golden(tmp_path, "test_no_domain_mode",
                          GOLDEN_FILES, update=update_goldens)


def test_versioned_unversioned_ids(
    binary, with_tags_index, tmp_path, update_goldens,
):
    """ENSP1 and ENSP1.5 must produce identical mapping intervals. Captured in
    the goldens so any future divergence between the two id forms is flagged."""
    bed = (
        "ENSP1\t1\t10\tV_PLAIN\n"
        "ENSP1.5\t1\t10\tV_VERSIONED\n"
    )
    run_mapping(binary, with_tags_index, bed, tmp_path)
    assert_matches_golden(tmp_path, "test_versioned_unversioned_ids",
                          GOLDEN_FILES, update=update_goldens)


def test_enst_input_equals_ensp(out_all):
    """An ENST query maps to the same genomic interval as the equivalent ENSP
    query — only input_id_type differs."""
    summary = summary_by_id(out_all)
    ensp, enst = summary["Q1_ENSP"], summary["Q2_ENST"]
    assert ensp["input_id_type"] == "ENSP"
    assert enst["input_id_type"] == "ENST"
    assert enst["protein_id"] == "ENSP1"
    assert enst["domain_genomic_start"] == ensp["domain_genomic_start"]
    assert enst["domain_genomic_end"] == ensp["domain_genomic_end"]


def test_selenoprotein_like_cds_mismatch(out_all):
    """A CDS whose length isn't a multiple of 3 (selenoprotein-like) is still
    mapped, but flagged: cds_length_mismatch=true and status carries
    `_cds_mismatch`."""
    q7 = summary_by_id(out_all)["Q7_SEC"]
    assert q7["cds_length_mismatch"] == "true"
    assert q7["cds_nt_remainder"] == "1"
    assert "_cds_mismatch" in q7["status"]
