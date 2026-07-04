"""BED12 output and output-mode equivalences.

Checks the BED12 block geometry, the `--bed12` add-on flag on a non-bed12
output mode, and that streaming (`--batch-size`) is byte-for-byte identical to
the one-shot run.
"""

from __future__ import annotations

import subprocess

from conftest import BED1


def test_bed12_block_geometry(out_all):
    """One BED12 row per successful domain query (7 for BED1), with correct
    chrom/start/end and block columns for the single- and split-exon cases."""
    bed12 = (out_all / "domain_blocks.bed").read_text().splitlines()
    assert len(bed12) == 7

    q1 = [l for l in bed12 if "Q1_ENSP" in l]
    assert len(q1) == 1
    parts = q1[0].split("\t")
    # ENSP1 aa1..10 = CDS_1, genomic 120..149 → chromStart 119, one 30-nt block.
    assert parts[0] == "chrA"
    assert parts[1] == "119"
    assert parts[2] == "149"
    assert parts[9] == "1"        # blockCount
    assert parts[10] == "30,"     # blockSizes
    assert parts[11] == "0,"      # blockStarts

    # Q8 is a 1+2 codon split → two blocks.
    q8 = [l for l in bed12 if "Q8_SPLIT12" in l][0].split("\t")
    assert q8[9] == "2"


def test_coding_plus_bed12_addon(binary, with_tags_index, out_all, tmp_path):
    """`--output coding --bed12` adds the BED12 alongside the coding TSV without
    emitting the isoform/intron tables, and the BED12 is identical to the one
    produced under `--output all`."""
    out_b12 = tmp_path / "out_coding_bed12"
    bed = tmp_path / "queries.bed"
    bed.write_text(BED1)
    subprocess.run(
        [str(binary), "map", "--index", str(with_tags_index), "--bed", str(bed),
         "--out-dir", str(out_b12), "--output", "coding", "--bed12"],
        check=True, capture_output=True, text=True,
    )
    assert (out_b12 / "domain_cds_segments.tsv").exists()
    assert (out_b12 / "domain_blocks.bed").exists()
    assert not (out_b12 / "isoform_structure.tsv").exists()
    assert ((out_b12 / "domain_blocks.bed").read_bytes()
            == (out_all / "domain_blocks.bed").read_bytes())


def test_batch_size_equivalence(binary, with_tags_index, out_all, tmp_path):
    """Streaming with a small `--batch-size` produces byte-identical output
    files (and an equivalent run_metadata.json modulo cli/timestamp)."""
    out_batched = tmp_path / "out_all_batched"
    bed = tmp_path / "queries.bed"
    bed.write_text(BED1)
    subprocess.run(
        [str(binary), "map", "--index", str(with_tags_index), "--bed", str(bed),
         "--out-dir", str(out_batched), "--output", "all", "--batch-size", "3"],
        check=True, capture_output=True, text=True,
    )
    for fname in (
        "domain_mapping_summary.tsv", "domain_cds_segments.tsv",
        "domain_cds_segments.bed", "domain_introns.tsv", "domain_introns.bed",
        "domain_span_with_introns.bed", "isoform_structure.tsv",
        "domain_blocks.bed", "unmapped_domains.tsv",
    ):
        a = (out_all / fname).read_bytes() if (out_all / fname).exists() else b""
        b = (out_batched / fname).read_bytes() if (out_batched / fname).exists() else b""
        assert a == b, f"{fname} differs between one-shot and batched runs"

    def canonical(p):
        return "\n".join(line for line in p.read_text().splitlines()
                         if '"timestamp_utc"' not in line and '"cli"' not in line)

    assert (canonical(out_all / "run_metadata.json")
            == canonical(out_batched / "run_metadata.json"))
