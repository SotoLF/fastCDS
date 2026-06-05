"""GTF-dialect and custom-protein compatibility.

Covers inputs that aren't the standard GENCODE/Ensembl-with-tags case:
  * a GTF that carries no tag attributes at all (flags must be NA, not false),
  * an NCBI RefSeq-style GTF (`gene` instead of `gene_name`, `locus_tag` that
    must not be mistaken for a real `tag`), and
  * proteins injected via scripts/append_custom_proteins.py.
"""

from __future__ import annotations

import subprocess
import sys

from conftest import REPO_ROOT, run_mapping, summary_by_id


def test_no_tags_gtf_reports_na_flags(binary, no_tags_index, tmp_path):
    """When the GTF has no tag attributes anywhere, is_mane_select and
    is_ensembl_canonical are NA — distinct from `false`, which means the GTF
    has tags but this transcript carries none."""
    run_mapping(binary, no_tags_index, "ENSP7\t1\t10\tQ_NOTAGS\n", tmp_path)
    qn = summary_by_id(tmp_path)["Q_NOTAGS"]
    assert qn["is_mane_select"] == "NA"
    assert qn["is_ensembl_canonical"] == "NA"


def test_refseq_style_gtf(binary, tmp_path):
    """NCBI RefSeq GTFs put the gene symbol in `gene "X"` (not `gene_name`) and
    carry `locus_tag "..."` but no real `tag "..."`. The parser must read the
    symbol from `gene`, strip the protein version, and NOT mistake `locus_tag`
    for a MANE/canonical tag (so the flags stay NA)."""
    refseq_gtf = tmp_path / "refseq_style.gtf"
    refseq_gtf.write_text(
        'chrZ\tRefSeq\ttranscript\t100\t199\t.\t+\t.\t'
        'gene_id "YEAST1"; transcript_id "NM_000001.1"; '
        'gbkey "mRNA"; gene "FOO"; locus_tag "YEAST1";\n'
        'chrZ\tRefSeq\texon\t100\t199\t.\t+\t.\t'
        'gene_id "YEAST1"; transcript_id "NM_000001.1"; '
        'gene "FOO"; locus_tag "YEAST1"; exon_number "1";\n'
        'chrZ\tRefSeq\tCDS\t100\t198\t.\t+\t0\t'
        'gene_id "YEAST1"; transcript_id "NM_000001.1"; '
        'gene "FOO"; locus_tag "YEAST1"; protein_id "NP_000001.1"; '
        'exon_number "1";\n'
    )
    refseq_idx = tmp_path / "refseq_style.idx"
    subprocess.run([str(binary), "index", "--gtf", str(refseq_gtf),
                    "--out", str(refseq_idx)],
                   check=True, capture_output=True, text=True)

    run_mapping(binary, refseq_idx, "NP_000001.1\t1\t10\tRSQ_TEST\n", tmp_path)
    rs = summary_by_id(tmp_path)["RSQ_TEST"]
    assert rs["protein_id"] == "NP_000001"          # version stripped
    assert rs["gene_name"] == "FOO"                 # symbol from `gene`
    assert rs["is_mane_select"] == "NA"             # locus_tag != real tag
    assert rs["is_ensembl_canonical"] == "NA"
    assert rs["status"] == "ok"


def test_custom_protein_injection(binary, synthetic_gtfs, tmp_path):
    """scripts/append_custom_proteins.py turns a TSV of genomic blocks into GTF
    rows with strand-aware exon numbering; after indexing the augmented GTF the
    custom protein maps like any reference one."""
    custom_tsv = tmp_path / "custom.tsv"
    custom_tsv.write_text(
        "protein_id\ttranscript_id\tgene_id\tgene_name\tchrom\tstrand\tblocks\tsource\n"
        "SYN_NEG\tSYN_NEG_tx\tSYN_NEG_g\tSYN_NEG\tchrSYN\t-\t500-700;800-1000\ttest\n"
    )
    custom_rows = tmp_path / "custom_rows.gtf"
    helper = REPO_ROOT / "scripts" / "append_custom_proteins.py"
    proc = subprocess.run(
        [sys.executable, str(helper), "--in", str(custom_tsv)],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stderr
    custom_rows.write_text(proc.stdout)   # the script emits the GTF rows on stdout

    rows = [l for l in custom_rows.read_text().splitlines() if not l.startswith("#")]
    assert len(rows) == 5  # transcript + 2 exons + 2 CDSs
    # On the '-' strand the high-coordinate exon (800-1000) is exon_number 1.
    exon1 = [l for l in rows if "\texon\t" in l and 'exon_number "1"' in l]
    assert sum("\t800\t" in l for l in exon1) == 1

    combined = tmp_path / "with_custom.gtf"
    combined.write_text(synthetic_gtfs["with_tags"].read_text() + custom_rows.read_text())
    custom_idx = tmp_path / "with_custom.idx"
    subprocess.run([str(binary), "index", "--gtf", str(combined),
                    "--out", str(custom_idx)],
                   check=True, capture_output=True, text=True)

    run_mapping(binary, custom_idx, "SYN_NEG\t10\t30\tSYN_DOM\n", tmp_path)
    cu = summary_by_id(tmp_path)["SYN_DOM"]
    assert cu["protein_id"] == "SYN_NEG"
    assert cu["chrom"] == "chrSYN"
    assert cu["strand"] == "-"
    assert cu["status"] == "ok"
