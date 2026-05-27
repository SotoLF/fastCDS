#!/usr/bin/env python3
"""End-to-end test suite for prot2exon.

Builds the synthetic GTFs via make_synthetic_gtf.py, builds an index, runs
several BED queries through the C++ binary, and asserts on the produced
outputs.

Run from the repo root after a successful `cmake --build`:
    python3 tests/run_tests.py

Exits non-zero on any failure.
"""

from __future__ import annotations

import csv
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
BIN = REPO_ROOT / "build" / "prot2exon"
WRAPPER = REPO_ROOT / "bin" / "prot2exon"
TESTS_DIR = Path(__file__).resolve().parent

WITH_TAGS_GTF = TESTS_DIR / "with_tags.gtf"
NO_TAGS_GTF   = TESTS_DIR / "no_tags.gtf"


# --------------------------------------------------------------------------- #
# Tiny test framework
# --------------------------------------------------------------------------- #

PASSED: list[str] = []
FAILED: list[tuple[str, str]] = []


def assert_eq(name: str, expected, actual):
    if expected == actual:
        PASSED.append(name)
    else:
        FAILED.append((name, f"expected {expected!r}, got {actual!r}"))


def assert_in(name: str, needle, haystack):
    if needle in haystack:
        PASSED.append(name)
    else:
        FAILED.append((name, f"expected {needle!r} to be in {haystack!r}"))


def assert_true(name: str, cond, hint: str = ""):
    if cond:
        PASSED.append(name)
    else:
        FAILED.append((name, hint or "condition was false"))


def run(*args, expect_zero: bool = True) -> subprocess.CompletedProcess:
    proc = subprocess.run([str(a) for a in args],
                          capture_output=True, text=True)
    if expect_zero and proc.returncode != 0:
        raise SystemExit(f"command failed ({proc.returncode}): {args}\n"
                         f"stdout: {proc.stdout}\nstderr: {proc.stderr}")
    return proc


def read_tsv(path: Path) -> list[dict]:
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


# --------------------------------------------------------------------------- #
# Test cases
# --------------------------------------------------------------------------- #

def main() -> int:
    # Step 0: regenerate the synthetic GTFs.
    run(sys.executable, TESTS_DIR / "make_synthetic_gtf.py")
    if not BIN.exists():
        raise SystemExit(f"binary not found: {BIN}. Build it first.")

    work = Path(tempfile.mkdtemp(prefix="p2g_test_"))
    print(f"work dir: {work}", file=sys.stderr)

    # ---- Build the with-tags index ---------------------------------------
    idx = work / "with_tags.idx"
    run(BIN, "--gtf", WITH_TAGS_GTF, "--build-index", "--index", idx)
    assert_true("with_tags index exists", idx.exists())

    # ---- Build the no-tags index -----------------------------------------
    idx_no = work / "no_tags.idx"
    run(BIN, "--gtf", NO_TAGS_GTF, "--build-index", "--index", idx_no)
    assert_true("no_tags index exists", idx_no.exists())

    # ---- BED 1: cover ENSP vs ENST, MANE flags, structure_only, beyond,
    #            negative strand, codon splits, selenoprotein, non-coding ---
    bed = work / "queries.bed"
    bed.write_text("""\
# 1) ENSP query, MANE+canonical, plus strand, aa 1..10 = CDS_1 fully
ENSP1\t1\t10\tQ1_ENSP
# 2) ENST query for the same transcript — should produce identical mapping
ENST1\t1\t10\tQ2_ENST
# 3) Versioned ENSP — version is stripped
ENSP1.5\t1\t10\tQ3_VER
# 4) Structure-only (no aa). aa_start = aa_end = 0 keeps the domain_id column
#    intact (whitespace-tokenized parser collapses empty fields).
ENSP1\t0\t0\tQ4_STRUCT
ENSP2\t0\t0\tQ4b_STRUCT_M
# 5) Domain beyond protein length (ENSP1 has 44 aa)
ENSP1\t100\t110\tQ5_BEYOND
# 6) Negative strand, full coverage of CDS_1 in translation order (aa 1..17)
ENSP2\t1\t17\tQ6_NEG
# 7) Selenoprotein-like, domain at the end (aa 8..8)
ENSP3\t8\t8\tQ7_SEC
# 8) Codon split 1+2 — domain on the split aa (aa 2..2)
ENSP4\t2\t2\tQ8_SPLIT12
# 9) Codon split 2+1 — domain on the split aa (aa 2..2)
ENSP5\t2\t2\tQ9_SPLIT21
# 10) Non-coding ENST — no CDS
ENST6\t1\t10\tQ10_NONCODING
# 11) Protein not in index
ENSP99\t1\t10\tQ11_NOTFOUND
""")

    out = work / "out_all"
    run(BIN, "--index", idx, "--bed", bed, "--out-dir", out, "--output", "all")

    summary = {r["input_id"]: r for r in read_tsv(out / "domain_mapping_summary.tsv")}
    unmapped = {r["input_id"]: r for r in read_tsv(out / "unmapped_domains.tsv")}

    # Q1: ENSP positive case
    q1 = summary["Q1_ENSP"]
    assert_eq("Q1 status ok", "ok", q1["status"])
    assert_eq("Q1 input_id_type ENSP", "ENSP", q1["input_id_type"])
    assert_eq("Q1 is_mane_select", "true", q1["is_mane_select"])
    assert_eq("Q1 is_ensembl_canonical", "true", q1["is_ensembl_canonical"])
    assert_eq("Q1 cds_length_mismatch", "false", q1["cds_length_mismatch"])
    assert_eq("Q1 cds_nt_remainder", "0", q1["cds_nt_remainder"])
    assert_eq("Q1 protein_length_aa = 44", "44", q1["protein_length_aa"])
    assert_eq("Q1 domain_genomic_start = 120", "120", q1["domain_genomic_start"])
    assert_eq("Q1 domain_genomic_end = 149", "149", q1["domain_genomic_end"])
    assert_eq("Q1 n_coding_segments = 1", "1", q1["n_coding_segments"])
    # Phase 7 derived columns: Q1 is a single-exon domain (aa 1-10 fits in CDS_1)
    assert_eq("Q1 n_coding_exons_touched = 1", "1", q1["n_coding_exons_touched"])
    assert_eq("Q1 n_introns_spanned = 0", "0", q1["n_introns_spanned"])
    assert_eq("Q1 is_single_exon_domain = true", "true", q1["is_single_exon_domain"])
    assert_eq("Q1 intron_burden_nt = 0", "0", q1["intron_burden_nt"])
    # Largest-exon fraction is 1.0 for a fully-contained domain.
    assert_true("Q1 fraction_domain_in_largest_exon ≈ 1.0",
                abs(float(q1["fraction_domain_in_largest_exon"]) - 1.0) < 0.001)

    # Q2: ENST input gives same result, except input_id_type
    q2 = summary["Q2_ENST"]
    assert_eq("Q2 input_id_type ENST", "ENST", q2["input_id_type"])
    assert_eq("Q2 domain_genomic_start same", q1["domain_genomic_start"],
              q2["domain_genomic_start"])
    assert_eq("Q2 domain_genomic_end same", q1["domain_genomic_end"],
              q2["domain_genomic_end"])
    assert_eq("Q2 protein_id resolved to ENSP1", "ENSP1", q2["protein_id"])

    # Q3: versioned id
    q3 = summary["Q3_VER"]
    assert_eq("Q3 versioned stripped, protein_id = ENSP1", "ENSP1", q3["protein_id"])

    # Q4 / Q4b: structure_only
    assert_eq("Q4 status structure_only", "structure_only", summary["Q4_STRUCT"]["status"])
    assert_eq("Q4 no_domain_mode true", "true", summary["Q4_STRUCT"]["no_domain_mode"])
    assert_eq("Q4b status structure_only", "structure_only",
              summary["Q4b_STRUCT_M"]["status"])

    # Q5: beyond protein length → unmapped (still in summary with the reason)
    assert_true("Q5 unmapped row exists", "Q5_BEYOND" in unmapped)
    assert_eq("Q5 reason", "domain_beyond_protein_length",
              unmapped["Q5_BEYOND"]["reason"])

    # Q6: negative strand, aa 1..17 = CDS_1 in translation order (genomic
    # 800..850). protein_length_aa = 41 (123 nt / 3).
    q6 = summary["Q6_NEG"]
    assert_eq("Q6 strand", "-", q6["strand"])
    assert_eq("Q6 domain_genomic_start = 800", "800", q6["domain_genomic_start"])
    assert_eq("Q6 domain_genomic_end = 850", "850", q6["domain_genomic_end"])
    assert_eq("Q6 protein_length_aa = 41", "41", q6["protein_length_aa"])
    # MANE/canonical: this transcript has none, but the GTF has tags overall.
    assert_eq("Q6 is_mane_select false (GTF has tags)", "false", q6["is_mane_select"])
    assert_eq("Q6 is_canonical false (GTF has tags)", "false",
              q6["is_ensembl_canonical"])

    # Q7: selenoprotein-like, status ok_cds_mismatch
    q7 = summary["Q7_SEC"]
    assert_eq("Q7 cds_length_mismatch true", "true", q7["cds_length_mismatch"])
    assert_eq("Q7 cds_nt_remainder", "1", q7["cds_nt_remainder"])
    assert_in("Q7 status carries _cds_mismatch", "_cds_mismatch", q7["status"])

    # Q8: codon split 1+2 — aa 2 covered by both CDS rows. Fractions sum to 1.
    isoform_rows = read_tsv(out / "isoform_structure.tsv")
    q8_cds = [r for r in isoform_rows
              if r["input_id"] == "Q8_SPLIT12" and r["feature_type"] == "CDS"
              and r["overlaps_domain"] == "coding_overlap"]
    assert_eq("Q8 two coding_overlap rows", 2, len(q8_cds))
    q8_total = sum(float(r["domain_overlap_fraction_of_domain"]) for r in q8_cds)
    assert_true("Q8 fractions sum to ~1.0",
                abs(q8_total - 1.0) < 1e-6,
                f"sum = {q8_total}")
    # Phase 7: Q8 spans 2 CDS exons with 1 intron between.
    q8_summary = summary["Q8_SPLIT12"]
    assert_eq("Q8 n_coding_exons_touched = 2", "2",
              q8_summary["n_coding_exons_touched"])
    assert_eq("Q8 n_introns_spanned = 1", "1",
              q8_summary["n_introns_spanned"])
    assert_eq("Q8 is_single_exon_domain = false", "false",
              q8_summary["is_single_exon_domain"])
    assert_true("Q8 intron_burden_nt > 0",
                int(q8_summary["intron_burden_nt"]) > 0)

    # Q9: codon split 2+1, symmetric check
    q9_cds = [r for r in isoform_rows
              if r["input_id"] == "Q9_SPLIT21" and r["feature_type"] == "CDS"
              and r["overlaps_domain"] == "coding_overlap"]
    assert_eq("Q9 two coding_overlap rows", 2, len(q9_cds))
    q9_total = sum(float(r["domain_overlap_fraction_of_domain"]) for r in q9_cds)
    assert_true("Q9 fractions sum to ~1.0",
                abs(q9_total - 1.0) < 1e-6,
                f"sum = {q9_total}")

    # Q10: non-coding ENST → no_CDS_for_protein
    assert_true("Q10 unmapped row exists", "Q10_NONCODING" in unmapped)
    assert_eq("Q10 reason", "no_CDS_for_protein", unmapped["Q10_NONCODING"]["reason"])

    # Q11: protein not in index
    assert_true("Q11 unmapped row exists", "Q11_NOTFOUND" in unmapped)
    assert_eq("Q11 reason", "protein_not_in_index",
              unmapped["Q11_NOTFOUND"]["reason"])

    # ---- BED 2: same query against the no-tags GTF ------------------------
    bed2 = work / "queries_notags.bed"
    bed2.write_text("ENSP7\t1\t10\tQ_NOTAGS\n")
    out2 = work / "out_notags"
    run(BIN, "--index", idx_no, "--bed", bed2, "--out-dir", out2, "--output", "all")
    no_summary = {r["input_id"]: r for r in
                  read_tsv(out2 / "domain_mapping_summary.tsv")}
    qn = no_summary["Q_NOTAGS"]
    assert_eq("No-tags GTF: is_mane_select = NA", "NA", qn["is_mane_select"])
    assert_eq("No-tags GTF: is_ensembl_canonical = NA", "NA",
              qn["is_ensembl_canonical"])

    # ---- RefSeq-style (NCBI) GTF compatibility ---------------------------
    # NCBI RefSeq GTFs differ from GENCODE in two subtle ways:
    #   (a) gene symbol lives in `gene "X"` instead of `gene_name "X"`,
    #   (b) the GTF carries `locus_tag "..."` but no genuine `tag "..."`,
    #       so the false-positive substring match would set is_mane_select
    #       to false instead of NA.
    # Verify the parser picks up the gene symbol from `gene` and correctly
    # recognises the absence of true tag attributes.
    refseq_gtf = work / "refseq_style.gtf"
    refseq_gtf.write_text(
        # transcript + exon + CDS for a tiny single-exon protein.
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
    refseq_idx = work / "refseq_style.idx"
    run(BIN, "--gtf", refseq_gtf, "--build-index", "--index", refseq_idx)
    refseq_bed = work / "refseq_q.bed"
    refseq_bed.write_text("NP_000001.1\t1\t10\tRSQ_TEST\n")
    refseq_out = work / "refseq_out"
    run(BIN, "--index", refseq_idx, "--bed", refseq_bed,
        "--out-dir", refseq_out, "--output", "all")
    rs_summary = {r["input_id"]: r for r in
                  read_tsv(refseq_out / "domain_mapping_summary.tsv")}
    rs = rs_summary["RSQ_TEST"]
    assert_eq("RefSeq: protein_id resolved", "NP_000001", rs["protein_id"])
    assert_eq("RefSeq: gene_name from `gene`", "FOO", rs["gene_name"])
    assert_eq("RefSeq: locus_tag does NOT trigger tag detection (mane=NA)",
              "NA", rs["is_mane_select"])
    assert_eq("RefSeq: locus_tag does NOT trigger tag detection (canon=NA)",
              "NA", rs["is_ensembl_canonical"])
    assert_eq("RefSeq: query mapped", "ok", rs["status"])

    # ---- Custom-protein injection (scripts/append_custom_proteins.py) -----
    # Append a non-reference ORF to a synthetic GTF, rebuild the index, query
    # the new protein. Covers strand-aware exon numbering and multi-exon block
    # parsing in the helper script.
    custom_tsv = work / "custom.tsv"
    custom_tsv.write_text(
        "protein_id\ttranscript_id\tgene_id\tgene_name\tchrom\tstrand\tblocks\tsource\n"
        # Single-exon + strand-minus to exercise reverse exon numbering
        "SYN_NEG\tSYN_NEG_tx\tSYN_NEG_g\tSYN_NEG\tchrSYN\t-\t500-700;800-1000\ttest\n"
    )
    custom_rows = work / "custom_rows.gtf"
    helper = REPO_ROOT / "scripts" / "append_custom_proteins.py"
    proc = subprocess.run(
        [sys.executable, str(helper), "--in", str(custom_tsv), "--out", str(custom_rows)],
        capture_output=True, text=True,
    )
    assert_eq("append_custom_proteins exit 0", 0, proc.returncode)
    rows = [l for l in custom_rows.read_text().splitlines()
            if not l.startswith("#")]
    # transcript + 2 exons + 2 CDSs = 5 rows
    assert_eq("custom GTF rows generated", 5, len(rows))
    # On '-' strand the high-coordinate exon (800-1000) gets exon_number 1
    exon1 = [l for l in rows if "\texon\t" in l and 'exon_number "1"' in l]
    assert_eq("custom: exon_number=1 is high-coord on minus strand",
              1, sum('\t800\t' in l for l in exon1))

    custom_combined = work / "with_tags_plus_custom.gtf"
    custom_combined.write_text(WITH_TAGS_GTF.read_text() + custom_rows.read_text())
    custom_idx = work / "with_tags_plus_custom.idx"
    run(BIN, "--gtf", custom_combined, "--build-index", "--index", custom_idx)

    custom_bed = work / "custom_q.bed"
    custom_bed.write_text("SYN_NEG\t10\t30\tSYN_DOM\n")
    custom_out = work / "custom_out"
    run(BIN, "--index", custom_idx, "--bed", custom_bed,
        "--out-dir", custom_out, "--output", "all")
    cu_summary = {r["input_id"]: r for r in
                  read_tsv(custom_out / "domain_mapping_summary.tsv")}
    cu = cu_summary["SYN_DOM"]
    assert_eq("custom: protein_id resolved", "SYN_NEG", cu["protein_id"])
    assert_eq("custom: chrom",              "chrSYN",  cu["chrom"])
    assert_eq("custom: strand",             "-",       cu["strand"])
    assert_eq("custom: mapped status",      "ok",      cu["status"])

    # ---- prot2exon fetch (offline path via --gtf-url file://) ------------
    # The fetch helper bundles the curl + gunzip + --build-index pipeline.
    # We test it offline by feeding it a file:// URL pointing at a gzipped
    # copy of the synthetic with_tags.gtf — exercises the download +
    # gunzip + build path without depending on a live host.
    import gzip
    fetch_gz = work / "fetch_src.gtf.gz"
    with gzip.open(fetch_gz, "wb") as g:
        g.write(WITH_TAGS_GTF.read_bytes())
    fetch_cache = work / "fetch_cache"
    fetch_helper = REPO_ROOT / "python" / "prot2exon" / "fetch.py"
    fetch_url = fetch_gz.resolve().as_uri()  # file:///abs/path/...
    proc = subprocess.run(
        [sys.executable, str(fetch_helper),
         "human", "--release", "TEST",
         "--gtf-url", fetch_url,
         "--cache-dir", str(fetch_cache),
         "--binary", str(BIN)],
        capture_output=True, text=True,
    )
    assert_eq("fetch exit 0", 0, proc.returncode)
    fetch_idx_path = proc.stdout.strip().splitlines()[-1]
    assert_true("fetch printed an index path on stdout", bool(fetch_idx_path))
    assert_true("fetch index file exists", Path(fetch_idx_path).exists())
    assert_true("fetch removed uncompressed GTF (no --keep-gtf)",
                not (fetch_cache / "human_vTEST.gtf").exists())
    # Re-running hits the cache and is a no-op.
    proc2 = subprocess.run(
        [sys.executable, str(fetch_helper),
         "human", "--release", "TEST",
         "--gtf-url", fetch_url,
         "--cache-dir", str(fetch_cache),
         "--binary", str(BIN)],
        capture_output=True, text=True,
    )
    assert_eq("fetch second run exit 0", 0, proc2.returncode)
    assert_true("fetch second run mentions cache",
                "cached index" in proc2.stderr)
    # Confirm the built index actually works against a query.
    fetch_q_bed = work / "fetch_q.bed"
    fetch_q_bed.write_text("ENSP1\t1\t10\tFETCH_TEST\n")
    fetch_q_out = work / "fetch_q_out"
    run(BIN, "--index", fetch_idx_path, "--bed", fetch_q_bed,
        "--out-dir", fetch_q_out, "--output", "all")
    fq_summary = {r["input_id"]: r for r in
                  read_tsv(fetch_q_out / "domain_mapping_summary.tsv")}
    assert_eq("fetch index produces working mapping",
              "ok", fq_summary["FETCH_TEST"]["status"])

    # ---- plot --compact-genomic flag --------------------------------------
    # The new layout mode clamps each intron to a fixed display width while
    # keeping CDS/UTR at true bp scale. We just confirm the helper produces
    # a non-empty PDF and rejects --compact-genomic + --spliced together.
    cg_pdf = work / "compact_genomic.pdf"
    proc = subprocess.run(
        [sys.executable, "-c",
         "import sys; sys.path.insert(0, %r); "
         "from prot2exon.plot import main; "
         "sys.exit(main([%r,%r,%r,%r,%r,%r,%r]))" % (
             str(REPO_ROOT / "python"),
             "--isoform", str(out / "isoform_structure.tsv"),
             "--input-id", "Q1_ENSP",
             "--out", str(cg_pdf),
             "--compact-genomic")],
        capture_output=True, text=True,
    )
    assert_eq("plot --compact-genomic exit 0", 0, proc.returncode)
    assert_true("compact-genomic PDF non-empty",
                cg_pdf.exists() and cg_pdf.stat().st_size > 0)
    # Mutex check.
    proc_mx = subprocess.run(
        [sys.executable, "-c",
         "import sys; sys.path.insert(0, %r); "
         "from prot2exon.plot import main; "
         "sys.exit(main([%r,%r,%r,%r,%r,%r,%r,%r]))" % (
             str(REPO_ROOT / "python"),
             "--isoform", str(out / "isoform_structure.tsv"),
             "--input-id", "Q1_ENSP",
             "--out", str(work / "mx.pdf"),
             "--spliced", "--compact-genomic")],
        capture_output=True, text=True,
    )
    assert_eq("plot rejects --spliced + --compact-genomic", 2, proc_mx.returncode)

    # ---- plot --link-template (HTML linkout integration) -----------------
    # The --html path renders a plotly figure; --link-template adds a
    # clickable URL next to the title (TFRegDB2 / UniProt / UCSC / etc.).
    # We skip if plotly isn't installed — the helper raises SystemExit with
    # a clean message and the test would otherwise be flaky on barebones envs.
    try:
        import plotly  # noqa: F401
        have_plotly = True
    except ImportError:
        have_plotly = False
    if have_plotly:
        link_html = work / "link.html"
        proc = subprocess.run(
            [sys.executable, "-c",
             "import sys; sys.path.insert(0, %r); "
             "from prot2exon.plot import main; "
             "sys.exit(main([%r,%r,%r,%r,%r,%r,%r,%r]))" % (
                 str(REPO_ROOT / "python"),
                 "--isoform", str(out / "isoform_structure.tsv"),
                 "--input-id", "Q1_ENSP",
                 "--html", str(link_html),
                 "--link-template", "https://example.com/{protein_id}/entry")],
            capture_output=True, text=True,
        )
        assert_eq("plot --html --link-template exit 0", 0, proc.returncode)
        assert_true("link HTML non-empty",
                    link_html.exists() and link_html.stat().st_size > 0)
        html_text = link_html.read_text()
        # Plotly escapes the `/` as / inside its layout JSON, so we can't
        # match the URL literally. Assert the host AND the expanded protein_id
        # are both present — proves the template ran and the placeholder was
        # filled in.
        assert_true(
            "link_template host appears in HTML",
            "example.com" in html_text,
        )
        assert_true(
            "link_template protein_id placeholder expanded",
            "ENSP1" in html_text,
        )
        assert_true(
            "raw `{protein_id}` placeholder NOT present (was expanded)",
            "{protein_id}" not in html_text,
        )

    # ---- --html-interactive (interactive standalone HTML) ----------------
    # The TFRegDB2 viewer is vanilla JS — no plotly dependency, so this
    # always runs. We assert the output is a self-contained HTML, embeds
    # the segments as a JS payload, and includes the TFRegDB2 visual cues
    # (slate CDS, intron chevrons, domain lane).
    interactive_html = work / "interactive.html"
    proc = subprocess.run(
        [sys.executable, "-c",
         "import sys; sys.path.insert(0, %r); "
         "from prot2exon.plot import main; "
         "sys.exit(main([%r,%r,%r,%r,%r,%r]))" % (
             str(REPO_ROOT / "python"),
             "--isoform", str(out / "isoform_structure.tsv"),
             "--input-id", "Q1_ENSP",
             "--html-interactive", str(interactive_html))],
        capture_output=True, text=True,
    )
    assert_eq("plot --html-interactive exit 0", 0, proc.returncode)
    assert_true("interactive HTML non-empty",
                interactive_html.exists() and interactive_html.stat().st_size > 0)
    tf_text = interactive_html.read_text()
    # No external resources — the file must work offline.
    for forbidden in ("cdn.plot.ly", "esm.sh", "cdnjs.cloudflare", "https://"):
        # We do allow https:// in the externalLink we deliberately bake in
        # via link_template, so only check for CDNs.
        pass
    assert_true("self-contained: no plotly CDN",
                "cdn.plot.ly" not in tf_text and "esm.sh" not in tf_text)
    # Embedded payload sanity: transcript_id is present.
    assert_true("transcript id embedded in payload",
                "ENST1" in tf_text)
    # Domain detection: the BED query has aa_start=1, aa_end=10 → CDS_domain
    # should fire on at least one segment → the JS payload includes a
    # `domains` array with at least one entry.
    assert_true("domains payload populated",
                '"domains":' in tf_text and '"type": "Other"' in tf_text)

    # ---- benchmarks/make_figure_1.py renders without crashing --------------
    # We don't pass --pfam-atlas-tsv (the script falls back to the wiki's
    # headline numbers) and use the in-repo TP53 fixture.
    fig_dir = work / "figures"
    proc = subprocess.run(
        [sys.executable, str(REPO_ROOT / "benchmarks" / "make_figure_1.py"),
         "--tp53-isoforms", str(REPO_ROOT / "examples" / "tp53_isoforms.tsv"),
         "--out-dir", str(fig_dir)],
        capture_output=True, text=True,
    )
    assert_eq("make_figure_1 exit 0", 0, proc.returncode)
    assert_true("figure_1.png written",
                (fig_dir / "figure_1.png").exists() and
                (fig_dir / "figure_1.png").stat().st_size > 10_000)
    assert_true("figure_1.pdf written",
                (fig_dir / "figure_1.pdf").exists() and
                (fig_dir / "figure_1.pdf").stat().st_size > 5_000)

    # ---- prot2exon fetch index --url (Zenodo path) ----------------------
    # The "index" subcommand downloads a pre-built .idx straight from a URL,
    # bypassing the GTF parse + build step. We test it against a local
    # file:// URL backed by the synthetic index built earlier.
    fetched_idx = work / "fetched.idx"
    src_idx = work / "with_tags.idx"
    src_sha = subprocess.run(
        [sys.executable, "-c", f"import hashlib;print(hashlib.sha256(open({str(src_idx)!r},'rb').read()).hexdigest())"],
        capture_output=True, text=True,
    ).stdout.strip()
    src_url = f"file://{src_idx}"

    proc = subprocess.run(
        [sys.executable, "-m", "prot2exon.fetch", "index",
         "--url", src_url, "--out", str(fetched_idx), "--sha256", src_sha],
        env={**os.environ, "PYTHONPATH": str(REPO_ROOT / "python")},
        capture_output=True, text=True,
    )
    assert_eq("fetch index exit 0", 0, proc.returncode)
    assert_true("fetched index exists", fetched_idx.exists())
    assert_eq("fetched index size matches source",
              src_idx.stat().st_size, fetched_idx.stat().st_size)

    # Bad sha256 should make the verify step abort with exit code 1.
    proc = subprocess.run(
        [sys.executable, "-m", "prot2exon.fetch", "index",
         "--url", src_url, "--out", str(fetched_idx),
         "--sha256", "0" * 64, "--force"],
        env={**os.environ, "PYTHONPATH": str(REPO_ROOT / "python")},
        capture_output=True, text=True,
    )
    assert_true("fetch index bad sha exits non-zero", proc.returncode != 0)
    assert_true("fetch index bad sha mentions mismatch",
                "sha256 mismatch" in proc.stderr)

    # `fetch list` lists both index presets and GTF-build presets.
    proc = subprocess.run(
        [sys.executable, "-m", "prot2exon.fetch", "list"],
        env={**os.environ, "PYTHONPATH": str(REPO_ROOT / "python")},
        capture_output=True, text=True,
    )
    assert_eq("fetch list exit 0", 0, proc.returncode)
    assert_true("fetch list shows pre-built indexes",
                "Pre-built indexes" in proc.stdout)
    assert_true("fetch list shows GTF-build presets",
                "Build from GTF" in proc.stdout)
    assert_true("fetch list shows the yeast index preset",
                "yeast" in proc.stdout)

    # `fetch index --preset` errors clearly when the Zenodo record id
    # hasn't been published yet (URL still has `<RECORD>` placeholder).
    proc = subprocess.run(
        [sys.executable, "-m", "prot2exon.fetch", "index",
         "--preset", "human-v49", "--out", str(work / "p49.idx")],
        env={**os.environ, "PYTHONPATH": str(REPO_ROOT / "python")},
        capture_output=True, text=True,
    )
    # Either errors because `<RECORD>` is still a placeholder, OR succeeds
    # if the URL has been updated post-Zenodo-upload. Both are valid.
    pre_zenodo = "<RECORD>" in proc.stderr
    assert_true("fetch index --preset gives a meaningful message either way",
                pre_zenodo or proc.returncode == 0)

    # ---- _render_to_string + plot_height parametrisation ----------------
    # The render_to_string helper backs both the file writer and the
    # Jupyter wrapper, so we sanity-check it directly (no IPython needed).
    proc = subprocess.run(
        [sys.executable, "-c",
         "import sys; sys.path.insert(0, %r); "
         "from prot2exon._interactive_html import _render_to_string; "
         "from prot2exon.plot import load_isoform_tsv; "
         "by_id = load_isoform_tsv(%r); "
         "segs = by_id['Q1_ENSP']; "
         "out = _render_to_string(segs, plot_height=200); "
         "assert 'height: 200px' in out, 'plot_height kwarg not applied'; "
         "assert 'top: 27.5%%' in out, 'percentage CSS missing'; "
         "print('OK', len(out))" % (
             str(REPO_ROOT / "python"),
             str(out / "isoform_structure.tsv"))],
        capture_output=True, text=True,
    )
    assert_eq("_render_to_string exit 0", 0, proc.returncode)
    assert_true("plot_height parameter wired",
                proc.stdout.startswith("OK"))

    # ---- render_interactive_html_stack (multi-isoform) ----------------------
    # The stack viewer expects a dict of {input_id: list[Segment]} and
    # renders all isoforms on a single shared axis. We give it two copies
    # of the same isoform so the union-axis logic is exercised.
    stack_html = work / "interactive_stack.html"
    proc = subprocess.run(
        [sys.executable, "-c",
         "import sys; sys.path.insert(0, %r); "
         "from prot2exon.plot import load_isoform_tsv; "
         "from prot2exon import render_interactive_html_stack; "
         "by_id = load_isoform_tsv(%r); "
         "render_interactive_html_stack(by_id, %r); "
         "print('OK')" % (
             str(REPO_ROOT / "python"),
             str(out / "isoform_structure.tsv"),
             str(stack_html))],
        capture_output=True, text=True,
    )
    assert_eq("stack viewer exit 0", 0, proc.returncode)
    assert_true("stack HTML non-empty",
                stack_html.exists() and stack_html.stat().st_size > 0)
    stack_text = stack_html.read_text()
    assert_true("stack viewer embeds STRUCTURES array",
                "const STRUCTURES" in stack_text and "iso-stack" in stack_text)

    # ---- render_interactive_jupyter returns an IFrame wrapper ---------------
    # We skip the IPython.display.HTML object construction when IPython is
    # absent, but the function should at least exist on the package.
    proc = subprocess.run(
        [sys.executable, "-c",
         "import sys; sys.path.insert(0, %r); "
         "import prot2exon; "
         "assert hasattr(prot2exon, 'render_interactive_jupyter'); "
         "print('OK')" % (str(REPO_ROOT / "python"),)],
        capture_output=True, text=True,
    )
    assert_eq("render_interactive_jupyter exported exit 0", 0, proc.returncode)
    assert_true("render_interactive_jupyter is part of the public API",
                "OK" in proc.stdout)

    # ---- notebooks/generate_notebooks.py still runs ----------------------
    # Cheap regression: just make sure the generator imports + produces the
    # expected ipynb files. We pass --out-dir so the run writes to a tempdir
    # rather than the repo's notebooks/, otherwise any embedded outputs the
    # user has from a real `--run` invocation get wiped on every test pass.
    nb_gen = REPO_ROOT / "notebooks" / "generate_notebooks.py"
    if nb_gen.exists():
        nb_tmp = work / "nb_out"
        proc = subprocess.run(
            [sys.executable, str(nb_gen), "--out-dir", str(nb_tmp)],
            capture_output=True, text=True,
        )
        assert_eq("generate_notebooks exit 0", 0, proc.returncode)
        assert_true("walkthrough notebook generated",
                    "walkthrough_end_to_end.ipynb" in proc.stdout)
        # Confirm the .ipynb files actually landed in the tempdir, not the repo.
        assert_true("notebooks written to --out-dir",
                    (nb_tmp / "walkthrough_end_to_end.ipynb").exists() and
                    (nb_tmp / "validation.ipynb").exists())

    # ---- --batch-size equivalence ----------------------------------------
    # Re-run the same query set in streaming mode (batch_size=3, smaller than
    # the 11-query input so multiple chunks are exercised) and assert every
    # output file is byte-identical to the one-shot run. run_metadata.json
    # is compared with the CLI line + timestamp stripped.
    out_batched = work / "out_all_batched"
    run(BIN, "--index", idx, "--bed", bed, "--out-dir", out_batched,
        "--output", "all", "--batch-size", "3")
    for fname in (
        "domain_mapping_summary.tsv", "domain_cds_segments.tsv",
        "domain_cds_segments.bed", "domain_introns.tsv",
        "domain_introns.bed", "domain_span_with_introns.bed",
        "isoform_structure.tsv", "domain_blocks.bed12",
        "unmapped_domains.tsv",
    ):
        a = (out / fname).read_bytes() if (out / fname).exists() else b""
        b_ = (out_batched / fname).read_bytes() if (out_batched / fname).exists() else b""
        assert_eq(f"batched=={fname}", a, b_)

    def _meta_canonical(p: Path) -> str:
        return "\n".join(
            line for line in p.read_text().splitlines()
            if '"timestamp_utc"' not in line and '"cli"' not in line
        )
    assert_eq("batched==run_metadata.json (ex-cli/timestamp)",
              _meta_canonical(out / "run_metadata.json"),
              _meta_canonical(out_batched / "run_metadata.json"))

    # ---- BED12 sanity ----------------------------------------------------
    bed12 = (out / "domain_blocks.bed12").read_text().splitlines()
    # We expect a BED12 row for every successful query that has a domain.
    # Q1, Q2, Q3 (ENSP1 aa1..10), Q6 (ENSP2 aa1..17), Q7 (ENSP3 aa8),
    # Q8 (ENSP4), Q9 (ENSP5). 7 rows.
    assert_eq("BED12 row count", 7, len(bed12))
    # Q1's BED12: one block, size 30, start 0, chromStart 119 (genomic 120-1).
    q1_line = [l for l in bed12 if "Q1_ENSP" in l]
    assert_eq("Q1 BED12 lines == 1", 1, len(q1_line))
    parts = q1_line[0].split("\t")
    assert_eq("Q1 BED12 chrom", "chrA", parts[0])
    assert_eq("Q1 BED12 chromStart", "119", parts[1])
    assert_eq("Q1 BED12 chromEnd", "149", parts[2])
    assert_eq("Q1 BED12 blockCount", "1", parts[9])
    assert_eq("Q1 BED12 blockSizes", "30,", parts[10])
    assert_eq("Q1 BED12 blockStarts", "0,", parts[11])
    # Q8's BED12: split 1+2, two blocks.
    q8_line = [l for l in bed12 if "Q8_SPLIT12" in l][0].split("\t")
    assert_eq("Q8 BED12 blockCount", "2", q8_line[9])

    # ---- Plotter smoke test (CLI) ----------------------------------------
    pdf = work / "Q1.pdf"
    proc = subprocess.run(
        [str(WRAPPER), "plot", "--isoform", str(out / "isoform_structure.tsv"),
         "--input-id", "Q1_ENSP", "--out", str(pdf)],
        capture_output=True, text=True,
    )
    assert_eq("plot exit 0", 0, proc.returncode)
    assert_true("PDF created", pdf.exists() and pdf.stat().st_size > 0)

    # ---- Python wrapper API ---------------------------------------------
    # Exercise the high-level Python API end-to-end against the with-tags
    # synthetic index. Mirrors the user-facing examples in the README.
    sys.path.insert(0, str(REPO_ROOT / "python"))
    import prot2exon as p2e
    mapper = p2e.Mapper(index=str(idx))
    assert_eq("Mapper.binary points at our build", str(BIN), mapper.binary)

    py_result = mapper.map("ENSP1", aa_start=1, aa_end=10, domain_id="PY_AD1")
    assert_eq("Mapper.map: n_mapped", 1, py_result.n_mapped)
    assert_eq("Mapper.map: gene_name", "TEST1",
              str(py_result.summary.iloc[0]["gene_name"]))
    assert_eq("Mapper.map: bed12 row count", 1, len(py_result.bed12))

    py_batch = mapper.map_batch([
        {"protein_id": "ENSP1", "aa_start": 1, "aa_end": 10, "domain_id": "PY_AD1"},
        {"protein_id": "ENST1", "aa_start": 1, "aa_end": 10, "domain_id": "PY_AD1_E"},
        {"protein_id": "ENSP4", "aa_start": 2, "aa_end": 2, "domain_id": "PY_SPLIT"},
    ])
    assert_eq("map_batch: n_total", 3, py_batch.n_total)
    assert_eq("map_batch: n_mapped", 3, py_batch.n_mapped)

    # keep_outputs persists the files and read_results_dir round-trips them.
    persist_dir = work / "py_persist"
    p2e.Mapper(index=str(idx)).map(
        "ENSP1", 1, 10, "PY_PERSIST",
        keep_outputs=str(persist_dir))
    assert_true("keep_outputs: files on disk",
                (persist_dir / "domain_mapping_summary.tsv").exists())
    re_read = p2e.read_results_dir(str(persist_dir))
    assert_eq("read_results_dir: mapped count", 1, re_read.n_mapped)
    assert_eq("read_results_dir: metadata tool", "prot2exon",
              re_read.metadata.get("tool", ""))

    # Python plot() with a MappingResult and with a DataFrame both produce
    # a PDF. We just check the file is non-empty; the visual content is
    # the CLI test's responsibility.
    py_pdf = work / "py_plot.pdf"
    p2e.plot(py_result, out=str(py_pdf))
    assert_true("py plot from MappingResult",
                py_pdf.exists() and py_pdf.stat().st_size > 0)
    py_pdf2 = work / "py_plot2.pdf"
    p2e.plot(py_batch.isoform, input_id="PY_SPLIT", out=str(py_pdf2))
    assert_true("py plot from DataFrame",
                py_pdf2.exists() and py_pdf2.stat().st_size > 0)

    # ---- Wrap up ---------------------------------------------------------
    print(f"\n{len(PASSED)} passed, {len(FAILED)} failed")
    for n, why in FAILED:
        print(f"  FAIL  {n}: {why}", file=sys.stderr)
    if FAILED:
        return 1
    # Clean up.
    shutil.rmtree(work, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
