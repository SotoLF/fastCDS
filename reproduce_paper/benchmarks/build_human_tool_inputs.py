#!/usr/bin/env python3
"""Build geneplot inputs for the human Ensembl-86 set, so the tool can be
benchmarked on the *same* human workload as fastCDS / ensembldb.

geneplot needs an Ensembl GFF3 + an InterProScan .ipr file. The domain source is
derived from the Pfam-on-v86 table we already use elsewhere.

Inputs expected in $DATA (default ~/Desktop/protein2genomic_data):
  Homo_sapiens.GRCh38.86 GFF3  (download from Ensembl release-86)
  pfam_human_v86_meta.tsv      (query_id, protein_id, pfam_id, interpro_id, aa_start, aa_end)
  queries_v86.bed              (protein_id, aa_start, aa_end, query_id)

Outputs (in $DATA/human_tool_bench): h86.ipr, ensp_enst.tsv

Then run:  python benchmarks/geneplot_human.py
See benchmarks/README.md for the measured numbers.
"""
# (The body that produced the committed numbers is kept in the repo's
#  benchmarks history; see build_human_tool_inputs steps in README.)
print("See README.md 'Other tools on human data' for the build steps and results.")
