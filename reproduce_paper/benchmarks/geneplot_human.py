"""geneplot on the human set (offline benchmark harness).

Times the geneplot mapper over the human Ensembl-86 query set to measure its
speed on real human data. Assumes the offline harness layout under
/tmp/p2gbench/ (geneplot checkout + inputs); see benchmarks/README.md.
"""
import time, sys, os, io, contextlib, logging
sys.path.insert(0, "/tmp/p2gbench/geneplot")
logging.getLogger("geneplot").setLevel(logging.ERROR)
import geneplot as gp

GFF = "h86.gff3"; IPR = "h86.ipr"
# ENSP -> ENST map
ensp2enst = {}
for line in open("ensp_enst.tsv"):
    p, t = line.split(); ensp2enst[p] = t
# first 1000 v86 query proteins (same set the other tools used)
queries = []
for i, line in enumerate(open("/home/goguxor/Desktop/protein2genomic_data/queries_v86.bed")):
    if i >= 1000: break
    queries.append(line.split("\t")[0])

# --- setup: build the gffutils SQLite db from the full human GFF3 (O(genome)) ---
if os.path.exists(GFF + ".db"): os.remove(GFF + ".db")
t0 = time.time()
gp.createGFFdb(GFF)
build_s = time.time() - t0
print(f"geneplot human gffutils db build: {build_s:.1f} s", flush=True)

with contextlib.redirect_stdout(io.StringIO()):
    g = gp.genome(GFF, iprfile=IPR, vcffiles="./")

# --- map each query protein's gene ---
mapped = 0; t0 = time.time()
for pid in queries:
    enst = ensp2enst.get(pid)
    if not enst: continue
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            gene = g.gene(mRNAid="transcript:"+enst, proteinid=pid)
            gene._proteindoms(IPR, pid)
            gene._transcriptpos_to_genomepos()
        mapped += 1
    except Exception:
        pass
map_s = time.time() - t0
print(f"geneplot human: mapped {mapped} genes in {map_s:.2f} s  ({mapped/map_s:.0f} genes/s mapping)", flush=True)
print(f"RESULT geneplot_human build_s={build_s:.1f} map_s={map_s:.2f} mapped={mapped} rate={mapped/map_s:.0f}", flush=True)
