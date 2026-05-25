"""Stratify the Pfam-A atlas (Phase 4 Notebook 1 output) by functional class.

Classifies each Pfam family into one of:
    dna_binding, catalytic, transmembrane, signaling, structural, other

Classification is keyword-based on the Pfam family NAME and DESCRIPTION (from
Pfam-A.clans.tsv). This is a deliberately simple, reviewable heuristic — for
the paper's purposes we only need a coarse split to ask whether encoding
architecture (single- vs multi-exon, fraction-in-largest, intron burden)
differs across functional classes.

Inputs:
  --atlas       pfam_atlas.tsv     (137K rows, output of Notebook 1)
  --pfam-clans  Pfam-A.clans.tsv   (5 cols: pfam_acc, clan, clan_name, name, desc)

Outputs:
  --out-tsv     pfam_atlas_stratified.tsv   (atlas + 'domain_class' column)
  --out-summary pfam_atlas_class_summary.tsv (per-class headline stats)
  --out-figure  pfam_atlas_class_figure.png  (3-panel: % single-exon, median exons,
                                              median intron burden, by class)
"""

from __future__ import annotations

import argparse
import csv
import statistics
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


CLASS_ORDER = ["dna_binding", "catalytic", "transmembrane",
               "signaling", "structural", "other"]

# Keywords are searched (case-insensitive) against both 'name' and 'description'.
# The classifier applies them in CLASS_ORDER; first match wins.
KEYWORDS: dict[str, list[str]] = {
    "dna_binding": [
        "zinc finger", "zinc_finger", "zinc-finger", "homeobox", "homeodomain",
        "helix-loop-helix", "hlh", "bzip", "leucine zipper", "leucine_zipper",
        "dna binding", "dna-binding", "dna_binding", "transcription factor",
        "ets-", " ets ", "t-box", "tbx", "mads", "pou domain", "sox-tcf",
        "paired domain", "forkhead", "winged helix", "winged-helix",
        " hmg ", "hmg-", "hmg_box", "basic helix", "hth_", "helix-turn-helix",
        "rrm", "kh domain", "kh_", " znf ", "znf-", "krab", "scan domain",
    ],
    "catalytic": [
        "kinase", "phosphatase", "hydrolase", "transferase", "oxidase",
        "reductase", "isomerase", "ligase", "peptidase", "protease",
        "proteinase", "dehydrogenase", "synthase", "synthetase", "esterase",
        " lyase", "polymerase", "helicase", "nuclease", "dehydratase",
        "deaminase", "decarboxylase", "carboxylase", "racemase", "methylase",
        "methyltransferase", "acetyltransferase", "glycosyltransferase",
        "glycosylase", "topoisomerase", "epimerase", "mutase", "atpase",
        "gtpase", "phospholipase", " sulfatase", "amidase", "aminotransferase",
        "ubiquitin ligase", "ubiquitin-conjugating", "ring finger", "ring_finger",
    ],
    "transmembrane": [
        "transmembrane", "7tm", "gpcr", "ion channel", "ion_channel",
        " transporter", "transporter ", "pore-forming", " porin", "abc transporter",
        " mfs ", "channel", "permease", "antiporter", "symporter", "uniporter",
        "aquaporin", "tetraspanin", " claudin", "connexin", "membrane attack",
    ],
    "signaling": [
        " sh2", " sh3", " pdz", " ph domain", "pleckstrin", " c2 domain",
        " ww domain", " fha ", "brct", " btb", " bro1", " mit ", "death domain",
        "death-effector", " tir ", "pyrin", "card domain", " card ", "irs",
        " ankyrin repeat", "armadillo", "g-protein", "ras-binding", "rho gtpase",
        "rhogef", "rhogap", "rasgef", "rasgap", "guanine nucleotide",
    ],
    "structural": [
        "collagen", "immunoglobulin", " ig domain", "fibronectin", " egf",
        "egf-like", "cadherin", "laminin", "ankyrin", " wd40", "wd_40",
        " kelch", "spectrin", "leucine-rich", "leucine rich", "lrr",
        "tetratricopeptide", " tpr ", "armadillo repeat", "annexin",
        "actin", "tubulin", "myosin", "tropomyosin", "keratin", "intermediate filament",
        "fibrillin", "elastin", "lectin", "c-type lectin", "fibrinogen",
        "extracellular", "scaffold", "structural",
    ],
}


def classify_one(name: str, desc: str) -> str:
    text = f" {name.lower()} {desc.lower()} "
    for cls in CLASS_ORDER[:-1]:
        for kw in KEYWORDS[cls]:
            if kw in text:
                return cls
    return "other"


def load_pfam_names(path: Path) -> dict[str, tuple[str, str]]:
    """pfam_acc -> (name, description)."""
    out: dict[str, tuple[str, str]] = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            acc, _clan, _clan_name, name, desc = parts[:5]
            out[acc] = (name, desc)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--atlas", required=True, type=Path)
    ap.add_argument("--pfam-clans", required=True, type=Path)
    ap.add_argument("--out-tsv", required=True, type=Path)
    ap.add_argument("--out-summary", required=True, type=Path)
    ap.add_argument("--out-figure", required=True, type=Path)
    args = ap.parse_args()

    names = load_pfam_names(args.pfam_clans)
    print(f"loaded {len(names):,} Pfam families")

    # Classify every Pfam family once (cache).
    pfam_class: dict[str, str] = {acc: classify_one(n, d)
                                  for acc, (n, d) in names.items()}
    cls_counts_families = Counter(pfam_class.values())
    print("Pfam families per class:")
    for c in CLASS_ORDER:
        print(f"  {c:14s}  {cls_counts_families[c]:6,}")

    # Walk the atlas, attach class, collect per-class metrics.
    per_class: dict[str, dict[str, list]] = defaultdict(
        lambda: {"exons": [], "fraction": [], "intron": [], "single": 0, "n": 0})
    n_total = 0
    n_unclassified = 0

    with open(args.atlas) as fin, open(args.out_tsv, "w") as fout:
        rdr = csv.DictReader(fin, delimiter="\t")
        fieldnames = rdr.fieldnames + ["domain_class", "pfam_name"]
        wtr = csv.DictWriter(fout, fieldnames=fieldnames, delimiter="\t")
        wtr.writeheader()
        for row in rdr:
            n_total += 1
            pf = row["pfam_id"]
            cls = pfam_class.get(pf, "other")
            name = names.get(pf, ("", ""))[0]
            if pf not in pfam_class:
                n_unclassified += 1
            row["domain_class"] = cls
            row["pfam_name"] = name
            wtr.writerow(row)
            bucket = per_class[cls]
            bucket["n"] += 1
            bucket["exons"].append(int(row["n_coding_exons_touched"]))
            bucket["fraction"].append(float(row["fraction_in_largest"]))
            bucket["intron"].append(int(row["intron_burden_nt"]))
            if row["is_single_exon"] == "True":
                bucket["single"] += 1

    print(f"\natlas rows: {n_total:,}  (unmapped to known Pfam name: {n_unclassified:,})")

    # Per-class summary.
    with open(args.out_summary, "w") as f:
        f.write("class\tn_domain_instances\tn_pfam_families\tpct_single_exon"
                "\tmedian_exons\tmean_exons\tmedian_fraction_in_largest"
                "\tmedian_intron_burden\tmean_intron_burden\n")
        for c in CLASS_ORDER:
            b = per_class[c]
            n = b["n"]
            if n == 0:
                continue
            pct_single = 100.0 * b["single"] / n
            med_exons = statistics.median(b["exons"])
            mean_exons = statistics.mean(b["exons"])
            med_frac = statistics.median(b["fraction"])
            med_intron = statistics.median(b["intron"])
            mean_intron = statistics.mean(b["intron"])
            f.write(f"{c}\t{n}\t{cls_counts_families[c]}\t{pct_single:.2f}"
                    f"\t{med_exons}\t{mean_exons:.2f}\t{med_frac:.3f}"
                    f"\t{med_intron}\t{mean_intron:.1f}\n")

    print(f"wrote {args.out_summary}")
    print(open(args.out_summary).read())

    # Figure: 3 panels — % single-exon, median exons, median intron burden, by class.
    classes_present = [c for c in CLASS_ORDER if per_class[c]["n"] > 0]
    labels = [c.replace("_", " ") for c in classes_present]
    pct_single = [100.0 * per_class[c]["single"] / per_class[c]["n"]
                  for c in classes_present]
    med_exons = [statistics.median(per_class[c]["exons"]) for c in classes_present]
    med_intron_log = [statistics.median(per_class[c]["intron"]) for c in classes_present]

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    colors = ["#3b7dd8", "#d83b3b", "#3bd87a", "#d8a83b", "#7a3bd8", "#888888"]
    colors = colors[:len(classes_present)]

    axes[0].bar(labels, pct_single, color=colors)
    axes[0].set_ylabel("% single-exon domains")
    axes[0].set_title("Encoded by exactly 1 CDS exon")
    axes[0].axhline(27.3, ls="--", color="black", lw=0.8, alpha=0.6,
                    label="overall: 27.3%")
    axes[0].legend(loc="upper right", fontsize=8)

    axes[1].bar(labels, med_exons, color=colors)
    axes[1].set_ylabel("Median CDS exons per domain")
    axes[1].set_title("Median exon count")
    axes[1].axhline(2, ls="--", color="black", lw=0.8, alpha=0.6,
                    label="overall: 2")
    axes[1].legend(loc="upper right", fontsize=8)

    axes[2].bar(labels, med_intron_log, color=colors)
    axes[2].set_ylabel("Median intron burden (nt)")
    axes[2].set_yscale("log")
    axes[2].set_title("Median intron burden (log)")
    axes[2].axhline(1833, ls="--", color="black", lw=0.8, alpha=0.6,
                    label="overall: 1,833")
    axes[2].legend(loc="upper right", fontsize=8)

    for ax in axes:
        for tick in ax.get_xticklabels():
            tick.set_rotation(35)
            tick.set_ha("right")

    fig.suptitle("Pfam-A proteome atlas — encoding architecture by domain class",
                 fontsize=12)
    fig.tight_layout()
    fig.savefig(args.out_figure, dpi=150, bbox_inches="tight")
    print(f"wrote {args.out_figure}")


if __name__ == "__main__":
    main()
