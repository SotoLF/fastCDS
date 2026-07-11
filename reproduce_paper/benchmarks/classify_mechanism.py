#!/usr/bin/env python3
"""Classify *why* a domain is not conserved in a target isoform.

Pure functions (no I/O) imported by isoform_conservation.py. Given a domain's
genomic coding intervals (from its source isoform) and a target isoform's full
feature structure (CDS / 5'UTR / 3'UTR / intron, genomic coords), decide the
splicing mechanism behind partial / lost conservation:

  conserved        inframe_coverage >= CONSERVED_T
  frame_disruption bases are present as CDS but in a shifted reading frame
  exon_skipping    domain bases fall in a target intron (spliced out)
  alt_N_terminal   domain bases fall in target 5'UTR or upstream of its CDS
  alt_C_terminal   domain bases fall in target 3'UTR / past its CDS (alt polyA)
  (terminal calls also cover bases beyond the target transcript bounds)

Classes/thresholds match the plan: conserved >=80%, partial 50-80%, lost <50%
on in-frame coverage.
"""
from __future__ import annotations

CONSERVED_T = 0.80   # inframe coverage at/above this -> conserved
PARTIAL_T = 0.50     # between PARTIAL_T and CONSERVED_T -> partial; below -> lost


def conservation_class(inframe_coverage: float) -> str:
    if inframe_coverage >= CONSERVED_T:
        return "conserved"
    if inframe_coverage >= PARTIAL_T:
        return "partial"
    return "lost"


def _overlap_bp(a0, a1, b0, b1) -> int:
    """Inclusive-coordinate overlap length of [a0,a1] and [b0,b1]."""
    lo, hi = max(a0, b0), min(a1, b1)
    return hi - lo + 1 if hi >= lo else 0


def classify_mechanism(domain_intervals, target, coverage, inframe_coverage):
    """Return the mechanism string for one (domain, target isoform) pair.

    domain_intervals: list of (g0, g1) genomic coding intervals of the domain
                      (inclusive, ascending), from the SOURCE isoform.
    target: dict with genomic feature lists of the TARGET isoform:
            cds   -> [(g0,g1), ...]   five_utr -> [...]   three_utr -> [...]
            intron-> [...]            tx_start, tx_end (transcript bounds)
            strand -> '+' / '-'
    """
    if inframe_coverage >= CONSERVED_T:
        return "conserved"
    # bases present but out of frame dominate -> frame disruption
    if coverage >= CONSERVED_T and inframe_coverage < CONSERVED_T:
        return "frame_disruption"

    # tally the domain bases that are NOT coding in the target, by context
    skip = nterm = cterm = 0
    strand = target["strand"]
    tx0, tx1 = target["tx_start"], target["tx_end"]
    for d0, d1 in domain_intervals:
        # remove the part covered by target CDS; classify the remainder
        coding = sum(_overlap_bp(d0, d1, c0, c1) for c0, c1 in target["cds"])
        missing = (d1 - d0 + 1) - coding
        if missing <= 0:
            continue
        intron = sum(_overlap_bp(d0, d1, i0, i1) for i0, i1 in target["intron"])
        five = sum(_overlap_bp(d0, d1, u0, u1) for u0, u1 in target["five_utr"])
        three = sum(_overlap_bp(d0, d1, u0, u1) for u0, u1 in target["three_utr"])
        # bases beyond the target transcript entirely -> terminal, side by strand
        outside_lo = max(0, tx0 - d0)            # genomic-left overhang
        outside_hi = max(0, d1 - tx1)            # genomic-right overhang
        # genomic-left is N-terminal on '+' strand, C-terminal on '-' strand
        if strand == "+":
            nterm += five + outside_lo
            cterm += three + outside_hi
        else:
            nterm += three + outside_hi
            cterm += five + outside_lo
        skip += intron

    tally = {"exon_skipping": skip, "alt_N_terminal": nterm, "alt_C_terminal": cterm}
    if max(tally.values()) == 0:
        return "frame_disruption" if coverage > inframe_coverage else "lost_other"
    return max(tally, key=tally.get)
