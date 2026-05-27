"""Parse HMMER `--domtblout` (Pfam scan) into prot2exon BED-like rows.

Extracted from `scripts/prepare_from_pfam.py` so the same logic is usable
from both the CLI and the Python API.
"""

from __future__ import annotations

import re
import sys

from ._mapping import (UniProtToEnsp, looks_like_ensp, looks_like_uniprot,
                       open_text, strip_version)


def parse(path: str, *, mode: str, min_score: float, min_length: int,
          id_type: str, mapping: UniProtToEnsp | None
          ) -> tuple[list, list]:
    """Return ``(rows, rejected)``.

    rows: ``[(ensp, aa_start, aa_end, domain_id, description), ...]``
    rejected: ``[(raw_line, reason), ...]``
    """
    rows: list[tuple[str, int, int, str, str]] = []
    rejected: list[tuple[str, str]] = []

    if id_type == "auto":
        samples: list[str] = []
        with open_text(path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip(): continue
                fields = line.split()
                protein_col = 3 if mode == "scan" else 0
                if len(fields) > protein_col:
                    samples.append(fields[protein_col])
                if len(samples) >= 50: break
        n_ensp = sum(1 for s in samples if looks_like_ensp(s))
        n_up = sum(1 for s in samples if looks_like_uniprot(s))
        id_type = "ensp" if n_ensp >= max(n_up, 1) else (
                  "uniprot" if n_up else "unknown")
        print(f"[pfam] auto-detected id_type = {id_type}", file=sys.stderr)

    if id_type == "uniprot" and mapping is None:
        raise ValueError(
            "domtblout uses UniProt accessions but no mapping given.")

    with open_text(path) as fh:
        for ln, line in enumerate(fh, start=1):
            if line.startswith("#") or not line.strip():
                continue
            fields = re.split(r"\s+", line.rstrip("\n"), maxsplit=22)
            if len(fields) < 22:
                rejected.append((line.rstrip("\n"), "too few columns"))
                continue
            if mode == "scan":
                hmm_name = fields[0]
                hmm_acc  = fields[1]
                protein  = fields[3]
            else:
                protein  = fields[0]
                hmm_name = fields[3]
                hmm_acc  = fields[4]
            try:
                this_score = float(fields[13])
                ali_from   = int(fields[17])
                ali_to     = int(fields[18])
            except (ValueError, IndexError):
                rejected.append((line.rstrip("\n"), "unparseable numerics"))
                continue
            description = fields[22] if len(fields) > 22 else ""

            if this_score < min_score:
                continue
            if ali_to < ali_from:
                rejected.append((line.rstrip("\n"), "ali_to<ali_from"))
                continue
            if (ali_to - ali_from + 1) < min_length:
                continue

            ensps: list[str]
            if id_type == "ensp":
                ensps = [strip_version(protein)] if looks_like_ensp(protein) else []
            elif id_type == "uniprot":
                ensps = mapping.lookup(protein) if mapping else []
            else:
                ensps = ([strip_version(protein)] if looks_like_ensp(protein)
                         else (mapping.lookup(protein) if mapping and looks_like_uniprot(protein) else []))
            if not ensps:
                rejected.append((line.rstrip("\n"), f"no ENSP for {protein}"))
                continue

            pfam_acc = strip_version(hmm_acc) if hmm_acc and hmm_acc != "-" else ""
            tag = pfam_acc or hmm_name or "Pfam_hit"
            domain_id = f"{hmm_name}_{tag}" if (hmm_name and pfam_acc) else f"{tag}_Pfam"
            desc_safe = re.sub(r"\s+", " ", description).strip()

            for ensp in ensps:
                rows.append((ensp, ali_from, ali_to, domain_id, desc_safe))

    return rows, rejected
