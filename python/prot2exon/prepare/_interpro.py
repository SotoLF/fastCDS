"""Parse InterProScan TSV into prot2exon BED-like rows.

Extracted from `scripts/prepare_from_interpro.py`.
"""

from __future__ import annotations

import re
import sys

from ._mapping import (UniProtToEnsp, looks_like_ensp, looks_like_uniprot,
                       open_text)


def detect_id_type(samples: list[str]) -> str:
    """Best-effort: vote on the dominant ID style in the file."""
    n_up = sum(1 for s in samples if looks_like_uniprot(s))
    n_ensp = sum(1 for s in samples if looks_like_ensp(s))
    if n_ensp >= max(n_up, 1):
        return "ensp"
    if n_up >= 1:
        return "uniprot"
    return "unknown"


def parse(path: str, *, analyses: set[str] | None, min_length: int,
          id_type: str, mapping: UniProtToEnsp | None,
          source_filter: set[str] | None) -> tuple[list, list]:
    """Return ``(rows, rejected)``."""
    rows: list[tuple[str, int, int, str, str]] = []
    rejected: list[tuple[str, str]] = []

    if id_type == "auto":
        samples: list[str] = []
        with open_text(path) as fh:
            for i, line in enumerate(fh):
                if line.startswith("#") or not line.strip(): continue
                samples.append(line.split("\t", 1)[0].strip())
                if i > 200: break
        id_type = detect_id_type(samples)
        print(f"[interpro] auto-detected id_type = {id_type}", file=sys.stderr)

    if id_type == "uniprot" and mapping is None:
        raise ValueError(
            "InterPro file uses UniProt accessions but no mapping given."
        )

    with open_text(path) as fh:
        for ln, line in enumerate(fh, start=1):
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                rejected.append((line.rstrip("\n"), "too few columns"))
                continue
            try:
                src_id      = parts[0].strip()
                analysis    = parts[3].strip()
                sig_acc     = parts[4].strip()
                sig_desc    = parts[5].strip()
                aa_start    = int(parts[6])
                aa_end      = int(parts[7])
            except (ValueError, IndexError):
                rejected.append((line.rstrip("\n"), "unparseable row"))
                continue
            if analyses and analysis not in analyses:
                continue
            if aa_end < aa_start:
                rejected.append((line.rstrip("\n"), "aa_end<aa_start"))
                continue
            if (aa_end - aa_start + 1) < min_length:
                continue
            ipr_acc = parts[11].strip() if len(parts) > 11 else "-"

            chosen_acc = ipr_acc if ipr_acc and ipr_acc != "-" else sig_acc
            if not chosen_acc:
                chosen_acc = "hit"
            domain_id = f"{chosen_acc}_{analysis}"
            sig_desc_safe = re.sub(r"\s+", " ", sig_desc).strip()

            ensps: list[str]
            if id_type == "ensp":
                ensps = [src_id.split(".", 1)[0]] if looks_like_ensp(src_id) else []
            elif id_type == "uniprot":
                ensps = mapping.lookup(src_id) if mapping else []
            else:
                ensps = [src_id.split(".", 1)[0]] if looks_like_ensp(src_id) else (
                    mapping.lookup(src_id) if mapping and looks_like_uniprot(src_id) else [])

            if not ensps:
                rejected.append((line.rstrip("\n"), f"no ENSP for {src_id}"))
                continue

            for ensp in ensps:
                rows.append((ensp, aa_start, aa_end, domain_id, sig_desc_safe))

    return rows, rejected
