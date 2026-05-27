"""Parse UniProt feature tables (.dat / REST .json) into prot2exon rows.

Extracted from `scripts/prepare_from_uniprot_features.py`. Supports both
the classic flat-file with ``FT`` records and the REST JSON output.
"""

from __future__ import annotations

import json
import re

from ._mapping import strip_version, open_text


DEFAULT_FEATURE_TYPES = {
    "DOMAIN", "REPEAT", "REGION", "ZN_FING", "DNA_BIND",
    "COILED", "TRANSMEM", "MOTIF", "TOPO_DOM",
}


# Match: "FT   DOMAIN          93..312" — fuzzy positions are skipped.
FT_HEADER_RE = re.compile(
    r"^FT\s+(\w+)\s+(?P<start><?\??\d+)\.\.(?P<end>>?\??\d+)"
)
FT_NOTE_RE = re.compile(r'/note="([^"]+)"')
FT_ID_RE   = re.compile(r'/id="([^"]+)"')
DR_ENSEMBL_RE = re.compile(r"^DR\s+Ensembl;\s+([^;\s]+);\s+([^;\s]+);\s+([^.\s]+)")
AC_RE   = re.compile(r"^AC\s+([^;]+);")


def parse_dat(path: str, *, feature_types: set[str]
              ) -> tuple[list, dict[str, list[str]], list]:
    """Parse a UniProt .dat / .txt / .dat.gz flat-file.

    Returns ``(rows, dr_ensp_by_acc, rejected)``.
    """
    rows: list[tuple[str, int, int, str, str]] = []
    dr_by_acc: dict[str, list[str]] = {}
    rejected: list[tuple[str, str]] = []

    cur_acc: str | None = None
    cur_dr: list[str] = []
    pending: dict | None = None

    def flush_pending():
        nonlocal pending
        if not pending or cur_acc is None:
            pending = None
            return
        rows.append((cur_acc, pending["s"], pending["e"], pending["domain"],
                     pending["desc"]))
        pending = None

    with open_text(path) as fh:
        for line in fh:
            if line.startswith("ID "):
                flush_pending()
                if cur_acc is not None:
                    dr_by_acc[cur_acc] = cur_dr
                cur_acc = None
                cur_dr = []
                continue
            if line.startswith("AC "):
                m = AC_RE.match(line)
                if m and cur_acc is None:
                    cur_acc = m.group(1).strip()
                continue
            if line.startswith("DR "):
                m = DR_ENSEMBL_RE.match(line)
                if m:
                    enst, ensp, ensg = m.group(1), m.group(2), m.group(3)
                    cur_dr.append(strip_version(ensp))
                continue
            if line.startswith("//"):
                flush_pending()
                if cur_acc is not None:
                    dr_by_acc[cur_acc] = list(cur_dr)
                cur_acc, cur_dr = None, []
                continue
            if not line.startswith("FT"):
                continue
            m = FT_HEADER_RE.match(line)
            if m:
                flush_pending()
                ftype = re.match(r"^FT\s+(\w+)", line).group(1)
                if ftype not in feature_types:
                    pending = None
                    continue
                s_str = m.group("start").lstrip("<>?")
                e_str = m.group("end").lstrip("<>?")
                if not s_str.isdigit() or not e_str.isdigit():
                    rejected.append((line.rstrip("\n"),
                                     "fuzzy location, skipped"))
                    pending = None
                    continue
                s, e = int(s_str), int(e_str)
                if e < s:
                    rejected.append((line.rstrip("\n"), "end<start"))
                    pending = None
                    continue
                pending = {"s": s, "e": e, "type": ftype,
                           "domain": ftype, "desc": ""}
                continue
            if pending is not None:
                m_note = FT_NOTE_RE.search(line)
                if m_note:
                    pending["desc"] = re.sub(r"\s+", " ", m_note.group(1)).strip()
                    if pending["desc"] and len(pending["desc"]) <= 40:
                        pending["domain"] = (
                            pending["type"] + "_" +
                            re.sub(r"[^A-Za-z0-9_.-]+", "_", pending["desc"]))
                m_id = FT_ID_RE.search(line)
                if m_id:
                    pending["domain"] = pending["type"] + "_" + m_id.group(1)

        flush_pending()
        if cur_acc is not None:
            dr_by_acc[cur_acc] = list(cur_dr)

    return rows, dr_by_acc, rejected


def parse_json(path: str, *, feature_types: set[str]
               ) -> tuple[list, dict[str, list[str]], list]:
    """Parse a UniProt REST .json file (single entry, list, or `{results: [...]}`)."""
    rows: list[tuple[str, int, int, str, str]] = []
    dr_by_acc: dict[str, list[str]] = {}
    rejected: list[tuple[str, str]] = []

    def consume_entry(entry: dict):
        acc = entry.get("primaryAccession", "")
        if not acc:
            return
        ensps: list[str] = []
        for xref in (entry.get("uniProtKBCrossReferences", []) or []):
            if xref.get("database") == "Ensembl":
                for prop in (xref.get("properties", []) or []):
                    if prop.get("key") in ("ProteinId", "EnsemblProteinId"):
                        v = prop.get("value", "")
                        if v: ensps.append(strip_version(v))
        if ensps:
            dr_by_acc[acc] = ensps
        for feat in (entry.get("features", []) or []):
            ftype = (feat.get("type") or "").upper().replace(" ", "_")
            if ftype not in feature_types and ftype.replace("-", "_") not in feature_types:
                continue
            loc = feat.get("location", {}) or {}
            s = (loc.get("start", {}) or {}).get("value")
            e = (loc.get("end", {}) or {}).get("value")
            if s is None or e is None or not isinstance(s, int) or not isinstance(e, int):
                rejected.append((json.dumps(feat), "fuzzy/missing location"))
                continue
            if e < s:
                rejected.append((json.dumps(feat), "end<start"))
                continue
            desc = (feat.get("description") or "").strip()
            domain_id = ftype
            fid = feat.get("featureId")
            if fid:
                domain_id = f"{ftype}_{fid}"
            elif desc:
                domain_id = (f"{ftype}_" +
                             re.sub(r"[^A-Za-z0-9_.-]+", "_", desc)[:40])
            rows.append((acc, s, e, domain_id,
                         re.sub(r"\s+", " ", desc).strip()))

    with open_text(path) as fh:
        first = fh.readline()
        fh.seek(0)
        if first.lstrip().startswith("{") and not first.lstrip().startswith('{"results"'):
            ok = True
            try:
                json.loads(first)
            except json.JSONDecodeError:
                ok = False
            if ok:
                for line in fh:
                    line = line.strip()
                    if not line: continue
                    try:
                        consume_entry(json.loads(line))
                    except json.JSONDecodeError as e:
                        rejected.append((line[:120], f"bad json: {e}"))
                return rows, dr_by_acc, rejected
        doc = json.load(fh)
        if isinstance(doc, dict) and "results" in doc:
            for entry in doc["results"]:
                consume_entry(entry)
        elif isinstance(doc, list):
            for entry in doc:
                consume_entry(entry)
        elif isinstance(doc, dict):
            consume_entry(doc)
        else:
            rejected.append((str(type(doc)), "unrecognized JSON shape"))
    return rows, dr_by_acc, rejected
