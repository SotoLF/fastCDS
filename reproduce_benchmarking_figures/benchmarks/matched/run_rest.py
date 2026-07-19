"""Ensembl REST route: one /map/translation call per query.

Serial and network-bound. The public endpoint allows 15 req/s, but the ~900 ms
HTTP round-trip is the real ceiling for a script written the way anyone
actually writes one, which is what this measures. That caps the ladder at
N = 1,000.

For the agreement runs (which need the mapped coordinates rather than the
timing) use ../run_ensembl_rest.py, which writes the classifier TSV and can
run a paced thread pool.

Usage: python run_rest.py <N> <ids.txt>
"""

from __future__ import annotations

import json
import sys
import urllib.request
from pathlib import Path

URL = ("https://rest.ensembl.org/map/translation/{pid}/1..50"
       "?content-type=application/json")


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("run_rest.py: need <N> <ids.txt>")
    n = int(sys.argv[1])
    ids_path = Path(sys.argv[2])
    if not ids_path.exists():
        raise SystemExit(f"run_rest.py: ids file not found: {ids_path}")

    ids = [l.strip() for l in open(ids_path) if l.strip()][:n]
    ok = 0
    for pid in ids:
        try:
            req = urllib.request.Request(
                URL.format(pid=pid), headers={"User-Agent": "fastCDS-bench"})
            with urllib.request.urlopen(req, timeout=30) as r:
                json.loads(r.read())
            ok += 1
        except Exception:
            pass
    print(f"rest mapped {ok} / {len(ids)}")


if __name__ == "__main__":
    main()
