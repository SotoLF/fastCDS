"""Query Ensembl REST API's /map/translation endpoint for a batch of
(protein_id, aa_start, aa_end, query_id) rows and emit a TSV in the same shape
as `ensembldb_query.R` so the comparison logic in validate_vs_ensembldb.py
can reuse the same classifier.

The public REST server is rate-limited (15 req/s, ~55K/hour). The script throttles
to stay just under that with a small safety margin.

Usage:
    python run_ensembl_rest.py <queries.bed> <out.tsv> [--limit N] [--qps 12]
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path


REST_URL = "https://rest.ensembl.org/map/translation/{pid}/{start}..{end}?content-type=application/json"


def query_one(pid: str, start: int, end: int, retries: int = 3) -> tuple[list, str]:
    """Return ([(chrom, start, end)], status). Status: 'ok', 'no_result', 'error', 'rate_limited'."""
    url = REST_URL.format(pid=pid, start=start, end=end)
    backoff = 1.0
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=30) as r:
                data = json.loads(r.read())
            mappings = data.get("mappings", [])
            if not mappings:
                return [], "no_result"
            out = []
            for m in mappings:
                if "gap" in m and m.get("gap", 0):
                    continue
                out.append((m["seq_region_name"], int(m["start"]), int(m["end"])))
            return out, "ok" if out else "no_result"
        except urllib.error.HTTPError as e:
            if e.code == 429:
                # Honor Retry-After if present, else backoff.
                ra = e.headers.get("Retry-After")
                wait = float(ra) if ra else backoff
                time.sleep(wait)
                backoff *= 2
                continue
            if e.code in (400, 404):
                # Bad request / not found - protein not in their annotation.
                return [], "no_result"
            return [], f"error_http{e.code}"
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as e:
            time.sleep(backoff)
            backoff *= 2
            continue
    return [], "error_retries_exhausted"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("bed", type=Path, help="BED-like queries (4 cols: protein_id, aa_start, aa_end, query_id)")
    ap.add_argument("out", type=Path, help="Output TSV (same schema as ensembldb_query.R)")
    ap.add_argument("--limit", type=int, default=1000,
                    help="Cap queries (REST is rate-limited; PLAN says 1000)")
    ap.add_argument("--qps", type=float, default=12.0,
                    help="Requests per second (REST limit is 15; we stay under)")
    args = ap.parse_args()

    sleep_per = 1.0 / args.qps
    rows = []
    with open(args.bed) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            rows.append((parts[0], int(parts[1]), int(parts[2]), parts[3]))
            if len(rows) >= args.limit:
                break

    print(f"querying {len(rows):,} ENSPs at {args.qps:.1f} req/s "
          f"(estimated {len(rows)/args.qps:.0f}s)", file=sys.stderr)
    t0 = time.perf_counter()
    with open(args.out, "w") as f:
        f.write("query_id\tchrom\tstart\tend\tstrand\tstatus\n")
        for i, (pid, s, e, qid) in enumerate(rows, 1):
            t_before = time.perf_counter()
            intervals, status = query_one(pid, s, e)
            if not intervals:
                f.write(f"{qid}\tNA\tNA\tNA\tNA\t{status}\n")
            else:
                for c, gs, ge in intervals:
                    f.write(f"{qid}\t{c}\t{gs}\t{ge}\t.\tok\n")
            # Throttle to qps.
            elapsed = time.perf_counter() - t_before
            if elapsed < sleep_per:
                time.sleep(sleep_per - elapsed)
            if i % 100 == 0:
                rate = i / (time.perf_counter() - t0)
                print(f"  {i:,}/{len(rows):,}  ({rate:.1f} q/s effective)", file=sys.stderr)

    total = time.perf_counter() - t0
    print(f"done in {total:.1f}s ({len(rows)/total:.2f} q/s overall)", file=sys.stderr)


if __name__ == "__main__":
    main()
