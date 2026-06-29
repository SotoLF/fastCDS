"""MappingResult dataclass + parsers for the files the C++ binary writes.

Reading an output directory back into Python is its own useful thing — for
example, when you already ran the binary on the command line and want to
analyze the TSVs interactively. `read_results_dir(path)` does that.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from typing import Any

import pandas as pd


# Column lists for the headerless BED outputs. We attach names when reading so
# downstream code can use `df.chrom` etc. instead of positional indexing.
BED6_COLS = ["chrom", "start", "end", "name", "score", "strand"]
BED12_COLS = ["chrom", "start", "end", "name", "score", "strand",
              "thickStart", "thickEnd", "itemRgb",
              "blockCount", "blockSizes", "blockStarts"]


def _read_tsv(path: str) -> pd.DataFrame:
    """Read one of the TSV outputs. Empty / missing files yield an empty
    DataFrame (with no columns)."""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return pd.DataFrame()
    return pd.read_csv(
        path, sep="\t",
        na_values=["NA"], keep_default_na=True,
        dtype=str,            # parse as string by default ...
        low_memory=False,
    ).pipe(_coerce_numerics)


# Columns that should be numeric where present. Keeping the conversion
# centralized means new C++ columns just need to land in one of these sets.
_INT_COLS = {
    "aa_start", "aa_end", "domain_length_aa", "domain_length_nt",
    "protein_length_aa", "domain_genomic_start", "domain_genomic_end",
    "n_coding_segments", "cds_nt_remainder",
    "feature_part", "exon_number",
    "feature_genomic_start", "feature_genomic_end", "feature_length_nt",
    "feature_order_genomic", "feature_order_transcript",
    "cds_nt_start", "cds_nt_end", "aa_start_encoded", "aa_end_encoded",
    "domain_overlap_genomic_start", "domain_overlap_genomic_end",
    "domain_overlap_cds_nt_start", "domain_overlap_cds_nt_end",
    "domain_overlap_aa_start", "domain_overlap_aa_end",
}
_FLOAT_COLS = {
    "domain_overlap_fraction_of_feature",
    "domain_overlap_fraction_of_domain",
}
_BOOL_COLS = {
    "fully_mapped", "no_domain_mode", "cds_length_mismatch",
}
_TRIBOOL_COLS = {
    "is_mane_select", "is_ensembl_canonical",
}


def _coerce_numerics(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    for c in df.columns:
        if c in _INT_COLS:
            df[c] = pd.to_numeric(df[c], errors="coerce").astype("Int64")
        elif c in _FLOAT_COLS:
            df[c] = pd.to_numeric(df[c], errors="coerce")
        elif c in _BOOL_COLS:
            df[c] = df[c].map({"true": True, "false": False}).astype("boolean")
        elif c in _TRIBOOL_COLS:
            # Keep as string ("true"/"false"/NaN) so downstream code can
            # distinguish "false" from "unknown" cleanly. pandas boolean dtype
            # collapses NaN onto something less obvious.
            df[c] = df[c].astype("string")
    return df


def _read_bed(path: str, *, ncols: int) -> pd.DataFrame:
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return pd.DataFrame()
    cols = BED12_COLS if ncols == 12 else BED6_COLS
    df = pd.read_csv(path, sep="\t", header=None, names=cols, dtype=str)
    for c in ("start", "end", "score", "thickStart", "thickEnd", "blockCount"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce").astype("Int64")
    return df


@dataclass
class MappingResult:
    """Everything the C++ binary produced for one invocation.

    Tables are pandas DataFrames, possibly empty. `metadata` is the parsed
    `run_metadata.json` (only present for `--output all`).
    """

    summary: pd.DataFrame
    isoform: pd.DataFrame
    cds_segments: pd.DataFrame
    introns: pd.DataFrame
    cds_bed: pd.DataFrame
    introns_bed: pd.DataFrame
    span_bed: pd.DataFrame
    bed12: pd.DataFrame
    unmapped: pd.DataFrame
    metadata: dict[str, Any] | None = None
    out_dir: str | None = None

    # ------------------------------------------------------------------ #
    # Stats
    # ------------------------------------------------------------------ #
    @property
    def n_total(self) -> int:
        return len(self.summary)

    @property
    def n_unmapped(self) -> int:
        return len(self.unmapped)

    @property
    def n_mapped(self) -> int:
        return self.n_total - self.n_unmapped

    def __repr__(self) -> str:
        return (f"MappingResult(n_total={self.n_total}, "
                f"n_mapped={self.n_mapped}, n_unmapped={self.n_unmapped}, "
                f"out_dir={self.out_dir!r})")

    # ------------------------------------------------------------------ #
    # Filters
    # ------------------------------------------------------------------ #
    def by_input_id(self, input_id: str) -> "MappingResult":
        """Return a new MappingResult containing only rows for one query."""
        keep = lambda df: (df[df["input_id"] == input_id]
                           if not df.empty and "input_id" in df.columns
                           else df)
        bed_keep = lambda df: (df[df["name"].str.contains(input_id, na=False)]
                               if not df.empty else df)
        return MappingResult(
            summary=keep(self.summary),
            isoform=keep(self.isoform),
            cds_segments=keep(self.cds_segments),
            introns=keep(self.introns),
            cds_bed=bed_keep(self.cds_bed),
            introns_bed=bed_keep(self.introns_bed),
            span_bed=bed_keep(self.span_bed),
            bed12=bed_keep(self.bed12),
            unmapped=keep(self.unmapped),
            metadata=self.metadata,
            out_dir=self.out_dir,
        )

    def write(self, out_dir: str) -> None:
        """Persist all tables back to a directory in the same layout the C++
        binary uses. Useful when you ran the wrapper with `keep_outputs=False`
        but later decide you want the files on disk."""
        os.makedirs(out_dir, exist_ok=True)

        def w_tsv(df: pd.DataFrame, name: str) -> None:
            if df.empty: return
            df.to_csv(os.path.join(out_dir, name), sep="\t",
                      index=False, na_rep="NA")

        def w_bed(df: pd.DataFrame, name: str) -> None:
            if df.empty: return
            df.to_csv(os.path.join(out_dir, name), sep="\t",
                      index=False, header=False, na_rep="NA")

        w_tsv(self.summary,      "domain_mapping_summary.tsv")
        w_tsv(self.isoform,      "isoform_structure.tsv")
        w_tsv(self.cds_segments, "domain_cds_segments.tsv")
        w_tsv(self.introns,      "domain_introns.tsv")
        w_tsv(self.unmapped,     "unmapped_domains.tsv")
        w_bed(self.cds_bed,      "domain_cds_segments.bed")
        w_bed(self.introns_bed,  "domain_introns.bed")
        w_bed(self.span_bed,     "domain_span_with_introns.bed")
        w_bed(self.bed12,        "domain_blocks.bed12")
        if self.metadata is not None:
            with open(os.path.join(out_dir, "run_metadata.json"), "w") as f:
                json.dump(self.metadata, f, indent=2)


def read_results_dir(out_dir: str) -> MappingResult:
    """Parse an existing fastCDS output directory into a MappingResult."""
    if not os.path.isdir(out_dir):
        raise FileNotFoundError(f"Not a directory: {out_dir}")
    j = lambda name: os.path.join(out_dir, name)
    meta = None
    meta_path = j("run_metadata.json")
    if os.path.exists(meta_path):
        with open(meta_path) as f:
            meta = json.load(f)
    return MappingResult(
        summary=_read_tsv(j("domain_mapping_summary.tsv")),
        isoform=_read_tsv(j("isoform_structure.tsv")),
        cds_segments=_read_tsv(j("domain_cds_segments.tsv")),
        introns=_read_tsv(j("domain_introns.tsv")),
        cds_bed=_read_bed(j("domain_cds_segments.bed"), ncols=6),
        introns_bed=_read_bed(j("domain_introns.bed"), ncols=6),
        span_bed=_read_bed(j("domain_span_with_introns.bed"), ncols=6),
        bed12=_read_bed(j("domain_blocks.bed12"), ncols=12),
        unmapped=_read_tsv(j("unmapped_domains.tsv")),
        metadata=meta,
        out_dir=out_dir,
    )
