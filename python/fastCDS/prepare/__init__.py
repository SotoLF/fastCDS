"""Read domain annotations from common sources into fastCDS-ready DataFrames.

Each ``from_*`` function returns a ``pandas.DataFrame`` with the canonical
columns ``protein_id, aa_start, aa_end, domain_id, description``. The
DataFrame is consumed directly by :func:`fastCDS.Mapper.map_batch` —
there's no need to write an intermediate BED file unless you want one.

Examples
--------

>>> import fastCDS as fc
>>> queries = fc.prepare.from_pfam("hits.dom", mode="scan", id_type="ensp")
>>> queries.head()
   protein_id  aa_start  aa_end                domain_id    description
0  ENSP000…         5      120  PF00069_protein_kinase   Protein kin…

>>> mapper = fc.Mapper(index="human.idx")
>>> result = mapper.map_batch(queries)
"""

from __future__ import annotations

from typing import Iterable

import pandas as pd

from ._mapping import UniProtToEnsp, dedup_and_sort_rows
from . import _interpro, _pfam, _uniprot


__all__ = [
    "from_pfam",
    "from_interproscan",
    "from_uniprot_features",
    "rows_to_dataframe",
    "load_uniprot_mapping",
]

_COLUMNS = ("protein_id", "aa_start", "aa_end", "domain_id", "description")


def rows_to_dataframe(rows: Iterable[tuple]) -> pd.DataFrame:
    """Convert ``(ensp, aa_start, aa_end, domain_id, description)`` rows
    to a DataFrame with the fastCDS-shaped column names."""
    return pd.DataFrame(list(rows), columns=list(_COLUMNS))


def load_uniprot_mapping(path: str, *, format: str = "ensembl_xref") -> UniProtToEnsp:
    """Load a UniProt → ENSP mapping for files that use UniProt accessions.

    ``format='ensembl_xref'`` (default) reads the Ensembl per-release UniProt
    xref TSV; ``format='simple'`` reads a two-column ``uniprot\\tensp`` TSV.
    """
    if format == "ensembl_xref":
        return UniProtToEnsp.from_ensembl_xref_tsv(path)
    if format == "simple":
        return UniProtToEnsp.from_simple_tsv(path)
    raise ValueError(f"unknown mapping format: {format!r} "
                     "(expected 'ensembl_xref' or 'simple')")


def from_pfam(path: str,
              *,
              mode: str = "scan",
              min_score: float = 0.0,
              min_length: int = 5,
              id_type: str = "auto",
              mapping: UniProtToEnsp | None = None,
              dedup: bool = True,
              return_rejected: bool = False,
              ) -> pd.DataFrame | tuple[pd.DataFrame, list]:
    """Read an HMMER ``--domtblout`` file into a DataFrame.

    Parameters
    ----------
    path : str
        Path to the HMMER ``--domtblout`` file.
    mode : {'scan', 'search'}
        ``'scan'`` for ``hmmscan`` output (target=HMM, query=protein — typical).
        ``'search'`` for ``hmmsearch`` output (target=protein, query=HMM).
    min_score : float
        Drop hits with this-domain bit score below this value.
    min_length : int
        Drop hits shorter than this many amino acids.
    id_type : {'auto', 'ensp', 'uniprot'}
        Protein-ID style. ``'auto'`` sniffs the first 50 rows.
    mapping : UniProtToEnsp, optional
        Required when ``id_type='uniprot'``. Build via :func:`load_uniprot_mapping`.
    dedup : bool
        Drop exact-duplicate rows and sort. Default True.
    return_rejected : bool
        If True, return ``(DataFrame, rejected_list)`` instead of just the
        DataFrame. Useful for QC.
    """
    rows, rejected = _pfam.parse(
        path, mode=mode, min_score=min_score, min_length=min_length,
        id_type=id_type, mapping=mapping,
    )
    if dedup:
        rows = dedup_and_sort_rows(rows)
    df = rows_to_dataframe(rows)
    if return_rejected:
        return df, rejected
    return df


def from_interproscan(path: str,
                      *,
                      analyses: set[str] | None = None,
                      min_length: int = 5,
                      id_type: str = "auto",
                      mapping: UniProtToEnsp | None = None,
                      source_filter: set[str] | None = None,
                      dedup: bool = True,
                      return_rejected: bool = False,
                      ) -> pd.DataFrame | tuple[pd.DataFrame, list]:
    """Read an InterProScan TSV (the ``-f TSV`` output) into a DataFrame.

    Parameters
    ----------
    path : str
        Path to the InterProScan ``.tsv`` (gzip OK).
    analyses : set[str], optional
        Restrict to specific analyses (e.g. ``{'Pfam', 'SMART'}``). Default:
        include all.
    min_length : int
        Drop hits shorter than this many amino acids.
    id_type : {'auto', 'ensp', 'uniprot'}
        Same sniffing semantics as :func:`from_pfam`.
    mapping : UniProtToEnsp, optional
        Required when ``id_type='uniprot'``.
    dedup : bool
        Drop exact-duplicate rows and sort. Default True.
    return_rejected : bool
        If True, return ``(DataFrame, rejected_list)``.
    """
    rows, rejected = _interpro.parse(
        path, analyses=analyses, min_length=min_length,
        id_type=id_type, mapping=mapping, source_filter=source_filter,
    )
    if dedup:
        rows = dedup_and_sort_rows(rows)
    df = rows_to_dataframe(rows)
    if return_rejected:
        return df, rejected
    return df


def from_uniprot_features(path: str,
                          *,
                          format: str = "auto",
                          feature_types: set[str] | None = None,
                          min_length: int = 5,
                          fallback_mapping: UniProtToEnsp | None = None,
                          dedup: bool = True,
                          return_rejected: bool = False,
                          ) -> pd.DataFrame | tuple[pd.DataFrame, list]:
    """Read a UniProt feature table (``.dat`` flat-file or REST ``.json``).

    UniProt entries usually carry their own Ensembl cross-references in
    ``DR Ensembl;`` lines — those are used first, with ``fallback_mapping``
    consulted only for entries that lack an Ensembl xref.

    Parameters
    ----------
    path : str
        Path to a UniProt ``.dat`` / ``.dat.gz`` / ``.txt`` / ``.json`` file.
    format : {'auto', 'dat', 'json'}
        ``'auto'`` decides from the extension.
    feature_types : set[str], optional
        Which ``FT`` types to keep. Default: ``{DOMAIN, REPEAT, REGION,
        ZN_FING, DNA_BIND, COILED, TRANSMEM, MOTIF, TOPO_DOM}``.
    min_length : int
        Drop features shorter than this many amino acids.
    fallback_mapping : UniProtToEnsp, optional
        Used only when an entry has no inline Ensembl xref.
    dedup : bool
        Drop exact-duplicate rows and sort. Default True.
    return_rejected : bool
        If True, return ``(DataFrame, rejected_list)``.
    """
    if format == "auto":
        low = path.lower().rstrip(".gz")
        format = "json" if low.endswith(".json") else "dat"
    feat_types = feature_types if feature_types is not None else _uniprot.DEFAULT_FEATURE_TYPES

    if format == "dat":
        raw, dr_by_acc, rejected = _uniprot.parse_dat(
            path, feature_types=feat_types)
    elif format == "json":
        raw, dr_by_acc, rejected = _uniprot.parse_json(
            path, feature_types=feat_types)
    else:
        raise ValueError(f"unknown format: {format!r}")

    rows: list[tuple[str, int, int, str, str]] = []
    for acc, s, e, did, desc in raw:
        if (e - s + 1) < min_length:
            continue
        ensps = dr_by_acc.get(acc, [])
        if not ensps and fallback_mapping is not None:
            ensps = fallback_mapping.lookup(acc)
        if not ensps:
            rejected.append((f"{acc} {s}-{e} {did}", "no ENSP"))
            continue
        for ensp in ensps:
            rows.append((ensp, s, e, did, desc))

    if dedup:
        rows = dedup_and_sort_rows(rows)
    df = rows_to_dataframe(rows)
    if return_rejected:
        return df, rejected
    return df
