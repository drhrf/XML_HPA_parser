"""Combine parsed genes into tables and write them to CSV, TSV, XLSX or Parquet."""

from __future__ import annotations

import os
import re
from collections.abc import Iterable, Sequence
from pathlib import Path

import pandas as pd

from ._version import __version__
from .parser import ParsedGene, _utc_iso
from .schema import SAMPLE_ASSAYS
from .validation import check_version_consistency

GeneInput = ParsedGene | Iterable[ParsedGene]

FORMATS = ("csv", "tsv", "xlsx", "parquet")


def _as_list(genes: GeneInput) -> list[ParsedGene]:
    if isinstance(genes, ParsedGene):
        return [genes]
    return list(genes)


def to_frame(
    genes: GeneInput,
    *,
    assays: Sequence[str] | None = None,
    kind: str = "samples",
) -> pd.DataFrame:
    """Concatenate the tables of one or more genes into one long DataFrame.

    Args:
        genes: A :class:`ParsedGene` or an iterable of them.
        assays: Keep only these ``assay_type`` values. ``None`` keeps all.
        kind: ``"samples"`` for per-sample values or ``"aggregates"`` for the
            aggregated ``level`` values.

    Returns:
        A DataFrame with the core columns first and any extra columns after.
    """
    if kind not in ("samples", "aggregates"):
        raise ValueError("kind must be 'samples' or 'aggregates'")
    frames = [getattr(g, kind) for g in _as_list(genes)]
    frames = [f for f in frames if len(f.columns)]
    if not frames:
        return pd.DataFrame()
    df = pd.concat(frames, ignore_index=True, sort=False)
    if assays is not None:
        df = df[df["assay_type"].isin(list(assays))].reset_index(drop=True)
    return df


def metadata_frame(genes: GeneInput) -> pd.DataFrame:
    """One row per gene with provenance (HPA/Ensembl versions, retrieval time, source)."""
    rows = [g.metadata() for g in _as_list(genes)]
    df = pd.DataFrame(rows)
    df.insert(0, "parser_version", __version__)
    df.insert(1, "exported_at", _utc_iso())
    return df


def _infer_format(path: Path, fmt: str | None) -> str:
    if fmt:
        fmt = fmt.lower()
    else:
        suffixes = [s.lower() for s in path.suffixes]
        fmt = next((s.lstrip(".") for s in reversed(suffixes) if s.lstrip(".") in FORMATS), "")
        if not fmt and suffixes and suffixes[-1] == ".xls":
            fmt = "xlsx"
    if fmt not in FORMATS:
        raise ValueError(f"Unsupported output format {fmt!r}; choose one of {FORMATS}")
    return fmt


def _sheet_name(name: str, used: set[str]) -> str:
    base = re.sub(r"[\[\]:*?/\\]", "_", name)[:31] or "sheet"
    candidate, i = base, 1
    while candidate.lower() in used:
        i += 1
        candidate = f"{base[: 31 - len(str(i)) - 1]}_{i}"
    used.add(candidate.lower())
    return candidate


def export(
    genes: GeneInput,
    out: str | os.PathLike[str],
    *,
    fmt: str | None = None,
    assays: Sequence[str] | None = None,
    include_aggregates: bool = True,
) -> list[Path]:
    """Write per-sample data, aggregates and metadata for one or more genes.

    * ``xlsx``: one sheet per assay with per-sample data, an ``aggregates``
      sheet and a ``metadata`` sheet (HPA version, dates, genes).
    * ``csv`` / ``tsv`` / ``parquet``: ``out`` holds the per-sample long table;
      ``<stem>_aggregates.<ext>`` and ``<stem>_metadata.<ext>`` are written
      next to it.

    Args:
        genes: Parsed gene(s).
        out: Output path. The format is taken from the extension unless ``fmt`` is given.
        fmt: One of ``csv``, ``tsv``, ``xlsx``, ``parquet``.
        assays: Assays to export. Defaults to all assays with per-sample data.
        include_aggregates: Also write the aggregated values.

    Returns:
        The paths written.
    """
    gene_list = _as_list(genes)
    check_version_consistency(
        {g.gene_name or g.ensembl_id or "?": g.hpa_version for g in gene_list}
    )
    path = Path(out)
    fmt = _infer_format(path, fmt)
    sample_assays = list(assays) if assays is not None else list(SAMPLE_ASSAYS)
    samples = to_frame(gene_list, assays=sample_assays, kind="samples")
    aggregates = to_frame(gene_list, assays=assays, kind="aggregates")
    meta = metadata_frame(gene_list)
    path.parent.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []

    if fmt == "xlsx":
        used: set[str] = set()
        with pd.ExcelWriter(path, engine="openpyxl") as writer:
            for assay in sample_assays:
                part = samples[samples["assay_type"] == assay] if len(samples) else samples
                if len(part):
                    part.to_excel(writer, sheet_name=_sheet_name(assay, used), index=False)
            if include_aggregates and len(aggregates):
                aggregates.to_excel(writer, sheet_name=_sheet_name("aggregates", used), index=False)
            meta.to_excel(writer, sheet_name=_sheet_name("metadata", used), index=False)
        return [path]

    stem = path.name[: -len("".join(path.suffixes))] if path.suffixes else path.name
    ext = "".join(path.suffixes) or f".{fmt}"

    def _write(df: pd.DataFrame, target: Path) -> None:
        if fmt == "csv":
            df.to_csv(target, index=False)
        elif fmt == "tsv":
            df.to_csv(target, index=False, sep="\t")
        else:
            try:
                df.to_parquet(target, index=False)
            except ImportError as exc:  # pragma: no cover - depends on environment
                raise ImportError(
                    "Parquet export needs pyarrow: pip install 'hpa-xml-parser[parquet]'"
                ) from exc
        written.append(target)

    _write(samples, path)
    if include_aggregates:
        _write(aggregates, path.with_name(f"{stem}_aggregates{ext}"))
    _write(meta, path.with_name(f"{stem}_metadata{ext}"))
    return written
