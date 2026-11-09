"""Parse a Human Protein Atlas per-gene XML file into tidy tables.

All attributes are read **by name**. Attributes that were not observed in HPA
v25 are preserved in extra columns (``extra_<name>`` for ``RNASample``,
``extra_label_<name>`` for the ``tissue``/``immuneCell`` label element) and
reported with :class:`~hpa_xml_parser.validation.HPAStructureWarning`.
"""

from __future__ import annotations

import datetime as _dt
import gzip
import io
import os
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path
from typing import IO, Any

import pandas as pd

from . import schema
from .validation import warn_structure

XMLSource = str | os.PathLike[str] | bytes | IO[bytes]

_GZIP_MAGIC = b"\x1f\x8b"
_STRING_SAMPLE_COLUMNS = (
    "gene_name",
    "ensembl_id",
    "assay_type",
    "source",
    "technology",
    "organ",
    "region",
    "ontology_terms",
    "sample_id",
    "unit",
    "sex",
    "hpa_version",
    "retrieved_at",
    "species",
    "label_attribute",
)


@dataclass
class ParsedGene:
    """Contents of one HPA per-gene XML file.

    Attributes:
        gene_name: HPA gene symbol (``entry/name``).
        ensembl_id: Ensembl gene ID (``entry/identifier/@id``).
        hpa_version: HPA release (``entry/@version``), e.g. ``"25"``.
        schema_version: XML schema version (``proteinAtlas/@schemaVersion``).
        ensembl_version: Ensembl release (``identifier/@version``).
        assembly: Genome assembly (``identifier/@assembly``).
        gencode_version: GENCODE release (``identifier/@gencodeVersion``).
        entry_url: Versioned HPA URL of the entry (``entry/@url``).
        retrieved_at: ISO 8601 UTC timestamp of when the XML was downloaded
            (for local files: the file's modification time unless given).
        origin: Path or URL the XML was read from.
        samples: Per-sample long table (see :data:`schema.SAMPLE_COLUMNS`).
        aggregates: Aggregated ``level`` values per tissue/region
            (see :data:`schema.AGGREGATE_COLUMNS`).
        assays: One row per ``rnaExpression`` element with its counts.
        issues: Structural warnings raised while parsing.
    """

    gene_name: str | None
    ensembl_id: str | None
    hpa_version: str | None
    schema_version: str | None
    ensembl_version: str | None
    assembly: str | None
    gencode_version: str | None
    entry_url: str | None
    retrieved_at: str | None
    origin: str | None
    samples: pd.DataFrame
    aggregates: pd.DataFrame
    assays: pd.DataFrame
    issues: list[str] = field(default_factory=list)

    def metadata(self) -> dict[str, Any]:
        """Return the provenance fields of this gene as a flat dictionary."""
        return {
            "gene_name": self.gene_name,
            "ensembl_id": self.ensembl_id,
            "hpa_version": self.hpa_version,
            "schema_version": self.schema_version,
            "ensembl_version": self.ensembl_version,
            "assembly": self.assembly,
            "gencode_version": self.gencode_version,
            "entry_url": self.entry_url,
            "retrieved_at": self.retrieved_at,
            "origin": self.origin,
            "n_samples": int(len(self.samples)),
            "issues": "; ".join(self.issues),
        }

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return (
            f"ParsedGene({self.gene_name!r}, {self.ensembl_id!r}, hpa_version={self.hpa_version!r},"
            f" samples={len(self.samples)}, aggregates={len(self.aggregates)})"
        )


def _utc_iso(ts: float | None = None) -> str:
    if ts is None:
        moment = _dt.datetime.now(_dt.timezone.utc)
    else:
        moment = _dt.datetime.fromtimestamp(ts, tz=_dt.timezone.utc)
    return moment.replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _read_bytes(source: XMLSource) -> tuple[bytes, str | None, float | None]:
    """Return raw (decompressed) bytes, an origin label and a file mtime."""
    origin: str | None = None
    mtime: float | None = None
    if isinstance(source, (bytes, bytearray)):
        raw = bytes(source)
    elif hasattr(source, "read"):
        raw = source.read()  # type: ignore[union-attr]
        origin = getattr(source, "name", None)
    else:
        path = Path(os.fspath(source))
        raw = path.read_bytes()
        origin = str(path)
        mtime = path.stat().st_mtime
    if raw[:2] == _GZIP_MAGIC:
        raw = gzip.decompress(raw)
    return raw, origin, mtime


def _num(value: str | None) -> float:
    if value is None or value == "":
        return float("nan")
    try:
        return float(value)
    except ValueError:
        return float("nan")


def _label_fields(
    label: ET.Element | None, assay: str, issues: list[str]
) -> tuple[dict[str, Any], dict[str, str]]:
    """Normalize a ``tissue``/``immuneCell`` element.

    Returns:
        A tuple of (core fields, extra columns).
    """
    core: dict[str, Any] = {
        "organ": None,
        "region": None,
        "ontology_terms": None,
        "label_attribute": None,
    }
    extras: dict[str, str] = {}
    if label is None:
        warn_structure(f"assay {assay!r}: <data> block without a label element", issues)
        return core, extras
    if label.tag not in schema.KNOWN_LABEL_TAGS:
        warn_structure(f"assay {assay!r}: unknown label element <{label.tag}>", issues)
    core["region"] = (label.text or "").strip() or None
    for attr in schema.GROUPING_ATTRIBUTES:
        if attr in label.attrib:
            core["organ"] = label.attrib[attr]
            core["label_attribute"] = f"{label.tag}@{attr}"
            break
    core["ontology_terms"] = label.attrib.get("ontologyTerms")
    for key, val in label.attrib.items():
        if key not in schema.KNOWN_LABEL_ATTRIBUTES:
            extras[schema.EXTRA_LABEL_PREFIX + key] = val
    return core, extras


def parse_xml(source: XMLSource, *, retrieved_at: str | None = None) -> ParsedGene:
    """Parse an HPA per-gene XML document.

    Args:
        source: Path to a ``.xml`` or ``.xml.gz`` file, raw bytes, or a binary
            file object. Gzip compression is detected from the content, not the
            file name.
        retrieved_at: ISO 8601 timestamp of when the XML was downloaded. If
            omitted, the file's modification time (UTC) is used for paths, and
            the current time otherwise.

    Returns:
        A :class:`ParsedGene` with per-sample and aggregated tables for every
        ``rnaExpression`` assay in the file.
    """
    raw, origin, mtime = _read_bytes(source)
    root = ET.parse(io.BytesIO(raw)).getroot()
    issues: list[str] = []

    entry = root.find("entry") if root.tag != "entry" else root
    if entry is None:
        raise ValueError("Not an HPA per-gene XML: no <entry> element found")
    if len(root.findall("entry")) > 1:
        warn_structure("XML contains more than one <entry>; only the first is parsed", issues)

    ident = entry.find("identifier")
    ident_attrib = ident.attrib if ident is not None else {}
    gene_name = (entry.findtext("name") or "").strip() or None
    ensembl_id = ident_attrib.get("id")
    hpa_version = entry.attrib.get("version")
    if hpa_version is None:
        warn_structure("entry/@version (HPA release) is missing", issues)
    if retrieved_at is None:
        retrieved_at = _utc_iso(mtime)

    common = {
        "gene_name": gene_name,
        "ensembl_id": ensembl_id,
        "hpa_version": hpa_version,
        "retrieved_at": retrieved_at,
    }

    sample_rows: list[dict[str, Any]] = []
    aggregate_rows: list[dict[str, Any]] = []
    assay_rows: list[dict[str, Any]] = []

    for rx in entry.findall("rnaExpression"):
        assay = rx.attrib.get("assayType")
        if assay is None:
            warn_structure("rnaExpression without @assayType", issues)
            assay = "unknown"
        elif assay not in schema.KNOWN_ASSAYS:
            warn_structure(f"unknown assayType {assay!r}", issues)
        assay_common = {
            **common,
            "assay_type": assay,
            "source": rx.attrib.get("source"),
            "technology": rx.attrib.get("technology"),
            "species": schema.ASSAY_SPECIES.get(assay),
        }
        n_data = n_samp = 0
        for data in rx.findall("data"):
            n_data += 1
            label = next((c for c in data if c.tag not in ("level", "RNASample")), None)
            label_core, label_extra = _label_fields(label, assay, issues)
            samples = data.findall("RNASample")
            n_samp += len(samples)
            for level in data.findall("level"):
                ltype = level.attrib.get("type")
                if ltype not in schema.KNOWN_LEVEL_TYPES:
                    warn_structure(f"assay {assay!r}: unknown level type {ltype!r}", issues)
                row = {
                    **assay_common,
                    **label_core,
                    "level_type": ltype,
                    "unit": level.attrib.get("unitRNA"),
                    "value": _num(level.attrib.get("expRNA")),
                    "n_samples": len(samples),
                    **label_extra,
                }
                for key, val in level.attrib.items():
                    if key not in schema.KNOWN_LEVEL_ATTRIBUTES:
                        row[schema.EXTRA_SAMPLE_PREFIX + key] = val
                aggregate_rows.append(row)
            for rs in samples:
                a = rs.attrib
                unknown = [k for k in a if k not in schema.KNOWN_SAMPLE_ATTRIBUTES]
                row = {
                    **assay_common,
                    **label_core,
                    "sample_id": a.get("sampleId"),
                    "unit": a.get("unitRNA"),
                    "value": _num(a.get("expRNA")),
                    "sex": a.get("sex"),
                    "age": _num(a.get("age")),
                    **label_extra,
                }
                for key in unknown:
                    row[schema.EXTRA_SAMPLE_PREFIX + key] = a[key]
                if unknown:
                    msg = (
                        f"assay {assay!r}: unknown RNASample attribute(s) "
                        f"{sorted(unknown)} kept as extra columns"
                    )
                    if msg not in issues:
                        warn_structure(msg, issues)
                if "expRNA" in a and row["value"] != row["value"]:  # NaN check
                    warn_structure(
                        f"assay {assay!r}: non-numeric expRNA {a['expRNA']!r} "
                        f"for sample {a.get('sampleId')!r}",
                        issues,
                    )
                sample_rows.append(row)
        if n_samp and assay in schema.AGGREGATE_ONLY_ASSAYS:
            warn_structure(f"assay {assay!r} unexpectedly contains RNASample data", issues)
        assay_rows.append(
            {
                "assay_type": assay,
                "source": assay_common["source"],
                "technology": assay_common["technology"],
                "species": assay_common["species"],
                "n_data_blocks": n_data,
                "n_samples": n_samp,
            }
        )

    return ParsedGene(
        gene_name=gene_name,
        ensembl_id=ensembl_id,
        hpa_version=hpa_version,
        schema_version=root.attrib.get("schemaVersion"),
        ensembl_version=ident_attrib.get("version"),
        assembly=ident_attrib.get("assembly"),
        gencode_version=ident_attrib.get("gencodeVersion"),
        entry_url=entry.attrib.get("url"),
        retrieved_at=retrieved_at,
        origin=origin,
        samples=_finalize(sample_rows, schema.SAMPLE_COLUMNS + schema.SAMPLE_EXTRA_COLUMNS),
        aggregates=_finalize(aggregate_rows, schema.AGGREGATE_COLUMNS),
        assays=pd.DataFrame(
            assay_rows,
            columns=[
                "assay_type",
                "source",
                "technology",
                "species",
                "n_data_blocks",
                "n_samples",
            ],
        ),
        issues=issues,
    )


def _finalize(rows: list[dict[str, Any]], core: tuple[str, ...]) -> pd.DataFrame:
    """Build a DataFrame with core columns first, extras after, and fixed dtypes."""
    extras: list[str] = []
    for row in rows:
        for key in row:
            if key not in core and key not in extras:
                extras.append(key)
    df = pd.DataFrame(rows, columns=list(core) + sorted(extras))
    for col in df.columns:
        if col in ("value", "age"):
            df[col] = df[col].astype("float64")
        elif col == "n_samples":
            df[col] = df[col].astype("int64")
        elif col in _STRING_SAMPLE_COLUMNS or col == "level_type" or col.startswith("extra_"):
            df[col] = df[col].astype("string")
    return df
