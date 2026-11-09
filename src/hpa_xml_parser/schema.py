"""Known structure of the HPA per-gene XML and the output column schema.

Everything listed here was observed in HPA release 25 (schema version 3.0);
see ``docs/hpa_xml_structure_v25.md``. Anything outside these sets triggers a
:class:`~hpa_xml_parser.validation.HPAStructureWarning` but is still parsed.
"""

from __future__ import annotations

#: ``rnaExpression/@assayType`` values that carry per-sample ``RNASample`` data.
SAMPLE_ASSAYS: tuple[str, ...] = (
    "tissue",
    "humanBrain",
    "pigBrain",
    "mouseBrain",
    "immuneCell",
)

#: ``rnaExpression/@assayType`` values that carry only aggregated ``level`` data
#: (or only summary elements).
AGGREGATE_ONLY_ASSAYS: tuple[str, ...] = (
    "consensusTissue",
    "humanBrainRegional",
    "pigBrainRegional",
    "mouseBrainRegional",
    "cellLine",
    "immuneCellLineage",
)

KNOWN_ASSAYS: frozenset[str] = frozenset(SAMPLE_ASSAYS + AGGREGATE_ONLY_ASSAYS)

#: Species of each assay (the XML does not state it; derived from the assay name).
ASSAY_SPECIES: dict[str, str] = {
    "tissue": "human",
    "humanBrain": "human",
    "immuneCell": "human",
    "consensusTissue": "human",
    "humanBrainRegional": "human",
    "cellLine": "human",
    "immuneCellLineage": "human",
    "pigBrain": "pig",
    "pigBrainRegional": "pig",
    "mouseBrain": "mouse",
    "mouseBrainRegional": "mouse",
}

#: Attributes observed on ``RNASample``.
KNOWN_SAMPLE_ATTRIBUTES: frozenset[str] = frozenset({"sampleId", "unitRNA", "expRNA", "sex", "age"})

#: Attributes of the label element (``tissue`` / ``immuneCell``) that hold the
#: grouping value written to the ``organ`` output column, in priority order.
GROUPING_ATTRIBUTES: tuple[str, ...] = ("organ", "region", "lineage")

#: All attributes observed on label elements.
KNOWN_LABEL_ATTRIBUTES: frozenset[str] = frozenset(GROUPING_ATTRIBUTES + ("ontologyTerms",))

#: Tags observed for the label element inside ``data``.
KNOWN_LABEL_TAGS: frozenset[str] = frozenset({"tissue", "immuneCell"})

#: Attributes observed on ``level``.
KNOWN_LEVEL_ATTRIBUTES: frozenset[str] = frozenset({"type", "unitRNA", "expRNA"})

#: ``level/@type`` values observed.
KNOWN_LEVEL_TYPES: frozenset[str] = frozenset(
    {"normalizedRNAExpression", "proteinCodingRNAExpression", "RNAExpression"}
)

#: Core columns of the per-sample (tidy, long) table, in output order.
SAMPLE_COLUMNS: tuple[str, ...] = (
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
    "value",
    "sex",
    "age",
    "hpa_version",
    "retrieved_at",
)

#: Columns appended after :data:`SAMPLE_COLUMNS` in every per-sample table.
SAMPLE_EXTRA_COLUMNS: tuple[str, ...] = ("species", "label_attribute")

#: Core columns of the aggregated (per tissue / region) table.
AGGREGATE_COLUMNS: tuple[str, ...] = (
    "gene_name",
    "ensembl_id",
    "assay_type",
    "source",
    "technology",
    "organ",
    "region",
    "ontology_terms",
    "level_type",
    "unit",
    "value",
    "n_samples",
    "hpa_version",
    "retrieved_at",
    "species",
    "label_attribute",
)

#: Prefix for unknown ``RNASample`` attributes preserved as extra columns.
EXTRA_SAMPLE_PREFIX = "extra_"
#: Prefix for unknown label-element attributes preserved as extra columns.
EXTRA_LABEL_PREFIX = "extra_label_"
