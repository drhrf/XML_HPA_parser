"""HPA XML Parser: per-sample RNA expression from Human Protein Atlas per-gene XML files.

Basic use::

    import hpa_xml_parser as hxp

    gene = hxp.fetch_gene("ACADM")          # or hxp.parse_xml("ENSG00000117054.xml.gz")
    df = hxp.to_frame(gene, assays=["humanBrain"])
    hxp.export(gene, "ACADM.xlsx")
"""

from ._version import __version__
from .export import export, metadata_frame, to_frame
from .fetch import (
    GeneResolutionError,
    download_xml,
    fetch_gene,
    fetch_genes,
    read_gene_list,
    resolve_gene,
)
from .parser import ParsedGene, parse_xml
from .schema import AGGREGATE_COLUMNS, SAMPLE_ASSAYS, SAMPLE_COLUMNS
from .validation import HPAStructureWarning, HPAVersionMismatchWarning

__all__ = [
    "AGGREGATE_COLUMNS",
    "GeneResolutionError",
    "HPAStructureWarning",
    "HPAVersionMismatchWarning",
    "ParsedGene",
    "SAMPLE_ASSAYS",
    "SAMPLE_COLUMNS",
    "__version__",
    "download_xml",
    "export",
    "fetch_gene",
    "fetch_genes",
    "metadata_frame",
    "parse_xml",
    "read_gene_list",
    "resolve_gene",
    "to_frame",
]
