"""Command-line interface: ``hpa-xml-parser``."""

from __future__ import annotations

import argparse
import logging
import sys
import warnings
from collections.abc import Sequence

from ._version import __version__
from .export import FORMATS, export
from .fetch import DEFAULT_BASE_URL, DEFAULT_RETRIES, DEFAULT_TIMEOUT, fetch_genes, read_gene_list
from .parser import ParsedGene, parse_xml
from .schema import KNOWN_ASSAYS, SAMPLE_ASSAYS
from .validation import check_version_consistency

log = logging.getLogger("hpa_xml_parser")


def _common_output_args(p: argparse.ArgumentParser) -> None:
    p.add_argument(
        "--assay",
        nargs="+",
        metavar="ASSAY",
        help=(
            "assay types to export (default: all with per-sample data: "
            + ", ".join(SAMPLE_ASSAYS)
            + ")"
        ),
    )
    p.add_argument("-o", "--out", required=True, help="output file (.csv, .tsv, .xlsx, .parquet)")
    p.add_argument("--format", choices=FORMATS, help="output format (default: from extension)")
    p.add_argument(
        "--no-aggregates", action="store_true", help="do not export aggregated level values"
    )


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser."""
    parser = argparse.ArgumentParser(
        prog="hpa-xml-parser",
        description=(
            "Extract per-sample RNA expression (nTPM, sex, age) from Human Protein Atlas "
            "per-gene XML files."
        ),
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    f = sub.add_parser("fetch", help="download genes from the HPA by Ensembl ID or symbol")
    f.add_argument("genes", nargs="*", help="Ensembl gene IDs or gene symbols")
    f.add_argument("--genes-file", help="text file with gene IDs/symbols (one per line)")
    _common_output_args(f)
    f.add_argument(
        "--base-url", default=DEFAULT_BASE_URL, help="HPA base URL (default: %(default)s)"
    )
    f.add_argument("--cache-dir", help="cache directory (default: ~/.cache/hpa-xml-parser)")
    f.add_argument("--no-cache", action="store_true", help="do not read or write the cache")
    f.add_argument("--refresh", action="store_true", help="re-download even if cached")
    f.add_argument("--timeout", type=float, default=DEFAULT_TIMEOUT, help="seconds per request")
    f.add_argument("--retries", type=int, default=DEFAULT_RETRIES, help="retries per request")
    f.add_argument(
        "--skip-errors", action="store_true", help="skip genes that fail instead of stopping"
    )

    p = sub.add_parser("parse", help="parse local .xml or .xml.gz files")
    p.add_argument("files", nargs="+", help="HPA per-gene XML files")
    _common_output_args(p)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Entry point. Returns the process exit code."""
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    warnings.simplefilter("always")
    logging.captureWarnings(True)

    if args.assay:
        unknown = sorted(set(args.assay) - KNOWN_ASSAYS)
        if unknown:
            log.warning("warning: assay type(s) not seen in HPA v25: %s", ", ".join(unknown))

    genes: list[ParsedGene]
    if args.command == "fetch":
        ids = list(args.genes)
        if args.genes_file:
            ids.extend(read_gene_list(args.genes_file))
        if not ids:
            log.error("error: give gene IDs/symbols or --genes-file")
            return 2
        genes = fetch_genes(
            ids,
            on_error="warn" if args.skip_errors else "raise",
            base_url=args.base_url,
            timeout=args.timeout,
            retries=args.retries,
            cache_dir=args.cache_dir,
            use_cache=not args.no_cache,
            refresh=args.refresh,
        )
    else:
        genes = [parse_xml(path) for path in args.files]
        check_version_consistency({g.gene_name or "?": g.hpa_version for g in genes})

    if not genes:
        log.error("error: no gene could be processed")
        return 1
    for g in genes:
        log.info(
            "%s (%s): HPA v%s, %d samples", g.gene_name, g.ensembl_id, g.hpa_version, len(g.samples)
        )
    for path in export(
        genes,
        args.out,
        fmt=args.format,
        assays=args.assay,
        include_aggregates=not args.no_aggregates,
    ):
        log.info("wrote %s", path)
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
