# Changelog

All notable changes to this project are documented in this file. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
- Zenodo DOI of release 2.0.0 ([10.5281/zenodo.24091185](https://doi.org/10.5281/zenodo.24091185)) in `CITATION.cff` and the
  README, with a DOI badge and a ready-to-use citation.

## [2.0.0] - 2026-11-09

Complete rewrite, validated against Human Protein Atlas release 25 (XML schema 3.0,
Ensembl 109).

### Added
- Installable Python package `hpa_xml_parser` (`pyproject.toml`, `src/` layout, Python ≥ 3.10)
  with the API `parse_xml`, `fetch_gene`, `fetch_genes`, `to_frame` and `export`.
- `hpa-xml-parser` command-line interface with `fetch` (Ensembl IDs, gene symbols or
  `--genes-file`) and `parse` (local `.xml`/`.xml.gz`) subcommands.
- Download by Ensembl ID or gene symbol (resolved through the HPA search API), with local
  cache, timeout, retries and a generic User-Agent. `--base-url` pins an archived release.
- Support for all assays with per-sample data: `tissue`, `humanBrain`, `pigBrain`,
  `mouseBrain` and `immuneCell`.
- Tidy long output with a fixed column schema, plus a table of HPA's aggregated values
  (nTPM, pTPM, TPM) for every assay, including `consensusTissue` and the `*Regional` assays.
- Provenance in every output: HPA release, XML schema, Ensembl, assembly and GENCODE
  versions, download time.
- Export to CSV, TSV, XLSX (one sheet per assay, plus `aggregates` and `metadata`) and
  Parquet.
- Warnings for unknown assay types, attributes and level types, and for genes from
  different HPA releases in one batch.
- pytest suite with HPA v25 fixtures (ACADM, GFAP, ACTB), mocked downloads, and a
  concordance test against HPA's aggregates; testthat tests for the R bridge; GitHub
  Actions CI (Python 3.10 and 3.12, R, ruff).
- `docs/hpa_xml_structure_v25.md` describing the XML structure the parser relies on.
- Example output `examples/ACADM_v25.*`.
- `CITATION.cff` and `.zenodo.json`.

### Changed
- Attributes are read by name. Version 1 assigned column names by position, so any
  added or missing attribute (such as the absent `age` on some `tissue` samples) could
  shift values into the wrong columns.
- The Shiny app now calls the Python package through reticulate instead of re-implementing
  the parser in R. It accepts uploads or gene IDs/symbols, shows the HPA release, exports
  the same files as the CLI, and adds boxplots by region, tissue, sex and age group.
- The notebook installs the package and uses the API; it no longer asks for a file path.
- README rewritten in English, with the output schema, validation and limitations.

### Removed
- `hpa_xml_parser.py` (Colab export with a Portuguese `input()` prompt and duplicated
  functions), replaced by the package.
- Theme selector and the shinydashboard, shinythemes and shinycssloaders dependencies.
- `README.html` and `README_files/`, a stale rendering of the old README.

### Fixed
- License statement: the project is licensed under GPL-3.0 (the README previously said MIT).
- README references to a non-existent `HPA_XML_Parser.py` and `example_data/` folder.

## [1.0.0]

Initial version: a Colab-exported Python script and an independent R/Shiny parser for the
`humanBrain` and `tissue` assays.

[2.0.0]: https://github.com/drhrf/XML_HPA_parser/releases/tag/v2.0.0
[1.0.0]: https://github.com/drhrf/XML_HPA_parser/releases/tag/v1.0.0
