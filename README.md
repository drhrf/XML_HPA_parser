# HPA XML Parser

[![CI](https://github.com/drhrf/XML_HPA_parser/actions/workflows/ci.yml/badge.svg)](https://github.com/drhrf/XML_HPA_parser/actions/workflows/ci.yml)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/drhrf/XML_HPA_parser/blob/main/HPA_XML_Parser.ipynb)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)

HPA XML Parser extracts **per-sample, gene-level RNA expression** (nTPM, with donor sex and
age) from the [Human Protein Atlas](https://www.proteinatlas.org/) (HPA) per-gene XML files.
It covers human tissues, human brain subregions, pig and mouse brain regions, and immune
cells. It also extracts HPA's own aggregated values per tissue or region, and records the
HPA and Ensembl releases in every output.

It is available as a Python package with a command-line interface, a Shiny app, and a
Google Colab notebook.

## Why this tool

* In the current HPA release (v25, Ensembl 109), the **bulk downloads** give **per-sample**
  values only at the **transcript** level. **Gene-level** values are given only
  **aggregated by tissue or region**.
* The per-gene XML (`https://www.proteinatlas.org/<ENSG>.xml`) is the HPA source that
  contains gene-level nTPM **for each sample**, together with sex and age.
* [HPAanalyze](https://bioconductor.org/packages/HPAanalyze/) (Bioconductor) reads the
  XML for antibody and immunohistochemistry data only, not for these RNA samples.

HPA XML Parser fills this gap. It reads every attribute **by name**, preserves attributes
it does not know instead of dropping or shifting them, warns when the XML structure
changes, and outputs one tidy table per run.

## Installation

Python ≥ 3.10:

```bash
pip install "hpa-xml-parser @ git+https://github.com/drhrf/XML_HPA_parser"
```

Optional Parquet export:

```bash
pip install "hpa-xml-parser[parquet] @ git+https://github.com/drhrf/XML_HPA_parser"
```

## Command line

```bash
# Download by gene symbol or Ensembl ID; one Excel sheet per assay
hpa-xml-parser fetch ACADM GFAP --assay humanBrain tissue --out results.xlsx

# Gene list from a text file (one per line; '#' starts a comment), CSV output
hpa-xml-parser fetch --genes-file genes.txt --out results.csv

# Parse XML files you already downloaded (.xml or .xml.gz)
hpa-xml-parser parse ENSG00000117054.xml.gz --out ACADM.tsv
```

Useful options: `--assay` (default: all assays with per-sample data), `--no-aggregates`,
`--format`, `--cache-dir`, `--no-cache`, `--refresh`, `--timeout`, `--retries`,
`--skip-errors`, and `--base-url`. For example, `--base-url https://v25.proteinatlas.org`
pins an archived release. Run `hpa-xml-parser fetch --help` for details.

Downloads use a timeout, retries with backoff, and the User-Agent
`hpa-xml-parser/<version>`. They are cached in `~/.cache/hpa-xml-parser` (or
`$HPA_XML_PARSER_CACHE`) together with the original download time.

## Python API

```python
import hpa_xml_parser as hxp

gene = hxp.fetch_gene("ACADM")              # or "ENSG00000117054"
gene = hxp.parse_xml("ENSG00000117054.xml.gz")  # local file

gene.hpa_version, gene.ensembl_version      # ('25', '109')
brain = hxp.to_frame(gene, assays=["humanBrain"])           # per-sample long table
regions = hxp.to_frame(gene, kind="aggregates")             # HPA aggregated values

genes = hxp.fetch_genes(["ACADM", "GFAP", "ACTB"])          # batch
hxp.export(genes, "results.xlsx")                           # or .csv / .tsv / .parquet
```

`parse_xml` and `fetch_gene` return a `ParsedGene` with `samples`, `aggregates`, `assays`
(counts per assay), version fields, and `issues` (structural warnings).

## Shiny app

The app is a graphical front end to the same Python package (through
[reticulate](https://rstudio.github.io/reticulate/)), so it produces the same columns and
files as the command line. It accepts uploaded XML files or gene IDs/symbols, shows the HPA
release, exports XLSX or CSV/TSV, and draws boxplots by region, tissue, sex and age group.

```r
install.packages(c("shiny", "bslib", "reticulate", "DT", "ggplot2"))
reticulate::py_install("hpa-xml-parser @ git+https://github.com/drhrf/XML_HPA_parser", pip = TRUE)
shiny::runApp("path/to/XML_HPA_parser")
```

To use a specific Python environment, set `HPA_XML_PARSER_PYTHON` (for example
`Sys.setenv(HPA_XML_PARSER_PYTHON = "~/venvs/hpa/bin/python")`) before starting the app.

## Notebook

[`HPA_XML_Parser.ipynb`](HPA_XML_Parser.ipynb) installs the package and walks through
download, the concordance check, a plot and export. Open it in Colab with the badge above.

## Output schema

### Per-sample table (`samples`; one sheet per assay in XLSX)

| Column | Content | XML source |
|---|---|---|
| `gene_name` | HPA gene symbol | `entry/name` |
| `ensembl_id` | Ensembl gene ID | `entry/identifier/@id` |
| `assay_type` | `tissue`, `humanBrain`, `pigBrain`, `mouseBrain`, `immuneCell` | `rnaExpression/@assayType` |
| `source`, `technology` | `HPA`, `RNAseq` | `rnaExpression/@source`, `@technology` |
| `organ` | grouping of the label: organ system (tissue), brain-region code such as `CTX` or `AMY` (brain assays), lineage (immune cells) | `tissue/@organ`, `tissue/@region` or `immuneCell/@lineage` |
| `region` | tissue, brain subregion or immune cell type | text of `tissue` / `immuneCell` |
| `ontology_terms` | UBERON term(s), when given | `tissue/@ontologyTerms` |
| `sample_id` | HPA sample ID (text) | `RNASample/@sampleId` |
| `unit` | always `nTPM` in v25 | `RNASample/@unitRNA` |
| `value` | expression value (float) | `RNASample/@expRNA` |
| `sex` | `Female` / `Male` | `RNASample/@sex` |
| `age` | age as reported (float; empty if absent) | `RNASample/@age` |
| `hpa_version` | HPA release, e.g. `25` | `entry/@version` |
| `retrieved_at` | UTC time the XML was downloaded (for local files: file modification time) | — |
| `species` | `human`, `pig` or `mouse` (derived from the assay name) | — |
| `label_attribute` | which attribute filled `organ`, e.g. `tissue@organ` | — |
| `extra_*`, `extra_label_*` | any attribute not seen in v25, kept verbatim | `RNASample/@*`, label `@*` |

### Aggregated table (`aggregates` sheet or `<name>_aggregates.<ext>`)

This table has the same identification columns plus `level_type` (`normalizedRNAExpression`
= nTPM, `proteinCodingRNAExpression` = pTPM, `RNAExpression` = TPM), `unit`, `value` and
`n_samples`. It covers every assay, including the aggregate-only ones (`consensusTissue`,
`humanBrainRegional`, `pigBrainRegional`, `mouseBrainRegional`).

### Metadata (`metadata` sheet or `<name>_metadata.<ext>`)

One row per gene with `parser_version`, `exported_at`, `hpa_version`, `schema_version`,
`ensembl_version`, `assembly`, `gencode_version`, `entry_url`, `retrieved_at`, `origin`,
`n_samples` and `issues`.

The structure behind this schema is documented in
[`docs/hpa_xml_structure_v25.md`](docs/hpa_xml_structure_v25.md).

## Validation

* **Concordance with HPA.** For every tissue, brain subregion and immune cell type in the
  test genes (ACADM, GFAP, ACTB), the mean of the extracted per-sample nTPM equals HPA's
  aggregated `normalizedRNAExpression` within 0.1 nTPM. That is the bound for two
  one-decimal roundings; the largest difference observed is 0.075 nTPM.
  See [`tests/test_concordance.py`](tests/test_concordance.py).
* **No column shift.** Every extracted row is checked against the raw `RNASample`
  attributes. A synthetic XML with an unknown attribute placed out of order is parsed
  correctly, and the unknown attribute is kept.
* CI runs the tests on Python 3.10 and 3.12, the R bridge tests, and ruff.

## Limitations

* **Post mortem tissue.** Human brain samples come from post mortem donors. The XML does
  not record sample origin, post mortem interval or cause of death.
* **Age range.** In v25, human brain donors are 18–94 years old, and 907 of the 966 brain
  samples come from donors aged 61 or older; younger ages are poorly represented. Human
  tissue ages range from 0.1 to 85, with some missing. Pig (`1`) and mouse (`0.2`) ages
  are constant, and immune cell samples have no age. The XML does not state age units.
* **No donor ID.** `sampleId` identifies a tissue sample, not a donor. Samples from the
  same donor cannot be linked, so they should not be treated as independent observations
  in statistical tests.
* **No hemisphere.** The XML has no hemisphere information (no attribute, region name or
  ID encoding).
* **nTPM only per sample.** pTPM and TPM are available only as aggregates.
* **Regional summaries are not means.** In v25, the 13-region human brain values
  (`humanBrainRegional`) equal the maximum subregion value, not the mean of the samples.
* **Pig and mouse data** are the values HPA reports for the orthologs of the human gene,
  as they appear in the human gene's XML. Their region codes differ from the human ones.
* **Dependence on HPA's XML format.** The parser is validated against HPA v25 (schema
  3.0). Unknown assay types, attributes or level types trigger a warning, and the data
  are kept, but future HPA releases may change the XML. The minor release (e.g. 25.1) is
  not in the XML; report it, with the download date, separately.

## Citation

If you use this tool, please cite it (see [`CITATION.cff`](CITATION.cff); a Zenodo DOI will
be added on release). Please also cite the Human Protein Atlas, as requested on its
website.

## License

GNU General Public License v3.0. See [`LICENSE`](LICENSE).
