# Example output: ACADM, HPA v25

Generated with hpa-xml-parser 2.0.0 from the XML of ACADM (ENSG00000117054) downloaded
from https://www.proteinatlas.org on 2026-11-09 (HPA release 25, Ensembl 109). The XML is
byte-identical to `tests/fixtures/ENSG00000117054_ACADM_v25.xml.gz` (SHA-256 of the
uncompressed file: `4839e056…dfebb03`).

```bash
hpa-xml-parser fetch ENSG00000117054 --out examples/ACADM_v25.xlsx
hpa-xml-parser fetch ENSG00000117054 --out examples/ACADM_v25.csv
```

| File | Content |
|---|---|
| `ACADM_v25.xlsx` | one sheet per assay (`tissue`, `humanBrain`, `pigBrain`, `mouseBrain`, `immuneCell`), plus `aggregates` and `metadata` |
| `ACADM_v25.csv` | all 1,447 samples in one long table |
| `ACADM_v25_aggregates.csv` | HPA aggregated nTPM/pTPM/TPM per tissue, region and subregion |
| `ACADM_v25_metadata.csv` | HPA/Ensembl versions, download and export times |

To reproduce from the fixture without network access:

```bash
hpa-xml-parser parse tests/fixtures/ENSG00000117054_ACADM_v25.xml.gz --out ACADM_v25.csv
```

The values are identical. Only `retrieved_at` differs: for local files it is the file's
modification time.
