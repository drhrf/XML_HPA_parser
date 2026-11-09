# Structure of the Human Protein Atlas per-gene XML (v25)

This document records what the per-gene XML served at
`https://www.proteinatlas.org/<ENSG>.xml` actually contains in HPA release 25.
Everything the parser does is based on these observations; nothing here is
inferred from the XSD or from documentation alone.

## Method

Files downloaded on 2026-11-09 with `curl -A "hpa-xml-parser/2.0.0"` and stored,
gzip-compressed, in `tests/fixtures/`:

| Gene  | Ensembl ID        | Fixture                                   | SHA-256 of uncompressed XML |
|-------|-------------------|-------------------------------------------|-----------------------------|
| ACADM | ENSG00000117054   | `ENSG00000117054_ACADM_v25.xml.gz`        | `4839e056076b4904498128cd2a61191ef1bf9c2162a1abda6aed32998dfebb03` |
| GFAP  | ENSG00000131095   | `ENSG00000131095_GFAP_v25.xml.gz`         | `d59f501e8f15e5b3d68c140188204253c4e8e4aa83646b24481dfc1ca9751149` |
| ACTB  | ENSG00000075624   | `ENSG00000075624_ACTB_v25.xml.gz`         | `b7ac4f255769220e832490a1389ce02d2bbd595051fa2486e95538d830b73bbe` |

Every `rnaExpression` element was enumerated with `xml.etree.ElementTree`, and
the attribute sets of all `data/*` children were tallied. The three genes have
an **identical** structure: same assay types, same attribute sets, and the same
number of `data` blocks and samples per assay.

## Version information

```xml
<proteinAtlas xsi:schemaLocation="https://v25.proteinatlas.org/download/proteinatlas.xsd" schemaVersion="3.0">
  <entry version="25" url="https://v25.proteinatlas.org/ENSG00000117054">
    <name>ACADM</name>
    <identifier id="ENSG00000117054" db="Ensembl" version="109" assembly="GRCh38.p13" gencodeVersion="43">
```

| Information        | Location                                   | Value (v25) |
|--------------------|--------------------------------------------|-------------|
| HPA release        | `entry/@version`                           | `25` |
| HPA release (URL)  | `entry/@url`, `proteinAtlas/@xsi:schemaLocation` (subdomain `v25.`) | `v25` |
| XML schema version | `proteinAtlas/@schemaVersion`              | `3.0` |
| Ensembl release    | `entry/identifier/@version`                | `109` |
| Genome assembly    | `entry/identifier/@assembly`               | `GRCh38.p13` |
| GENCODE release    | `entry/identifier/@gencodeVersion`         | `43` |

The **minor** release number (the "25.1" shown on the website) does **not**
appear anywhere in the XML, and the HTTP response headers carry no version
either. The parser therefore records `hpa_version = "25"` plus the schema,
Ensembl, assembly and GENCODE versions, and the retrieval timestamp. If the
minor release matters, it must be stated by the user (for example in the
manuscript's Methods) together with the download date.

## `rnaExpression` assay types

All `rnaExpression` elements have `source="HPA"` and `technology="RNAseq"`.
Counts are identical for ACADM, GFAP and ACTB.

| `assayType`          | `data` blocks | `RNASample` | Label element & attributes          | Per-sample data |
|----------------------|--------------:|------------:|-------------------------------------|-----------------|
| `consensusTissue`    | 51            | 0           | `tissue` (`organ`, `ontologyTerms`*) | no (aggregated only) |
| `tissue`             | 40            | 186         | `tissue` (`organ`, `ontologyTerms`)  | **yes** |
| `humanBrain`         | 193           | 966         | `tissue` (`organ`, `ontologyTerms`*) | **yes** |
| `humanBrainRegional` | 13            | 0           | `tissue` (`region`, `ontologyTerms`*) | no (aggregated only) |
| `mouseBrain`         | 17            | 67          | `tissue` (`region`, `ontologyTerms`) | **yes** |
| `mouseBrainRegional` | 11            | 0           | `tissue` (`region`, `ontologyTerms`*) | no (aggregated only) |
| `pigBrain`           | 30            | 119         | `tissue` (`region`, `ontologyTerms`*) | **yes** |
| `pigBrainRegional`   | 13            | 0           | `tissue` (`region`, `ontologyTerms`*) | no (aggregated only) |
| `cellLine`           | 0             | 0           | none              | no (only summary elements) |
| `immuneCell`         | 19            | 109         | `immuneCell` (`lineage`)            | **yes** |
| `immuneCellLineage`  | 0             | 0           | none              | no (only summary elements) |

\* `ontologyTerms` is missing on some elements (for example, 183 of the 193
`humanBrain` subregions have no UBERON term).

Besides `data`, some `rnaExpression` elements contain summary children:
`rnaSpecificity` (`@specificity`, `@description`), `rnaDistribution`
(`@description`, text) and `rnaExpressionCluster` (`@clusterID`, text).
`cellTypeExpression` (scRNA-seq, `singleCellTypeExpression` with nCPM) is a
separate element. It holds aggregated values per cell type only, with no
samples.

### Layout of a `data` block

```xml
<rnaExpression source="HPA" technology="RNAseq" assayType="humanBrain">
  <data>
    <tissue organ="AMY" ontologyTerms="UBERON:0001876">Amygdala</tissue>
    <level type="normalizedRNAExpression"   unitRNA="nTPM" expRNA="2.9" />
    <level type="proteinCodingRNAExpression" unitRNA="pTPM" expRNA="3.0" />
    <level type="RNAExpression"             unitRNA="TPM"  expRNA="2.2" />
    <RNASample sampleId="1987" unitRNA="nTPM" expRNA="2.4" sex="Female" age="91" />
    ...
  </data>
</rnaExpression>
```

* `level`: always the three types above (`nTPM`, `pTPM`, `TPM`), one of each per
  `data` block. These are the aggregated values.
* `RNASample`: attributes observed anywhere in the three files:
  `sampleId, unitRNA, expRNA, sex, age`. Per-sample values are **nTPM only**
  (`unitRNA="nTPM"` in all 1,447 samples per gene). There are no per-sample
  pTPM or TPM values.
* `age` is **absent** (not empty) on 3 `tissue` samples (bone marrow 97,
  placenta 200, skeletal muscle 341) and on **all** `immuneCell` samples. Any
  positional parser (such as v1 of this tool) shifts columns when an attribute
  is missing or added.

### Meaning of `@organ` and `@region` differs between assays

| Assay | Label attribute | Meaning | Example |
|-------|-----------------|---------|---------|
| `tissue`, `consensusTissue` | `@organ` | organ system | `organ="Connective & Soft tissue"`, text `Adipose tissue` |
| `humanBrain` | `@organ` | **brain-region code** (13 codes) | `organ="CTX"`, text `Angular gyrus` |
| `pigBrain`, `mouseBrain`, `*Regional` | `@region` | brain-region code | `region="AMY"`, text `Amygdala` |
| `immuneCell` | `@lineage` | immune lineage | `lineage="T-cells"`, text `MAIT T-cell` |

The `humanBrain` assay stores the region code in `@organ`, while the pig and
mouse assays use `@region` for the same purpose. The parser normalizes this.

Region codes observed:

* human (`humanBrain` / `humanBrainRegional`): AMY, BG, CB, CTX, CP, HPF, HY, M,
  MB, P, SC, TH, WM (13). `humanBrain` has 193 subregions.
* pig (`pigBrain`): AMY, BG, CB, CC, CTX, HPF, HY, MB, OLF, PM, SC, TH (30
  subregions). `pigBrainRegional` uses AMY, BG, CB, CTX, HPF, HY, **M, P, WM**,
  MB, OLF, SC, TH. The codes PM and CC in the per-sample assay do not match
  M/P/WM in the regional assay.
* mouse (`mouseBrain`): AMY, BG, CB, CC, CTX, HPF, HY, MB, OLF, PM, TH (17
  subregions). `mouseBrainRegional` uses WM instead of CC.

Two pig subregion names carry a trailing space in the XML
(`"Dorsal hippocampus "`, `"Ventral hippocampus "`).

### Sample metadata

| Assay | Samples | Sex | Age values |
|-------|--------:|-----|-----------|
| `tissue` | 186 | 105 F / 81 M | 0.1 – 85 (8 samples with `0.1`; 3 missing) |
| `humanBrain` | 966 | 398 F / 568 M | 18, 42, 47 … 94 (21 distinct values) |
| `pigBrain` | 119 | 59 F / 60 M | all `1` |
| `mouseBrain` | 67 | 33 F / 34 M | all `0.2` |
| `immuneCell` | 109 | 56 F / 53 M | absent |

Units of `age` are not stated in the XML. The pig (`1`) and mouse (`0.2`)
values are consistent with years, but this is an interpretation.

`sampleId` is unique per `RNASample` within an assay (966 distinct IDs for 966
`humanBrain` samples). **There is no donor/individual identifier.** Samples
from the same donor cannot be linked, except through the (`sex`, `age`) pair,
which is not guaranteed to be unique. `humanBrain` has 23 distinct (sex, age)
pairs.

## Agreement between per-sample and aggregated values

For every `data` block with samples (tissue, humanBrain, pigBrain,
mouseBrain, immuneCell; all three genes), the **arithmetic mean** of the
`RNASample/@expRNA` values matches the block's
`level[@type="normalizedRNAExpression"]` value with a maximum absolute
difference of **0.075 nTPM**. This is consistent with both values being
rounded to one decimal place (worst case 0.05 + 0.05 = 0.1). The median does
**not** match (differences of up to 55 nTPM).

The 13-region values in `*BrainRegional` are **not** the mean of the samples
in that region. For human brain, in all three genes, the regional value
equals the **maximum** of the subregion nTPM values (for example, ACADM CTX:
regional 5.9, mean of samples 0.88, maximum subregion 5.9). For pig, it is
mostly the maximum, with exceptions: for ACADM, regional MB = 33.0 =
substantia nigra, while superior colliculi = 35.6. For pig and mouse, regions
without a per-sample counterpart (M, P, WM) take the value of a single
subregion (PM *Medulla oblongata*, PM *Pons*, CC *Corpus callosum*).

Consequence: concordance tests must compare samples with the aggregate of
**their own `data` block** (subregion or tissue), not with the 13-region
summary.

## Answers to the Stage 0 questions

**(a) Pig and mouse data.** Yes. The *human* gene XML already contains
`pigBrain` (119 samples, 30 subregions) and `mouseBrain` (67 samples, 17
subregions) with per-sample nTPM, sex and age, plus `pigBrainRegional` and
`mouseBrainRegional` aggregates. The species values are given for the human
gene (via HPA's ortholog mapping). No separate endpoint or ortholog XML is
needed. Caveats: pig age is `1` and mouse age is `0.2` for every sample, so
age is not an analyzable variable for these species; and region codes differ
between species and between per-sample and regional assays (see above).

**(b) Hemisphere.** **No.** There is no hemisphere attribute on `RNASample`
or `tissue`. No subregion name contains left/right/hemisphere (the matches for
"lateral" are anatomical names such as *Lateral amygdala*). `sampleId`s are
bare integers with no encoded side. Hemisphere is therefore not available
from the XML.

**(c) Version.** `entry/@version="25"` (HPA release), `proteinAtlas/@schemaVersion="3.0"`,
`identifier/@version="109"` (Ensembl), `@assembly="GRCh38.p13"`,
`@gencodeVersion="43"`, and the `v25.` subdomain in `entry/@url`. The minor
release (25.1) is not present.

## Other observations relevant to users

* Per-sample values are nTPM only. pTPM and TPM exist only as aggregates.
* `tissue` has 40 tissues with samples, while `consensusTissue` has 51
  (consensus with GTEx, aggregated only).
* Tissue names may carry a numeric suffix (`Endometrium 1`, `Skin 1`,
  `Stomach 1`), kept as published.
* The `cellLine` and `immuneCellLineage` assays contain only summary elements
  (specificity, distribution, cluster) and no per-entity values.
