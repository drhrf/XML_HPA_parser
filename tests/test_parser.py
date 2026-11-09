"""Parsing of real HPA v25 fixtures."""

from __future__ import annotations

import gzip
import io
import math
import warnings
import xml.etree.ElementTree as ET

import pandas as pd
import pytest
from conftest import ACADM, ALL_FIXTURES, EXPECTED_SAMPLES, RETRIEVED_AT

import hpa_xml_parser as hxp
from hpa_xml_parser.schema import SAMPLE_COLUMNS, SAMPLE_EXTRA_COLUMNS


def test_metadata(acadm):
    assert acadm.gene_name == "ACADM"
    assert acadm.ensembl_id == "ENSG00000117054"
    assert acadm.hpa_version == "25"
    assert acadm.schema_version == "3.0"
    assert acadm.ensembl_version == "109"
    assert acadm.assembly == "GRCh38.p13"
    assert acadm.gencode_version == "43"
    assert acadm.entry_url == "https://v25.proteinatlas.org/ENSG00000117054"
    assert acadm.retrieved_at == RETRIEVED_AT
    assert acadm.issues == []


@pytest.mark.parametrize("path", ALL_FIXTURES, ids=lambda p: p.name.split("_")[1])
def test_sample_counts_per_assay(path):
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # real v25 data must not trigger structure warnings
        gene = hxp.parse_xml(path)
    counts = gene.samples.groupby("assay_type").size().to_dict()
    assert counts == EXPECTED_SAMPLES
    assert len(gene.samples) == sum(EXPECTED_SAMPLES.values()) == 1447


def test_column_order_and_dtypes(acadm):
    df = acadm.samples
    assert tuple(df.columns) == SAMPLE_COLUMNS + SAMPLE_EXTRA_COLUMNS
    assert df["value"].dtype == "float64"
    assert df["age"].dtype == "float64"
    for col in ("sample_id", "sex", "unit", "region", "organ", "hpa_version"):
        assert isinstance(df[col].dtype, pd.StringDtype), col


def test_value_ranges(acadm):
    df = acadm.samples
    assert set(df["sex"].dropna()) == {"Female", "Male"}
    assert set(df["unit"]) == {"nTPM"}
    assert (df["value"] >= 0).all()
    assert df["value"].notna().all()
    assert set(df["species"]) == {"human", "pig", "mouse"}
    brain = df[df["assay_type"] == "humanBrain"]
    assert brain["age"].between(18, 94).all()
    assert set(df.loc[df["assay_type"] == "pigBrain", "age"]) == {1.0}
    assert set(df.loc[df["assay_type"] == "mouseBrain", "age"]) == {0.2}
    assert df.loc[df["assay_type"] == "immuneCell", "age"].isna().all()


def test_known_rows(acadm):
    df = acadm.samples.set_index(["assay_type", "sample_id"])
    amy = df.loc[("humanBrain", "1987")]
    assert (amy["region"], amy["organ"], amy["value"], amy["sex"], amy["age"]) == (
        "Amygdala",
        "AMY",
        2.4,
        "Female",
        91.0,
    )
    assert amy["ontology_terms"] == "UBERON:0001876"
    assert amy["label_attribute"] == "tissue@organ"
    # Sample without an age attribute: sex must not be shifted into age.
    bm = df.loc[("tissue", "97")]
    assert (bm["region"], bm["value"], bm["sex"]) == ("Bone marrow", 5.2, "Female")
    assert math.isnan(bm["age"])
    pig = df.loc[("pigBrain", "1052")]
    assert (pig["organ"], pig["region"], pig["label_attribute"]) == (
        "AMY",
        "Amygdala",
        "tissue@region",
    )
    imm = df.loc[("immuneCell", "976")]
    assert (imm["organ"], imm["region"], imm["value"]) == ("Granulocytes", "Basophil", 18.1)
    assert imm["label_attribute"] == "immuneCell@lineage"


def test_no_column_shift_against_raw_xml(acadm, acadm_bytes):
    """Every parsed row must equal the attributes of its RNASample, read by name."""
    root = ET.fromstring(acadm_bytes)
    parsed = acadm.samples.set_index(["assay_type", "sample_id"])
    n = 0
    for rx in root.iter("rnaExpression"):
        for rs in rx.iter("RNASample"):
            row = parsed.loc[(rx.get("assayType"), rs.get("sampleId"))]
            assert row["value"] == float(rs.get("expRNA"))
            assert row["unit"] == rs.get("unitRNA")
            assert row["sex"] == rs.get("sex")
            if rs.get("age") is None:
                assert math.isnan(row["age"])
            else:
                assert row["age"] == float(rs.get("age"))
            n += 1
    assert n == len(acadm.samples)


def test_aggregates(acadm):
    agg = acadm.aggregates
    assert set(agg["level_type"]) == {
        "normalizedRNAExpression",
        "proteinCodingRNAExpression",
        "RNAExpression",
    }
    by_unit = agg.groupby("level_type")["unit"].unique().to_dict()
    assert list(by_unit["normalizedRNAExpression"]) == ["nTPM"]
    regional = agg[(agg["assay_type"] == "humanBrainRegional")]
    assert regional["region"].nunique() == 13
    assert (regional["n_samples"] == 0).all()
    nt = agg[(agg["assay_type"] == "humanBrain") & (agg["level_type"] == "normalizedRNAExpression")]
    assert len(nt) == 193
    assert nt["n_samples"].sum() == EXPECTED_SAMPLES["humanBrain"]


def test_assay_summary(acadm):
    assays = acadm.assays.set_index("assay_type")
    assert set(assays.index) == {
        "consensusTissue", "tissue", "humanBrain", "humanBrainRegional", "mouseBrainRegional",
        "mouseBrain", "pigBrainRegional", "pigBrain", "cellLine", "immuneCell",
        "immuneCellLineage",
    }  # fmt: skip
    assert (assays["source"] == "HPA").all()
    assert (assays["technology"] == "RNAseq").all()


def test_gz_plain_bytes_and_fileobj_identical(tmp_path, acadm, acadm_bytes):
    plain = tmp_path / "ACADM.xml"
    plain.write_bytes(acadm_bytes)
    renamed_gz = tmp_path / "ACADM_no_suffix"
    renamed_gz.write_bytes(gzip.compress(acadm_bytes))
    for src in (plain, str(plain), acadm_bytes, io.BytesIO(acadm_bytes), renamed_gz, ACADM):
        gene = hxp.parse_xml(src, retrieved_at=RETRIEVED_AT)
        pd.testing.assert_frame_equal(gene.samples, acadm.samples)
        pd.testing.assert_frame_equal(gene.aggregates, acadm.aggregates)


def test_retrieved_at_defaults_to_file_mtime(tmp_path, acadm_bytes):
    import os

    p = tmp_path / "x.xml"
    p.write_bytes(acadm_bytes)
    os.utime(p, (1_700_000_000, 1_700_000_000))
    assert hxp.parse_xml(p).retrieved_at == "2023-11-14T22:13:20Z"


def test_not_hpa_xml():
    with pytest.raises(ValueError, match="no <entry>"):
        hxp.parse_xml(b"<root><foo/></root>")
