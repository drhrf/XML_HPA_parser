"""Robustness against XML structure not seen in HPA v25."""

from __future__ import annotations

import math
import warnings

import pandas as pd
import pytest

import hpa_xml_parser as hxp
from hpa_xml_parser import HPAStructureWarning, HPAVersionMismatchWarning

SYNTHETIC = b"""<?xml version="1.0" encoding="UTF-8"?>
<proteinAtlas schemaVersion="9.9">
  <entry version="%(version)s" url="https://v99.proteinatlas.org/ENSG00000000001">
    <name>TEST1</name>
    <identifier id="ENSG00000000001" db="Ensembl" version="999" assembly="GRCh99" gencodeVersion="99"/>
    <rnaExpression source="HPA" technology="RNAseq" assayType="humanBrain">
      <data>
        <tissue hemisphereCode="L" ontologyTerms="UBERON:1" organ="CTX">Test cortex</tissue>
        <level expRNA="2.0" unitRNA="nTPM" type="normalizedRNAExpression"/>
        <RNASample age="30" hemisphere="left" sex="Male" expRNA="1.5" unitRNA="nTPM" sampleId="S1"/>
        <RNASample sampleId="S2" sex="Female" unitRNA="nTPM" expRNA="2.5"/>
      </data>
    </rnaExpression>
    <rnaExpression source="HPA" technology="snRNAseq" assayType="newAssay">
      <data>
        <tissue region="XX">New region</tissue>
        <level type="newLevelType" unitRNA="nTPM" expRNA="7"/>
        <RNASample sampleId="N1" unitRNA="nTPM" expRNA="7" sex="Male" age="40"/>
      </data>
    </rnaExpression>
  </entry>
</proteinAtlas>
"""


def synthetic(version: str = "99") -> bytes:
    return SYNTHETIC % {b"version": version.encode()}


def test_extra_attribute_out_of_order_is_kept_not_shifted():
    with pytest.warns(HPAStructureWarning) as record:
        gene = hxp.parse_xml(synthetic())
    messages = " | ".join(str(w.message) for w in record)
    assert "hemisphere" in messages
    assert "newAssay" in messages
    assert "newLevelType" in messages
    df = gene.samples.set_index("sample_id")
    s1 = df.loc["S1"]
    assert (s1["value"], s1["sex"], s1["age"], s1["unit"]) == (1.5, "Male", 30.0, "nTPM")
    assert s1["extra_hemisphere"] == "left"
    assert s1["extra_label_hemisphereCode"] == "L"
    assert (s1["organ"], s1["region"], s1["ontology_terms"]) == ("CTX", "Test cortex", "UBERON:1")
    s2 = df.loc["S2"]
    assert (s2["value"], s2["sex"]) == (2.5, "Female")
    assert math.isnan(s2["age"])
    assert pd.isna(s2["extra_hemisphere"])
    # Unknown assay is still parsed.
    assert df.loc["N1", "assay_type"] == "newAssay"
    assert df.loc["N1", "technology"] == "snRNAseq"
    # Core columns come first, extras after.
    assert list(gene.samples.columns[:15]) == list(hxp.SAMPLE_COLUMNS)
    assert gene.issues


def test_version_fields_from_synthetic():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gene = hxp.parse_xml(synthetic())
    assert gene.hpa_version == "99"
    assert gene.schema_version == "9.9"
    assert gene.ensembl_version == "999"


def test_version_mismatch_in_batch_warns(acadm, tmp_path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", HPAStructureWarning)
        other = hxp.parse_xml(synthetic("26"))
    with pytest.warns(HPAVersionMismatchWarning, match="different HPA releases"):
        hxp.export([acadm, other], tmp_path / "mixed.csv")


@pytest.mark.filterwarnings("ignore::hpa_xml_parser.HPAStructureWarning")
def test_missing_entry_version_warns():
    xml = synthetic().replace(b'version="99" url', b"url")
    with pytest.warns(HPAStructureWarning, match="entry/@version"):
        gene = hxp.parse_xml(xml)
    assert gene.hpa_version is None
