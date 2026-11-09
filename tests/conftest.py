"""Shared fixtures: HPA v25 per-gene XML files downloaded on 2026-11-09."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

import hpa_xml_parser as hxp

FIXTURES = Path(__file__).parent / "fixtures"
ACADM = FIXTURES / "ENSG00000117054_ACADM_v25.xml.gz"
GFAP = FIXTURES / "ENSG00000131095_GFAP_v25.xml.gz"
ACTB = FIXTURES / "ENSG00000075624_ACTB_v25.xml.gz"
ALL_FIXTURES = (ACADM, GFAP, ACTB)
RETRIEVED_AT = "2026-11-09T16:23:00Z"

#: Samples per assay observed in HPA v25 (identical for all three genes).
EXPECTED_SAMPLES = {
    "tissue": 186,
    "humanBrain": 966,
    "pigBrain": 119,
    "mouseBrain": 67,
    "immuneCell": 109,
}


@pytest.fixture(scope="session")
def acadm() -> hxp.ParsedGene:
    return hxp.parse_xml(ACADM, retrieved_at=RETRIEVED_AT)


@pytest.fixture(scope="session")
def all_genes() -> list[hxp.ParsedGene]:
    return [hxp.parse_xml(p, retrieved_at=RETRIEVED_AT) for p in ALL_FIXTURES]


@pytest.fixture(scope="session")
def acadm_bytes() -> bytes:
    return gzip.decompress(ACADM.read_bytes())
