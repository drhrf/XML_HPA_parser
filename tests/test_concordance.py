"""Internal validation: per-sample values reproduce HPA's own aggregates.

For every ``<data>`` block (one tissue, brain subregion or immune cell type)
the arithmetic mean of the per-sample nTPM values must equal the block's
``level[@type="normalizedRNAExpression"]`` value.

Tolerance: both the per-sample values and the aggregate are published with one
decimal place, so each carries a rounding error of at most 0.05 nTPM. The mean
of rounded values therefore differs from the rounded mean by at most
0.05 + 0.05 = 0.1 nTPM. The largest difference observed in HPA v25 for ACADM,
GFAP and ACTB is 0.075 nTPM.

Note: the 13-region values of ``humanBrainRegional`` are NOT means of the
samples. In HPA v25 they equal the maximum subregion value (tested below), so
concordance is checked at the level of the ``<data>`` block.
"""

from __future__ import annotations

import pytest
from conftest import ALL_FIXTURES

import hpa_xml_parser as hxp
from hpa_xml_parser.schema import SAMPLE_ASSAYS

TOLERANCE_NTPM = 0.1
KEYS = ["assay_type", "organ", "region"]


def _concordance(gene: hxp.ParsedGene):
    means = gene.samples.groupby(KEYS)["value"].agg(["mean", "size"]).reset_index()
    agg = gene.aggregates
    agg = agg[(agg["level_type"] == "normalizedRNAExpression") & (agg["n_samples"] > 0)]
    merged = means.merge(agg[KEYS + ["value", "n_samples"]], on=KEYS, how="outer", indicator=True)
    return merged


@pytest.mark.parametrize("path", ALL_FIXTURES, ids=lambda p: p.name.split("_")[1])
def test_sample_mean_matches_hpa_aggregate(path):
    merged = _concordance(hxp.parse_xml(path))
    assert (merged["_merge"] == "both").all(), (
        "every data block must pair samples with an aggregate"
    )
    assert (merged["size"] == merged["n_samples"]).all()
    diff = (merged["mean"] - merged["value"]).abs()
    assert diff.max() <= TOLERANCE_NTPM, merged.loc[diff.idxmax()].to_dict()
    assert set(merged["assay_type"]) == set(SAMPLE_ASSAYS)


def test_acadm_concordance_summary(acadm):
    merged = _concordance(acadm)
    diff = (merged["mean"] - merged["value"]).abs()
    assert len(merged) == 40 + 193 + 30 + 17 + 19  # tissues + subregions + pig + mouse + immune
    assert round(diff.max(), 3) == 0.075


@pytest.mark.parametrize("path", ALL_FIXTURES, ids=lambda p: p.name.split("_")[1])
def test_human_regional_value_is_max_of_subregions(path):
    """Documents HPA behavior the manuscript must not misdescribe as a mean."""
    agg = hxp.parse_xml(path).aggregates
    nt = agg[agg["level_type"] == "normalizedRNAExpression"]
    sub = nt[nt["assay_type"] == "humanBrain"].groupby("organ")["value"].max()
    reg = nt[nt["assay_type"] == "humanBrainRegional"].set_index("organ")["value"]
    assert len(reg) == 13
    assert (sub.reindex(reg.index) == reg).all()
