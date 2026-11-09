"""Export formats and the command-line interface."""

from __future__ import annotations

import pandas as pd
import pytest
from conftest import ACADM, ACTB, GFAP

import hpa_xml_parser as hxp
from hpa_xml_parser import cli
from hpa_xml_parser.schema import SAMPLE_ASSAYS, SAMPLE_COLUMNS


def test_to_frame_batch_and_filter(all_genes):
    df = hxp.to_frame(all_genes)
    assert len(df) == 3 * 1447
    assert df.groupby("gene_name").size().to_dict() == {"ACADM": 1447, "ACTB": 1447, "GFAP": 1447}
    brain = hxp.to_frame(all_genes, assays=["humanBrain"])
    assert set(brain["assay_type"]) == {"humanBrain"} and len(brain) == 3 * 966
    agg = hxp.to_frame(all_genes, kind="aggregates")
    assert set(agg["gene_name"]) == {"ACADM", "GFAP", "ACTB"}
    with pytest.raises(ValueError):
        hxp.to_frame(all_genes, kind="other")


def test_xlsx_sheets(tmp_path, acadm):
    out = tmp_path / "ACADM.xlsx"
    assert hxp.export(acadm, out) == [out]
    xl = pd.ExcelFile(out)
    assert xl.sheet_names == list(SAMPLE_ASSAYS) + ["aggregates", "metadata"]
    brain = pd.read_excel(xl, "humanBrain", dtype={"sample_id": str})
    assert len(brain) == 966
    assert list(brain.columns[: len(SAMPLE_COLUMNS)]) == list(SAMPLE_COLUMNS)
    meta = pd.read_excel(xl, "metadata", dtype=str)
    assert meta.loc[0, "hpa_version"] == "25"
    assert meta.loc[0, "parser_version"] == hxp.__version__
    assert meta.loc[0, "ensembl_version"] == "109"


def test_xlsx_assay_selection(tmp_path, acadm):
    out = tmp_path / "sel.xlsx"
    hxp.export(acadm, out, assays=["humanBrain", "tissue"], include_aggregates=False)
    assert pd.ExcelFile(out).sheet_names == ["humanBrain", "tissue", "metadata"]


@pytest.mark.parametrize("ext,sep", [(".csv", ","), (".tsv", "\t")])
def test_flat_files(tmp_path, acadm, ext, sep):
    out = tmp_path / f"ACADM{ext}"
    written = hxp.export(acadm, out)
    assert [p.name for p in written] == [
        f"ACADM{ext}",
        f"ACADM_aggregates{ext}",
        f"ACADM_metadata{ext}",
    ]
    back = pd.read_csv(out, sep=sep, dtype={"sample_id": str, "hpa_version": str})
    assert len(back) == 1447
    pd.testing.assert_series_equal(back["value"], acadm.samples["value"].reset_index(drop=True))
    assert (back["hpa_version"] == "25").all()


def test_parquet_roundtrip(tmp_path, acadm):
    pytest.importorskip("pyarrow")
    out = tmp_path / "ACADM.parquet"
    hxp.export(acadm, out)
    back = pd.read_parquet(out)
    assert len(back) == 1447
    assert (tmp_path / "ACADM_aggregates.parquet").exists()


def test_unknown_format(tmp_path, acadm):
    with pytest.raises(ValueError, match="Unsupported"):
        hxp.export(acadm, tmp_path / "x.json")


def test_cli_parse(tmp_path, caplog):
    out = tmp_path / "res.xlsx"
    code = cli.main(
        ["parse", str(ACADM), str(GFAP), str(ACTB), "--assay", "humanBrain", "--out", str(out)]
    )
    assert code == 0
    df = pd.read_excel(out, "humanBrain")
    assert len(df) == 3 * 966


def test_cli_fetch_uses_api(tmp_path, monkeypatch, acadm):
    seen = {}

    def fake_fetch_genes(ids, **kwargs):
        seen["ids"] = ids
        seen.update(kwargs)
        return [acadm]

    monkeypatch.setattr(cli, "fetch_genes", fake_fetch_genes)
    genes_file = tmp_path / "genes.txt"
    genes_file.write_text("GFAP\n")
    out = tmp_path / "r.csv"
    code = cli.main(
        ["fetch", "ACADM", "--genes-file", str(genes_file), "--out", str(out), "--no-cache",
         "--timeout", "5", "--retries", "1", "--skip-errors"]
    )  # fmt: skip
    assert code == 0
    assert seen["ids"] == ["ACADM", "GFAP"]
    assert seen["use_cache"] is False and seen["timeout"] == 5 and seen["retries"] == 1
    assert seen["on_error"] == "warn"
    assert out.exists()


def test_cli_errors(tmp_path):
    assert cli.main(["fetch", "--out", str(tmp_path / "x.csv")]) == 2
    assert cli.main(["parse", str(tmp_path / "missing.xml"), "--out", str(tmp_path / "x.csv")]) == 1


def test_cli_version(capsys):
    with pytest.raises(SystemExit):
        cli.main(["--version"])
    assert hxp.__version__ in capsys.readouterr().out
