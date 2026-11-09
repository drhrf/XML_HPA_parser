"""Download, symbol resolution, caching and batches, with the network mocked."""

from __future__ import annotations

import gzip
import json

import pytest
import requests
from conftest import ACADM, ACTB, GFAP

import hpa_xml_parser as hxp
from hpa_xml_parser import fetch as fetch_mod

XML = {
    "ENSG00000117054": gzip.decompress(ACADM.read_bytes()),
    "ENSG00000131095": gzip.decompress(GFAP.read_bytes()),
    "ENSG00000075624": gzip.decompress(ACTB.read_bytes()),
}
SEARCH = {
    "ACADM": [
        {"Gene": "ACADM", "Gene synonym": ["ACAD1", "MCAD", "MCADH"], "Ensembl": "ENSG00000117054"}
    ],
    "MCAD": [
        {"Gene": "ACADM", "Gene synonym": ["ACAD1", "MCAD", "MCADH"], "Ensembl": "ENSG00000117054"}
    ],
    "GFAP": [
        {"Gene": "GFAP", "Gene synonym": [], "Ensembl": "ENSG00000131095"},
        {"Gene": "GFAPX", "Gene synonym": ["GFAP"], "Ensembl": "ENSG99999999999"},
    ],
    "AMBIG": [
        {"Gene": "A1", "Gene synonym": ["AMBIG"], "Ensembl": "ENSG00000000001"},
        {"Gene": "A2", "Gene synonym": ["AMBIG"], "Ensembl": "ENSG00000000002"},
    ],
}


class FakeResponse:
    def __init__(self, status: int, content: bytes):
        self.status_code = status
        self.content = content

    def json(self):
        return json.loads(self.content)

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError(f"{self.status_code}")


class FakeSession:
    """Stands in for requests.Session; records every call."""

    def __init__(self):
        self.calls: list[tuple[str, dict]] = []
        self.headers = {"User-Agent": fetch_mod.USER_AGENT}

    def get(self, url, params=None, timeout=None):
        assert timeout is not None, "every request must set a timeout"
        self.calls.append((url, params or {}))
        if url.endswith("/api/search_download.php"):
            return FakeResponse(200, json.dumps(SEARCH.get(params["search"].upper(), [])).encode())
        ensg = url.rsplit("/", 1)[-1].removesuffix(".xml")
        if ensg in XML:
            return FakeResponse(200, XML[ensg])
        return FakeResponse(404, b"Not found")


@pytest.fixture
def session():
    return FakeSession()


@pytest.fixture(autouse=True)
def no_real_network(monkeypatch):
    def _blocked(*a, **k):
        raise AssertionError("tests must not use the network")

    monkeypatch.setattr(requests.Session, "get", _blocked)


def test_make_session_sets_user_agent_and_retries():
    s = fetch_mod.make_session(retries=5)
    assert s.headers["User-Agent"] == f"hpa-xml-parser/{hxp.__version__}"
    assert s.get_adapter("https://www.proteinatlas.org").max_retries.total == 5


def test_resolve_ensembl_id_without_network(session):
    assert hxp.resolve_gene("ensg00000117054.12", session=session) == "ENSG00000117054"
    assert session.calls == []


def test_resolve_symbol_exact_beats_synonym(session):
    assert hxp.resolve_gene("gfap", session=session) == "ENSG00000131095"
    url, params = session.calls[0]
    assert url == "https://www.proteinatlas.org/api/search_download.php"
    assert params["search"] == "gfap"


def test_resolve_symbol_via_unique_synonym_warns(session):
    with pytest.warns(UserWarning, match="synonym"):
        assert hxp.resolve_gene("MCAD", session=session) == "ENSG00000117054"


def test_resolve_ambiguous_and_missing(session):
    with pytest.raises(hxp.GeneResolutionError, match="ambiguous"):
        hxp.resolve_gene("AMBIG", session=session)
    with pytest.raises(hxp.GeneResolutionError, match="No HPA gene"):
        hxp.resolve_gene("NOPE", session=session)


def test_fetch_gene_uses_cache(tmp_path, session):
    g1 = hxp.fetch_gene("ACADM", session=session, cache_dir=tmp_path)
    n_calls = len(session.calls)
    assert n_calls == 2  # search + download
    assert g1.origin == "https://www.proteinatlas.org/ENSG00000117054.xml"
    assert len(g1.samples) == 1447
    g2 = hxp.fetch_gene("ENSG00000117054", session=session, cache_dir=tmp_path)
    assert len(session.calls) == n_calls  # served from cache
    assert g2.retrieved_at == g1.retrieved_at  # original download time preserved
    cached = list(tmp_path.rglob("*.xml.gz"))
    assert [p.name for p in cached] == ["ENSG00000117054.xml.gz"]
    hxp.fetch_gene("ENSG00000117054", session=session, cache_dir=tmp_path, refresh=True)
    assert len(session.calls) == n_calls + 1


def test_fetch_gene_no_cache(tmp_path, session):
    hxp.fetch_gene("ENSG00000117054", session=session, cache_dir=tmp_path, use_cache=False)
    assert not any(tmp_path.iterdir())


def test_base_url_pins_release(tmp_path, session):
    g = hxp.fetch_gene(
        "ENSG00000117054",
        session=session,
        cache_dir=tmp_path,
        base_url="https://v25.proteinatlas.org",
    )
    assert g.origin == "https://v25.proteinatlas.org/ENSG00000117054.xml"
    assert (tmp_path / "v25.proteinatlas.org" / "ENSG00000117054.xml.gz").exists()


def test_http_error_and_non_xml(tmp_path, session):
    with pytest.raises(requests.HTTPError):
        hxp.fetch_gene("ENSG00000000404", session=session, cache_dir=tmp_path)
    XML["ENSG00000000405"] = b"<html>maintenance</html>"
    try:
        with pytest.raises(ValueError, match="not an HPA XML"):
            hxp.fetch_gene("ENSG00000000405", session=session, cache_dir=tmp_path)
    finally:
        del XML["ENSG00000000405"]


def test_batch_dedupes_and_skips(tmp_path, session):
    with pytest.warns(UserWarning) as record:
        genes = hxp.fetch_genes(
            ["ACADM", "GFAP", "MCAD", "acadm", "NOPE", "ENSG00000075624"],
            session=session,
            cache_dir=tmp_path,
            on_error="warn",
        )
    assert [g.gene_name for g in genes] == ["ACADM", "GFAP", "ACTB"]
    text = " ".join(str(w.message) for w in record)
    assert "Skipping 'NOPE'" in text
    assert "same gene" in text
    with pytest.raises(hxp.GeneResolutionError):
        hxp.fetch_genes(["NOPE"], session=session, cache_dir=tmp_path)


def test_read_gene_list(tmp_path):
    p = tmp_path / "genes.txt"
    p.write_text("# my genes\nACADM\nGFAP, ACTB  # controls\n\nENSG00000117054\tMCAD\n")
    assert hxp.read_gene_list(p) == ["ACADM", "GFAP", "ACTB", "ENSG00000117054", "MCAD"]
