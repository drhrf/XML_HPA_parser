"""Download HPA per-gene XML files, with gene-symbol resolution and a local cache."""

from __future__ import annotations

import gzip
import json
import os
import re
import warnings
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from ._version import __version__
from .parser import ParsedGene, _utc_iso, parse_xml
from .validation import check_version_consistency

#: Default HPA base URL. Use e.g. ``https://v25.proteinatlas.org`` to pin a release.
DEFAULT_BASE_URL = "https://www.proteinatlas.org"
#: Generic User-Agent sent with every request.
USER_AGENT = f"hpa-xml-parser/{__version__}"
DEFAULT_TIMEOUT = 60.0
DEFAULT_RETRIES = 3

_ENSG_RE = re.compile(r"^ENSG\d{11}(\.\d+)?$", re.IGNORECASE)


class GeneResolutionError(LookupError):
    """A gene symbol could not be resolved to exactly one Ensembl ID."""


def default_cache_dir() -> Path:
    """Return the cache directory (``$HPA_XML_PARSER_CACHE`` or ``~/.cache/hpa-xml-parser``)."""
    env = os.environ.get("HPA_XML_PARSER_CACHE")
    if env:
        return Path(env).expanduser()
    xdg = os.environ.get("XDG_CACHE_HOME")
    base = Path(xdg).expanduser() if xdg else Path.home() / ".cache"
    return base / "hpa-xml-parser"


def make_session(retries: int = DEFAULT_RETRIES) -> requests.Session:
    """Create a :class:`requests.Session` with retries and the package User-Agent."""
    session = requests.Session()
    retry = Retry(
        total=retries,
        connect=retries,
        read=retries,
        status=retries,
        backoff_factor=1.0,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset({"GET"}),
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    session.headers["User-Agent"] = USER_AGENT
    return session


def is_ensembl_id(identifier: str) -> bool:
    """Return ``True`` if ``identifier`` looks like a human Ensembl gene ID."""
    return bool(_ENSG_RE.match(identifier.strip()))


def resolve_gene(
    identifier: str,
    *,
    session: requests.Session | None = None,
    base_url: str = DEFAULT_BASE_URL,
    timeout: float = DEFAULT_TIMEOUT,
) -> str:
    """Resolve a gene symbol (or Ensembl ID) to an Ensembl gene ID via the HPA search API.

    An exact (case-insensitive) match on the HPA gene name wins. Otherwise a
    unique match on a gene synonym is accepted.

    Args:
        identifier: Ensembl ID (returned unchanged, version suffix removed) or gene symbol.
        session: Optional session (see :func:`make_session`).
        base_url: HPA base URL.
        timeout: Request timeout in seconds.

    Raises:
        GeneResolutionError: If no gene or more than one gene matches.
    """
    ident = identifier.strip()
    if is_ensembl_id(ident):
        return ident.split(".")[0].upper()
    session = session or make_session()
    resp = session.get(
        f"{base_url}/api/search_download.php",
        params={"search": ident, "format": "json", "columns": "g,gs,eg", "compress": "no"},
        timeout=timeout,
    )
    resp.raise_for_status()
    hits: list[dict[str, Any]] = resp.json() or []
    exact = [h for h in hits if str(h.get("Gene", "")).upper() == ident.upper()]
    if len(exact) == 1:
        return str(exact[0]["Ensembl"])
    if not exact:
        syn = [
            h
            for h in hits
            if ident.upper() in {str(s).upper() for s in (h.get("Gene synonym") or [])}
        ]
        if len(syn) == 1:
            warnings.warn(
                f"{ident!r} resolved via synonym to {syn[0].get('Gene')} ({syn[0]['Ensembl']})",
                UserWarning,
                stacklevel=2,
            )
            return str(syn[0]["Ensembl"])
        exact = syn
    if not exact:
        raise GeneResolutionError(f"No HPA gene matches {ident!r}")
    found = ", ".join(f"{h.get('Gene')} ({h.get('Ensembl')})" for h in exact)
    raise GeneResolutionError(f"{ident!r} is ambiguous in HPA: {found}; use an Ensembl ID")


def _cache_paths(cache_dir: Path, ensembl_id: str, base_url: str) -> tuple[Path, Path]:
    host = re.sub(r"^https?://", "", base_url).strip("/").replace("/", "_")
    folder = cache_dir / host
    return folder / f"{ensembl_id}.xml.gz", folder / f"{ensembl_id}.json"


def download_xml(
    ensembl_id: str,
    *,
    session: requests.Session | None = None,
    base_url: str = DEFAULT_BASE_URL,
    timeout: float = DEFAULT_TIMEOUT,
    cache_dir: str | os.PathLike[str] | None = None,
    use_cache: bool = True,
    refresh: bool = False,
) -> tuple[bytes, dict[str, Any]]:
    """Download (or read from cache) the XML of one gene.

    Returns:
        A tuple ``(xml_bytes, info)`` where ``info`` holds ``url``,
        ``retrieved_at`` and ``from_cache``.
    """
    url = f"{base_url}/{ensembl_id}.xml"
    cdir = Path(cache_dir) if cache_dir is not None else default_cache_dir()
    xml_path, meta_path = _cache_paths(cdir, ensembl_id, base_url)
    if use_cache and not refresh and xml_path.exists() and meta_path.exists():
        info = json.loads(meta_path.read_text(encoding="utf-8"))
        info["from_cache"] = True
        return gzip.decompress(xml_path.read_bytes()), info

    session = session or make_session()
    resp = session.get(url, timeout=timeout)
    resp.raise_for_status()
    content = resp.content
    if b"<proteinAtlas" not in content[:2000] and b"<entry" not in content[:2000]:
        raise ValueError(f"Response from {url} is not an HPA XML document")
    info = {"url": url, "retrieved_at": _utc_iso(), "from_cache": False}
    if use_cache:
        xml_path.parent.mkdir(parents=True, exist_ok=True)
        xml_path.write_bytes(gzip.compress(content))
        meta_path.write_text(
            json.dumps({k: info[k] for k in ("url", "retrieved_at")}), encoding="utf-8"
        )
    return content, info


def fetch_gene(
    identifier: str,
    *,
    session: requests.Session | None = None,
    base_url: str = DEFAULT_BASE_URL,
    timeout: float = DEFAULT_TIMEOUT,
    retries: int = DEFAULT_RETRIES,
    cache_dir: str | os.PathLike[str] | None = None,
    use_cache: bool = True,
    refresh: bool = False,
) -> ParsedGene:
    """Fetch and parse one gene given an Ensembl ID or gene symbol.

    ``retrieved_at`` is the time of the original download, also when the XML
    is served from the cache.
    """
    session = session or make_session(retries)
    ensembl_id = resolve_gene(identifier, session=session, base_url=base_url, timeout=timeout)
    content, info = download_xml(
        ensembl_id,
        session=session,
        base_url=base_url,
        timeout=timeout,
        cache_dir=cache_dir,
        use_cache=use_cache,
        refresh=refresh,
    )
    gene = parse_xml(content, retrieved_at=info["retrieved_at"])
    gene.origin = info["url"]
    return gene


def read_gene_list(path: str | os.PathLike[str]) -> list[str]:
    """Read gene identifiers from a text file.

    Identifiers may be separated by newlines, commas, tabs or spaces. Blank
    lines and text after ``#`` are ignored.
    """
    genes: list[str] = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        line = line.split("#", 1)[0]
        genes.extend(tok for tok in re.split(r"[\s,;]+", line) if tok)
    return genes


def fetch_genes(
    identifiers: Iterable[str],
    *,
    on_error: str = "raise",
    **kwargs: Any,
) -> list[ParsedGene]:
    """Fetch several genes, sharing one HTTP session.

    Args:
        identifiers: Ensembl IDs and/or gene symbols. Duplicates are fetched once.
        on_error: ``"raise"`` (default) or ``"warn"`` to skip failing genes.
        **kwargs: Passed to :func:`fetch_gene`.

    Warns:
        HPAVersionMismatchWarning: If genes come from different HPA releases.
    """
    if on_error not in ("raise", "warn"):
        raise ValueError("on_error must be 'raise' or 'warn'")
    if kwargs.get("session") is None:
        kwargs["session"] = make_session(kwargs.get("retries", DEFAULT_RETRIES))
    genes: list[ParsedGene] = []
    seen: set[str] = set()
    for ident in identifiers:
        key = ident.strip().upper()
        if not key or key in seen:
            continue
        seen.add(key)
        try:
            gene = fetch_gene(ident, **kwargs)
        except Exception as exc:  # noqa: BLE001 - reported to the user
            if on_error == "raise":
                raise
            warnings.warn(f"Skipping {ident!r}: {exc}", UserWarning, stacklevel=2)
            continue
        if any(g.ensembl_id == gene.ensembl_id for g in genes):
            warnings.warn(
                f"{ident!r} is the same gene as an earlier entry ({gene.ensembl_id}); skipped",
                UserWarning,
                stacklevel=2,
            )
            continue
        genes.append(gene)
    check_version_consistency({g.gene_name or g.ensembl_id or "?": g.hpa_version for g in genes})
    return genes
