"""Structural validation warnings.

The parser never fails silently on unexpected XML structure: anything not
observed in the reference release is reported through
:class:`HPAStructureWarning` (and recorded in ``ParsedGene.issues``) and the
data are kept, not dropped or shifted.
"""

from __future__ import annotations

import warnings


class HPAStructureWarning(UserWarning):
    """The XML contains structure that was not observed in HPA v25."""


class HPAVersionMismatchWarning(UserWarning):
    """Genes in one batch come from different HPA releases."""


def warn_structure(message: str, issues: list[str] | None = None) -> None:
    """Emit an :class:`HPAStructureWarning` and record it in ``issues``."""
    if issues is not None:
        issues.append(message)
    warnings.warn(message, HPAStructureWarning, stacklevel=3)


def check_version_consistency(versions: dict[str, str | None]) -> bool:
    """Warn if the genes in a batch do not share one HPA release.

    Args:
        versions: Mapping of gene identifier to its ``hpa_version``.

    Returns:
        ``True`` if all versions agree, ``False`` otherwise.
    """
    distinct = {v for v in versions.values()}
    if len(distinct) <= 1:
        return True
    detail = ", ".join(f"{g}={v}" for g, v in sorted(versions.items()))
    warnings.warn(
        f"Genes in this batch come from different HPA releases: {detail}",
        HPAVersionMismatchWarning,
        stacklevel=2,
    )
    return False
