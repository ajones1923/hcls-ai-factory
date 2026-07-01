"""
Dependency-pin + license gate (A7).

Two checks that keep the Apache-2.0 platform redistributable and reproducible:
  * every dependency is exact-pinned (``==``) — for 21 CFR Part 11 reproducibility;
  * every bundled package is permissively licensed — non-commercial / academic-only licenses
    are blockers, copyleft / ShareAlike are flagged for derivative-redistribution review.

Stdlib-only; usable as a CI gate (`audit(...)["ok"]`) or interactively.
"""
from __future__ import annotations

import re
from typing import Any

_PERMISSIVE = ("mit", "bsd", "apache", "isc", "psf", "python software foundation",
               "zlib", "unlicense", "cc0", "public domain", "mpl-2.0", "mpl 2.0")
_NON_COMMERCIAL = ("noncommercial", "non-commercial", "cc-by-nc", "-nc-", "research only",
                   "research-only", "academic only", "academic-only", "non commercial")
_COPYLEFT = ("agpl", "lgpl", "gpl", "cc-by-sa", "sharealike", "share-alike", "-sa-", "eupl")


def classify_license(license_str: str | None) -> str:
    """-> 'permissive' | 'copyleft' | 'non-commercial' | 'unknown'. Non-commercial wins (most restrictive)."""
    if not license_str:
        return "unknown"
    t = license_str.lower()
    if any(x in t for x in _NON_COMMERCIAL):
        return "non-commercial"
    if any(x in t for x in _COPYLEFT):
        return "copyleft"
    if any(x in t for x in _PERMISSIVE):
        return "permissive"
    return "unknown"


def is_pinned(requirement: str) -> bool:
    """A requirement is reproducible only if exact-pinned (==) or hash-pinned."""
    r = requirement.split("#", 1)[0].strip()
    if not r:
        return True                       # blank/comment lines don't count
    return "==" in r or "@" in r          # 'pkg==1.2.3' or 'pkg @ url' (vcs/url pin)


def check_requirements(requirements: list[str]) -> list[str]:
    """Return UNPINNED issues for a requirements list."""
    issues = []
    for line in requirements:
        s = line.split("#", 1)[0].strip()
        if not s:
            continue
        if not is_pinned(s):
            issues.append(f"UNPINNED: {s}")
    return issues


def audit(packages: list[dict[str, Any]], requirements: list[str] | None = None) -> dict[str, Any]:
    """Audit packages [{name, version, license}] (+ optional requirements) for the A7 gate.

    ``ok`` is False if any package is non-commercial-licensed or any requirement is unpinned.
    Copyleft/unknown are surfaced as warnings, not hard blocks.
    """
    buckets: dict[str, list[str]] = {"permissive": [], "copyleft": [], "non-commercial": [], "unknown": []}
    for p in packages:
        buckets[classify_license(p.get("license"))].append(p.get("name", "?"))
    pin_issues = check_requirements(requirements or [])
    blocking = list(buckets["non-commercial"]) + pin_issues
    return {
        "ok": not blocking,
        "blocking": blocking,
        "warnings": {"copyleft": buckets["copyleft"], "unknown": buckets["unknown"]},
        "by_license": buckets,
        "n_packages": len(packages),
    }
