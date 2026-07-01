"""
Security utilities for the HCLS AI Factory precision medicine platform.

Consolidated from:
  - core/engines/genomic-foundation/web-portal/app/security.py   (headers, rate limiting)
  - core/engines/genomic-foundation/web-portal/app/validation.py  (input validation)

All functions are framework-agnostic except ``add_security_headers`` which
operates on a Flask ``Response`` object.

Provides:
  - sanitize_search_query         -- strip dangerous chars from free-text search
  - validate_milvus_filter        -- token-whitelist validation for Milvus boolean exprs
  - sanitize_gene_name            -- HGNC gene name validation
  - sanitize_chromosome           -- chromosome identifier validation
  - validate_patient_id           -- alphanumeric patient/sample ID
  - validate_smiles               -- SMILES string validation (RDKit optional)
  - validate_pdb_id               -- PDB identifier validation
  - rate_limit_check              -- sliding-window rate limiter
  - add_security_headers          -- Flask response security headers
  - sanitize_filename             -- safe filename for file outputs
"""

from __future__ import annotations

import logging
import os
import re
import threading
import time
from collections import defaultdict, deque
from pathlib import Path
from typing import Any, Deque, Dict, Optional, Tuple

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Regex patterns
# ---------------------------------------------------------------------------

# Gene names: standard HGNC symbols (letters, digits, hyphens, underscores)
_GENE_RE = re.compile(r"^[A-Za-z0-9_\-]+$")

# Chromosome: optional "chr" prefix + 1-22/X/Y/M/MT
_CHROM_RE = re.compile(r"^(chr)?([0-9]{1,2}|[XYMT]{1,2})$", re.IGNORECASE)

# Patient / sample IDs: start with alphanum, then alphanum/underscore/hyphen
_PATIENT_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_\-]*$")

# PDB identifiers: 4 alphanumeric characters
_PDB_ID_RE = re.compile(r"^[A-Za-z0-9]{4}$")

# Safe filename characters
_SAFE_FILENAME_RE = re.compile(r"[^A-Za-z0-9_\-.]")

# Dangerous shell metacharacters
_SHELL_META_RE = re.compile(r"[;&|`$\\]")

# Path traversal
_PATH_TRAVERSAL_RE = re.compile(r"\.\.")

# Milvus filter token whitelist (same as milvus_client.py -- kept in sync)
_MILVUS_TOKEN_RE = re.compile(
    r"""
    (?:
      \d+(?:\.\d+)?
    | "[^"]*"
    | '[^']*'
    | ==|!=|>=|<=|>|<
    | \b(?:and|or|not|in|AND|OR|NOT|IN)\b
    | \b(?:gene|chrom|pos|impact|consequence|clinical_significance
          |rsid|am_class|am_pathogenicity|qual)\b
    | [()\[\],]
    | \s+
    )
    """,
    re.VERBOSE,
)

# ---------------------------------------------------------------------------
# Free-text search sanitization
# ---------------------------------------------------------------------------


def sanitize_search_query(query: str, max_length: int = 500) -> str:
    """
    Sanitize a free-text search query.

    - Strips leading/trailing whitespace
    - Removes shell metacharacters
    - Truncates to ``max_length``
    - Rejects empty queries after sanitization
    """
    if not query:
        raise ValueError("Search query must not be empty")

    query = query.strip()
    query = _SHELL_META_RE.sub("", query)
    query = _PATH_TRAVERSAL_RE.sub("", query)
    query = query[:max_length]

    if not query.strip():
        raise ValueError("Search query contains only invalid characters")

    return query


# ---------------------------------------------------------------------------
# Milvus filter validation
# ---------------------------------------------------------------------------


def validate_milvus_filter(expr: str) -> str:
    """
    Validate a Milvus boolean filter expression using a token whitelist.

    Only allows known field names (gene, chrom, pos, ...), comparison
    operators, logical keywords, numeric and string literals, and grouping
    punctuation.

    Raises ``ValueError`` on disallowed tokens.
    """
    if not expr or not expr.strip():
        raise ValueError("Filter expression is empty")

    residual = _MILVUS_TOKEN_RE.sub("", expr).strip()
    if residual:
        raise ValueError(
            f"Milvus filter expression contains disallowed tokens: {residual!r}"
        )
    return expr


# ---------------------------------------------------------------------------
# Gene / chromosome / patient validation
# ---------------------------------------------------------------------------


def sanitize_gene_name(gene: str) -> str:
    """Validate and return a gene name. Raises ``ValueError`` on bad input."""
    gene = gene.strip()
    if not gene:
        raise ValueError("Gene name must not be empty")
    if len(gene) > 50:
        raise ValueError(f"Gene name too long ({len(gene)} chars, max 50)")
    if not _GENE_RE.match(gene):
        raise ValueError(f"Invalid gene name: {gene!r}")
    return gene


def sanitize_chromosome(chrom: str) -> str:
    """Validate and return a chromosome identifier."""
    chrom = chrom.strip()
    if not chrom:
        raise ValueError("Chromosome must not be empty")
    if not _CHROM_RE.match(chrom):
        raise ValueError(f"Invalid chromosome: {chrom!r}")
    return chrom


def validate_patient_id(patient_id: str) -> Tuple[bool, Optional[str]]:
    """
    Validate a patient/sample ID.

    Returns ``(True, None)`` if valid, ``(False, error_message)`` if not.
    """
    if not patient_id:
        return False, "Patient ID is required"

    patient_id = patient_id.strip()

    if not _PATIENT_ID_RE.match(patient_id):
        return False, "Invalid patient ID format (alphanumeric, underscores, hyphens)"

    if len(patient_id) > 50:
        return False, "Patient ID too long (max 50 characters)"

    return True, None


# ---------------------------------------------------------------------------
# SMILES validation (RDKit optional)
# ---------------------------------------------------------------------------


def validate_smiles(smiles: str) -> Tuple[bool, Optional[str]]:
    """
    Validate a SMILES string.

    If RDKit is installed, uses ``Chem.MolFromSmiles`` for chemical
    validation.  Otherwise performs basic character-set checks.

    Returns ``(True, None)`` if valid, ``(False, error_message)`` if not.
    """
    if not smiles or not smiles.strip():
        return False, "SMILES string must not be empty"

    smiles = smiles.strip()

    if len(smiles) > 5000:
        return False, "SMILES string too long (max 5000 characters)"

    # Reject obvious injection attempts
    if _SHELL_META_RE.search(smiles):
        return False, "SMILES contains disallowed characters"

    try:
        from rdkit import Chem

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, f"RDKit could not parse SMILES: {smiles!r}"
        return True, None
    except ImportError:
        # Basic character-set check (no RDKit available)
        allowed = set(
            "abcdefghijklmnopqrstuvwxyz"
            "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
            "0123456789"
            "()[]=#@+-/\\.%"
        )
        invalid = set(smiles) - allowed
        if invalid:
            return False, f"SMILES contains unexpected characters: {''.join(sorted(invalid))}"
        return True, None


# ---------------------------------------------------------------------------
# PDB ID validation
# ---------------------------------------------------------------------------


def validate_pdb_id(pdb_id: str) -> Tuple[bool, Optional[str]]:
    """
    Validate a PDB identifier (4 alphanumeric characters).

    Returns ``(True, None)`` if valid, ``(False, error_message)`` if not.
    """
    if not pdb_id or not pdb_id.strip():
        return False, "PDB ID must not be empty"

    pdb_id = pdb_id.strip().upper()

    if not _PDB_ID_RE.match(pdb_id):
        return False, f"Invalid PDB ID: {pdb_id!r} (must be 4 alphanumeric characters)"

    return True, None


# ---------------------------------------------------------------------------
# Rate limiting (sliding window)
# ---------------------------------------------------------------------------

_rate_lock = threading.Lock()
_rate_windows: Dict[str, Deque[float]] = defaultdict(deque)


def rate_limit_check(
    key: str,
    max_requests: int = 100,
    window_seconds: int = 60,
) -> bool:
    """
    Sliding-window rate limiter.

    Parameters
    ----------
    key : str
        Identifier to rate-limit (e.g. client IP, user ID).
    max_requests : int
        Maximum requests allowed within the window.
    window_seconds : int
        Window duration in seconds.

    Returns
    -------
    bool
        ``True`` if the request is allowed, ``False`` if rate-limited.
    """
    now = time.monotonic()
    cutoff = now - window_seconds

    with _rate_lock:
        window = _rate_windows[key]

        # Purge expired entries
        while window and window[0] < cutoff:
            window.popleft()

        if len(window) >= max_requests:
            return False

        window.append(now)
        return True


# ---------------------------------------------------------------------------
# Flask security headers
# ---------------------------------------------------------------------------


def add_security_headers(response: Any) -> Any:
    """
    Add standard security headers to a Flask ``Response`` object.

    - X-Frame-Options: DENY
    - X-Content-Type-Options: nosniff
    - X-XSS-Protection: 1; mode=block
    - Referrer-Policy: strict-origin-when-cross-origin
    - Content-Security-Policy (restrictive baseline)
    - Permissions-Policy (disable geolocation, mic, camera)

    Can also be used as a Flask ``after_request`` handler::

        app.after_request(add_security_headers)
    """
    headers = response.headers
    headers["X-Frame-Options"] = "DENY"
    headers["X-Content-Type-Options"] = "nosniff"
    headers["X-XSS-Protection"] = "1; mode=block"
    headers["Referrer-Policy"] = "strict-origin-when-cross-origin"
    headers["Content-Security-Policy"] = (
        "default-src 'self'; "
        "script-src 'self' 'unsafe-inline' 'unsafe-eval' https://cdn.jsdelivr.net; "
        "style-src 'self' 'unsafe-inline' https://cdn.jsdelivr.net; "
        "font-src 'self' https://cdn.jsdelivr.net; "
        "img-src 'self' data:; "
        "connect-src 'self'"
    )
    headers["Permissions-Policy"] = "geolocation=(), microphone=(), camera=()"
    return response


# ---------------------------------------------------------------------------
# Filename sanitization
# ---------------------------------------------------------------------------


def sanitize_filename(
    filename: str,
    max_length: int = 255,
    replacement: str = "_",
) -> str:
    """
    Sanitize a filename for safe filesystem use.

    - Strips path components (prevents directory traversal)
    - Replaces unsafe characters with ``replacement``
    - Truncates to ``max_length``
    - Raises ``ValueError`` if the result is empty
    """
    if not filename:
        raise ValueError("Filename must not be empty")

    # Remove directory components
    filename = Path(filename).name

    # Replace unsafe characters
    filename = _SAFE_FILENAME_RE.sub(replacement, filename)

    # Collapse multiple replacements
    if replacement:
        filename = re.sub(re.escape(replacement) + "+", replacement, filename)
        filename = filename.strip(replacement)

    # Truncate
    filename = filename[:max_length]

    if not filename:
        raise ValueError("Filename contains only unsafe characters")

    return filename
