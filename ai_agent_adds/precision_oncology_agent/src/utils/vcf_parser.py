"""
VCF (Variant Call Format) file parsing utilities.

Provides functions for reading VCF files, extracting gene and consequence
annotations from INFO fields, filtering variants, and generating summary
statistics.  Supports standard VCF 4.x format as produced by GATK,
DeepVariant, Strelka, Mutect2, and Parabricks pipelines.

Author: Adam Jones
Date: February 2026
"""

import logging
import re
from typing import Dict, List

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# INFO field annotation patterns
# ---------------------------------------------------------------------------

# SnpEff ANN= format: ANN=Allele|Annotation|Impact|Gene|...
# Standard SnpEff has gene at position 4 (0-indexed field 3)
_ANN_GENE_PATTERN = re.compile(
    r"ANN=[^|]*\|[^|]*\|[^|]*\|([^|]+)\|", re.IGNORECASE
)
# Simplified ANN format: ANN=Allele|Gene|Consequence|Impact (gene at position 2)
_ANN_GENE_SIMPLE_PATTERN = re.compile(
    r"ANN=[^|]*\|([A-Z][A-Z0-9_.-]*)\|", re.IGNORECASE
)
_ANN_CONSEQUENCE_PATTERN = re.compile(
    r"ANN=[^|]*\|([^|]+)\|", re.IGNORECASE
)

# VEP CSQ= format: CSQ=Allele|Consequence|...|SYMBOL|...
# Position of SYMBOL depends on the VEP fields header; we look for
# a SYMBOL= key-value or fall back to positional extraction.
_CSQ_GENE_PATTERN = re.compile(
    r"CSQ=[^;]*?(?:\|){3}([A-Za-z0-9_.-]+)", re.IGNORECASE
)
_CSQ_CONSEQUENCE_PATTERN = re.compile(
    r"CSQ=[^|]*\|([^|]+)\|", re.IGNORECASE
)

# Simple GENE= or GENEINFO= tag
_GENE_TAG_PATTERN = re.compile(
    r"(?:GENE|GENEINFO)=([A-Za-z0-9_.-]+)", re.IGNORECASE
)

# Simple consequence / effect tags
_EFFECT_TAG_PATTERN = re.compile(
    r"(?:EFFECT|CONSEQUENCE|IMPACT)=([^;]+)", re.IGNORECASE
)

# Variant type classification
_INDEL_PATTERN = re.compile(r"^[ACGTNacgtn]{2,}$")


# ===================================================================
# Core Parsers
# ===================================================================

def parse_vcf_file(filepath: str) -> List[Dict]:
    """
    Read a VCF file and parse data lines into variant dictionaries.

    Skips meta-information lines (starting with ``##``) and uses the
    header line (starting with ``#CHROM``) to identify columns.

    Parameters
    ----------
    filepath : str
        Path to the VCF file.

    Returns
    -------
    list of dict
        Parsed variant records with keys: chrom, pos, id, ref, alt,
        qual, filter, info (and any sample columns if present).

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    """
    logger.info("Parsing VCF file: %s", filepath)

    with open(filepath, "r", encoding="utf-8") as f:
        text = f.read()

    variants = parse_vcf_text(text)
    logger.info("Parsed %d variants from %s", len(variants), filepath)
    return variants


def parse_vcf_text(text: str) -> List[Dict]:
    """
    Parse VCF-formatted text into variant dictionaries.

    Parameters
    ----------
    text : str
        VCF content as a string.

    Returns
    -------
    list of dict
        Parsed variant records with keys: chrom, pos, id, ref, alt,
        qual, filter, info. Additional sample columns are included
        under their header names.
    """
    variants: List[Dict] = []
    header_columns: List[str] = []

    for line in text.splitlines():
        line = line.strip()

        # Skip empty lines
        if not line:
            continue

        # Skip meta-information lines
        if line.startswith("##"):
            continue

        # Parse header line
        if line.startswith("#CHROM") or line.startswith("#chrom"):
            header_columns = line.lstrip("#").split("\t")
            continue

        # Parse data lines
        fields = line.split("\t")

        # Use standard VCF columns if no header was found
        if not header_columns:
            header_columns = [
                "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"
            ]

        variant: Dict = {}
        for idx, col_name in enumerate(header_columns):
            if idx < len(fields):
                variant[col_name.lower()] = fields[idx]

        # Normalise common field names
        if "chrom" not in variant and "#chrom" in variant:
            variant["chrom"] = variant.pop("#chrom")

        # Convert pos to int if possible
        if "pos" in variant:
            try:
                variant["pos"] = int(variant["pos"])
            except (ValueError, TypeError):
                pass

        # Convert qual to float if possible
        if "qual" in variant and variant["qual"] != ".":
            try:
                variant["qual"] = float(variant["qual"])
            except (ValueError, TypeError):
                pass

        variants.append(variant)

    return variants


# ===================================================================
# INFO Field Extraction
# ===================================================================

def extract_gene_from_info(info: str) -> str:
    """
    Extract gene symbol from a VCF INFO field.

    Searches for gene annotations in the following order:
    1. SnpEff ANN= annotation
    2. VEP CSQ= annotation
    3. GENE= or GENEINFO= tag

    Parameters
    ----------
    info : str
        The INFO field string from a VCF record.

    Returns
    -------
    str
        Gene symbol if found, empty string otherwise.
    """
    if not info:
        return ""

    # Try ANN= (SnpEff standard format — gene at position 4)
    match = _ANN_GENE_PATTERN.search(info)
    if match:
        return match.group(1).strip()

    # Try ANN= (simplified format — gene at position 2)
    match = _ANN_GENE_SIMPLE_PATTERN.search(info)
    if match:
        return match.group(1).strip()

    # Try CSQ= (VEP)
    match = _CSQ_GENE_PATTERN.search(info)
    if match:
        return match.group(1).strip()

    # Try GENE= / GENEINFO= tag
    match = _GENE_TAG_PATTERN.search(info)
    if match:
        gene = match.group(1).strip()
        # GENEINFO may have format GENE:ID — take just the gene symbol
        if ":" in gene:
            gene = gene.split(":")[0]
        return gene

    return ""


def extract_consequence_from_info(info: str) -> str:
    """
    Extract variant consequence / effect from a VCF INFO field.

    Searches for consequence annotations in the following order:
    1. SnpEff ANN= annotation (second pipe-delimited field)
    2. VEP CSQ= annotation (second pipe-delimited field)
    3. EFFECT= or CONSEQUENCE= or IMPACT= tag

    Parameters
    ----------
    info : str
        The INFO field string from a VCF record.

    Returns
    -------
    str
        Consequence description if found, empty string otherwise.
    """
    if not info:
        return ""

    # Try ANN= (SnpEff) — consequence is in the second pipe field
    match = _ANN_CONSEQUENCE_PATTERN.search(info)
    if match:
        return match.group(1).strip()

    # Try CSQ= (VEP)
    match = _CSQ_CONSEQUENCE_PATTERN.search(info)
    if match:
        return match.group(1).strip()

    # Try explicit tags
    match = _EFFECT_TAG_PATTERN.search(info)
    if match:
        return match.group(1).strip()

    return ""


# ===================================================================
# Filtering
# ===================================================================

def filter_pass_variants(variants: List[Dict]) -> List[Dict]:
    """
    Filter variants to keep only those with PASS in the FILTER field.

    Also includes variants with an empty filter or "." (treated as PASS
    by many callers).

    Parameters
    ----------
    variants : list of dict
        Parsed variant records.

    Returns
    -------
    list of dict
        Variants where FILTER is PASS, ".", or empty.
    """
    pass_variants = []
    for v in variants:
        filt = str(v.get("filter", "")).strip().upper()
        if filt in ("PASS", ".", ""):
            pass_variants.append(v)

    logger.info(
        "Filtered to %d PASS variants from %d total",
        len(pass_variants), len(variants),
    )
    return pass_variants


# ===================================================================
# Summary Statistics
# ===================================================================

def summarize_variants(variants: List[Dict]) -> Dict:
    """
    Generate summary statistics for a list of parsed variants.

    Counts total variants, SNVs, indels (insertions + deletions),
    and unique genes affected.

    Parameters
    ----------
    variants : list of dict
        Parsed variant records with at least ref, alt, and info fields.

    Returns
    -------
    dict
        Summary with keys:
            - total: Total variant count
            - snvs: Single nucleotide variant count
            - insertions: Insertion count
            - deletions: Deletion count
            - indels: Total indel count (insertions + deletions)
            - genes_affected: List of unique gene symbols found
            - genes_count: Number of unique genes
    """
    total = len(variants)
    snvs = 0
    insertions = 0
    deletions = 0
    genes: set = set()

    for v in variants:
        ref = str(v.get("ref", ""))
        alt = str(v.get("alt", ""))
        info = str(v.get("info", ""))

        # Handle multi-allelic sites (comma-separated ALT)
        alt_alleles = alt.split(",")

        for allele in alt_alleles:
            allele = allele.strip()
            if len(ref) == 1 and len(allele) == 1:
                snvs += 1
            elif len(ref) < len(allele):
                insertions += 1
            elif len(ref) > len(allele):
                deletions += 1
            else:
                # MNV or complex — count as SNV-like
                snvs += 1

        # Extract gene
        gene = extract_gene_from_info(info)
        if gene:
            genes.add(gene)

    genes_sorted = sorted(genes)

    summary = {
        "total": total,
        "snvs": snvs,
        "insertions": insertions,
        "deletions": deletions,
        "indels": insertions + deletions,
        "genes_affected": genes_sorted,
        "genes_count": len(genes_sorted),
    }

    logger.info(
        "Variant summary: %d total, %d SNVs, %d indels, %d genes",
        total, snvs, insertions + deletions, len(genes_sorted),
    )
    return summary
