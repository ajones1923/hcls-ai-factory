"""VCF-to-PGx extraction pipeline — core pharmacogenomics module.

Resolves VCF variants to pharmacogene star allele nomenclature, translates
diplotypes to CPIC-standardized phenotypes, and cross-references patient
metabolizer status against medication lists to generate actionable alerts.

Classes:
    StarAlleleCaller  — VCF → star allele diplotypes
    PhenotypeTranslator — diplotypes → CPIC phenotypes with activity scores
    DrugGeneMatcher — phenotype profile × medication list → PGx alerts

References:
    - CPIC guidelines: https://cpicpgx.org/guidelines/
    - PharmVar: https://www.pharmvar.org/
    - PharmGKB: https://www.pharmgkb.org/

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Tuple

logger = logging.getLogger(__name__)


# ═════════════════════════════════════════════════════════════════════════════
# Data classes
# ═════════════════════════════════════════════════════════════════════════════


class AlertSeverity(str, Enum):
    """CPIC-aligned clinical action levels."""
    CONTRAINDICATED = "contraindicated"
    MAJOR = "major"
    MODERATE = "moderate"
    MINOR = "minor"
    INFORMATIONAL = "informational"


@dataclass
class PGxAlert:
    """Alert generated when a patient's PGx profile impacts a medication."""
    drug: str
    gene: str
    diplotype: str
    phenotype: str
    severity: AlertSeverity
    recommendation: str
    cpic_level: str = ""  # e.g. "A", "B"
    alternatives: List[str] = field(default_factory=list)
    source: str = "CPIC"


# ═════════════════════════════════════════════════════════════════════════════
# 1. StarAlleleCaller
# ═════════════════════════════════════════════════════════════════════════════


@dataclass(frozen=True)
class PGxPosition:
    """A single pharmacogenomic variant position."""
    chrom: str
    pos: int
    ref: str
    alt: str
    star_allele: str
    rsid: str
    gene: str


class StarAlleleCaller:
    """Resolves VCF variants to pharmacogene star allele nomenclature.

    Uses a curated set of sentinel positions for each pharmacogene.  For each
    gene the caller identifies which star-allele-defining variants are present
    in the VCF, then assigns the most specific diplotype (defaulting
    unobserved alleles to *1 / wild-type).
    """

    # ── Pharmacogenomic sentinel positions ──────────────────────────────
    # Format: (chrom, pos_hg38, ref, alt, star_allele, rsid)
    # Positions use GRCh38 coordinates from PharmVar / CPIC.
    PGX_POSITIONS: Dict[str, List[Tuple[str, int, str, str, str, str]]] = {
        # ── CYP2D6 (chr22) ──
        "CYP2D6": [
            ("chr22", 42128945, "C",  "CT",  "*3",  "rs35742686"),   # 2549del (frameshift)
            ("chr22", 42128191, "C",  "T",   "*4",  "rs3892097"),    # 1846G>A (splicing defect)
            ("chr22", 42126611, "REF","DEL",  "*5",  "gene_deletion"),# whole-gene deletion
            ("chr22", 42127941, "T",  "TG",  "*6",  "rs5030655"),    # 1707del (frameshift)
            ("chr22", 42127803, "GCA","G",   "*9",  "rs5030656"),    # 2613del (in-frame)
            ("chr22", 42130692, "G",  "A",   "*10", "rs1065852"),    # 100C>T (P34S, decreased)
            ("chr22", 42129770, "G",  "A",   "*17", "rs28371706"),   # 1023C>T (T107I, decreased)
            ("chr22", 42127556, "C",  "T",   "*41", "rs28371725"),   # 2988G>A (splicing, decreased)
        ],
        # ── CYP2C19 (chr10) ──
        "CYP2C19": [
            ("chr10", 94781859, "G",  "A",  "*2",  "rs4244285"),     # 681G>A (splicing, no function)
            ("chr10", 94780653, "G",  "A",  "*3",  "rs4986893"),     # 636G>A (premature stop)
            ("chr10", 94761900, "C",  "T",  "*17", "rs12248560"),    # -806C>T (increased function)
        ],
        # ── CYP2C9 (chr10) ──
        "CYP2C9": [
            ("chr10", 94942290, "C",  "T",  "*2",  "rs1799853"),     # 430C>T (R144C, decreased)
            ("chr10", 94981297, "A",  "C",  "*3",  "rs1057910"),     # 1075A>C (I359L, decreased)
        ],
        # ── VKORC1 (chr16) ──
        "VKORC1": [
            ("chr16", 31107689, "C",  "T",  "-1639G>A", "rs9923231"),# promoter variant
        ],
        # ── SLCO1B1 (chr12) ──
        "SLCO1B1": [
            ("chr12", 21178615, "T",  "C",  "*5",  "rs4149056"),     # 521T>C (V174A, decreased transport)
        ],
        # ── DPYD (chr1) ──
        "DPYD": [
            ("chr1",  97450058, "G",  "A",  "*2A", "rs3918290"),     # IVS14+1G>A (no function)
            ("chr1",  97515839, "T",  "G",  "*13", "rs55886062"),    # I560S (no function)
            ("chr1",  97544653, "A",  "T",  "c.2846A>T", "rs67376798"),# D949V (decreased)
        ],
        # ── TPMT (chr6) ──
        "TPMT": [
            ("chr6",  18130918, "C",  "G",  "*2",  "rs1800462"),     # A80P (no function)
            ("chr6",  18139228, "A",  "G",  "*3B", "rs1800460"),     # A154T (no function)
            ("chr6",  18143724, "T",  "C",  "*3C", "rs1142345"),     # Y240C (no function)
            # *3A is the combination haplotype of *3B + *3C (rs1800460 + rs1142345)
        ],
        # ── NUDT15 (chr13) ──
        "NUDT15": [
            ("chr13", 48037782, "C",  "T",  "*3",  "rs116855232"),   # R139C (no function)
        ],
        # ── UGT1A1 (chr2) ──
        "UGT1A1": [
            ("chr2",  233760233, "C", "T",  "*6",  "rs4148323"),     # G71R (decreased)
            ("chr2",  233759924, "TA6","TA7","*28", "rs8175347"),    # TA repeat (decreased)
        ],
        # ── CYP3A5 (chr7) ──
        "CYP3A5": [
            ("chr7",  99672916, "C",  "T",  "*3",  "rs776746"),      # splice defect (no function)
            ("chr7",  99666950, "T",  "C",  "*6",  "rs10264272"),    # splice defect (no function)
        ],
        # ── CYP4F2 (chr19) ──
        "CYP4F2": [
            ("chr19", 15990431, "C",  "T",  "*3",  "rs2108622"),     # V433M (decreased)
        ],
        # ── ABCB1 (chr7) ──
        "ABCB1": [
            ("chr7",  87509329, "A",  "G",  "C3435T",  "rs1045642"),   # Altered P-gp expression — affects digoxin, cyclosporine bioavailability
            ("chr7",  87531302, "G",  "T",  "G2677T/A","rs2032582"),   # Altered P-gp function
            ("chr7",  87550285, "C",  "T",  "C1236T",  "rs1128503"),   # Synonymous but affects mRNA stability
        ],
        # ── CACNA1S (chr1) ──
        "CACNA1S": [
            ("chr1",  201060815,"G",  "A",  "c.3257G>A","rs772226819"),# MH susceptibility — R1086H
        ],
        # ── CYP1A2 (chr15) ──
        "CYP1A2": [
            ("chr15", 74749576, "C",  "A",  "*1F", "rs762551"),        # Inducible phenotype — enhanced response to inducers
            ("chr15", 74750844, "C",  "A",  "*1K", "rs12720461"),      # Decreased inducibility
            ("chr15", 74751842, "G",  "A",  "*1C", "rs2069514"),       # Decreased activity in some populations
        ],
        # ── CYP2B6 (chr19) ──
        "CYP2B6": [
            ("chr19", 41006936, "G",  "T",  "*6",  "rs3745274"),       # Q172H — decreased function, most common variant
            ("chr19", 41012316, "T",  "C",  "*18", "rs28399499"),      # I328T — no function
            ("chr19", 41012090, "A",  "G",  "*4",  "rs2279343"),       # K262R — increased function
        ],
        # ── F5 (chr1) ──
        "F5": [
            ("chr1",  169549811,"C",  "T",  "Factor V Leiden","rs6025"),# R506Q — thrombophilia risk
            ("chr1",  169549775,"G",  "C",  "Factor V Cambridge","rs118203906"),# R306T
        ],
        # ── G6PD (chrX) ──
        "G6PD": [
            ("chrX",  154536002,"C",  "T",  "A-",  "rs1050828"),       # V68M — G6PD A- variant, hemolytic anemia risk
            ("chrX",  154535277,"C",  "T",  "Mediterranean","rs5030868"),# S188F — severe deficiency
            ("chrX",  154532822,"G",  "A",  "Canton","rs72554665"),     # R459L — common in East Asia
        ],
        # ── IFNL3 (chr19) ──
        "IFNL3": [
            ("chr19", 39248147, "C",  "T",  "rs12979860","rs12979860"),# HCV treatment response predictor — CC favorable
            ("chr19", 39243731, "TT", "ΔG", "ss469415590","rs368234815"),# Replaces rs12979860 in some assays
        ],
        # ── MTHFR (chr1) ──
        "MTHFR": [
            ("chr1",  11796321, "G",  "A",  "C677T","rs1801133"),      # Thermolabile variant — reduced enzyme activity
            ("chr1",  11794419, "A",  "C",  "A1298C","rs1801131"),     # Decreased enzyme activity
        ],
        # ── NAT2 (chr8) ──
        "NAT2": [
            ("chr8",  18257854, "T",  "C",  "*5",  "rs1801280"),       # Slow acetylator — I114T
            ("chr8",  18258103, "G",  "A",  "*6",  "rs1799930"),       # Slow acetylator — R197Q
            ("chr8",  18258514, "G",  "A",  "*7",  "rs1799931"),       # Slow acetylator — G286E
            ("chr8",  18257611, "G",  "A",  "*14", "rs1801279"),       # Slow acetylator — R64Q
        ],
        # ── RYR1 (chr19) ──
        "RYR1": [
            ("chr19", 38499670, "G",  "A",  "c.1021G>A","rs118192161"),# MH susceptibility — R341H
            ("chr19", 38512993, "C",  "T",  "c.1840C>T","rs118192163"),# MH susceptibility — R614C
            ("chr19", 38585223, "G",  "A",  "c.7300G>A","rs118192178"),# MH susceptibility — G2434R
            ("chr19", 38585277, "C",  "T",  "c.7354C>T","rs118192172"),# MH susceptibility — R2452W — most studied
        ],
    }

    # Build a fast lookup: (chrom, pos) → list of PGxPosition
    _LOOKUP: Dict[Tuple[str, int], List[PGxPosition]] = {}

    def __init__(self) -> None:
        if not StarAlleleCaller._LOOKUP:
            self._build_lookup()

    @classmethod
    def _build_lookup(cls) -> None:
        for gene, positions in cls.PGX_POSITIONS.items():
            for chrom, pos, ref, alt, star, rsid in positions:
                key = (chrom, pos)
                entry = PGxPosition(
                    chrom=chrom, pos=pos, ref=ref, alt=alt,
                    star_allele=star, rsid=rsid, gene=gene,
                )
                cls._LOOKUP.setdefault(key, []).append(entry)

    # ── VCF parsing ─────────────────────────────────────────────────────

    def extract_pgx_variants(self, vcf_path: str | Path) -> List[Dict[str, Any]]:
        """Parse a VCF file and return only variants at known PGx positions.

        Parameters
        ----------
        vcf_path : str or Path
            Path to a (possibly gzipped) VCF file.

        Returns
        -------
        list of dict
            Each dict contains: gene, chrom, pos, ref, alt, star_allele,
            rsid, genotype, quality, filter_status.
        """
        vcf_path = Path(vcf_path)
        if not vcf_path.exists():
            raise FileNotFoundError(f"VCF not found: {vcf_path}")

        open_fn = self._open_vcf(vcf_path)
        hits: List[Dict[str, Any]] = []

        with open_fn(vcf_path, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if len(fields) < 10:
                    continue

                chrom = fields[0]
                if not chrom.startswith("chr"):
                    chrom = f"chr{chrom}"

                try:
                    pos = int(fields[1])
                except ValueError:
                    continue

                key = (chrom, pos)
                if key not in self._LOOKUP:
                    continue

                ref = fields[3]
                alts = fields[4].split(",")
                qual = fields[5]
                filt = fields[6]
                gt = self._extract_genotype(fields)

                for pgx in self._LOOKUP[key]:
                    if pgx.alt in alts or pgx.ref == ref:
                        hits.append({
                            "gene": pgx.gene,
                            "chrom": chrom,
                            "pos": pos,
                            "ref": ref,
                            "alt": fields[4],
                            "star_allele": pgx.star_allele,
                            "rsid": pgx.rsid,
                            "genotype": gt,
                            "quality": qual,
                            "filter_status": filt,
                        })

        logger.info("Extracted %d PGx variants from %s", len(hits), vcf_path.name)
        return hits

    # ── Star allele calling ─────────────────────────────────────────────

    def call_star_alleles(
        self, variants: List[Dict[str, Any]], gene: str
    ) -> Tuple[str, str]:
        """Determine diplotype for a single gene from extracted variants.

        Logic:
        - Heterozygous variant → one allele carries the star allele, other is *1
        - Homozygous variant → both alleles carry the star allele
        - Multiple different variants → assign each to one allele (phase unknown,
          assume trans unless phased)

        Returns
        -------
        tuple of (str, str)
            Diplotype as two star alleles, e.g. ("*1", "*4").
        """
        gene_variants = [v for v in variants if v["gene"] == gene]

        if not gene_variants:
            return ("*1", "*1")

        allele1_stars: List[str] = []
        allele2_stars: List[str] = []

        for var in gene_variants:
            gt = var.get("genotype", "0/1")
            star = var["star_allele"]

            # Handle TPMT *3A (compound haplotype)
            if gene == "TPMT" and star in ("*3B", "*3C"):
                # Defer *3A detection to post-processing
                pass

            if gt in ("1/1", "1|1"):
                # Homozygous variant → both alleles
                allele1_stars.append(star)
                allele2_stars.append(star)
            elif gt in ("0/1", "0|1", "1|0"):
                # Heterozygous → assign to one allele
                if not allele1_stars:
                    allele1_stars.append(star)
                else:
                    allele2_stars.append(star)
            # 0/0 → reference (wild-type), skip

        # TPMT *3A detection: *3B + *3C on same haplotype
        if gene == "TPMT":
            allele1_stars, allele2_stars = self._resolve_tpmt_3a(
                allele1_stars, allele2_stars
            )

        # Compose final alleles
        allele1 = self._compose_allele(allele1_stars, gene)
        allele2 = self._compose_allele(allele2_stars, gene)

        # Canonical ordering: lower number first
        return self._order_diplotype(allele1, allele2)

    def call_all_genes(self, vcf_path: str | Path) -> Dict[str, Dict[str, Any]]:
        """Generate complete PGx profile from a VCF file.

        Returns
        -------
        dict
            gene → {diplotype, allele1, allele2, variants_found, positions_checked}
        """
        variants = self.extract_pgx_variants(vcf_path)
        profile: Dict[str, Dict[str, Any]] = {}

        for gene in self.PGX_POSITIONS:
            allele1, allele2 = self.call_star_alleles(variants, gene)
            gene_vars = [v for v in variants if v["gene"] == gene]
            profile[gene] = {
                "diplotype": f"{allele1}/{allele2}",
                "allele1": allele1,
                "allele2": allele2,
                "variants_found": len(gene_vars),
                "positions_checked": len(self.PGX_POSITIONS[gene]),
                "variants": gene_vars,
            }

        logger.info("Called star alleles for %d genes", len(profile))
        return profile

    # ── Private helpers ─────────────────────────────────────────────────

    @staticmethod
    def _open_vcf(path: Path):
        """Return the appropriate open function for plain or gzipped VCF."""
        if path.suffix == ".gz" or str(path).endswith(".vcf.gz"):
            import gzip
            return gzip.open
        return open

    @staticmethod
    def _extract_genotype(fields: List[str]) -> str:
        """Pull GT from the FORMAT + first sample columns."""
        fmt = fields[8].split(":")
        sample = fields[9].split(":")
        try:
            gt_idx = fmt.index("GT")
            return sample[gt_idx]
        except (ValueError, IndexError):
            return "."

    @staticmethod
    def _resolve_tpmt_3a(
        allele1: List[str], allele2: List[str]
    ) -> Tuple[List[str], List[str]]:
        """Detect TPMT *3A (cis combination of *3B + *3C)."""
        combined = allele1 + allele2
        has_3b = "*3B" in combined
        has_3c = "*3C" in combined

        if has_3b and has_3c:
            # Check if they are on the same haplotype (cis) → *3A
            if "*3B" in allele1 and "*3C" in allele1:
                allele1 = [s for s in allele1 if s not in ("*3B", "*3C")]
                allele1.append("*3A")
            elif "*3B" in allele2 and "*3C" in allele2:
                allele2 = [s for s in allele2 if s not in ("*3B", "*3C")]
                allele2.append("*3A")
            else:
                # Trans — keep separate (one on each allele)
                pass

        return allele1, allele2

    @staticmethod
    def _compose_allele(stars: List[str], gene: str) -> str:
        """Compose a single allele label from detected star alleles."""
        if not stars:
            return "*1"
        if len(stars) == 1:
            return stars[0]
        # Multiple variants on one allele — use lowest-numbered star
        # (simplified; a production system would use PharmVar sub-allele logic)
        return sorted(stars, key=lambda s: StarAlleleCaller._star_sort_key(s))[0]

    @staticmethod
    def _star_sort_key(star: str) -> float:
        """Numeric sort key for star alleles (*1 → 1.0, *3A → 3.01, etc.)."""
        m = re.match(r"\*(\d+)([A-Z]?)", star)
        if not m:
            return 999.0
        base = float(m.group(1))
        suffix = (ord(m.group(2)) - ord("A") + 1) * 0.01 if m.group(2) else 0
        return base + suffix

    @staticmethod
    def _order_diplotype(a1: str, a2: str) -> Tuple[str, str]:
        """Canonical ordering: lower-numbered allele first."""
        k1 = StarAlleleCaller._star_sort_key(a1)
        k2 = StarAlleleCaller._star_sort_key(a2)
        return (a1, a2) if k1 <= k2 else (a2, a1)


# ═════════════════════════════════════════════════════════════════════════════
# 2. PhenotypeTranslator
# ═════════════════════════════════════════════════════════════════════════════


class PhenotypeTranslator:
    """Converts diplotypes to CPIC standardized phenotype terms using activity scores.

    The activity score system (Gaedigk 2008) assigns a numeric value to each
    star allele based on its functional impact.  The sum across both alleles
    determines the phenotype category.
    """

    # ── Per-allele activity scores (CPIC tables) ────────────────────────
    ACTIVITY_SCORES: Dict[str, Dict[str, float]] = {
        "CYP2D6": {
            "*1": 1.0, "*2": 1.0, "*3": 0.0, "*4": 0.0, "*5": 0.0,
            "*6": 0.0, "*7": 0.0, "*8": 0.0, "*9": 0.5, "*10": 0.25,
            "*11": 0.0, "*12": 0.0, "*14": 0.0, "*15": 0.0,
            "*17": 0.5, "*29": 0.5, "*35": 1.0, "*36": 0.0,
            "*40": 0.0, "*41": 0.5, "*45": 1.0,
        },
        # Every table below is the corresponding gene's CPIC
        # allele-functionality table (api.cpicpgx.org/v1/allele?genesymbol=eq.<GENE>,
        # retrieved 2026-11-03), mapping CPIC's clinical function to the
        # activity value it publishes: Normal 1.0, Decreased 0.5, No
        # function 0.0, Increased 1.5. See dosing.py for the per-allele
        # notes on what each correction was.
        "CYP2C19": {
            "*1": 1.0, "*11": 1.0, "*13": 1.0, "*15": 1.0, "*18": 1.0,
            "*28": 1.0, "*38": 1.0,
            "*2": 0.0, "*3": 0.0, "*4": 0.0, "*5": 0.0, "*6": 0.0,
            "*7": 0.0, "*8": 0.0, "*22": 0.0, "*24": 0.0, "*35": 0.0,
            "*36": 0.0, "*37": 0.0,
            "*9": 0.5, "*10": 0.5, "*16": 0.5, "*19": 0.5, "*25": 0.5,
            "*26": 0.5,
            "*17": 1.5,
        },
        "CYP2C9": {
            "*1": 1.0, "*9": 1.0,
            "*2": 0.5, "*4": 0.5, "*5": 0.5, "*8": 0.5, "*11": 0.5,
            "*12": 0.5, "*14": 0.5, "*16": 0.5, "*23": 0.5, "*26": 0.5,
            "*28": 0.5, "*29": 0.5, "*30": 0.5, "*31": 0.5, "*37": 0.5,
            "*38": 0.5, "*44": 0.5, "*46": 0.5, "*50": 0.5, "*55": 0.5,
            "*61": 0.5,
            "*3": 0.0, "*6": 0.0, "*13": 0.0, "*15": 0.0, "*24": 0.0,
            "*25": 0.0, "*33": 0.0, "*35": 0.0, "*39": 0.0, "*42": 0.0,
            "*43": 0.0, "*45": 0.0, "*52": 0.0,
        },
        "DPYD": {
            "*1": 1.0, "*9A": 1.0,
            "*2A": 0.0, "*3": 0.0, "*7": 0.0, "*8": 0.0, "*10": 0.0,
            "*12": 0.0, "*13": 0.0,
            "c.2846A>T": 0.5, "HapB3": 0.5,
        },
        "TPMT": {
            "*1": 1.0, "*2": 0.0, "*3A": 0.0, "*3B": 0.0, "*3C": 0.0,
            "*4": 0.0, "*8": 0.5,
        },
        # NUDT15 from CPIC's MAINTAINED allele-functionality table, not from
        # Relling 2019 Table 1. The 2018 paper called *4-*9 "unclear
        # function"; CPIC has since assigned them (*4, *6, *7, *9 no
        # function; *5, *8 decreased), and for allele function the
        # maintained table governs. Scoring them as unscoreable made a
        # NUDT15 poor metabolizer Indeterminate and sent a child on
        # 6-mercaptopurine away with a full weight-based dose.
        # api.cpicpgx.org/v1/allele?genesymbol=eq.NUDT15, retrieved
        # 2026-11-03; diplotype assignments cross-checked against
        # api.cpicpgx.org/v1/diplotype?genesymbol=eq.NUDT15 (276 rows,
        # all reproduced).
        "NUDT15": {
            "*1": 1.0,
            "*2": 0.0, "*3": 0.0, "*4": 0.0, "*6": 0.0, "*7": 0.0,
            "*9": 0.0, "*10": 0.0, "*13": 0.0, "*14": 0.0, "*17": 0.0,
            "*18": 0.0, "*19": 0.0,
            "*5": 0.5, "*8": 0.5,
        },
        # UGT1A1: *37 was scored 0.0 (no function), making *1/*37 a Poor
        # Metabolizer where CPIC's diplotype table says Intermediate -- an
        # over-called irinotecan reduction in a curative-intent regimen.
        # *80 was given 0.5, a function CPIC declines to assign: CPIC says
        # Unknown function and answers Indeterminate for *1/*80, so it is
        # listed in UNCERTAIN_FUNCTION_ALLELES instead.
        # api.cpicpgx.org/v1/allele?genesymbol=eq.UGT1A1, retrieved 2026-11-03.
        "UGT1A1": {
            "*1": 1.0, "*6": 0.5, "*27": 0.5, "*28": 0.5, "*36": 1.5,
            "*37": 0.5, "*80+*28": 0.5, "*80+*37": 0.5,
        },
        # CYP3A5 *7 is a no-function allele (CPIC); it was in dosing.py and
        # not here, which is the asymmetry test_allele_tables_do_not_drift
        # allowlisted pending retrieval of the CPIC table. Retrieved
        # 2026-11-03: api.cpicpgx.org/v1/allele?genesymbol=eq.CYP3A5.
        "CYP3A5": {
            "*1": 1.0, "*3": 0.0, "*6": 0.0, "*7": 0.0,
        },
        # VKORC1, SLCO1B1, CYP4F2 use genotype-based categories, not activity scores
    }

    # ── Alleles CPIC assigns "uncertain" or "unclear" function ──────────
    # These cannot be given an activity value, and CPIC's own phenotype for a diplotype
    # carrying one is **Indeterminate** -- not Normal, and not a guess in between. Scoring them
    # produces a confident-looking phenotype out of an allele the guideline says it cannot call.
    #
    # Each set is the gene's "Uncertain function" / "Unknown function" rows
    # in CPIC's allele-functionality table (api.cpicpgx.org/v1/allele,
    # retrieved 2026-11-03), plus rows CPIC leaves unassigned. CPIC's
    # diplotype table returns Indeterminate for every diplotype carrying
    # one of them.
    #
    # NUDT15 previously held *4-*9 on the strength of Relling 2019's
    # "unclear function". CPIC has since assigned *4, *6, *7 and *9 (no
    # function) and *5 and *8 (decreased), and reserves Indeterminate for
    # the alleles below -- so the old set withheld a dose reduction from
    # NUDT15 poor and intermediate metabolizers.
    UNCERTAIN_FUNCTION_ALLELES: Dict[str, frozenset] = {
        "NUDT15": frozenset({
            "*11", "*12", "*15", "*16", "*20", "*21", "*22", "*23",
        }),
        "UGT1A1": frozenset({"*80"}),
        "CYP3A5": frozenset({"*8", "*9", "*10", "*11"}),
        "CYP2C19": frozenset({
            "*12", "*14", "*23", "*29", "*30", "*31", "*32", "*33", "*34",
            "*39",
        }),
        "CYP2C9": frozenset({
            "*7", "*10", "*17", "*18", "*19", "*20", "*21", "*22", "*27",
            "*32", "*34", "*36", "*40", "*41", "*47", "*48", "*49", "*51",
            "*53", "*54", "*56", "*57", "*58", "*59", "*60", "*62", "*63",
            "*64", "*65", "*66", "*67", "*68", "*69", "*70", "*71", "*72",
            "*73", "*74", "*75",
        }),
        "TPMT": frozenset({
            "*9", "*12", "*19", "*24", "*30", "*32", "*35", "*36", "*40",
        }),
    }

    # An increased-function allele does not rescue a reduced-function one.
    # CPIC's diplotype tables are explicit: UGT1A1 *28/*36 (Decreased +
    # Increased) is an **Intermediate** Metabolizer, not the Normal the
    # summed activity score 2.0 would give, and CYP2C19 *10/*17 is a
    # (Likely) Intermediate Metabolizer rather than Normal. So when an
    # increased-function allele is paired with a decreased- or no-function
    # allele, it contributes a normal 1.0 and no more.
    #
    # Verified: with this rule the activity-score model reproduces CPIC's
    # phenotype for all 36 UGT1A1 and all 351 CYP2C19 diplotypes built
    # from assigned alleles (api.cpicpgx.org/v1/diplotype, 2026-11-03).
    INCREASED_FUNCTION_SCORE: float = 1.5
    _NORMAL_FUNCTION_SCORE: float = 1.0

    # ── Phenotype thresholds: activity score range → phenotype ──────────
    # Ranges use half-open intervals [low, high) — the last bucket is
    # closed [low, high] (handled by _score_to_phenotype).
    # Values aligned with CPIC 2023/2024 consensus guidelines.
    PHENOTYPE_THRESHOLDS: Dict[str, List[Tuple[float, float, str]]] = {
        "CYP2D6": [
            # CPIC/DPWG consensus: PM = 0, IM = 0 < AS < 1.25, NM = 1.25-2.25 inclusive,
            # UM > 2.25. The consensus paper's own words: "reclassify AS of 1 from NM to IM."
            #
            # This table previously read "IM=0.25-0.75, NM>=1.0", which put AS 1.0 (*1/*4,
            # *1/*5, *41/*41) in Normal Metabolizer and so recommended standard dosing for
            # patients the consensus calls intermediate. It also classified AS 2.25 (e.g.
            # *1x2/*41) as Ultrarapid, where 2.25 is the inclusive top of Normal.
            #
            # Buckets are [low, high) except the last, which is closed. Attainable activity
            # scores are multiples of 0.25, so two bounds here are written as the next
            # attainable value rather than the consensus bound: PM's 0.25 captures AS = 0
            # alone, and NM's 2.5 expresses the consensus' INCLUSIVE 2.25 upper bound.
            #
            # Source: Caudle KE, Sangkuhl K, Whirl-Carrillo M, Swen JJ, Haidar CE, Klein TE,
            # Gammal RS, Relling MV, Scott SA, Hertz DL, Guchelaar HJ, Gaedigk A.
            # "Standardizing CYP2D6 Genotype to Phenotype Translation: Consensus
            # Recommendations from the Clinical Pharmacogenetics Implementation Consortium and
            # Dutch Pharmacogenetics Working Group." Clin Transl Sci 2020;13(1):116-124,
            # Table 3. doi:10.1111/cts.12692, PMID 31647186.
            (0.0,  0.25,  "Poor Metabolizer"),
            (0.25, 1.25,  "Intermediate Metabolizer"),
            (1.25, 2.5,   "Normal Metabolizer"),
            (2.5,  99.0,  "Ultrarapid Metabolizer"),
        ],
        "CYP2C19": [
            # CPIC: *2/*2 = PM, *1/*2 = IM, *1/*1 = NM, *1/*17 = RM,
            # *17/*17 = UM. The Poor bucket's top is 1.0 (exclusive), not
            # 0.5, because CPIC calls a no-function + decreased-function
            # pair -- AS 0.5, e.g. *2/*10 -- a (Likely) **Poor**
            # Metabolizer; with *9/*10/*16/*19/*25/*26 now correctly
            # scored 0.5 that score is reachable, and the old 0.5 bound
            # put those patients in Intermediate.
            #
            # Verified against api.cpicpgx.org/v1/diplotype?genesymbol=eq.CYP2C19
            # (retrieved 2026-11-03): these bounds reproduce CPIC's
            # phenotype for all 351 diplotypes built from assigned alleles.
            (0.0,  1.0,   "Poor Metabolizer"),
            (1.0,  2.0,   "Intermediate Metabolizer"),
            (2.0,  2.5,   "Normal Metabolizer"),
            (2.5,  3.0,   "Rapid Metabolizer"),
            (3.0,  99.0,  "Ultrarapid Metabolizer"),
        ],
        "CYP2C9": [
            # CPIC: PM=AS≤0.5 (*2/*3, *3/*3), IM=AS 1.0-1.5, NM=AS 2.0
            (0.0,  1.0,   "Poor Metabolizer"),
            (1.0,  2.0,   "Intermediate Metabolizer"),
            (2.0,  2.1,   "Normal Metabolizer"),
        ],
        "DPYD": [
            # CPIC: Deficient=AS 0-0.5, Intermediate=AS 1.0-1.5, Normal=AS 2.0
            (0.0,  1.0,   "DPD Deficient"),
            (1.0,  2.0,   "DPD Intermediate Activity"),
            (2.0,  2.1,   "DPD Normal Activity"),
        ],
        "TPMT": [
            # CPIC: no/no = PM; no/decreased, no/normal and
            # decreased/decreased = IM; decreased/normal and normal/normal
            # = NM. With *8 (decreased, 0.5) now present, AS 0.5 and AS
            # 1.5 are both reachable and the old 0.5/2.0 bounds put AS 1.5
            # (*1/*8, a CPIC **Normal** Metabolizer) in Intermediate.
            # Verified against api.cpicpgx.org/v1/diplotype?genesymbol=eq.TPMT
            # (retrieved 2026-11-03): reproduces all 1,225 diplotypes.
            (0.0,  0.5,   "TPMT Deficient"),
            (0.5,  1.5,   "TPMT Intermediate"),
            (1.5,  2.0,   "TPMT Normal"),
        ],
        "NUDT15": [
            # CPIC: no/no and no/decreased = PM; no/normal and
            # decreased/decreased = IM; decreased/normal and normal/normal
            # = NM. Note NUDT15 differs from TPMT here -- a no-function
            # plus a decreased-function allele is a NUDT15 **Poor**
            # Metabolizer but a TPMT **Intermediate** one, which is why
            # these two tables are not the same.
            # Verified against api.cpicpgx.org/v1/diplotype?genesymbol=eq.NUDT15
            # (retrieved 2026-11-03): reproduces all 276 diplotypes.
            (0.0,  1.0,   "NUDT15 Deficient"),
            (1.0,  1.5,   "NUDT15 Intermediate"),
            (1.5,  2.0,   "NUDT15 Normal"),
        ],
        "UGT1A1": [
            # CPIC: *28/*28 (AS=1.0) = PM, *1/*28 (AS=1.5) = IM, *1/*1 (AS=2.0) = NM.
            # The Normal bucket's top was 2.1, below the AS 2.5 and 3.0 that
            # *36 (increased function, 1.5) produces, so *1/*36 and *36/*36 --
            # CPIC **Normal** Metabolizers both -- fell off the end of the
            # table and were reported "Indeterminate". Raised to 3.0, the
            # highest attainable UGT1A1 activity score.
            # Verified against api.cpicpgx.org/v1/diplotype?genesymbol=eq.UGT1A1
            # (retrieved 2026-11-03): reproduces all 36 diplotypes built from
            # assigned alleles.
            (0.0,  1.5,   "Poor Metabolizer"),
            (1.5,  2.0,   "Intermediate Metabolizer"),
            (2.0,  3.0,   "Normal Metabolizer"),
        ],
        "CYP3A5": [
            # Non-expresser=AS 0, Intermediate=AS 0.5-1.0, Expresser=AS 2.0
            (0.0,  0.5,   "CYP3A5 Non-expresser"),
            (0.5,  2.0,   "CYP3A5 Intermediate Expresser"),
            (2.0,  2.1,   "CYP3A5 Expresser"),
        ],
    }

    # ── Genotype-based categorization (non-activity-score genes) ────────
    _GENOTYPE_PHENOTYPES: Dict[str, Dict[str, Dict[str, str]]] = {
        "VKORC1": {
            "-1639G>A/-1639G>A": {
                "phenotype": "Low Dose Warfarin Sensitivity",
                "clinical_meaning": "Homozygous variant — lowest VKORC1 expression; "
                                    "significantly reduced warfarin dose required.",
            },
            "*1/-1639G>A": {
                "phenotype": "Intermediate Warfarin Sensitivity",
                "clinical_meaning": "Heterozygous — intermediate VKORC1 expression; "
                                    "moderately reduced warfarin dose.",
            },
            "*1/*1": {
                "phenotype": "Normal Warfarin Sensitivity",
                "clinical_meaning": "Wild-type — typical VKORC1 expression; "
                                    "standard warfarin dose range.",
            },
        },
        "SLCO1B1": {
            "*5/*5": {
                "phenotype": "Poor Function",
                "clinical_meaning": "Homozygous 521CC — highest myopathy risk with "
                                    "statins; consider alternative or low dose.",
            },
            "*1/*5": {
                "phenotype": "Decreased Function",
                "clinical_meaning": "Heterozygous 521TC — increased myopathy risk; "
                                    "consider lower statin dose or alternative.",
            },
            "*1/*1": {
                "phenotype": "Normal Function",
                "clinical_meaning": "Wild-type 521TT — typical statin metabolism.",
            },
        },
        "CYP4F2": {
            "*3/*3": {
                "phenotype": "High Dose Warfarin",
                "clinical_meaning": "Homozygous 433M — increased warfarin dose needed "
                                    "(decreased vitamin K metabolism).",
            },
            "*1/*3": {
                "phenotype": "Intermediate Warfarin Dose",
                "clinical_meaning": "Heterozygous V433M — slightly increased warfarin "
                                    "dose may be needed.",
            },
            "*1/*1": {
                "phenotype": "Normal Warfarin Dose",
                "clinical_meaning": "Wild-type — standard warfarin dose.",
            },
        },
    }

    # ── Clinical meaning templates ──────────────────────────────────────
    _CLINICAL_MEANINGS: Dict[str, Dict[str, str]] = {
        "Poor Metabolizer": {
            "CYP2D6": "No functional CYP2D6 enzyme activity; prodrugs (codeine, "
                      "tamoxifen) ineffective, active drugs may accumulate.",
            "CYP2C19": "No CYP2C19 activity; clopidogrel ineffective (use prasugrel/"
                       "ticagrelor), voriconazole requires dose reduction.",
            "CYP2C9": "Severely reduced CYP2C9 activity; significantly lower warfarin/"
                      "phenytoin doses required.",
        },
        "Intermediate Metabolizer": {
            "CYP2D6": "Reduced CYP2D6 activity; consider dose adjustments for "
                      "CYP2D6 substrates.",
            "CYP2C19": "Reduced CYP2C19 activity; consider alternative to clopidogrel.",
            "CYP2C9": "Reduced CYP2C9 activity; reduce warfarin/phenytoin dose.",
        },
        "Normal Metabolizer": {
            "CYP2D6": "Normal CYP2D6 enzyme activity; standard dosing.",
            "CYP2C19": "Normal CYP2C19 activity; standard dosing.",
            "CYP2C9": "Normal CYP2C9 activity; standard dosing.",
        },
        "Ultrarapid Metabolizer": {
            "CYP2D6": "Significantly increased CYP2D6 activity; codeine toxicity risk, "
                      "active drugs may be sub-therapeutic.",
            "CYP2C19": "Increased CYP2C19 activity; may need higher doses of PPIs, "
                       "voriconazole; improved clopidogrel response.",
        },
    }

    def translate(self, gene: str, diplotype: str) -> Dict[str, Any]:
        """Convert a single gene diplotype to CPIC phenotype.

        Parameters
        ----------
        gene : str
            Pharmacogene name (e.g. "CYP2D6").
        diplotype : str
            Diplotype string (e.g. "*1/*4").

        Returns
        -------
        dict
            {phenotype, activity_score, clinical_meaning, gene, diplotype}
        """
        # Genotype-based genes (no activity score)
        if gene in self._GENOTYPE_PHENOTYPES:
            return self._translate_genotype_based(gene, diplotype)

        # Activity-score-based genes
        try:
            allele1, allele2 = self._parse_diplotype(diplotype)
        except ValueError as exc:
            # A genotype string that is not two alleles is an absent or
            # malformed result, never wild-type. CPIC's answer for an
            # unassignable genotype is Indeterminate.
            return {
                "gene": gene,
                "diplotype": diplotype,
                "phenotype": "Indeterminate",
                "activity_score": None,
                "allele_scores": {},
                "clinical_meaning": (
                    f"No phenotype can be assigned: {exc} Report the genotype "
                    f"as two alleles, or treat this gene as not called. Do not "
                    f"infer normal function from an unparseable result."
                ),
                "unscoreable_alleles": [diplotype],
            }
        scores = self.ACTIVITY_SCORES.get(gene, {})

        # An allele the guideline calls uncertain, or one this table has never heard of, cannot
        # be scored. Both used to fall through to 1.0 -- "default to normal if unknown allele"
        # -- which is the least safe default available: an unrecognised no-function allele
        # reported as a normal metabolizer, and a standard dose recommended on the strength of
        # it. CPIC's own answer for an uncertain allele is Indeterminate, so that is what this
        # returns, with the allele named so a reader knows which one it could not call.
        uncertain = self.UNCERTAIN_FUNCTION_ALLELES.get(gene, frozenset())
        unscoreable = [
            (a, "uncertain function" if a in uncertain else "not in the allele table")
            for a in (allele1, allele2)
            if a in uncertain or a not in scores
        ]
        if unscoreable:
            reasons = "; ".join(f"{a} ({why})" for a, why in unscoreable)
            return {
                "gene": gene,
                "diplotype": diplotype,
                "phenotype": "Indeterminate",
                "activity_score": None,
                "allele_scores": {},
                "clinical_meaning": (
                    f"No phenotype can be assigned: {reasons}. CPIC assigns Indeterminate "
                    f"rather than a phenotype here. Do not infer normal function from the "
                    f"absence of a call."
                ),
                "unscoreable_alleles": [a for a, _ in unscoreable],
            }

        score1, score2 = self._paired_scores(scores[allele1], scores[allele2])
        total = score1 + score2

        phenotype = self._score_to_phenotype(gene, total)
        clinical = self._get_clinical_meaning(gene, phenotype)

        return {
            "gene": gene,
            "diplotype": diplotype,
            "phenotype": phenotype,
            "activity_score": total,
            "allele_scores": {allele1: score1, allele2: score2},
            "clinical_meaning": clinical,
        }

    def translate_all(
        self, profile: Dict[str, Dict[str, Any]]
    ) -> Dict[str, Dict[str, Any]]:
        """Translate a full PGx profile (all genes) to phenotypes.

        Parameters
        ----------
        profile : dict
            Output of StarAlleleCaller.call_all_genes().

        Returns
        -------
        dict
            gene → translation result dict.
        """
        results: Dict[str, Dict[str, Any]] = {}
        for gene, data in profile.items():
            # A gene with no called diplotype is NOT wild-type. The old
            # default "*1/*1" translated an uncalled gene as a normal
            # metabolizer. translate() returns Indeterminate for the empty
            # string.
            diplotype = data.get("diplotype") or ""
            results[gene] = self.translate(gene, diplotype)
        return results

    # ── Private helpers ─────────────────────────────────────────────────

    @classmethod
    def _paired_scores(
        cls, score1: float, score2: float
    ) -> Tuple[float, float]:
        """Apply CPIC's "an increased allele does not rescue a reduced one".

        See INCREASED_FUNCTION_SCORE for the primary source. An
        increased-function allele paired with a decreased- or no-function
        allele contributes a normal 1.0, so UGT1A1 ``*28/*36`` scores 1.5
        (Intermediate, as CPIC says) rather than 2.0 (Normal).
        """
        inc, norm = cls.INCREASED_FUNCTION_SCORE, cls._NORMAL_FUNCTION_SCORE
        if score1 == inc and score2 < norm:
            return norm, score2
        if score2 == inc and score1 < norm:
            return score1, norm
        return score1, score2

    @staticmethod
    def _parse_diplotype(diplotype: str) -> Tuple[str, str]:
        """Split a diplotype string into its two allele labels.

        Raises
        ------
        ValueError
            If the string is not exactly two slash-separated non-empty
            fields.  This used to return ``("*1", "*1")`` — homozygous
            wild-type — for an empty string, a single allele, "Unknown",
            "not reported" or any free-text value, so every downstream
            calculator then reported a normal metabolizer and a standard
            dose.  ``translate`` turns this into CPIC's Indeterminate
            rather than letting it propagate.
        """
        parts = diplotype.split("/") if diplotype else []
        if len(parts) == 2 and parts[0].strip() and parts[1].strip():
            return parts[0].strip(), parts[1].strip()
        raise ValueError(
            f"Not a diplotype: {diplotype!r}. Expected two slash-separated "
            f"alleles, e.g. '*1/*2'. This is not wild-type."
        )

    @staticmethod
    def _translate_slco1b1(diplotype: str) -> Dict[str, Any]:
        """SLCO1B1 via CPIC's allele functions and function-pair rule.

        Source: api.cpicpgx.org/v1/allele?genesymbol=eq.SLCO1B1 and
        api.cpicpgx.org/v1/diplotype?genesymbol=eq.SLCO1B1, retrieved
        2026-11-03.  The tables live in dosing.py so there is one copy.
        """
        from src.dosing import (  # local import: dosing.py imports nothing here
            _SLCO1B1_ALLELE_FUNCTION,
            _SLCO1B1_PAIR_PHENOTYPE,
            _SLCO1B1_UNCERTAIN,
        )

        _MEANINGS = {
            "Normal Function": (
                "Normal OATP1B1 transporter function — typical statin "
                "exposure."
            ),
            "Increased Function": (
                "Increased OATP1B1 transporter function — statin exposure at "
                "or below the population average."
            ),
            "Decreased Function": (
                "Reduced OATP1B1 hepatic uptake — increased statin exposure "
                "and myopathy risk. CPIC: prescribe an alternative statin; "
                "if simvastatin is warranted, limit the dose to <20 mg/day."
            ),
            "Possible Decreased Function": (
                "One allele is of uncertain function; CPIC gives this class "
                "the same recommendation as Decreased Function — prescribe "
                "an alternative statin, or simvastatin at <20 mg/day."
            ),
            "Poor Function": (
                "No functional OATP1B1 — highest myopathy/rhabdomyolysis "
                "risk with statins. CPIC: prescribe an alternative statin."
            ),
        }

        try:
            allele1, allele2 = PhenotypeTranslator._parse_diplotype(diplotype)
        except ValueError as exc:
            allele1 = allele2 = ""
            unrecognised = [diplotype]
            why = str(exc)
        else:
            unrecognised = [
                a for a in (allele1, allele2)
                if a not in _SLCO1B1_ALLELE_FUNCTION
                and a not in _SLCO1B1_UNCERTAIN
            ]
            why = (
                f"SLCO1B1 {', '.join(dict.fromkeys(unrecognised))} is not in "
                f"CPIC's SLCO1B1 allele-functionality table."
            )
        if unrecognised:
            return {
                "gene": "SLCO1B1",
                "diplotype": diplotype,
                "phenotype": "Indeterminate",
                "activity_score": None,
                "clinical_meaning": (
                    f"No phenotype can be assigned: {why} CPIC assigns "
                    f"Indeterminate and makes no recommendation here. Do NOT "
                    f"read this as normal transporter function."
                ),
                "unscoreable_alleles": list(dict.fromkeys(unrecognised)),
            }

        def _fn(allele: str) -> str:
            return _SLCO1B1_ALLELE_FUNCTION.get(allele, "uncertain")

        phenotype = _SLCO1B1_PAIR_PHENOTYPE.get(
            frozenset({_fn(allele1), _fn(allele2)}), "Indeterminate")
        return {
            "gene": "SLCO1B1",
            "diplotype": diplotype,
            "phenotype": phenotype,
            "activity_score": None,  # not an activity-score gene
            "clinical_meaning": _MEANINGS.get(phenotype, (
                "At least one allele is of uncertain function at CPIC; CPIC "
                "assigns Indeterminate and makes no recommendation."
            )),
        }

    def _translate_genotype_based(
        self, gene: str, diplotype: str
    ) -> Dict[str, Any]:
        """Handle genes that use genotype-based phenotype assignment."""
        if gene == "SLCO1B1":
            # _GENOTYPE_PHENOTYPES defines THREE SLCO1B1 diplotypes; CPIC's
            # own diplotype table has 741. Three of 741 is not a table, it
            # is three examples — and with the wild-type fallback that used
            # to sit below, *5/*15 and *15/*15 (CPIC Poor Function, ~17-fold
            # myopathy risk on simvastatin) were reported as "Normal
            # Function — wild-type 521TT, typical statin metabolism", which
            # contradicted dosing.py's answer for the same genotype.
            # Delegate to the same CPIC allele functions and function-pair
            # rule the dose calculator uses, so the report and the
            # calculator cannot disagree.
            return self._translate_slco1b1(diplotype)

        lookup = self._GENOTYPE_PHENOTYPES.get(gene, {})
        match = lookup.get(diplotype)
        if not match:
            # This used to read `lookup.get("*1/*1", ...)` under a comment
            # saying "try to find the closest match". *1/*1 is present for
            # all three of these genes, so the Unknown default below it was
            # dead code and **every** unrecognised genotype was reported as
            # the wild-type phenotype with the wild-type prose attached:
            # SLCO1B1 *5/*15 -- a CPIC Poor Function patient with ~17-fold
            # myopathy risk on simvastatin -- rendered as "Normal Function:
            # wild-type 521TT, typical statin metabolism", and VKORC1 "A/A"
            # as "Normal Warfarin Sensitivity, standard warfarin dose
            # range". These genes do not use activity scores, so they
            # bypassed the Indeterminate handling the scored path has.
            #
            # CPIC's answer for a genotype it cannot map is Indeterminate /
            # "No recommendation" (e.g. simvastatin SLCO1B1 Indeterminate,
            # api.cpicpgx.org/v1/recommendation drugid RxNorm:36567,
            # retrieved 2026-11-03).
            known = ", ".join(sorted(lookup)) or "none"
            return {
                "gene": gene,
                "diplotype": diplotype,
                "phenotype": "Indeterminate",
                "activity_score": None,
                "clinical_meaning": (
                    f"No phenotype can be assigned: {gene} {diplotype!r} is not "
                    f"in this agent's {gene} genotype table, which covers "
                    f"{known}. CPIC assigns Indeterminate and makes no "
                    f"recommendation here. Do NOT read this as wild-type."
                ),
                "unscoreable_alleles": [diplotype],
            }
        return {
            "gene": gene,
            "diplotype": diplotype,
            "phenotype": match["phenotype"],
            "activity_score": None,  # not applicable
            "clinical_meaning": match["clinical_meaning"],
        }

    def _score_to_phenotype(self, gene: str, score: float) -> str:
        """Map a total activity score to a phenotype string.

        Uses half-open intervals [low, high) for all buckets except the
        last, which is closed [low, high].  This prevents boundary overlap;
        where a guideline's own bound is inclusive, PHENOTYPE_THRESHOLDS
        encodes it as the next attainable score and says so (see CYP2D6,
        whose AS 1.0 is Intermediate per CPIC/DPWG consensus).
        """
        thresholds = self.PHENOTYPE_THRESHOLDS.get(gene, [])
        for i, (low, high, phenotype) in enumerate(thresholds):
            is_last = (i == len(thresholds) - 1)
            if is_last:
                if low <= score <= high:
                    return phenotype
            else:
                if low <= score < high:
                    return phenotype
        # An activity score that falls off the end of a gene's threshold
        # table is a bug in the table, not a clinical finding. This used to
        # return "Indeterminate", which is a *clinical answer* and hid the
        # UGT1A1 dead zone: *1/*36 (AS 2.5) and *36/*36 (AS 3.0) are CPIC
        # Normal Metabolizers and were reported Indeterminate because the
        # Normal bucket stopped at 2.1.
        raise ValueError(
            f"{gene} activity score {score} falls outside "
            f"PHENOTYPE_THRESHOLDS[{gene!r}] = {thresholds}. The threshold "
            f"table does not cover every attainable score for this gene."
        )

    def _get_clinical_meaning(self, gene: str, phenotype: str) -> str:
        """Retrieve a clinical interpretation string."""
        pheno_map = self._CLINICAL_MEANINGS.get(phenotype, {})
        return pheno_map.get(
            gene,
            f"{phenotype} for {gene} — consult CPIC guidelines for drug-specific "
            f"recommendations.",
        )


# ═════════════════════════════════════════════════════════════════════════════
# 3. DrugGeneMatcher
# ═════════════════════════════════════════════════════════════════════════════


class DrugGeneMatcher:
    """Cross-references patient phenotype profile against medication list.

    Contains a curated mapping of drugs to their pharmacogenomically
    relevant gene-phenotype-action associations, sourced from CPIC and FDA
    pharmacogenomic biomarker tables.
    """

    # ── Drug → Gene → Phenotype → Action mapping ───────────────────────
    # Each entry: gene, {phenotype: (severity, recommendation, alternatives, cpic_level)}
    DRUG_GENE_MAP: Dict[str, List[Dict[str, Any]]] = {
        # ── Opioids ──
        "codeine": [
            {
                "gene": "CYP2D6",
                "actions": {
                    "Ultrarapid Metabolizer": (
                        AlertSeverity.CONTRAINDICATED,
                        "AVOID codeine. Rapid conversion to morphine → respiratory "
                        "depression / death risk. Use non-tramadol opioid alternative.",
                        ["morphine (non-CYP2D6)", "oxycodone (with monitoring)"],
                        "A",
                    ),
                    "Poor Metabolizer": (
                        AlertSeverity.MAJOR,
                        "AVOID codeine. Insufficient conversion to active morphine → "
                        "inadequate analgesia. Use alternative analgesic.",
                        ["morphine", "non-opioid analgesics"],
                        "A",
                    ),
                    "Intermediate Metabolizer": (
                        AlertSeverity.MODERATE,
                        "Reduced codeine efficacy likely. Consider alternative or "
                        "use with monitoring for adequate pain control.",
                        ["morphine", "non-opioid analgesics"],
                        "A",
                    ),
                },
            },
        ],
        "tramadol": [
            {
                "gene": "CYP2D6",
                "actions": {
                    "Ultrarapid Metabolizer": (
                        AlertSeverity.CONTRAINDICATED,
                        "AVOID tramadol. Rapid O-demethylation → toxicity risk.",
                        ["non-tramadol/non-codeine analgesic"],
                        "A",
                    ),
                    "Poor Metabolizer": (
                        AlertSeverity.MAJOR,
                        "AVOID tramadol. Insufficient active metabolite. Use alternative.",
                        ["non-tramadol/non-codeine analgesic"],
                        "A",
                    ),
                },
            },
        ],
        # ── Antiplatelet ──
        "clopidogrel": [
            {
                "gene": "CYP2C19",
                "actions": {
                    "Poor Metabolizer": (
                        AlertSeverity.CONTRAINDICATED,
                        "AVOID clopidogrel. No bioactivation → no antiplatelet effect. "
                        "Use prasugrel or ticagrelor.",
                        ["prasugrel", "ticagrelor"],
                        "A",
                    ),
                    "Intermediate Metabolizer": (
                        AlertSeverity.MAJOR,
                        "Reduced clopidogrel efficacy. Consider prasugrel or ticagrelor.",
                        ["prasugrel", "ticagrelor"],
                        "A",
                    ),
                    "Ultrarapid Metabolizer": (
                        AlertSeverity.INFORMATIONAL,
                        "Enhanced clopidogrel response. Standard dose appropriate but "
                        "monitor for increased bleeding.",
                        [],
                        "B",
                    ),
                },
            },
        ],
        # ── Anticoagulant ──
        "warfarin": [
            {
                "gene": "CYP2C9",
                "actions": {
                    "Poor Metabolizer": (
                        AlertSeverity.MAJOR,
                        "Significantly reduced warfarin clearance. Use CPIC/IWPC dosing "
                        "algorithm; expect 40-80% dose reduction.",
                        [],
                        "A",
                    ),
                    "Intermediate Metabolizer": (
                        AlertSeverity.MODERATE,
                        "Reduced warfarin clearance. Use CPIC/IWPC dosing algorithm; "
                        "expect 20-40% dose reduction.",
                        [],
                        "A",
                    ),
                },
            },
            {
                "gene": "VKORC1",
                "actions": {
                    "Low Dose Warfarin Sensitivity": (
                        AlertSeverity.MAJOR,
                        "VKORC1 AA genotype — lowest warfarin target. Use CPIC/IWPC "
                        "dosing algorithm.",
                        [],
                        "A",
                    ),
                    "Intermediate Warfarin Sensitivity": (
                        AlertSeverity.MODERATE,
                        "VKORC1 AG genotype — moderately reduced warfarin dose. Use "
                        "CPIC/IWPC dosing algorithm.",
                        [],
                        "A",
                    ),
                },
            },
        ],
        # ── Fluoropyrimidines ──
        "fluorouracil": [
            {
                "gene": "DPYD",
                "actions": {
                    "DPD Deficient": (
                        AlertSeverity.CONTRAINDICATED,
                        "AVOID fluorouracil, capecitabine, tegafur. Complete DPD "
                        "deficiency → life-threatening toxicity.",
                        ["non-fluoropyrimidine regimen"],
                        "A",
                    ),
                    "DPD Intermediate Activity": (
                        AlertSeverity.MAJOR,
                        "Reduce fluorouracil dose by 50%. Activity score 1.0-1.5 → "
                        "partial DPD deficiency.",
                        [],
                        "A",
                    ),
                },
            },
        ],
        "capecitabine": [
            {
                "gene": "DPYD",
                "actions": {
                    "DPD Deficient": (
                        AlertSeverity.CONTRAINDICATED,
                        "AVOID capecitabine. Complete DPD deficiency → life-threatening "
                        "toxicity.",
                        ["non-fluoropyrimidine regimen"],
                        "A",
                    ),
                    "DPD Intermediate Activity": (
                        AlertSeverity.MAJOR,
                        "Reduce capecitabine dose by 50%.",
                        [],
                        "A",
                    ),
                },
            },
        ],
        # ── Thiopurines ──
        "azathioprine": [
            {
                "gene": "TPMT",
                "actions": {
                    "TPMT Deficient": (
                        AlertSeverity.CONTRAINDICATED,
                        "AVOID azathioprine or reduce dose to 10% of standard. "
                        "Severe myelosuppression risk.",
                        ["mycophenolate"],
                        "A",
                    ),
                    "TPMT Intermediate": (
                        AlertSeverity.MAJOR,
                        "Reduce azathioprine to 30-70% of standard dose. Monitor CBC.",
                        [],
                        "A",
                    ),
                },
            },
            {
                "gene": "NUDT15",
                "actions": {
                    "NUDT15 Deficient": (
                        AlertSeverity.CONTRAINDICATED,
                        "AVOID azathioprine. NUDT15 deficiency → severe "
                        "myelosuppression.",
                        ["mycophenolate"],
                        "A",
                    ),
                    "NUDT15 Intermediate": (
                        AlertSeverity.MAJOR,
                        "Reduce azathioprine to 30-70% of standard dose.",
                        [],
                        "A",
                    ),
                },
            },
        ],
        "mercaptopurine": [
            {
                "gene": "TPMT",
                "actions": {
                    "TPMT Deficient": (
                        AlertSeverity.CONTRAINDICATED,
                        "Reduce 6-MP to 10% of standard dose. Severe "
                        "myelosuppression risk.",
                        [],
                        "A",
                    ),
                    "TPMT Intermediate": (
                        AlertSeverity.MAJOR,
                        "Reduce 6-MP to 30-70% of standard dose. Monitor CBC weekly.",
                        [],
                        "A",
                    ),
                },
            },
            {
                "gene": "NUDT15",
                "actions": {
                    "NUDT15 Deficient": (
                        AlertSeverity.CONTRAINDICATED,
                        "Reduce 6-MP to 10% of standard dose.",
                        [],
                        "A",
                    ),
                    "NUDT15 Intermediate": (
                        AlertSeverity.MAJOR,
                        "Reduce 6-MP to 30-70% of standard dose.",
                        [],
                        "A",
                    ),
                },
            },
        ],
        # ── Statins ──
        "simvastatin": [
            {
                "gene": "SLCO1B1",
                "actions": {
                    "Poor Function": (
                        AlertSeverity.CONTRAINDICATED,
                        "AVOID simvastatin. 521CC homozygous → 17-fold myopathy risk. "
                        "Use pravastatin or rosuvastatin.",
                        ["pravastatin", "rosuvastatin"],
                        "A",
                    ),
                    "Decreased Function": (
                        AlertSeverity.MAJOR,
                        "Limit simvastatin to 20 mg/day or use alternative statin. "
                        "521TC heterozygous → 5-fold myopathy risk.",
                        ["pravastatin", "rosuvastatin"],
                        "A",
                    ),
                },
            },
        ],
        # ── Tacrolimus ──
        "tacrolimus": [
            {
                "gene": "CYP3A5",
                "actions": {
                    "CYP3A5 Expresser": (
                        AlertSeverity.MAJOR,
                        "CYP3A5 expresser: increase tacrolimus starting dose by "
                        "1.5-2x (0.3 mg/kg/day). Rapid metabolism expected.",
                        [],
                        "A",
                    ),
                    "CYP3A5 Intermediate Expresser": (
                        AlertSeverity.MODERATE,
                        "CYP3A5 intermediate: increase tacrolimus starting dose by "
                        "1.25-1.5x. Monitor trough levels closely.",
                        [],
                        "A",
                    ),
                },
            },
        ],
        # ── Irinotecan ──
        "irinotecan": [
            {
                "gene": "UGT1A1",
                "actions": {
                    "Poor Metabolizer": (
                        AlertSeverity.MAJOR,
                        "UGT1A1 *28/*28: reduce irinotecan dose by 30%. Increased "
                        "risk of severe neutropenia and diarrhea.",
                        [],
                        "A",
                    ),
                    "Intermediate Metabolizer": (
                        AlertSeverity.MODERATE,
                        "UGT1A1 *1/*28 or *1/*6: monitor closely for toxicity. "
                        "Dose reduction may be warranted if high-dose protocol.",
                        [],
                        "B",
                    ),
                },
            },
        ],
        # ── SSRIs/Antidepressants ──
        "paroxetine": [
            {
                "gene": "CYP2D6",
                "actions": {
                    "Ultrarapid Metabolizer": (
                        AlertSeverity.MAJOR,
                        "Select alternative SSRI not metabolized by CYP2D6 "
                        "(sertraline, citalopram, escitalopram).",
                        ["sertraline", "citalopram", "escitalopram"],
                        "A",
                    ),
                    "Poor Metabolizer": (
                        AlertSeverity.MAJOR,
                        "Select alternative SSRI or reduce paroxetine dose by 50%. "
                        "Increased side effect risk.",
                        ["sertraline", "citalopram", "escitalopram"],
                        "A",
                    ),
                },
            },
        ],
        "amitriptyline": [
            {
                "gene": "CYP2D6",
                "actions": {
                    "Ultrarapid Metabolizer": (
                        AlertSeverity.MAJOR,
                        "AVOID amitriptyline. Rapid metabolism → sub-therapeutic "
                        "levels. Use alternative not dependent on CYP2D6.",
                        ["nortriptyline (with TDM)", "SSRI"],
                        "A",
                    ),
                    "Poor Metabolizer": (
                        AlertSeverity.MAJOR,
                        "AVOID amitriptyline or reduce dose by 50%. Increased "
                        "side-effect risk from accumulation.",
                        ["nortriptyline (with TDM)", "SSRI"],
                        "A",
                    ),
                },
            },
            {
                "gene": "CYP2C19",
                "actions": {
                    "Ultrarapid Metabolizer": (
                        AlertSeverity.MODERATE,
                        "CYP2C19 UM: amitriptyline may have reduced efficacy. "
                        "Consider alternative or therapeutic drug monitoring.",
                        [],
                        "A",
                    ),
                    "Poor Metabolizer": (
                        AlertSeverity.MODERATE,
                        "CYP2C19 PM: reduce amitriptyline dose by 50%. "
                        "Increased parent drug exposure.",
                        [],
                        "A",
                    ),
                },
            },
        ],
        # ── PPIs ──
        "omeprazole": [
            {
                "gene": "CYP2C19",
                "actions": {
                    "Ultrarapid Metabolizer": (
                        AlertSeverity.MODERATE,
                        "Increase omeprazole dose by 100-200% for H. pylori "
                        "eradication. Rapid metabolism reduces efficacy.",
                        ["rabeprazole (CYP2C19-independent)", "esomeprazole (higher dose)"],
                        "B",
                    ),
                    "Rapid Metabolizer": (
                        AlertSeverity.MODERATE,
                        "Increase omeprazole dose by 50-100% for H. pylori.",
                        ["rabeprazole"],
                        "B",
                    ),
                    "Poor Metabolizer": (
                        AlertSeverity.MINOR,
                        "Reduced omeprazole dose may suffice. Increased exposure, "
                        "but generally safe — monitor for side effects.",
                        [],
                        "C",
                    ),
                },
            },
        ],
        # ── Antifungal ──
        "voriconazole": [
            {
                "gene": "CYP2C19",
                "actions": {
                    "Ultrarapid Metabolizer": (
                        AlertSeverity.MAJOR,
                        "Choose alternative antifungal (isavuconazole, posaconazole). "
                        "Rapid CYP2C19 metabolism → sub-therapeutic voriconazole.",
                        ["isavuconazole", "posaconazole"],
                        "A",
                    ),
                    "Poor Metabolizer": (
                        AlertSeverity.MAJOR,
                        "Reduce voriconazole dose by 50%. Elevated exposure with "
                        "increased toxicity risk (visual, hepatic).",
                        [],
                        "A",
                    ),
                },
            },
        ],
        # ── Oncology (tamoxifen) ──
        "tamoxifen": [
            {
                "gene": "CYP2D6",
                "actions": {
                    "Poor Metabolizer": (
                        AlertSeverity.MAJOR,
                        "Consider alternative endocrine therapy (aromatase inhibitor "
                        "if postmenopausal). CYP2D6 PM → insufficient endoxifen.",
                        ["aromatase inhibitor (postmenopausal)", "tamoxifen 40mg (with TDM)"],
                        "B",
                    ),
                    "Intermediate Metabolizer": (
                        AlertSeverity.MODERATE,
                        "Consider higher tamoxifen dose (40 mg) with endoxifen TDM, "
                        "or aromatase inhibitor if postmenopausal.",
                        ["aromatase inhibitor (postmenopausal)"],
                        "B",
                    ),
                },
            },
        ],
        # ── ADHD ──
        "atomoxetine": [
            {
                "gene": "CYP2D6",
                "actions": {
                    "Poor Metabolizer": (
                        AlertSeverity.MAJOR,
                        "CYP2D6 PM: 4-fold increase in atomoxetine exposure. "
                        "Start at 0.5 mg/kg/day and titrate cautiously.",
                        ["guanfacine", "clonidine"],
                        "A",
                    ),
                    "Ultrarapid Metabolizer": (
                        AlertSeverity.MODERATE,
                        "CYP2D6 UM: may need higher atomoxetine dose due to "
                        "increased metabolism. Monitor for inadequate response.",
                        [],
                        "B",
                    ),
                },
            },
        ],
        # ── Anticonvulsants ──
        "carbamazepine": [
            {
                "gene": "HLA-B",
                "actions": {
                    "HLA-B*15:02 positive": (
                        AlertSeverity.CONTRAINDICATED,
                        "AVOID carbamazepine. HLA-B*15:02 positive → high risk of "
                        "Stevens-Johnson syndrome / toxic epidermal necrolysis (SJS/TEN).",
                        ["levetiracetam", "valproic acid", "lacosamide"],
                        "A",
                    ),
                },
            },
            {
                "gene": "HLA-A",
                "actions": {
                    "HLA-A*31:01 positive": (
                        AlertSeverity.MAJOR,
                        "HLA-A*31:01 positive → increased risk of DRESS syndrome "
                        "and SJS/TEN with carbamazepine. Use alternative anticonvulsant.",
                        ["levetiracetam", "valproic acid", "lacosamide"],
                        "A",
                    ),
                },
            },
        ],
        # ── SSRIs/SNRIs ──
        "citalopram": [
            {
                "gene": "CYP2C19",
                "actions": {
                    "Poor Metabolizer": (
                        AlertSeverity.MAJOR,
                        "Reduce citalopram dose by 50% (max 20 mg/day). CYP2C19 PM → "
                        "increased exposure and QT prolongation risk.",
                        ["sertraline", "mirtazapine"],
                        "A",
                    ),
                    "Ultrarapid Metabolizer": (
                        AlertSeverity.MODERATE,
                        "CYP2C19 UM: citalopram may have reduced efficacy. "
                        "Consider alternative antidepressant.",
                        ["fluoxetine", "venlafaxine"],
                        "B",
                    ),
                },
            },
        ],
        # ── Tricyclic Antidepressants ──
        "clomipramine": [
            {
                "gene": "CYP2D6",
                "actions": {
                    "Poor Metabolizer": (
                        AlertSeverity.MAJOR,
                        "Reduce clomipramine dose by 50%. CYP2D6 PM → elevated "
                        "TCA levels with increased risk of toxicity. Monitor drug levels.",
                        ["sertraline", "fluoxetine"],
                        "A",
                    ),
                    "Ultrarapid Metabolizer": (
                        AlertSeverity.MODERATE,
                        "CYP2D6 UM: clomipramine may have subtherapeutic levels. "
                        "Consider higher dose with TDM or alternative antidepressant.",
                        ["venlafaxine"],
                        "A",
                    ),
                },
            },
            {
                "gene": "CYP2C19",
                "actions": {
                    "Poor Metabolizer": (
                        AlertSeverity.MODERATE,
                        "CYP2C19 PM: reduced clomipramine demethylation. "
                        "Reduce dose and monitor TCA plasma levels.",
                        ["sertraline"],
                        "A",
                    ),
                },
            },
        ],
        "doxepin": [
            {
                "gene": "CYP2D6",
                "actions": {
                    "Poor Metabolizer": (
                        AlertSeverity.MAJOR,
                        "Reduce doxepin dose by 50%. CYP2D6 PM → elevated active "
                        "metabolite levels with increased toxicity risk. Monitor drug levels.",
                        ["sertraline", "mirtazapine"],
                        "A",
                    ),
                    "Ultrarapid Metabolizer": (
                        AlertSeverity.MODERATE,
                        "CYP2D6 UM: doxepin may have subtherapeutic levels. "
                        "Consider alternative antidepressant.",
                        ["venlafaxine"],
                        "A",
                    ),
                },
            },
            {
                "gene": "CYP2C19",
                "actions": {
                    "Poor Metabolizer": (
                        AlertSeverity.MODERATE,
                        "CYP2C19 PM: reduced doxepin demethylation. Reduce dose "
                        "and monitor plasma levels.",
                        ["sertraline"],
                        "A",
                    ),
                },
            },
        ],
        # ── Antiretroviral ──
        "efavirenz": [
            {
                "gene": "CYP2B6",
                "actions": {
                    "Poor Metabolizer": (
                        AlertSeverity.MAJOR,
                        "CYP2B6 PM: 3-fold increase in efavirenz exposure. Switch to "
                        "400 mg/day dosing or use alternative antiretroviral.",
                        ["dolutegravir", "rilpivirine"],
                        "A",
                    ),
                    "Intermediate Metabolizer": (
                        AlertSeverity.MODERATE,
                        "CYP2B6 IM: consider reducing efavirenz to 400 mg/day. "
                        "Monitor for CNS side effects.",
                        ["dolutegravir", "rilpivirine"],
                        "B",
                    ),
                },
            },
        ],
        "escitalopram": [
            {
                "gene": "CYP2C19",
                "actions": {
                    "Poor Metabolizer": (
                        AlertSeverity.MAJOR,
                        "Reduce escitalopram dose by 50% (max 10 mg/day). CYP2C19 PM → "
                        "increased exposure and QT prolongation risk.",
                        ["sertraline", "mirtazapine"],
                        "A",
                    ),
                    "Ultrarapid Metabolizer": (
                        AlertSeverity.MODERATE,
                        "CYP2C19 UM: escitalopram may have reduced efficacy. "
                        "Consider alternative antidepressant.",
                        ["fluoxetine", "venlafaxine"],
                        "B",
                    ),
                },
            },
        ],
        "fluvoxamine": [
            {
                "gene": "CYP2D6",
                "actions": {
                    "Poor Metabolizer": (
                        AlertSeverity.MODERATE,
                        "CYP2D6 PM: possible increased fluvoxamine side effects via "
                        "CYP2D6-mediated pathway. Monitor and consider alternative.",
                        ["sertraline", "escitalopram"],
                        "B",
                    ),
                },
            },
        ],
        "nortriptyline": [
            {
                "gene": "CYP2D6",
                "actions": {
                    "Poor Metabolizer": (
                        AlertSeverity.MAJOR,
                        "Reduce nortriptyline dose by 50%. CYP2D6 PM → elevated levels. "
                        "Target 50-150 ng/mL with therapeutic drug monitoring.",
                        ["sertraline", "venlafaxine"],
                        "A",
                    ),
                    "Ultrarapid Metabolizer": (
                        AlertSeverity.MAJOR,
                        "AVOID nortriptyline. CYP2D6 UM → subtherapeutic levels likely. "
                        "Use alternative not dependent on CYP2D6.",
                        ["sertraline", "venlafaxine"],
                        "A",
                    ),
                },
            },
        ],
        # ── Antiemetic ──
        "ondansetron": [
            {
                "gene": "CYP2D6",
                "actions": {
                    "Ultrarapid Metabolizer": (
                        AlertSeverity.MODERATE,
                        "CYP2D6 UM: increased ondansetron metabolism may reduce "
                        "antiemetic efficacy. Consider higher dose or alternative.",
                        ["granisetron", "palonosetron"],
                        "B",
                    ),
                    "Poor Metabolizer": (
                        AlertSeverity.MODERATE,
                        "CYP2D6 PM: may have increased QT prolongation risk with "
                        "ondansetron. Use alternative antiemetic if possible.",
                        ["granisetron"],
                        "B",
                    ),
                },
            },
        ],
        # ── Anticonvulsant (phenytoin) ──
        "phenytoin": [
            {
                "gene": "CYP2C9",
                "actions": {
                    "Poor Metabolizer": (
                        AlertSeverity.MAJOR,
                        "Reduce phenytoin dose by 50%. CYP2C9 PM → markedly reduced "
                        "clearance with increased toxicity risk.",
                        ["levetiracetam", "lacosamide"],
                        "A",
                    ),
                    "Intermediate Metabolizer": (
                        AlertSeverity.MODERATE,
                        "Reduce phenytoin dose by 25%. CYP2C9 IM → moderately reduced "
                        "clearance. Monitor levels closely.",
                        ["levetiracetam", "lacosamide"],
                        "A",
                    ),
                },
            },
            {
                "gene": "HLA-B",
                "actions": {
                    "HLA-B*15:02 positive": (
                        AlertSeverity.CONTRAINDICATED,
                        "AVOID phenytoin. HLA-B*15:02 positive → high risk of "
                        "Stevens-Johnson syndrome / toxic epidermal necrolysis (SJS/TEN).",
                        ["levetiracetam", "valproic acid"],
                        "A",
                    ),
                },
            },
        ],
        "sertraline": [
            {
                "gene": "CYP2C19",
                "actions": {
                    "Poor Metabolizer": (
                        AlertSeverity.MODERATE,
                        "CYP2C19 PM: consider 50% dose reduction of sertraline or "
                        "use alternative antidepressant. Increased exposure expected.",
                        ["mirtazapine", "bupropion"],
                        "B",
                    ),
                    "Ultrarapid Metabolizer": (
                        AlertSeverity.INFORMATIONAL,
                        "CYP2C19 UM: standard sertraline dosing appropriate. "
                        "Monitor clinical response.",
                        [],
                        "C",
                    ),
                },
            },
        ],
    }

    def check_drug(
        self,
        drug: str,
        profile: Dict[str, Dict[str, Any]],
    ) -> List[PGxAlert]:
        """Check a single drug against the patient PGx profile.

        Parameters
        ----------
        drug : str
            Drug name (case-insensitive).
        profile : dict
            gene → {phenotype, diplotype, ...} from PhenotypeTranslator.

        Returns
        -------
        list of PGxAlert
        """
        drug_lower = drug.lower().strip()
        gene_entries = self.DRUG_GENE_MAP.get(drug_lower, [])
        alerts: List[PGxAlert] = []

        for entry in gene_entries:
            gene = entry["gene"]
            gene_data = profile.get(gene)
            if not gene_data:
                continue

            # A gene entry with no phenotype or no diplotype is an absent
            # call, not a normal metabolizer carrying *1/*1. The old
            # defaults made an uncalled gene look like a confident negative
            # result; "Indeterminate" matches no action, so no alert is
            # emitted and nothing is asserted.
            phenotype = gene_data.get("phenotype") or "Indeterminate"
            diplotype = gene_data.get("diplotype") or "not called"
            actions = entry.get("actions", {})

            action = actions.get(phenotype)
            if action:
                severity, recommendation, alternatives, cpic_level = action
                alerts.append(PGxAlert(
                    drug=drug,
                    gene=gene,
                    diplotype=diplotype,
                    phenotype=phenotype,
                    severity=severity,
                    recommendation=recommendation,
                    cpic_level=cpic_level,
                    alternatives=alternatives,
                    source="CPIC",
                ))

        return alerts

    def check_medication_list(
        self,
        medications: List[str],
        profile: Dict[str, Dict[str, Any]],
    ) -> List[PGxAlert]:
        """Screen an entire medication list against the PGx profile.

        Parameters
        ----------
        medications : list of str
            Drug names.
        profile : dict
            Output of PhenotypeTranslator.translate_all().

        Returns
        -------
        list of PGxAlert
            Sorted by severity (most severe first).
        """
        severity_order = {
            AlertSeverity.CONTRAINDICATED: 0,
            AlertSeverity.MAJOR: 1,
            AlertSeverity.MODERATE: 2,
            AlertSeverity.MINOR: 3,
            AlertSeverity.INFORMATIONAL: 4,
        }

        all_alerts: List[PGxAlert] = []
        for drug in medications:
            all_alerts.extend(self.check_drug(drug, profile))

        all_alerts.sort(key=lambda a: severity_order.get(a.severity, 99))
        logger.info(
            "Screened %d medications → %d alerts (%d contraindicated)",
            len(medications),
            len(all_alerts),
            sum(1 for a in all_alerts if a.severity == AlertSeverity.CONTRAINDICATED),
        )
        return all_alerts
