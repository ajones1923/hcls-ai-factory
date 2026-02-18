"""
VCF Parser - Extract variants from VCF files into evidence objects.
"""
import gzip
from collections.abc import Generator
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from loguru import logger

try:
    from cyvcf2 import VCF
    HAS_CYVCF2 = True
except ImportError:
    HAS_CYVCF2 = False
    logger.warning("cyvcf2 not installed. Using fallback parser.")


@dataclass
class VariantEvidence:
    """Represents a single variant as an evidence object."""

    # Core variant information
    variant_id: str
    chrom: str
    pos: int
    ref: str
    alt: str
    qual: float
    filter_status: str

    # Genotype information
    genotype: str = "."
    depth: int = 0

    # Annotation (populated by annotator)
    gene: str | None = None
    consequence: str | None = None
    impact: str | None = None
    hgvs_c: str | None = None
    hgvs_p: str | None = None

    # Clinical information
    clinical_significance: str | None = None
    rsid: str | None = None
    disease_associations: str | None = None

    # AlphaMissense predictions
    am_pathogenicity: float | None = None  # Score 0-1
    am_class: str | None = None  # likely_benign, ambiguous, likely_pathogenic

    # Metadata
    source: str = "VCF"
    tags: list = field(default_factory=list)

    @property
    def text_summary(self) -> str:
        """Generate human-readable text summary for embedding."""
        parts = [
            f"Variant at {self.chrom}:{self.pos}",
            f"{self.ref}>{self.alt}",
        ]

        if self.gene:
            parts.append(f"in gene {self.gene}")

        if self.consequence:
            parts.append(f"({self.consequence})")

        if self.impact:
            parts.append(f"with {self.impact} impact")

        if self.hgvs_p:
            parts.append(f"protein change: {self.hgvs_p}")

        if self.rsid and self.rsid != ".":
            parts.append(f"dbSNP: {self.rsid}")

        if self.clinical_significance:
            parts.append(f"clinical: {self.clinical_significance}")

        if self.disease_associations:
            parts.append(f"associated with: {self.disease_associations}")

        if self.am_pathogenicity is not None:
            parts.append(f"AlphaMissense: {self.am_pathogenicity:.3f} ({self.am_class})")

        # Include disease tags if present
        disease_tags = [t.replace("disease:", "") for t in self.tags if t.startswith("disease:")]
        if disease_tags:
            parts.append(f"linked diseases: {'; '.join(disease_tags)}")

        parts.append(f"quality score {self.qual:.1f}")

        if self.genotype and self.genotype != ".":
            gt_desc = {
                "0/0": "homozygous reference",
                "0/1": "heterozygous",
                "1/0": "heterozygous",
                "1/1": "homozygous alternate",
            }.get(self.genotype, self.genotype)
            parts.append(f"genotype: {gt_desc}")

        return " ".join(parts)

    def to_dict(self) -> dict:
        """Convert to dictionary for storage."""
        return {
            "variant_id": self.variant_id,
            "chrom": self.chrom,
            "pos": self.pos,
            "ref": self.ref,
            "alt": self.alt,
            "qual": self.qual,
            "filter_status": self.filter_status,
            "genotype": self.genotype,
            "depth": self.depth,
            "gene": self.gene,
            "consequence": self.consequence,
            "impact": self.impact,
            "hgvs_c": self.hgvs_c,
            "hgvs_p": self.hgvs_p,
            "clinical_significance": self.clinical_significance,
            "rsid": self.rsid,
            "am_pathogenicity": self.am_pathogenicity,
            "am_class": self.am_class,
            "source": self.source,
            "tags": self.tags,
            "text_summary": self.text_summary,
        }


class VCFParser:
    """Parse VCF files and yield VariantEvidence objects."""

    def __init__(
        self,
        vcf_path: Path,
        min_qual: float = 20.0,
        include_chromosomes: list | None = None,
        max_variants: int | None = None,
        exclude_ref_calls: bool = True,
    ):
        self.vcf_path = Path(vcf_path)
        self.min_qual = min_qual
        self.include_chromosomes = include_chromosomes
        self.max_variants = max_variants
        self.exclude_ref_calls = exclude_ref_calls

        if not self.vcf_path.exists():
            raise FileNotFoundError(f"VCF file not found: {self.vcf_path}")

        logger.info(f"Initialized VCF parser for: {self.vcf_path}")
        if exclude_ref_calls:
            logger.info("Filtering: excluding homozygous reference (0/0) calls")

    def parse(self) -> Generator[VariantEvidence, None, None]:
        """Parse VCF and yield VariantEvidence objects."""
        if HAS_CYVCF2:
            yield from self._parse_cyvcf2()
        else:
            yield from self._parse_fallback()

    def _parse_cyvcf2(self) -> Generator[VariantEvidence, None, None]:
        """Parse using cyvcf2 (fast, recommended)."""
        vcf = VCF(str(self.vcf_path))
        count = 0

        for variant in vcf:
            # Apply filters
            if self.max_variants and count >= self.max_variants:
                break

            if variant.QUAL and self.min_qual > variant.QUAL:
                continue

            chrom = variant.CHROM
            if self.include_chromosomes and chrom not in self.include_chromosomes:
                continue

            # Handle multi-allelic variants
            for _i, alt in enumerate(variant.ALT):
                if alt is None:
                    continue

                # Truncate long alleles in ID to keep it reasonable
                ref_short = variant.REF[:50] if len(variant.REF) > 50 else variant.REF
                alt_short = alt[:50] if len(alt) > 50 else alt
                variant_id = f"{chrom}_{variant.POS}_{ref_short}_{alt_short}"

                # Extract genotype
                gt = "."
                depth = 0
                if variant.genotypes:
                    gt_tuple = variant.genotypes[0][:2]  # First sample
                    gt = f"{gt_tuple[0]}/{gt_tuple[1]}"

                    # Skip homozygous reference calls (0/0) if filtering enabled
                    if self.exclude_ref_calls and gt == "0/0":
                        continue

                    # Also skip no-call genotypes (./.  or -1/-1)
                    if gt_tuple[0] == -1 or gt_tuple[1] == -1:
                        continue

                    try:
                        depth = variant.format("DP")[0] if variant.format("DP") is not None else 0
                    except (KeyError, TypeError, IndexError):
                        depth = 0

                # Extract rsID
                rsid = variant.ID if variant.ID else None

                evidence = VariantEvidence(
                    variant_id=variant_id,
                    chrom=chrom,
                    pos=variant.POS,
                    ref=variant.REF,
                    alt=alt,
                    qual=variant.QUAL or 0.0,
                    filter_status=variant.FILTER or "PASS",
                    genotype=gt,
                    depth=depth,
                    rsid=rsid,
                    source="GIAB WGS VCF (GRCh38)",
                    tags=["giab", "wgs", "hg002", "grch38"],
                )

                yield evidence
                count += 1

        vcf.close()
        logger.info(f"Parsed {count} variants from VCF")

    def _parse_fallback(self) -> Generator[VariantEvidence, None, None]:
        """Fallback parser using standard library (slower)."""
        logger.warning("Using fallback VCF parser. Install cyvcf2 for better performance.")

        open_func = gzip.open if str(self.vcf_path).endswith('.gz') else open
        count = 0

        with open_func(self.vcf_path, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue

                if self.max_variants and count >= self.max_variants:
                    break

                fields = line.strip().split('\t')
                if len(fields) < 8:
                    continue

                chrom, pos, rsid, ref, alt, qual, filt, info = fields[:8]

                # Parse quality
                try:
                    qual_float = float(qual) if qual != '.' else 0.0
                except ValueError:
                    qual_float = 0.0

                if qual_float < self.min_qual:
                    continue

                if self.include_chromosomes and chrom not in self.include_chromosomes:
                    continue

                # Handle multi-allelic
                for alt_allele in alt.split(','):
                    # Truncate long alleles in ID to keep it reasonable
                    ref_short = ref[:50] if len(ref) > 50 else ref
                    alt_short = alt_allele[:50] if len(alt_allele) > 50 else alt_allele
                    variant_id = f"{chrom}_{pos}_{ref_short}_{alt_short}"

                    # Parse genotype if present
                    gt = "."
                    depth = 0
                    if len(fields) > 9:
                        format_fields = fields[8].split(':')
                        sample_fields = fields[9].split(':')
                        if 'GT' in format_fields:
                            gt_idx = format_fields.index('GT')
                            gt = sample_fields[gt_idx].replace('|', '/')

                            # Skip homozygous reference calls (0/0) if filtering enabled
                            if self.exclude_ref_calls and gt == "0/0":
                                continue

                            # Also skip no-call genotypes
                            if gt == "./." or gt == ".|.":
                                continue

                        if 'DP' in format_fields:
                            dp_idx = format_fields.index('DP')
                            try:
                                depth = int(sample_fields[dp_idx])
                            except (ValueError, IndexError):
                                depth = 0

                    evidence = VariantEvidence(
                        variant_id=variant_id,
                        chrom=chrom,
                        pos=int(pos),
                        ref=ref,
                        alt=alt_allele,
                        qual=qual_float,
                        filter_status=filt,
                        genotype=gt,
                        depth=depth,
                        rsid=rsid if rsid != '.' else None,
                        source="GIAB WGS VCF (GRCh38)",
                        tags=["giab", "wgs", "hg002", "grch38"],
                    )

                    yield evidence
                    count += 1

        logger.info(f"Parsed {count} variants from VCF (fallback parser)")

    def count_variants(self) -> int:
        """Count total variants in VCF (for progress bars)."""
        if HAS_CYVCF2:
            vcf = VCF(str(self.vcf_path))
            # Use index if available for fast counting
            count = sum(1 for _ in vcf)
            vcf.close()
            return count
        else:
            open_func = gzip.open if str(self.vcf_path).endswith('.gz') else open
            with open_func(self.vcf_path, 'rt') as f:
                return sum(1 for line in f if not line.startswith('#'))
