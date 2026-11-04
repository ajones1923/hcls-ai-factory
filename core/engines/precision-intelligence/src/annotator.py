"""
Variant Annotator - Add gene names, consequences, and clinical information to variants.
"""
import subprocess
import tempfile
import json
import gzip
from pathlib import Path
from typing import Optional, Dict, List, Tuple
from dataclasses import dataclass
from loguru import logger

from .vcf_parser import VariantEvidence


@dataclass
class AnnotationResult:
    """Result of variant annotation."""
    gene: Optional[str] = None
    consequence: Optional[str] = None
    impact: Optional[str] = None
    hgvs_c: Optional[str] = None
    hgvs_p: Optional[str] = None
    clinical_significance: Optional[str] = None
    biotype: Optional[str] = None
    sift: Optional[str] = None
    polyphen: Optional[str] = None


class VariantAnnotator:
    """
    Annotate variants with gene names, consequences, and clinical information.

    Supports multiple annotation strategies:
    1. VEP (Ensembl Variant Effect Predictor) - most comprehensive
    2. BED file overlap - lightweight gene name annotation
    3. Pre-computed annotation database lookup
    """

    def __init__(
        self,
        use_vep: bool = True,
        vep_cache_dir: Optional[Path] = None,
        gene_bed_file: Optional[Path] = None,
        species: str = "homo_sapiens",
        assembly: str = "GRCh38",
    ):
        self.use_vep = use_vep
        self.vep_cache_dir = vep_cache_dir
        self.gene_bed_file = gene_bed_file
        self.species = species
        self.assembly = assembly

        # Cache for gene lookups
        self._gene_cache: Dict[str, str] = {}
        self._gene_intervals: Optional[List] = None

        if gene_bed_file and Path(gene_bed_file).exists():
            self._load_gene_bed(gene_bed_file)

    def _load_gene_bed(self, bed_file: Path) -> None:
        """Load gene regions from BED file for fast lookup."""
        logger.info(f"Loading gene regions from: {bed_file}")
        self._gene_intervals = []

        with open(bed_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    chrom, start, end, gene = parts[:4]
                    self._gene_intervals.append({
                        'chrom': chrom,
                        'start': int(start),
                        'end': int(end),
                        'gene': gene,
                    })

        logger.info(f"Loaded {len(self._gene_intervals)} gene regions")

    def annotate(self, variant: VariantEvidence) -> VariantEvidence:
        """Annotate a single variant."""

        # Try VEP first if enabled
        if self.use_vep:
            annotation = self._annotate_with_vep_api(variant)
            if annotation and annotation.gene:
                return self._apply_annotation(variant, annotation)

        # Fallback to BED file lookup
        if self._gene_intervals:
            gene = self._lookup_gene_bed(variant.chrom, variant.pos)
            if gene:
                variant.gene = gene
                return variant

        return variant

    def annotate_batch(self, variants: List[VariantEvidence]) -> List[VariantEvidence]:
        """Annotate a batch of variants (more efficient for VEP)."""
        if self.use_vep and len(variants) > 1:
            return self._annotate_batch_vep(variants)
        else:
            return [self.annotate(v) for v in variants]

    def _annotate_with_vep_api(self, variant: VariantEvidence) -> Optional[AnnotationResult]:
        """
        Annotate using Ensembl VEP REST API.
        Note: For production, use local VEP installation for speed.
        """
        try:
            import requests
        except ImportError:
            logger.warning("requests not installed. Cannot use VEP API.")
            return None

        # Build VEP query
        query = f"{variant.chrom}:{variant.pos}:{variant.ref}:{variant.alt}"

        url = f"https://rest.ensembl.org/vep/{self.species}/hgvs/{query}"
        headers = {"Content-Type": "application/json"}

        try:
            response = requests.get(
                url,
                headers=headers,
                params={"content-type": "application/json"},
                timeout=10,
            )

            if response.status_code == 200:
                data = response.json()
                if data and len(data) > 0:
                    return self._parse_vep_response(data[0])
        except Exception as e:
            logger.debug(f"VEP API error for {variant.variant_id}: {e}")

        return None

    def _annotate_batch_vep(self, variants: List[VariantEvidence]) -> List[VariantEvidence]:
        """Annotate batch using VEP REST API POST endpoint."""
        try:
            import requests
        except ImportError:
            return [self.annotate(v) for v in variants]

        # Build batch query
        queries = []
        for v in variants:
            queries.append(f"{v.chrom} {v.pos} . {v.ref} {v.alt} . . .")

        url = f"https://rest.ensembl.org/vep/{self.species}/region"
        headers = {"Content-Type": "application/json", "Accept": "application/json"}

        try:
            response = requests.post(
                url,
                headers=headers,
                json={"variants": queries},
                timeout=60,
            )

            if response.status_code == 200:
                results = response.json()
                # Map results back to variants
                result_map = {}
                for r in results:
                    key = f"{r.get('seq_region_name')}_{r.get('start')}_{r.get('allele_string', '').replace('/', '_')}"
                    result_map[key] = self._parse_vep_response(r)

                for v in variants:
                    key = f"{v.chrom}_{v.pos}_{v.ref}_{v.alt}"
                    if key in result_map:
                        self._apply_annotation(v, result_map[key])

        except Exception as e:
            logger.warning(f"VEP batch API error: {e}")

        return variants

    def _parse_vep_response(self, data: dict) -> AnnotationResult:
        """Parse VEP API response into AnnotationResult."""
        result = AnnotationResult()

        # Get most severe consequence
        if 'most_severe_consequence' in data:
            result.consequence = data['most_severe_consequence'].replace('_', ' ')

        # Parse transcript consequences
        consequences = data.get('transcript_consequences', [])
        if consequences:
            # Use first consequence (usually most severe)
            tc = consequences[0]
            result.gene = tc.get('gene_symbol')
            result.impact = tc.get('impact')
            result.biotype = tc.get('biotype')
            result.hgvs_c = tc.get('hgvsc')
            result.hgvs_p = tc.get('hgvsp')

            # SIFT/PolyPhen predictions
            if 'sift_prediction' in tc:
                result.sift = tc['sift_prediction']
            if 'polyphen_prediction' in tc:
                result.polyphen = tc['polyphen_prediction']

        # Clinical significance from ClinVar
        colocated = data.get('colocated_variants', [])
        for cv in colocated:
            if 'clin_sig' in cv:
                result.clinical_significance = ', '.join(cv['clin_sig'])
                break

        return result

    def _apply_annotation(self, variant: VariantEvidence, annotation: AnnotationResult) -> VariantEvidence:
        """Apply annotation result to variant."""
        variant.gene = annotation.gene
        variant.consequence = annotation.consequence
        variant.impact = annotation.impact
        variant.hgvs_c = annotation.hgvs_c
        variant.hgvs_p = annotation.hgvs_p
        variant.clinical_significance = annotation.clinical_significance
        return variant

    def _lookup_gene_bed(self, chrom: str, pos: int) -> Optional[str]:
        """Look up gene name from BED intervals."""
        # Check cache first
        cache_key = f"{chrom}:{pos}"
        if cache_key in self._gene_cache:
            return self._gene_cache[cache_key]

        # Search intervals
        if self._gene_intervals:
            for interval in self._gene_intervals:
                if (interval['chrom'] == chrom and
                    interval['start'] <= pos <= interval['end']):
                    self._gene_cache[cache_key] = interval['gene']
                    return interval['gene']

        return None


class LocalVEPAnnotator(VariantAnnotator):
    """
    Annotate using local VEP installation (Docker or native).
    Much faster than REST API for large datasets.
    """

    def __init__(
        self,
        vep_docker_image: str = "ensemblorg/ensembl-vep:release_110.1",
        cache_dir: Path = Path("/data/vep_cache"),
        **kwargs
    ):
        super().__init__(**kwargs)
        self.vep_docker_image = vep_docker_image
        self.cache_dir = cache_dir

    def annotate_vcf_file(
        self,
        input_vcf: Path,
        output_vcf: Path,
        extra_args: Optional[List[str]] = None,
    ) -> Path:
        """
        Annotate entire VCF file using local VEP.
        Returns path to annotated VCF.
        """
        cmd = [
            "docker", "run", "--rm",
            "-v", f"{input_vcf.parent}:/data/input:ro",
            "-v", f"{output_vcf.parent}:/data/output",
            "-v", f"{self.cache_dir}:/opt/vep/.vep",
            self.vep_docker_image,
            "vep",
            "-i", f"/data/input/{input_vcf.name}",
            "-o", f"/data/output/{output_vcf.name}",
            "--cache",
            "--offline",
            "--assembly", self.assembly,
            "--species", self.species,
            "--format", "vcf",
            "--vcf",
            "--symbol",
            "--biotype",
            "--numbers",
            "--canonical",
            "--protein",
            "--hgvs",
            "--af",
            "--sift", "b",
            "--polyphen", "b",
        ]

        if extra_args:
            cmd.extend(extra_args)

        logger.info(f"Running VEP annotation: {' '.join(cmd)}")

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
            if result.returncode != 0:
                logger.error(f"VEP failed: {result.stderr}")
                raise RuntimeError(f"VEP annotation failed: {result.stderr}")
            logger.info("VEP annotation complete")
            return output_vcf
        except subprocess.TimeoutExpired:
            raise RuntimeError("VEP annotation timed out after 1 hour")


# ---------------------------------------------------------------------------
# ClinVar aggregate germline classification -> internal impact tier
#
# Primary sources (retrieved 2026-11-03):
#  * ClinVar's own release VCF header,
#    https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
#    (##fileDate=2026-11-01):
#      CLNSIG     "Aggregate germline classification for this single variant;
#                  multiple values are separated by a vertical bar"
#      CLNSIGCONF "Conflicting germline classification for this single variant"
#    The distinct CLNSIG values in that release include `Pathogenic`,
#    `Likely_pathogenic`, `Pathogenic/Likely_pathogenic` and
#    `Conflicting_classifications_of_pathogenicity` as four SEPARATE aggregate
#    terms; the conflicting term travels with
#    `CLNREVSTAT=criteria_provided,_conflicting_classifications` and a
#    CLNSIGCONF field enumerating the submissions that disagree. It is
#    therefore not a pathogenic call.
#  * NCBI, ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/README (updated
#    2026-01-13): "ClinicalSignificance  character, semicolon-separated list of
#    aggregate values of germline classification calculated for this variant."
#  * NCBI ClinVar_VCV.xsd, typeSubmitterReviewStatusValue, enumerates
#    "criteria provided, conflicting classifications".
#  * NCBI, ncbi/clinvar `somatic_classifications/ClassificationOnClinVar.md`:
#    a "Conflicting classifications of ..." aggregate is what ClinVar reports
#    when submitters' terms disagree across the pathogenic / uncertain / benign
#    boundary, with the conflicting terms and submission counts given.
#
# Mapping policy, and what it deliberately does NOT claim:
#  * Pathogenic, Likely pathogenic and their aggregates stay in the HIGH tier
#    (the tier the chat UI's triage filter uses) but carry distinct
#    `consequence` labels, so the ACMG/AMP difference between a reportable
#    Pathogenic call and a Likely pathogenic one survives into the report.
#  * Conflicting classifications of pathogenicity sits with Uncertain
#    significance (MODERATE) and is tagged `clinvar:conflicting`. ClinVar's own
#    aggregation rule makes it a disagreement, not a consensus pathogenic call.
#  * Risk-allele and non-classification terms (drug response, risk factor,
#    association, Affects, protective, other, not provided) are recognised but
#    given NO impact tier: assigning one would be a severity claim no retrieved
#    document supports. See CLINICAL_REVIEW.md §13.
#  * An unrecognised term leaves `impact` untouched and logs a warning, rather
#    than being guessed into a tier.
#
# Keys are lower-cased exact ClinVar terms (underscores in the VCF are spaces
# in the tab-delimited release).
#
# Each value is (impact tier, consequence label, tag). The tag is what makes the
# tier reachable downstream: `consequence` is overwritten by ClinVar's molecular
# `Type` column whenever that is present, so without the tag the
# Pathogenic / Likely pathogenic / Conflicting distinction would again be
# invisible to any consumer that reads only `impact`.
_CLINVAR_IMPACT_BY_TERM: Dict[str, Tuple[Optional[str], Optional[str], Optional[str]]] = {
    # --- pathogenic family -------------------------------------------------
    'pathogenic': ('HIGH', 'pathogenic variant', 'clinvar:pathogenic'),
    'pathogenic, low penetrance': (
        'HIGH', 'pathogenic variant, low penetrance', 'clinvar:pathogenic_low_penetrance'),
    'pathogenic/likely pathogenic': (
        'HIGH', 'pathogenic/likely pathogenic variant', 'clinvar:pathogenic_likely_pathogenic'),
    'pathogenic/likely pathogenic/pathogenic, low penetrance': (
        'HIGH', 'pathogenic/likely pathogenic variant, low penetrance',
        'clinvar:pathogenic_likely_pathogenic'),
    'pathogenic/likely pathogenic/established risk allele': (
        'HIGH', 'pathogenic/likely pathogenic variant', 'clinvar:pathogenic_likely_pathogenic'),
    'pathogenic/pathogenic, low penetrance': (
        'HIGH', 'pathogenic variant, low penetrance', 'clinvar:pathogenic_low_penetrance'),
    'likely pathogenic': ('HIGH', 'likely pathogenic variant', 'clinvar:likely_pathogenic'),
    'likely pathogenic, low penetrance': (
        'HIGH', 'likely pathogenic variant, low penetrance',
        'clinvar:likely_pathogenic_low_penetrance'),
    'likely pathogenic/likely risk allele': (
        'HIGH', 'likely pathogenic variant', 'clinvar:likely_pathogenic'),
    # --- disagreement between submitters ----------------------------------
    'conflicting classifications of pathogenicity': (
        'MODERATE', 'variant with conflicting classifications of pathogenicity',
        'clinvar:conflicting'),
    # --- uncertain ---------------------------------------------------------
    'uncertain significance': ('MODERATE', 'VUS', 'clinvar:uncertain_significance'),
    'uncertain significance/uncertain risk allele': (
        'MODERATE', 'VUS', 'clinvar:uncertain_significance'),
    'uncertain risk allele': (
        'MODERATE', 'uncertain risk allele', 'clinvar:uncertain_risk_allele'),
    # --- benign family -----------------------------------------------------
    'benign': ('LOW', 'benign variant', 'clinvar:benign'),
    'likely benign': ('LOW', 'likely benign variant', 'clinvar:likely_benign'),
    'benign/likely benign': ('LOW', 'benign/likely benign variant', 'clinvar:benign_likely_benign'),
    # --- recognised, but no severity tier asserted ------------------------
    'established risk allele': (None, 'established risk allele', 'clinvar:established_risk_allele'),
    'likely risk allele': (None, 'likely risk allele', 'clinvar:likely_risk_allele'),
    'risk allele': (None, 'risk allele', 'clinvar:risk_allele'),
    'risk factor': (None, 'risk factor', 'clinvar:risk_factor'),
    'drug response': (None, 'drug response', 'clinvar:drug_response'),
    'association': (None, 'association', 'clinvar:association'),
    'association not found': (None, 'association not found', 'clinvar:association_not_found'),
    'affects': (None, 'affects', 'clinvar:affects'),
    'protective': (None, 'protective', 'clinvar:protective'),
    'confers sensitivity': (None, 'confers sensitivity', 'clinvar:confers_sensitivity'),
    'other': (None, None, 'clinvar:other'),
    'not provided': (None, None, 'clinvar:not_provided'),
    'no classification for the single variant': (None, None, 'clinvar:no_classification'),
    'no classifications from unflagged records': (None, None, 'clinvar:no_classification'),
}

# Most-severe-wins when ClinVar reports several aggregate values for one variant
# (e.g. `Pathogenic|other`, `Benign/Likely_benign|drug_response`).
_IMPACT_SEVERITY = {'HIGH': 3, 'MODERATE': 2, 'LOW': 1}

_warned_clinvar_terms: set = set()


def map_clinvar_significance(
    clinical_significance: Optional[str],
) -> Tuple[Optional[str], Optional[str], List[str]]:
    """Map a ClinVar aggregate germline classification to (impact, consequence, tags).

    `clinical_significance` is the raw ClinVar string: one or more aggregate
    terms separated by ';' (tab-delimited release) or '|' (VCF). Returns
    ``(None, None, [])`` for an empty or wholly unrecognised value, so the
    caller leaves any existing impact alone rather than inventing one.

    See `_CLINVAR_IMPACT_BY_TERM` for the primary sources and the mapping policy.
    """
    if not clinical_significance:
        return None, None, []

    raw_terms = [
        t.strip()
        for part in str(clinical_significance).split(';')
        for t in part.split('|')
        if t.strip()
    ]

    best_impact: Optional[str] = None
    best_consequence: Optional[str] = None
    tags: List[str] = []
    recognised = False

    for term in raw_terms:
        key = term.lower().replace('_', ' ')
        if key not in _CLINVAR_IMPACT_BY_TERM:
            if key not in _warned_clinvar_terms:
                _warned_clinvar_terms.add(key)
                logger.warning(
                    f"Unrecognised ClinVar germline classification {term!r}: "
                    "no impact tier assigned"
                )
            continue
        recognised = True
        impact, consequence, tag = _CLINVAR_IMPACT_BY_TERM[key]
        if tag and tag not in tags:
            tags.append(tag)
        if impact is not None and (
            best_impact is None
            or _IMPACT_SEVERITY[impact] > _IMPACT_SEVERITY[best_impact]
        ):
            best_impact = impact
            best_consequence = consequence
        elif best_consequence is None and consequence is not None:
            best_consequence = consequence

    if not recognised:
        return None, None, []
    return best_impact, best_consequence, tags


class ClinVarAnnotator:
    """
    Fast local ClinVar annotation using pre-downloaded variant_summary.txt.gz.
    Much faster than VEP API calls for matching known pathogenic variants.
    """

    def __init__(
        self,
        clinvar_file: Path,
        assembly: str = "GRCh38",
    ):
        self.clinvar_file = Path(clinvar_file)
        self.assembly = assembly
        self._variant_db: Dict[str, dict] = {}
        self._loaded = False

        if not self.clinvar_file.exists():
            raise FileNotFoundError(f"ClinVar file not found: {self.clinvar_file}")

    def load(self, progress_callback=None) -> int:
        """
        Load ClinVar database into memory for fast lookups.
        Returns number of variants loaded.
        """
        if self._loaded:
            return len(self._variant_db)

        logger.info(f"Loading ClinVar database from: {self.clinvar_file}")
        count = 0
        open_func = gzip.open if str(self.clinvar_file).endswith('.gz') else open

        with open_func(self.clinvar_file, 'rt') as f:
            header = None
            for line_num, line in enumerate(f):
                if line_num == 0:
                    header = line.strip().split('\t')
                    continue

                fields = line.strip().split('\t')
                if len(fields) < 33:
                    continue

                # Only load variants for our assembly
                assembly_idx = header.index('Assembly') if 'Assembly' in header else 16
                if fields[assembly_idx] != self.assembly:
                    continue

                try:
                    # Extract key fields
                    chrom = fields[header.index('Chromosome')] if 'Chromosome' in header else fields[18]
                    pos_vcf = fields[header.index('PositionVCF')] if 'PositionVCF' in header else fields[30]
                    ref_vcf = fields[header.index('ReferenceAlleleVCF')] if 'ReferenceAlleleVCF' in header else fields[31]
                    alt_vcf = fields[header.index('AlternateAlleleVCF')] if 'AlternateAlleleVCF' in header else fields[32]

                    # Skip entries without valid position
                    if not pos_vcf or pos_vcf == '-1' or pos_vcf == 'na':
                        continue

                    # Normalize chromosome
                    if not chrom.startswith('chr'):
                        chrom = f"chr{chrom}"

                    # Create lookup key
                    key = f"{chrom}_{pos_vcf}_{ref_vcf}_{alt_vcf}"

                    # Extract annotation data
                    gene_symbol = fields[header.index('GeneSymbol')] if 'GeneSymbol' in header else fields[4]
                    clin_sig = fields[header.index('ClinicalSignificance')] if 'ClinicalSignificance' in header else fields[6]
                    phenotypes = fields[header.index('PhenotypeList')] if 'PhenotypeList' in header else fields[13]
                    rs_id = fields[header.index('RS# (dbSNP)')] if 'RS# (dbSNP)' in header else fields[9]
                    review_status = fields[header.index('ReviewStatus')] if 'ReviewStatus' in header else fields[24]
                    variant_type = fields[header.index('Type')] if 'Type' in header else fields[1]
                    name = fields[header.index('Name')] if 'Name' in header else fields[2]

                    self._variant_db[key] = {
                        'gene': gene_symbol if gene_symbol and gene_symbol != '-' else None,
                        'clinical_significance': clin_sig if clin_sig and clin_sig != '-' else None,
                        'phenotypes': phenotypes if phenotypes and phenotypes != '-' else None,
                        'rsid': f"rs{rs_id}" if rs_id and rs_id != '-' else None,
                        'review_status': review_status if review_status and review_status != '-' else None,
                        'variant_type': variant_type if variant_type and variant_type != '-' else None,
                        'hgvs': name if name and name != '-' else None,
                    }
                    count += 1

                    if progress_callback and count % 100000 == 0:
                        progress_callback(count)

                except (IndexError, ValueError) as e:
                    continue

        self._loaded = True
        logger.info(f"Loaded {count} ClinVar variants for {self.assembly}")
        return count

    def annotate(self, variant: VariantEvidence) -> VariantEvidence:
        """Annotate a single variant with ClinVar data."""
        if not self._loaded:
            self.load()

        key = f"{variant.chrom}_{variant.pos}_{variant.ref}_{variant.alt}"

        if key in self._variant_db:
            data = self._variant_db[key]
            variant.gene = data['gene']
            variant.clinical_significance = data['clinical_significance']
            variant.rsid = data['rsid']

            # Store phenotypes/disease associations
            if data['phenotypes']:
                # Clean up phenotypes - remove internal IDs
                phenotypes = data['phenotypes'].split('|')
                clean_phenotypes = []
                for p in phenotypes:
                    # Remove MONDO:, MedGen:, OMIM:, etc. prefixes from internal parts
                    if not p.startswith('MONDO:') and not p.startswith('MedGen:') and not p.startswith('OMIM:'):
                        clean_phenotypes.append(p)
                if clean_phenotypes:
                    # Set disease associations directly on the variant
                    variant.disease_associations = '; '.join(clean_phenotypes[:3])  # Limit to 3 diseases
                    variant.tags.append(f"disease:{variant.disease_associations}")

            # Set consequence/impact from the ClinVar aggregate germline
            # classification, matching exact terms (see
            # map_clinvar_significance for the primary sources). The previous
            # substring chain stamped "Conflicting classifications of
            # pathogenicity" as HIGH impact, identical to an unambiguous
            # Pathogenic call, and its `elif 'likely pathogenic'` branch was
            # unreachable because every such string also contains
            # 'pathogenic'.
            if data['clinical_significance']:
                impact, consequence, sig_tags = map_clinvar_significance(
                    data['clinical_significance']
                )
                if impact is not None:
                    variant.impact = impact
                if consequence is not None:
                    variant.consequence = data['variant_type'] or consequence
                for tag in sig_tags:
                    if tag not in variant.tags:
                        variant.tags.append(tag)

            # Add HGVS notation if available
            if data['hgvs']:
                # Parse HGVS - typically format is "NM_xxx:c.xxx (p.xxx)"
                hgvs = data['hgvs']
                if ':c.' in hgvs:
                    variant.hgvs_c = hgvs.split('(')[0].strip() if '(' in hgvs else hgvs
                if '(p.' in hgvs:
                    start = hgvs.find('(p.')
                    end = hgvs.find(')', start)
                    if start != -1 and end != -1:
                        variant.hgvs_p = hgvs[start+1:end]

        return variant

    def annotate_batch(self, variants: List[VariantEvidence]) -> List[VariantEvidence]:
        """Annotate a batch of variants."""
        if not self._loaded:
            self.load()
        return [self.annotate(v) for v in variants]

    def get_stats(self) -> dict:
        """Get statistics about loaded ClinVar data."""
        if not self._loaded:
            return {"loaded": False, "count": 0}

        # Count by clinical significance
        sig_counts = {}
        for data in self._variant_db.values():
            sig = data.get('clinical_significance', 'Unknown')
            if sig:
                # Simplify significance
                sig_simple = sig.split('/')[0].strip()
                sig_counts[sig_simple] = sig_counts.get(sig_simple, 0) + 1

        return {
            "loaded": True,
            "total_variants": len(self._variant_db),
            "by_significance": sig_counts,
        }


class AlphaMissenseAnnotator:
    """
    AlphaMissense annotation for missense variants.
    Uses Google DeepMind's AlphaMissense predictions for pathogenicity.

    AlphaMissense provides:
    - Pathogenicity score (0-1): Higher = more likely pathogenic
    - Classification: likely_benign (<0.34), ambiguous (0.34-0.564), likely_pathogenic (>0.564)
    """

    def __init__(
        self,
        alphamissense_file: Path,
    ):
        self.alphamissense_file = Path(alphamissense_file)
        self._variant_db: Dict[str, dict] = {}
        self._loaded = False

        if not self.alphamissense_file.exists():
            raise FileNotFoundError(f"AlphaMissense file not found: {self.alphamissense_file}")

    def load(self, progress_callback=None) -> int:
        """
        Load AlphaMissense database into memory for fast lookups.
        Returns number of variants loaded.

        Note: This loads ~71M variants and requires ~8-10GB RAM.
        """
        if self._loaded:
            return len(self._variant_db)

        logger.info(f"Loading AlphaMissense database from: {self.alphamissense_file}")
        logger.info("This may take a few minutes and use ~8-10GB RAM...")
        count = 0
        open_func = gzip.open if str(self.alphamissense_file).endswith('.gz') else open

        with open_func(self.alphamissense_file, 'rt') as f:
            for line_num, line in enumerate(f):
                # Skip header/comment lines
                if line.startswith('#'):
                    continue

                fields = line.strip().split('\t')
                if len(fields) < 10:
                    continue

                try:
                    # Format: CHROM, POS, REF, ALT, genome, uniprot_id, transcript_id,
                    #         protein_variant, am_pathogenicity, am_class
                    chrom = fields[0]
                    pos = fields[1]
                    ref = fields[2]
                    alt = fields[3]
                    am_pathogenicity = float(fields[8])
                    am_class = fields[9]

                    # Create lookup key (same format as ClinVar)
                    key = f"{chrom}_{pos}_{ref}_{alt}"

                    self._variant_db[key] = {
                        'am_pathogenicity': am_pathogenicity,
                        'am_class': am_class,
                        'protein_variant': fields[7] if len(fields) > 7 else None,
                        'uniprot_id': fields[5] if len(fields) > 5 else None,
                    }
                    count += 1

                    if progress_callback and count % 1000000 == 0:
                        progress_callback(count)

                except (IndexError, ValueError) as e:
                    continue

        self._loaded = True
        logger.info(f"Loaded {count:,} AlphaMissense predictions")
        return count

    def annotate(self, variant: VariantEvidence) -> VariantEvidence:
        """Annotate a single variant with AlphaMissense prediction."""
        if not self._loaded:
            self.load()

        key = f"{variant.chrom}_{variant.pos}_{variant.ref}_{variant.alt}"

        if key in self._variant_db:
            data = self._variant_db[key]
            variant.am_pathogenicity = data['am_pathogenicity']
            variant.am_class = data['am_class']

            # Add tag for high-confidence pathogenic predictions
            if data['am_pathogenicity'] >= 0.8:
                variant.tags.append("alphamissense:high_pathogenic")
            elif data['am_pathogenicity'] <= 0.2:
                variant.tags.append("alphamissense:high_benign")

        return variant

    def annotate_batch(self, variants: List[VariantEvidence]) -> List[VariantEvidence]:
        """Annotate a batch of variants."""
        if not self._loaded:
            self.load()
        return [self.annotate(v) for v in variants]

    def get_stats(self) -> dict:
        """Get statistics about loaded AlphaMissense data."""
        if not self._loaded:
            return {"loaded": False, "count": 0}

        # Count by classification
        class_counts = {"likely_benign": 0, "ambiguous": 0, "likely_pathogenic": 0}
        for data in self._variant_db.values():
            am_class = data.get('am_class', 'unknown')
            if am_class in class_counts:
                class_counts[am_class] += 1

        return {
            "loaded": True,
            "total_variants": len(self._variant_db),
            "by_class": class_counts,
        }
