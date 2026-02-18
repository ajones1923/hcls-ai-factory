"""
Variant Annotator - Add gene names, consequences, and clinical information to variants.
"""
import gzip
import subprocess
from dataclasses import dataclass
from pathlib import Path

from loguru import logger

from .vcf_parser import VariantEvidence


@dataclass
class AnnotationResult:
    """Result of variant annotation."""
    gene: str | None = None
    consequence: str | None = None
    impact: str | None = None
    hgvs_c: str | None = None
    hgvs_p: str | None = None
    clinical_significance: str | None = None
    biotype: str | None = None
    sift: str | None = None
    polyphen: str | None = None


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
        vep_cache_dir: Path | None = None,
        gene_bed_file: Path | None = None,
        species: str = "homo_sapiens",
        assembly: str = "GRCh38",
    ):
        self.use_vep = use_vep
        self.vep_cache_dir = vep_cache_dir
        self.gene_bed_file = gene_bed_file
        self.species = species
        self.assembly = assembly

        # Cache for gene lookups
        self._gene_cache: dict[str, str] = {}
        self._gene_intervals: list | None = None

        if gene_bed_file and Path(gene_bed_file).exists():
            self._load_gene_bed(gene_bed_file)

    def _load_gene_bed(self, bed_file: Path) -> None:
        """Load gene regions from BED file for fast lookup."""
        logger.info(f"Loading gene regions from: {bed_file}")
        self._gene_intervals = []

        with open(bed_file) as f:
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

    def annotate_batch(self, variants: list[VariantEvidence]) -> list[VariantEvidence]:
        """Annotate a batch of variants (more efficient for VEP)."""
        if self.use_vep and len(variants) > 1:
            return self._annotate_batch_vep(variants)
        else:
            return [self.annotate(v) for v in variants]

    def _annotate_with_vep_api(self, variant: VariantEvidence) -> AnnotationResult | None:
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

    def _annotate_batch_vep(self, variants: list[VariantEvidence]) -> list[VariantEvidence]:
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

    def _lookup_gene_bed(self, chrom: str, pos: int) -> str | None:
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
        extra_args: list[str] | None = None,
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
        except subprocess.TimeoutExpired as e:
            raise RuntimeError("VEP annotation timed out after 1 hour") from e


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
        self._variant_db: dict[str, dict] = {}
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

            # Set consequence/impact based on clinical significance
            if data['clinical_significance']:
                sig_lower = data['clinical_significance'].lower()
                if 'pathogenic' in sig_lower:
                    variant.impact = 'HIGH'
                    variant.consequence = data['variant_type'] or 'pathogenic variant'
                elif 'likely pathogenic' in sig_lower:
                    variant.impact = 'HIGH'
                    variant.consequence = data['variant_type'] or 'likely pathogenic variant'
                elif 'benign' in sig_lower:
                    variant.impact = 'LOW'
                    variant.consequence = data['variant_type'] or 'benign variant'
                elif 'uncertain' in sig_lower:
                    variant.impact = 'MODERATE'
                    variant.consequence = data['variant_type'] or 'VUS'

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

    def annotate_batch(self, variants: list[VariantEvidence]) -> list[VariantEvidence]:
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
        self._variant_db: dict[str, dict] = {}
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
            for _line_num, line in enumerate(f):
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

    def annotate_batch(self, variants: list[VariantEvidence]) -> list[VariantEvidence]:
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
