"""
Reference cases for the ClinVar germline-classification mapping.

Every expected value below is a term ClinVar itself reports in the
`CLNSIG` field of its release VCF
(https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz,
##fileDate=2026-11-01, "Aggregate germline classification for this single
variant"). `Pathogenic`, `Likely_pathogenic`, `Pathogenic/Likely_pathogenic`
and `Conflicting_classifications_of_pathogenicity` are four separate aggregate
terms in that release; the conflicting term travels with
`CLNREVSTAT=criteria_provided,_conflicting_classifications` and a `CLNSIGCONF`
field enumerating the submissions that disagree.
"""
import pytest

from src.annotator import map_clinvar_significance


class TestConflictingIsNotPathogenic:
    """A conflicting aggregate must not be reported as HIGH impact."""

    def test_conflicting_is_not_high_impact(self):
        impact, _, tags = map_clinvar_significance(
            "Conflicting classifications of pathogenicity"
        )
        assert impact == 'MODERATE'
        assert 'clinvar:conflicting' in tags

    def test_conflicting_matches_vus_tier_not_pathogenic_tier(self):
        conflicting, _, _ = map_clinvar_significance(
            "Conflicting classifications of pathogenicity"
        )
        vus, _, _ = map_clinvar_significance("Uncertain significance")
        pathogenic, _, _ = map_clinvar_significance("Pathogenic")
        assert conflicting == vus
        assert conflicting != pathogenic

    def test_conflicting_in_vcf_underscore_spelling(self):
        """The VCF spells aggregate terms with underscores; both must map."""
        impact, _, tags = map_clinvar_significance(
            "Conflicting_classifications_of_pathogenicity"
        )
        assert impact == 'MODERATE'
        assert 'clinvar:conflicting' in tags

    def test_conflicting_with_secondary_term(self):
        """`Conflicting_classifications_of_pathogenicity|risk_factor` is real data."""
        impact, _, tags = map_clinvar_significance(
            "Conflicting_classifications_of_pathogenicity|risk_factor"
        )
        assert impact == 'MODERATE'
        assert 'clinvar:conflicting' in tags


class TestPathogenicTiersAreDistinguishable:
    """The Likely-pathogenic branch used to be unreachable dead code."""

    def test_likely_pathogenic_is_reachable_and_labelled(self):
        impact, consequence, tags = map_clinvar_significance("Likely pathogenic")
        assert impact == 'HIGH'
        assert consequence == 'likely pathogenic variant'
        assert 'clinvar:likely_pathogenic' in tags

    def test_pathogenic_and_likely_pathogenic_are_separable(self):
        _, p_consequence, p_tags = map_clinvar_significance("Pathogenic")
        _, lp_consequence, lp_tags = map_clinvar_significance("Likely pathogenic")
        assert p_consequence != lp_consequence
        assert p_tags != lp_tags

    def test_combined_aggregate_term(self):
        impact, consequence, tags = map_clinvar_significance(
            "Pathogenic/Likely pathogenic"
        )
        assert impact == 'HIGH'
        assert 'clinvar:pathogenic_likely_pathogenic' in tags
        assert consequence == 'pathogenic/likely pathogenic variant'


class TestBenignAndUncertain:

    @pytest.mark.parametrize("term,expected", [
        ("Benign", 'LOW'),
        ("Likely benign", 'LOW'),
        ("Benign/Likely benign", 'LOW'),
        ("Uncertain significance", 'MODERATE'),
        ("Uncertain risk allele", 'MODERATE'),
    ])
    def test_tiers(self, term, expected):
        impact, _, _ = map_clinvar_significance(term)
        assert impact == expected

    def test_likely_benign_is_not_swallowed_by_benign_substring(self):
        _, consequence, tags = map_clinvar_significance("Likely benign")
        assert consequence == 'likely benign variant'
        assert 'clinvar:likely_benign' in tags


class TestNoTierAsserted:
    """Terms with no severity tier must not be guessed into one."""

    @pytest.mark.parametrize("term", [
        "drug response",
        "risk factor",
        "association",
        "Affects",
        "protective",
        "other",
        "not provided",
        "no classification for the single variant",
        "Likely risk allele",
        "Established risk allele",
    ])
    def test_recognised_without_impact(self, term):
        impact, _, _ = map_clinvar_significance(term)
        assert impact is None

    def test_unrecognised_term_yields_nothing(self):
        impact, consequence, tags = map_clinvar_significance("Tier I - Strong")
        assert impact is None
        assert consequence is None
        assert tags == []

    @pytest.mark.parametrize("value", [None, "", "-"])
    def test_empty_values(self, value):
        assert map_clinvar_significance(value) == (None, None, [])


class TestMultiValueMostSevereWins:

    @pytest.mark.parametrize("term,expected", [
        ("Pathogenic|other", 'HIGH'),
        ("Pathogenic|risk_factor", 'HIGH'),
        ("Benign/Likely_benign|drug_response", 'LOW'),
        ("Likely_pathogenic|association", 'HIGH'),
        ("Uncertain_significance|drug_response", 'MODERATE'),
        ("Pathogenic;Benign", 'HIGH'),
    ])
    def test_most_severe_recognised_term_wins(self, term, expected):
        impact, _, _ = map_clinvar_significance(term)
        assert impact == expected


class TestAnnotatorAppliesMapping:
    """End-to-end through ClinVarAnnotator.annotate, with a stubbed variant DB."""

    def _annotator(self, tmp_path, clin_sig, variant_type='single nucleotide variant'):
        from src.annotator import ClinVarAnnotator
        clinvar_file = tmp_path / 'variant_summary.txt.gz'
        clinvar_file.write_bytes(b'')
        ann = ClinVarAnnotator(clinvar_file)
        ann._loaded = True
        ann._variant_db = {
            'chr9_123456_A_G': {
                'gene': 'VCP',
                'clinical_significance': clin_sig,
                'phenotypes': None,
                'rsid': 'rs12345',
                'review_status': 'criteria provided, conflicting classifications',
                'variant_type': variant_type,
                'hgvs': None,
            }
        }
        return ann

    def test_conflicting_variant_is_not_high(self, tmp_path, sample_variant):
        ann = self._annotator(tmp_path, 'Conflicting classifications of pathogenicity')
        sample_variant.impact = None
        out = ann.annotate(sample_variant)
        assert out.impact == 'MODERATE'
        assert 'clinvar:conflicting' in out.tags

    def test_pathogenic_variant_stays_high(self, tmp_path, sample_variant):
        ann = self._annotator(tmp_path, 'Pathogenic')
        sample_variant.impact = None
        out = ann.annotate(sample_variant)
        assert out.impact == 'HIGH'
        assert 'clinvar:pathogenic' in out.tags

    def test_unrecognised_term_leaves_impact_untouched(self, tmp_path, sample_variant):
        ann = self._annotator(tmp_path, 'Tier II - Potential')
        sample_variant.impact = 'MODIFIER'
        out = ann.annotate(sample_variant)
        assert out.impact == 'MODIFIER'
        assert not any(t.startswith('clinvar:') for t in out.tags)
