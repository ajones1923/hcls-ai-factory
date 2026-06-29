"""Tests for src/models.py — enums, collection models, and agent models.

Covers all 14 enums (valid values, str representation), all 14 collection
models (creation, to_embedding_text, field validation), PGxQuery, PGxResponse,
PGxAlert, SearchHit, CrossCollectionResult, ComparativeResult.

Author: Adam Jones
Date: March 2026
"""

from src.models import (
    MetabolizerPhenotype, TransporterFunction, HLAStatus, EnzymeDeficiency,
    GuidelineBody, CPICLevel, ClinicalAction, AlertLevel, InteractionType,
    EvidenceLevel, ReactionSeverity, DrugCategory, PGxWorkflowType,
    InhibitorStrength,
    GeneReference, DrugGuideline, DrugInteraction, HLAHypersensitivity,
    Phenoconversion, DosingAlgorithm, ClinicalEvidence, PopulationData,
    PGxClinicalTrial, FDALabel, DrugAlternative, PatientProfile,
    ImplementationProtocol, EducationMaterial,
    SearchHit, CrossCollectionResult, ComparativeResult,
    PGxAlert, PGxQuery, PGxResponse,
)


# ═══════════════════════════════════════════════════════════════════════
# ENUM TESTS
# ═══════════════════════════════════════════════════════════════════════

class TestMetabolizerPhenotype:
    def test_values(self):
        assert MetabolizerPhenotype.ULTRA_RAPID == "ultra_rapid"
        assert MetabolizerPhenotype.POOR == "poor"
        assert MetabolizerPhenotype.NORMAL == "normal"
        assert MetabolizerPhenotype.INTERMEDIATE == "intermediate"
        assert MetabolizerPhenotype.RAPID == "rapid"

    def test_count(self):
        assert len(MetabolizerPhenotype) == 5

    def test_str(self):
        assert str(MetabolizerPhenotype.POOR) == "MetabolizerPhenotype.POOR"


class TestTransporterFunction:
    def test_values(self):
        assert TransporterFunction.INCREASED == "increased"
        assert TransporterFunction.NORMAL == "normal"
        assert TransporterFunction.DECREASED == "decreased"
        assert TransporterFunction.POOR == "poor"

    def test_count(self):
        assert len(TransporterFunction) == 4


class TestHLAStatus:
    def test_values(self):
        assert HLAStatus.POSITIVE == "positive"
        assert HLAStatus.NEGATIVE == "negative"

    def test_count(self):
        assert len(HLAStatus) == 2


class TestEnzymeDeficiency:
    def test_values(self):
        assert EnzymeDeficiency.NORMAL == "normal"
        assert EnzymeDeficiency.INTERMEDIATE == "intermediate"
        assert EnzymeDeficiency.DEFICIENT == "deficient"

    def test_count(self):
        assert len(EnzymeDeficiency) == 3


class TestGuidelineBody:
    def test_values(self):
        assert GuidelineBody.CPIC == "CPIC"
        assert GuidelineBody.DPWG == "DPWG"
        assert GuidelineBody.FDA == "FDA"
        assert GuidelineBody.CPNDS == "CPNDS"

    def test_count(self):
        assert len(GuidelineBody) == 4


class TestCPICLevel:
    def test_values(self):
        assert CPICLevel.A == "A"
        assert CPICLevel.A_B == "A/B"
        assert CPICLevel.B == "B"
        assert CPICLevel.C == "C"
        assert CPICLevel.D == "D"

    def test_count(self):
        assert len(CPICLevel) == 5


class TestClinicalAction:
    def test_values(self):
        assert ClinicalAction.STANDARD == "standard"
        assert ClinicalAction.DOSE_ADJUST == "dose_adjust"
        assert ClinicalAction.ALTERNATIVE == "alternative"
        assert ClinicalAction.AVOID == "avoid"
        assert ClinicalAction.CONTRAINDICATED == "contraindicated"

    def test_count(self):
        assert len(ClinicalAction) == 5


class TestAlertLevel:
    def test_values(self):
        assert AlertLevel.INFO == "info"
        assert AlertLevel.WARNING == "warning"
        assert AlertLevel.CRITICAL == "critical"

    def test_count(self):
        assert len(AlertLevel) == 3


class TestInteractionType:
    def test_values(self):
        assert InteractionType.PK == "pharmacokinetic"
        assert InteractionType.PD == "pharmacodynamic"
        assert InteractionType.EFFICACY == "efficacy"
        assert InteractionType.TOXICITY == "toxicity"

    def test_count(self):
        assert len(InteractionType) == 4


class TestEvidenceLevel:
    def test_values(self):
        assert EvidenceLevel.LEVEL_1A == "1A"
        assert EvidenceLevel.LEVEL_2B == "2B"
        assert EvidenceLevel.LEVEL_4 == "4"

    def test_count(self):
        assert len(EvidenceLevel) == 6


class TestReactionSeverity:
    def test_values(self):
        assert ReactionSeverity.LIFE_THREATENING == "life_threatening"
        assert ReactionSeverity.SEVERE == "severe"
        assert ReactionSeverity.MODERATE == "moderate"
        assert ReactionSeverity.MILD == "mild"

    def test_count(self):
        assert len(ReactionSeverity) == 4


class TestDrugCategory:
    def test_opioids(self):
        assert DrugCategory.OPIOIDS == "opioids"

    def test_ppi(self):
        assert DrugCategory.PPI == "proton_pump_inhibitors"

    def test_count(self):
        assert len(DrugCategory) == 12


class TestPGxWorkflowType:
    def test_values(self):
        assert PGxWorkflowType.GENE_QUERY == "gene_query"
        assert PGxWorkflowType.HLA_SCREEN == "hla_screen"
        assert PGxWorkflowType.DOSING_QUERY == "dosing_query"

    def test_count(self):
        assert len(PGxWorkflowType) == 6


class TestInhibitorStrength:
    def test_values(self):
        assert InhibitorStrength.STRONG == "strong"
        assert InhibitorStrength.MODERATE == "moderate"
        assert InhibitorStrength.WEAK == "weak"

    def test_count(self):
        assert len(InhibitorStrength) == 3


# ═══════════════════════════════════════════════════════════════════════
# COLLECTION MODEL TESTS
# ═══════════════════════════════════════════════════════════════════════

class TestGeneReference:
    def test_create(self):
        g = GeneReference(id="CYP2D6_star4", gene="CYP2D6", star_allele="*4",
                          text_chunk="CYP2D6 *4 is a no-function allele")
        assert g.gene == "CYP2D6"
        assert g.star_allele == "*4"

    def test_to_embedding_text(self):
        g = GeneReference(id="test", gene="CYP2D6", star_allele="*4",
                          function_status="no function",
                          text_chunk="splicing defect",
                          defining_variants="rs3892097")
        text = g.to_embedding_text()
        assert "CYP2D6" in text
        assert "*4" in text
        assert "no function" in text
        assert "rs3892097" in text

    def test_optional_fields(self):
        g = GeneReference(id="t", gene="G", star_allele="*1", text_chunk="t")
        assert g.activity_score is None
        assert g.allele_frequency_global is None


class TestDrugGuideline:
    def test_create(self):
        dg = DrugGuideline(id="g1", gene="CYP2D6", drug="codeine",
                           phenotype="PM", recommendation="avoid",
                           text_chunk="test chunk")
        assert dg.guideline_body == GuidelineBody.CPIC
        assert dg.cpic_level == CPICLevel.A

    def test_to_embedding_text(self):
        dg = DrugGuideline(id="g1", gene="CYP2D6", drug="codeine",
                           phenotype="PM", recommendation="Avoid codeine",
                           text_chunk="chunk", dose_adjustment="reduce 50%")
        text = dg.to_embedding_text()
        assert "codeine" in text
        assert "Avoid codeine" in text
        assert "reduce 50%" in text


class TestDrugInteraction:
    def test_create_and_embed(self):
        di = DrugInteraction(id="i1", drug="codeine", gene="CYP2D6",
                             effect_description="poor metabolism",
                             text_chunk="interaction chunk",
                             clinical_significance="major")
        text = di.to_embedding_text()
        assert "codeine" in text
        assert "CYP2D6" in text
        assert "major" in text


class TestHLAHypersensitivity:
    def test_create_and_embed(self):
        h = HLAHypersensitivity(id="h1", hla_allele="HLA-B*57:01",
                                drug="abacavir", reaction_type="AHS",
                                risk_if_positive="48% risk",
                                recommendation="contraindicated",
                                text_chunk="chunk")
        text = h.to_embedding_text()
        assert "HLA-B*57:01" in text
        assert "abacavir" in text


class TestPhenoconversion:
    def test_create_and_embed(self):
        p = Phenoconversion(id="p1", affected_enzyme="CYP2D6",
                            precipitant_drug="fluoxetine",
                            effect_on_phenotype="NM to PM",
                            text_chunk="chunk",
                            affected_substrate_drugs="codeine, tramadol")
        text = p.to_embedding_text()
        assert "fluoxetine" in text
        assert "CYP2D6" in text
        assert "codeine" in text


class TestDosingAlgorithm:
    def test_create_and_embed(self):
        d = DosingAlgorithm(id="d1", drug="warfarin",
                            genes_involved="CYP2C9, VKORC1",
                            algorithm_name="IWPC",
                            text_chunk="chunk",
                            formula_description="sqrt(dose) = ...")
        text = d.to_embedding_text()
        assert "IWPC" in text
        assert "warfarin" in text
        assert "CYP2C9" in text


class TestClinicalEvidence:
    def test_create_and_embed(self):
        c = ClinicalEvidence(id="pmid123", title="PGx RCT",
                             text_chunk="outcomes data",
                             gene="CYP2D6", drug="codeine",
                             outcome_measure="pain relief",
                             outcome_value="NNT=3")
        text = c.to_embedding_text()
        assert "PGx RCT" in text
        assert "CYP2D6" in text
        assert "NNT=3" in text


class TestPopulationData:
    def test_create_and_embed(self):
        p = PopulationData(id="pop1", gene="CYP2D6", star_allele="*4",
                           population="European", allele_frequency=0.12,
                           text_chunk="chunk", source_study="1000 Genomes")
        text = p.to_embedding_text()
        assert "European" in text
        assert "0.12" in text


class TestPGxClinicalTrial:
    def test_create_and_embed(self):
        t = PGxClinicalTrial(id="nct1", title="PGx Trial",
                             text_summary="summary", gene="DPYD",
                             drug="capecitabine",
                             outcome_summary="reduced toxicity")
        text = t.to_embedding_text()
        assert "PGx Trial" in text
        assert "DPYD" in text


class TestFDALabel:
    def test_create_and_embed(self):
        f = FDALabel(id="fda1", drug="abacavir", gene="HLA-B",
                     labeling_section="Warnings", text_chunk="chunk",
                     label_type="testing required")
        text = f.to_embedding_text()
        assert "abacavir" in text
        assert "Warnings" in text
        assert "testing required" in text


class TestDrugAlternative:
    def test_create_and_embed(self):
        da = DrugAlternative(id="alt1", primary_drug="codeine",
                             gene="CYP2D6", phenotype="PM",
                             alternative_drug="morphine",
                             rationale="non-CYP2D6 pathway",
                             text_chunk="chunk")
        text = da.to_embedding_text()
        assert "codeine" in text
        assert "morphine" in text
        assert "CYP2D6" in text


class TestPatientProfile:
    def test_create_and_embed(self):
        pp = PatientProfile(id="pp1", patient_id="P001", gene="CYP2D6",
                            diplotype="*1/*4", phenotype="IM",
                            activity_score=1.0, text_chunk="chunk")
        text = pp.to_embedding_text()
        assert "*1/*4" in text
        assert "1.0" in text

    def test_no_activity_score(self):
        pp = PatientProfile(id="pp2", patient_id="P002", gene="VKORC1",
                            diplotype="*1/*1", text_chunk="chunk")
        text = pp.to_embedding_text()
        assert "Activity score" not in text


class TestImplementationProtocol:
    def test_create_and_embed(self):
        ip = ImplementationProtocol(id="imp1", institution="VUMC",
                                    title="PREDICT", text_chunk="chunk",
                                    workflow_type="preemptive",
                                    outcomes="50% dose changes")
        text = ip.to_embedding_text()
        assert "VUMC" in text
        assert "preemptive" in text


class TestEducationMaterial:
    def test_create_and_embed(self):
        em = EducationMaterial(id="edu1", title="PGx 101",
                               text_chunk="chunk", audience="clinician",
                               topic="CYP2D6 basics")
        text = em.to_embedding_text()
        assert "PGx 101" in text
        assert "clinician" in text


# ═══════════════════════════════════════════════════════════════════════
# SEARCH RESULT MODELS
# ═══════════════════════════════════════════════════════════════════════

class TestSearchHit:
    def test_create(self):
        hit = SearchHit(collection="Evidence", id="123", score=0.95,
                        text="evidence text")
        assert hit.score == 0.95
        assert hit.metadata == {}


class TestCrossCollectionResult:
    def test_hit_count(self, sample_search_hits):
        ccr = CrossCollectionResult(query="test", hits=sample_search_hits)
        assert ccr.hit_count == 3

    def test_hits_by_collection(self, sample_search_hits):
        ccr = CrossCollectionResult(query="test", hits=sample_search_hits)
        grouped = ccr.hits_by_collection()
        assert "GeneReference" in grouped
        assert "DrugGuideline" in grouped
        assert "Evidence" in grouped

    def test_empty(self):
        ccr = CrossCollectionResult(query="empty")
        assert ccr.hit_count == 0
        assert ccr.hits_by_collection() == {}


class TestComparativeResult:
    def test_total_hits(self, sample_search_hits):
        ea = CrossCollectionResult(query="a", hits=sample_search_hits[:2])
        eb = CrossCollectionResult(query="b", hits=sample_search_hits[2:])
        comp = ComparativeResult(query="compare", entity_a="CYP2D6",
                                 entity_b="CYP2C19", evidence_a=ea,
                                 evidence_b=eb)
        assert comp.total_hits == 3


# ═══════════════════════════════════════════════════════════════════════
# AGENT MODELS
# ═══════════════════════════════════════════════════════════════════════

class TestPGxAlert:
    def test_create(self):
        alert = PGxAlert(alert_level=AlertLevel.CRITICAL, gene="CYP2D6",
                         drug="codeine", phenotype="PM",
                         recommendation="avoid codeine")
        assert alert.alert_level == AlertLevel.CRITICAL
        assert alert.gene == "CYP2D6"


class TestPGxQuery:
    def test_create_minimal(self):
        q = PGxQuery(question="Is codeine safe?")
        assert q.patient_id is None
        assert q.medication_list is None

    def test_create_full(self):
        q = PGxQuery(question="test", patient_id="P001",
                     medication_list=["codeine", "fluoxetine"],
                     gene_filter="CYP2D6", drug_filter="codeine")
        assert len(q.medication_list) == 2


class TestPGxResponse:
    def test_create(self, sample_evidence, sample_alerts):
        r = PGxResponse(question="test", answer="response",
                        evidence=sample_evidence, alerts=sample_alerts,
                        knowledge_used=["CYP2D6"])
        assert r.question == "test"
        assert len(r.alerts) == 2
        assert r.timestamp  # auto-generated
