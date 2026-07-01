"""Tests for domain knowledge base in src/knowledge.py.

Covers domain counts, drug data, gene data, clinical scale data.

Author: Adam Jones
Date: March 2026
"""


from src.knowledge import (
    CLINICAL_SCALES,
    KNOWLEDGE_VERSION,
    NEURO_DOMAINS,
    NEURO_DRUGS,
    NEURO_GENES,
    get_domain_count,
    get_drug_count,
    get_drugs_by_domain,
    get_gene_count,
    get_genes_by_domain,
    get_scale_count,
)


class TestKnowledgeVersion:
    """Tests for KNOWLEDGE_VERSION metadata."""

    def test_version_exists(self):
        assert KNOWLEDGE_VERSION["version"]

    def test_has_counts(self):
        counts = KNOWLEDGE_VERSION["counts"]
        assert counts["disease_domains"] >= 8
        assert counts["clinical_scales"] == 10

    def test_has_sources(self):
        assert len(KNOWLEDGE_VERSION["sources"]) > 0


class TestNeuroDomains:
    """Tests for NEURO_DOMAINS data."""

    def test_domain_count(self):
        assert get_domain_count() >= 10

    def test_expected_domains(self):
        expected = [
            "cerebrovascular", "degenerative", "epilepsy", "movement",
            "ms", "headache", "neuromuscular", "neuro_oncology",
        ]
        for domain in expected:
            assert domain in NEURO_DOMAINS, f"Missing domain: {domain}"

    def test_domains_have_descriptions(self):
        for name, domain in NEURO_DOMAINS.items():
            assert domain["description"], f"Domain {name} missing description"

    def test_domains_have_key_conditions(self):
        for name, domain in NEURO_DOMAINS.items():
            assert len(domain["key_conditions"]) > 0, f"Domain {name} has no key conditions"

    def test_cerebrovascular_has_stroke(self):
        conditions = NEURO_DOMAINS["cerebrovascular"]["key_conditions"]
        stroke_conditions = [c for c in conditions if "stroke" in c.lower()]
        assert len(stroke_conditions) > 0


class TestNeuroDrugs:
    """Tests for NEURO_DRUGS data."""

    def test_drug_count(self):
        assert get_drug_count() >= 40

    def test_drugs_have_required_fields(self):
        for drug in NEURO_DRUGS:
            assert drug["name"], "Drug missing name"
            assert drug["class"], f"Drug {drug['name']} missing class"
            assert drug["domain"], f"Drug {drug['name']} missing domain"
            assert drug["mechanism"], f"Drug {drug['name']} missing mechanism"

    def test_drugs_by_domain(self):
        stroke_drugs = get_drugs_by_domain("cerebrovascular")
        assert len(stroke_drugs) > 0
        assert any(d["name"] == "alteplase" for d in stroke_drugs)

    def test_epilepsy_drugs(self):
        epilepsy_drugs = get_drugs_by_domain("epilepsy")
        assert len(epilepsy_drugs) > 0
        drug_names = [d["name"] for d in epilepsy_drugs]
        assert "levetiracetam" in drug_names

    def test_ms_drugs(self):
        ms_drugs = get_drugs_by_domain("ms")
        assert len(ms_drugs) > 0
        drug_names = [d["name"] for d in ms_drugs]
        assert "ocrelizumab" in drug_names

    def test_parkinson_drugs(self):
        pd_drugs = get_drugs_by_domain("movement")
        assert len(pd_drugs) > 0
        drug_names = [d["name"] for d in pd_drugs]
        assert "levodopa" in drug_names

    def test_headache_drugs(self):
        headache_drugs = get_drugs_by_domain("headache")
        assert len(headache_drugs) > 0
        drug_names = [d["name"] for d in headache_drugs]
        assert "erenumab" in drug_names

    def test_unique_drug_names(self):
        names = [d["name"] for d in NEURO_DRUGS]
        assert len(names) == len(set(names))


class TestNeuroGenes:
    """Tests for NEURO_GENES data."""

    def test_gene_count(self):
        assert get_gene_count() >= 35

    def test_genes_have_required_fields(self):
        for gene in NEURO_GENES:
            assert gene["gene"], "Gene missing symbol"
            assert gene["domain"], f"Gene {gene['gene']} missing domain"
            assert gene["condition"], f"Gene {gene['gene']} missing condition"

    def test_genes_by_domain(self):
        epilepsy_genes = get_genes_by_domain("epilepsy")
        assert len(epilepsy_genes) > 0
        gene_names = [g["gene"] for g in epilepsy_genes]
        assert "SCN1A" in gene_names

    def test_movement_genes(self):
        movement_genes = get_genes_by_domain("movement")
        gene_names = [g["gene"] for g in movement_genes]
        assert "LRRK2" in gene_names
        assert "GBA1" in gene_names

    def test_degenerative_genes(self):
        deg_genes = get_genes_by_domain("degenerative")
        gene_names = [g["gene"] for g in deg_genes]
        assert "APP" in gene_names
        assert "APOE" in gene_names


class TestClinicalScales:
    """Tests for CLINICAL_SCALES data."""

    def test_scale_count(self):
        assert get_scale_count() == 10

    def test_expected_scales(self):
        expected = ["nihss", "gcs", "moca", "updrs", "edss", "mrs", "hit6", "alsfrs", "aspects", "hoehn_yahr"]
        for scale in expected:
            assert scale in CLINICAL_SCALES, f"Missing scale: {scale}"

    def test_scales_have_full_name(self):
        for name, scale in CLINICAL_SCALES.items():
            assert scale["full_name"], f"Scale {name} missing full_name"

    def test_scales_have_max_score(self):
        for name, scale in CLINICAL_SCALES.items():
            assert "max_score" in scale, f"Scale {name} missing max_score"

    def test_scales_have_interpretation(self):
        for name, scale in CLINICAL_SCALES.items():
            assert len(scale["interpretation"]) > 0, f"Scale {name} has no interpretations"

    def test_nihss_max_score(self):
        assert CLINICAL_SCALES["nihss"]["max_score"] == 42

    def test_gcs_max_score(self):
        assert CLINICAL_SCALES["gcs"]["max_score"] == 15

    def test_moca_max_score(self):
        assert CLINICAL_SCALES["moca"]["max_score"] == 30

    def test_edss_max_score(self):
        assert CLINICAL_SCALES["edss"]["max_score"] == 10.0
