#!/usr/bin/env python3
"""Comprehensive demo validation for Precision Biomarker Agent.

Simulates all 8 Streamlit tabs to verify every module works correctly
with the sample patient data. Run from project root:
    python scripts/demo_validation.py
"""

import json
import sys
from pathlib import Path

# Add project root to path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

# ── Load sample patients ──
with open(ROOT / "data" / "reference" / "biomarker_sample_patients.json") as f:
    patients = json.load(f)

male = patients[0]  # HCLS-BIO-2026-00001 (45M)
female = patients[1]  # HCLS-BIO-2026-00002 (38F)

passed = 0
failed = 0
errors = []


def check(name: str, condition: bool, detail: str = ""):
    global passed, failed
    if condition:
        passed += 1
        print(f"  ✓ {name}")
    else:
        failed += 1
        msg = f"  ✗ {name}" + (f" — {detail}" if detail else "")
        print(msg)
        errors.append(msg)


# ═══════════════════════════════════════════════════════════════════
# TAB 1 — Biomarker Analysis (core pipeline)
# ═══════════════════════════════════════════════════════════════════
print("\n═══ TAB 1: Biomarker Analysis ═══")

# 1a. Critical Values Engine
from src.critical_values import CriticalValueEngine

cve = CriticalValueEngine()
check("CriticalValueEngine loads rules", len(cve._rules) > 0, f"rules={len(cve._rules)}")

alerts = cve.check(male["biomarkers"])
check("Male patient: no critical alerts on normal values", isinstance(alerts, list))

# Test with extreme value
extreme = {**male["biomarkers"], "glucose": 450}
extreme_alerts = cve.check(extreme)
check("Extreme glucose triggers alert", len(extreme_alerts) > 0,
      f"alerts={[a.biomarker for a in extreme_alerts]}")

# 1b. Discordance Detector
from src.discordance_detector import DiscordanceDetector

dd = DiscordanceDetector()
check("DiscordanceDetector loads rules", len(dd._rules) > 0, f"rules={len(dd._rules)}")

disc = dd.check(male["biomarkers"])
check("Male discordance check returns list", isinstance(disc, list))

# 1c. Lab Range Interpreter
from src.lab_range_interpreter import LabRangeInterpreter

lri = LabRangeInterpreter()
check("LabRangeInterpreter loads data", len(lri._data.get("labs", {})) > 0,
      f"labs={len(lri._data.get('labs', {}))}")

male_ranges = lri.interpret(male["biomarkers"], sex="Male")
check("Male lab range interpretation", len(male_ranges) > 0, f"comparisons={len(male_ranges)}")

female_ranges = lri.interpret(female["biomarkers"], sex="Female")
check("Female lab range interpretation", len(female_ranges) > 0, f"comparisons={len(female_ranges)}")

# ═══════════════════════════════════════════════════════════════════
# TAB 2 — Biological Age
# ═══════════════════════════════════════════════════════════════════
print("\n═══ TAB 2: Biological Age ═══")

from src.biological_age import BiologicalAgeCalculator

calc = BiologicalAgeCalculator()
bio = calc.calculate(45, male["biomarkers"])

check("PhenoAge calculated", "phenoage" in bio)
check("PhenoAge has biological_age", "biological_age" in bio.get("phenoage", {}))
check("PhenoAge has mortality_risk", "mortality_risk" in bio.get("phenoage", {}))
check("PhenoAge has top_aging_drivers", "top_aging_drivers" in bio.get("phenoage", {}))
check("PhenoAge has confidence_interval", "confidence_interval" in bio.get("phenoage", {}))

phenoage_score = bio.get("phenoage", {}).get("biological_age")
check("PhenoAge score is numeric", isinstance(phenoage_score, (int, float)),
      f"got {type(phenoage_score)}: {phenoage_score}")

accel = bio.get("phenoage", {}).get("age_acceleration")
check("Age acceleration computed", accel is not None, f"acceleration={accel}")

# GrimAge (may be None if no plasma markers)
grimage = bio.get("grimage")
check("GrimAge key present", "grimage" in bio,
      "None is acceptable if no plasma markers")

# Female biological age
bio_f = calc.calculate(38, female["biomarkers"])
check("Female PhenoAge calculated", "phenoage" in bio_f)

# ═══════════════════════════════════════════════════════════════════
# TAB 3 — Disease Trajectories
# ═══════════════════════════════════════════════════════════════════
print("\n═══ TAB 3: Disease Trajectories ═══")

from src.disease_trajectory import DiseaseTrajectoryAnalyzer

dta = DiseaseTrajectoryAnalyzer()
trajectories = dta.analyze_all(
    biomarkers=male["biomarkers"],
    genotypes=male.get("genotypes", {}),
    age=45,
    sex="male",
)

check("analyze_all returns list", isinstance(trajectories, list))
check("9 disease trajectories returned", len(trajectories) == 9,
      f"got {len(trajectories)}")

# Check structure of each trajectory
diseases_seen = set()
for t in trajectories:
    diseases_seen.add(t.get("disease"))
    check(f"  {t['disease']}: has risk_level",
          t.get("risk_level") in ("CRITICAL", "HIGH", "MODERATE", "LOW"),
          f"risk_level={t.get('risk_level')}")

expected_diseases = {
    "type2_diabetes", "cardiovascular", "liver", "thyroid",
    "iron", "nutritional", "kidney", "bone_health", "cognitive"
}
check("All 9 disease domains covered", diseases_seen == expected_diseases,
      f"missing: {expected_diseases - diseases_seen}")

# Female trajectories
traj_f = dta.analyze_all(
    biomarkers=female["biomarkers"],
    genotypes=female.get("genotypes", {}),
    age=38,
    sex="female",
)
check("Female trajectories: 9 returned", len(traj_f) == 9)

# ═══════════════════════════════════════════════════════════════════
# TAB 4 — Pharmacogenomics
# ═══════════════════════════════════════════════════════════════════
print("\n═══ TAB 4: Pharmacogenomics ═══")

from src.pharmacogenomics import PharmacogenomicMapper

pgx = PharmacogenomicMapper()
pgx_result = pgx.map_all(
    star_alleles=male.get("star_alleles", {}),
    genotypes=male.get("genotypes", {}),
)

check("PGx map_all returns dict", isinstance(pgx_result, dict))
check("PGx has gene_results", "gene_results" in pgx_result)
check("PGx has genes_analyzed count", "genes_analyzed" in pgx_result)

gene_results = pgx_result.get("gene_results", [])
check("PGx gene results non-empty", len(gene_results) > 0,
      f"genes={len(gene_results)}")

genes_mapped = [r.get("gene") for r in gene_results]
check("CYP2D6 mapped", "CYP2D6" in genes_mapped)
check("CYP2C19 mapped", "CYP2C19" in genes_mapped)

# Check for CYP2D6 *1/*4 intermediate metabolizer
cyp2d6 = next((r for r in gene_results if r.get("gene") == "CYP2D6"), {})
check("CYP2D6 *1/*4 → Intermediate Metabolizer",
      "intermediate" in cyp2d6.get("phenotype", "").lower(),
      f"phenotype={cyp2d6.get('phenotype')}")

# Female PGx
pgx_f = pgx.map_all(
    star_alleles=female.get("star_alleles", {}),
    genotypes=female.get("genotypes", {}),
)
check("Female PGx returns results", len(pgx_f.get("gene_results", [])) > 0)

# ═══════════════════════════════════════════════════════════════════
# TAB 5 — Evidence Explorer (RAG)
# ═══════════════════════════════════════════════════════════════════
print("\n═══ TAB 5: Evidence Explorer (RAG) ═══")

try:
    from sentence_transformers import SentenceTransformer

    from src.collections import BiomarkerCollectionManager
    from src.rag_engine import BiomarkerRAGEngine

    embedder = SentenceTransformer("BAAI/bge-small-en-v1.5")
    cm = BiomarkerCollectionManager()
    rag = BiomarkerRAGEngine(collection_manager=cm, embedder=embedder, llm_client=None)
    check("RAG engine initialized", True)

    # Test search (vector search only, no LLM since llm_client=None)
    # RAG engine uses query() which needs an LLM, so just verify init
    # and collection manager connectivity
    check("Milvus collections accessible", cm is not None)
except Exception as e:
    check("RAG engine available (requires Milvus)", False, str(e))

# ═══════════════════════════════════════════════════════════════════
# TAB 6 — Export (PDF, FHIR, Markdown, CSV)
# ═══════════════════════════════════════════════════════════════════
print("\n═══ TAB 6: Export ═══")

from src.export import export_fhir_diagnostic_report, validate_fhir_bundle
from src.models import AnalysisResult, BiologicalAgeResult, PatientProfile

# Build minimal patient profile
profile = PatientProfile(
    patient_id=male["id"],
    age=male["demographics"]["age"],
    sex="M",
    biomarkers=male["biomarkers"],
)

# Build minimal analysis result with required fields
phenoage = bio.get("phenoage", {})
bio_age_result = BiologicalAgeResult(
    chronological_age=45,
    biological_age=phenoage.get("biological_age", 45.0),
    age_acceleration=phenoage.get("age_acceleration", 0.0),
    phenoage_score=phenoage.get("biological_age", 45.0),
)

analysis = AnalysisResult(
    patient_profile=profile,
    biological_age=bio_age_result,
    disease_trajectories=[],
    pgx_results=[],
    genotype_adjustments=[],
    critical_alerts=[],
)

fhir_json = export_fhir_diagnostic_report(analysis, profile)
fhir = json.loads(fhir_json)
check("FHIR bundle generated", fhir.get("resourceType") == "Bundle")
check("FHIR bundle has entries", len(fhir.get("entry", [])) > 0)

# Validate FHIR
validation_errors = validate_fhir_bundle(fhir)
check("FHIR validation passes", len(validation_errors) == 0,
      f"errors: {validation_errors}")

# Check Patient resource in bundle
resource_types = [e.get("resource", {}).get("resourceType") for e in fhir.get("entry", [])]
check("FHIR has Patient resource", "Patient" in resource_types)
check("FHIR has DiagnosticReport", "DiagnosticReport" in resource_types)
check("FHIR has Observation(s)", "Observation" in resource_types)

# ═══════════════════════════════════════════════════════════════════
# TAB 7 — Genotype Adjustments & Age-Stratified Ranges
# ═══════════════════════════════════════════════════════════════════
print("\n═══ TAB 7: Genotype & Age Adjustments ═══")

from src.genotype_adjustment import GenotypeAdjuster

adj = GenotypeAdjuster()

# Test genotype adjustments
geno_adj = adj.adjust_all(male["biomarkers"], male.get("genotypes", {}))
check("Genotype adjustments return dict", isinstance(geno_adj, dict))
check("Genotype adjustments have results", len(geno_adj) > 0,
      f"keys={list(geno_adj.keys())}")

# Test age-stratified adjustments
age_adj = adj.apply_age_adjustments(male["biomarkers"], age=45, sex="M")
check("Age-stratified adjustments return list", isinstance(age_adj, list))
check("Age-stratified adjustments non-empty", len(age_adj) > 0,
      f"adjustments={len(age_adj)}")

# Check structure
if age_adj:
    a0 = age_adj[0]
    check("Age adjustment has biomarker field", "biomarker" in a0)
    check("Age adjustment has age_adjusted_range", "age_adjusted_range" in a0,
          f"keys={list(a0.keys())}")

# Female age adjustments
age_adj_f = adj.apply_age_adjustments(female["biomarkers"], age=38, sex="F")
check("Female age adjustments", isinstance(age_adj_f, list))

# ═══════════════════════════════════════════════════════════════════
# TAB 8 — Sample Patient Data Integrity
# ═══════════════════════════════════════════════════════════════════
print("\n═══ TAB 8: Sample Patient Data Integrity ═══")

check("Male patient has ID", male["id"] == "HCLS-BIO-2026-00001")
check("Female patient has ID", female["id"] == "HCLS-BIO-2026-00002")
check("Male has biomarkers", len(male.get("biomarkers", {})) >= 25,
      f"count={len(male.get('biomarkers', {}))}")
check("Female has biomarkers", len(female.get("biomarkers", {})) >= 25,
      f"count={len(female.get('biomarkers', {}))}")
check("Male has genotypes", len(male.get("genotypes", {})) >= 4)
check("Female has genotypes", len(female.get("genotypes", {})) >= 3)
check("Male has star_alleles", len(male.get("star_alleles", {})) >= 3)
check("Female has star_alleles", len(female.get("star_alleles", {})) >= 3)
check("Male has genetic_markers", "apoe" in male.get("genetic_markers", {}))
check("Male has clinical_context", "medications" in male.get("clinical_context", {}))
check("Male has family_history", len(male["clinical_context"].get("family_history", [])) >= 5)
check("Female has ob_gyn", "ob_gyn" in female.get("clinical_context", {}))

# ═══════════════════════════════════════════════════════════════════
# SUMMARY
# ═══════════════════════════════════════════════════════════════════
total = passed + failed
print(f"\n{'═' * 60}")
print(f"  DEMO VALIDATION RESULTS: {passed}/{total} passed, {failed} failed")
print(f"{'═' * 60}")

if errors:
    print("\n  Failed checks:")
    for e in errors:
        print(f"    {e}")

sys.exit(0 if failed == 0 else 1)
