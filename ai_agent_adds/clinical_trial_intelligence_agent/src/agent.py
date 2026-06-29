"""Clinical Trial Intelligence Agent -- autonomous reasoning across clinical trial data.

Implements the plan -> search -> evaluate -> synthesize -> report pattern from the
VAST AI OS AgentEngine model. The agent can:

1. Parse complex multi-part questions about clinical trial design and strategy
2. Plan a search strategy across 14 domain-specific collections
3. Execute multi-collection retrieval via the TrialRAGEngine
4. Evaluate evidence quality and completeness
5. Synthesize cross-functional insights with regulatory alerts
6. Generate structured reports with trial-specific formatting

Mapping to VAST AI OS:
  - AgentEngine entry point: TrialIntelligenceAgent.run()
  - Plan -> search_plan()
  - Execute -> rag_engine.query()
  - Reflect -> evaluate_evidence()
  - Report -> generate_report()

Author: Adam Jones
Date: March 2026
"""

from dataclasses import dataclass, field
from datetime import datetime, timezone
from enum import Enum
from typing import Dict, List, Optional

from src.models import TrialWorkflowType


class EvidenceLevel(str, Enum):
    """Clinical evidence hierarchy for trial recommendations."""
    LEVEL_1A = "1a"       # Systematic review of RCTs
    LEVEL_1B = "1b"       # Individual RCT
    LEVEL_2A = "2a"       # Systematic review of cohort studies
    LEVEL_2B = "2b"       # Individual cohort study
    LEVEL_3 = "3"         # Case-control study
    LEVEL_4 = "4"         # Case series
    LEVEL_5 = "5"         # Expert opinion
    REGULATORY = "reg"    # Regulatory guidance (FDA/EMA/ICH)


class TrialPhase(str, Enum):
    """Clinical trial phases."""
    EARLY_PHASE_1 = "early_phase_1"
    PHASE_1 = "phase_1"
    PHASE_1_2 = "phase_1_2"
    PHASE_2 = "phase_2"
    PHASE_2_3 = "phase_2_3"
    PHASE_3 = "phase_3"
    PHASE_4 = "phase_4"


class RegulatoryBody(str, Enum):
    """Regulatory bodies referenced in trial intelligence."""
    FDA = "fda"
    EMA = "ema"
    PMDA = "pmda"
    ICH = "ich"
    WHO = "who"


class SeverityLevel(str, Enum):
    """Finding severity classification."""
    CRITICAL = "critical"
    HIGH = "high"
    MODERATE = "moderate"
    LOW = "low"
    INFORMATIONAL = "informational"


# =====================================================================
# RESPONSE DATACLASS
# =====================================================================

@dataclass
class TrialResponse:
    """Complete response from the trial intelligence agent.

    Attributes:
        question: Original query text.
        answer: LLM-synthesised answer text.
        results: Ranked search results used for synthesis.
        workflow: Clinical trial workflow that was applied.
        confidence: Overall confidence score (0.0 - 1.0).
        citations: Formatted citation list.
        search_time_ms: Total search time in milliseconds.
        collections_searched: Number of collections queried.
        patient_context_used: Whether patient context was injected.
        timestamp: ISO 8601 timestamp of response generation.
    """
    question: str = ""
    answer: str = ""
    results: list = field(default_factory=list)
    workflow: Optional[TrialWorkflowType] = None
    confidence: float = 0.0
    citations: List[Dict[str, str]] = field(default_factory=list)
    search_time_ms: float = 0.0
    collections_searched: int = 0
    patient_context_used: bool = False
    timestamp: str = field(
        default_factory=lambda: datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    )


# =====================================================================
# CLINICAL TRIAL SYSTEM PROMPT
# =====================================================================

TRIAL_SYSTEM_PROMPT = """\
You are a clinical trial intelligence system embedded within the HCLS AI Factory \
precision medicine platform. You have deep expertise in clinical trial design, \
patient recruitment, regulatory strategy, and competitive analysis across all \
therapeutic areas. You provide evidence-based recommendations for protocol \
optimization, patient-trial matching, site selection, eligibility criteria \
analysis, adaptive trial design, safety signal detection, regulatory document \
preparation, and competitive intelligence. Always cite specific trials \
(NCT numbers), regulatory guidance (ICH, FDA), and evidence levels. Never \
fabricate trial data or regulatory precedents.

Your responses must adhere to the following standards:

1. **Trial Citations** -- Always cite clinical trials using NCT identifiers with \
   clickable links: [NCT01234567](https://clinicaltrials.gov/study/NCT01234567). \
   Include trial phase, sponsor, status (recruiting, completed, terminated), and \
   primary endpoint when available.

2. **Regulatory References** -- Cite FDA guidance documents, ICH guidelines \
   (E6(R3), E8(R1), E9(R1), E17, E20), EMA scientific advice, and PMDA guidance \
   with document identifiers. Format regulatory references as: \
   [ICH E9(R1)] Estimands and Sensitivity Analysis in Clinical Trials.

3. **CRITICAL Findings** -- Flag the following as CRITICAL with prominent visual \
   markers and immediate action recommendations:
   - Safety signals requiring DSMB notification or trial hold
   - Protocol deviations affecting data integrity
   - Regulatory deficiency letters or clinical hold orders
   - Enrollment below 50% of target at midpoint
   - Significant imbalance in baseline characteristics
   - Unblinding events or protocol amendments affecting primary endpoint
   - GCP violations or site-level compliance failures
   - Serious unexpected adverse reactions (SUSARs)
   - Futility analysis triggering early stopping criteria
   - Data integrity concerns from site monitoring

4. **Severity Badges** -- Classify all findings using standardised severity levels: \
   [CRITICAL], [HIGH], [MODERATE], [LOW], [INFORMATIONAL]. Place the badge at the \
   start of each finding or recommendation line.

5. **Statistical Context** -- Include relevant statistical parameters: sample size \
   calculations, power analysis, effect sizes (Cohen's d, hazard ratios, odds ratios), \
   confidence intervals, p-values with appropriate multiplicity adjustments (Bonferroni, \
   Hochberg, Holm). Reference the estimand framework (ICH E9(R1)) when discussing \
   treatment effects.

6. **Structured Formatting** -- Organise responses with clear sections: \
   Trial Summary, Design Rationale, Regulatory Considerations, Statistical Plan, \
   Risk Assessment, Recommendations, and Timeline. Use bullet points and numbered \
   lists for actionable items.

7. **Genomic Cross-Reference** -- When biomarker or genetic stratification is \
   relevant (e.g., companion diagnostics, enrichment designs, umbrella/basket \
   trials, precision oncology), cross-reference with the genomic_evidence \
   collection and recommend appropriate biomarker validation strategies.

8. **Regulatory Frameworks** -- Reference and interpret validated regulatory \
   frameworks: 21 CFR Parts 11/50/56/312, ICH E6(R3) GCP, ICH E8(R1) general \
   considerations, ICH E9(R1) estimands, ICH E17 multi-regional trials, \
   ICH E20 adaptive designs, FDA Breakthrough/Fast Track/Accelerated Approval \
   pathways, EMA PRIME designation, Orphan Drug designations.

9. **Competitive Analysis** -- For competitive intelligence queries, provide \
   structured comparison across: mechanism of action, trial design, enrollment \
   targets, primary/secondary endpoints, interim results, projected timelines, \
   and regulatory strategy differentiation.

10. **Limitations** -- You are a clinical trial intelligence support tool. You \
    do NOT replace regulatory affairs professionals, biostatisticians, or clinical \
    development leadership. All recommendations require review by qualified \
    professionals with knowledge of the specific therapeutic context, regulatory \
    jurisdiction, and sponsor development strategy. Explicitly state when evidence \
    is limited or when specialist consultation is recommended."""


# =====================================================================
# WORKFLOW-SPECIFIC COLLECTION BOOST WEIGHTS
# =====================================================================
# Maps each TrialWorkflowType to collection weight overrides (multipliers).
# Collections not listed retain their base weight (1.0x). Values > 1.0
# boost the collection; values < 1.0 would suppress it.

WORKFLOW_COLLECTION_BOOST: Dict[TrialWorkflowType, Dict[str, float]] = {

    # ── Protocol Design ──────────────────────────────────────────────
    TrialWorkflowType.PROTOCOL_DESIGN: {
        "trial_protocols": 2.5,
        "trial_endpoints": 2.0,
        "trial_regulatory": 1.5,
        "trial_guidelines": 1.5,
        "trial_adaptive": 1.3,
        "trial_literature": 1.2,
        "trial_results": 1.1,
        "trial_biomarkers": 1.0,
    },

    # ── Patient-Trial Matching ───────────────────────────────────────
    TrialWorkflowType.PATIENT_MATCHING: {
        "trial_eligibility": 2.5,
        "trial_protocols": 1.8,
        "trial_biomarkers": 1.5,
        "trial_sites": 1.3,
        "trial_safety": 1.2,
        "genomic_evidence": 1.5,
        "trial_guidelines": 1.0,
        "trial_literature": 1.0,
    },

    # ── Site Selection ───────────────────────────────────────────────
    TrialWorkflowType.SITE_SELECTION: {
        "trial_sites": 2.5,
        "trial_investigators": 2.0,
        "trial_protocols": 1.3,
        "trial_results": 1.2,
        "trial_regulatory": 1.1,
        "trial_literature": 1.0,
        "trial_rwe": 1.0,
    },

    # ── Eligibility Criteria Analysis ────────────────────────────────
    TrialWorkflowType.ELIGIBILITY_ANALYSIS: {
        "trial_eligibility": 2.5,
        "trial_protocols": 2.0,
        "trial_guidelines": 1.5,
        "trial_biomarkers": 1.3,
        "trial_safety": 1.2,
        "trial_literature": 1.2,
        "trial_rwe": 1.1,
        "genomic_evidence": 1.0,
    },

    # ── Endpoint Strategy ────────────────────────────────────────────
    TrialWorkflowType.ENDPOINT_STRATEGY: {
        "trial_endpoints": 2.5,
        "trial_results": 2.0,
        "trial_protocols": 1.5,
        "trial_regulatory": 1.5,
        "trial_guidelines": 1.3,
        "trial_biomarkers": 1.2,
        "trial_literature": 1.1,
        "trial_adaptive": 1.0,
    },

    # ── Regulatory Strategy ──────────────────────────────────────────
    TrialWorkflowType.REGULATORY_STRATEGY: {
        "trial_regulatory": 2.5,
        "trial_guidelines": 2.0,
        "trial_protocols": 1.5,
        "trial_endpoints": 1.3,
        "trial_results": 1.3,
        "trial_safety": 1.2,
        "trial_literature": 1.1,
        "trial_adaptive": 1.0,
    },

    # ── Competitive Intelligence ─────────────────────────────────────
    TrialWorkflowType.COMPETITIVE_INTELLIGENCE: {
        "trial_results": 2.5,
        "trial_protocols": 2.0,
        "trial_endpoints": 1.8,
        "trial_literature": 1.5,
        "trial_regulatory": 1.3,
        "trial_sites": 1.1,
        "trial_biomarkers": 1.0,
        "trial_safety": 1.0,
    },

    # ── Safety Monitoring ────────────────────────────────────────────
    TrialWorkflowType.SAFETY_MONITORING: {
        "trial_safety": 2.5,
        "trial_regulatory": 1.8,
        "trial_protocols": 1.5,
        "trial_results": 1.3,
        "trial_guidelines": 1.3,
        "trial_biomarkers": 1.2,
        "trial_literature": 1.1,
        "trial_endpoints": 1.0,
    },

    # ── Adaptive Design ─────────────────────────────────────────────
    TrialWorkflowType.ADAPTIVE_DESIGN: {
        "trial_adaptive": 2.5,
        "trial_endpoints": 2.0,
        "trial_protocols": 1.8,
        "trial_regulatory": 1.5,
        "trial_guidelines": 1.3,
        "trial_results": 1.2,
        "trial_literature": 1.1,
        "trial_biomarkers": 1.0,
    },

    # ── Biomarker Strategy ───────────────────────────────────────────
    TrialWorkflowType.BIOMARKER_STRATEGY: {
        "trial_biomarkers": 2.5,
        "genomic_evidence": 2.0,
        "trial_eligibility": 1.5,
        "trial_endpoints": 1.5,
        "trial_protocols": 1.3,
        "trial_literature": 1.2,
        "trial_regulatory": 1.1,
        "trial_results": 1.0,
    },

    # ── Real-World Evidence Analysis ─────────────────────────────────
    TrialWorkflowType.RWE_ANALYSIS: {
        "trial_rwe": 2.5,
        "trial_results": 2.0,
        "trial_literature": 1.5,
        "trial_regulatory": 1.3,
        "trial_endpoints": 1.2,
        "trial_protocols": 1.1,
        "trial_safety": 1.0,
        "trial_guidelines": 1.0,
    },

    # ── Recruitment Optimization ─────────────────────────────────────
    TrialWorkflowType.RECRUITMENT_OPTIMIZATION: {
        "trial_sites": 2.5,
        "trial_eligibility": 2.0,
        "trial_investigators": 1.8,
        "trial_rwe": 1.5,
        "trial_protocols": 1.3,
        "trial_literature": 1.0,
        "trial_biomarkers": 1.0,
    },

    # ── General (balanced across all collections) ────────────────────
    TrialWorkflowType.GENERAL: {
        "trial_protocols": 1.2,
        "trial_eligibility": 1.1,
        "trial_endpoints": 1.1,
        "trial_sites": 1.0,
        "trial_investigators": 0.9,
        "trial_results": 1.2,
        "trial_regulatory": 1.1,
        "trial_literature": 1.2,
        "trial_biomarkers": 1.0,
        "trial_safety": 1.1,
        "trial_rwe": 1.0,
        "trial_adaptive": 0.9,
        "trial_guidelines": 1.1,
        "genomic_evidence": 0.8,
    },
}


# =====================================================================
# KNOWLEDGE DOMAIN DICTIONARIES
# =====================================================================
# Comprehensive clinical knowledge for entity detection and context
# enrichment. Used by the agent's search_plan() to identify entities
# in user queries and map them to workflows.

TRIAL_CONDITIONS: Dict[str, Dict[str, object]] = {
    # ── Oncology ─────────────────────────────────────────────────────
    "non-small cell lung cancer": {
        "aliases": ["nsclc", "lung cancer", "lung adenocarcinoma",
                    "lung squamous cell carcinoma", "non-small cell"],
        "workflows": [TrialWorkflowType.PROTOCOL_DESIGN, TrialWorkflowType.BIOMARKER_STRATEGY],
        "search_terms": ["EGFR", "ALK", "PD-L1", "KRAS G12C", "immunotherapy lung"],
    },
    "breast cancer": {
        "aliases": ["bc", "her2+ breast cancer", "triple negative breast cancer",
                    "tnbc", "hr+ breast cancer", "er+ breast cancer"],
        "workflows": [TrialWorkflowType.PROTOCOL_DESIGN, TrialWorkflowType.BIOMARKER_STRATEGY],
        "search_terms": ["HER2", "ER", "PR", "CDK4/6", "PARP inhibitor"],
    },
    "colorectal cancer": {
        "aliases": ["crc", "colon cancer", "rectal cancer",
                    "metastatic colorectal cancer", "mcrc"],
        "workflows": [TrialWorkflowType.PROTOCOL_DESIGN, TrialWorkflowType.ENDPOINT_STRATEGY],
        "search_terms": ["MSI-H", "KRAS", "BRAF V600E", "cetuximab"],
    },
    "acute myeloid leukemia": {
        "aliases": ["aml", "acute myelogenous leukemia"],
        "workflows": [TrialWorkflowType.PROTOCOL_DESIGN, TrialWorkflowType.BIOMARKER_STRATEGY],
        "search_terms": ["FLT3", "IDH1", "IDH2", "NPM1", "venetoclax"],
    },
    "multiple myeloma": {
        "aliases": ["mm", "myeloma", "plasma cell myeloma"],
        "workflows": [TrialWorkflowType.PROTOCOL_DESIGN, TrialWorkflowType.COMPETITIVE_INTELLIGENCE],
        "search_terms": ["CAR-T BCMA", "bispecific", "daratumumab"],
    },

    # ── Cardiovascular ───────────────────────────────────────────────
    "heart failure": {
        "aliases": ["hf", "chf", "hfref", "hfpef", "congestive heart failure"],
        "workflows": [TrialWorkflowType.PROTOCOL_DESIGN, TrialWorkflowType.ENDPOINT_STRATEGY],
        "search_terms": ["LVEF", "NT-proBNP", "SGLT2i", "ARNI"],
    },
    "atrial fibrillation": {
        "aliases": ["afib", "af", "atrial fib"],
        "workflows": [TrialWorkflowType.PROTOCOL_DESIGN, TrialWorkflowType.SAFETY_MONITORING],
        "search_terms": ["anticoagulation", "ablation", "DOAC", "stroke prevention"],
    },

    # ── Neurology ────────────────────────────────────────────────────
    "alzheimer disease": {
        "aliases": ["alzheimers", "ad", "alzheimer's disease",
                    "alzheimer's", "mild cognitive impairment", "mci"],
        "workflows": [TrialWorkflowType.PROTOCOL_DESIGN, TrialWorkflowType.ENDPOINT_STRATEGY],
        "search_terms": ["amyloid", "tau", "ADAS-Cog", "CDR-SB", "lecanemab"],
    },
    "parkinson disease": {
        "aliases": ["parkinsons", "pd", "parkinson's disease", "parkinson's"],
        "workflows": [TrialWorkflowType.PROTOCOL_DESIGN, TrialWorkflowType.BIOMARKER_STRATEGY],
        "search_terms": ["alpha-synuclein", "UPDRS", "levodopa", "dopamine"],
    },

    # ── Immunology / Autoimmune ──────────────────────────────────────
    "rheumatoid arthritis": {
        "aliases": ["ra", "rheumatoid"],
        "workflows": [TrialWorkflowType.PROTOCOL_DESIGN, TrialWorkflowType.COMPETITIVE_INTELLIGENCE],
        "search_terms": ["TNF inhibitor", "JAK inhibitor", "ACR20", "DAS28"],
    },
    "systemic lupus erythematosus": {
        "aliases": ["sle", "lupus"],
        "workflows": [TrialWorkflowType.PROTOCOL_DESIGN, TrialWorkflowType.ELIGIBILITY_ANALYSIS],
        "search_terms": ["anti-dsDNA", "complement", "SRI-4", "SLEDAI"],
    },
    "crohn disease": {
        "aliases": ["crohns", "crohn's disease", "crohn's"],
        "workflows": [TrialWorkflowType.PROTOCOL_DESIGN, TrialWorkflowType.ENDPOINT_STRATEGY],
        "search_terms": ["anti-TNF", "vedolizumab", "ustekinumab", "CDAI"],
    },
    "ulcerative colitis": {
        "aliases": ["uc", "colitis"],
        "workflows": [TrialWorkflowType.PROTOCOL_DESIGN, TrialWorkflowType.ENDPOINT_STRATEGY],
        "search_terms": ["Mayo score", "mucosal healing", "JAK inhibitor"],
    },

    # ── Metabolic ────────────────────────────────────────────────────
    "type 2 diabetes": {
        "aliases": ["t2dm", "t2d", "type 2 diabetes mellitus", "diabetes mellitus"],
        "workflows": [TrialWorkflowType.PROTOCOL_DESIGN, TrialWorkflowType.ENDPOINT_STRATEGY],
        "search_terms": ["HbA1c", "GLP-1", "SGLT2", "insulin", "CVOT"],
    },
    "obesity": {
        "aliases": ["obese", "overweight", "weight management", "metabolic syndrome"],
        "workflows": [TrialWorkflowType.PROTOCOL_DESIGN, TrialWorkflowType.COMPETITIVE_INTELLIGENCE],
        "search_terms": ["GLP-1 RA", "semaglutide", "tirzepatide", "BMI"],
    },
    "nash": {
        "aliases": ["nonalcoholic steatohepatitis", "non-alcoholic steatohepatitis",
                    "mash", "metabolic dysfunction-associated steatohepatitis",
                    "fatty liver disease", "nafld"],
        "workflows": [TrialWorkflowType.PROTOCOL_DESIGN, TrialWorkflowType.ENDPOINT_STRATEGY],
        "search_terms": ["fibrosis", "NAS score", "liver biopsy", "FIB-4"],
    },

    # ── Infectious Disease ───────────────────────────────────────────
    "hiv": {
        "aliases": ["human immunodeficiency virus", "hiv-1", "hiv/aids"],
        "workflows": [TrialWorkflowType.PROTOCOL_DESIGN, TrialWorkflowType.REGULATORY_STRATEGY],
        "search_terms": ["antiretroviral", "viral load", "CD4", "broadly neutralizing antibody"],
    },
    "respiratory syncytial virus": {
        "aliases": ["rsv", "rsv infection"],
        "workflows": [TrialWorkflowType.PROTOCOL_DESIGN, TrialWorkflowType.COMPETITIVE_INTELLIGENCE],
        "search_terms": ["nirsevimab", "vaccine", "prophylaxis", "infant"],
    },

    # ── Rare Disease ─────────────────────────────────────────────────
    "spinal muscular atrophy": {
        "aliases": ["sma", "sma type 1", "sma type 2"],
        "workflows": [TrialWorkflowType.REGULATORY_STRATEGY, TrialWorkflowType.ELIGIBILITY_ANALYSIS],
        "search_terms": ["SMN1", "SMN2", "nusinersen", "risdiplam", "gene therapy"],
    },
    "duchenne muscular dystrophy": {
        "aliases": ["dmd", "duchenne"],
        "workflows": [TrialWorkflowType.REGULATORY_STRATEGY, TrialWorkflowType.PROTOCOL_DESIGN],
        "search_terms": ["dystrophin", "exon skipping", "gene therapy", "6MWT"],
    },
    "cystic fibrosis": {
        "aliases": ["cf"],
        "workflows": [TrialWorkflowType.PROTOCOL_DESIGN, TrialWorkflowType.BIOMARKER_STRATEGY],
        "search_terms": ["CFTR", "trikafta", "elexacaftor", "FEV1", "sweat chloride"],
    },

    # ── Dermatology ──────────────────────────────────────────────────
    "atopic dermatitis": {
        "aliases": ["ad", "eczema", "atopic eczema"],
        "workflows": [TrialWorkflowType.PROTOCOL_DESIGN, TrialWorkflowType.COMPETITIVE_INTELLIGENCE],
        "search_terms": ["dupilumab", "JAK inhibitor", "EASI", "IGA"],
    },
    "psoriasis": {
        "aliases": ["plaque psoriasis", "psoriatic disease"],
        "workflows": [TrialWorkflowType.PROTOCOL_DESIGN, TrialWorkflowType.ENDPOINT_STRATEGY],
        "search_terms": ["IL-17", "IL-23", "PASI 90", "bimekizumab"],
    },

    # ── Neuro-Oncology ────────────────────────────────────────────────
    "glioblastoma": {
        "aliases": ["gbm", "glioblastoma multiforme", "idh-mutant glioma",
                    "idh-wild-type glioblastoma", "high-grade glioma"],
        "workflows": [TrialWorkflowType.PROTOCOL_DESIGN, TrialWorkflowType.BIOMARKER_STRATEGY],
        "search_terms": ["MGMT methylation", "temozolomide", "radiation", "IDH mutation",
                        "bevacizumab", "TTFields", "tumor treating fields"],
    },

    # ── Hepatology / GI Oncology ──────────────────────────────────────
    "hepatocellular carcinoma": {
        "aliases": ["hcc", "liver cancer", "hepatoma",
                    "child-pugh", "bclc staging"],
        "workflows": [TrialWorkflowType.PROTOCOL_DESIGN, TrialWorkflowType.COMPETITIVE_INTELLIGENCE],
        "search_terms": ["AFP", "sorafenib", "lenvatinib", "atezolizumab + bevacizumab",
                        "TACE", "Child-Pugh classification", "BCLC"],
    },
    "pancreatic cancer": {
        "aliases": ["pdac", "pancreatic ductal adenocarcinoma",
                    "pancreatic adenocarcinoma", "exocrine pancreatic cancer"],
        "workflows": [TrialWorkflowType.PROTOCOL_DESIGN, TrialWorkflowType.ENDPOINT_STRATEGY],
        "search_terms": ["CA 19-9", "FOLFIRINOX", "gemcitabine", "nab-paclitaxel",
                        "KRAS", "stromal targeting", "pancreatic stellate cells"],
    },

    # ── Gynecologic Oncology ──────────────────────────────────────────
    "ovarian cancer": {
        "aliases": ["ovarian carcinoma", "epithelial ovarian cancer",
                    "high-grade serous ovarian cancer", "hgsoc",
                    "brca-mutant ovarian cancer"],
        "workflows": [TrialWorkflowType.PROTOCOL_DESIGN, TrialWorkflowType.BIOMARKER_STRATEGY],
        "search_terms": ["BRCA mutation", "HRD", "PARP inhibitor", "platinum sensitivity",
                        "olaparib", "niraparib", "bevacizumab maintenance"],
    },

    # ── Urology ───────────────────────────────────────────────────────
    "bladder cancer": {
        "aliases": ["urothelial carcinoma", "urothelial cancer",
                    "transitional cell carcinoma", "muscle-invasive bladder cancer",
                    "non-muscle-invasive bladder cancer", "mibc", "nmibc"],
        "workflows": [TrialWorkflowType.PROTOCOL_DESIGN, TrialWorkflowType.BIOMARKER_STRATEGY],
        "search_terms": ["PD-L1", "enfortumab vedotin", "erdafitinib", "FGFR",
                        "avelumab maintenance", "BCG", "Nectin-4"],
    },

    # ── Hematology (non-malignant) ────────────────────────────────────
    "myelodysplastic syndrome": {
        "aliases": ["mds", "myelodysplastic syndromes", "myelodysplasia",
                    "refractory anemia", "refractory cytopenia"],
        "workflows": [TrialWorkflowType.PROTOCOL_DESIGN, TrialWorkflowType.BIOMARKER_STRATEGY],
        "search_terms": ["IPSS-R", "luspatercept", "azacitidine", "decitabine",
                        "imetelstat", "SF3B1", "TP53", "ring sideroblasts"],
    },

    # ── Pulmonology ───────────────────────────────────────────────────
    "idiopathic pulmonary fibrosis": {
        "aliases": ["ipf", "pulmonary fibrosis", "usual interstitial pneumonia",
                    "uip pattern"],
        "workflows": [TrialWorkflowType.PROTOCOL_DESIGN, TrialWorkflowType.ENDPOINT_STRATEGY],
        "search_terms": ["FVC", "nintedanib", "pirfenidone", "HRCT", "honeycombing",
                        "GAP index", "anti-fibrotic"],
    },

    # ── Nephrology ────────────────────────────────────────────────────
    "chronic kidney disease": {
        "aliases": ["ckd", "chronic renal disease", "renal insufficiency",
                    "diabetic kidney disease", "dkd", "ckd stage 3-5"],
        "workflows": [TrialWorkflowType.PROTOCOL_DESIGN, TrialWorkflowType.ENDPOINT_STRATEGY],
        "search_terms": ["eGFR", "SGLT2 inhibitor", "finerenone", "albuminuria",
                        "UACR", "kidney failure", "dialysis"],
    },

    # ── Neurology (headache) ──────────────────────────────────────────
    "migraine": {
        "aliases": ["episodic migraine", "chronic migraine",
                    "migraine with aura", "migraine without aura"],
        "workflows": [TrialWorkflowType.PROTOCOL_DESIGN, TrialWorkflowType.COMPETITIVE_INTELLIGENCE],
        "search_terms": ["CGRP", "erenumab", "galcanezumab", "fremanezumab",
                        "rimegepant", "monthly migraine days", "acute treatment"],
    },

    # ── Endocrinology ─────────────────────────────────────────────────
    "type 1 diabetes": {
        "aliases": ["t1dm", "t1d", "type 1 diabetes mellitus",
                    "juvenile diabetes", "insulin-dependent diabetes"],
        "workflows": [TrialWorkflowType.PROTOCOL_DESIGN, TrialWorkflowType.REGULATORY_STRATEGY],
        "search_terms": ["teplizumab", "islet cell autoimmunity", "C-peptide",
                        "anti-CD3", "GAD65 antibody", "beta-cell preservation"],
    },

    # ── Hepatology (metabolic) ────────────────────────────────────────
    "nash_mafld": {
        "aliases": ["nash", "mafld", "metabolic associated fatty liver disease",
                    "nonalcoholic steatohepatitis", "steatohepatitis",
                    "mash", "masld"],
        "workflows": [TrialWorkflowType.PROTOCOL_DESIGN, TrialWorkflowType.ENDPOINT_STRATEGY],
        "search_terms": ["NAS score", "fibrosis stage", "resmetirom", "liver biopsy",
                        "FIB-4", "ELF score", "THR-beta agonist"],
    },

    # ── Rheumatology ──────────────────────────────────────────────────
    "systemic sclerosis": {
        "aliases": ["ssc", "scleroderma", "systemic scleroderma",
                    "diffuse cutaneous systemic sclerosis", "limited cutaneous ssc"],
        "workflows": [TrialWorkflowType.PROTOCOL_DESIGN, TrialWorkflowType.ELIGIBILITY_ANALYSIS],
        "search_terms": ["mRSS", "modified Rodnan skin score", "nintedanib for ILD",
                        "tocilizumab", "interstitial lung disease", "Scl-70", "anti-centromere"],
    },
}


TRIAL_DRUGS: Dict[str, Dict[str, object]] = {
    # ── Oncology ─────────────────────────────────────────────────────
    "pembrolizumab": {
        "aliases": ["keytruda"],
        "mechanism": "Anti-PD-1 monoclonal antibody",
        "indications": ["NSCLC", "melanoma", "head and neck", "MSI-H solid tumors"],
        "workflows": ["protocol_design", "competitive_intelligence"],
    },
    "nivolumab": {
        "aliases": ["opdivo"],
        "mechanism": "Anti-PD-1 monoclonal antibody",
        "indications": ["NSCLC", "melanoma", "renal cell carcinoma", "hepatocellular carcinoma"],
        "workflows": ["protocol_design", "competitive_intelligence"],
    },
    "trastuzumab": {
        "aliases": ["herceptin"],
        "mechanism": "Anti-HER2 monoclonal antibody",
        "indications": ["HER2+ breast cancer", "HER2+ gastric cancer"],
        "workflows": ["protocol_design", "biomarker_strategy"],
    },
    "trastuzumab deruxtecan": {
        "aliases": ["t-dxd", "enhertu"],
        "mechanism": "Anti-HER2 antibody-drug conjugate",
        "indications": ["HER2+ breast cancer", "HER2-low breast cancer", "NSCLC", "gastric"],
        "workflows": ["protocol_design", "competitive_intelligence"],
    },
    "osimertinib": {
        "aliases": ["tagrisso"],
        "mechanism": "Third-generation EGFR TKI",
        "indications": ["EGFR-mutant NSCLC"],
        "workflows": ["protocol_design", "biomarker_strategy"],
    },
    "sotorasib": {
        "aliases": ["lumakras"],
        "mechanism": "KRAS G12C inhibitor",
        "indications": ["KRAS G12C NSCLC"],
        "workflows": ["protocol_design", "biomarker_strategy"],
    },
    "venetoclax": {
        "aliases": ["venclexta"],
        "mechanism": "BCL-2 inhibitor",
        "indications": ["CLL", "AML"],
        "workflows": ["protocol_design", "competitive_intelligence"],
    },

    # ── Metabolic / Cardio ───────────────────────────────────────────
    "semaglutide": {
        "aliases": ["ozempic", "wegovy", "rybelsus"],
        "mechanism": "GLP-1 receptor agonist",
        "indications": ["T2DM", "obesity", "NASH", "cardiovascular risk"],
        "workflows": ["protocol_design", "competitive_intelligence"],
    },
    "tirzepatide": {
        "aliases": ["mounjaro", "zepbound"],
        "mechanism": "Dual GIP/GLP-1 receptor agonist",
        "indications": ["T2DM", "obesity"],
        "workflows": ["protocol_design", "competitive_intelligence"],
    },
    "empagliflozin": {
        "aliases": ["jardiance"],
        "mechanism": "SGLT2 inhibitor",
        "indications": ["T2DM", "heart failure", "CKD"],
        "workflows": ["protocol_design", "endpoint_strategy"],
    },
    "dapagliflozin": {
        "aliases": ["farxiga"],
        "mechanism": "SGLT2 inhibitor",
        "indications": ["T2DM", "heart failure", "CKD"],
        "workflows": ["protocol_design", "endpoint_strategy"],
    },

    # ── Immunology ───────────────────────────────────────────────────
    "dupilumab": {
        "aliases": ["dupixent"],
        "mechanism": "Anti-IL-4/IL-13 monoclonal antibody",
        "indications": ["atopic dermatitis", "asthma", "CRSwNP"],
        "workflows": ["protocol_design", "competitive_intelligence"],
    },
    "upadacitinib": {
        "aliases": ["rinvoq"],
        "mechanism": "Selective JAK1 inhibitor",
        "indications": ["RA", "atopic dermatitis", "ulcerative colitis", "Crohn's"],
        "workflows": ["protocol_design", "competitive_intelligence"],
    },
    "risankizumab": {
        "aliases": ["skyrizi"],
        "mechanism": "Anti-IL-23 p19 monoclonal antibody",
        "indications": ["psoriasis", "Crohn's disease", "ulcerative colitis"],
        "workflows": ["protocol_design", "competitive_intelligence"],
    },
    "secukinumab": {
        "aliases": ["cosentyx"],
        "mechanism": "Anti-IL-17A monoclonal antibody",
        "indications": ["psoriasis", "psoriatic arthritis", "ankylosing spondylitis"],
        "workflows": ["protocol_design", "competitive_intelligence"],
    },

    # ── Neurology ────────────────────────────────────────────────────
    "lecanemab": {
        "aliases": ["leqembi"],
        "mechanism": "Anti-amyloid beta monoclonal antibody",
        "indications": ["Alzheimer's disease early stage"],
        "workflows": ["protocol_design", "safety_monitoring"],
    },
    "donanemab": {
        "aliases": ["kisunla"],
        "mechanism": "Anti-amyloid beta monoclonal antibody (N3pG)",
        "indications": ["Alzheimer's disease early stage"],
        "workflows": ["protocol_design", "competitive_intelligence"],
    },

    # ── Gene / Cell Therapy ──────────────────────────────────────────
    "axicabtagene ciloleucel": {
        "aliases": ["yescarta", "axi-cel"],
        "mechanism": "Anti-CD19 CAR-T cell therapy",
        "indications": ["DLBCL", "follicular lymphoma"],
        "workflows": ["protocol_design", "safety_monitoring"],
    },
    "tisagenlecleucel": {
        "aliases": ["kymriah", "tisa-cel"],
        "mechanism": "Anti-CD19 CAR-T cell therapy",
        "indications": ["B-ALL", "DLBCL"],
        "workflows": ["protocol_design", "safety_monitoring"],
    },
    "onasemnogene abeparvovec": {
        "aliases": ["zolgensma"],
        "mechanism": "AAV9-based SMN1 gene therapy",
        "indications": ["SMA type 1"],
        "workflows": ["protocol_design", "regulatory_strategy"],
    },

    # ── Infectious Disease ───────────────────────────────────────────
    "lenacapavir": {
        "aliases": ["sunlenca"],
        "mechanism": "Capsid inhibitor",
        "indications": ["HIV-1 infection", "HIV PrEP"],
        "workflows": ["protocol_design", "regulatory_strategy"],
    },
    "nirmatrelvir/ritonavir": {
        "aliases": ["paxlovid"],
        "mechanism": "SARS-CoV-2 3CLpro inhibitor + CYP3A4 inhibitor",
        "indications": ["COVID-19"],
        "workflows": ["protocol_design", "regulatory_strategy"],
    },

    # ── ADCs (Antibody-Drug Conjugates) ───────────────────────────────
    "enfortumab vedotin": {
        "aliases": ["padcev", "ev"],
        "mechanism": "Anti-nectin-4 antibody-drug conjugate",
        "indications": ["urothelial carcinoma", "bladder cancer"],
        "workflows": ["protocol_design", "competitive_intelligence"],
    },
    "sacituzumab govitecan": {
        "aliases": ["trodelvy"],
        "mechanism": "Anti-Trop-2 antibody-drug conjugate",
        "indications": ["TNBC", "urothelial carcinoma", "HR+/HER2- breast cancer"],
        "workflows": ["protocol_design", "competitive_intelligence"],
    },
    "datopotamab deruxtecan": {
        "aliases": ["dato-dxd", "ds-1062"],
        "mechanism": "Anti-Trop-2 antibody-drug conjugate",
        "indications": ["NSCLC", "HR+/HER2- breast cancer"],
        "workflows": ["protocol_design", "competitive_intelligence"],
    },

    # ── Bispecifics / TCEs ────────────────────────────────────────────
    "elranatamab": {
        "aliases": ["elrexfio"],
        "mechanism": "BCMA x CD3 bispecific T-cell engager",
        "indications": ["multiple myeloma"],
        "workflows": ["protocol_design", "safety_monitoring"],
    },

    # ── NASH / Metabolic ──────────────────────────────────────────────
    "resmetirom": {
        "aliases": ["rezdiffra", "mgl-3196"],
        "mechanism": "Thyroid hormone receptor-beta (THR-beta) agonist",
        "indications": ["NASH", "MASH", "MAFLD"],
        "workflows": ["protocol_design", "regulatory_strategy"],
    },

    # ── Immunology / Autoimmune ───────────────────────────────────────
    "teplizumab": {
        "aliases": ["tzield"],
        "mechanism": "Anti-CD3 monoclonal antibody",
        "indications": ["type 1 diabetes delay of onset", "T1DM stage 2"],
        "workflows": ["protocol_design", "regulatory_strategy"],
    },

    # ── Hematology ────────────────────────────────────────────────────
    "luspatercept": {
        "aliases": ["reblozyl"],
        "mechanism": "TGF-beta superfamily ligand trap (modified activin receptor IIB-Fc)",
        "indications": ["MDS with ring sideroblasts", "beta-thalassemia"],
        "workflows": ["protocol_design", "competitive_intelligence"],
    },
    "imetelstat": {
        "aliases": ["imetelstat sodium"],
        "mechanism": "Telomerase inhibitor (oligonucleotide)",
        "indications": ["lower-risk MDS", "myelofibrosis"],
        "workflows": ["protocol_design", "biomarker_strategy"],
    },

    # ── Targeted Oncology ─────────────────────────────────────────────
    "capivasertib": {
        "aliases": ["truqap", "azd5363"],
        "mechanism": "AKT inhibitor (pan-AKT kinase inhibitor)",
        "indications": ["HR+/HER2- breast cancer with PIK3CA/AKT1/PTEN alteration"],
        "workflows": ["protocol_design", "biomarker_strategy"],
    },
    "zolbetuximab": {
        "aliases": ["zolbetuximab-clzb", "imab362"],
        "mechanism": "Anti-CLDN18.2 monoclonal antibody",
        "indications": ["CLDN18.2-positive gastric cancer", "GEJ adenocarcinoma"],
        "workflows": ["protocol_design", "biomarker_strategy"],
    },

    # ── Complement Pathway ────────────────────────────────────────────
    "crovalimab": {
        "aliases": ["piasky", "sky59"],
        "mechanism": "Anti-complement C5 recycling antibody (long-acting)",
        "indications": ["paroxysmal nocturnal hemoglobinuria", "PNH"],
        "workflows": ["protocol_design", "competitive_intelligence"],
    },
    "iptacopan": {
        "aliases": ["fabhalta", "lnp023"],
        "mechanism": "Complement factor B inhibitor (oral small molecule)",
        "indications": ["PNH", "C3 glomerulopathy", "IgA nephropathy"],
        "workflows": ["protocol_design", "competitive_intelligence"],
    },
}


TRIAL_BIOMARKERS: Dict[str, Dict[str, str]] = {
    "pd-l1": {
        "full_name": "Programmed Death-Ligand 1",
        "assay": "IHC (22C3, SP263, SP142, 28-8 clones)",
        "significance": "Predictive biomarker for anti-PD-1/PD-L1 immunotherapy response",
        "workflows": "biomarker_strategy,protocol_design",
    },
    "her2": {
        "full_name": "Human Epidermal Growth Factor Receptor 2",
        "assay": "IHC + FISH/ISH",
        "significance": "Predictive biomarker for HER2-targeted therapy; enrichment stratifier",
        "workflows": "biomarker_strategy,protocol_design",
    },
    "egfr": {
        "full_name": "Epidermal Growth Factor Receptor",
        "assay": "NGS, PCR-based mutation testing",
        "significance": "Predictive biomarker for EGFR TKI therapy in NSCLC",
        "workflows": "biomarker_strategy,eligibility_analysis",
    },
    "msi-h": {
        "full_name": "Microsatellite Instability-High / dMMR",
        "assay": "IHC for MMR proteins, PCR, NGS",
        "significance": "Predictive for checkpoint inhibitor response; tumor-agnostic indication",
        "workflows": "biomarker_strategy,protocol_design",
    },
    "ctdna": {
        "full_name": "Circulating Tumor DNA",
        "assay": "Liquid biopsy (ddPCR, NGS-based)",
        "significance": "MRD detection, treatment response monitoring, resistance mechanism ID",
        "workflows": "biomarker_strategy,endpoint_strategy",
    },
    "nt-probnp": {
        "full_name": "N-terminal pro-B-type Natriuretic Peptide",
        "assay": "Immunoassay",
        "significance": "Heart failure diagnosis, prognosis, and outcome endpoint in CV trials",
        "workflows": "endpoint_strategy,safety_monitoring",
    },
    "hba1c": {
        "full_name": "Glycated Hemoglobin A1c",
        "assay": "HPLC, immunoassay",
        "significance": "Primary efficacy endpoint in diabetes trials; glycemic control marker",
        "workflows": "endpoint_strategy,eligibility_analysis",
    },
    "cdr-sb": {
        "full_name": "Clinical Dementia Rating Sum of Boxes",
        "assay": "Clinical assessment scale",
        "significance": "Primary efficacy endpoint in Alzheimer's disease clinical trials",
        "workflows": "endpoint_strategy,protocol_design",
    },
    "amyloid pet": {
        "full_name": "Amyloid PET Imaging",
        "assay": "PET with florbetapir, florbetaben, flutemetamol",
        "significance": "Enrollment biomarker and treatment response endpoint in AD trials",
        "workflows": "biomarker_strategy,eligibility_analysis",
    },
    "fev1": {
        "full_name": "Forced Expiratory Volume in 1 Second",
        "assay": "Spirometry",
        "significance": "Primary endpoint in respiratory/CF trials; lung function measure",
        "workflows": "endpoint_strategy,protocol_design",
    },
    "acr20": {
        "full_name": "ACR 20% Response Criteria",
        "assay": "Composite clinical assessment",
        "significance": "Primary endpoint in rheumatoid arthritis clinical trials",
        "workflows": "endpoint_strategy,protocol_design",
    },
    "pasi 90": {
        "full_name": "Psoriasis Area and Severity Index 90% Response",
        "assay": "Clinical assessment scale",
        "significance": "Primary/co-primary endpoint in psoriasis clinical trials",
        "workflows": "endpoint_strategy,protocol_design",
    },
    "tmb": {
        "full_name": "Tumor Mutational Burden",
        "assay": "Whole exome sequencing, targeted panels (F1CDx)",
        "significance": "Predictive biomarker for immunotherapy; FDA-approved companion diagnostic",
        "workflows": "biomarker_strategy,eligibility_analysis",
    },
    "brca1/2": {
        "full_name": "BRCA1 and BRCA2 Germline/Somatic Mutations",
        "assay": "NGS, targeted sequencing",
        "significance": "Predictive for PARP inhibitor response; enrichment biomarker",
        "workflows": "biomarker_strategy,protocol_design",
    },
    "alk": {
        "full_name": "Anaplastic Lymphoma Kinase Fusion",
        "assay": "FISH, IHC, NGS",
        "significance": "Predictive biomarker for ALK TKI therapy in NSCLC",
        "workflows": "biomarker_strategy,eligibility_analysis",
    },
    "troponin": {
        "full_name": "High-Sensitivity Cardiac Troponin",
        "assay": "hs-cTnI / hs-cTnT immunoassay",
        "significance": "Cardiac safety monitoring in oncology trials; cardiotoxicity detection",
        "workflows": "safety_monitoring,endpoint_strategy",
    },
    "fgfr": {
        "full_name": "Fibroblast Growth Factor Receptor",
        "assay": "NGS, FISH, IHC",
        "significance": "Predictive biomarker for FGFR-targeted therapy (erdafitinib); urothelial, cholangiocarcinoma",
        "workflows": "biomarker_strategy,eligibility_analysis",
    },
    "ros1": {
        "full_name": "ROS Proto-Oncogene 1 Receptor Tyrosine Kinase",
        "assay": "FISH, IHC, NGS (fusion detection)",
        "significance": "Predictive biomarker for ROS1 TKI therapy (crizotinib, entrectinib) in NSCLC",
        "workflows": "biomarker_strategy,eligibility_analysis",
    },
    "cldn18.2": {
        "full_name": "Claudin 18 Isoform 2",
        "assay": "IHC (VENTANA CLDN18 (43-14A) RxDx Assay)",
        "significance": "Predictive biomarker for zolbetuximab eligibility in gastric/GEJ cancer",
        "workflows": "biomarker_strategy,protocol_design",
    },
    "mrd": {
        "full_name": "Minimal Residual Disease",
        "assay": "Flow cytometry, NGS-based (clonoSEQ), PCR",
        "significance": "Response depth biomarker in hematologic malignancies; emerging endpoint in AML, ALL, MM",
        "workflows": "biomarker_strategy,endpoint_strategy",
    },
    "nfl": {
        "full_name": "Neurofilament Light Chain",
        "assay": "Simoa (serum/CSF), ELISA",
        "significance": "Neurodegeneration biomarker; axonal injury marker in MS, ALS, AD trials",
        "workflows": "biomarker_strategy,endpoint_strategy",
    },
    "ca_19_9": {
        "full_name": "Carbohydrate Antigen 19-9",
        "assay": "Immunoassay (serum)",
        "significance": "Pancreatic cancer monitoring; response assessment and surveillance biomarker",
        "workflows": "endpoint_strategy,protocol_design",
    },
    "afp": {
        "full_name": "Alpha-Fetoprotein",
        "assay": "Immunoassay (serum)",
        "significance": "HCC screening and monitoring; prognostic biomarker and treatment response marker",
        "workflows": "biomarker_strategy,endpoint_strategy",
    },
    "trop2": {
        "full_name": "Trophoblast Cell Surface Antigen 2 (TACSTD2)",
        "assay": "IHC",
        "significance": "ADC target biomarker; sacituzumab govitecan and datopotamab deruxtecan eligibility",
        "workflows": "biomarker_strategy,protocol_design",
    },
}


# =====================================================================
# SEARCH PLAN DATACLASS
# =====================================================================

@dataclass
class SearchPlan:
    """Agent's plan for answering a clinical trial intelligence question.

    The search plan captures all entities detected in the user's question
    and the strategy the agent will use to retrieve evidence from the
    14 trial-specific Milvus collections.
    """
    question: str
    conditions: List[str] = field(default_factory=list)
    drugs: List[str] = field(default_factory=list)
    biomarkers: List[str] = field(default_factory=list)
    relevant_workflows: List[TrialWorkflowType] = field(default_factory=list)
    search_strategy: str = "broad"  # broad, targeted, comparative, regulatory
    sub_questions: List[str] = field(default_factory=list)
    identified_topics: List[str] = field(default_factory=list)


# =====================================================================
# CLINICAL TRIAL INTELLIGENCE AGENT
# =====================================================================

class TrialIntelligenceAgent:
    """Autonomous Clinical Trial Intelligence Agent.

    Wraps the multi-collection TrialRAGEngine with planning and reasoning
    capabilities. Designed to answer complex cross-functional questions
    about clinical trial design, regulatory strategy, patient recruitment,
    and competitive intelligence.

    Example queries this agent handles:
    - "Design a Phase 3 protocol for KRAS G12C NSCLC with adaptive interim analysis"
    - "What eligibility criteria optimize enrollment for HFrEF trials?"
    - "Compare pembrolizumab vs nivolumab in first-line NSCLC trial designs"
    - "What FDA accelerated approval pathway requirements apply to our AML program?"
    - "Identify optimal sites for a decentralized atopic dermatitis trial"
    - "What safety monitoring plan is needed for a CAR-T trial?"
    - "What biomarker enrichment strategy maximizes response rate in TNBC?"
    - "Analyze competing Alzheimer's trials and differentiation opportunities"

    Usage:
        agent = TrialIntelligenceAgent(rag_engine)
        plan = agent.search_plan("Design a Phase 3 NSCLC protocol")
        response = agent.run("Design a Phase 3 NSCLC protocol")
    """

    def __init__(self, rag_engine):
        """Initialize agent with a configured RAG engine.

        Args:
            rag_engine: TrialRAGEngine instance with Milvus collections connected.
        """
        self.rag = rag_engine
        self.knowledge = {
            "conditions": TRIAL_CONDITIONS,
            "biomarkers": TRIAL_BIOMARKERS,
            "drugs": TRIAL_DRUGS,
        }

    # ── Public API ──────────────────────────────────────────────────

    def run(
        self,
        query: str,
        workflow_type: Optional[TrialWorkflowType] = None,
        patient_context: Optional[dict] = None,
        **kwargs,
    ) -> TrialResponse:
        """Execute the full agent pipeline: plan -> search -> evaluate -> synthesize.

        Args:
            query: Natural language question about clinical trials.
            workflow_type: Optional workflow override for collection boosting.
            patient_context: Optional patient data for patient-trial matching.
            **kwargs: Additional query parameters (top_k, collection_filter).

        Returns:
            TrialResponse with findings, recommendations, and metadata.
        """
        # Phase 1: Plan
        plan = self.search_plan(query)

        # Phase 2: Determine workflow (allow override)
        workflow = workflow_type or (
            plan.relevant_workflows[0] if plan.relevant_workflows else None
        )

        # Phase 3: Search via RAG engine
        top_k = kwargs.get("top_k", 5)

        response = self.rag.query(
            question=query,
            workflow=workflow,
            top_k=top_k,
            patient_context=patient_context,
        )

        # Phase 4: Evaluate and potentially expand
        if hasattr(response, "results") and response.results is not None:
            quality = self.evaluate_evidence(response.results)
            if quality == "insufficient" and plan.sub_questions:
                for sub_q in plan.sub_questions[:2]:
                    sub_response = self.rag.search(sub_q, top_k=top_k)
                    if sub_response:
                        response.results.extend(sub_response)

        return response

    def search_plan(self, question: str) -> SearchPlan:
        """Analyze a question and create an optimised search plan.

        Detects clinical conditions, drugs, and biomarkers in the question
        text. Determines relevant trial workflows, chooses a search
        strategy, and generates sub-questions for comprehensive retrieval
        across collections.

        Args:
            question: The user's natural language question.

        Returns:
            SearchPlan with all detected entities and retrieval strategy.
        """
        plan = SearchPlan(question=question)

        # Step 1: Detect entities
        entities = self._detect_entities(question)
        plan.conditions = entities.get("conditions", [])
        plan.drugs = entities.get("drugs", [])
        plan.biomarkers = entities.get("biomarkers", [])

        # Step 2: Determine relevant workflows
        plan.relevant_workflows = [self._detect_workflow(question)]
        # Add entity-derived workflows
        for condition in plan.conditions:
            info = TRIAL_CONDITIONS.get(condition, {})
            for wf in info.get("workflows", []):
                if wf not in plan.relevant_workflows:
                    plan.relevant_workflows.append(wf)

        # Step 3: Choose search strategy
        plan.search_strategy = self._choose_strategy(
            question, plan.conditions, plan.drugs,
        )

        # Step 4: Generate sub-questions
        plan.sub_questions = self._generate_sub_questions(plan)

        # Step 5: Compile identified topics
        plan.identified_topics = (
            plan.conditions + plan.drugs + plan.biomarkers
        )

        return plan

    def evaluate_evidence(self, results) -> str:
        """Evaluate the quality and coverage of retrieved evidence.

        Uses collection diversity and hit count to assess whether
        the retrieved evidence is sufficient for a comprehensive answer.

        Args:
            results: List of search results from the RAG engine.

        Returns:
            "sufficient", "partial", or "insufficient".
        """
        if not results:
            return "insufficient"

        total_hits = len(results)
        collections_seen = set()

        for result in results:
            if hasattr(result, "collection"):
                collections_seen.add(result.collection)
            elif isinstance(result, dict):
                collections_seen.add(result.get("collection", "unknown"))

        num_collections = len(collections_seen)

        if num_collections >= 3 and total_hits >= 10:
            return "sufficient"
        elif num_collections >= 2 and total_hits >= 5:
            return "partial"
        else:
            return "insufficient"

    def generate_report(self, results, workflow: Optional[TrialWorkflowType] = None) -> str:
        """Generate a structured clinical trial intelligence report.

        Args:
            results: Response object from run() or rag.query().
            workflow: Optional workflow type for section customisation.

        Returns:
            Formatted markdown report string.
        """
        timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        question = results.question if hasattr(results, "question") else ""
        plan = self.search_plan(question) if question else SearchPlan(question="")

        report_lines = [
            "# Clinical Trial Intelligence Report",
            f"**Query:** {question}",
            f"**Generated:** {timestamp}",
            f"**Workflows:** {', '.join(wf.value for wf in plan.relevant_workflows)}",
            f"**Strategy:** {plan.search_strategy}",
            "",
        ]

        # Detected entities
        if plan.conditions or plan.drugs or plan.biomarkers:
            report_lines.extend([
                "---",
                "",
                "## Detected Clinical Entities",
                "",
            ])
            if plan.conditions:
                report_lines.append(
                    f"- **Conditions/Indications:** {', '.join(plan.conditions)}"
                )
            if plan.drugs:
                report_lines.append(
                    f"- **Investigational Agents:** {', '.join(plan.drugs)}"
                )
            if plan.biomarkers:
                report_lines.append(
                    f"- **Biomarkers:** {', '.join(plan.biomarkers)}"
                )
            report_lines.append("")

        # Critical findings check
        critical_flags = []
        if hasattr(results, "results") and results.results:
            for r in results.results:
                meta = r.metadata if hasattr(r, "metadata") else {}
                if meta.get("urgency") == "critical" or meta.get("safety_signal"):
                    critical_flags.append(r)

        if critical_flags:
            report_lines.extend([
                "---",
                "",
                "## [CRITICAL] Safety / Regulatory Alerts",
                "",
            ])
            for flag in critical_flags:
                text = flag.text if hasattr(flag, "text") else str(flag)
                report_lines.append(
                    f"- **[CRITICAL]** {text[:200]} -- "
                    f"immediate review required."
                )
            report_lines.append("")

        # Analysis section
        report_lines.extend([
            "---",
            "",
            "## Analysis",
            "",
        ])

        if hasattr(results, "answer"):
            report_lines.append(results.answer)
        elif hasattr(results, "summary"):
            report_lines.append(results.summary)
        elif isinstance(results, str):
            report_lines.append(results)
        else:
            report_lines.append("No analysis generated.")

        report_lines.append("")

        # Regulatory considerations
        if workflow in (TrialWorkflowType.REGULATORY_STRATEGY,
                        TrialWorkflowType.PROTOCOL_DESIGN,
                        TrialWorkflowType.ADAPTIVE_DESIGN):
            report_lines.extend([
                "---",
                "",
                "## Regulatory Considerations",
                "",
                "- ICH E6(R3): Good Clinical Practice",
                "- ICH E8(R1): General Considerations for Clinical Studies",
                "- ICH E9(R1): Statistical Principles -- Estimands and Sensitivity Analysis",
                "- 21 CFR Part 312: Investigational New Drug Application",
                "",
            ])

        # Confidence and metadata
        confidence = results.confidence if hasattr(results, "confidence") else 0.0
        report_lines.extend([
            "---",
            "",
            "## Metadata",
            "",
            f"- **Confidence Score:** {confidence:.3f}",
            f"- **Collections Searched:** {results.collections_searched if hasattr(results, 'collections_searched') else 'N/A'}",
            f"- **Search Time:** {results.search_time_ms if hasattr(results, 'search_time_ms') else 'N/A'} ms",
            "",
            "---",
            "",
            "*This report is generated by the Clinical Trial Intelligence Agent "
            "within the HCLS AI Factory. All recommendations require review by "
            "qualified regulatory, clinical, and statistical professionals.*",
        ])

        return "\n".join(report_lines)

    # ── Workflow Detection ──────────────────────────────────────────

    def _detect_workflow(self, question: str) -> TrialWorkflowType:
        """Detect the most relevant workflow from a question.

        Uses keyword-based heuristics to identify which of the 13 trial
        workflows is most relevant to the query.

        Args:
            question: The user's natural language question.

        Returns:
            Most relevant TrialWorkflowType.
        """
        text_upper = question.upper()

        workflow_scores: Dict[TrialWorkflowType, float] = {}

        keyword_workflow_map = {
            TrialWorkflowType.PROTOCOL_DESIGN: [
                "PROTOCOL", "TRIAL DESIGN", "STUDY DESIGN", "PHASE 1",
                "PHASE 2", "PHASE 3", "PHASE 4", "RANDOMIZED", "RANDOMISED",
                "DOUBLE-BLIND", "PLACEBO-CONTROLLED", "OPEN-LABEL",
                "CROSSOVER", "PARALLEL", "SAMPLE SIZE", "POWER",
            ],
            TrialWorkflowType.PATIENT_MATCHING: [
                "PATIENT MATCHING", "PATIENT SELECTION", "PATIENT-TRIAL",
                "MATCH PATIENT", "FIND TRIAL", "ELIGIBLE PATIENT",
                "CANDIDATE PATIENT", "ENROLL PATIENT",
            ],
            TrialWorkflowType.SITE_SELECTION: [
                "SITE SELECTION", "SITE IDENTIFICATION", "INVESTIGATOR SITE",
                "CLINICAL SITE", "ENROLLMENT SITE", "SITE PERFORMANCE",
                "SITE FEASIBILITY", "DECENTRALIZED",
            ],
            TrialWorkflowType.ELIGIBILITY_ANALYSIS: [
                "ELIGIBILITY", "INCLUSION CRITERIA", "EXCLUSION CRITERIA",
                "INCLUSION/EXCLUSION", "ENROLL", "ELIGIBLE",
                "CRITERIA ANALYSIS", "SCREEN FAILURE",
            ],
            TrialWorkflowType.ENDPOINT_STRATEGY: [
                "ENDPOINT", "PRIMARY ENDPOINT", "SECONDARY ENDPOINT",
                "SURROGATE ENDPOINT", "COMPOSITE ENDPOINT", "PRO",
                "PATIENT REPORTED OUTCOME", "OVERALL SURVIVAL", "OS",
                "PROGRESSION FREE SURVIVAL", "PFS", "ORR", "DOR",
                "RESPONSE RATE",
            ],
            TrialWorkflowType.REGULATORY_STRATEGY: [
                "REGULATORY", "FDA", "EMA", "PMDA", "IND", "NDA", "BLA",
                "MAA", "BREAKTHROUGH", "ACCELERATED APPROVAL", "FAST TRACK",
                "PRIORITY REVIEW", "ORPHAN", "REMS", "ICH",
                "GUIDANCE", "ADVISORY COMMITTEE",
            ],
            TrialWorkflowType.COMPETITIVE_INTELLIGENCE: [
                "COMPETITIVE", "COMPETITOR", "LANDSCAPE", "COMPARE",
                "VS ", "VERSUS", "HEAD-TO-HEAD", "DIFFERENTIATION",
                "MARKET", "PIPELINE", "COMPETING",
            ],
            TrialWorkflowType.SAFETY_MONITORING: [
                "SAFETY", "ADVERSE EVENT", "SAE", "SUSAR", "DSMB",
                "DMC", "SAFETY SIGNAL", "DOSE LIMITING TOXICITY",
                "DLT", "MAXIMUM TOLERATED DOSE", "MTD", "CRS",
                "CYTOKINE RELEASE", "NEUROTOXICITY", "ICANS",
            ],
            TrialWorkflowType.ADAPTIVE_DESIGN: [
                "ADAPTIVE", "BAYESIAN", "INTERIM ANALYSIS", "FUTILITY",
                "DOSE ESCALATION", "BASKET", "UMBRELLA", "PLATFORM",
                "SEAMLESS", "RESPONSE ADAPTIVE", "ENRICHMENT DESIGN",
                "MASTER PROTOCOL",
            ],
            TrialWorkflowType.BIOMARKER_STRATEGY: [
                "BIOMARKER", "COMPANION DIAGNOSTIC", "CDX", "PREDICTIVE",
                "PROGNOSTIC", "CTDNA", "LIQUID BIOPSY", "GENOMIC",
                "MOLECULAR", "PRECISION", "STRATIFICATION",
            ],
            TrialWorkflowType.RWE_ANALYSIS: [
                "REAL WORLD", "RWE", "RWD", "REGISTRY", "CLAIMS DATA",
                "EHR", "ELECTRONIC HEALTH RECORD", "PRAGMATIC",
                "OBSERVATIONAL", "POST-MARKETING",
            ],
            TrialWorkflowType.RECRUITMENT_OPTIMIZATION: [
                "RECRUITMENT", "ENROLLMENT", "RETENTION", "DIVERSITY",
                "MINORITY", "UNDERREPRESENTED", "ACCRUAL",
                "RECRUITMENT RATE", "SCREEN-TO-ENROLL",
            ],
        }

        for wf, keywords in keyword_workflow_map.items():
            for kw in keywords:
                if kw in text_upper:
                    workflow_scores[wf] = workflow_scores.get(wf, 0) + 1.0

        if not workflow_scores:
            return TrialWorkflowType.GENERAL

        sorted_workflows = sorted(
            workflow_scores.items(), key=lambda x: x[1], reverse=True,
        )

        return sorted_workflows[0][0]

    # ── Entity Detection ────────────────────────────────────────────

    def _detect_entities(self, question: str) -> Dict[str, List[str]]:
        """Detect clinical trial entities in the question text.

        Scans for conditions, drugs, and biomarkers using the knowledge
        dictionaries. Performs case-insensitive matching against canonical
        names and aliases.

        Args:
            question: The user's natural language question.

        Returns:
            Dict with keys 'conditions', 'drugs', 'biomarkers' mapping
            to lists of detected entity names.
        """
        import re

        entities: Dict[str, List[str]] = {
            "conditions": [],
            "drugs": [],
            "biomarkers": [],
        }

        text_lower = question.lower()

        # Detect conditions
        for condition, info in TRIAL_CONDITIONS.items():
            if condition in text_lower:
                if condition not in entities["conditions"]:
                    entities["conditions"].append(condition)
                continue
            aliases = info.get("aliases", [])
            for alias in aliases:
                if len(alias) <= 3:
                    pattern = r'\b' + re.escape(alias) + r'\b'
                    if re.search(pattern, text_lower):
                        if condition not in entities["conditions"]:
                            entities["conditions"].append(condition)
                        break
                else:
                    if alias.lower() in text_lower:
                        if condition not in entities["conditions"]:
                            entities["conditions"].append(condition)
                        break

        # Detect drugs
        for drug, info in TRIAL_DRUGS.items():
            if drug.lower() in text_lower:
                if drug not in entities["drugs"]:
                    entities["drugs"].append(drug)
                continue
            aliases = info.get("aliases", [])
            for alias in aliases:
                if alias.lower() in text_lower:
                    if drug not in entities["drugs"]:
                        entities["drugs"].append(drug)
                    break

        # Detect biomarkers
        for biomarker, info in TRIAL_BIOMARKERS.items():
            if biomarker.lower() in text_lower:
                if biomarker not in entities["biomarkers"]:
                    entities["biomarkers"].append(biomarker)
                continue
            full_name = info.get("full_name", "")
            if full_name and full_name.lower() in text_lower:
                if biomarker not in entities["biomarkers"]:
                    entities["biomarkers"].append(biomarker)

        return entities

    # ── Search Strategy ─────────────────────────────────────────────

    def _build_search_strategy(
        self,
        entities: Dict[str, List[str]],
        workflow: TrialWorkflowType,
    ) -> str:
        """Build a descriptive search strategy based on entities and workflow.

        Args:
            entities: Detected entities dict from _detect_entities.
            workflow: Determined workflow type.

        Returns:
            Strategy description string for logging/debugging.
        """
        parts = [f"Workflow: {workflow.value}"]

        if entities.get("conditions"):
            parts.append(f"Conditions: {', '.join(entities['conditions'])}")
        if entities.get("drugs"):
            parts.append(f"Drugs: {', '.join(entities['drugs'])}")
        if entities.get("biomarkers"):
            parts.append(f"Biomarkers: {', '.join(entities['biomarkers'])}")

        # Determine collection priorities
        boosts = WORKFLOW_COLLECTION_BOOST.get(workflow, {})
        top_collections = sorted(
            boosts.items(), key=lambda x: x[1], reverse=True,
        )[:5]
        if top_collections:
            parts.append(
                "Priority collections: "
                + ", ".join(f"{c}({w:.1f}x)" for c, w in top_collections)
            )

        return " | ".join(parts)

    def _choose_strategy(
        self,
        text: str,
        conditions: List[str],
        drugs: List[str],
    ) -> str:
        """Choose search strategy: broad, targeted, comparative, or regulatory.

        Args:
            text: Original query text.
            conditions: Detected conditions.
            drugs: Detected drugs.

        Returns:
            Strategy name string.
        """
        text_upper = text.upper()

        # Comparative queries
        if ("COMPARE" in text_upper or " VS " in text_upper
                or "VERSUS" in text_upper or "DIFFERENCE BETWEEN" in text_upper
                or "COMPETING" in text_upper or "HEAD-TO-HEAD" in text_upper):
            return "comparative"

        # Regulatory-focused queries
        regulatory_keywords = [
            "FDA", "EMA", "IND", "NDA", "BLA", "ICH", "REGULATORY",
            "GUIDANCE", "ACCELERATED APPROVAL", "BREAKTHROUGH",
        ]
        if any(kw in text_upper for kw in regulatory_keywords):
            return "regulatory"

        # Targeted: specific condition + drug or single focused entity
        if (len(conditions) == 1 and len(drugs) <= 1) or (
            len(conditions) <= 1 and len(drugs) == 1
        ):
            if conditions or drugs:
                return "targeted"

        return "broad"

    # ── Sub-Question Generation ────────────────────────────────────

    def _generate_sub_questions(self, plan: SearchPlan) -> List[str]:
        """Generate sub-questions for comprehensive retrieval.

        Decomposes the main question into focused sub-queries based on
        the detected entities and workflow type. Enables multi-hop
        retrieval across different aspects of the clinical trial question.

        Args:
            plan: SearchPlan with detected entities and workflows.

        Returns:
            List of sub-question strings (typically 2-4 questions).
        """
        sub_questions: List[str] = []

        condition_label = plan.conditions[0] if plan.conditions else "the indication"
        drug_label = plan.drugs[0] if plan.drugs else "the investigational agent"
        biomarker_label = plan.biomarkers[0] if plan.biomarkers else "the biomarker"

        primary_wf = plan.relevant_workflows[0] if plan.relevant_workflows else TrialWorkflowType.GENERAL

        # ── Pattern 1: Protocol Design ────────────────────────────
        if primary_wf == TrialWorkflowType.PROTOCOL_DESIGN:
            sub_questions = [
                f"What are the key design considerations for clinical trials in {condition_label}?",
                f"What endpoints have been used in successful {condition_label} trials?",
                f"What regulatory precedents exist for {condition_label} trial approval?",
                f"What sample size and statistical approaches are typical for {condition_label}?",
            ]

        # ── Pattern 2: Eligibility Analysis ──────────────────────
        elif primary_wf == TrialWorkflowType.ELIGIBILITY_ANALYSIS:
            sub_questions = [
                f"What are standard inclusion criteria for {condition_label} trials?",
                f"What common exclusion criteria cause screen failures in {condition_label}?",
                f"How do eligibility criteria impact enrollment diversity in {condition_label}?",
                f"What biomarker requirements apply to {condition_label} trials?",
            ]

        # ── Pattern 3: Regulatory Strategy ───────────────────────
        elif primary_wf == TrialWorkflowType.REGULATORY_STRATEGY:
            sub_questions = [
                f"What FDA regulatory pathways apply to {drug_label} for {condition_label}?",
                f"What ICH guidelines are most relevant to {condition_label} trials?",
                f"What precedent decisions exist for {condition_label} approvals?",
                f"What are the key regulatory risks for {drug_label}?",
            ]

        # ── Pattern 4: Competitive Intelligence ──────────────────
        elif primary_wf == TrialWorkflowType.COMPETITIVE_INTELLIGENCE:
            sub_questions = [
                f"What competing trials are active in {condition_label}?",
                f"How does {drug_label} compare to other agents in {condition_label}?",
                f"What are the key differentiators for {condition_label} development programs?",
                f"What are the projected timelines for competing {condition_label} programs?",
            ]

        # ── Pattern 5: Safety Monitoring ─────────────────────────
        elif primary_wf == TrialWorkflowType.SAFETY_MONITORING:
            sub_questions = [
                f"What are the known safety signals for {drug_label}?",
                f"What DSMB charter provisions apply to {condition_label} trials?",
                f"What safety monitoring schedule is recommended for {drug_label}?",
                f"What stopping rules apply to {condition_label} studies?",
            ]

        # ── Pattern 6: Endpoint Strategy ─────────────────────────
        elif primary_wf == TrialWorkflowType.ENDPOINT_STRATEGY:
            sub_questions = [
                f"What primary endpoints have FDA accepted for {condition_label}?",
                f"What surrogate endpoints are validated for {condition_label}?",
                f"What patient-reported outcomes are used in {condition_label} trials?",
                f"What composite endpoints have been successful in {condition_label}?",
            ]

        # ── Pattern 7: Biomarker Strategy ────────────────────────
        elif primary_wf == TrialWorkflowType.BIOMARKER_STRATEGY:
            sub_questions = [
                f"What is the validation status of {biomarker_label} in {condition_label}?",
                f"What companion diagnostic requirements exist for {biomarker_label}?",
                f"How should {biomarker_label} be incorporated into trial enrichment?",
                f"What regulatory guidance applies to {biomarker_label} as a trial endpoint?",
            ]

        # ── Pattern 8: Adaptive Design ───────────────────────────
        elif primary_wf == TrialWorkflowType.ADAPTIVE_DESIGN:
            sub_questions = [
                f"What adaptive design options are suitable for {condition_label}?",
                f"What interim analysis timing is optimal for {condition_label} trials?",
                f"What regulatory guidance applies to adaptive designs in {condition_label}?",
                f"What operational challenges exist for adaptive {condition_label} trials?",
            ]

        # ── Pattern 9: Comparative ───────────────────────────────
        elif plan.search_strategy == "comparative":
            entities = plan.conditions + plan.drugs
            if len(entities) >= 2:
                sub_questions = [
                    f"What is the evidence for {entities[0]} in clinical trials?",
                    f"What is the evidence for {entities[1]} in clinical trials?",
                    f"What head-to-head trials compare {entities[0]} and {entities[1]}?",
                ]
            else:
                sub_questions = [
                    f"What are the treatment options being studied for {condition_label}?",
                    f"What comparative evidence exists for {condition_label} therapies?",
                ]

        # ── Default ──────────────────────────────────────────────
        else:
            sub_questions = [
                f"What are the current clinical trials for {condition_label}?",
                f"What is the regulatory landscape for {condition_label}?",
                f"What evidence supports clinical development in {condition_label}?",
            ]

        return sub_questions
