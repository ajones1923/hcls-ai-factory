"""Rare Disease Diagnostic Agent -- autonomous reasoning for rare disease diagnostics.

Implements the plan -> search -> evaluate -> synthesize -> report pattern from the
VAST AI OS AgentEngine model. The agent can:

1. Parse complex phenotype descriptions coded in HPO terms
2. Plan a search strategy across 14 domain-specific collections
3. Execute multi-collection retrieval via the RareDiseaseRAGEngine
4. Evaluate evidence quality using ACMG criteria and diagnostic certainty
5. Synthesize differential diagnoses with gene-disease-phenotype linking
6. Generate structured diagnostic reports with OMIM/Orphanet references

Mapping to VAST AI OS:
  - AgentEngine entry point: RareDiseaseAgent.run()
  - Plan -> search_plan()
  - Execute -> rag_engine.query()
  - Reflect -> evaluate_evidence()
  - Report -> generate_report()

Author: Adam Jones
Date: March 2026
"""

import re
from dataclasses import dataclass, field
from datetime import datetime, timezone
from enum import Enum
from typing import Dict, List, Optional

from src.models import DiagnosticWorkflowType


# =====================================================================
# ENUMS
# =====================================================================


class EvidenceLevel(str, Enum):
    """Clinical evidence hierarchy for rare disease diagnostics."""
    LEVEL_1 = "1"         # Systematic review / meta-analysis
    LEVEL_2 = "2"         # Well-designed controlled study
    LEVEL_3 = "3"         # Cohort / case-control study
    LEVEL_4 = "4"         # Case series / case report
    LEVEL_5 = "5"         # Expert opinion / consensus
    GUIDELINE = "guideline"   # Practice guideline (ACMG, ACOG, NORD)
    REGISTRY = "registry"     # Patient registry / natural history study


class ACMGClassification(str, Enum):
    """ACMG/AMP variant classification."""
    PATHOGENIC = "pathogenic"
    LIKELY_PATHOGENIC = "likely_pathogenic"
    VUS = "uncertain_significance"
    LIKELY_BENIGN = "likely_benign"
    BENIGN = "benign"


class InheritancePattern(str, Enum):
    """Mendelian inheritance patterns."""
    AUTOSOMAL_DOMINANT = "autosomal_dominant"
    AUTOSOMAL_RECESSIVE = "autosomal_recessive"
    X_LINKED_DOMINANT = "x_linked_dominant"
    X_LINKED_RECESSIVE = "x_linked_recessive"
    MITOCHONDRIAL = "mitochondrial"
    DIGENIC = "digenic"
    MULTIFACTORIAL = "multifactorial"
    DE_NOVO = "de_novo"


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
class DiagnosticResult:
    """Complete response from the rare disease diagnostic agent.

    Attributes:
        question: Original query text.
        answer: LLM-synthesised diagnostic assessment.
        results: Ranked search results used for synthesis.
        workflow: Diagnostic workflow that was applied.
        confidence: Overall confidence score (0.0 - 1.0).
        citations: Formatted citation list.
        search_time_ms: Total search time in milliseconds.
        collections_searched: Number of collections queried.
        patient_context_used: Whether patient context was injected.
        differential_diagnoses: Ordered list of candidate diagnoses.
        timestamp: ISO 8601 timestamp of response generation.
    """
    question: str = ""
    answer: str = ""
    results: list = field(default_factory=list)
    workflow: Optional[DiagnosticWorkflowType] = None
    confidence: float = 0.0
    citations: List[Dict[str, str]] = field(default_factory=list)
    search_time_ms: float = 0.0
    collections_searched: int = 0
    patient_context_used: bool = False
    differential_diagnoses: List[str] = field(default_factory=list)
    timestamp: str = field(
        default_factory=lambda: datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    )


# =====================================================================
# RARE DISEASE SYSTEM PROMPT
# =====================================================================

RARE_DISEASE_SYSTEM_PROMPT = """\
You are a rare disease diagnostic intelligence system within the HCLS AI Factory. \
You have deep expertise in clinical genetics, dysmorphology, metabolic medicine, \
and the diagnostic odyssey faced by rare disease patients. You analyze patient \
phenotypes coded in HPO terms, interpret genomic variants using ACMG criteria, \
identify therapeutic options including gene therapies and orphan drugs, and \
connect patients with clinical trials. You reference OMIM, Orphanet, \
GeneReviews, and ClinVar. You are aware that 7,000+ rare diseases exist, \
affecting 300M people globally, with an average diagnostic odyssey of 5-7 years. \
Your goal is to shorten that odyssey to hours. Always cite specific OMIM \
numbers, gene symbols, HPO terms, and evidence levels. Never fabricate \
genetic data.

Your responses must adhere to the following standards:

1. **HPO Term Usage** -- Always use Human Phenotype Ontology terms with their \
   identifiers when describing phenotypic features. Format as: \
   [HP:0001250](https://hpo.jax.org/browse/term/HP:0001250) Seizures. \
   Group related phenotypes by organ system (neurological, musculoskeletal, \
   craniofacial, metabolic, ophthalmologic, cardiac, dermatologic, etc.).

2. **OMIM / Orphanet References** -- Cite diseases using OMIM MIM numbers \
   and Orphanet identifiers with clickable links: \
   [OMIM:219700](https://omim.org/entry/219700) Cystic Fibrosis; \
   [ORPHA:586](https://www.orpha.net/consor/cgi-bin/OC_Exp.php?lng=en&Expert=586) \
   Cystic Fibrosis. Include gene symbol, inheritance pattern, and \
   prevalence when available.

3. **ACMG Variant Classification** -- For variant interpretation queries, \
   apply the ACMG/AMP 2015 guidelines (Richards et al., Genet Med 17:405-424). \
   Classify variants as Pathogenic, Likely Pathogenic, VUS, Likely Benign, \
   or Benign. List supporting evidence criteria (PVS1, PS1-PS4, PM1-PM6, \
   PP1-PP5, BS1-BS4, BP1-BP7, BA1). Reference ClinVar assertion counts \
   and conflicting interpretations.

4. **CRITICAL Findings** -- Flag the following as CRITICAL with prominent \
   visual markers and immediate action recommendations:
   - Metabolic emergencies (hyperammonemia, hypoglycemia, acidosis)
   - Variants in actionable genes with available therapies
   - Newborn screening follow-up requiring urgent confirmatory testing
   - Pathogenic variants in genes with gene therapy trials enrolling
   - Carrier status for X-linked lethal conditions in pregnant patients
   - Progressive conditions where early intervention alters prognosis
   - Pharmacogenomic interactions with orphan drugs
   - Consanguinity with autosomal recessive findings

5. **Severity Badges** -- Classify all findings using standardised severity \
   levels: [CRITICAL], [HIGH], [MODERATE], [LOW], [INFORMATIONAL]. Place \
   the badge at the start of each finding or recommendation line.

6. **Differential Diagnosis Format** -- Organise differential diagnoses as \
   a ranked list with: disease name, OMIM number, causative gene(s), \
   inheritance pattern, overlapping phenotypes (HPO terms), distinguishing \
   features, and diagnostic confirmation method. Include both common and \
   rare differentials. Assign a phenotype match score (0-100%) based on \
   HPO term overlap.

7. **Structured Formatting** -- Organise responses with clear sections: \
   Clinical Summary, Phenotype Analysis, Differential Diagnosis, Variant \
   Interpretation, Recommended Testing, Therapeutic Options, Clinical \
   Trials, Natural History, Genetic Counselling Considerations, and \
   References. Use bullet points and numbered lists for actionable items.

8. **Gene-Disease Validity** -- Reference ClinGen gene-disease validity \
   classifications (Definitive, Strong, Moderate, Limited, Disputed, \
   Refuted) when discussing gene-phenotype associations. Include \
   ClinGen actionability scores where relevant.

9. **Therapeutic Landscape** -- For therapy queries, cover: FDA/EMA-approved \
   orphan drugs, gene therapies (AAV, lentiviral, mRNA, base editing, \
   prime editing), enzyme replacement therapies (ERT), substrate reduction \
   therapies (SRT), antisense oligonucleotides (ASOs), read-through agents, \
   chaperone therapies, and symptomatic management. Include regulatory \
   designations (Orphan Drug, Breakthrough Therapy, PRIME, RMAT).

10. **Limitations** -- You are a rare disease diagnostic support tool. You \
    do NOT replace clinical geneticists, genetic counselors, or metabolic \
    specialists. All diagnostic suggestions require validation through \
    appropriate clinical and laboratory testing. Variant interpretations \
    should be confirmed by CLIA/CAP-certified laboratories. Explicitly \
    state when evidence is limited or when specialist consultation is \
    recommended. Never provide a definitive diagnosis -- always present \
    findings as differential diagnoses requiring clinical correlation."""


# =====================================================================
# WORKFLOW-SPECIFIC COLLECTION BOOST WEIGHTS
# =====================================================================
# Maps each DiagnosticWorkflowType to collection weight overrides (multipliers).
# Collections not listed retain their base weight (1.0x). Values > 1.0
# boost the collection; values < 1.0 would suppress it.

WORKFLOW_COLLECTION_BOOST: Dict[DiagnosticWorkflowType, Dict[str, float]] = {

    # -- Phenotype-Driven Diagnosis ------------------------------------
    DiagnosticWorkflowType.PHENOTYPE_DRIVEN: {
        "rd_phenotypes": 2.5,
        "rd_diseases": 2.0,
        "rd_genes": 1.5,
        "rd_case_reports": 1.5,
        "rd_literature": 1.3,
        "rd_natural_history": 1.2,
        "rd_guidelines": 1.1,
        "rd_variants": 1.0,
    },

    # -- Variant Interpretation ----------------------------------------
    DiagnosticWorkflowType.VARIANT_INTERPRETATION: {
        "rd_variants": 2.5,
        "rd_genes": 2.0,
        "genomic_evidence": 2.0,
        "rd_diseases": 1.5,
        "rd_literature": 1.3,
        "rd_guidelines": 1.3,
        "rd_phenotypes": 1.1,
        "rd_case_reports": 1.0,
    },

    # -- Differential Diagnosis ----------------------------------------
    DiagnosticWorkflowType.DIFFERENTIAL_DIAGNOSIS: {
        "rd_diseases": 2.5,
        "rd_phenotypes": 2.0,
        "rd_genes": 1.8,
        "rd_case_reports": 1.5,
        "rd_literature": 1.3,
        "rd_natural_history": 1.2,
        "rd_variants": 1.1,
        "rd_guidelines": 1.0,
    },

    # -- Gene Therapy Eligibility --------------------------------------
    DiagnosticWorkflowType.GENE_THERAPY_ELIGIBILITY: {
        "rd_therapies": 2.5,
        "rd_trials": 2.0,
        "rd_genes": 1.8,
        "rd_diseases": 1.5,
        "rd_guidelines": 1.3,
        "rd_variants": 1.2,
        "rd_natural_history": 1.1,
        "rd_literature": 1.0,
    },

    # -- Newborn Screening Triage --------------------------------------
    DiagnosticWorkflowType.NEWBORN_SCREENING: {
        "rd_newborn_screening": 2.5,
        "rd_diseases": 2.0,
        "rd_therapies": 1.8,
        "rd_genes": 1.5,
        "rd_natural_history": 1.5,
        "rd_guidelines": 1.3,
        "rd_pathways": 1.2,
        "rd_variants": 1.0,
    },

    # -- Metabolic Crisis Workup ---------------------------------------
    DiagnosticWorkflowType.METABOLIC_WORKUP: {
        "rd_pathways": 2.5,
        "rd_diseases": 2.0,
        "rd_phenotypes": 1.8,
        "rd_therapies": 1.5,
        "rd_guidelines": 1.5,
        "rd_genes": 1.3,
        "rd_newborn_screening": 1.2,
        "rd_natural_history": 1.1,
    },

    # -- Carrier Screening ---------------------------------------------
    DiagnosticWorkflowType.CARRIER_SCREENING: {
        "rd_genes": 2.5,
        "rd_variants": 2.0,
        "rd_diseases": 1.8,
        "genomic_evidence": 1.5,
        "rd_guidelines": 1.5,
        "rd_literature": 1.2,
        "rd_registries": 1.1,
        "rd_phenotypes": 1.0,
    },

    # -- Prenatal Diagnosis --------------------------------------------
    DiagnosticWorkflowType.PRENATAL_DIAGNOSIS: {
        "rd_genes": 2.5,
        "rd_variants": 2.0,
        "rd_diseases": 1.8,
        "rd_guidelines": 1.5,
        "rd_natural_history": 1.5,
        "rd_therapies": 1.3,
        "rd_literature": 1.2,
        "genomic_evidence": 1.0,
    },

    # -- Natural History -----------------------------------------------
    DiagnosticWorkflowType.NATURAL_HISTORY: {
        "rd_natural_history": 2.5,
        "rd_diseases": 2.0,
        "rd_registries": 1.8,
        "rd_literature": 1.5,
        "rd_phenotypes": 1.3,
        "rd_case_reports": 1.2,
        "rd_guidelines": 1.1,
        "rd_therapies": 1.0,
    },

    # -- Therapy Selection ---------------------------------------------
    DiagnosticWorkflowType.THERAPY_SELECTION: {
        "rd_therapies": 2.5,
        "rd_guidelines": 2.0,
        "rd_trials": 1.8,
        "rd_diseases": 1.5,
        "rd_genes": 1.3,
        "rd_literature": 1.2,
        "rd_natural_history": 1.1,
        "rd_variants": 1.0,
    },

    # -- Clinical Trial Matching ---------------------------------------
    DiagnosticWorkflowType.CLINICAL_TRIAL_MATCHING: {
        "rd_trials": 2.5,
        "rd_therapies": 2.0,
        "rd_diseases": 1.5,
        "rd_genes": 1.5,
        "rd_variants": 1.3,
        "rd_registries": 1.2,
        "rd_guidelines": 1.0,
        "rd_literature": 1.0,
    },

    # -- Genetic Counseling --------------------------------------------
    DiagnosticWorkflowType.GENETIC_COUNSELING: {
        "rd_genes": 2.5,
        "rd_diseases": 2.0,
        "rd_variants": 1.8,
        "rd_natural_history": 1.5,
        "rd_guidelines": 1.5,
        "rd_registries": 1.3,
        "rd_literature": 1.2,
        "rd_phenotypes": 1.0,
    },

    # -- General (balanced across all collections) ---------------------
    DiagnosticWorkflowType.GENERAL: {
        "rd_phenotypes": 1.2,
        "rd_diseases": 1.2,
        "rd_genes": 1.1,
        "rd_variants": 1.1,
        "rd_literature": 1.2,
        "rd_case_reports": 1.0,
        "rd_therapies": 1.0,
        "rd_trials": 0.9,
        "rd_guidelines": 1.1,
        "rd_pathways": 0.9,
        "rd_registries": 0.8,
        "rd_natural_history": 1.0,
        "rd_newborn_screening": 0.8,
        "genomic_evidence": 0.8,
    },
}


# =====================================================================
# KNOWLEDGE DOMAIN DICTIONARIES
# =====================================================================
# Comprehensive rare disease knowledge for entity detection and context
# enrichment. Used by the agent's search_plan() to identify entities
# in user queries and map them to workflows.

RARE_DISEASE_CONDITIONS: Dict[str, Dict[str, object]] = {

    # -- Lysosomal Storage Disorders -----------------------------------
    "gaucher disease": {
        "omim": "OMIM:230800",
        "gene": "GBA1",
        "inheritance": "autosomal_recessive",
        "aliases": ["gaucher", "glucocerebrosidase deficiency", "gd type 1",
                    "gd type 2", "gd type 3"],
        "workflows": [DiagnosticWorkflowType.THERAPY_SELECTION,
                       DiagnosticWorkflowType.METABOLIC_WORKUP],
        "search_terms": ["GBA1", "glucocerebrosidase", "ERT", "imiglucerase",
                         "eliglustat", "splenomegaly", "Erlenmeyer flask deformity"],
    },
    "fabry disease": {
        "omim": "OMIM:301500",
        "gene": "GLA",
        "inheritance": "x_linked",
        "aliases": ["fabry", "anderson-fabry disease",
                    "alpha-galactosidase a deficiency"],
        "workflows": [DiagnosticWorkflowType.THERAPY_SELECTION,
                       DiagnosticWorkflowType.CARRIER_SCREENING],
        "search_terms": ["GLA", "alpha-galactosidase", "agalsidase",
                         "migalastat", "acroparesthesia", "angiokeratoma"],
    },
    "pompe disease": {
        "omim": "OMIM:232300",
        "gene": "GAA",
        "inheritance": "autosomal_recessive",
        "aliases": ["pompe", "glycogen storage disease type ii",
                    "gsd ii", "acid maltase deficiency"],
        "workflows": [DiagnosticWorkflowType.NEWBORN_SCREENING,
                       DiagnosticWorkflowType.THERAPY_SELECTION],
        "search_terms": ["GAA", "acid alpha-glucosidase", "alglucosidase alfa",
                         "avalglucosidase alfa", "cardiomyopathy", "hypotonia"],
    },
    "mucopolysaccharidosis type i": {
        "omim": "OMIM:607014",
        "gene": "IDUA",
        "inheritance": "autosomal_recessive",
        "aliases": ["mps i", "hurler syndrome", "hurler-scheie syndrome",
                    "scheie syndrome", "mps1"],
        "workflows": [DiagnosticWorkflowType.GENE_THERAPY_ELIGIBILITY,
                       DiagnosticWorkflowType.NEWBORN_SCREENING],
        "search_terms": ["IDUA", "alpha-L-iduronidase", "laronidase",
                         "HSCT", "dysostosis multiplex", "corneal clouding"],
    },
    "niemann-pick disease type c": {
        "omim": "OMIM:257220",
        "gene": "NPC1",
        "inheritance": "autosomal_recessive",
        "aliases": ["npc", "niemann-pick c", "niemann pick type c"],
        "workflows": [DiagnosticWorkflowType.DIFFERENTIAL_DIAGNOSIS,
                       DiagnosticWorkflowType.METABOLIC_WORKUP],
        "search_terms": ["NPC1", "NPC2", "miglustat", "filipin staining",
                         "vertical supranuclear gaze palsy", "hepatosplenomegaly"],
    },

    # -- Neuromuscular Disorders ---------------------------------------
    "spinal muscular atrophy": {
        "omim": "OMIM:253300",
        "gene": "SMN1",
        "inheritance": "autosomal_recessive",
        "aliases": ["sma", "sma type 1", "sma type 2", "sma type 3",
                    "werdnig-hoffmann disease"],
        "workflows": [DiagnosticWorkflowType.GENE_THERAPY_ELIGIBILITY,
                       DiagnosticWorkflowType.NEWBORN_SCREENING],
        "search_terms": ["SMN1", "SMN2", "nusinersen", "onasemnogene",
                         "risdiplam", "motor neuron", "hypotonia"],
    },
    "duchenne muscular dystrophy": {
        "omim": "OMIM:310200",
        "gene": "DMD",
        "inheritance": "x_linked_recessive",
        "aliases": ["dmd", "duchenne", "becker muscular dystrophy",
                    "dystrophinopathy"],
        "workflows": [DiagnosticWorkflowType.GENE_THERAPY_ELIGIBILITY,
                       DiagnosticWorkflowType.VARIANT_INTERPRETATION],
        "search_terms": ["DMD", "dystrophin", "exon skipping", "eteplirsen",
                         "delandistrogene", "creatine kinase", "Gowers sign"],
    },
    "myotonic dystrophy type 1": {
        "omim": "OMIM:160900",
        "gene": "DMPK",
        "inheritance": "autosomal_dominant",
        "aliases": ["dm1", "steinert disease", "myotonic dystrophy",
                    "dystrophia myotonica"],
        "workflows": [DiagnosticWorkflowType.GENETIC_COUNSELING,
                       DiagnosticWorkflowType.NATURAL_HISTORY],
        "search_terms": ["DMPK", "CTG repeat", "myotonia", "anticipation",
                         "cataracts", "cardiac conduction"],
    },
    "charcot-marie-tooth disease": {
        "omim": "OMIM:118220",
        "gene": "PMP22",
        "inheritance": "autosomal_dominant",
        "aliases": ["cmt", "cmt1a", "hereditary motor sensory neuropathy",
                    "hmsn", "peroneal muscular atrophy"],
        "workflows": [DiagnosticWorkflowType.DIFFERENTIAL_DIAGNOSIS,
                       DiagnosticWorkflowType.VARIANT_INTERPRETATION],
        "search_terms": ["PMP22", "MPZ", "GJB1", "MFN2", "nerve conduction",
                         "pes cavus", "distal weakness"],
    },

    # -- Connective Tissue Disorders -----------------------------------
    "marfan syndrome": {
        "omim": "OMIM:154700",
        "gene": "FBN1",
        "inheritance": "autosomal_dominant",
        "aliases": ["marfan", "mfs"],
        "workflows": [DiagnosticWorkflowType.PHENOTYPE_DRIVEN,
                       DiagnosticWorkflowType.GENETIC_COUNSELING],
        "search_terms": ["FBN1", "fibrillin-1", "aortic root dilation",
                         "lens subluxation", "Ghent criteria", "losartan"],
    },
    "ehlers-danlos syndrome": {
        "omim": "OMIM:130000",
        "gene": "COL5A1",
        "inheritance": "autosomal_dominant",
        "aliases": ["eds", "eds hypermobile", "eds vascular",
                    "heds", "veds", "ehlers danlos"],
        "workflows": [DiagnosticWorkflowType.PHENOTYPE_DRIVEN,
                       DiagnosticWorkflowType.DIFFERENTIAL_DIAGNOSIS],
        "search_terms": ["COL5A1", "COL3A1", "TNXB", "joint hypermobility",
                         "skin hyperextensibility", "Beighton score"],
    },
    "osteogenesis imperfecta": {
        "omim": "OMIM:166200",
        "gene": "COL1A1",
        "inheritance": "autosomal_dominant",
        "aliases": ["oi", "brittle bone disease",
                    "osteogenesis imperfecta type i"],
        "workflows": [DiagnosticWorkflowType.PHENOTYPE_DRIVEN,
                       DiagnosticWorkflowType.THERAPY_SELECTION],
        "search_terms": ["COL1A1", "COL1A2", "collagen type I", "blue sclerae",
                         "fractures", "bisphosphonates", "dentinogenesis imperfecta"],
    },

    # -- Metabolic Disorders -------------------------------------------
    "phenylketonuria": {
        "omim": "OMIM:261600",
        "gene": "PAH",
        "inheritance": "autosomal_recessive",
        "aliases": ["pku", "hyperphenylalaninemia", "pah deficiency"],
        "workflows": [DiagnosticWorkflowType.NEWBORN_SCREENING,
                       DiagnosticWorkflowType.METABOLIC_WORKUP],
        "search_terms": ["PAH", "phenylalanine hydroxylase", "phenylalanine",
                         "sapropterin", "pegvaliase", "BH4 responsive"],
    },
    "maple syrup urine disease": {
        "omim": "OMIM:248600",
        "gene": "BCKDHA",
        "inheritance": "autosomal_recessive",
        "aliases": ["msud", "branched-chain ketoaciduria"],
        "workflows": [DiagnosticWorkflowType.NEWBORN_SCREENING,
                       DiagnosticWorkflowType.METABOLIC_WORKUP],
        "search_terms": ["BCKDHA", "BCKDHB", "DBT", "branched-chain amino acids",
                         "leucine", "encephalopathy", "metabolic crisis"],
    },
    "ornithine transcarbamylase deficiency": {
        "omim": "OMIM:311250",
        "gene": "OTC",
        "inheritance": "x_linked",
        "aliases": ["otc deficiency", "otcd",
                    "ornithine carbamoyltransferase deficiency"],
        "workflows": [DiagnosticWorkflowType.METABOLIC_WORKUP,
                       DiagnosticWorkflowType.GENE_THERAPY_ELIGIBILITY],
        "search_terms": ["OTC", "hyperammonemia", "urea cycle", "citrulline",
                         "orotic acid", "nitrogen scavengers", "liver transplant"],
    },
    "medium-chain acyl-coa dehydrogenase deficiency": {
        "omim": "OMIM:201450",
        "gene": "ACADM",
        "inheritance": "autosomal_recessive",
        "aliases": ["mcadd", "mcad deficiency"],
        "workflows": [DiagnosticWorkflowType.NEWBORN_SCREENING,
                       DiagnosticWorkflowType.METABOLIC_WORKUP],
        "search_terms": ["ACADM", "acylcarnitine", "C8", "hypoketotic hypoglycemia",
                         "fasting intolerance", "newborn screening"],
    },
    "galactosemia": {
        "omim": "OMIM:230400",
        "gene": "GALT",
        "inheritance": "autosomal_recessive",
        "aliases": ["classic galactosemia", "galactose-1-phosphate uridylyltransferase deficiency"],
        "workflows": [DiagnosticWorkflowType.NEWBORN_SCREENING,
                       DiagnosticWorkflowType.METABOLIC_WORKUP],
        "search_terms": ["GALT", "galactose", "Beutler test", "E. coli sepsis",
                         "cataracts", "galactose-restricted diet"],
    },

    # -- Neurogenetic Disorders ----------------------------------------
    "rett syndrome": {
        "omim": "OMIM:312750",
        "gene": "MECP2",
        "inheritance": "x_linked_dominant",
        "aliases": ["rett", "mecp2 disorder"],
        "workflows": [DiagnosticWorkflowType.PHENOTYPE_DRIVEN,
                       DiagnosticWorkflowType.GENE_THERAPY_ELIGIBILITY],
        "search_terms": ["MECP2", "hand stereotypies", "regression",
                         "trofinetide", "breathing irregularities", "microcephaly"],
    },
    "fragile x syndrome": {
        "omim": "OMIM:300624",
        "gene": "FMR1",
        "inheritance": "x_linked_dominant",
        "aliases": ["fragile x", "fra(x)", "fxs",
                    "martin-bell syndrome"],
        "workflows": [DiagnosticWorkflowType.GENETIC_COUNSELING,
                       DiagnosticWorkflowType.CARRIER_SCREENING],
        "search_terms": ["FMR1", "CGG repeat", "FMRP", "intellectual disability",
                         "macroorchidism", "premutation", "FXTAS"],
    },
    "huntington disease": {
        "omim": "OMIM:143100",
        "gene": "HTT",
        "inheritance": "autosomal_dominant",
        "aliases": ["huntington", "huntington's disease", "hd",
                    "huntington chorea"],
        "workflows": [DiagnosticWorkflowType.GENETIC_COUNSELING,
                       DiagnosticWorkflowType.NATURAL_HISTORY],
        "search_terms": ["HTT", "huntingtin", "CAG repeat", "chorea",
                         "anticipation", "presymptomatic testing", "tominersen"],
    },
    "tuberous sclerosis complex": {
        "omim": "OMIM:191100",
        "gene": "TSC1",
        "inheritance": "autosomal_dominant",
        "aliases": ["tsc", "tuberous sclerosis", "bourneville disease"],
        "workflows": [DiagnosticWorkflowType.PHENOTYPE_DRIVEN,
                       DiagnosticWorkflowType.THERAPY_SELECTION],
        "search_terms": ["TSC1", "TSC2", "mTOR", "everolimus",
                         "angiomyolipoma", "cortical tubers", "subependymal nodules"],
    },

    # -- Hematologic Disorders -----------------------------------------
    "sickle cell disease": {
        "omim": "OMIM:603903",
        "gene": "HBB",
        "inheritance": "autosomal_recessive",
        "aliases": ["sickle cell", "scd", "sickle cell anemia",
                    "hb ss disease"],
        "workflows": [DiagnosticWorkflowType.GENE_THERAPY_ELIGIBILITY,
                       DiagnosticWorkflowType.NEWBORN_SCREENING],
        "search_terms": ["HBB", "hemoglobin S", "vaso-occlusive", "exagamglogene",
                         "lovotibeglogene", "hydroxyurea", "crizanlizumab"],
    },
    "beta-thalassemia": {
        "omim": "OMIM:613985",
        "gene": "HBB",
        "inheritance": "autosomal_recessive",
        "aliases": ["thalassemia major", "cooley anemia",
                    "beta-thal", "transfusion-dependent thalassemia"],
        "workflows": [DiagnosticWorkflowType.GENE_THERAPY_ELIGIBILITY,
                       DiagnosticWorkflowType.CARRIER_SCREENING],
        "search_terms": ["HBB", "hemoglobin electrophoresis", "betibeglogene",
                         "exagamglogene", "luspatercept", "iron chelation"],
    },
    "hemophilia a": {
        "omim": "OMIM:306700",
        "gene": "F8",
        "inheritance": "x_linked_recessive",
        "aliases": ["hemophilia", "factor viii deficiency",
                    "classic hemophilia"],
        "workflows": [DiagnosticWorkflowType.GENE_THERAPY_ELIGIBILITY,
                       DiagnosticWorkflowType.THERAPY_SELECTION],
        "search_terms": ["F8", "factor VIII", "valoctocogene", "emicizumab",
                         "fitusiran", "inhibitor development", "Bethesda assay"],
    },

    # -- Immunodeficiency Disorders ------------------------------------
    "severe combined immunodeficiency": {
        "omim": "OMIM:602450",
        "gene": "IL2RG",
        "inheritance": "x_linked_recessive",
        "aliases": ["scid", "bubble boy disease", "x-linked scid",
                    "x-scid"],
        "workflows": [DiagnosticWorkflowType.NEWBORN_SCREENING,
                       DiagnosticWorkflowType.GENE_THERAPY_ELIGIBILITY],
        "search_terms": ["IL2RG", "ADA", "JAK3", "TREC", "lymphopenia",
                         "HSCT", "gene therapy", "adenosine deaminase"],
    },

    # -- Pulmonary Disorders -------------------------------------------
    "cystic fibrosis": {
        "omim": "OMIM:219700",
        "gene": "CFTR",
        "inheritance": "autosomal_recessive",
        "aliases": ["cf", "mucoviscidosis"],
        "workflows": [DiagnosticWorkflowType.NEWBORN_SCREENING,
                       DiagnosticWorkflowType.THERAPY_SELECTION],
        "search_terms": ["CFTR", "F508del", "sweat chloride", "trikafta",
                         "elexacaftor", "ivacaftor", "pancreatic insufficiency"],
    },

    # -- Skeletal Disorders --------------------------------------------
    "achondroplasia": {
        "omim": "OMIM:100800",
        "gene": "FGFR3",
        "inheritance": "autosomal_dominant",
        "aliases": ["achondroplasia", "acp"],
        "workflows": [DiagnosticWorkflowType.PHENOTYPE_DRIVEN,
                       DiagnosticWorkflowType.THERAPY_SELECTION],
        "search_terms": ["FGFR3", "G380R", "vosoritide", "rhizomelic shortening",
                         "foramen magnum stenosis", "CNP analog"],
    },

    # -- Dermatologic Disorders ----------------------------------------
    "epidermolysis bullosa": {
        "omim": "OMIM:226700",
        "gene": "COL7A1",
        "inheritance": "autosomal_recessive",
        "aliases": ["eb", "dystrophic eb", "junctional eb",
                    "eb simplex"],
        "workflows": [DiagnosticWorkflowType.GENE_THERAPY_ELIGIBILITY,
                       DiagnosticWorkflowType.THERAPY_SELECTION],
        "search_terms": ["COL7A1", "LAMB3", "KRT5", "beremagene geperpavec",
                         "blistering", "wound healing", "skin fragility"],
    },

    # -- Neurometabolic Disorders --------------------------------------
    "adrenoleukodystrophy": {
        "omim": "OMIM:300100",
        "gene": "ABCD1",
        "inheritance": "x_linked_recessive",
        "aliases": ["ald", "x-ald", "x-linked adrenoleukodystrophy",
                    "cerebral ald"],
        "workflows": [DiagnosticWorkflowType.NEWBORN_SCREENING,
                       DiagnosticWorkflowType.GENE_THERAPY_ELIGIBILITY],
        "search_terms": ["ABCD1", "VLCFA", "elivaldogene", "Lorenzo's oil",
                         "HSCT", "adrenal insufficiency", "demyelination"],
    },
    "krabbe disease": {
        "omim": "OMIM:245200",
        "gene": "GALC",
        "inheritance": "autosomal_recessive",
        "aliases": ["krabbe", "globoid cell leukodystrophy",
                    "galactosylceramidase deficiency"],
        "workflows": [DiagnosticWorkflowType.NEWBORN_SCREENING,
                       DiagnosticWorkflowType.GENE_THERAPY_ELIGIBILITY],
        "search_terms": ["GALC", "galactosylceramidase", "psychosine",
                         "HSCT", "irritability", "spasticity", "optic atrophy"],
    },

    # -- Ciliopathies --------------------------------------------------
    "polycystic kidney disease": {
        "omim": "OMIM:173900",
        "gene": "PKD1",
        "inheritance": "autosomal_dominant",
        "aliases": ["pkd", "adpkd", "autosomal dominant pkd"],
        "workflows": [DiagnosticWorkflowType.NATURAL_HISTORY,
                       DiagnosticWorkflowType.THERAPY_SELECTION],
        "search_terms": ["PKD1", "PKD2", "polycystin", "tolvaptan",
                         "kidney volume", "hypertension", "cyst growth"],
    },

    # -- Cardiac Genetic Disorders -------------------------------------
    "long qt syndrome": {
        "omim": "OMIM:192500",
        "gene": "KCNQ1",
        "inheritance": "autosomal_dominant",
        "aliases": ["lqts", "long qt", "lqt1", "lqt2", "lqt3",
                    "romano-ward syndrome"],
        "workflows": [DiagnosticWorkflowType.VARIANT_INTERPRETATION,
                       DiagnosticWorkflowType.GENETIC_COUNSELING],
        "search_terms": ["KCNQ1", "KCNH2", "SCN5A", "QTc prolongation",
                         "syncope", "beta-blocker", "ICD", "Schwartz score"],
    },
    "hypertrophic cardiomyopathy": {
        "omim": "OMIM:192600",
        "gene": "MYH7",
        "inheritance": "autosomal_dominant",
        "aliases": ["hcm", "hypertrophic cm", "hocm",
                    "familial hypertrophic cardiomyopathy"],
        "workflows": [DiagnosticWorkflowType.VARIANT_INTERPRETATION,
                       DiagnosticWorkflowType.GENETIC_COUNSELING],
        "search_terms": ["MYH7", "MYBPC3", "mavacamten", "septal hypertrophy",
                         "LVOT obstruction", "sudden cardiac death", "cascade testing"],
    },

    # -- Metabolic Disorders (expanded) -----------------------------------
    "homocystinuria": {
        "omim": "OMIM:236200",
        "gene": "CBS",
        "inheritance": "autosomal_recessive",
        "aliases": ["cystathionine beta-synthase deficiency", "cbs deficiency",
                    "classic homocystinuria"],
        "workflows": [DiagnosticWorkflowType.NEWBORN_SCREENING,
                       DiagnosticWorkflowType.METABOLIC_WORKUP],
        "search_terms": ["CBS", "cystathionine beta-synthase", "homocysteine",
                         "methionine", "lens subluxation", "pyridoxine responsive",
                         "betaine", "thromboembolism"],
    },
    "biotinidase deficiency": {
        "omim": "OMIM:253260",
        "gene": "BTD",
        "inheritance": "autosomal_recessive",
        "aliases": ["biotinidase def", "late-onset multiple carboxylase deficiency"],
        "workflows": [DiagnosticWorkflowType.NEWBORN_SCREENING,
                       DiagnosticWorkflowType.METABOLIC_WORKUP],
        "search_terms": ["BTD", "biotinidase", "biotin", "carboxylase",
                         "seizures", "hearing loss", "alopecia", "skin rash"],
    },
    "glutaric aciduria type 1": {
        "omim": "OMIM:231670",
        "gene": "GCDH",
        "inheritance": "autosomal_recessive",
        "aliases": ["glutaric acidemia type 1", "ga1", "gai",
                    "glutaryl-coa dehydrogenase deficiency"],
        "workflows": [DiagnosticWorkflowType.NEWBORN_SCREENING,
                       DiagnosticWorkflowType.METABOLIC_WORKUP],
        "search_terms": ["GCDH", "glutaryl-CoA dehydrogenase", "glutaric acid",
                         "macrocephaly", "striatal necrosis", "dystonia",
                         "encephalopathic crisis", "carnitine supplementation"],
    },
    "tyrosinemia type 1": {
        "omim": "OMIM:276700",
        "gene": "FAH",
        "inheritance": "autosomal_recessive",
        "aliases": ["tyrosinemia", "hepatorenal tyrosinemia",
                    "fumarylacetoacetate hydrolase deficiency"],
        "workflows": [DiagnosticWorkflowType.NEWBORN_SCREENING,
                       DiagnosticWorkflowType.THERAPY_SELECTION],
        "search_terms": ["FAH", "fumarylacetoacetate hydrolase", "nitisinone",
                         "succinylacetone", "liver failure", "renal tubular dysfunction",
                         "hepatocellular carcinoma", "NTBC"],
    },

    # -- Neurogenetic Disorders (expanded) --------------------------------
    "cdkl5 deficiency disorder": {
        "omim": "OMIM:300672",
        "gene": "CDKL5",
        "inheritance": "x_linked_dominant",
        "aliases": ["cdkl5 deficiency", "cdkl5 epileptic encephalopathy",
                    "early infantile epileptic encephalopathy type 2"],
        "workflows": [DiagnosticWorkflowType.PHENOTYPE_DRIVEN,
                       DiagnosticWorkflowType.GENE_THERAPY_ELIGIBILITY],
        "search_terms": ["CDKL5", "cyclin-dependent kinase-like 5", "ganaxolone",
                         "infantile spasms", "cortical visual impairment",
                         "hand stereotypies", "Rett-like"],
    },
    "canavan disease": {
        "omim": "OMIM:271900",
        "gene": "ASPA",
        "inheritance": "autosomal_recessive",
        "aliases": ["canavan", "aspartoacylase deficiency",
                    "spongy degeneration of the brain"],
        "workflows": [DiagnosticWorkflowType.CARRIER_SCREENING,
                       DiagnosticWorkflowType.GENE_THERAPY_ELIGIBILITY],
        "search_terms": ["ASPA", "aspartoacylase", "N-acetylaspartic acid",
                         "NAA", "macrocephaly", "leukodystrophy",
                         "Ashkenazi Jewish", "white matter disease"],
    },
    "alexander disease": {
        "omim": "OMIM:203450",
        "gene": "GFAP",
        "inheritance": "autosomal_dominant",
        "aliases": ["alexander", "fibrinoid leukodystrophy",
                    "dysmyelinogenic leukodystrophy"],
        "workflows": [DiagnosticWorkflowType.DIFFERENTIAL_DIAGNOSIS,
                       DiagnosticWorkflowType.PHENOTYPE_DRIVEN],
        "search_terms": ["GFAP", "glial fibrillary acidic protein", "Rosenthal fibers",
                         "macrocephaly", "frontal predominant leukodystrophy",
                         "seizures", "bulbar dysfunction"],
    },
    "neuronal ceroid lipofuscinosis": {
        "omim": "OMIM:204200",
        "gene": "CLN3",
        "inheritance": "autosomal_recessive",
        "aliases": ["ncl", "batten disease", "cln3 disease",
                    "juvenile neuronal ceroid lipofuscinosis"],
        "workflows": [DiagnosticWorkflowType.DIFFERENTIAL_DIAGNOSIS,
                       DiagnosticWorkflowType.THERAPY_SELECTION],
        "search_terms": ["CLN3", "CLN2", "TPP1", "cerliponase alfa", "Brineura",
                         "vision loss", "seizures", "motor decline",
                         "curvilinear profiles", "lipofuscin"],
    },

    # -- Hematologic Disorders (expanded) ---------------------------------
    "pyruvate kinase deficiency": {
        "omim": "OMIM:266200",
        "gene": "PKLR",
        "inheritance": "autosomal_recessive",
        "aliases": ["pk deficiency", "pyruvate kinase deficiency of red cells"],
        "workflows": [DiagnosticWorkflowType.DIFFERENTIAL_DIAGNOSIS,
                       DiagnosticWorkflowType.THERAPY_SELECTION],
        "search_terms": ["PKLR", "pyruvate kinase", "mitapivat", "hemolytic anemia",
                         "reticulocytosis", "iron overload", "splenectomy",
                         "chronic hemolysis"],
    },
    "congenital dyserythropoietic anemia": {
        "omim": "OMIM:224120",
        "gene": "SEC23B",
        "inheritance": "autosomal_recessive",
        "aliases": ["cda", "cda type ii", "hempas"],
        "workflows": [DiagnosticWorkflowType.DIFFERENTIAL_DIAGNOSIS,
                       DiagnosticWorkflowType.PHENOTYPE_DRIVEN],
        "search_terms": ["SEC23B", "CDAN1", "KLF1", "dyserythropoiesis",
                         "ineffective erythropoiesis", "iron overload",
                         "binucleated erythroblasts", "jaundice"],
    },

    # -- Immunodeficiency Disorders (expanded) ----------------------------
    "dock8 deficiency": {
        "omim": "OMIM:243700",
        "gene": "DOCK8",
        "inheritance": "autosomal_recessive",
        "aliases": ["dock8 immunodeficiency", "autosomal recessive hyper-ige syndrome",
                    "dock8 deficient hies"],
        "workflows": [DiagnosticWorkflowType.DIFFERENTIAL_DIAGNOSIS,
                       DiagnosticWorkflowType.GENE_THERAPY_ELIGIBILITY],
        "search_terms": ["DOCK8", "dedicator of cytokinesis 8", "hyper-IgE",
                         "eczema", "recurrent viral infections", "skin abscesses",
                         "HSCT", "T cell lymphopenia"],
    },
    "activated pi3k delta syndrome": {
        "omim": "OMIM:615513",
        "gene": "PIK3CD",
        "inheritance": "autosomal_dominant",
        "aliases": ["apds", "apds1", "pasli disease",
                    "p110 delta activating mutation"],
        "workflows": [DiagnosticWorkflowType.THERAPY_SELECTION,
                       DiagnosticWorkflowType.DIFFERENTIAL_DIAGNOSIS],
        "search_terms": ["PIK3CD", "PIK3R1", "PI3K delta", "leniolisib",
                         "lymphoproliferation", "immunodeficiency",
                         "recurrent sinopulmonary infections", "lymphoma risk"],
    },

    # -- Dermatologic / Renal / Syndromic (expanded) ----------------------
    "dystrophic epidermolysis bullosa": {
        "omim": "OMIM:226600",
        "gene": "COL7A1",
        "inheritance": "autosomal_recessive",
        "aliases": ["deb", "rdeb", "recessive dystrophic eb",
                    "dominant dystrophic eb"],
        "workflows": [DiagnosticWorkflowType.GENE_THERAPY_ELIGIBILITY,
                       DiagnosticWorkflowType.THERAPY_SELECTION],
        "search_terms": ["COL7A1", "collagen VII", "beremagene geperpavec", "Vyjuvek",
                         "blistering", "esophageal stricture", "mitten deformity",
                         "squamous cell carcinoma risk", "wound care"],
    },
    "alport syndrome": {
        "omim": "OMIM:301050",
        "gene": "COL4A5",
        "inheritance": "x_linked",
        "aliases": ["alport", "hereditary nephritis",
                    "progressive hereditary nephritis"],
        "workflows": [DiagnosticWorkflowType.NATURAL_HISTORY,
                       DiagnosticWorkflowType.DIFFERENTIAL_DIAGNOSIS],
        "search_terms": ["COL4A5", "COL4A3", "COL4A4", "collagen IV",
                         "hematuria", "sensorineural hearing loss",
                         "lenticonus", "ESRD", "ACE inhibitor"],
    },
    "bardet-biedl syndrome": {
        "omim": "OMIM:209900",
        "gene": "BBS1",
        "inheritance": "autosomal_recessive",
        "aliases": ["bbs", "laurence-moon-bardet-biedl syndrome"],
        "workflows": [DiagnosticWorkflowType.PHENOTYPE_DRIVEN,
                       DiagnosticWorkflowType.DIFFERENTIAL_DIAGNOSIS],
        "search_terms": ["BBS1", "BBS10", "BBSome", "setmelanotide",
                         "retinitis pigmentosa", "obesity", "polydactyly",
                         "renal anomalies", "intellectual disability", "ciliopathy"],
    },
}


# =====================================================================
# RARE DISEASE GENES
# =====================================================================

RARE_DISEASE_GENES: Dict[str, Dict[str, str]] = {
    "CFTR": {
        "full_name": "Cystic Fibrosis Transmembrane Conductance Regulator",
        "omim_gene": "602421",
        "disease": "Cystic fibrosis",
        "inheritance": "autosomal_recessive",
        "significance": "Most common AR lethal disease in Caucasians; CFTR modulator therapy available",
        "workflows": "newborn_screening,therapy_selection",
    },
    "SMN1": {
        "full_name": "Survival Motor Neuron 1",
        "omim_gene": "600354",
        "disease": "Spinal muscular atrophy",
        "inheritance": "autosomal_recessive",
        "significance": "Gene therapy (onasemnogene), ASO (nusinersen), SMN2 modifier (risdiplam) available",
        "workflows": "gene_therapy_eligibility,newborn_screening",
    },
    "DMD": {
        "full_name": "Dystrophin",
        "omim_gene": "300377",
        "disease": "Duchenne/Becker muscular dystrophy",
        "inheritance": "x_linked_recessive",
        "significance": "Exon skipping therapies; gene therapy under development; reading frame determines severity",
        "workflows": "gene_therapy_eligibility,variant_interpretation",
    },
    "GBA1": {
        "full_name": "Glucosylceramidase Beta 1",
        "omim_gene": "606463",
        "disease": "Gaucher disease",
        "inheritance": "autosomal_recessive",
        "significance": "ERT and SRT available; GBA1 variants also risk factor for Parkinson disease",
        "workflows": "therapy_selection,metabolic_workup",
    },
    "HBB": {
        "full_name": "Hemoglobin Subunit Beta",
        "omim_gene": "141900",
        "disease": "Sickle cell disease, Beta-thalassemia",
        "inheritance": "autosomal_recessive",
        "significance": "Gene therapy and gene editing (exagamglogene) approved; most common hemoglobinopathy",
        "workflows": "gene_therapy_eligibility,newborn_screening",
    },
    "MECP2": {
        "full_name": "Methyl-CpG Binding Protein 2",
        "omim_gene": "300005",
        "disease": "Rett syndrome",
        "inheritance": "x_linked_dominant",
        "significance": "Trofinetide approved; gene therapy in clinical trials; almost exclusively affects females",
        "workflows": "gene_therapy_eligibility,phenotype_driven",
    },
    "FMR1": {
        "full_name": "Fragile X Messenger Ribonucleoprotein 1",
        "omim_gene": "309550",
        "disease": "Fragile X syndrome",
        "inheritance": "x_linked_dominant",
        "significance": "Most common inherited cause of intellectual disability; premutation causes FXTAS and FXPOI",
        "workflows": "carrier_screening,genetic_counseling",
    },
    "HTT": {
        "full_name": "Huntingtin",
        "omim_gene": "613004",
        "disease": "Huntington disease",
        "inheritance": "autosomal_dominant",
        "significance": "CAG repeat expansion; anticipation; presymptomatic testing ethical considerations",
        "workflows": "genetic_counseling,natural_history",
    },
    "PAH": {
        "full_name": "Phenylalanine Hydroxylase",
        "omim_gene": "612349",
        "disease": "Phenylketonuria",
        "inheritance": "autosomal_recessive",
        "significance": "Paradigm for newborn screening; BH4-responsive variants; pegvaliase for adults",
        "workflows": "newborn_screening,metabolic_workup",
    },
    "F8": {
        "full_name": "Coagulation Factor VIII",
        "omim_gene": "300841",
        "disease": "Hemophilia A",
        "inheritance": "x_linked_recessive",
        "significance": "Gene therapy (valoctocogene) approved; bispecific antibody (emicizumab) for inhibitor patients",
        "workflows": "gene_therapy_eligibility,therapy_selection",
    },
    "FGFR3": {
        "full_name": "Fibroblast Growth Factor Receptor 3",
        "omim_gene": "134934",
        "disease": "Achondroplasia, Thanatophoric dysplasia",
        "inheritance": "autosomal_dominant",
        "significance": "Vosoritide (CNP analog) approved; most common skeletal dysplasia; G380R hotspot",
        "workflows": "phenotype_driven,therapy_selection",
    },
    "FBN1": {
        "full_name": "Fibrillin 1",
        "omim_gene": "134797",
        "disease": "Marfan syndrome",
        "inheritance": "autosomal_dominant",
        "significance": "Ghent nosology criteria; aortic surveillance critical; losartan/atenolol management",
        "workflows": "phenotype_driven,genetic_counseling",
    },
    "COL7A1": {
        "full_name": "Collagen Type VII Alpha 1 Chain",
        "omim_gene": "120120",
        "disease": "Dystrophic epidermolysis bullosa",
        "inheritance": "autosomal_recessive",
        "significance": "Gene therapy (beremagene geperpavec) approved for topical application",
        "workflows": "gene_therapy_eligibility,therapy_selection",
    },
    "TSC1": {
        "full_name": "TSC Complex Subunit 1 (Hamartin)",
        "omim_gene": "605284",
        "disease": "Tuberous sclerosis complex",
        "inheritance": "autosomal_dominant",
        "significance": "mTOR inhibitors (everolimus) for SEGA and renal AML; high de novo rate",
        "workflows": "phenotype_driven,therapy_selection",
    },
    "TSC2": {
        "full_name": "TSC Complex Subunit 2 (Tuberin)",
        "omim_gene": "191092",
        "disease": "Tuberous sclerosis complex",
        "inheritance": "autosomal_dominant",
        "significance": "More severe phenotype than TSC1; contiguous gene syndrome with PKD1",
        "workflows": "phenotype_driven,therapy_selection",
    },
    "COL1A1": {
        "full_name": "Collagen Type I Alpha 1 Chain",
        "omim_gene": "120150",
        "disease": "Osteogenesis imperfecta",
        "inheritance": "autosomal_dominant",
        "significance": "Null alleles milder (quantitative); glycine substitutions more severe (qualitative)",
        "workflows": "phenotype_driven,variant_interpretation",
    },
    "IDUA": {
        "full_name": "Alpha-L-Iduronidase",
        "omim_gene": "252800",
        "disease": "MPS I (Hurler/Scheie)",
        "inheritance": "autosomal_recessive",
        "significance": "ERT (laronidase) and HSCT available; severity spectrum from Hurler to Scheie",
        "workflows": "gene_therapy_eligibility,newborn_screening",
    },
    "ABCD1": {
        "full_name": "ATP Binding Cassette Subfamily D Member 1",
        "omim_gene": "300371",
        "disease": "X-linked adrenoleukodystrophy",
        "inheritance": "x_linked_recessive",
        "significance": "Gene therapy (elivaldogene) approved; HSCT for cerebral ALD; NBS with VLCFA",
        "workflows": "newborn_screening,gene_therapy_eligibility",
    },
    "GLA": {
        "full_name": "Galactosidase Alpha",
        "omim_gene": "300644",
        "disease": "Fabry disease",
        "inheritance": "x_linked",
        "significance": "ERT (agalsidase) and chaperone (migalastat) available; late-onset variants common",
        "workflows": "therapy_selection,carrier_screening",
    },
    "GAA": {
        "full_name": "Acid Alpha-Glucosidase",
        "omim_gene": "606800",
        "disease": "Pompe disease",
        "inheritance": "autosomal_recessive",
        "significance": "ERT (avalglucosidase alfa) available; CRIM status affects ERT response",
        "workflows": "newborn_screening,therapy_selection",
    },
    "GALC": {
        "full_name": "Galactosylceramidase",
        "omim_gene": "606890",
        "disease": "Krabbe disease",
        "inheritance": "autosomal_recessive",
        "significance": "Pre-symptomatic HSCT is time-critical; NBS enables early intervention",
        "workflows": "newborn_screening,gene_therapy_eligibility",
    },
    "PKD1": {
        "full_name": "Polycystin 1",
        "omim_gene": "601313",
        "disease": "Autosomal dominant polycystic kidney disease",
        "inheritance": "autosomal_dominant",
        "significance": "Tolvaptan slows progression; PKD1 truncating more severe than PKD2",
        "workflows": "natural_history,therapy_selection",
    },
    "KCNQ1": {
        "full_name": "Potassium Voltage-Gated Channel Subfamily Q Member 1",
        "omim_gene": "607542",
        "disease": "Long QT syndrome type 1",
        "inheritance": "autosomal_dominant",
        "significance": "Exercise-triggered arrhythmia; beta-blockers highly effective; homozygous = Jervell-Lange-Nielsen",
        "workflows": "variant_interpretation,genetic_counseling",
    },
    "MYH7": {
        "full_name": "Myosin Heavy Chain 7",
        "omim_gene": "160760",
        "disease": "Hypertrophic cardiomyopathy",
        "inheritance": "autosomal_dominant",
        "significance": "Mavacamten approved; genotype-negative HCM has better prognosis; cascade testing essential",
        "workflows": "variant_interpretation,genetic_counseling",
    },
    "MYBPC3": {
        "full_name": "Myosin Binding Protein C3",
        "omim_gene": "600958",
        "disease": "Hypertrophic cardiomyopathy",
        "inheritance": "autosomal_dominant",
        "significance": "Most common HCM gene; truncating variants; late-onset penetrance common",
        "workflows": "variant_interpretation,genetic_counseling",
    },
    "IL2RG": {
        "full_name": "Interleukin 2 Receptor Subunit Gamma",
        "omim_gene": "308380",
        "disease": "X-linked severe combined immunodeficiency",
        "inheritance": "x_linked_recessive",
        "significance": "Most common SCID type; gene therapy successful; NBS with TREC assay",
        "workflows": "newborn_screening,gene_therapy_eligibility",
    },
    "NPC1": {
        "full_name": "NPC Intracellular Cholesterol Transporter 1",
        "omim_gene": "607623",
        "disease": "Niemann-Pick disease type C",
        "inheritance": "autosomal_recessive",
        "significance": "Miglustat (off-label) only treatment; arimoclomol failed Phase 3; diagnostic delay common",
        "workflows": "differential_diagnosis,metabolic_workup",
    },
    "OTC": {
        "full_name": "Ornithine Transcarbamylase",
        "omim_gene": "300461",
        "disease": "OTC deficiency",
        "inheritance": "x_linked",
        "significance": "Most common urea cycle disorder; neonatal hyperammonemia; carrier females may decompensate",
        "workflows": "metabolic_workup,gene_therapy_eligibility",
    },
    "PMP22": {
        "full_name": "Peripheral Myelin Protein 22",
        "omim_gene": "601097",
        "disease": "Charcot-Marie-Tooth type 1A",
        "inheritance": "autosomal_dominant",
        "significance": "1.5 Mb duplication on 17p11.2; most common inherited neuropathy",
        "workflows": "variant_interpretation,differential_diagnosis",
    },
    "DMPK": {
        "full_name": "DM1 Protein Kinase",
        "omim_gene": "605377",
        "disease": "Myotonic dystrophy type 1",
        "inheritance": "autosomal_dominant",
        "significance": "CTG repeat expansion; anticipation through maternal transmission; multisystem disease",
        "workflows": "genetic_counseling,natural_history",
    },
    "CBS": {
        "full_name": "Cystathionine Beta-Synthase",
        "omim_gene": "613381",
        "disease": "Homocystinuria",
        "inheritance": "autosomal_recessive",
        "significance": "Pyridoxine-responsive variants have better prognosis; betaine and methionine-restricted diet",
        "workflows": "newborn_screening,metabolic_workup",
    },
    "BTD": {
        "full_name": "Biotinidase",
        "omim_gene": "609019",
        "disease": "Biotinidase deficiency",
        "inheritance": "autosomal_recessive",
        "significance": "Treatable IEM on NBS; lifelong oral biotin prevents symptoms; profound vs partial deficiency",
        "workflows": "newborn_screening,metabolic_workup",
    },
    "GCDH": {
        "full_name": "Glutaryl-CoA Dehydrogenase",
        "omim_gene": "608801",
        "disease": "Glutaric aciduria type 1",
        "inheritance": "autosomal_recessive",
        "significance": "Encephalopathic crises triggered by illness; macrocephaly as early sign; NBS with C5DC",
        "workflows": "newborn_screening,metabolic_workup",
    },
    "FAH": {
        "full_name": "Fumarylacetoacetate Hydrolase",
        "omim_gene": "613871",
        "disease": "Tyrosinemia type 1",
        "inheritance": "autosomal_recessive",
        "significance": "Nitisinone (NTBC) transformed prognosis; succinylacetone on NBS is pathognomonic",
        "workflows": "newborn_screening,therapy_selection",
    },
    "CDKL5": {
        "full_name": "Cyclin-Dependent Kinase-Like 5",
        "omim_gene": "300203",
        "disease": "CDKL5 deficiency disorder",
        "inheritance": "x_linked_dominant",
        "significance": "Ganaxolone approved for seizures; distinct from Rett syndrome; early-onset epilepsy",
        "workflows": "phenotype_driven,gene_therapy_eligibility",
    },
    "ASPA": {
        "full_name": "Aspartoacylase",
        "omim_gene": "608034",
        "disease": "Canavan disease",
        "inheritance": "autosomal_recessive",
        "significance": "Elevated NAA on MRS is diagnostic; Ashkenazi Jewish carrier screening; gene therapy in trials",
        "workflows": "carrier_screening,gene_therapy_eligibility",
    },
    "GFAP": {
        "full_name": "Glial Fibrillary Acidic Protein",
        "omim_gene": "137780",
        "disease": "Alexander disease",
        "inheritance": "autosomal_dominant",
        "significance": "De novo dominant mutations; Rosenthal fibers on biopsy; frontal-predominant leukodystrophy",
        "workflows": "differential_diagnosis,phenotype_driven",
    },
    "CLN3": {
        "full_name": "CLN3 Lysosomal/Endosomal Transmembrane Protein",
        "omim_gene": "607042",
        "disease": "Neuronal ceroid lipofuscinosis (juvenile)",
        "inheritance": "autosomal_recessive",
        "significance": "1 kb common deletion accounts for most cases; vision loss then seizures then cognitive decline",
        "workflows": "differential_diagnosis,therapy_selection",
    },
    "PKLR": {
        "full_name": "Pyruvate Kinase L/R",
        "omim_gene": "609712",
        "disease": "Pyruvate kinase deficiency",
        "inheritance": "autosomal_recessive",
        "significance": "Mitapivat (PKR activator) approved; most common enzyme defect of RBC glycolysis",
        "workflows": "differential_diagnosis,therapy_selection",
    },
    "ADAMTS13": {
        "full_name": "ADAM Metallopeptidase with Thrombospondin Type 1 Motif 13",
        "omim_gene": "604134",
        "disease": "Thrombotic thrombocytopenic purpura (congenital)",
        "inheritance": "autosomal_recessive",
        "significance": "Severe ADAMTS13 deficiency (<10%); caplacizumab for acute TTP; recombinant ADAMTS13 in development",
        "workflows": "differential_diagnosis,therapy_selection",
    },
    "DOCK8": {
        "full_name": "Dedicator of Cytokinesis 8",
        "omim_gene": "611432",
        "disease": "DOCK8 immunodeficiency syndrome",
        "inheritance": "autosomal_recessive",
        "significance": "Curative HSCT recommended; severe viral skin infections; distinguished from STAT3-HIES",
        "workflows": "differential_diagnosis,gene_therapy_eligibility",
    },
    "PIK3CD": {
        "full_name": "Phosphatidylinositol-4,5-Bisphosphate 3-Kinase Catalytic Subunit Delta",
        "omim_gene": "602839",
        "disease": "Activated PI3K delta syndrome (APDS)",
        "inheritance": "autosomal_dominant",
        "significance": "Leniolisib approved as first targeted therapy for APDS; gain-of-function mutations",
        "workflows": "therapy_selection,differential_diagnosis",
    },
    "CTLA4": {
        "full_name": "Cytotoxic T-Lymphocyte Associated Protein 4",
        "omim_gene": "123890",
        "disease": "CTLA-4 haploinsufficiency with autoimmune infiltration",
        "inheritance": "autosomal_dominant",
        "significance": "Abatacept (CTLA4-Ig) effective; variable penetrance; immune dysregulation and lymphocytic infiltration",
        "workflows": "therapy_selection,differential_diagnosis",
    },
    "COL4A3": {
        "full_name": "Collagen Type IV Alpha 3 Chain",
        "omim_gene": "120070",
        "disease": "Alport syndrome (autosomal form)",
        "inheritance": "autosomal_recessive",
        "significance": "With COL4A4 causes autosomal Alport; thin basement membrane nephropathy in heterozygotes; ACE inhibitors slow progression",
        "workflows": "natural_history,differential_diagnosis",
    },
    "BBS1": {
        "full_name": "Bardet-Biedl Syndrome 1",
        "omim_gene": "209901",
        "disease": "Bardet-Biedl syndrome",
        "inheritance": "autosomal_recessive",
        "significance": "Most common BBS gene; M390R founder mutation; ciliopathy with retinal, renal, and obesity phenotype; setmelanotide for obesity",
        "workflows": "phenotype_driven,differential_diagnosis",
    },
}


# =====================================================================
# RARE DISEASE PHENOTYPE PATTERNS
# =====================================================================

RARE_DISEASE_PHENOTYPES: Dict[str, Dict[str, object]] = {
    "hypotonia": {
        "hpo_id": "HP:0001252",
        "aliases": ["floppy infant", "decreased muscle tone", "low tone"],
        "associated_diseases": ["SMA", "Pompe disease", "Prader-Willi syndrome",
                                "Down syndrome", "congenital myopathy"],
        "workflows": [DiagnosticWorkflowType.PHENOTYPE_DRIVEN,
                       DiagnosticWorkflowType.DIFFERENTIAL_DIAGNOSIS],
    },
    "seizures": {
        "hpo_id": "HP:0001250",
        "aliases": ["epilepsy", "convulsions", "fits",
                    "epileptic encephalopathy"],
        "associated_diseases": ["Dravet syndrome", "Tuberous sclerosis",
                                "Rett syndrome", "Angelman syndrome",
                                "pyridoxine-dependent epilepsy"],
        "workflows": [DiagnosticWorkflowType.PHENOTYPE_DRIVEN,
                       DiagnosticWorkflowType.DIFFERENTIAL_DIAGNOSIS],
    },
    "intellectual disability": {
        "hpo_id": "HP:0001249",
        "aliases": ["developmental delay", "cognitive impairment",
                    "global developmental delay", "learning disability"],
        "associated_diseases": ["Fragile X syndrome", "Down syndrome",
                                "Rett syndrome", "PKU untreated",
                                "mucopolysaccharidoses"],
        "workflows": [DiagnosticWorkflowType.PHENOTYPE_DRIVEN,
                       DiagnosticWorkflowType.DIFFERENTIAL_DIAGNOSIS],
    },
    "hepatosplenomegaly": {
        "hpo_id": "HP:0001433",
        "aliases": ["enlarged liver and spleen", "hepatomegaly with splenomegaly"],
        "associated_diseases": ["Gaucher disease", "Niemann-Pick disease",
                                "MPS disorders", "Wolman disease",
                                "glycogen storage diseases"],
        "workflows": [DiagnosticWorkflowType.METABOLIC_WORKUP,
                       DiagnosticWorkflowType.DIFFERENTIAL_DIAGNOSIS],
    },
    "hyperammonemia": {
        "hpo_id": "HP:0001987",
        "aliases": ["elevated ammonia", "high ammonia", "ammonia elevation"],
        "associated_diseases": ["OTC deficiency", "CPS1 deficiency",
                                "citrullinemia", "argininosuccinate lyase deficiency",
                                "organic acidemias"],
        "workflows": [DiagnosticWorkflowType.METABOLIC_WORKUP,
                       DiagnosticWorkflowType.NEWBORN_SCREENING],
    },
    "short stature": {
        "hpo_id": "HP:0004322",
        "aliases": ["growth failure", "dwarfism", "growth retardation",
                    "proportionate short stature", "disproportionate short stature"],
        "associated_diseases": ["Achondroplasia", "Turner syndrome",
                                "Noonan syndrome", "Silver-Russell syndrome",
                                "MPS disorders"],
        "workflows": [DiagnosticWorkflowType.PHENOTYPE_DRIVEN,
                       DiagnosticWorkflowType.DIFFERENTIAL_DIAGNOSIS],
    },
    "cardiomyopathy": {
        "hpo_id": "HP:0001638",
        "aliases": ["dilated cardiomyopathy", "hypertrophic cardiomyopathy",
                    "restrictive cardiomyopathy", "cardiac hypertrophy"],
        "associated_diseases": ["Pompe disease", "Fabry disease",
                                "Friedreich ataxia", "Barth syndrome",
                                "Danon disease"],
        "workflows": [DiagnosticWorkflowType.PHENOTYPE_DRIVEN,
                       DiagnosticWorkflowType.DIFFERENTIAL_DIAGNOSIS],
    },
    "corneal clouding": {
        "hpo_id": "HP:0007957",
        "aliases": ["corneal opacity", "cloudy cornea"],
        "associated_diseases": ["MPS I", "MPS VI", "Fabry disease",
                                "cystinosis", "mucolipidosis"],
        "workflows": [DiagnosticWorkflowType.PHENOTYPE_DRIVEN,
                       DiagnosticWorkflowType.METABOLIC_WORKUP],
    },
    "coarse facial features": {
        "hpo_id": "HP:0000280",
        "aliases": ["coarse facies", "gargoylism", "coarsened features"],
        "associated_diseases": ["MPS disorders", "mucolipidosis",
                                "mannosidosis", "fucosidosis",
                                "GM1 gangliosidosis"],
        "workflows": [DiagnosticWorkflowType.PHENOTYPE_DRIVEN,
                       DiagnosticWorkflowType.METABOLIC_WORKUP],
    },
    "failure to thrive": {
        "hpo_id": "HP:0001508",
        "aliases": ["poor weight gain", "faltering growth",
                    "feeding difficulties"],
        "associated_diseases": ["Cystic fibrosis", "celiac disease",
                                "organic acidemias", "mitochondrial disorders",
                                "Noonan syndrome"],
        "workflows": [DiagnosticWorkflowType.PHENOTYPE_DRIVEN,
                       DiagnosticWorkflowType.METABOLIC_WORKUP],
    },
    "dystonia": {
        "hpo_id": "HP:0001332",
        "aliases": ["dystonic posturing", "torsion dystonia"],
        "associated_diseases": ["Glutaric aciduria type 1", "Wilson disease",
                                "DYT1 dystonia", "Lesch-Nyhan syndrome",
                                "neurodegeneration with brain iron accumulation"],
        "workflows": [DiagnosticWorkflowType.DIFFERENTIAL_DIAGNOSIS,
                       DiagnosticWorkflowType.PHENOTYPE_DRIVEN],
    },
    "ataxia": {
        "hpo_id": "HP:0001251",
        "aliases": ["cerebellar ataxia", "gait ataxia", "truncal ataxia",
                    "spinocerebellar ataxia"],
        "associated_diseases": ["Friedreich ataxia", "Ataxia-telangiectasia",
                                "Niemann-Pick C", "SCA types",
                                "mitochondrial disorders"],
        "workflows": [DiagnosticWorkflowType.DIFFERENTIAL_DIAGNOSIS,
                       DiagnosticWorkflowType.PHENOTYPE_DRIVEN],
    },
    "skin blistering": {
        "hpo_id": "HP:0008066",
        "aliases": ["blistering", "bullae", "vesicles",
                    "skin fragility"],
        "associated_diseases": ["Epidermolysis bullosa subtypes",
                                "pemphigus", "dermatitis herpetiformis"],
        "workflows": [DiagnosticWorkflowType.PHENOTYPE_DRIVEN,
                       DiagnosticWorkflowType.DIFFERENTIAL_DIAGNOSIS],
    },
    "metabolic acidosis": {
        "hpo_id": "HP:0001942",
        "aliases": ["acidemia", "low bicarbonate", "anion gap acidosis"],
        "associated_diseases": ["Organic acidemias", "MSUD",
                                "mitochondrial disorders", "renal tubular acidosis",
                                "lactic acidosis"],
        "workflows": [DiagnosticWorkflowType.METABOLIC_WORKUP,
                       DiagnosticWorkflowType.NEWBORN_SCREENING],
    },
    "hypoglycemia": {
        "hpo_id": "HP:0001943",
        "aliases": ["low blood sugar", "hypoketotic hypoglycemia",
                    "ketotic hypoglycemia"],
        "associated_diseases": ["MCADD", "Glycogen storage diseases",
                                "hyperinsulinism", "fatty acid oxidation defects",
                                "Beckwith-Wiedemann syndrome"],
        "workflows": [DiagnosticWorkflowType.METABOLIC_WORKUP,
                       DiagnosticWorkflowType.NEWBORN_SCREENING],
    },
    "lens subluxation": {
        "hpo_id": "HP:0001083",
        "aliases": ["ectopia lentis", "dislocated lens", "lens dislocation"],
        "associated_diseases": ["Marfan syndrome", "Homocystinuria",
                                "Weill-Marchesani syndrome",
                                "sulfite oxidase deficiency"],
        "workflows": [DiagnosticWorkflowType.PHENOTYPE_DRIVEN,
                       DiagnosticWorkflowType.DIFFERENTIAL_DIAGNOSIS],
    },
    "recurrent infections": {
        "hpo_id": "HP:0002719",
        "aliases": ["frequent infections", "immunodeficiency",
                    "susceptibility to infections"],
        "associated_diseases": ["SCID", "Common variable immunodeficiency",
                                "Chronic granulomatous disease",
                                "Wiskott-Aldrich syndrome",
                                "Cystic fibrosis"],
        "workflows": [DiagnosticWorkflowType.DIFFERENTIAL_DIAGNOSIS,
                       DiagnosticWorkflowType.PHENOTYPE_DRIVEN],
    },
    "hearing loss": {
        "hpo_id": "HP:0000365",
        "aliases": ["deafness", "sensorineural hearing loss",
                    "conductive hearing loss"],
        "associated_diseases": ["Connexin 26 (GJB2)", "Usher syndrome",
                                "Pendred syndrome", "Waardenburg syndrome",
                                "Alport syndrome"],
        "workflows": [DiagnosticWorkflowType.DIFFERENTIAL_DIAGNOSIS,
                       DiagnosticWorkflowType.PHENOTYPE_DRIVEN],
    },
    "retinitis pigmentosa": {
        "hpo_id": "HP:0000510",
        "aliases": ["rod-cone dystrophy", "night blindness",
                    "progressive vision loss"],
        "associated_diseases": ["Usher syndrome", "Bardet-Biedl syndrome",
                                "Refsum disease", "Leber congenital amaurosis",
                                "Stargardt disease"],
        "workflows": [DiagnosticWorkflowType.GENE_THERAPY_ELIGIBILITY,
                       DiagnosticWorkflowType.DIFFERENTIAL_DIAGNOSIS],
    },
    "aortic root dilation": {
        "hpo_id": "HP:0002616",
        "aliases": ["aortic aneurysm", "aortic root enlargement",
                    "aortic dilation"],
        "associated_diseases": ["Marfan syndrome", "Loeys-Dietz syndrome",
                                "Vascular EDS", "Turner syndrome",
                                "bicuspid aortic valve"],
        "workflows": [DiagnosticWorkflowType.PHENOTYPE_DRIVEN,
                       DiagnosticWorkflowType.GENETIC_COUNSELING],
    },
    "cherry red spot": {
        "hpo_id": "HP:0010729",
        "aliases": ["macular cherry red spot", "cherry-red macula"],
        "associated_diseases": ["Tay-Sachs disease", "Niemann-Pick type A",
                                "GM1 gangliosidosis", "Sandhoff disease",
                                "sialidosis"],
        "workflows": [DiagnosticWorkflowType.PHENOTYPE_DRIVEN,
                       DiagnosticWorkflowType.METABOLIC_WORKUP],
    },
    "lactic acidosis": {
        "hpo_id": "HP:0003128",
        "aliases": ["elevated lactate", "hyperlactatemia",
                    "lactic acid elevation"],
        "associated_diseases": ["Mitochondrial disorders", "MELAS",
                                "Leigh syndrome", "pyruvate dehydrogenase deficiency",
                                "respiratory chain defects"],
        "workflows": [DiagnosticWorkflowType.METABOLIC_WORKUP,
                       DiagnosticWorkflowType.DIFFERENTIAL_DIAGNOSIS],
    },
    "photosensitivity": {
        "hpo_id": "HP:0000992",
        "aliases": ["sun sensitivity", "cutaneous photosensitivity",
                    "light sensitivity", "UV sensitivity"],
        "associated_diseases": ["Xeroderma pigmentosum", "Cockayne syndrome",
                                "erythropoietic protoporphyria",
                                "Bloom syndrome", "Rothmund-Thomson syndrome"],
        "workflows": [DiagnosticWorkflowType.PHENOTYPE_DRIVEN,
                       DiagnosticWorkflowType.DIFFERENTIAL_DIAGNOSIS],
    },
    "sensorineural hearing loss": {
        "hpo_id": "HP:0000407",
        "aliases": ["snhl", "nerve deafness", "cochlear hearing loss",
                    "bilateral sensorineural hearing loss"],
        "associated_diseases": ["Usher syndrome", "Pendred syndrome",
                                "mitochondrial hearing loss (m.1555A>G)",
                                "Waardenburg syndrome", "Alport syndrome"],
        "workflows": [DiagnosticWorkflowType.DIFFERENTIAL_DIAGNOSIS,
                       DiagnosticWorkflowType.PHENOTYPE_DRIVEN],
    },
    "progressive external ophthalmoplegia": {
        "hpo_id": "HP:0000590",
        "aliases": ["peo", "chronic progressive external ophthalmoplegia",
                    "cpeo", "ophthalmoparesis"],
        "associated_diseases": ["Mitochondrial DNA deletions", "POLG mutations",
                                "Kearns-Sayre syndrome", "ANT1 deficiency",
                                "Twinkle helicase deficiency"],
        "workflows": [DiagnosticWorkflowType.DIFFERENTIAL_DIAGNOSIS,
                       DiagnosticWorkflowType.PHENOTYPE_DRIVEN],
    },
    "angiokeratoma": {
        "hpo_id": "HP:0001014",
        "aliases": ["angiokeratoma corporis diffusum",
                    "angiokeratomas", "vascular skin lesion"],
        "associated_diseases": ["Fabry disease", "fucosidosis",
                                "beta-mannosidosis", "galactosialidosis",
                                "Schindler disease"],
        "workflows": [DiagnosticWorkflowType.PHENOTYPE_DRIVEN,
                       DiagnosticWorkflowType.METABOLIC_WORKUP],
    },
    "microcephaly": {
        "hpo_id": "HP:0000252",
        "aliases": ["small head circumference", "reduced head circumference",
                    "primary microcephaly", "secondary microcephaly"],
        "associated_diseases": ["Rett syndrome", "Angelman syndrome",
                                "Smith-Lemli-Opitz syndrome", "Seckel syndrome",
                                "congenital Zika syndrome"],
        "workflows": [DiagnosticWorkflowType.PHENOTYPE_DRIVEN,
                       DiagnosticWorkflowType.DIFFERENTIAL_DIAGNOSIS],
    },
    "macrocephaly": {
        "hpo_id": "HP:0000256",
        "aliases": ["large head", "megalencephaly", "increased head circumference",
                    "relative macrocephaly"],
        "associated_diseases": ["Canavan disease", "Glutaric aciduria type 1",
                                "Alexander disease", "Sotos syndrome",
                                "PTEN hamartoma syndrome"],
        "workflows": [DiagnosticWorkflowType.PHENOTYPE_DRIVEN,
                       DiagnosticWorkflowType.DIFFERENTIAL_DIAGNOSIS],
    },
    "developmental regression": {
        "hpo_id": "HP:0002376",
        "aliases": ["regression", "loss of skills", "neurodevelopmental regression",
                    "psychomotor regression"],
        "associated_diseases": ["Rett syndrome", "neuronal ceroid lipofuscinosis",
                                "Krabbe disease", "metachromatic leukodystrophy",
                                "Sanfilippo syndrome"],
        "workflows": [DiagnosticWorkflowType.DIFFERENTIAL_DIAGNOSIS,
                       DiagnosticWorkflowType.PHENOTYPE_DRIVEN],
    },
    "spasticity": {
        "hpo_id": "HP:0001257",
        "aliases": ["spastic paraplegia", "upper motor neuron signs",
                    "increased muscle tone", "spastic gait"],
        "associated_diseases": ["Hereditary spastic paraplegia", "Krabbe disease",
                                "adrenoleukodystrophy", "metachromatic leukodystrophy",
                                "Canavan disease"],
        "workflows": [DiagnosticWorkflowType.DIFFERENTIAL_DIAGNOSIS,
                       DiagnosticWorkflowType.PHENOTYPE_DRIVEN],
    },
}


# =====================================================================
# SEARCH PLAN DATACLASS
# =====================================================================

@dataclass
class SearchPlan:
    """Agent's plan for answering a rare disease diagnostic question.

    The search plan captures all entities detected in the user's question
    and the strategy the agent will use to retrieve evidence from the
    14 rare-disease-specific Milvus collections.
    """
    question: str
    conditions: List[str] = field(default_factory=list)
    genes: List[str] = field(default_factory=list)
    phenotypes: List[str] = field(default_factory=list)
    relevant_workflows: List[DiagnosticWorkflowType] = field(default_factory=list)
    search_strategy: str = "broad"  # broad, targeted, phenotype_match, variant_focused
    sub_questions: List[str] = field(default_factory=list)
    identified_topics: List[str] = field(default_factory=list)


# =====================================================================
# RARE DISEASE DIAGNOSTIC AGENT
# =====================================================================

class RareDiseaseAgent:
    """Autonomous Rare Disease Diagnostic Agent.

    Wraps the multi-collection RareDiseaseRAGEngine with planning and reasoning
    capabilities. Designed to answer complex questions about rare disease
    diagnosis, variant interpretation, therapeutic options, and genetic
    counseling.

    Example queries this agent handles:
    - "Infant with hypotonia, hepatomegaly, and cardiomyopathy -- differential diagnosis?"
    - "Interpret NM_000492.4(CFTR):c.1521_1523del (p.Phe508del) using ACMG criteria"
    - "Is this patient eligible for gene therapy for SMA type 1?"
    - "Positive newborn screen for MCADD -- what is the workup?"
    - "Carrier risk for Tay-Sachs in Ashkenazi Jewish couple"
    - "Natural history of Fabry disease in heterozygous females"
    - "Available clinical trials for Duchenne muscular dystrophy gene therapy"
    - "Metabolic crisis workup for infant with hyperammonemia"

    Usage:
        agent = RareDiseaseAgent(rag_engine)
        plan = agent.search_plan("Infant with hypotonia and cardiomyopathy")
        response = agent.run("Infant with hypotonia and cardiomyopathy")
    """

    def __init__(self, rag_engine):
        """Initialize agent with a configured RAG engine.

        Args:
            rag_engine: RareDiseaseRAGEngine instance with Milvus collections connected.
        """
        self.rag = rag_engine
        self.knowledge = {
            "conditions": RARE_DISEASE_CONDITIONS,
            "genes": RARE_DISEASE_GENES,
            "phenotypes": RARE_DISEASE_PHENOTYPES,
        }

    # -- Public API ----------------------------------------------------

    def run(
        self,
        query: str,
        workflow_type: Optional[DiagnosticWorkflowType] = None,
        patient_context: Optional[dict] = None,
        **kwargs,
    ) -> DiagnosticResult:
        """Execute the full agent pipeline: plan -> search -> evaluate -> synthesize.

        Args:
            query: Natural language question about rare disease diagnostics.
            workflow_type: Optional workflow override for collection boosting.
            patient_context: Optional patient data for personalised diagnosis.
            **kwargs: Additional query parameters (top_k, collection_filter).

        Returns:
            DiagnosticResult with findings, differential diagnoses, and metadata.
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

        Detects rare diseases, genes, and phenotypes in the question text.
        Determines relevant diagnostic workflows, chooses a search strategy,
        and generates sub-questions for comprehensive retrieval across
        collections.

        Args:
            question: The user's natural language question.

        Returns:
            SearchPlan with all detected entities and retrieval strategy.
        """
        plan = SearchPlan(question=question)

        # Step 1: Detect entities
        entities = self._detect_entities(question)
        plan.conditions = entities.get("conditions", [])
        plan.genes = entities.get("genes", [])
        plan.phenotypes = entities.get("phenotypes", [])

        # Step 2: Determine relevant workflows
        plan.relevant_workflows = [self._detect_workflow(question)]
        # Add entity-derived workflows
        for condition in plan.conditions:
            info = RARE_DISEASE_CONDITIONS.get(condition, {})
            for wf in info.get("workflows", []):
                if wf not in plan.relevant_workflows:
                    plan.relevant_workflows.append(wf)

        # Step 3: Choose search strategy
        plan.search_strategy = self._choose_strategy(
            question, plan.conditions, plan.genes, plan.phenotypes,
        )

        # Step 4: Generate sub-questions
        plan.sub_questions = self._generate_sub_questions(plan)

        # Step 5: Compile identified topics
        plan.identified_topics = (
            plan.conditions + plan.genes + plan.phenotypes
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

    def generate_report(
        self,
        results,
        workflow: Optional[DiagnosticWorkflowType] = None,
    ) -> str:
        """Generate a structured rare disease diagnostic report.

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
            "# Rare Disease Diagnostic Intelligence Report",
            f"**Query:** {question}",
            f"**Generated:** {timestamp}",
            f"**Workflows:** {', '.join(wf.value for wf in plan.relevant_workflows)}",
            f"**Strategy:** {plan.search_strategy}",
            "",
        ]

        # Detected entities
        if plan.conditions or plan.genes or plan.phenotypes:
            report_lines.extend([
                "---",
                "",
                "## Detected Clinical Entities",
                "",
            ])
            if plan.conditions:
                report_lines.append(
                    f"- **Conditions:** {', '.join(plan.conditions)}"
                )
            if plan.genes:
                report_lines.append(
                    f"- **Genes:** {', '.join(plan.genes)}"
                )
            if plan.phenotypes:
                report_lines.append(
                    f"- **Phenotypes (HPO):** {', '.join(plan.phenotypes)}"
                )
            report_lines.append("")

        # Critical findings check
        critical_flags = []
        if hasattr(results, "results") and results.results:
            for r in results.results:
                meta = r.metadata if hasattr(r, "metadata") else {}
                if meta.get("urgency") == "critical" or meta.get("actionable"):
                    critical_flags.append(r)

        if critical_flags:
            report_lines.extend([
                "---",
                "",
                "## [CRITICAL] Urgent Diagnostic / Therapeutic Alerts",
                "",
            ])
            for flag in critical_flags:
                text = flag.text if hasattr(flag, "text") else str(flag)
                report_lines.append(
                    f"- **[CRITICAL]** {text[:200]} -- "
                    f"immediate clinical action required."
                )
            report_lines.append("")

        # Analysis section
        report_lines.extend([
            "---",
            "",
            "## Diagnostic Assessment",
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

        # Differential diagnosis section for phenotype-driven workflows
        if workflow in (DiagnosticWorkflowType.PHENOTYPE_DRIVEN,
                        DiagnosticWorkflowType.DIFFERENTIAL_DIAGNOSIS,
                        DiagnosticWorkflowType.METABOLIC_WORKUP):
            report_lines.extend([
                "---",
                "",
                "## Diagnostic Considerations",
                "",
                "- Prioritize diagnoses with available therapies (treatability bias)",
                "- Consider consanguinity and ethnicity in carrier frequency estimation",
                "- Trio whole-exome/genome sequencing recommended for undiagnosed cases",
                "- Re-analysis of negative genomic data at 12-month intervals",
                "",
            ])

        # Gene therapy / clinical trial section for relevant workflows
        if workflow in (DiagnosticWorkflowType.GENE_THERAPY_ELIGIBILITY,
                        DiagnosticWorkflowType.CLINICAL_TRIAL_MATCHING,
                        DiagnosticWorkflowType.THERAPY_SELECTION):
            report_lines.extend([
                "---",
                "",
                "## Therapeutic Landscape Considerations",
                "",
                "- Confirm molecular diagnosis before gene therapy eligibility assessment",
                "- Check pre-existing AAV neutralizing antibody titers for AAV-based therapies",
                "- Review FDA Orphan Drug and Breakthrough Therapy designations",
                "- Consider compassionate use / expanded access programs",
                "",
            ])

        # Genetic counseling note for relevant workflows
        if workflow in (DiagnosticWorkflowType.GENETIC_COUNSELING,
                        DiagnosticWorkflowType.CARRIER_SCREENING,
                        DiagnosticWorkflowType.PRENATAL_DIAGNOSIS):
            report_lines.extend([
                "---",
                "",
                "## Genetic Counseling Considerations",
                "",
                "- Confirm variant segregation with parental testing",
                "- Address recurrence risk based on inheritance pattern",
                "- Discuss reproductive options (PGT-M, prenatal diagnosis)",
                "- Refer to genetic counselor for psychosocial support",
                "- Connect family with disease-specific patient advocacy organizations",
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
            "*This report is generated by the Rare Disease Diagnostic Agent "
            "within the HCLS AI Factory. All diagnostic suggestions require "
            "validation by clinical geneticists, genetic counselors, and "
            "appropriate laboratory testing. Variant interpretations should "
            "be confirmed by CLIA/CAP-certified laboratories. This tool does "
            "NOT provide definitive diagnoses.*",
        ])

        return "\n".join(report_lines)

    # -- Workflow Detection --------------------------------------------

    def _detect_workflow(self, question: str) -> DiagnosticWorkflowType:
        """Detect the most relevant workflow from a question.

        Uses keyword-based heuristics to identify which of the 13 diagnostic
        workflows is most relevant to the query.

        Args:
            question: The user's natural language question.

        Returns:
            Most relevant DiagnosticWorkflowType.
        """
        text_upper = question.upper()

        workflow_scores: Dict[DiagnosticWorkflowType, float] = {}

        keyword_workflow_map = {
            DiagnosticWorkflowType.PHENOTYPE_DRIVEN: [
                "PHENOTYPE", "HPO", "HP:", "CLINICAL FEATURES",
                "DYSMORPHIC", "SYNDROMIC", "FACIAL FEATURES",
                "PRESENTING WITH", "SYMPTOMS", "SIGNS",
            ],
            DiagnosticWorkflowType.VARIANT_INTERPRETATION: [
                "VARIANT", "MUTATION", "ACMG", "PATHOGENIC",
                "VUS", "BENIGN", "CLASSIFY", "INTERPRETATION",
                "NM_", "C.", "P.", "MISSENSE", "NONSENSE",
                "FRAMESHIFT", "SPLICE", "CLINVAR",
            ],
            DiagnosticWorkflowType.DIFFERENTIAL_DIAGNOSIS: [
                "DIFFERENTIAL", "DIAGNOSIS", "DIAGNOSE",
                "WHAT DISEASE", "WHICH CONDITION", "RULE OUT",
                "CONSIDER", "DDX", "WORKUP",
            ],
            DiagnosticWorkflowType.GENE_THERAPY_ELIGIBILITY: [
                "GENE THERAPY", "AAV", "LENTIVIRAL", "ELIGIBLE",
                "ELIGIBILITY", "ZOLGENSMA", "LUXTURNA",
                "BASE EDITING", "PRIME EDITING", "CRISPR",
                "GENE REPLACEMENT", "GENE TRANSFER",
            ],
            DiagnosticWorkflowType.NEWBORN_SCREENING: [
                "NEWBORN SCREEN", "NBS", "NEWBORN", "DRIED BLOOD SPOT",
                "TANDEM MASS", "ACYLCARNITINE", "RUSP",
                "RECOMMENDED UNIFORM SCREENING PANEL",
                "POSITIVE SCREEN", "SCREENING FOLLOW-UP",
            ],
            DiagnosticWorkflowType.METABOLIC_WORKUP: [
                "METABOLIC", "INBORN ERROR", "IEM", "AMMONIA",
                "LACTIC ACID", "LACTATE", "ORGANIC ACID",
                "AMINO ACID", "ACYLCARNITINE", "METABOLIC CRISIS",
                "ACIDOSIS", "HYPERAMMONEMIA", "HYPOGLYCEMIA",
            ],
            DiagnosticWorkflowType.CARRIER_SCREENING: [
                "CARRIER", "CARRIER SCREENING", "CARRIER STATUS",
                "HETEROZYGOUS", "CARRIER RISK", "CARRIER TESTING",
                "EXPANDED CARRIER SCREEN", "CARRIER FREQUENCY",
            ],
            DiagnosticWorkflowType.PRENATAL_DIAGNOSIS: [
                "PRENATAL", "FETAL", "AMNIOCENTESIS", "CVS",
                "CHORIONIC VILLUS", "PREGNANCY", "PREIMPLANTATION",
                "PGT-M", "NIPT", "CELL-FREE DNA",
            ],
            DiagnosticWorkflowType.NATURAL_HISTORY: [
                "NATURAL HISTORY", "PROGNOSIS", "PROGRESSION",
                "LIFE EXPECTANCY", "DISEASE COURSE", "OUTCOMES",
                "LONG-TERM", "LONGITUDINAL", "SURVIVAL",
            ],
            DiagnosticWorkflowType.THERAPY_SELECTION: [
                "TREATMENT", "THERAPY", "DRUG", "ERT", "SRT",
                "ORPHAN DRUG", "MEDICATION", "ENZYME REPLACEMENT",
                "SUBSTRATE REDUCTION", "ASO", "CHAPERONE",
            ],
            DiagnosticWorkflowType.CLINICAL_TRIAL_MATCHING: [
                "CLINICAL TRIAL", "NCT", "TRIAL", "ENROLL",
                "STUDY", "RECRUITING", "INVESTIGATIONAL",
                "EXPERIMENTAL", "PHASE 1", "PHASE 2", "PHASE 3",
            ],
            DiagnosticWorkflowType.GENETIC_COUNSELING: [
                "COUNSELING", "COUNSELLING", "RECURRENCE RISK",
                "INHERITANCE", "FAMILY", "AUTOSOMAL",
                "X-LINKED", "ANTICIPATION", "PENETRANCE",
                "EXPRESSIVITY", "CASCADE", "PRESYMPTOMATIC",
            ],
        }

        for wf, keywords in keyword_workflow_map.items():
            for kw in keywords:
                if kw in text_upper:
                    workflow_scores[wf] = workflow_scores.get(wf, 0) + 1.0

        if not workflow_scores:
            return DiagnosticWorkflowType.GENERAL

        sorted_workflows = sorted(
            workflow_scores.items(), key=lambda x: x[1], reverse=True,
        )

        return sorted_workflows[0][0]

    # -- Entity Detection ----------------------------------------------

    def _detect_entities(self, question: str) -> Dict[str, List[str]]:
        """Detect rare disease entities in the question text.

        Scans for conditions, genes, and phenotypes using the knowledge
        dictionaries. Performs case-insensitive matching against canonical
        names and aliases.

        Args:
            question: The user's natural language question.

        Returns:
            Dict with keys 'conditions', 'genes', 'phenotypes' mapping
            to lists of detected entity names.
        """
        entities: Dict[str, List[str]] = {
            "conditions": [],
            "genes": [],
            "phenotypes": [],
        }

        text_lower = question.lower()

        # Detect conditions
        for condition, info in RARE_DISEASE_CONDITIONS.items():
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

        # Detect genes (case-sensitive for gene symbols)
        for gene, info in RARE_DISEASE_GENES.items():
            # Gene symbols are uppercase -- search in original text
            pattern = r'\b' + re.escape(gene) + r'\b'
            if re.search(pattern, question):
                if gene not in entities["genes"]:
                    entities["genes"].append(gene)
                continue
            # Also check full name (case-insensitive)
            full_name = info.get("full_name", "")
            if full_name and full_name.lower() in text_lower:
                if gene not in entities["genes"]:
                    entities["genes"].append(gene)

        # Detect phenotypes
        for phenotype, info in RARE_DISEASE_PHENOTYPES.items():
            if phenotype in text_lower:
                if phenotype not in entities["phenotypes"]:
                    entities["phenotypes"].append(phenotype)
                continue
            aliases = info.get("aliases", [])
            for alias in aliases:
                if alias.lower() in text_lower:
                    if phenotype not in entities["phenotypes"]:
                        entities["phenotypes"].append(phenotype)
                    break
            # Check HPO ID
            hpo_id = info.get("hpo_id", "")
            if hpo_id and hpo_id in question:
                if phenotype not in entities["phenotypes"]:
                    entities["phenotypes"].append(phenotype)

        return entities

    # -- Search Strategy -----------------------------------------------

    def _choose_strategy(
        self,
        text: str,
        conditions: List[str],
        genes: List[str],
        phenotypes: List[str],
    ) -> str:
        """Choose search strategy based on query type.

        Args:
            text: Original query text.
            conditions: Detected conditions.
            genes: Detected genes.
            phenotypes: Detected phenotypes.

        Returns:
            Strategy name string.
        """
        text_upper = text.upper()

        # Variant-focused queries
        variant_keywords = [
            "VARIANT", "MUTATION", "NM_", "C.", "P.",
            "ACMG", "PATHOGENIC", "VUS", "CLASSIFY",
        ]
        if any(kw in text_upper for kw in variant_keywords):
            return "variant_focused"

        # Phenotype matching queries (multiple phenotypes, no specific condition)
        if len(phenotypes) >= 2 and len(conditions) == 0:
            return "phenotype_match"

        # Targeted: specific condition or gene
        if (len(conditions) == 1 and len(genes) <= 1) or (
            len(conditions) <= 1 and len(genes) == 1
        ):
            if conditions or genes:
                return "targeted"

        return "broad"

    # -- Sub-Question Generation ----------------------------------------

    def _generate_sub_questions(self, plan: SearchPlan) -> List[str]:
        """Generate sub-questions for comprehensive retrieval.

        Decomposes the main question into focused sub-queries based on
        the detected entities and workflow type. Enables multi-hop
        retrieval across different aspects of the diagnostic question.

        Args:
            plan: SearchPlan with detected entities and workflows.

        Returns:
            List of sub-question strings (typically 2-4 questions).
        """
        sub_questions: List[str] = []

        condition_label = plan.conditions[0] if plan.conditions else "the rare disease"
        gene_label = plan.genes[0] if plan.genes else "the causative gene"
        phenotype_label = plan.phenotypes[0] if plan.phenotypes else "the phenotype"

        primary_wf = (
            plan.relevant_workflows[0]
            if plan.relevant_workflows
            else DiagnosticWorkflowType.GENERAL
        )

        # -- Pattern 1: Phenotype-Driven --------------------------------
        if primary_wf == DiagnosticWorkflowType.PHENOTYPE_DRIVEN:
            sub_questions = [
                f"What rare diseases are associated with {phenotype_label}?",
                f"What genes cause {phenotype_label} in rare disease?",
                f"What is the differential diagnosis for {phenotype_label}?",
                f"What diagnostic tests confirm diseases presenting with {phenotype_label}?",
            ]

        # -- Pattern 2: Variant Interpretation --------------------------
        elif primary_wf == DiagnosticWorkflowType.VARIANT_INTERPRETATION:
            sub_questions = [
                f"What is the ACMG classification of variants in {gene_label}?",
                f"What diseases are caused by pathogenic variants in {gene_label}?",
                f"What is the ClinVar evidence for variants in {gene_label}?",
                f"What functional studies support pathogenicity of {gene_label} variants?",
            ]

        # -- Pattern 3: Differential Diagnosis --------------------------
        elif primary_wf == DiagnosticWorkflowType.DIFFERENTIAL_DIAGNOSIS:
            sub_questions = [
                f"What is the differential diagnosis for {condition_label}?",
                f"What phenotypic features distinguish {condition_label} from similar conditions?",
                f"What genetic tests are needed to confirm {condition_label}?",
                f"What is the prevalence and natural history of {condition_label}?",
            ]

        # -- Pattern 4: Gene Therapy Eligibility ------------------------
        elif primary_wf == DiagnosticWorkflowType.GENE_THERAPY_ELIGIBILITY:
            sub_questions = [
                f"What gene therapies are available or in trials for {condition_label}?",
                f"What are the eligibility criteria for {condition_label} gene therapy?",
                f"What is the safety profile of gene therapy for {condition_label}?",
                f"What outcomes have been reported for {condition_label} gene therapy?",
            ]

        # -- Pattern 5: Newborn Screening -------------------------------
        elif primary_wf == DiagnosticWorkflowType.NEWBORN_SCREENING:
            sub_questions = [
                f"What is the newborn screening follow-up protocol for {condition_label}?",
                f"What confirmatory tests are needed after positive NBS for {condition_label}?",
                f"What is the false positive rate for {condition_label} newborn screening?",
                f"What early interventions improve outcomes in {condition_label}?",
            ]

        # -- Pattern 6: Metabolic Workup --------------------------------
        elif primary_wf == DiagnosticWorkflowType.METABOLIC_WORKUP:
            sub_questions = [
                f"What metabolic conditions present with {phenotype_label}?",
                f"What biochemical tests are diagnostic for {condition_label}?",
                f"What is the acute management of metabolic crisis in {condition_label}?",
                f"What long-term dietary and pharmacologic management exists for {condition_label}?",
            ]

        # -- Pattern 7: Carrier Screening -------------------------------
        elif primary_wf == DiagnosticWorkflowType.CARRIER_SCREENING:
            sub_questions = [
                f"What is the carrier frequency of {condition_label} by ethnicity?",
                f"What testing methodology detects carriers for {gene_label}?",
                f"What reproductive options exist for carriers of {condition_label}?",
                f"What genetic counseling considerations apply to {condition_label} carriers?",
            ]

        # -- Pattern 8: Genetic Counseling ------------------------------
        elif primary_wf == DiagnosticWorkflowType.GENETIC_COUNSELING:
            sub_questions = [
                f"What is the inheritance pattern and recurrence risk for {condition_label}?",
                f"What is the penetrance and expressivity of {gene_label} variants?",
                f"What cascade testing strategy is recommended for {condition_label}?",
                f"What psychosocial considerations apply to {condition_label} families?",
            ]

        # -- Pattern 9: Natural History ---------------------------------
        elif primary_wf == DiagnosticWorkflowType.NATURAL_HISTORY:
            sub_questions = [
                f"What is the natural history and prognosis of {condition_label}?",
                f"What are the major milestones and complications of {condition_label}?",
                f"What surveillance protocol is recommended for {condition_label}?",
                f"What quality of life measures are used in {condition_label}?",
            ]

        # -- Pattern 10: Therapy Selection ------------------------------
        elif primary_wf == DiagnosticWorkflowType.THERAPY_SELECTION:
            sub_questions = [
                f"What approved therapies exist for {condition_label}?",
                f"What is the evidence for treatment efficacy in {condition_label}?",
                f"What monitoring is required during treatment of {condition_label}?",
                f"What new therapies are in development for {condition_label}?",
            ]

        # -- Pattern 11: Clinical Trial Matching ------------------------
        elif primary_wf == DiagnosticWorkflowType.CLINICAL_TRIAL_MATCHING:
            sub_questions = [
                f"What clinical trials are recruiting for {condition_label}?",
                f"What are the eligibility criteria for {condition_label} trials?",
                f"What therapeutic modalities are being tested for {condition_label}?",
                f"What trial sites are available for {condition_label} studies?",
            ]

        # -- Default ----------------------------------------------------
        else:
            sub_questions = [
                f"What is known about {condition_label} in the rare disease literature?",
                f"What genes and variants are associated with {condition_label}?",
                f"What therapeutic options exist for {condition_label}?",
            ]

        return sub_questions

    # -- Build Search Strategy -----------------------------------------

    def _build_search_strategy(
        self,
        entities: Dict[str, List[str]],
        workflow: DiagnosticWorkflowType,
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
        if entities.get("genes"):
            parts.append(f"Genes: {', '.join(entities['genes'])}")
        if entities.get("phenotypes"):
            parts.append(f"Phenotypes: {', '.join(entities['phenotypes'])}")

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
