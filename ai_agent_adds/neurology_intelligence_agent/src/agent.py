"""Neurology Intelligence Agent -- autonomous reasoning across neurological data.

Implements the plan -> search -> evaluate -> synthesize -> report pattern from the
VAST AI OS AgentEngine model. The agent can:

1. Parse complex multi-part questions about neurological conditions and management
2. Plan a search strategy across 14 domain-specific collections
3. Execute multi-collection retrieval via the NeuroRAGEngine
4. Evaluate evidence quality and completeness
5. Synthesize cross-functional insights with clinical alerts
6. Generate structured reports with neurology-specific formatting

Mapping to VAST AI OS:
  - AgentEngine entry point: NeurologyAgent.run()
  - Plan -> search_plan()
  - Execute -> rag_engine.query()
  - Reflect -> evaluate_evidence()
  - Report -> generate_report()

Domain coverage:
  - Cerebrovascular disease (stroke, TIA, SAH, cavernous malformations)
  - Neurodegenerative disorders (Alzheimer's, Parkinson's, FTD, DLB, ALS, MSA, PSP, CBD)
  - Epilepsy and seizure disorders (focal, generalized, ILAE 2017 classification)
  - Multiple sclerosis and neuroimmunology (RRMS, PPMS, SPMS, NMOSD, MOGAD)
  - Brain tumors and neuro-oncology (glioblastoma, meningioma, CNS lymphoma)
  - Movement disorders (PD, ET, dystonia, Huntington's, atypical parkinsonism)
  - Headache and migraine (episodic, chronic, cluster, medication overuse)
  - Neuromuscular disease (MG, ALS, GBS, CIDP, SMA, muscular dystrophies)

Clinical scales computed: NIHSS, GCS, MoCA, UPDRS, EDSS, mRS, HIT-6, ALSFRS-R
Guideline bodies referenced: AAN, EAN, ILAE, MDS, IHS

Author: Adam Jones
Date: March 2026
"""

from dataclasses import dataclass, field
from datetime import datetime, timezone
from enum import Enum
from typing import Dict, List, Optional


# =====================================================================
# ENUMS
# =====================================================================

class NeuroWorkflowType(str, Enum):
    """Types of neurology query workflows."""
    STROKE_ACUTE = "stroke_acute"
    STROKE_PREVENTION = "stroke_prevention"
    DEMENTIA_EVALUATION = "dementia_evaluation"
    EPILEPSY_CLASSIFICATION = "epilepsy_classification"
    MS_MANAGEMENT = "ms_management"
    MOVEMENT_DISORDER = "movement_disorder"
    HEADACHE_DIAGNOSIS = "headache_diagnosis"
    NEUROMUSCULAR_EVAL = "neuromuscular_eval"
    NEURO_ONCOLOGY = "neuro_oncology"
    GENERAL = "general"


class EvidenceLevel(str, Enum):
    """Clinical evidence hierarchy for neurological recommendations."""
    CLASS_I = "I"          # Prospective, randomized, controlled clinical trial
    CLASS_II = "II"        # Prospective matched group cohort study
    CLASS_III = "III"      # Retrospective cohort, case-control, cross-sectional
    CLASS_IV = "IV"        # Case series, case reports, expert opinion
    GUIDELINE = "guideline"  # Published guideline recommendation (AAN, EAN, ILAE, MDS)
    META_ANALYSIS = "meta"   # Systematic review / meta-analysis


class SeverityLevel(str, Enum):
    """Finding severity classification."""
    CRITICAL = "critical"
    HIGH = "high"
    MODERATE = "moderate"
    LOW = "low"
    INFORMATIONAL = "informational"


class ImagingModality(str, Enum):
    """Neuroimaging modalities."""
    MRI_BRAIN = "mri_brain"
    MRI_SPINE = "mri_spine"
    CT_HEAD = "ct_head"
    CT_ANGIOGRAPHY = "ct_angiography"
    MR_ANGIOGRAPHY = "mr_angiography"
    PET_FDG = "pet_fdg"
    PET_AMYLOID = "pet_amyloid"
    PET_TAU = "pet_tau"
    DAT_SPECT = "dat_spect"
    CONVENTIONAL_ANGIOGRAPHY = "conventional_angiography"


class ElectrophysiologyType(str, Enum):
    """Electrophysiology study types."""
    EEG_ROUTINE = "eeg_routine"
    EEG_CONTINUOUS = "eeg_continuous"
    EEG_VIDEO = "eeg_video"
    EMG = "emg"
    NCS = "ncs"
    EVOKED_POTENTIALS = "evoked_potentials"
    POLYSOMNOGRAPHY = "polysomnography"


# =====================================================================
# RESPONSE DATACLASS
# =====================================================================

@dataclass
class NeuroResponse:
    """Complete response from the neurology intelligence agent.

    Attributes:
        question: Original query text.
        answer: LLM-synthesised answer text.
        results: Ranked search results used for synthesis.
        workflow: Neuro workflow that was applied.
        confidence: Overall confidence score (0.0 - 1.0).
        citations: Formatted citation list.
        search_time_ms: Total search time in milliseconds.
        collections_searched: Number of collections queried.
        patient_context_used: Whether patient context was injected.
        clinical_alerts: Any critical clinical findings flagged.
        timestamp: ISO 8601 timestamp of response generation.
    """
    question: str = ""
    answer: str = ""
    results: list = field(default_factory=list)
    workflow: Optional[NeuroWorkflowType] = None
    confidence: float = 0.0
    citations: List[Dict[str, str]] = field(default_factory=list)
    search_time_ms: float = 0.0
    collections_searched: int = 0
    patient_context_used: bool = False
    clinical_alerts: List[str] = field(default_factory=list)
    timestamp: str = field(
        default_factory=lambda: datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    )


# =====================================================================
# NEUROLOGY SYSTEM PROMPT
# =====================================================================

NEURO_SYSTEM_PROMPT = """\
You are a neurology clinical decision support system within the HCLS AI Factory. \
You have deep expertise in cerebrovascular disease, neurodegenerative disorders, \
epilepsy, multiple sclerosis, brain tumors, movement disorders, headache/migraine, \
and neuromuscular disease. You analyze neuroimaging (MRI, CT, PET, DaT-SPECT), \
electrophysiology (EEG, EMG, NCS), CSF biomarkers (A\u03b242, p-tau181, NfL), \
and neurogenetic data. You calculate validated clinical scales (NIHSS, GCS, MoCA, \
UPDRS, EDSS, mRS, HIT-6, ALSFRS-R). You reference AAN, EAN, ILAE, and MDS \
guidelines. You identify stroke tPA/thrombectomy eligibility, classify seizures \
per ILAE 2017, stage dementia using NIA-AA ATN framework, and detect atypical \
parkinsonism red flags. Never fabricate clinical data.

Your responses must adhere to the following standards:

1. **Clinical Citations** -- Always cite clinical evidence using PubMed identifiers \
   with clickable links: [PMID:12345678](https://pubmed.ncbi.nlm.nih.gov/12345678/). \
   Include study type (RCT, cohort, meta-analysis), sample size, and key endpoints \
   when available. For clinical trials, cite NCT identifiers: \
   [NCT01234567](https://clinicaltrials.gov/study/NCT01234567).

2. **Guideline References** -- Cite AAN practice parameters and practice advisories \
   with evidence classification (Class I-IV) and recommendation level (Level A-U). \
   Reference EAN guidelines, ILAE classification and treatment guidelines, MDS \
   diagnostic criteria, and IHS classification (ICHD-3). Format as: \
   [AAN 2019] Evidence-based guideline: Treatment of convulsive status epilepticus \
   (Level A recommendation, Class I evidence).

3. **CRITICAL Findings** -- Flag the following as CRITICAL with prominent visual \
   markers and immediate action recommendations:
   - Acute stroke within tPA window (<4.5 hours) or thrombectomy window (<24 hours)
   - Status epilepticus (convulsive or non-convulsive)
   - Acute intracranial hemorrhage or expanding hematoma
   - Signs of herniation (pupil asymmetry, Cushing triad, posturing)
   - Acute cord compression with progressive neurological deficit
   - Guillain-Barre syndrome with respiratory compromise (FVC <20 mL/kg)
   - Myasthenic crisis with respiratory failure
   - Malignant cerebral edema requiring decompressive craniectomy evaluation
   - Acute hydrocephalus requiring emergent CSF diversion
   - Brain death evaluation criteria met

4. **Severity Badges** -- Classify all findings using standardised severity levels: \
   [CRITICAL], [HIGH], [MODERATE], [LOW], [INFORMATIONAL]. Place the badge at the \
   start of each finding or recommendation line.

5. **Clinical Scale Scoring** -- When calculating neurological assessment scales, \
   show the component-level breakdown and total score with interpretation:
   - NIHSS: 0-42 (0=no deficit, 1-4=minor, 5-15=moderate, 16-20=moderate-severe, 21-42=severe)
   - GCS: 3-15 (Eye 1-4, Verbal 1-5, Motor 1-6; <=8=coma, 9-12=moderate, 13-15=mild)
   - MoCA: 0-30 (>=26=normal, 18-25=mild cognitive impairment, <18=dementia range)
   - UPDRS Part III: 0-132 (motor examination; higher=worse)
   - EDSS: 0-10 (0=normal, 1-3.5=mild, 4-5.5=moderate, 6-7.5=severe, 8-9.5=restricted)
   - mRS: 0-6 (0=no symptoms, 1-2=favorable outcome, 3-5=unfavorable, 6=dead)
   - HIT-6: 36-78 (<=49=little/no impact, 50-55=some, 56-59=substantial, 60+=severe)
   - ALSFRS-R: 0-48 (higher=better function; rate of decline predicts prognosis)

6. **Structured Formatting** -- Organise responses with clear sections: \
   Clinical Summary, Diagnostic Findings, Differential Diagnosis, Imaging \
   Interpretation, Electrophysiology, Biomarker Analysis, Treatment \
   Recommendations, Follow-up Plan, and Guideline References. Use bullet \
   points and numbered lists for actionable items.

7. **Genomic Cross-Reference** -- When neurogenetic conditions are relevant \
   (e.g., APOE genotyping in AD, GBA/LRRK2 in PD, SCN1A in Dravet syndrome, \
   NOTCH3 in CADASIL, HTT in Huntington's, SOD1/C9orf72 in ALS, SMN1 in SMA), \
   cross-reference with the genomic_evidence collection and discuss genetic \
   counseling implications and targeted therapy eligibility.

8. **Neuroimaging Interpretation** -- For imaging-related queries, specify the \
   modality, sequence (T1, T2, FLAIR, DWI, ADC, SWI, GRE, MRA, CTA), relevant \
   findings, and differential diagnosis. Reference standardized scoring systems \
   (ASPECTS for stroke, McDonald criteria lesion requirements for MS, Fazekas \
   scale for white matter disease, Scheltens visual rating for atrophy).

9. **Treatment Decision Trees** -- For treatment questions, provide structured \
   decision trees considering: first-line vs. second-line therapy, mechanism \
   of action, contraindications, drug interactions, monitoring requirements, \
   insurance/access considerations, and evidence quality for each option. \
   Reference AAN/EAN treatment guidelines and FDA-approved indications.

10. **Limitations** -- You are a neurology clinical decision support tool. You \
    do NOT replace board-certified neurologists, neuroradiologists, or \
    neurosurgeons. All recommendations require review by qualified clinicians \
    with access to the full clinical picture, examination findings, and \
    imaging studies. Explicitly state when evidence is limited, when urgent \
    in-person evaluation is needed, or when specialist consultation \
    (epileptology, neuro-oncology, neuromuscular, movement disorders) \
    is recommended."""


# =====================================================================
# WORKFLOW-SPECIFIC COLLECTION BOOST WEIGHTS
# =====================================================================
# Maps each NeuroWorkflowType to collection weight overrides (multipliers).
# Collections not listed retain their base weight (1.0x). Values > 1.0
# boost the collection; values < 1.0 would suppress it.

WORKFLOW_COLLECTION_BOOST: Dict[NeuroWorkflowType, Dict[str, float]] = {

    # -- Acute Stroke ---------------------------------------------------
    NeuroWorkflowType.STROKE_ACUTE: {
        "neuro_cerebrovascular": 2.5,
        "neuro_imaging": 2.0,
        "neuro_guidelines": 1.8,
        "neuro_trials": 1.5,
        "neuro_literature": 1.3,
        "neuro_electrophysiology": 0.8,
        "neuro_degenerative": 0.5,
        "neuro_epilepsy": 0.7,
        "neuro_oncology": 0.5,
        "neuro_ms": 0.5,
        "neuro_movement": 0.5,
        "neuro_headache": 0.7,
        "neuro_neuromuscular": 0.5,
        "genomic_evidence": 0.8,
    },

    # -- Stroke Prevention -----------------------------------------------
    NeuroWorkflowType.STROKE_PREVENTION: {
        "neuro_cerebrovascular": 2.5,
        "neuro_guidelines": 2.0,
        "neuro_trials": 1.8,
        "neuro_literature": 1.5,
        "neuro_imaging": 1.3,
        "genomic_evidence": 1.2,
        "neuro_electrophysiology": 0.8,
        "neuro_degenerative": 0.6,
        "neuro_epilepsy": 0.5,
        "neuro_oncology": 0.5,
        "neuro_ms": 0.5,
        "neuro_movement": 0.6,
        "neuro_headache": 0.8,
        "neuro_neuromuscular": 0.5,
    },

    # -- Dementia Evaluation ---------------------------------------------
    NeuroWorkflowType.DEMENTIA_EVALUATION: {
        "neuro_degenerative": 2.5,
        "neuro_imaging": 2.0,
        "neuro_guidelines": 1.8,
        "neuro_literature": 1.5,
        "genomic_evidence": 1.5,
        "neuro_trials": 1.3,
        "neuro_electrophysiology": 1.0,
        "neuro_movement": 1.2,
        "neuro_cerebrovascular": 1.0,
        "neuro_epilepsy": 0.7,
        "neuro_oncology": 0.6,
        "neuro_ms": 0.6,
        "neuro_headache": 0.5,
        "neuro_neuromuscular": 0.5,
    },

    # -- Epilepsy Classification -----------------------------------------
    NeuroWorkflowType.EPILEPSY_CLASSIFICATION: {
        "neuro_epilepsy": 2.5,
        "neuro_electrophysiology": 2.0,
        "neuro_imaging": 1.8,
        "neuro_guidelines": 1.8,
        "neuro_literature": 1.3,
        "neuro_trials": 1.2,
        "genomic_evidence": 1.5,
        "neuro_degenerative": 0.7,
        "neuro_cerebrovascular": 0.8,
        "neuro_oncology": 0.8,
        "neuro_ms": 0.5,
        "neuro_movement": 0.5,
        "neuro_headache": 0.5,
        "neuro_neuromuscular": 0.5,
    },

    # -- MS Management ---------------------------------------------------
    NeuroWorkflowType.MS_MANAGEMENT: {
        "neuro_ms": 2.5,
        "neuro_imaging": 2.0,
        "neuro_guidelines": 1.8,
        "neuro_trials": 1.5,
        "neuro_literature": 1.5,
        "neuro_electrophysiology": 1.2,
        "genomic_evidence": 1.0,
        "neuro_neuromuscular": 0.8,
        "neuro_cerebrovascular": 0.7,
        "neuro_degenerative": 0.7,
        "neuro_epilepsy": 0.5,
        "neuro_oncology": 0.5,
        "neuro_movement": 0.5,
        "neuro_headache": 0.6,
    },

    # -- Movement Disorder -----------------------------------------------
    NeuroWorkflowType.MOVEMENT_DISORDER: {
        "neuro_movement": 2.5,
        "neuro_imaging": 1.8,
        "neuro_guidelines": 1.8,
        "neuro_degenerative": 1.5,
        "neuro_literature": 1.3,
        "neuro_trials": 1.3,
        "genomic_evidence": 1.5,
        "neuro_electrophysiology": 1.0,
        "neuro_cerebrovascular": 0.7,
        "neuro_epilepsy": 0.6,
        "neuro_oncology": 0.5,
        "neuro_ms": 0.5,
        "neuro_headache": 0.5,
        "neuro_neuromuscular": 0.8,
    },

    # -- Headache Diagnosis -----------------------------------------------
    NeuroWorkflowType.HEADACHE_DIAGNOSIS: {
        "neuro_headache": 2.5,
        "neuro_guidelines": 2.0,
        "neuro_imaging": 1.5,
        "neuro_literature": 1.5,
        "neuro_trials": 1.3,
        "neuro_cerebrovascular": 1.2,
        "neuro_electrophysiology": 0.8,
        "neuro_degenerative": 0.5,
        "neuro_epilepsy": 0.6,
        "neuro_oncology": 0.8,
        "neuro_ms": 0.6,
        "neuro_movement": 0.5,
        "neuro_neuromuscular": 0.5,
        "genomic_evidence": 0.8,
    },

    # -- Neuromuscular Evaluation ----------------------------------------
    NeuroWorkflowType.NEUROMUSCULAR_EVAL: {
        "neuro_neuromuscular": 2.5,
        "neuro_electrophysiology": 2.0,
        "neuro_guidelines": 1.8,
        "neuro_literature": 1.5,
        "neuro_trials": 1.3,
        "genomic_evidence": 1.5,
        "neuro_imaging": 1.2,
        "neuro_degenerative": 1.0,
        "neuro_cerebrovascular": 0.5,
        "neuro_epilepsy": 0.5,
        "neuro_oncology": 0.6,
        "neuro_ms": 0.6,
        "neuro_movement": 0.8,
        "neuro_headache": 0.5,
    },

    # -- Neuro-Oncology --------------------------------------------------
    NeuroWorkflowType.NEURO_ONCOLOGY: {
        "neuro_oncology": 2.5,
        "neuro_imaging": 2.0,
        "neuro_guidelines": 1.5,
        "neuro_trials": 1.8,
        "neuro_literature": 1.5,
        "genomic_evidence": 1.5,
        "neuro_electrophysiology": 0.8,
        "neuro_epilepsy": 1.0,
        "neuro_cerebrovascular": 0.7,
        "neuro_degenerative": 0.5,
        "neuro_ms": 0.5,
        "neuro_movement": 0.5,
        "neuro_headache": 0.8,
        "neuro_neuromuscular": 0.5,
    },

    # -- General (balanced across all collections) -----------------------
    NeuroWorkflowType.GENERAL: {
        "neuro_literature": 1.2,
        "neuro_trials": 1.1,
        "neuro_imaging": 1.1,
        "neuro_electrophysiology": 1.0,
        "neuro_degenerative": 1.0,
        "neuro_cerebrovascular": 1.0,
        "neuro_epilepsy": 1.0,
        "neuro_oncology": 0.9,
        "neuro_ms": 1.0,
        "neuro_movement": 1.0,
        "neuro_headache": 0.9,
        "neuro_neuromuscular": 0.9,
        "neuro_guidelines": 1.2,
        "genomic_evidence": 0.8,
    },
}


# =====================================================================
# KNOWLEDGE DOMAIN DICTIONARIES
# =====================================================================
# Comprehensive neurological knowledge for entity detection and context
# enrichment. Used by the agent's search_plan() to identify entities
# in user queries and map them to workflows.

NEURO_CONDITIONS: Dict[str, Dict[str, object]] = {

    # -- Cerebrovascular Disease ----------------------------------------
    "ischemic stroke": {
        "aliases": ["acute ischemic stroke", "ais", "cerebral infarction",
                    "thromboembolic stroke", "large vessel occlusion", "lvo",
                    "lacunar infarct", "watershed infarct"],
        "workflows": [NeuroWorkflowType.STROKE_ACUTE, NeuroWorkflowType.STROKE_PREVENTION],
        "search_terms": ["tPA", "alteplase", "tenecteplase", "thrombectomy",
                        "NIHSS", "ASPECTS", "penumbra", "DWI-FLAIR mismatch"],
    },
    "hemorrhagic stroke": {
        "aliases": ["intracerebral hemorrhage", "ich", "brain bleed",
                    "hypertensive hemorrhage", "lobar hemorrhage",
                    "intraparenchymal hemorrhage"],
        "workflows": [NeuroWorkflowType.STROKE_ACUTE],
        "search_terms": ["hematoma expansion", "ICH score", "spot sign",
                        "surgical evacuation", "anticoagulant reversal",
                        "blood pressure management"],
    },
    "subarachnoid hemorrhage": {
        "aliases": ["sah", "aneurysmal sah", "ruptured aneurysm",
                    "thunderclap headache"],
        "workflows": [NeuroWorkflowType.STROKE_ACUTE, NeuroWorkflowType.HEADACHE_DIAGNOSIS],
        "search_terms": ["Hunt-Hess", "Fisher grade", "vasospasm", "nimodipine",
                        "coiling", "clipping", "delayed cerebral ischemia",
                        "external ventricular drain"],
    },
    "transient ischemic attack": {
        "aliases": ["tia", "mini stroke", "transient neurological deficit"],
        "workflows": [NeuroWorkflowType.STROKE_PREVENTION],
        "search_terms": ["ABCD2 score", "dual antiplatelet", "carotid stenosis",
                        "secondary prevention", "atrial fibrillation"],
    },
    "carotid stenosis": {
        "aliases": ["carotid artery stenosis", "carotid artery disease",
                    "internal carotid stenosis"],
        "workflows": [NeuroWorkflowType.STROKE_PREVENTION],
        "search_terms": ["carotid endarterectomy", "carotid stenting",
                        "NASCET criteria", "duplex ultrasound"],
    },
    "cerebral venous thrombosis": {
        "aliases": ["cvt", "cerebral venous sinus thrombosis", "cvst",
                    "dural sinus thrombosis", "sagittal sinus thrombosis"],
        "workflows": [NeuroWorkflowType.STROKE_ACUTE],
        "search_terms": ["anticoagulation", "MR venography", "headache",
                        "papilledema", "empty delta sign"],
    },

    # -- Neurodegenerative Disorders ------------------------------------
    "alzheimer disease": {
        "aliases": ["alzheimers", "alzheimer's disease", "alzheimer's",
                    "ad", "early onset alzheimer", "familial alzheimer"],
        "workflows": [NeuroWorkflowType.DEMENTIA_EVALUATION],
        "search_terms": ["amyloid beta", "A\u03b242", "p-tau181", "NfL",
                        "MoCA", "CDR-SB", "ADAS-Cog", "ATN framework",
                        "lecanemab", "donanemab", "APOE", "hippocampal atrophy"],
    },
    "mild cognitive impairment": {
        "aliases": ["mci", "amnestic mci", "non-amnestic mci",
                    "prodromal alzheimer"],
        "workflows": [NeuroWorkflowType.DEMENTIA_EVALUATION],
        "search_terms": ["MoCA", "neuropsychological testing", "amyloid PET",
                        "CSF biomarkers", "conversion rate", "cholinesterase inhibitor"],
    },
    "frontotemporal dementia": {
        "aliases": ["ftd", "behavioral variant ftd", "bvftd",
                    "primary progressive aphasia", "ppa",
                    "semantic dementia", "nonfluent aphasia",
                    "frontotemporal lobar degeneration", "ftld"],
        "workflows": [NeuroWorkflowType.DEMENTIA_EVALUATION],
        "search_terms": ["frontal atrophy", "MAPT", "GRN", "C9orf72",
                        "TDP-43", "behavioral disinhibition", "executive dysfunction"],
    },
    "dementia with lewy bodies": {
        "aliases": ["dlb", "lewy body dementia", "lbd"],
        "workflows": [NeuroWorkflowType.DEMENTIA_EVALUATION, NeuroWorkflowType.MOVEMENT_DISORDER],
        "search_terms": ["visual hallucinations", "REM sleep behavior disorder",
                        "parkinsonism", "fluctuating cognition", "DaT-SPECT",
                        "cholinesterase inhibitor", "neuroleptic sensitivity"],
    },
    "vascular dementia": {
        "aliases": ["vascular cognitive impairment", "vci",
                    "multi-infarct dementia", "subcortical vascular dementia"],
        "workflows": [NeuroWorkflowType.DEMENTIA_EVALUATION, NeuroWorkflowType.STROKE_PREVENTION],
        "search_terms": ["white matter hyperintensities", "Fazekas", "lacunar infarcts",
                        "strategic infarct", "cerebral small vessel disease"],
    },

    # -- Parkinsonism & Movement Disorders ------------------------------
    "parkinson disease": {
        "aliases": ["parkinsons", "parkinson's disease", "parkinson's",
                    "pd", "idiopathic parkinson"],
        "workflows": [NeuroWorkflowType.MOVEMENT_DISORDER],
        "search_terms": ["bradykinesia", "rigidity", "resting tremor",
                        "UPDRS", "Hoehn and Yahr", "levodopa",
                        "dopamine agonist", "DaT-SPECT", "alpha-synuclein",
                        "GBA", "LRRK2", "deep brain stimulation"],
    },
    "essential tremor": {
        "aliases": ["et", "benign essential tremor", "familial tremor",
                    "action tremor"],
        "workflows": [NeuroWorkflowType.MOVEMENT_DISORDER],
        "search_terms": ["postural tremor", "kinetic tremor", "propranolol",
                        "primidone", "focused ultrasound", "DBS thalamus"],
    },
    "dystonia": {
        "aliases": ["cervical dystonia", "torticollis", "blepharospasm",
                    "writer's cramp", "generalized dystonia", "dyt1",
                    "dopa-responsive dystonia"],
        "workflows": [NeuroWorkflowType.MOVEMENT_DISORDER],
        "search_terms": ["botulinum toxin", "DBS GPi", "anticholinergic",
                        "TOR1A", "DYT1", "Fahn-Marsden scale"],
    },
    "huntington disease": {
        "aliases": ["huntingtons", "huntington's disease", "huntington's",
                    "hd", "huntington chorea"],
        "workflows": [NeuroWorkflowType.MOVEMENT_DISORDER, NeuroWorkflowType.DEMENTIA_EVALUATION],
        "search_terms": ["HTT", "CAG repeat", "chorea", "tetrabenazine",
                        "deutetrabenazine", "genetic counseling",
                        "UHDRS", "presymptomatic testing"],
    },
    "progressive supranuclear palsy": {
        "aliases": ["psp", "steele-richardson-olszewski syndrome"],
        "workflows": [NeuroWorkflowType.MOVEMENT_DISORDER],
        "search_terms": ["vertical gaze palsy", "postural instability",
                        "hummingbird sign", "midbrain atrophy", "tau",
                        "Richardson syndrome", "falls"],
    },
    "multiple system atrophy": {
        "aliases": ["msa", "msa-p", "msa-c", "shy-drager syndrome",
                    "striatonigral degeneration", "olivopontocerebellar atrophy"],
        "workflows": [NeuroWorkflowType.MOVEMENT_DISORDER],
        "search_terms": ["autonomic failure", "cerebellar ataxia",
                        "orthostatic hypotension", "hot cross bun sign",
                        "alpha-synuclein", "putaminal rim sign"],
    },
    "corticobasal degeneration": {
        "aliases": ["cbd", "corticobasal syndrome", "cbs"],
        "workflows": [NeuroWorkflowType.MOVEMENT_DISORDER],
        "search_terms": ["asymmetric parkinsonism", "alien limb",
                        "cortical sensory loss", "myoclonus", "apraxia",
                        "tau", "asymmetric cortical atrophy"],
    },

    # -- Epilepsy & Seizure Disorders -----------------------------------
    "epilepsy": {
        "aliases": ["seizure disorder", "recurrent seizures"],
        "workflows": [NeuroWorkflowType.EPILEPSY_CLASSIFICATION],
        "search_terms": ["EEG", "antiseizure medication", "ASM",
                        "ILAE classification", "drug-resistant epilepsy",
                        "epilepsy surgery", "vagus nerve stimulator"],
    },
    "focal epilepsy": {
        "aliases": ["partial epilepsy", "temporal lobe epilepsy", "tle",
                    "frontal lobe epilepsy", "focal aware seizure",
                    "focal impaired awareness seizure", "complex partial seizure"],
        "workflows": [NeuroWorkflowType.EPILEPSY_CLASSIFICATION],
        "search_terms": ["focal onset", "temporal sclerosis", "mesial temporal",
                        "hippocampal sclerosis", "carbamazepine", "lamotrigine",
                        "levetiracetam", "epilepsy surgery", "SEEG"],
    },
    "generalized epilepsy": {
        "aliases": ["idiopathic generalized epilepsy", "ige",
                    "juvenile myoclonic epilepsy", "jme",
                    "childhood absence epilepsy", "juvenile absence epilepsy",
                    "generalized tonic-clonic seizure", "gtcs"],
        "workflows": [NeuroWorkflowType.EPILEPSY_CLASSIFICATION],
        "search_terms": ["generalized spike-wave", "valproate", "lamotrigine",
                        "levetiracetam", "ethosuximide", "absence seizure",
                        "myoclonic seizure", "3 Hz spike-wave"],
    },
    "status epilepticus": {
        "aliases": ["se", "convulsive status epilepticus",
                    "nonconvulsive status epilepticus", "ncse",
                    "refractory status epilepticus", "super-refractory se"],
        "workflows": [NeuroWorkflowType.EPILEPSY_CLASSIFICATION],
        "search_terms": ["benzodiazepine", "lorazepam", "midazolam",
                        "fosphenytoin", "levetiracetam", "valproate",
                        "continuous EEG", "burst suppression", "pentobarbital"],
    },
    "dravet syndrome": {
        "aliases": ["severe myoclonic epilepsy of infancy", "smei"],
        "workflows": [NeuroWorkflowType.EPILEPSY_CLASSIFICATION],
        "search_terms": ["SCN1A", "febrile seizure", "fenfluramine",
                        "stiripentol", "cannabidiol", "sodium channel"],
    },

    # -- Multiple Sclerosis & Neuroimmunology ---------------------------
    "multiple sclerosis": {
        "aliases": ["ms", "relapsing-remitting ms", "rrms",
                    "primary progressive ms", "ppms",
                    "secondary progressive ms", "spms",
                    "clinically isolated syndrome", "cis"],
        "workflows": [NeuroWorkflowType.MS_MANAGEMENT],
        "search_terms": ["oligoclonal bands", "McDonald criteria", "EDSS",
                        "FLAIR lesions", "Dawson fingers", "DMT",
                        "disease-modifying therapy", "natalizumab",
                        "ocrelizumab", "ofatumumab", "JCV antibody"],
    },
    "neuromyelitis optica spectrum disorder": {
        "aliases": ["nmosd", "devic disease", "neuromyelitis optica",
                    "nmo", "aqp4 antibody disease"],
        "workflows": [NeuroWorkflowType.MS_MANAGEMENT],
        "search_terms": ["AQP4-IgG", "aquaporin-4", "longitudinally extensive",
                        "transverse myelitis", "optic neuritis",
                        "eculizumab", "inebilizumab", "satralizumab"],
    },
    "mog antibody disease": {
        "aliases": ["mogad", "mog-ad", "mog antibody associated disease",
                    "mog-igg disease"],
        "workflows": [NeuroWorkflowType.MS_MANAGEMENT],
        "search_terms": ["MOG-IgG", "optic neuritis", "ADEM",
                        "transverse myelitis", "cortical encephalitis",
                        "steroid-responsive"],
    },

    # -- Brain Tumors / Neuro-Oncology ----------------------------------
    "glioblastoma": {
        "aliases": ["gbm", "glioblastoma multiforme", "grade 4 glioma",
                    "idh-wildtype glioblastoma", "high-grade glioma"],
        "workflows": [NeuroWorkflowType.NEURO_ONCOLOGY],
        "search_terms": ["MGMT methylation", "temozolomide", "radiation",
                        "IDH mutation", "bevacizumab", "TTFields",
                        "tumor treating fields", "pseudoprogression",
                        "RANO criteria"],
    },
    "low-grade glioma": {
        "aliases": ["lgg", "diffuse astrocytoma", "oligodendroglioma",
                    "grade 2 glioma", "idh-mutant glioma",
                    "idh-mutant astrocytoma"],
        "workflows": [NeuroWorkflowType.NEURO_ONCOLOGY],
        "search_terms": ["IDH1", "IDH2", "1p/19q codeletion", "ATRX",
                        "observation", "early surgery", "radiation timing",
                        "seizure management"],
    },
    "meningioma": {
        "aliases": ["meningeal tumor", "convexity meningioma",
                    "parasagittal meningioma", "skull base meningioma",
                    "atypical meningioma"],
        "workflows": [NeuroWorkflowType.NEURO_ONCOLOGY],
        "search_terms": ["dural tail", "Simpson grade", "WHO grade",
                        "observation", "surgery", "stereotactic radiosurgery",
                        "NF2", "somatostatin receptor"],
    },
    "brain metastases": {
        "aliases": ["cerebral metastases", "brain mets",
                    "metastatic brain tumor", "leptomeningeal disease",
                    "leptomeningeal carcinomatosis"],
        "workflows": [NeuroWorkflowType.NEURO_ONCOLOGY],
        "search_terms": ["SRS", "stereotactic radiosurgery", "WBRT",
                        "GPA score", "immunotherapy", "targeted therapy",
                        "blood-brain barrier", "leptomeningeal"],
    },
    "cns lymphoma": {
        "aliases": ["primary cns lymphoma", "pcnsl",
                    "primary brain lymphoma"],
        "workflows": [NeuroWorkflowType.NEURO_ONCOLOGY],
        "search_terms": ["high-dose methotrexate", "rituximab",
                        "whole-brain radiation", "immunocompromised",
                        "EBV", "stereotactic biopsy", "diffusion restriction"],
    },

    # -- Headache & Migraine --------------------------------------------
    "migraine": {
        "aliases": ["migraine with aura", "migraine without aura",
                    "episodic migraine", "chronic migraine",
                    "hemiplegic migraine", "vestibular migraine",
                    "menstrual migraine", "retinal migraine"],
        "workflows": [NeuroWorkflowType.HEADACHE_DIAGNOSIS],
        "search_terms": ["CGRP", "triptan", "HIT-6", "MIDAS",
                        "erenumab", "galcanezumab", "fremanezumab",
                        "rimegepant", "gepant", "ditan",
                        "monthly migraine days", "preventive therapy"],
    },
    "cluster headache": {
        "aliases": ["cluster", "episodic cluster headache",
                    "chronic cluster headache", "trigeminal autonomic cephalalgia"],
        "workflows": [NeuroWorkflowType.HEADACHE_DIAGNOSIS],
        "search_terms": ["oxygen therapy", "sumatriptan injection",
                        "verapamil", "galcanezumab", "circadian",
                        "autonomic symptoms", "trigeminal autonomic cephalalgia"],
    },
    "tension-type headache": {
        "aliases": ["tension headache", "tth", "muscle contraction headache",
                    "chronic tension-type headache"],
        "workflows": [NeuroWorkflowType.HEADACHE_DIAGNOSIS],
        "search_terms": ["bilateral", "pressing quality", "pericranial tenderness",
                        "NSAID", "amitriptyline", "stress management"],
    },
    "medication overuse headache": {
        "aliases": ["moh", "analgesic rebound headache",
                    "medication adaptation headache"],
        "workflows": [NeuroWorkflowType.HEADACHE_DIAGNOSIS],
        "search_terms": ["withdrawal", "detoxification", "triptan overuse",
                        "opioid overuse", "bridging therapy", "preventive initiation"],
    },
    "idiopathic intracranial hypertension": {
        "aliases": ["iih", "pseudotumor cerebri", "benign intracranial hypertension"],
        "workflows": [NeuroWorkflowType.HEADACHE_DIAGNOSIS],
        "search_terms": ["papilledema", "opening pressure", "acetazolamide",
                        "optic nerve sheath fenestration", "CSF shunting",
                        "transverse sinus stenosis", "visual field loss"],
    },

    # -- Neuromuscular Disease ------------------------------------------
    "myasthenia gravis": {
        "aliases": ["mg", "myasthenia", "ocular myasthenia",
                    "generalized myasthenia gravis",
                    "acetylcholine receptor antibody positive mg",
                    "musk myasthenia"],
        "workflows": [NeuroWorkflowType.NEUROMUSCULAR_EVAL],
        "search_terms": ["AChR antibody", "MuSK antibody", "LRP4",
                        "pyridostigmine", "thymectomy", "IVIG", "PLEX",
                        "eculizumab", "efgartigimod", "FcRn inhibitor",
                        "myasthenic crisis", "decrement on RNS"],
    },
    "amyotrophic lateral sclerosis": {
        "aliases": ["als", "lou gehrig disease", "motor neuron disease",
                    "mnd", "upper motor neuron disease",
                    "lower motor neuron disease"],
        "workflows": [NeuroWorkflowType.NEUROMUSCULAR_EVAL],
        "search_terms": ["ALSFRS-R", "riluzole", "edaravone", "tofersen",
                        "SOD1", "C9orf72", "FUS", "TARDBP", "fasciculations",
                        "El Escorial criteria", "Awaji criteria",
                        "FVC", "respiratory failure", "NfL"],
    },
    "guillain-barre syndrome": {
        "aliases": ["gbs", "acute inflammatory demyelinating polyneuropathy",
                    "aidp", "miller fisher syndrome", "mfs"],
        "workflows": [NeuroWorkflowType.NEUROMUSCULAR_EVAL],
        "search_terms": ["ascending weakness", "areflexia", "albuminocytologic dissociation",
                        "IVIG", "plasmapheresis", "FVC monitoring",
                        "Campylobacter", "anti-ganglioside antibody",
                        "Hughes disability scale"],
    },
    "chronic inflammatory demyelinating polyneuropathy": {
        "aliases": ["cidp", "chronic inflammatory demyelinating polyradiculoneuropathy"],
        "workflows": [NeuroWorkflowType.NEUROMUSCULAR_EVAL],
        "search_terms": ["IVIG", "corticosteroids", "PLEX",
                        "conduction block", "temporal dispersion",
                        "EFNS/PNS criteria", "anti-NF155", "anti-CNTN1"],
    },
    "spinal muscular atrophy": {
        "aliases": ["sma", "sma type 1", "sma type 2", "sma type 3",
                    "werdnig-hoffmann disease", "kugelberg-welander disease"],
        "workflows": [NeuroWorkflowType.NEUROMUSCULAR_EVAL],
        "search_terms": ["SMN1", "SMN2", "nusinersen", "risdiplam",
                        "onasemnogene abeparvovec", "gene therapy",
                        "motor milestone", "newborn screening"],
    },
    "duchenne muscular dystrophy": {
        "aliases": ["dmd", "duchenne", "becker muscular dystrophy", "bmd"],
        "workflows": [NeuroWorkflowType.NEUROMUSCULAR_EVAL],
        "search_terms": ["dystrophin", "exon skipping", "gene therapy",
                        "corticosteroids", "6MWT", "CK elevation",
                        "Gowers sign", "cardiomyopathy"],
    },
    "peripheral neuropathy": {
        "aliases": ["diabetic neuropathy", "polyneuropathy",
                    "small fiber neuropathy", "large fiber neuropathy",
                    "charcot-marie-tooth", "cmt"],
        "workflows": [NeuroWorkflowType.NEUROMUSCULAR_EVAL],
        "search_terms": ["nerve conduction study", "EMG", "distal symmetric",
                        "neuropathic pain", "gabapentin", "pregabalin",
                        "duloxetine", "skin biopsy", "IENFD"],
    },

    # -- Autoimmune / Neuroimmunology (expanded) ----------------------------
    "autoimmune_encephalitis": {
        "aliases": ["autoimmune encephalitis", "anti-nmda receptor encephalitis",
                    "nmdar encephalitis", "lgi1 encephalitis", "caspr2 encephalitis",
                    "limbic encephalitis"],
        "workflows": [NeuroWorkflowType.EPILEPSY_CLASSIFICATION, NeuroWorkflowType.MS_MANAGEMENT],
        "search_terms": ["NMDA receptor antibody", "LGI1 antibody", "CASPR2 antibody",
                        "immunotherapy", "IVIG", "rituximab", "plasma exchange",
                        "faciobrachial dystonic seizures", "extreme delta brush"],
    },
    "neuromyelitis_optica": {
        "aliases": ["nmo", "neuromyelitis optica spectrum", "devic disease",
                    "aqp4-igg nmosd", "mog-igg nmosd"],
        "workflows": [NeuroWorkflowType.MS_MANAGEMENT],
        "search_terms": ["AQP4-IgG", "MOG-IgG", "longitudinally extensive transverse myelitis",
                        "optic neuritis", "area postrema syndrome", "eculizumab",
                        "inebilizumab", "satralizumab"],
    },
    "narcolepsy": {
        "aliases": ["narcolepsy type 1", "narcolepsy type 2", "narcolepsy with cataplexy",
                    "orexin deficiency", "hypocretin deficiency"],
        "workflows": [NeuroWorkflowType.GENERAL],
        "search_terms": ["orexin", "hypocretin", "cataplexy", "MSLT",
                        "sleep onset REM", "sodium oxybate", "pitolisant",
                        "modafinil", "CSF hypocretin", "HLA-DQB1*06:02"],
    },
    "restless_legs_syndrome": {
        "aliases": ["rls", "willis-ekbom disease", "restless legs",
                    "periodic limb movement disorder", "plmd"],
        "workflows": [NeuroWorkflowType.MOVEMENT_DISORDER],
        "search_terms": ["dopaminergic", "iron deficiency", "ferritin",
                        "augmentation", "gabapentin enacarbil", "pramipexole",
                        "ropinirole", "iron infusion", "periodic limb movements"],
    },
    "normal_pressure_hydrocephalus": {
        "aliases": ["nph", "normal pressure hydrocephalus", "idiopathic nph",
                    "hakim syndrome"],
        "workflows": [NeuroWorkflowType.DEMENTIA_EVALUATION, NeuroWorkflowType.MOVEMENT_DISORDER],
        "search_terms": ["gait apraxia", "cognitive impairment", "urinary incontinence",
                        "ventriculomegaly", "Evans index", "CSF tap test",
                        "ventriculoperitoneal shunt", "lumbar drain trial"],
    },
    "creutzfeldt_jakob_disease": {
        "aliases": ["cjd", "prion disease", "spongiform encephalopathy",
                    "variant cjd", "sporadic cjd", "familial cjd"],
        "workflows": [NeuroWorkflowType.DEMENTIA_EVALUATION],
        "search_terms": ["prion", "14-3-3 protein", "RT-QuIC", "rapidly progressive dementia",
                        "myoclonus", "cortical ribboning", "periodic sharp wave complexes",
                        "DWI restriction", "PRNP gene"],
    },
    "chiari_malformation": {
        "aliases": ["chiari", "chiari type 1", "chiari malformation type 1",
                    "arnold-chiari malformation", "tonsillar ectopia"],
        "workflows": [NeuroWorkflowType.HEADACHE_DIAGNOSIS],
        "search_terms": ["tonsillar herniation", "syringomyelia", "foramen magnum",
                        "posterior fossa decompression", "CSF flow study",
                        "Valsalva headache", "cine MRI"],
    },
    "idiopathic_intracranial_hypertension_expanded": {
        "aliases": ["iih expanded", "pseudotumor cerebri syndrome",
                    "raised intracranial pressure idiopathic"],
        "workflows": [NeuroWorkflowType.HEADACHE_DIAGNOSIS],
        "search_terms": ["papilledema", "elevated opening pressure", "acetazolamide",
                        "optic nerve sheath fenestration", "venous sinus stenting",
                        "visual field loss", "transverse sinus stenosis",
                        "empty sella", "obesity"],
    },
    "bells_palsy": {
        "aliases": ["bell's palsy", "bell palsy", "idiopathic facial palsy",
                    "facial nerve palsy", "lower motor neuron facial weakness"],
        "workflows": [NeuroWorkflowType.NEUROMUSCULAR_EVAL],
        "search_terms": ["facial nerve", "LMN pattern", "House-Brackmann scale",
                        "prednisone", "valacyclovir", "eye protection",
                        "synkinesis", "electrodiagnostics"],
    },
    "trigeminal_neuralgia": {
        "aliases": ["tic douloureux", "trigeminal neuralgia type 1",
                    "trigeminal neuralgia type 2", "tn"],
        "workflows": [NeuroWorkflowType.HEADACHE_DIAGNOSIS],
        "search_terms": ["carbamazepine", "oxcarbazepine", "microvascular decompression",
                        "gamma knife", "trigeminal nerve", "neurovascular conflict",
                        "trigger zone", "lancinating pain"],
    },
    "cerebral_venous_thrombosis_expanded": {
        "aliases": ["cvt expanded", "dural sinus thrombosis syndrome",
                    "cortical vein thrombosis"],
        "workflows": [NeuroWorkflowType.STROKE_ACUTE],
        "search_terms": ["anticoagulation", "heparin", "MR venography",
                        "empty delta sign", "venous infarct", "hemorrhagic transformation",
                        "oral contraceptives", "prothrombotic states"],
    },
    "cavernous_malformation": {
        "aliases": ["cavernoma", "cavernous angioma", "cerebral cavernous malformation",
                    "ccm"],
        "workflows": [NeuroWorkflowType.NEURO_ONCOLOGY, NeuroWorkflowType.EPILEPSY_CLASSIFICATION],
        "search_terms": ["seizures", "hemorrhage", "popcorn appearance", "hemosiderin ring",
                        "surgical resection", "CCM1", "CCM2", "CCM3",
                        "familial cavernous malformation"],
    },
    "spinal_cord_injury": {
        "aliases": ["sci", "spinal cord injury", "traumatic myelopathy",
                    "spinal cord compression", "acute spinal cord injury"],
        "workflows": [NeuroWorkflowType.NEUROMUSCULAR_EVAL],
        "search_terms": ["ASIA classification", "ASIA impairment scale",
                        "neurogenic shock", "spinal shock", "methylprednisolone",
                        "decompressive surgery", "rehabilitation", "autonomic dysreflexia"],
    },
    "traumatic_brain_injury": {
        "aliases": ["tbi", "concussion", "mild tbi", "moderate tbi", "severe tbi",
                    "diffuse axonal injury", "dai"],
        "workflows": [NeuroWorkflowType.GENERAL],
        "search_terms": ["GCS", "Glasgow Coma Scale", "intracranial pressure monitoring",
                        "cerebral perfusion pressure", "decompressive craniectomy",
                        "post-traumatic epilepsy", "chronic traumatic encephalopathy",
                        "post-concussion syndrome"],
    },
    "status_epilepticus_expanded": {
        "aliases": ["refractory se", "super-refractory status epilepticus",
                    "established status epilepticus"],
        "workflows": [NeuroWorkflowType.EPILEPSY_CLASSIFICATION],
        "search_terms": ["benzodiazepine protocol", "lorazepam", "midazolam",
                        "fosphenytoin", "levetiracetam IV", "valproate IV",
                        "continuous EEG monitoring", "burst suppression target",
                        "pentobarbital coma", "treatment algorithm",
                        "convulsive SE", "non-convulsive SE"],
    },
}


NEURO_DRUGS: Dict[str, Dict[str, object]] = {

    # -- Stroke / Cerebrovascular ----------------------------------------
    "alteplase": {
        "aliases": ["tpa", "activase", "rt-pa", "tissue plasminogen activator"],
        "mechanism": "Tissue plasminogen activator (thrombolytic)",
        "indications": ["acute ischemic stroke within 4.5 hours"],
        "workflows": ["stroke_acute"],
    },
    "tenecteplase": {
        "aliases": ["tnkase"],
        "mechanism": "Modified tissue plasminogen activator (thrombolytic)",
        "indications": ["acute ischemic stroke", "large vessel occlusion"],
        "workflows": ["stroke_acute"],
    },
    "clopidogrel": {
        "aliases": ["plavix"],
        "mechanism": "P2Y12 receptor antagonist (antiplatelet)",
        "indications": ["secondary stroke prevention", "TIA"],
        "workflows": ["stroke_prevention"],
    },
    "apixaban": {
        "aliases": ["eliquis"],
        "mechanism": "Direct factor Xa inhibitor (DOAC)",
        "indications": ["stroke prevention in atrial fibrillation"],
        "workflows": ["stroke_prevention"],
    },

    # -- Dementia / Neurodegeneration ------------------------------------
    "lecanemab": {
        "aliases": ["leqembi"],
        "mechanism": "Anti-amyloid beta monoclonal antibody (protofibrils)",
        "indications": ["early Alzheimer's disease with confirmed amyloid"],
        "workflows": ["dementia_evaluation"],
    },
    "donanemab": {
        "aliases": ["kisunla"],
        "mechanism": "Anti-amyloid beta monoclonal antibody (N3pG-modified Abeta)",
        "indications": ["early Alzheimer's disease with confirmed amyloid"],
        "workflows": ["dementia_evaluation"],
    },
    "donepezil": {
        "aliases": ["aricept"],
        "mechanism": "Cholinesterase inhibitor",
        "indications": ["Alzheimer's disease", "dementia with Lewy bodies"],
        "workflows": ["dementia_evaluation"],
    },
    "memantine": {
        "aliases": ["namenda"],
        "mechanism": "NMDA receptor antagonist",
        "indications": ["moderate-to-severe Alzheimer's disease"],
        "workflows": ["dementia_evaluation"],
    },
    "levodopa-carbidopa": {
        "aliases": ["sinemet", "levodopa", "l-dopa", "carbidopa-levodopa",
                    "rytary", "duopa"],
        "mechanism": "Dopamine precursor + decarboxylase inhibitor",
        "indications": ["Parkinson disease"],
        "workflows": ["movement_disorder"],
    },
    "pramipexole": {
        "aliases": ["mirapex"],
        "mechanism": "Dopamine D2/D3 receptor agonist",
        "indications": ["Parkinson disease", "restless legs syndrome"],
        "workflows": ["movement_disorder"],
    },
    "ropinirole": {
        "aliases": ["requip"],
        "mechanism": "Dopamine D2/D3 receptor agonist",
        "indications": ["Parkinson disease", "restless legs syndrome"],
        "workflows": ["movement_disorder"],
    },

    # -- Epilepsy --------------------------------------------------------
    "levetiracetam": {
        "aliases": ["keppra"],
        "mechanism": "SV2A modulator",
        "indications": ["focal epilepsy", "generalized epilepsy"],
        "workflows": ["epilepsy_classification"],
    },
    "lamotrigine": {
        "aliases": ["lamictal"],
        "mechanism": "Sodium channel blocker, glutamate release inhibitor",
        "indications": ["focal epilepsy", "generalized epilepsy", "bipolar disorder"],
        "workflows": ["epilepsy_classification"],
    },
    "valproate": {
        "aliases": ["valproic acid", "depakote", "depakene", "divalproex"],
        "mechanism": "Multiple (GABA, sodium channels, T-type calcium channels)",
        "indications": ["generalized epilepsy", "focal epilepsy", "status epilepticus"],
        "workflows": ["epilepsy_classification"],
    },
    "carbamazepine": {
        "aliases": ["tegretol"],
        "mechanism": "Sodium channel blocker",
        "indications": ["focal epilepsy", "trigeminal neuralgia"],
        "workflows": ["epilepsy_classification"],
    },
    "cenobamate": {
        "aliases": ["xcopri"],
        "mechanism": "Sodium channel blocker + GABA-A modulator",
        "indications": ["focal epilepsy (drug-resistant)"],
        "workflows": ["epilepsy_classification"],
    },
    "fenfluramine": {
        "aliases": ["fintepla"],
        "mechanism": "Serotonin releasing agent (sigma-1 agonist)",
        "indications": ["Dravet syndrome", "Lennox-Gastaut syndrome"],
        "workflows": ["epilepsy_classification"],
    },
    "cannabidiol": {
        "aliases": ["epidiolex", "cbd"],
        "mechanism": "Multiple (GPR55 antagonist, TRPV1 agonist, adenosine reuptake inhibitor)",
        "indications": ["Dravet syndrome", "Lennox-Gastaut syndrome", "tuberous sclerosis"],
        "workflows": ["epilepsy_classification"],
    },

    # -- Multiple Sclerosis ----------------------------------------------
    "ocrelizumab": {
        "aliases": ["ocrevus"],
        "mechanism": "Anti-CD20 monoclonal antibody",
        "indications": ["RRMS", "PPMS"],
        "workflows": ["ms_management"],
    },
    "natalizumab": {
        "aliases": ["tysabri"],
        "mechanism": "Anti-alpha-4 integrin monoclonal antibody",
        "indications": ["RRMS (highly active)"],
        "workflows": ["ms_management"],
    },
    "ofatumumab": {
        "aliases": ["kesimpta"],
        "mechanism": "Anti-CD20 monoclonal antibody (subcutaneous)",
        "indications": ["RRMS", "active SPMS"],
        "workflows": ["ms_management"],
    },
    "dimethyl fumarate": {
        "aliases": ["tecfidera", "dmf"],
        "mechanism": "Nrf2 activator, NF-kB modulator",
        "indications": ["RRMS"],
        "workflows": ["ms_management"],
    },
    "fingolimod": {
        "aliases": ["gilenya"],
        "mechanism": "Sphingosine-1-phosphate receptor modulator",
        "indications": ["RRMS"],
        "workflows": ["ms_management"],
    },

    # -- Headache / Migraine ---------------------------------------------
    "erenumab": {
        "aliases": ["aimovig"],
        "mechanism": "CGRP receptor monoclonal antibody",
        "indications": ["episodic migraine prevention", "chronic migraine prevention"],
        "workflows": ["headache_diagnosis"],
    },
    "galcanezumab": {
        "aliases": ["emgality"],
        "mechanism": "Anti-CGRP monoclonal antibody",
        "indications": ["migraine prevention", "episodic cluster headache"],
        "workflows": ["headache_diagnosis"],
    },
    "fremanezumab": {
        "aliases": ["ajovy"],
        "mechanism": "Anti-CGRP monoclonal antibody",
        "indications": ["episodic migraine prevention", "chronic migraine prevention"],
        "workflows": ["headache_diagnosis"],
    },
    "rimegepant": {
        "aliases": ["nurtec"],
        "mechanism": "CGRP receptor antagonist (gepant, oral)",
        "indications": ["acute migraine treatment", "migraine prevention"],
        "workflows": ["headache_diagnosis"],
    },

    # -- Neuromuscular ---------------------------------------------------
    "riluzole": {
        "aliases": ["rilutek"],
        "mechanism": "Glutamate release inhibitor / sodium channel blocker",
        "indications": ["amyotrophic lateral sclerosis"],
        "workflows": ["neuromuscular_eval"],
    },
    "tofersen": {
        "aliases": ["qalsody"],
        "mechanism": "Antisense oligonucleotide targeting SOD1 mRNA",
        "indications": ["SOD1-ALS"],
        "workflows": ["neuromuscular_eval"],
    },
    "eculizumab": {
        "aliases": ["soliris"],
        "mechanism": "Anti-complement C5 monoclonal antibody",
        "indications": ["generalized myasthenia gravis (AChR+)", "NMOSD (AQP4+)"],
        "workflows": ["neuromuscular_eval", "ms_management"],
    },
    "efgartigimod": {
        "aliases": ["vyvgart"],
        "mechanism": "Neonatal Fc receptor (FcRn) blocker",
        "indications": ["generalized myasthenia gravis (AChR+)", "CIDP"],
        "workflows": ["neuromuscular_eval"],
    },
    "nusinersen": {
        "aliases": ["spinraza"],
        "mechanism": "Antisense oligonucleotide modulating SMN2 splicing",
        "indications": ["spinal muscular atrophy"],
        "workflows": ["neuromuscular_eval"],
    },
    "risdiplam": {
        "aliases": ["evrysdi"],
        "mechanism": "SMN2 splicing modifier (oral small molecule)",
        "indications": ["spinal muscular atrophy"],
        "workflows": ["neuromuscular_eval"],
    },

    # -- Expanded drugs (10 new entries) ------------------------------------
    "rimegepant_expanded": {
        "aliases": ["nurtec odt", "rimegepant"],
        "mechanism": "CGRP receptor antagonist (gepant, oral dissolving tablet)",
        "indications": ["acute migraine treatment", "episodic migraine prevention"],
        "workflows": ["headache_diagnosis"],
    },
    "ubrogepant": {
        "aliases": ["ubrelvy"],
        "mechanism": "CGRP receptor antagonist (gepant, oral)",
        "indications": ["acute migraine treatment"],
        "workflows": ["headache_diagnosis"],
    },
    "tolebrutinib": {
        "aliases": [],
        "mechanism": "Bruton's tyrosine kinase (BTK) inhibitor (CNS-penetrant)",
        "indications": ["relapsing MS (Phase III)", "progressive MS (Phase III)"],
        "workflows": ["ms_management"],
    },
    "ublituximab": {
        "aliases": ["briumvi"],
        "mechanism": "Anti-CD20 monoclonal antibody (glycoengineered)",
        "indications": ["relapsing forms of MS"],
        "workflows": ["ms_management"],
    },
    "ganaxolone": {
        "aliases": ["ztalmy"],
        "mechanism": "Positive allosteric modulator of GABA-A receptors (neurosteroid)",
        "indications": ["CDKL5 deficiency disorder seizures"],
        "workflows": ["epilepsy_classification"],
    },
    "valiltramiprosate": {
        "aliases": ["alz-801"],
        "mechanism": "Anti-amyloid aggregation inhibitor (oral, prodrug of tramiprosate)",
        "indications": ["Alzheimer's disease in APOE4/4 homozygous carriers"],
        "workflows": ["dementia_evaluation"],
    },
    "zuranolone": {
        "aliases": ["zurzuvae"],
        "mechanism": "Neuroactive steroid, positive allosteric modulator of GABA-A receptors",
        "indications": ["postpartum depression", "major depressive disorder"],
        "workflows": ["general"],
    },
    "atogepant": {
        "aliases": ["qulipta"],
        "mechanism": "CGRP receptor antagonist (gepant, oral daily)",
        "indications": ["episodic migraine prevention", "chronic migraine prevention"],
        "workflows": ["headache_diagnosis"],
    },
    "safinamide": {
        "aliases": ["xadago"],
        "mechanism": "Selective reversible MAO-B inhibitor + sodium channel modulator",
        "indications": ["Parkinson disease adjunct to levodopa for motor fluctuations"],
        "workflows": ["movement_disorder"],
    },
    "tofersen_expanded": {
        "aliases": ["qalsody tofersen"],
        "mechanism": "Antisense oligonucleotide targeting SOD1 mRNA",
        "indications": ["SOD1-ALS (superoxide dismutase 1 amyotrophic lateral sclerosis)"],
        "workflows": ["neuromuscular_eval"],
    },
}


NEURO_BIOMARKERS: Dict[str, Dict[str, str]] = {
    "amyloid-beta-42": {
        "full_name": "Amyloid Beta 1-42 (A\u03b242)",
        "assay": "CSF immunoassay (Lumipulse, Elecsys), plasma (Simoa, IP-MS)",
        "significance": "Core AD biomarker (A in ATN framework); decreased CSF A\u03b242 or "
                        "elevated plasma A\u03b242/40 ratio indicates amyloid pathology",
        "workflows": "dementia_evaluation",
    },
    "p-tau181": {
        "full_name": "Phosphorylated Tau 181",
        "assay": "CSF immunoassay, plasma (Simoa, Lumipulse)",
        "significance": "Core AD biomarker (T in ATN); elevated in Alzheimer's disease, "
                        "correlates with tau tangle pathology and disease stage",
        "workflows": "dementia_evaluation",
    },
    "p-tau217": {
        "full_name": "Phosphorylated Tau 217",
        "assay": "Plasma immunoassay (Lumipulse, ALZpath pTau217)",
        "significance": "Highest-performing blood-based AD biomarker; "
                        "discriminates AD from non-AD with >95% accuracy in validation studies",
        "workflows": "dementia_evaluation",
    },
    "nfl": {
        "full_name": "Neurofilament Light Chain (NfL)",
        "assay": "Serum/CSF Simoa, ELISA",
        "significance": "Non-specific neuronal injury marker; elevated in MS, ALS, "
                        "FTD, stroke, traumatic brain injury; tracks disease activity in MS "
                        "and prognosis in ALS; FDA cleared as MS biomarker",
        "workflows": "dementia_evaluation,ms_management,neuromuscular_eval",
    },
    "gfap": {
        "full_name": "Glial Fibrillary Acidic Protein (GFAP)",
        "assay": "Serum/CSF Simoa, ELISA",
        "significance": "Astrocytic injury marker; elevated in AD (correlates with amyloid), "
                        "traumatic brain injury, NMOSD attacks, glioblastoma; "
                        "emerging plasma biomarker for AD screening",
        "workflows": "dementia_evaluation,ms_management",
    },
    "csf-ocb": {
        "full_name": "CSF Oligoclonal Bands",
        "assay": "Isoelectric focusing with immunofixation",
        "significance": "Intrathecal IgG synthesis; present in >95% of MS patients; "
                        "McDonald 2017 criteria allows substitution for DIT; "
                        "also seen in NMOSD, CNS infections, autoimmune encephalitis",
        "workflows": "ms_management",
    },
    "aqp4-igg": {
        "full_name": "Aquaporin-4 IgG Antibody",
        "assay": "Cell-based assay (CBA), ELISA, flow cytometry",
        "significance": "Highly specific biomarker for NMOSD; seropositive in ~75% of NMOSD; "
                        "distinguishes NMOSD from MS; required for many DMT approvals",
        "workflows": "ms_management",
    },
    "mog-igg": {
        "full_name": "Myelin Oligodendrocyte Glycoprotein IgG Antibody",
        "assay": "Live cell-based assay (CBA)",
        "significance": "Defines MOG antibody disease (MOGAD); distinguishes from MS and NMOSD; "
                        "seropositivity guides treatment decisions; may fluctuate",
        "workflows": "ms_management",
    },
    "achr-antibody": {
        "full_name": "Acetylcholine Receptor Antibody",
        "assay": "Radioimmunoassay (RIA), cell-based assay",
        "significance": "Diagnostic biomarker for myasthenia gravis; seropositive in ~85% "
                        "generalized MG and ~50% ocular MG; modulating, binding, blocking subtypes",
        "workflows": "neuromuscular_eval",
    },
    "musk-antibody": {
        "full_name": "Muscle-Specific Kinase Antibody",
        "assay": "Cell-based assay (CBA), RIA",
        "significance": "Diagnostic biomarker for MuSK-MG (~5-8% of MG); "
                        "bulbar predominance, poor response to cholinesterase inhibitors; "
                        "guides treatment (rituximab preferred)",
        "workflows": "neuromuscular_eval",
    },
    "anti-ganglioside": {
        "full_name": "Anti-Ganglioside Antibodies (GM1, GQ1b, GD1a)",
        "assay": "ELISA, line immunoassay",
        "significance": "GBS subtype classification: anti-GM1 in AMAN, "
                        "anti-GQ1b in Miller Fisher syndrome; "
                        "helps predict clinical course and outcome",
        "workflows": "neuromuscular_eval",
    },
    "csf-protein": {
        "full_name": "CSF Total Protein / Albuminocytologic Dissociation",
        "assay": "Lumbar puncture with protein measurement",
        "significance": "Elevated protein with normal cell count (albuminocytologic dissociation) "
                        "in GBS and CIDP; also elevated in spinal cord tumors and CNS infections",
        "workflows": "neuromuscular_eval",
    },
    "dat-spect": {
        "full_name": "Dopamine Transporter SPECT (DaT-SPECT / DaTscan)",
        "assay": "I-123 ioflupane SPECT imaging",
        "significance": "Reduced striatal uptake distinguishes neurodegenerative parkinsonism "
                        "(PD, MSA, PSP, DLB) from essential tremor, drug-induced parkinsonism, "
                        "and psychogenic tremor; supports DLB diagnosis",
        "workflows": "movement_disorder,dementia_evaluation",
    },
    "amyloid-pet": {
        "full_name": "Amyloid PET Imaging",
        "assay": "PET with florbetapir (Amyvid), florbetaben (Neuraceq), "
                "or flutemetamol (Vizamyl)",
        "significance": "Enrollment and outcome biomarker in AD trials; confirms amyloid "
                        "pathology for anti-amyloid therapy eligibility (lecanemab, donanemab); "
                        "negative scan essentially excludes AD",
        "workflows": "dementia_evaluation",
    },
    "tau-pet": {
        "full_name": "Tau PET Imaging",
        "assay": "PET with flortaucipir (Tauvid) or MK-6240",
        "significance": "Visualizes tau tangle distribution; correlates with cognitive decline "
                        "better than amyloid PET; Braak staging in vivo; helps distinguish "
                        "AD from non-AD tauopathies; emerging trial endpoint",
        "workflows": "dementia_evaluation",
    },
    "eeg-findings": {
        "full_name": "Electroencephalography Findings",
        "assay": "Scalp EEG (routine, prolonged, continuous, video-EEG)",
        "significance": "Seizure classification (focal vs generalized), epileptiform discharge "
                        "localization, status epilepticus diagnosis, encephalopathy grading, "
                        "periodic patterns (NCSE, CJD, SSPE), sleep staging",
        "workflows": "epilepsy_classification",
    },

    # -- Expanded biomarkers (5 new entries) --------------------------------
    "14-3-3-protein": {
        "full_name": "14-3-3 Protein (CSF)",
        "assay": "CSF Western blot, ELISA",
        "significance": "Marker of rapid neuronal destruction; elevated in sporadic CJD "
                        "(sensitivity ~90%); also elevated in herpes encephalitis, stroke, "
                        "and other conditions causing rapid neuronal death; best used with "
                        "RT-QuIC for CJD diagnosis",
        "workflows": "dementia_evaluation",
    },
    "rt-quic": {
        "full_name": "Real-Time Quaking-Induced Conversion (RT-QuIC)",
        "assay": "CSF prion seed amplification assay",
        "significance": "Gold-standard antemortem prion detection; sensitivity >92% and "
                        "specificity ~100% for sporadic CJD; detects misfolded prion protein "
                        "seed amplification; has largely replaced brain biopsy for CJD diagnosis",
        "workflows": "dementia_evaluation",
    },
    "flair-hyperintensity": {
        "full_name": "FLAIR Hyperintensity Pattern (MRI)",
        "assay": "MRI FLAIR sequence (fluid-attenuated inversion recovery)",
        "significance": "Non-specific marker seen across multiple neurological conditions: "
                        "MS periventricular lesions (Dawson fingers), autoimmune encephalitis "
                        "(mesial temporal), status epilepticus (cortical), PRES (posterior), "
                        "glioma infiltration; pattern and distribution guide differential diagnosis",
        "workflows": "ms_management,epilepsy_classification,dementia_evaluation",
    },
    "dat-binding-ratio": {
        "full_name": "DaT Binding Ratio (DaT-SPECT Quantification)",
        "assay": "Semi-quantitative analysis of 123I-ioflupane SPECT",
        "significance": "Quantifies striatal dopamine transporter density; specific binding ratio "
                        "(SBR) below age-adjusted normal indicates presynaptic dopaminergic "
                        "deficit; putamen-to-caudate ratio helps distinguish PD from atypical "
                        "parkinsonism; serial measurements track disease progression",
        "workflows": "movement_disorder",
    },
    "orexin-level": {
        "full_name": "CSF Orexin (Hypocretin-1) Level",
        "assay": "CSF immunoassay (radioimmunoassay or ELISA)",
        "significance": "Diagnostic biomarker for narcolepsy type 1; CSF orexin-A level "
                        "<110 pg/mL is diagnostic with >87% sensitivity and ~99% specificity; "
                        "reflects hypothalamic orexin neuron loss; normal levels in narcolepsy "
                        "type 2 and other hypersomnias",
        "workflows": "general",
    },
}


# =====================================================================
# SEARCH PLAN DATACLASS
# =====================================================================

@dataclass
class SearchPlan:
    """Agent's plan for answering a neurology intelligence question.

    The search plan captures all entities detected in the user's question
    and the strategy the agent will use to retrieve evidence from the
    14 neurology-specific Milvus collections.
    """
    question: str
    conditions: List[str] = field(default_factory=list)
    drugs: List[str] = field(default_factory=list)
    biomarkers: List[str] = field(default_factory=list)
    relevant_workflows: List[NeuroWorkflowType] = field(default_factory=list)
    search_strategy: str = "broad"  # broad, targeted, differential, emergency
    sub_questions: List[str] = field(default_factory=list)
    identified_topics: List[str] = field(default_factory=list)


# =====================================================================
# NEUROLOGY INTELLIGENCE AGENT
# =====================================================================

class NeurologyAgent:
    """Autonomous Neurology Intelligence Agent.

    Wraps the multi-collection NeuroRAGEngine with planning and reasoning
    capabilities. Designed to answer complex cross-functional questions
    about neurological diagnosis, treatment, and management.

    Example queries this agent handles:
    - "68-year-old with acute left hemiparesis, NIHSS 14, last known well 2 hours ago"
    - "Classify this EEG pattern: 3 Hz generalized spike-wave in a 7-year-old"
    - "Differentiate Parkinson disease from MSA given autonomic failure and poor levodopa response"
    - "MRI shows new T2/FLAIR lesions with gadolinium enhancement -- assess for MS relapse"
    - "CSF shows elevated protein, normal glucose, 2 WBCs -- differential for polyneuropathy"
    - "APOE e4/e4, MoCA 22/30, hippocampal atrophy on MRI -- stage dementia per NIA-AA ATN"
    - "Design a treatment plan for drug-resistant focal epilepsy after two ASM failures"
    - "Evaluate headache red flags: thunderclap onset, papilledema, CN VI palsy"

    Usage:
        agent = NeurologyAgent(rag_engine)
        plan = agent.search_plan("68yo acute stroke, NIHSS 14, 2 hours from onset")
        response = agent.run("68yo acute stroke, NIHSS 14, 2 hours from onset")
    """

    def __init__(self, rag_engine):
        """Initialize agent with a configured RAG engine.

        Args:
            rag_engine: NeuroRAGEngine instance with Milvus collections connected.
        """
        self.rag = rag_engine
        self.knowledge = {
            "conditions": NEURO_CONDITIONS,
            "biomarkers": NEURO_BIOMARKERS,
            "drugs": NEURO_DRUGS,
        }

    # -- Public API --------------------------------------------------------

    def run(
        self,
        query: str,
        workflow_type: Optional[NeuroWorkflowType] = None,
        patient_context: Optional[dict] = None,
        **kwargs,
    ) -> NeuroResponse:
        """Execute the full agent pipeline: plan -> search -> evaluate -> synthesize.

        Args:
            query: Natural language question about neurology.
            workflow_type: Optional workflow override for collection boosting.
            patient_context: Optional patient data for clinical context.
            **kwargs: Additional query parameters (top_k, collection_filter).

        Returns:
            NeuroResponse with findings, recommendations, and metadata.
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

        # Phase 5: Check for clinical alerts
        if hasattr(response, "clinical_alerts"):
            response.clinical_alerts = self._check_clinical_alerts(query, plan)

        return response

    def search_plan(self, question: str) -> SearchPlan:
        """Analyze a question and create an optimised search plan.

        Detects neurological conditions, drugs, and biomarkers in the question
        text. Determines relevant neuro workflows, chooses a search
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
            info = NEURO_CONDITIONS.get(condition, {})
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

    def generate_report(
        self,
        results,
        workflow: Optional[NeuroWorkflowType] = None,
    ) -> str:
        """Generate a structured neurology intelligence report.

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
            "# Neurology Intelligence Report",
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
                    f"- **Conditions:** {', '.join(plan.conditions)}"
                )
            if plan.drugs:
                report_lines.append(
                    f"- **Medications:** {', '.join(plan.drugs)}"
                )
            if plan.biomarkers:
                report_lines.append(
                    f"- **Biomarkers/Imaging:** {', '.join(plan.biomarkers)}"
                )
            report_lines.append("")

        # Clinical alerts check
        alerts = self._check_clinical_alerts(question, plan)
        if alerts:
            report_lines.extend([
                "---",
                "",
                "## [CRITICAL] Clinical Alerts",
                "",
            ])
            for alert in alerts:
                report_lines.append(f"- **[CRITICAL]** {alert}")
            report_lines.append("")

        # Critical findings from results
        critical_flags = []
        if hasattr(results, "results") and results.results:
            for r in results.results:
                meta = r.metadata if hasattr(r, "metadata") else {}
                if meta.get("urgency") == "critical" or meta.get("safety_alert"):
                    critical_flags.append(r)

        if critical_flags:
            if not alerts:
                report_lines.extend([
                    "---",
                    "",
                    "## [CRITICAL] Safety Alerts",
                    "",
                ])
            for flag in critical_flags:
                text = flag.text if hasattr(flag, "text") else str(flag)
                report_lines.append(
                    f"- **[CRITICAL]** {text[:200]} -- "
                    f"immediate neurological review required."
                )
            report_lines.append("")

        # Analysis section
        report_lines.extend([
            "---",
            "",
            "## Clinical Analysis",
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

        # Guideline references section
        if workflow in (NeuroWorkflowType.STROKE_ACUTE,
                        NeuroWorkflowType.STROKE_PREVENTION):
            report_lines.extend([
                "---",
                "",
                "## Guideline References",
                "",
                "- AHA/ASA 2019: Guidelines for Early Management of Acute Ischemic Stroke",
                "- AHA/ASA 2022: Guideline for Management of Spontaneous Intracerebral Hemorrhage",
                "- AHA/ASA 2021: Secondary Prevention of Stroke and TIA",
                "- ESO 2022: European Stroke Organisation Guidelines on Intravenous Thrombolysis",
                "",
            ])
        elif workflow == NeuroWorkflowType.EPILEPSY_CLASSIFICATION:
            report_lines.extend([
                "---",
                "",
                "## Guideline References",
                "",
                "- ILAE 2017: Operational Classification of Seizure Types",
                "- ILAE 2017: Classification of the Epilepsies",
                "- AAN/AES 2022: Treatment of New-Onset Epilepsy Practice Guideline",
                "- AES/NAEC 2016: Drug-Resistant Epilepsy Evaluation",
                "",
            ])
        elif workflow == NeuroWorkflowType.DEMENTIA_EVALUATION:
            report_lines.extend([
                "---",
                "",
                "## Guideline References",
                "",
                "- NIA-AA 2018: Research Framework for AD (ATN Classification)",
                "- NIA-AA 2024: Revised Clinical Criteria for Alzheimer's Disease",
                "- AAN 2018: Practice Parameter: Mild Cognitive Impairment",
                "- DLB Consortium 2017: Diagnosis and Management of DLB (Fourth Report)",
                "",
            ])
        elif workflow == NeuroWorkflowType.MS_MANAGEMENT:
            report_lines.extend([
                "---",
                "",
                "## Guideline References",
                "",
                "- AAN 2018: Practice Guideline for Disease-Modifying Therapies in MS",
                "- EAN/ECTRIMS 2024: Treatment Guidelines for RRMS and Progressive MS",
                "- McDonald 2017: Diagnostic Criteria for Multiple Sclerosis",
                "- IPND 2015: International Consensus Diagnostic Criteria for NMOSD",
                "",
            ])
        elif workflow == NeuroWorkflowType.MOVEMENT_DISORDER:
            report_lines.extend([
                "---",
                "",
                "## Guideline References",
                "",
                "- MDS 2015: Clinical Diagnostic Criteria for Parkinson's Disease",
                "- AAN 2019: Practice Parameter: Treatment of PD with Motor Fluctuations",
                "- MDS 2017: Criteria for Prodromal Parkinson's Disease",
                "- AAN 2011: Practice Parameter: Essential Tremor Update",
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
            "*This report is generated by the Neurology Intelligence Agent "
            "within the HCLS AI Factory. All recommendations require review by "
            "board-certified neurologists with access to the full clinical picture, "
            "examination findings, and imaging studies.*",
        ])

        return "\n".join(report_lines)

    # -- Clinical Alert Detection ----------------------------------------

    def _check_clinical_alerts(
        self,
        question: str,
        plan: SearchPlan,
    ) -> List[str]:
        """Check for clinical situations requiring urgent alerts.

        Scans the query text and detected entities for patterns that
        indicate time-critical neurological emergencies.

        Args:
            question: Original query text.
            plan: SearchPlan with detected entities.

        Returns:
            List of alert strings for critical findings.
        """
        alerts: List[str] = []
        text_upper = question.upper()

        # Acute stroke within treatment window
        stroke_keywords = ["ACUTE STROKE", "ISCHEMIC STROKE", "LVO",
                           "LARGE VESSEL OCCLUSION", "HEMIPARESIS",
                           "HEMIPLEGIA", "APHASIA"]
        time_keywords = ["HOUR", "MINUTES", "ONSET", "LAST KNOWN WELL",
                         "LKW", "WITHIN"]
        if any(kw in text_upper for kw in stroke_keywords):
            if any(kw in text_upper for kw in time_keywords):
                alerts.append(
                    "Possible acute stroke within treatment window detected. "
                    "Evaluate for IV thrombolysis (< 4.5h) and/or mechanical "
                    "thrombectomy (< 24h with favorable perfusion imaging). "
                    "Activate stroke code immediately."
                )

        # Status epilepticus
        if "STATUS EPILEPTICUS" in text_upper or (
            "CONTINUOUS SEIZURE" in text_upper or "PROLONGED SEIZURE" in text_upper
        ):
            alerts.append(
                "Status epilepticus suspected. Initiate benzodiazepine protocol "
                "immediately (lorazepam 0.1 mg/kg IV or midazolam 0.2 mg/kg IM). "
                "If refractory after second-line ASM, consider continuous infusion "
                "and continuous EEG monitoring."
            )

        # Intracranial hemorrhage
        if any(kw in text_upper for kw in [
            "INTRACEREBRAL HEMORRHAGE", "ICH", "BRAIN BLEED",
            "HEMATOMA EXPANSION", "INTRAPARENCHYMAL",
        ]):
            if any(kw in text_upper for kw in [
                "ANTICOAGULANT", "WARFARIN", "DOAC", "COUMADIN",
            ]):
                alerts.append(
                    "Anticoagulant-associated intracerebral hemorrhage. "
                    "Immediate reversal required: idarucizumab for dabigatran, "
                    "andexanet alfa or 4F-PCC for factor Xa inhibitors, "
                    "4F-PCC + vitamin K for warfarin. Blood pressure target < 140 mmHg."
                )

        # Herniation signs
        if any(kw in text_upper for kw in [
            "HERNIATION", "UNCAL HERNIATION", "TONSILLAR HERNIATION",
            "CUSHING TRIAD", "BLOWN PUPIL", "FIXED DILATED PUPIL",
            "POSTURING", "DECEREBRATE", "DECORTICATE",
        ]):
            alerts.append(
                "Signs of cerebral herniation detected. Immediate intervention "
                "required: head of bed elevation 30 degrees, hyperventilation "
                "to pCO2 30-35 mmHg, osmotic therapy (mannitol 1 g/kg or "
                "hypertonic saline 23.4%), and emergent neurosurgical consultation."
            )

        # Respiratory failure in neuromuscular disease
        if any(kw in text_upper for kw in [
            "MYASTHENIC CRISIS", "GBS", "GUILLAIN-BARRE",
        ]):
            if any(kw in text_upper for kw in [
                "FVC", "RESPIRATORY", "INTUBAT", "VENTILAT",
                "DYSPNEA", "ORTHOPNEA",
            ]):
                alerts.append(
                    "Neuromuscular respiratory failure risk. Monitor FVC and NIF "
                    "serially (q4-6h). Intubation criteria: FVC < 15-20 mL/kg, "
                    "NIF > -20 to -30 cmH2O, or clinical deterioration. "
                    "Initiate IVIG or plasmapheresis urgently."
                )

        # Thunderclap headache (SAH concern)
        if any(kw in text_upper for kw in [
            "THUNDERCLAP", "WORST HEADACHE", "SUDDEN ONSET HEADACHE",
        ]):
            alerts.append(
                "Thunderclap headache raises concern for subarachnoid hemorrhage. "
                "Emergent non-contrast CT head required. If CT negative, lumbar "
                "puncture for xanthochromia is indicated (after 6-12 hours from "
                "onset). Consider CTA for aneurysm evaluation."
            )

        return alerts

    # -- Workflow Detection -----------------------------------------------

    def _detect_workflow(self, question: str) -> NeuroWorkflowType:
        """Detect the most relevant workflow from a question.

        Uses keyword-based heuristics to identify which of the 9 neurology
        workflows is most relevant to the query.

        Args:
            question: The user's natural language question.

        Returns:
            Most relevant NeuroWorkflowType.
        """
        text_upper = question.upper()

        workflow_scores: Dict[NeuroWorkflowType, float] = {}

        keyword_workflow_map = {
            NeuroWorkflowType.STROKE_ACUTE: [
                "ACUTE STROKE", "ISCHEMIC STROKE", "TPA", "ALTEPLASE",
                "TENECTEPLASE", "THROMBECTOMY", "NIHSS", "ASPECTS",
                "LARGE VESSEL OCCLUSION", "LVO", "PENUMBRA",
                "HEMORRHAGIC STROKE", "INTRACEREBRAL HEMORRHAGE", "ICH",
                "SUBARACHNOID", "SAH", "HEMIPARESIS", "HEMIPLEGIA",
                "CEREBRAL VENOUS THROMBOSIS", "CVT", "STROKE CODE",
                "HEMATOMA", "DWI RESTRICTION", "PERFUSION MISMATCH",
            ],
            NeuroWorkflowType.STROKE_PREVENTION: [
                "SECONDARY PREVENTION", "ANTIPLATELET", "ANTICOAGULATION",
                "CAROTID STENOSIS", "ENDARTERECTOMY", "CAROTID STENTING",
                "TIA", "TRANSIENT ISCHEMIC", "ABCD2", "AFIB",
                "ATRIAL FIBRILLATION", "PFO", "PATENT FORAMEN OVALE",
                "INTRACRANIAL STENOSIS", "STROKE RISK",
            ],
            NeuroWorkflowType.DEMENTIA_EVALUATION: [
                "DEMENTIA", "ALZHEIMER", "COGNITIVE IMPAIRMENT", "MCI",
                "MOCA", "MMSE", "CDR", "FRONTOTEMPORAL", "FTD",
                "LEWY BODY", "DLB", "VASCULAR DEMENTIA",
                "AMYLOID", "TAU", "ATN", "MEMORY LOSS", "COGNITIVE DECLINE",
                "NEUROPSYCHOLOGICAL", "A\u03b242", "P-TAU", "HIPPOCAMPAL ATROPHY",
                "LECANEMAB", "DONANEMAB", "CHOLINESTERASE",
            ],
            NeuroWorkflowType.EPILEPSY_CLASSIFICATION: [
                "EPILEPSY", "SEIZURE", "EEG", "SPIKE-WAVE", "SPIKE AND WAVE",
                "STATUS EPILEPTICUS", "ANTISEIZURE", "ASM", "AED",
                "ANTICONVULSANT", "ILAE", "FOCAL SEIZURE", "TONIC-CLONIC",
                "ABSENCE", "MYOCLONIC", "TEMPORAL LOBE EPILEPSY",
                "DRUG-RESISTANT EPILEPSY", "EPILEPSY SURGERY",
                "VAGUS NERVE", "VNS", "DRAVET", "LENNOX-GASTAUT",
                "CONVULSION", "ELECTROENCEPHALOGRAPHY",
            ],
            NeuroWorkflowType.MS_MANAGEMENT: [
                "MULTIPLE SCLEROSIS", " MS ", "RRMS", "PPMS", "SPMS",
                "RELAPSING-REMITTING", "CLINICALLY ISOLATED", "CIS",
                "OLIGOCLONAL", "MCDONALD CRITERIA", "EDSS",
                "DISEASE-MODIFYING", "DMT", "NATALIZUMAB", "OCRELIZUMAB",
                "OFATUMUMAB", "FINGOLIMOD", "TECFIDERA", "NMOSD",
                "NEUROMYELITIS OPTICA", "MOGAD", "AQP4",
                "OPTIC NEURITIS", "TRANSVERSE MYELITIS",
                "DAWSON FINGERS", "PERIVENTRICULAR",
            ],
            NeuroWorkflowType.MOVEMENT_DISORDER: [
                "PARKINSON", "TREMOR", "BRADYKINESIA", "RIGIDITY",
                "LEVODOPA", "L-DOPA", "DOPAMINE AGONIST", "DBS",
                "DEEP BRAIN STIMULATION", "UPDRS", "HOEHN AND YAHR",
                "ESSENTIAL TREMOR", "DYSTONIA", "HUNTINGTON", "CHOREA",
                "ATYPICAL PARKINSONISM", "PSP", "MSA", "CBD",
                "CORTICOBASAL", "SUPRANUCLEAR PALSY",
                "MULTIPLE SYSTEM ATROPHY", "GAIT DISORDER",
                "DAT-SPECT", "DATSCAN", "MOVEMENT DISORDER",
                "TICS", "TOURETTE", "RESTLESS LEGS", "ATAXIA",
            ],
            NeuroWorkflowType.HEADACHE_DIAGNOSIS: [
                "HEADACHE", "MIGRAINE", "CLUSTER HEADACHE",
                "TENSION HEADACHE", "CGRP", "TRIPTAN", "HIT-6",
                "MIDAS", "THUNDERCLAP", "AURA", "PHOTOPHOBIA",
                "MEDICATION OVERUSE", "ANALGESIC REBOUND",
                "INTRACRANIAL HYPERTENSION", "PSEUDOTUMOR CEREBRI",
                "IIH", "PAPILLEDEMA", "CEPHALALGIA",
                "ERENUMAB", "GALCANEZUMAB", "RIMEGEPANT",
            ],
            NeuroWorkflowType.NEUROMUSCULAR_EVAL: [
                "MYASTHENIA", "ALS", "MOTOR NEURON", "AMYOTROPHIC",
                "GUILLAIN-BARRE", "GBS", "CIDP", "NEUROPATHY",
                "EMG", "NERVE CONDUCTION", "NCS", "FASCICULATION",
                "MUSCLE WEAKNESS", "PTOSIS", "DIPLOPIA",
                "ALSFRS", "RILUZOLE", "TOFERSEN",
                "SMA", "SPINAL MUSCULAR ATROPHY",
                "DUCHENNE", "MUSCULAR DYSTROPHY",
                "NEUROMUSCULAR JUNCTION", "ACETYLCHOLINE RECEPTOR",
                "MUSK", "EFGARTIGIMOD", "ECULIZUMAB",
                "PERIPHERAL NEUROPATHY", "CHARCOT-MARIE-TOOTH",
            ],
            NeuroWorkflowType.NEURO_ONCOLOGY: [
                "GLIOBLASTOMA", "GBM", "GLIOMA", "BRAIN TUMOR",
                "MENINGIOMA", "BRAIN METASTASIS", "BRAIN METS",
                "TEMOZOLOMIDE", "MGMT", "IDH MUTATION",
                "RANO CRITERIA", "PSEUDOPROGRESSION",
                "CNS LYMPHOMA", "PCNSL", "TTFIELDS",
                "STEREOTACTIC RADIOSURGERY", "GAMMA KNIFE",
                "LEPTOMENINGEAL", "NEURO-ONCOLOGY",
                "RADIATION NECROSIS", "OLIGODENDROGLIOMA",
            ],
        }

        for wf, keywords in keyword_workflow_map.items():
            for kw in keywords:
                if kw in text_upper:
                    workflow_scores[wf] = workflow_scores.get(wf, 0) + 1.0

        if not workflow_scores:
            return NeuroWorkflowType.GENERAL

        sorted_workflows = sorted(
            workflow_scores.items(), key=lambda x: x[1], reverse=True,
        )

        return sorted_workflows[0][0]

    # -- Entity Detection -------------------------------------------------

    def _detect_entities(self, question: str) -> Dict[str, List[str]]:
        """Detect neurological entities in the question text.

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
        for condition, info in NEURO_CONDITIONS.items():
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
        for drug, info in NEURO_DRUGS.items():
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
        for biomarker, info in NEURO_BIOMARKERS.items():
            if biomarker.lower() in text_lower:
                if biomarker not in entities["biomarkers"]:
                    entities["biomarkers"].append(biomarker)
                continue
            full_name = info.get("full_name", "")
            if full_name and full_name.lower() in text_lower:
                if biomarker not in entities["biomarkers"]:
                    entities["biomarkers"].append(biomarker)

        return entities

    # -- Search Strategy ---------------------------------------------------

    def _build_search_strategy(
        self,
        entities: Dict[str, List[str]],
        workflow: NeuroWorkflowType,
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
        """Choose search strategy: broad, targeted, differential, or emergency.

        Args:
            text: Original query text.
            conditions: Detected conditions.
            drugs: Detected drugs.

        Returns:
            Strategy name string.
        """
        text_upper = text.upper()

        # Emergency queries
        emergency_keywords = [
            "ACUTE STROKE", "STATUS EPILEPTICUS", "HERNIATION",
            "MYASTHENIC CRISIS", "THUNDERCLAP", "BRAIN DEATH",
            "CORD COMPRESSION", "EMERGENT", "EMERGENCY",
            "RESPIRATORY FAILURE", "INTUBAT",
        ]
        if any(kw in text_upper for kw in emergency_keywords):
            return "emergency"

        # Differential diagnosis queries
        if ("DIFFERENTIAL" in text_upper or "DIFFERENTIATE" in text_upper
                or "DISTINGUISH" in text_upper or "COMPARE" in text_upper
                or " VS " in text_upper or "VERSUS" in text_upper
                or "RULE OUT" in text_upper or "WORKUP" in text_upper):
            return "differential"

        # Targeted: specific condition + drug or single focused entity
        if (len(conditions) == 1 and len(drugs) <= 1) or (
            len(conditions) <= 1 and len(drugs) == 1
        ):
            if conditions or drugs:
                return "targeted"

        return "broad"

    # -- Sub-Question Generation -------------------------------------------

    def _generate_sub_questions(self, plan: SearchPlan) -> List[str]:
        """Generate sub-questions for comprehensive retrieval.

        Decomposes the main question into focused sub-queries based on
        the detected entities and workflow type. Enables multi-hop
        retrieval across different aspects of the neurology question.

        Args:
            plan: SearchPlan with detected entities and workflows.

        Returns:
            List of sub-question strings (typically 2-4 questions).
        """
        sub_questions: List[str] = []

        condition_label = plan.conditions[0] if plan.conditions else "the neurological condition"
        drug_label = plan.drugs[0] if plan.drugs else "the medication"
        biomarker_label = plan.biomarkers[0] if plan.biomarkers else "the biomarker"

        primary_wf = (
            plan.relevant_workflows[0] if plan.relevant_workflows
            else NeuroWorkflowType.GENERAL
        )

        # -- Pattern 1: Acute Stroke ---------------------------------
        if primary_wf == NeuroWorkflowType.STROKE_ACUTE:
            sub_questions = [
                "What are the inclusion and exclusion criteria for IV thrombolysis in acute ischemic stroke?",
                "What imaging criteria determine thrombectomy eligibility in large vessel occlusion?",
                f"What is the acute management protocol for {condition_label}?",
                "What are the blood pressure targets and monitoring requirements for acute stroke?",
            ]

        # -- Pattern 2: Stroke Prevention ----------------------------
        elif primary_wf == NeuroWorkflowType.STROKE_PREVENTION:
            sub_questions = [
                f"What secondary prevention strategies are recommended for {condition_label}?",
                "What antiplatelet or anticoagulation regimen is indicated based on stroke etiology?",
                "What vascular risk factor management targets apply for stroke prevention?",
                "What is the role of carotid intervention in stroke prevention?",
            ]

        # -- Pattern 3: Dementia Evaluation --------------------------
        elif primary_wf == NeuroWorkflowType.DEMENTIA_EVALUATION:
            sub_questions = [
                f"What diagnostic criteria apply to {condition_label}?",
                f"What biomarker and imaging workup is recommended for {condition_label}?",
                f"What is the NIA-AA ATN staging for {condition_label}?",
                f"What treatment options are available for {condition_label}?",
            ]

        # -- Pattern 4: Epilepsy Classification ----------------------
        elif primary_wf == NeuroWorkflowType.EPILEPSY_CLASSIFICATION:
            sub_questions = [
                f"How should seizures in {condition_label} be classified per ILAE 2017?",
                f"What EEG findings are expected in {condition_label}?",
                f"What antiseizure medications are first-line for {condition_label}?",
                "When should epilepsy surgery be considered for drug-resistant epilepsy?",
            ]

        # -- Pattern 5: MS Management --------------------------------
        elif primary_wf == NeuroWorkflowType.MS_MANAGEMENT:
            sub_questions = [
                f"Does the presentation meet McDonald 2017 criteria for {condition_label}?",
                f"What disease-modifying therapies are recommended for {condition_label}?",
                f"What is the monitoring plan for {drug_label} in MS?",
                "What distinguishes MS from NMOSD and MOGAD?",
            ]

        # -- Pattern 6: Movement Disorder ----------------------------
        elif primary_wf == NeuroWorkflowType.MOVEMENT_DISORDER:
            sub_questions = [
                f"What are the MDS diagnostic criteria for {condition_label}?",
                f"What red flags differentiate atypical parkinsonism from {condition_label}?",
                f"What is the evidence-based treatment algorithm for {condition_label}?",
                f"What role does {biomarker_label} play in diagnosing {condition_label}?",
            ]

        # -- Pattern 7: Headache Diagnosis ---------------------------
        elif primary_wf == NeuroWorkflowType.HEADACHE_DIAGNOSIS:
            sub_questions = [
                f"What ICHD-3 criteria apply to {condition_label}?",
                "What headache red flags require urgent neuroimaging?",
                f"What preventive therapy options exist for {condition_label}?",
                f"What acute treatment is recommended for {condition_label}?",
            ]

        # -- Pattern 8: Neuromuscular Evaluation ---------------------
        elif primary_wf == NeuroWorkflowType.NEUROMUSCULAR_EVAL:
            sub_questions = [
                f"What diagnostic criteria and workup apply to {condition_label}?",
                f"What electrodiagnostic findings (EMG/NCS) are expected in {condition_label}?",
                f"What is the evidence-based treatment for {condition_label}?",
                f"What is the role of {biomarker_label} in {condition_label} diagnosis?",
            ]

        # -- Pattern 9: Neuro-Oncology --------------------------------
        elif primary_wf == NeuroWorkflowType.NEURO_ONCOLOGY:
            sub_questions = [
                f"What is the standard of care treatment for {condition_label}?",
                f"What molecular markers guide treatment decisions in {condition_label}?",
                f"What clinical trials are available for {condition_label}?",
                "How should pseudoprogression be distinguished from true progression on imaging?",
            ]

        # -- Default ---------------------------------------------------
        else:
            sub_questions = [
                f"What are the diagnostic criteria for {condition_label}?",
                f"What is the recommended management approach for {condition_label}?",
                f"What is the current evidence base for {condition_label}?",
            ]

        return sub_questions
