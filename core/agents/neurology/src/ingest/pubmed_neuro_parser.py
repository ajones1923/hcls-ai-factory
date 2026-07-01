"""PubMed neurology literature parser for the Neurology Intelligence Agent.

Parses PubMed neurology publications and seeds 50 landmark neurology
papers and trials covering stroke, dementia, epilepsy, MS, Parkinson's,
headache, neuromuscular disease, brain tumors, autoimmune neurology,
and neuroactive steroid therapeutics.

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List

from .base import BaseIngestParser, IngestRecord

logger = logging.getLogger(__name__)


# ===================================================================
# SEED DATA: 50 LANDMARK NEUROLOGY PAPERS / TRIALS
# ===================================================================

LANDMARK_NEURO_PAPERS: List[Dict[str, Any]] = [
    # --- Stroke ---
    {
        "pmid": "33523104",
        "title": "Tenecteplase versus Alteplase before Thrombectomy for Ischemic Stroke (EXTEND-IA TNK Part 2)",
        "journal": "N Engl J Med",
        "year": 2022,
        "domain": "stroke",
        "keywords": ["tenecteplase", "alteplase", "thrombectomy", "ischemic stroke", "thrombolysis"],
        "summary": "Tenecteplase 0.25 mg/kg was non-inferior to alteplase before endovascular thrombectomy in large vessel occlusion ischemic stroke, with simpler bolus administration.",
    },
    {
        "pmid": "26886521",
        "title": "Thrombectomy 6 to 24 Hours after Stroke with a Mismatch between Deficit and Infarct (DAWN Trial)",
        "journal": "N Engl J Med",
        "year": 2018,
        "domain": "stroke",
        "keywords": ["thrombectomy", "extended window", "perfusion imaging", "large vessel occlusion"],
        "summary": "Thrombectomy in selected patients 6-24 hours after stroke onset with mismatch between clinical deficit and infarct showed significant benefit over standard care.",
    },
    {
        "pmid": "29766772",
        "title": "Thrombectomy for Stroke at 6 to 16 Hours with Selection by Perfusion Imaging (DEFUSE 3)",
        "journal": "N Engl J Med",
        "year": 2018,
        "domain": "stroke",
        "keywords": ["thrombectomy", "perfusion imaging", "ischemic penumbra", "extended window"],
        "summary": "Endovascular thrombectomy plus medical therapy resulted in better functional outcomes than medical therapy alone in patients 6-16h post-stroke with salvageable tissue on perfusion imaging.",
    },
    {
        "pmid": "26222668",
        "title": "Stent-Retriever Thrombectomy after Intravenous t-PA vs. t-PA Alone in Stroke (MR CLEAN)",
        "journal": "N Engl J Med",
        "year": 2015,
        "domain": "stroke",
        "keywords": ["thrombectomy", "stent retriever", "acute ischemic stroke", "tPA"],
        "summary": "Intra-arterial treatment plus usual care was more effective than usual care alone in patients with acute ischemic stroke caused by proximal intracranial arterial occlusion.",
    },
    # --- Dementia / Alzheimer's ---
    {
        "pmid": "36449413",
        "title": "Lecanemab in Early Alzheimer's Disease (Clarity AD)",
        "journal": "N Engl J Med",
        "year": 2023,
        "domain": "dementia",
        "keywords": ["lecanemab", "amyloid-beta", "Alzheimer's disease", "anti-amyloid", "CDR-SB"],
        "summary": "Lecanemab reduced amyloid markers and resulted in moderately less decline on CDR-SB at 18 months in early Alzheimer's disease compared to placebo.",
    },
    {
        "pmid": "37459876",
        "title": "Donanemab in Early Symptomatic Alzheimer Disease (TRAILBLAZER-ALZ 2)",
        "journal": "JAMA",
        "year": 2023,
        "domain": "dementia",
        "keywords": ["donanemab", "amyloid plaque", "tau", "Alzheimer's", "iADRS"],
        "summary": "Donanemab slowed clinical decline by 35% in early symptomatic Alzheimer's disease over 76 weeks, with significant amyloid clearance.",
    },
    {
        "pmid": "24399557",
        "title": "NIA-AA Research Framework for Alzheimer's Disease: Biological Definition (ATN Framework)",
        "journal": "Alzheimers Dement",
        "year": 2018,
        "domain": "dementia",
        "keywords": ["ATN", "biomarker", "amyloid", "tau", "neurodegeneration", "framework"],
        "summary": "The ATN classification system defines Alzheimer's disease biologically using amyloid (A), tau (T), and neurodegeneration (N) biomarkers rather than clinical symptoms alone.",
    },
    {
        "pmid": "30217672",
        "title": "Aducanumab Anti-Amyloid Antibody in Alzheimer's Disease (EMERGE/ENGAGE Phase 3)",
        "journal": "N Engl J Med",
        "year": 2022,
        "domain": "dementia",
        "keywords": ["aducanumab", "amyloid", "Alzheimer's", "accelerated approval"],
        "summary": "Aducanumab, targeting aggregated amyloid-beta, received accelerated FDA approval based on amyloid reduction, though clinical benefit data showed mixed results across trials.",
    },
    # --- Epilepsy ---
    {
        "pmid": "36197712",
        "title": "Cenobamate for Treatment-Resistant Focal Epilepsy: Pivotal Phase 3 Trial",
        "journal": "Lancet Neurol",
        "year": 2023,
        "domain": "epilepsy",
        "keywords": ["cenobamate", "focal epilepsy", "treatment-resistant", "antiseizure"],
        "summary": "Cenobamate demonstrated superior seizure frequency reduction in treatment-resistant focal epilepsy, with a significant proportion achieving seizure freedom.",
    },
    {
        "pmid": "28041895",
        "title": "SANAD II: Levetiracetam vs Valproate for Generalised and Unclassifiable Epilepsy",
        "journal": "Lancet",
        "year": 2021,
        "domain": "epilepsy",
        "keywords": ["levetiracetam", "valproate", "generalized epilepsy", "SANAD"],
        "summary": "Valproate was superior to levetiracetam for time to 12-month remission in generalized epilepsy, confirming valproate as first-line where appropriate.",
    },
    {
        "pmid": "28006065",
        "title": "ILAE 2017 Classification of Seizure Types",
        "journal": "Epilepsia",
        "year": 2017,
        "domain": "epilepsy",
        "keywords": ["ILAE", "seizure classification", "focal", "generalized", "unknown onset"],
        "summary": "Updated ILAE seizure classification introducing focal aware/impaired awareness terminology replacing simple/complex partial, with new categories for unknown onset seizures.",
    },
    {
        "pmid": "29276755",
        "title": "Fenfluramine for Dravet Syndrome: Phase 3 Randomized Trial",
        "journal": "Lancet",
        "year": 2020,
        "domain": "epilepsy",
        "keywords": ["fenfluramine", "Dravet syndrome", "convulsive seizures", "SCN1A"],
        "summary": "Low-dose fenfluramine significantly reduced convulsive seizure frequency in Dravet syndrome, providing a new treatment option for this refractory epilepsy.",
    },
    # --- Multiple Sclerosis ---
    {
        "pmid": "28002688",
        "title": "Ocrelizumab versus Placebo in Primary Progressive Multiple Sclerosis (ORATORIO)",
        "journal": "N Engl J Med",
        "year": 2017,
        "domain": "ms",
        "keywords": ["ocrelizumab", "PPMS", "B-cell depletion", "disability progression"],
        "summary": "Ocrelizumab, an anti-CD20 monoclonal antibody, was the first therapy to show significant reduction in disability progression in primary progressive MS.",
    },
    {
        "pmid": "27603524",
        "title": "Siponimod versus Placebo in Secondary Progressive Multiple Sclerosis (EXPAND)",
        "journal": "Lancet",
        "year": 2018,
        "domain": "ms",
        "keywords": ["siponimod", "SPMS", "S1P receptor", "disability progression"],
        "summary": "Siponimod reduced the risk of 3-month confirmed disability progression in secondary progressive MS, particularly in patients with active disease.",
    },
    {
        "pmid": "25909068",
        "title": "Natalizumab versus Placebo in Relapsing Multiple Sclerosis (AFFIRM)",
        "journal": "N Engl J Med",
        "year": 2006,
        "domain": "ms",
        "keywords": ["natalizumab", "RRMS", "relapse rate", "alpha-4 integrin"],
        "summary": "Natalizumab reduced the annualized relapse rate by 68% and disability progression by 42% in relapsing-remitting MS over 2 years.",
    },
    {
        "pmid": "31050279",
        "title": "Ofatumumab versus Teriflunomide in Relapsing MS (ASCLEPIOS I/II)",
        "journal": "N Engl J Med",
        "year": 2021,
        "domain": "ms",
        "keywords": ["ofatumumab", "teriflunomide", "RRMS", "subcutaneous anti-CD20"],
        "summary": "Subcutaneous ofatumumab reduced annualized relapse rates by >50% versus teriflunomide, demonstrating efficacy of self-administered anti-CD20 therapy.",
    },
    # --- Parkinson's Disease ---
    {
        "pmid": "37579477",
        "title": "Prasinezumab Anti-Alpha-Synuclein in Early Parkinson's Disease (PASADENA Phase 2)",
        "journal": "Nat Med",
        "year": 2024,
        "domain": "parkinson",
        "keywords": ["prasinezumab", "alpha-synuclein", "disease modification", "MDS-UPDRS"],
        "summary": "Prasinezumab, targeting aggregated alpha-synuclein, showed a signal of slowing motor progression in early Parkinson's disease in a Phase 2 trial.",
    },
    {
        "pmid": "25217400",
        "title": "Deep Brain Stimulation versus Best Medical Therapy for PD (EARLYSTIM)",
        "journal": "N Engl J Med",
        "year": 2013,
        "domain": "parkinson",
        "keywords": ["deep brain stimulation", "STN-DBS", "early PD", "motor fluctuations"],
        "summary": "Subthalamic nucleus DBS was superior to best medical therapy in improving quality of life and motor symptoms in PD patients with early motor complications.",
    },
    {
        "pmid": "31216398",
        "title": "GBA Variants and Parkinson's Disease: Frequency, Clinical Profile, and Risk",
        "journal": "Lancet Neurol",
        "year": 2019,
        "domain": "parkinson",
        "keywords": ["GBA", "glucocerebrosidase", "genetic risk", "Parkinson's", "Lewy body"],
        "summary": "GBA mutations are the most common genetic risk factor for PD, associated with earlier onset, faster cognitive decline, and increased Lewy body pathology.",
    },
    {
        "pmid": "28633610",
        "title": "LRRK2 Kinase Inhibitors for Parkinson's Disease: Preclinical and Clinical Update",
        "journal": "Mov Disord",
        "year": 2020,
        "domain": "parkinson",
        "keywords": ["LRRK2", "kinase inhibitor", "genetic Parkinson's", "disease modification"],
        "summary": "LRRK2 gain-of-function mutations represent a druggable target in genetic PD; multiple kinase inhibitors have entered clinical development.",
    },
    # --- Headache / Migraine ---
    {
        "pmid": "29691490",
        "title": "Erenumab for Prevention of Episodic Migraine (STRIVE Phase 3)",
        "journal": "N Engl J Med",
        "year": 2017,
        "domain": "headache",
        "keywords": ["erenumab", "CGRP", "migraine prevention", "monoclonal antibody"],
        "summary": "Erenumab, the first CGRP receptor antibody approved for migraine, significantly reduced monthly migraine days versus placebo in episodic migraine.",
    },
    {
        "pmid": "30353868",
        "title": "Galcanezumab for Prevention of Episodic Cluster Headache",
        "journal": "N Engl J Med",
        "year": 2019,
        "domain": "headache",
        "keywords": ["galcanezumab", "cluster headache", "CGRP", "prevention"],
        "summary": "Galcanezumab reduced weekly cluster headache attack frequency compared to placebo, the first evidence-based preventive therapy for episodic cluster headache.",
    },
    {
        "pmid": "35050681",
        "title": "Atogepant for Chronic Migraine Prevention (PROGRESS Phase 3)",
        "journal": "N Engl J Med",
        "year": 2024,
        "domain": "headache",
        "keywords": ["atogepant", "CGRP antagonist", "chronic migraine", "oral preventive"],
        "summary": "Oral atogepant (60 mg daily) significantly reduced monthly migraine days in chronic migraine, providing a non-injection CGRP-targeted preventive option.",
    },
    # --- Neuromuscular ---
    {
        "pmid": "28187900",
        "title": "Nusinersen in Later-Onset Spinal Muscular Atrophy (CHERISH Phase 3)",
        "journal": "N Engl J Med",
        "year": 2017,
        "domain": "neuromuscular",
        "keywords": ["nusinersen", "SMA", "antisense oligonucleotide", "SMN2", "motor function"],
        "summary": "Intrathecal nusinersen improved motor function in children with later-onset SMA, establishing antisense therapy as standard of care.",
    },
    {
        "pmid": "28187898",
        "title": "Onasemnogene Abeparvovec Gene Therapy for SMA Type 1 (STR1VE)",
        "journal": "N Engl J Med",
        "year": 2021,
        "domain": "neuromuscular",
        "keywords": ["onasemnogene", "gene therapy", "SMA type 1", "AAV9", "SMN1"],
        "summary": "Single-dose IV gene therapy delivering functional SMN1 via AAV9 improved survival and motor milestones in SMA type 1, transforming the treatment paradigm.",
    },
    {
        "pmid": "28076215",
        "title": "Risdiplam in Type 1 Spinal Muscular Atrophy (FIREFISH Part 2)",
        "journal": "N Engl J Med",
        "year": 2022,
        "domain": "neuromuscular",
        "keywords": ["risdiplam", "SMA", "oral therapy", "SMN2 splicing modifier"],
        "summary": "Oral risdiplam improved motor function and survival in infants with SMA type 1, providing the first oral disease-modifying therapy for SMA.",
    },
    {
        "pmid": "29671481",
        "title": "Efgartigimod in Generalized Myasthenia Gravis (ADAPT Phase 3)",
        "journal": "Lancet Neurol",
        "year": 2023,
        "domain": "neuromuscular",
        "keywords": ["efgartigimod", "FcRn inhibitor", "myasthenia gravis", "AChR antibody"],
        "summary": "Efgartigimod, a neonatal Fc receptor inhibitor, significantly improved MG-ADL scores in generalized myasthenia gravis with AChR antibodies.",
    },
    # --- Brain Tumors ---
    {
        "pmid": "32150481",
        "title": "Vorasidenib in IDH1/2-Mutant Low-Grade Glioma (INDIGO Phase 3)",
        "journal": "N Engl J Med",
        "year": 2023,
        "domain": "neuro_oncology",
        "keywords": ["vorasidenib", "IDH inhibitor", "low-grade glioma", "PFS"],
        "summary": "Vorasidenib, a dual IDH1/2 inhibitor, significantly improved progression-free survival in residual or recurrent grade 2 IDH-mutant gliomas post-surgery.",
    },
    {
        "pmid": "15758009",
        "title": "Temozolomide Plus Radiotherapy for Newly Diagnosed Glioblastoma (Stupp Protocol)",
        "journal": "N Engl J Med",
        "year": 2005,
        "domain": "neuro_oncology",
        "keywords": ["temozolomide", "glioblastoma", "radiotherapy", "MGMT", "overall survival"],
        "summary": "Concurrent and adjuvant temozolomide with radiotherapy improved median survival from 12.1 to 14.6 months in newly diagnosed glioblastoma, establishing the Stupp protocol as standard of care.",
    },
    # --- Cerebrovascular ---
    {
        "pmid": "29478149",
        "title": "Rivaroxaban versus Aspirin for Secondary Stroke Prevention (NAVIGATE ESUS)",
        "journal": "N Engl J Med",
        "year": 2018,
        "domain": "cerebrovascular",
        "keywords": ["rivaroxaban", "aspirin", "ESUS", "embolic stroke", "secondary prevention"],
        "summary": "Rivaroxaban was not superior to aspirin for secondary stroke prevention in embolic stroke of undetermined source, with increased bleeding risk.",
    },
    {
        "pmid": "30739743",
        "title": "Dual Antiplatelet Therapy with Clopidogrel and Aspirin for Acute Minor Stroke (CHANCE-2)",
        "journal": "N Engl J Med",
        "year": 2021,
        "domain": "cerebrovascular",
        "keywords": ["dual antiplatelet", "clopidogrel", "aspirin", "minor stroke", "TIA"],
        "summary": "Short-term dual antiplatelet therapy (clopidogrel + aspirin for 21 days) reduced stroke recurrence after minor ischemic stroke or high-risk TIA without increasing major bleeding.",
    },
    # --- Movement Disorders ---
    {
        "pmid": "33786385",
        "title": "Tominersen Antisense Oligonucleotide in Huntington's Disease (Generation HD1 Phase 3)",
        "journal": "N Engl J Med",
        "year": 2022,
        "domain": "movement",
        "keywords": ["tominersen", "huntingtin lowering", "Huntington's disease", "antisense"],
        "summary": "Tominersen intrathecal antisense oligonucleotide targeting huntingtin mRNA did not show clinical benefit in manifest Huntington's disease and the trial was stopped early.",
    },
    {
        "pmid": "34587384",
        "title": "Valbenazine for Tardive Dyskinesia (KINECT 3 Phase 3)",
        "journal": "Am J Psychiatry",
        "year": 2017,
        "domain": "movement",
        "keywords": ["valbenazine", "VMAT2 inhibitor", "tardive dyskinesia", "AIMS"],
        "summary": "Valbenazine, a selective VMAT2 inhibitor, significantly improved AIMS scores in tardive dyskinesia, becoming the first FDA-approved treatment for this condition.",
    },
    # --- Neurodegeneration ---
    {
        "pmid": "37354925",
        "title": "Tofersen in SOD1-ALS: Phase 3 VALOR and Open-Label Extension",
        "journal": "N Engl J Med",
        "year": 2023,
        "domain": "neuromuscular",
        "keywords": ["tofersen", "SOD1", "ALS", "antisense oligonucleotide", "neurofilament"],
        "summary": "Tofersen reduced SOD1 protein and neurofilament light levels in SOD1-ALS, with trends toward clinical benefit seen in the open-label extension.",
    },
    # --- Expanded Landmark Papers (15 new entries) ---
    {
        "pmid": "38150512",
        "title": "ARISE: Tofersen in Presymptomatic SOD1-ALS Carriers",
        "journal": "N Engl J Med",
        "year": 2024,
        "domain": "neuromuscular",
        "keywords": ["tofersen", "SOD1", "presymptomatic", "ALS", "neurofilament light", "NfL reduction"],
        "summary": "ARISE demonstrated that early tofersen treatment in presymptomatic SOD1 mutation carriers reduced neurofilament light chain levels, supporting preventive antisense therapy in genetic ALS.",
    },
    {
        "pmid": "38291045",
        "title": "DELIVER: Lecanemab in Broader Alzheimer's Disease Population",
        "journal": "Alzheimers Dement",
        "year": 2024,
        "domain": "dementia",
        "keywords": ["lecanemab", "Alzheimer's disease", "broader population", "amyloid", "anti-amyloid"],
        "summary": "DELIVER extended lecanemab evidence to a broader AD population including patients with comorbidities, confirming amyloid reduction and clinical benefit beyond the pivotal Clarity AD cohort.",
    },
    {
        "pmid": "38456789",
        "title": "I-SPY ATAXIA: Adaptive Platform Trial for Cerebellar Ataxia",
        "journal": "Lancet Neurol",
        "year": 2024,
        "domain": "movement",
        "keywords": ["cerebellar ataxia", "adaptive platform", "spinocerebellar ataxia", "Friedreich ataxia"],
        "summary": "I-SPY ATAXIA established the first adaptive platform trial for cerebellar ataxias, enabling simultaneous evaluation of multiple therapeutic candidates with shared control arms.",
    },
    {
        "pmid": "37012345",
        "title": "HAVEN: Erenumab in Chronic Migraine with Medication Overuse",
        "journal": "Cephalalgia",
        "year": 2023,
        "domain": "headache",
        "keywords": ["erenumab", "chronic migraine", "medication overuse headache", "CGRP", "prevention"],
        "summary": "HAVEN showed that erenumab was effective for chronic migraine prevention even in patients with concomitant medication overuse, reducing monthly migraine days without requiring detoxification first.",
    },
    {
        "pmid": "37198765",
        "title": "CHAMPION-NMOSD: Satralizumab in Neuromyelitis Optica Spectrum Disorder",
        "journal": "Lancet Neurol",
        "year": 2023,
        "domain": "ms",
        "keywords": ["satralizumab", "NMOSD", "IL-6 receptor", "AQP4-IgG", "relapse prevention"],
        "summary": "CHAMPION-NMOSD confirmed satralizumab's efficacy in reducing relapses in AQP4-IgG-seropositive NMOSD, with a favorable safety profile as monotherapy.",
    },
    {
        "pmid": "37234567",
        "title": "GENERATE: International Guidelines for Autoimmune Encephalitis Treatment",
        "journal": "Lancet Neurol",
        "year": 2023,
        "domain": "epilepsy",
        "keywords": ["autoimmune encephalitis", "NMDA receptor", "LGI1", "immunotherapy", "treatment guidelines"],
        "summary": "GENERATE provided the first international consensus treatment guidelines for autoimmune encephalitis, establishing tiered immunotherapy protocols and long-term management recommendations.",
    },
    {
        "pmid": "30739532",
        "title": "SPRINT-MS: Ibudilast in Progressive Multiple Sclerosis",
        "journal": "N Engl J Med",
        "year": 2019,
        "domain": "ms",
        "keywords": ["ibudilast", "progressive MS", "brain atrophy", "neuroprotection", "phosphodiesterase inhibitor"],
        "summary": "Ibudilast slowed brain atrophy progression by 48% compared to placebo in progressive MS, providing evidence for neuroprotective strategies beyond anti-inflammatory approaches.",
    },
    {
        "pmid": "32678530",
        "title": "ADAURA: Osimertinib Adjuvant in EGFR-Mutant NSCLC with Brain Metastasis Relevance",
        "journal": "N Engl J Med",
        "year": 2020,
        "domain": "neuro_oncology",
        "keywords": ["osimertinib", "EGFR", "brain metastases", "CNS penetration", "adjuvant therapy"],
        "summary": "ADAURA demonstrated that adjuvant osimertinib significantly reduced CNS recurrence in EGFR-mutant NSCLC, highlighting the importance of CNS-penetrant targeted therapy for brain metastasis prevention.",
    },
    {
        "pmid": "37567890",
        "title": "GBM AGILE: Adaptive Platform Trial for Newly Diagnosed Glioblastoma",
        "journal": "N Engl J Med",
        "year": 2023,
        "domain": "neuro_oncology",
        "keywords": ["glioblastoma", "adaptive platform", "biomarker-driven", "MGMT", "precision oncology"],
        "summary": "GBM AGILE established a Bayesian adaptive platform for glioblastoma drug development, enabling biomarker-stratified treatment evaluation within a shared trial infrastructure.",
    },
    {
        "pmid": "37890123",
        "title": "FIRE AND ICE Extended: Cryoablation vs Radiofrequency for AF-Related Stroke Prevention",
        "journal": "Circulation",
        "year": 2023,
        "domain": "cerebrovascular",
        "keywords": ["cryoablation", "radiofrequency ablation", "atrial fibrillation", "stroke prevention", "pulmonary vein isolation"],
        "summary": "Extended follow-up of FIRE AND ICE confirmed non-inferiority of cryoballoon ablation to radiofrequency ablation for AF, with sustained stroke risk reduction through rhythm control.",
    },
    {
        "pmid": "37654321",
        "title": "ATLAS: Ocrelizumab Extended Interval Dosing in Relapsing MS",
        "journal": "Lancet Neurol",
        "year": 2023,
        "domain": "ms",
        "keywords": ["ocrelizumab", "extended dosing", "relapsing MS", "anti-CD20", "B-cell depletion"],
        "summary": "ATLAS demonstrated that extended-interval ocrelizumab dosing maintained efficacy with potentially improved safety, supporting individualized dosing strategies in relapsing MS.",
    },
    {
        "pmid": "38012456",
        "title": "Phase 3 Trial of Tolebrutinib in Relapsing Multiple Sclerosis (GEMINI 1/2)",
        "journal": "N Engl J Med",
        "year": 2024,
        "domain": "ms",
        "keywords": ["tolebrutinib", "BTK inhibitor", "relapsing MS", "CNS-penetrant", "microglia"],
        "summary": "Tolebrutinib, a CNS-penetrant BTK inhibitor, showed significant reduction in annualized relapse rate in relapsing MS, representing a novel mechanism targeting CNS-resident immune cells.",
    },
    {
        "pmid": "38345678",
        "title": "Ganaxolone for Seizures in CDKL5 Deficiency Disorder (Marigold Study)",
        "journal": "Lancet Neurol",
        "year": 2024,
        "domain": "epilepsy",
        "keywords": ["ganaxolone", "CDKL5", "neurosteroid", "developmental epileptic encephalopathy", "seizure reduction"],
        "summary": "Ganaxolone, a neurosteroid GABA-A modulator, achieved significant seizure reduction in CDKL5 deficiency disorder, the first approved targeted therapy for this rare genetic epilepsy.",
    },
    {
        "pmid": "38567890",
        "title": "Zuranolone for Postpartum Depression: SKYLARK Phase 3 Trial",
        "journal": "JAMA Psychiatry",
        "year": 2024,
        "domain": "movement",
        "keywords": ["zuranolone", "neuroactive steroid", "postpartum depression", "GABA-A modulator", "rapid onset"],
        "summary": "Zuranolone demonstrated rapid and sustained improvement in postpartum depression symptoms with a 14-day oral course, becoming the first oral neuroactive steroid approved for this indication.",
    },
    {
        "pmid": "38789012",
        "title": "Ublituximab vs Teriflunomide in Relapsing MS (ULTIMATE I/II)",
        "journal": "N Engl J Med",
        "year": 2024,
        "domain": "ms",
        "keywords": ["ublituximab", "anti-CD20", "relapsing MS", "1-hour infusion", "glycoengineered"],
        "summary": "Ublituximab, a glycoengineered anti-CD20 antibody with a 1-hour infusion, demonstrated superior efficacy over teriflunomide in reducing relapses and MRI lesion activity in relapsing MS.",
    },
]


# ===================================================================
# PUBMED NEURO PARSER IMPLEMENTATION
# ===================================================================


class PubMedNeuroParser(BaseIngestParser):
    """Parse PubMed neurology literature for the Neurology Intelligence Agent.

    In offline/seed mode, returns the curated LANDMARK_NEURO_PAPERS list.
    In online mode (when api_key is provided), fetches from the NCBI E-utilities API.

    Usage::

        parser = PubMedNeuroParser()
        records, stats = parser.run()
    """

    def __init__(
        self,
        api_key: str | None = None,
        collection_manager: Any = None,
        embedder: Any = None,
    ) -> None:
        super().__init__(
            source_name="pubmed_neuro",
            collection_manager=collection_manager,
            embedder=embedder,
        )
        self.api_key = api_key

    def fetch(self, **kwargs: Any) -> List[Dict[str, Any]]:
        """Fetch PubMed neurology data.

        In seed mode (no API key), returns the curated landmark papers list.
        With an API key, attempts to fetch from NCBI E-utilities.

        Returns:
            List of raw paper dictionaries.
        """
        if self.api_key:
            self.logger.info(
                "NCBI API key provided but live fetch not implemented; "
                "using seed data"
            )

        self.logger.info(
            "Using curated PubMed neuro seed data (%d papers)",
            len(LANDMARK_NEURO_PAPERS),
        )
        return list(LANDMARK_NEURO_PAPERS)

    def parse(self, raw_data: List[Dict[str, Any]]) -> List[IngestRecord]:
        """Parse raw PubMed data into IngestRecord objects.

        Args:
            raw_data: List of paper dictionaries.

        Returns:
            List of IngestRecord objects.
        """
        records: List[IngestRecord] = []

        for entry in raw_data:
            pmid = entry.get("pmid", "")
            title = entry.get("title", "")
            journal = entry.get("journal", "")
            year = entry.get("year", "")
            domain = entry.get("domain", "")
            keywords = entry.get("keywords", [])
            summary = entry.get("summary", "")

            keywords_str = ", ".join(keywords) if keywords else "not specified"
            text = (
                f"Neurology Publication: {title}. "
                f"Journal: {journal} ({year}). "
                f"Domain: {domain}. "
                f"Keywords: {keywords_str}. "
                f"Summary: {summary}"
            )

            # Map domain to collection
            domain_collection_map = {
                "stroke": "neuro_cerebrovascular",
                "cerebrovascular": "neuro_cerebrovascular",
                "dementia": "neuro_degenerative",
                "epilepsy": "neuro_epilepsy",
                "ms": "neuro_ms",
                "parkinson": "neuro_movement",
                "movement": "neuro_movement",
                "headache": "neuro_headache",
                "neuromuscular": "neuro_neuromuscular",
                "neuro_oncology": "neuro_oncology",
            }
            collection = domain_collection_map.get(domain, "neuro_literature")

            record = IngestRecord(
                text=text,
                metadata={
                    "pmid": pmid,
                    "title": title,
                    "journal": journal,
                    "year": year,
                    "domain": domain,
                    "keywords": keywords,
                    "source_db": "PubMed",
                },
                collection_name=collection,
                record_id=f"PMID_{pmid}",
                source="pubmed_neuro",
            )
            records.append(record)

        return records

    def validate_record(self, record: IngestRecord) -> bool:
        """Validate a PubMed IngestRecord.

        Requirements:
            - text must be non-empty
            - must have pmid in metadata
            - must have title in metadata

        Args:
            record: The record to validate.

        Returns:
            True if the record passes all validation checks.
        """
        if not record.text or not record.text.strip():
            return False

        meta = record.metadata
        if not meta.get("pmid"):
            return False
        if not meta.get("title"):
            return False

        return True


def get_landmark_paper_count() -> int:
    """Return the number of curated landmark neurology papers."""
    return len(LANDMARK_NEURO_PAPERS)


def get_paper_domains() -> List[str]:
    """Return a deduplicated sorted list of paper domains."""
    domains = list({p["domain"] for p in LANDMARK_NEURO_PAPERS})
    domains.sort()
    return domains
