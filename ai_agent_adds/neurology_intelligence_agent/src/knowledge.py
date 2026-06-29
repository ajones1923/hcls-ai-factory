"""Neurology Intelligence Agent -- Domain Knowledge Base.

Comprehensive neurology knowledge covering 10+ disease domains, key drugs,
genes, clinical scales, neurodegenerative diseases, epilepsy syndromes,
stroke protocols, headache classifications, MS DMT tiers, and
domain-specific terminology for neurology RAG-based clinical decision support.

Author: Adam Jones
Date: March 2026
"""

from typing import Any, Dict, List


# ===================================================================
# KNOWLEDGE BASE VERSION
# ===================================================================

KNOWLEDGE_VERSION: Dict[str, Any] = {
    "version": "2.0.0",
    "last_updated": "2026-03-22",
    "revision_notes": "Major expansion -- 10+ neurological domains, 42 drugs, "
                      "38+ genes, 10 clinical scales, 15 neurodegenerative diseases, "
                      "12 epilepsy syndromes, 6 stroke protocols, 8 headache "
                      "classifications, 3 MS DMT tiers, sleep and neuroimmunology domains.",
    "sources": [
        "AAN Practice Guidelines",
        "AHA/ASA Stroke Guidelines",
        "ILAE Classification and Treatment Guidelines",
        "McDonald Criteria for MS (2017 revision)",
        "MDS-UPDRS for Parkinson's Disease",
        "ICHD-3 Headache Classification",
        "WHO 2021 CNS Tumor Classification",
        "NIA-AA ATN Framework for Alzheimer's Disease",
        "ACNS EEG Terminology Standards",
        "PubMed / MEDLINE",
        "AASM Sleep Disorder Classification",
        "Autoimmune Encephalitis Consensus Criteria",
    ],
    "counts": {
        "disease_domains": 10,
        "drugs": 42,
        "genes": 38,
        "clinical_scales": 10,
        "imaging_protocols": 55,
        "eeg_patterns": 35,
        "landmark_papers": 35,
        "neurodegenerative_diseases": 15,
        "epilepsy_syndromes": 12,
        "stroke_protocols": 6,
        "headache_classifications": 8,
        "ms_dmt_tiers": 3,
    },
}


# ===================================================================
# NEUROLOGICAL DISEASE DOMAINS
# ===================================================================

NEURO_DOMAINS: Dict[str, Dict[str, Any]] = {
    "cerebrovascular": {
        "description": "Cerebrovascular disease including ischemic stroke, "
                       "hemorrhagic stroke, TIA, and cerebral venous thrombosis.",
        "key_conditions": [
            "Large vessel occlusion stroke",
            "Small vessel ischemic disease",
            "Cardioembolic stroke",
            "Intracerebral hemorrhage",
            "Subarachnoid hemorrhage",
            "Transient ischemic attack",
            "Cerebral venous sinus thrombosis",
            "Moyamoya disease",
            "Carotid artery stenosis",
        ],
        "key_scales": ["NIHSS", "ASPECTS", "mRS", "GCS"],
    },
    "degenerative": {
        "description": "Neurodegenerative diseases including Alzheimer's, "
                       "frontotemporal dementia, Lewy body disease, and prion diseases.",
        "key_conditions": [
            "Alzheimer's disease",
            "Frontotemporal dementia (bvFTD, PPA)",
            "Lewy body dementia",
            "Vascular dementia",
            "Progressive supranuclear palsy",
            "Corticobasal degeneration",
            "Normal pressure hydrocephalus",
            "Creutzfeldt-Jakob disease",
            "Posterior cortical atrophy",
        ],
        "key_scales": ["MoCA", "MMSE", "CDR"],
    },
    "epilepsy": {
        "description": "Epilepsy and seizure disorders including focal, "
                       "generalized, and epileptic encephalopathies.",
        "key_conditions": [
            "Temporal lobe epilepsy",
            "Juvenile myoclonic epilepsy",
            "Childhood absence epilepsy",
            "Dravet syndrome",
            "Lennox-Gastaut syndrome",
            "West syndrome (infantile spasms)",
            "Status epilepticus",
            "Autoimmune encephalitis with seizures",
            "PNES (psychogenic non-epileptic seizures)",
        ],
        "key_scales": ["seizure_frequency", "Engel_class"],
    },
    "movement": {
        "description": "Movement disorders including Parkinson's disease, "
                       "essential tremor, dystonia, and Huntington disease.",
        "key_conditions": [
            "Parkinson's disease",
            "Essential tremor",
            "Dystonia (focal, generalized)",
            "Huntington disease",
            "Multiple system atrophy",
            "Progressive supranuclear palsy",
            "Wilson disease",
            "Tardive dyskinesia",
            "Restless legs syndrome",
        ],
        "key_scales": ["MDS-UPDRS Part III", "Hoehn & Yahr"],
    },
    "ms": {
        "description": "Multiple sclerosis and related CNS inflammatory demyelinating diseases.",
        "key_conditions": [
            "Relapsing-remitting MS",
            "Secondary progressive MS",
            "Primary progressive MS",
            "Clinically isolated syndrome",
            "Neuromyelitis optica spectrum disorder",
            "MOG antibody-associated disease",
            "Acute disseminated encephalomyelitis",
        ],
        "key_scales": ["EDSS"],
    },
    "headache": {
        "description": "Primary and secondary headache disorders per ICHD-3.",
        "key_conditions": [
            "Migraine without aura",
            "Migraine with aura",
            "Chronic migraine",
            "Tension-type headache",
            "Cluster headache",
            "Medication overuse headache",
            "New daily persistent headache",
            "Trigeminal neuralgia",
            "Idiopathic intracranial hypertension",
        ],
        "key_scales": ["HIT-6", "MIDAS"],
    },
    "neuromuscular": {
        "description": "Neuromuscular diseases including ALS, myasthenia gravis, "
                       "neuropathies, and myopathies.",
        "key_conditions": [
            "Amyotrophic lateral sclerosis",
            "Myasthenia gravis",
            "Guillain-Barre syndrome",
            "Chronic inflammatory demyelinating polyneuropathy",
            "Spinal muscular atrophy",
            "Duchenne muscular dystrophy",
            "Charcot-Marie-Tooth disease",
            "Inclusion body myositis",
            "Lambert-Eaton myasthenic syndrome",
        ],
        "key_scales": ["ALSFRS-R", "MG-ADL"],
    },
    "neuro_oncology": {
        "description": "Primary and metastatic CNS tumors with molecular classification.",
        "key_conditions": [
            "Glioblastoma (IDH-wildtype)",
            "Astrocytoma (IDH-mutant)",
            "Oligodendroglioma (IDH-mutant, 1p19q-codeleted)",
            "Meningioma",
            "Brain metastases",
            "Primary CNS lymphoma",
            "Medulloblastoma",
            "Ependymoma",
            "Vestibular schwannoma",
        ],
        "key_scales": ["KPS", "RANO"],
    },
    "sleep_neurology": {
        "description": "Sleep-related neurological disorders including narcolepsy, "
                       "REM sleep behavior disorder, parasomnias, and circadian rhythm disorders.",
        "key_conditions": [
            "Narcolepsy type 1 (with cataplexy)",
            "Narcolepsy type 2 (without cataplexy)",
            "REM sleep behavior disorder",
            "Obstructive sleep apnea with neurological sequelae",
            "Central sleep apnea",
            "Restless legs syndrome / Willis-Ekbom disease",
            "Periodic limb movement disorder",
            "Idiopathic hypersomnia",
            "Kleine-Levin syndrome",
            "Fatal familial insomnia",
            "Circadian rhythm sleep-wake disorders",
        ],
        "key_scales": ["ESS", "MSLT", "PSG"],
    },
    "neuroimmunology": {
        "description": "Autoimmune and antibody-mediated neurological disorders "
                       "including autoimmune encephalitis, neurosarcoidosis, and "
                       "CNS vasculitis.",
        "key_conditions": [
            "Anti-NMDA receptor encephalitis",
            "LGI1 antibody encephalitis",
            "CASPR2 antibody encephalitis",
            "MOG antibody-associated disease",
            "AQP4-IgG neuromyelitis optica spectrum disorder",
            "Autoimmune GFAP astrocytopathy",
            "Paraneoplastic cerebellar degeneration",
            "Stiff-person syndrome",
            "Neurosarcoidosis",
            "Primary angiitis of the CNS",
            "Hashimoto encephalopathy (SREAT)",
            "Bickerstaff brainstem encephalitis",
        ],
        "key_scales": ["mRS", "CASE_score"],
    },
}


# ===================================================================
# NEURODEGENERATIVE DISEASES (15 diseases)
# ===================================================================

NEURODEGENERATIVE_DISEASES: Dict[str, Dict[str, Any]] = {
    "alzheimers_early_onset": {
        "name": "Early-Onset Alzheimer's Disease",
        "genes": ["APP", "PSEN1", "PSEN2"],
        "biomarkers": ["amyloid-beta 42", "phospho-tau 181", "phospho-tau 217",
                       "amyloid PET", "tau PET", "NfL"],
        "clinical_scales": ["MoCA", "MMSE", "CDR", "ADAS-Cog"],
        "staging": {
            "preclinical": "Biomarker-positive, cognitively normal",
            "prodromal_MCI": "Mild cognitive impairment with AD biomarkers",
            "mild_dementia": "CDR 0.5-1, functional decline in complex tasks",
            "moderate_dementia": "CDR 2, requires assistance with ADLs",
            "severe_dementia": "CDR 3, fully dependent",
        },
        "treatments": ["lecanemab", "donanemab", "donepezil", "memantine"],
        "key_trials": ["CLARITY AD", "TRAILBLAZER-ALZ 2", "DIAN-TU"],
    },
    "alzheimers_late_onset": {
        "name": "Late-Onset Alzheimer's Disease",
        "genes": ["APOE", "TREM2", "CLU", "BIN1", "ABCA7"],
        "biomarkers": ["amyloid-beta 42/40 ratio", "phospho-tau 181",
                       "phospho-tau 217", "GFAP", "NfL", "amyloid PET"],
        "clinical_scales": ["MoCA", "MMSE", "CDR", "ADAS-Cog", "FAQ"],
        "staging": {
            "preclinical": "A+T-N-, asymptomatic",
            "early_MCI": "A+T+N-, subjective cognitive decline",
            "MCI_due_to_AD": "A+T+N+/-, objective impairment, preserved ADLs",
            "mild_dementia": "A+T+N+, functional decline",
            "moderate_dementia": "Requires significant support",
            "severe_dementia": "Fully dependent, non-verbal",
        },
        "treatments": ["lecanemab", "donanemab", "donepezil", "rivastigmine",
                       "galantamine", "memantine"],
        "key_trials": ["CLARITY AD", "TRAILBLAZER-ALZ 2", "AHEAD 3-45"],
    },
    "frontotemporal_bvFTD": {
        "name": "Behavioral Variant Frontotemporal Dementia",
        "genes": ["MAPT", "GRN", "C9orf72"],
        "biomarkers": ["NfL", "progranulin (GRN mutations)", "TDP-43 (CSF/plasma)"],
        "clinical_scales": ["FBI", "CBI-R", "FTLD-CDR"],
        "staging": {
            "possible": "3 of 6 behavioral features",
            "probable": "Functional decline + frontal/temporal atrophy on imaging",
            "definite": "Histopathological confirmation or known pathogenic mutation",
        },
        "treatments": ["SSRIs for behavioral symptoms", "trazodone",
                       "no disease-modifying therapy approved"],
        "key_trials": ["ALLFTD", "GENFI", "ARTFL-LEFFTDS"],
    },
    "semantic_PPA": {
        "name": "Semantic Variant Primary Progressive Aphasia",
        "genes": ["GRN", "C9orf72", "MAPT"],
        "biomarkers": ["NfL", "anterior temporal atrophy on MRI"],
        "clinical_scales": ["WAB", "BNT", "PPVT"],
        "staging": {
            "early": "Anomia, impaired single-word comprehension",
            "moderate": "Surface dyslexia, prosopagnosia, behavioral changes",
            "severe": "Global semantic loss, mutism",
        },
        "treatments": ["speech-language therapy", "SSRIs for behavioral symptoms"],
        "key_trials": ["ALLFTD", "LEFFTDS"],
    },
    "nonfluent_PPA": {
        "name": "Nonfluent/Agrammatic Variant Primary Progressive Aphasia",
        "genes": ["GRN", "MAPT"],
        "biomarkers": ["NfL", "left posterior fronto-insular atrophy"],
        "clinical_scales": ["WAB", "NAT"],
        "staging": {
            "early": "Effortful speech, agrammatism, preserved comprehension",
            "moderate": "Apraxia of speech worsens, limited verbal output",
            "severe": "Mutism, may develop PSP or CBD features",
        },
        "treatments": ["speech-language therapy",
                       "no disease-modifying therapy approved"],
        "key_trials": ["ALLFTD"],
    },
    "dementia_lewy_bodies": {
        "name": "Dementia with Lewy Bodies",
        "genes": ["GBA1", "SNCA", "APOE"],
        "biomarkers": ["DAT scan (reduced striatal uptake)",
                       "MIBG myocardial scintigraphy", "polysomnography (RBD)",
                       "alpha-synuclein seed amplification assay"],
        "clinical_scales": ["MoCA", "MMSE", "CDR", "NPI", "DLB Consortium criteria"],
        "staging": {
            "prodromal": "Isolated RBD, mild cognitive impairment",
            "mild": "Fluctuating cognition, visual hallucinations, parkinsonism",
            "moderate": "Significant cognitive and motor impairment",
            "severe": "Fully dependent, recurrent falls",
        },
        "treatments": ["donepezil", "rivastigmine", "carbidopa-levodopa (cautious)",
                       "pimavanserin", "melatonin for RBD"],
        "key_trials": ["DIAMOND-Lewy"],
    },
    "parkinsons_disease": {
        "name": "Parkinson's Disease",
        "genes": ["LRRK2", "GBA1", "SNCA", "PARK2", "PINK1", "PARK7"],
        "biomarkers": ["DAT scan", "alpha-synuclein seed amplification assay",
                       "NfL", "skin biopsy (phospho-alpha-synuclein)"],
        "clinical_scales": ["MDS-UPDRS", "Hoehn & Yahr", "PDQ-39", "MoCA"],
        "staging": {
            "preclinical": "Prodromal features (RBD, hyposmia, constipation)",
            "stage_1": "Unilateral motor symptoms",
            "stage_2": "Bilateral motor symptoms, no balance impairment",
            "stage_3": "Bilateral disease with postural instability",
            "stage_4": "Severe disability, able to walk or stand unassisted",
            "stage_5": "Wheelchair or bed-bound",
        },
        "treatments": ["levodopa/carbidopa", "pramipexole", "ropinirole",
                       "rasagiline", "safinamide", "amantadine", "DBS"],
        "key_trials": ["LEAP", "STEADY-PD III", "PASADENA", "SPARK"],
    },
    "als_sporadic": {
        "name": "Sporadic Amyotrophic Lateral Sclerosis",
        "genes": ["no single causative gene", "risk loci: UNC13A, C9orf72 expansions in ~10%"],
        "biomarkers": ["NfL (CSF and serum)", "phospho-NfH",
                       "EMG (active denervation)", "TDP-43"],
        "clinical_scales": ["ALSFRS-R", "SVC (slow vital capacity)", "MUNIX"],
        "staging": {
            "stage_1": "Symptom onset in one region",
            "stage_2": "Spread to second region",
            "stage_3": "Spread to third region",
            "stage_4": "Respiratory failure requiring ventilation",
        },
        "treatments": ["riluzole", "edaravone", "AMX0035 (sodium phenylbutyrate/taurursodiol)",
                       "multidisciplinary care", "NIV", "PEG"],
        "key_trials": ["CENTAUR", "PHOENIX", "ATLAS"],
    },
    "als_familial": {
        "name": "Familial Amyotrophic Lateral Sclerosis",
        "genes": ["SOD1", "C9orf72", "TARDBP", "FUS", "TBK1", "NEK1"],
        "biomarkers": ["NfL", "phospho-NfH", "genetic testing",
                       "EMG (active denervation)"],
        "clinical_scales": ["ALSFRS-R", "SVC", "MUNIX"],
        "staging": {
            "presymptomatic": "Gene carrier, no symptoms",
            "stage_1": "Symptom onset in one region",
            "stage_2": "Spread to second region",
            "stage_3": "Spread to third region",
            "stage_4": "Respiratory failure",
        },
        "treatments": ["tofersen (SOD1)", "riluzole", "edaravone",
                       "genetic counseling"],
        "key_trials": ["VALOR (tofersen)", "ATLAS (presymptomatic SOD1)"],
    },
    "huntingtons": {
        "name": "Huntington Disease",
        "genes": ["HTT (CAG repeat expansion >= 36)"],
        "biomarkers": ["NfL", "mutant huntingtin protein",
                       "caudate atrophy on MRI"],
        "clinical_scales": ["UHDRS", "TFC", "TMS", "MoCA"],
        "staging": {
            "premanifest": "Gene-positive, no motor symptoms (TFC 13)",
            "stage_1": "Early motor signs, functionally independent (TFC 11-13)",
            "stage_2": "Functional decline, may work with accommodations (TFC 7-10)",
            "stage_3": "Cannot work or manage finances (TFC 3-6)",
            "stage_4": "Requires substantial care (TFC 1-2)",
            "stage_5": "Total care required (TFC 0)",
        },
        "treatments": ["tetrabenazine", "deutetrabenazine", "valbenazine",
                       "SSRIs", "antipsychotics for chorea"],
        "key_trials": ["GENERATION HD1", "PRECISION-HD"],
    },
    "msa_c": {
        "name": "Multiple System Atrophy - Cerebellar Type (MSA-C)",
        "genes": ["COQ2 (rare)", "GBA1 (risk factor)"],
        "biomarkers": ["NfL (elevated)", "hot-cross-bun sign on MRI",
                       "cerebellar/pontine atrophy"],
        "clinical_scales": ["UMSARS", "ICARS", "COMPASS-31"],
        "staging": {
            "possible": "Sporadic cerebellar ataxia with autonomic dysfunction",
            "probable": "Cerebellar ataxia + urogenital dysfunction + cerebellar atrophy",
            "definite": "Neuropathological confirmation (GCI with alpha-synuclein)",
        },
        "treatments": ["midodrine/droxidopa (orthostatic hypotension)",
                       "fludrocortisone", "physical therapy"],
        "key_trials": ["MOVEMENT-MSA"],
    },
    "msa_p": {
        "name": "Multiple System Atrophy - Parkinsonian Type (MSA-P)",
        "genes": ["COQ2 (rare)", "GBA1 (risk factor)"],
        "biomarkers": ["NfL (elevated)", "putaminal rim sign on MRI",
                       "putaminal atrophy"],
        "clinical_scales": ["UMSARS", "MDS-UPDRS", "COMPASS-31"],
        "staging": {
            "possible": "Sporadic parkinsonism (poor levodopa response) + autonomic dysfunction",
            "probable": "Parkinsonism + urogenital dysfunction + putaminal abnormalities",
            "definite": "Neuropathological confirmation (GCI with alpha-synuclein)",
        },
        "treatments": ["levodopa trial (typically poor response)",
                       "midodrine/droxidopa", "physical therapy"],
        "key_trials": ["MOVEMENT-MSA"],
    },
    "psp": {
        "name": "Progressive Supranuclear Palsy",
        "genes": ["MAPT (H1 haplotype risk)", "LRRK2 (rare)"],
        "biomarkers": ["NfL (elevated)", "midbrain atrophy (hummingbird sign)",
                       "tau PET (limited sensitivity)"],
        "clinical_scales": ["PSP-RS", "SEADL", "MoCA"],
        "staging": {
            "suggestive": "Gradual onset, predominantly tau-related features",
            "possible": "Vertical supranuclear gaze palsy or slowing + postural instability",
            "probable": "Vertical supranuclear gaze palsy + falls in first year",
            "definite": "Neuropathological confirmation (4R tauopathy)",
        },
        "treatments": ["levodopa trial (usually poor response)",
                       "coenzyme Q10 (insufficient evidence)", "physical therapy",
                       "no approved disease-modifying therapy"],
        "key_trials": ["PROSPECT-PSP", "PASSPORT"],
    },
    "cbd": {
        "name": "Corticobasal Degeneration",
        "genes": ["MAPT (H1 haplotype risk)"],
        "biomarkers": ["NfL", "asymmetric cortical atrophy on MRI",
                       "tau PET (limited availability)"],
        "clinical_scales": ["CBD-FRS", "MDS-UPDRS", "MoCA"],
        "staging": {
            "possible": "Asymmetric limb rigidity/akinesia + cortical features (apraxia, alien limb)",
            "probable": "Above + imaging support (asymmetric atrophy, hypometabolism)",
            "definite": "Neuropathological confirmation (4R tauopathy with astrocytic plaques)",
        },
        "treatments": ["levodopa trial (usually poor response)",
                       "clonazepam for myoclonus", "botulinum toxin for dystonia",
                       "occupational therapy"],
        "key_trials": ["TauRx LUCIDITY (general tauopathy)"],
    },
    "prion_disease": {
        "name": "Creutzfeldt-Jakob Disease (Prion Disease)",
        "genes": ["PRNP"],
        "biomarkers": ["14-3-3 protein (CSF)", "RT-QuIC (CSF)",
                       "total tau (CSF, markedly elevated)", "cortical ribboning on DWI MRI"],
        "clinical_scales": ["MRC Prion Disease Rating Scale"],
        "staging": {
            "prodromal": "Psychiatric symptoms, fatigue, cognitive changes",
            "early": "Rapidly progressive dementia, myoclonus",
            "advanced": "Akinetic mutism, periodic sharp wave complexes on EEG",
            "terminal": "Decerebrate state",
        },
        "treatments": ["no effective treatment", "palliative care",
                       "doxycycline (investigational)", "quinacrine (negative trials)"],
        "key_trials": ["PRION-1", "UK National Prion Monitoring Cohort"],
    },
}


# ===================================================================
# EPILEPSY SYNDROMES (12 syndromes)
# ===================================================================

EPILEPSY_SYNDROMES: Dict[str, Dict[str, Any]] = {
    "dravet": {
        "name": "Dravet Syndrome",
        "gene": "SCN1A (>80% have pathogenic variant)",
        "seizure_types": ["febrile seizures (prolonged)", "generalized tonic-clonic",
                          "myoclonic", "focal", "absence", "status epilepticus"],
        "eeg_pattern": "Normal early; later generalized spike-wave, photosensitivity",
        "first_line_asm": ["stiripentol", "cannabidiol (Epidiolex)", "fenfluramine (Fintepla)"],
        "contraindicated_asm": ["carbamazepine", "oxcarbazepine", "phenytoin",
                                 "lamotrigine", "vigabatrin"],
        "surgical_candidate": False,
    },
    "lennox_gastaut": {
        "name": "Lennox-Gastaut Syndrome",
        "gene": "Multiple (often structural/unknown etiology)",
        "seizure_types": ["tonic", "atonic (drop attacks)", "atypical absence",
                          "myoclonic", "generalized tonic-clonic"],
        "eeg_pattern": "Slow spike-and-wave (<2.5 Hz), paroxysmal fast activity in sleep",
        "first_line_asm": ["valproate", "lamotrigine", "rufinamide",
                           "cannabidiol", "clobazam"],
        "contraindicated_asm": ["carbamazepine (may worsen atypical absence/myoclonic)"],
        "surgical_candidate": False,  # Callosotomy considered for drop attacks
    },
    "west_syndrome": {
        "name": "West Syndrome (Infantile Spasms)",
        "gene": "Multiple (TSC1/TSC2, ARX, CDKL5, STXBP1, others)",
        "seizure_types": ["epileptic (infantile) spasms", "clusters on awakening"],
        "eeg_pattern": "Hypsarrhythmia (chaotic high-amplitude slow waves with multifocal spikes)",
        "first_line_asm": ["ACTH (adrenocorticotropic hormone)", "vigabatrin (if TSC)",
                           "prednisolone"],
        "contraindicated_asm": ["carbamazepine", "oxcarbazepine"],
        "surgical_candidate": True,  # If focal etiology identified
    },
    "jme": {
        "name": "Juvenile Myoclonic Epilepsy",
        "gene": "EFHC1, GABRA1, CLCN2 (polygenic)",
        "seizure_types": ["myoclonic jerks (morning predominance)",
                          "generalized tonic-clonic", "absence (30%)"],
        "eeg_pattern": "4-6 Hz generalized polyspike-and-wave, photosensitivity",
        "first_line_asm": ["valproate", "levetiracetam", "lamotrigine"],
        "contraindicated_asm": ["carbamazepine", "oxcarbazepine", "phenytoin",
                                 "vigabatrin", "tiagabine"],
        "surgical_candidate": False,
    },
    "cae": {
        "name": "Childhood Absence Epilepsy",
        "gene": "GABRG2, CACNA1H (polygenic)",
        "seizure_types": ["typical absence seizures (many per day)"],
        "eeg_pattern": "3 Hz generalized spike-and-wave, hyperventilation-provoked",
        "first_line_asm": ["ethosuximide", "valproate", "lamotrigine"],
        "contraindicated_asm": ["carbamazepine", "oxcarbazepine", "phenytoin",
                                 "vigabatrin", "tiagabine"],
        "surgical_candidate": False,
    },
    "tle_hippocampal": {
        "name": "Temporal Lobe Epilepsy with Hippocampal Sclerosis",
        "gene": "Usually acquired (febrile seizures, infection, trauma); rare genetic forms",
        "seizure_types": ["focal impaired awareness (deja vu, epigastric rising, automatisms)",
                          "focal to bilateral tonic-clonic"],
        "eeg_pattern": "Temporal sharp waves/spikes, temporal intermittent rhythmic delta activity (TIRDA)",
        "first_line_asm": ["carbamazepine", "oxcarbazepine", "lamotrigine",
                           "levetiracetam", "lacosamide"],
        "contraindicated_asm": [],
        "surgical_candidate": True,  # Anterior temporal lobectomy - 60-80% seizure freedom
    },
    "bects": {
        "name": "Benign Epilepsy with Centrotemporal Spikes (BECTS/Rolandic Epilepsy)",
        "gene": "GRIN2A (rare), ELP4 (susceptibility)",
        "seizure_types": ["focal motor (face/arm), often nocturnal",
                          "secondary generalization (rare)"],
        "eeg_pattern": "Centrotemporal (rolandic) spikes, activated by drowsiness/sleep",
        "first_line_asm": ["often no treatment needed (self-limited)",
                           "levetiracetam if treatment required",
                           "carbamazepine", "oxcarbazepine"],
        "contraindicated_asm": [],
        "surgical_candidate": False,
    },
    "focal_cortical_dysplasia": {
        "name": "Focal Cortical Dysplasia-associated Epilepsy",
        "gene": "MTOR, DEPDC5, NPRL2, NPRL3, SLC35A2 (somatic/germline)",
        "seizure_types": ["focal aware", "focal impaired awareness",
                          "focal to bilateral tonic-clonic"],
        "eeg_pattern": "Focal spikes/sharp waves, continuous spikes in slow sleep (some types)",
        "first_line_asm": ["carbamazepine", "oxcarbazepine", "lamotrigine",
                           "levetiracetam", "lacosamide"],
        "contraindicated_asm": [],
        "surgical_candidate": True,  # Often best outcome with complete resection
    },
    "tuberous_sclerosis_epilepsy": {
        "name": "Tuberous Sclerosis Complex - Epilepsy",
        "gene": "TSC1 (hamartin), TSC2 (tuberin)",
        "seizure_types": ["infantile spasms", "focal seizures", "generalized seizures"],
        "eeg_pattern": "Multifocal spikes, hypsarrhythmia (in infantile spasms)",
        "first_line_asm": ["vigabatrin (first-line for infantile spasms in TSC)",
                           "everolimus (mTOR inhibitor)", "cannabidiol"],
        "contraindicated_asm": [],
        "surgical_candidate": True,  # If a single tuber is the epileptogenic focus
    },
    "progressive_myoclonic": {
        "name": "Progressive Myoclonic Epilepsies",
        "gene": "CSTB (Unverricht-Lundborg), EPM2A/NHLRC1 (Lafora), CLN genes (NCL), MERRF (mitochondrial)",
        "seizure_types": ["action myoclonus", "generalized tonic-clonic",
                          "absence (some forms)"],
        "eeg_pattern": "Generalized spike-wave, photoparoxysmal response, progressive background slowing",
        "first_line_asm": ["valproate", "levetiracetam", "clonazepam",
                           "piracetam (for myoclonus)"],
        "contraindicated_asm": ["phenytoin (may worsen cerebellar symptoms)",
                                 "carbamazepine", "lamotrigine (Lafora)"],
        "surgical_candidate": False,
    },
    "cdkl5_epilepsy": {
        "name": "CDKL5 Deficiency Disorder",
        "gene": "CDKL5",
        "seizure_types": ["infantile-onset tonic seizures", "epileptic spasms",
                          "hypermotor seizures", "myoclonic seizures"],
        "eeg_pattern": "Variable; may show hypsarrhythmia or multifocal discharges",
        "first_line_asm": ["ganaxolone", "cannabidiol", "clobazam", "vigabatrin"],
        "contraindicated_asm": [],
        "surgical_candidate": False,
    },
    "glut1_deficiency": {
        "name": "GLUT1 Deficiency Syndrome",
        "gene": "SLC2A1",
        "seizure_types": ["generalized tonic-clonic", "absence", "myoclonic",
                          "focal seizures"],
        "eeg_pattern": "Generalized 2.5-4 Hz spike-wave, worsened by fasting",
        "first_line_asm": ["ketogenic diet (first-line, treats underlying metabolic defect)",
                           "triheptanoin"],
        "contraindicated_asm": ["valproate (inhibits fatty acid oxidation)",
                                 "phenobarbital (impairs GLUT1 function)"],
        "surgical_candidate": False,
    },
}


# ===================================================================
# STROKE PROTOCOLS (6 protocols)
# ===================================================================

STROKE_PROTOCOLS: Dict[str, Dict[str, Any]] = {
    "tpa_eligibility": {
        "time_window": "0-4.5 hours from symptom onset",
        "criteria": [
            "Age >= 18 years",
            "Clinical diagnosis of ischemic stroke causing measurable neurological deficit",
            "Onset of symptoms < 4.5 hours before treatment",
            "CT/MRI excluding hemorrhage",
        ],
        "contraindications": [
            "Active internal bleeding",
            "History of intracranial hemorrhage",
            "Recent intracranial/spinal surgery (< 3 months)",
            "Intracranial neoplasm, AVM, or aneurysm",
            "Platelet count < 100,000",
            "INR > 1.7 or PT > 15 seconds",
            "Blood glucose < 50 mg/dL",
            "Severe uncontrolled hypertension (>185/110 despite treatment)",
            "CT showing multilobar infarction (hypodensity > 1/3 cerebral hemisphere)",
            "3-4.5 hour extension exclusions: age > 80, NIHSS > 25, oral anticoagulant use, diabetes + prior stroke",
        ],
        "key_trials": ["NINDS", "ECASS III", "IST-3", "WAKE-UP"],
    },
    "thrombectomy_dawn": {
        "time_window": "6-24 hours from last known well",
        "criteria": [
            "ICA or M1 MCA occlusion",
            "Age >= 18 years",
            "NIHSS >= 10",
            "Pre-stroke mRS 0-1",
            "Clinical-imaging mismatch per DAWN criteria:",
            "  - Age >= 80: NIHSS >= 10 and infarct core < 21 mL",
            "  - Age < 80: NIHSS >= 10 and infarct core < 31 mL",
            "  - Age < 80: NIHSS >= 20 and infarct core 31-51 mL",
        ],
        "contraindications": [
            "Large established infarct core exceeding mismatch criteria",
            "Pre-stroke mRS >= 2",
            "No target vessel occlusion",
            "Rapidly improving symptoms",
        ],
        "key_trials": ["DAWN"],
    },
    "thrombectomy_defuse3": {
        "time_window": "6-16 hours from last known well",
        "criteria": [
            "ICA or M1 MCA occlusion",
            "Age 18-90 years",
            "NIHSS >= 6",
            "Pre-stroke mRS 0-2",
            "Ischemic core volume < 70 mL",
            "Mismatch ratio >= 1.8",
            "Mismatch volume >= 15 mL",
            "Tmax > 6s volume (penumbra) identified on CTP or MR perfusion",
        ],
        "contraindications": [
            "Infarct core >= 70 mL",
            "Mismatch ratio < 1.8",
            "No target vessel occlusion",
            "Pre-stroke mRS >= 3",
        ],
        "key_trials": ["DEFUSE 3"],
    },
    "hemorrhagic_management": {
        "time_window": "Immediate (hyperacute phase 0-6 hours critical)",
        "criteria": [
            "CT-confirmed intracerebral hemorrhage",
            "Blood pressure management: target SBP < 140 mmHg (INTERACT2/ATACH-2)",
            "Reversal of anticoagulation if applicable",
            "ICH volume estimation (ABC/2 method)",
            "ICH Score calculation for prognosis",
            "Assessment for surgical evacuation (posterior fossa, lobar with deterioration)",
        ],
        "contraindications": [
            "Surgery contraindicated: deep basal ganglia hemorrhage (STICH trial)",
            "Comfort care only if ICH Score suggests futility",
        ],
        "key_trials": ["INTERACT2", "ATACH-2", "MISTIE III", "STICH", "STICH II",
                       "FASTEST", "TICH-2"],
    },
    "sah_management": {
        "time_window": "Immediate; vasospasm window days 3-14",
        "criteria": [
            "CT-confirmed subarachnoid hemorrhage (or xanthochromia on LP)",
            "Hunt & Hess / WFNS grading",
            "CTA or DSA to identify aneurysm",
            "Secure aneurysm within 24 hours (coil vs clip)",
            "Nimodipine 60 mg q4h for 21 days (vasospasm prophylaxis)",
            "EVD if hydrocephalus",
            "Target euvolemia, avoid hyponatremia",
            "Transcranial Doppler monitoring for vasospasm",
        ],
        "contraindications": [
            "Avoid hypotension before aneurysm secured",
            "Avoid antifibrinolytics beyond 72 hours (thrombotic risk)",
        ],
        "key_trials": ["ISAT", "BRAT", "CONSCIOUS-1", "NEWTON"],
    },
    "tia_workup": {
        "time_window": "Urgent (within 24-48 hours)",
        "criteria": [
            "Transient neurological symptoms resolved within 24 hours",
            "ABCD2 score for risk stratification",
            "Brain MRI with DWI (to rule out acute infarct)",
            "Neurovascular imaging (CTA or MRA of head and neck)",
            "Cardiac evaluation: ECG, telemetry, echocardiography",
            "Extended cardiac monitoring (14-30 day) if no AF detected",
            "Lab work: CBC, BMP, lipid panel, HbA1c, coagulation",
            "Dual antiplatelet therapy (aspirin + clopidogrel) for 21 days if high risk",
            "Statin initiation (high-intensity)",
        ],
        "contraindications": [
            "Dual antiplatelet > 90 days (bleeding risk outweighs benefit)",
            "Anticoagulation not indicated unless AF or cardioembolic source",
        ],
        "key_trials": ["CHANCE", "POINT", "SOCRATES", "NAVIGATE ESUS",
                       "CRYSTAL AF"],
    },
}


# ===================================================================
# HEADACHE CLASSIFICATIONS (8 types per ICHD-3)
# ===================================================================

HEADACHE_CLASSIFICATIONS: Dict[str, Dict[str, Any]] = {
    "migraine_without_aura": {
        "ichd3_code": "1.1",
        "criteria": [
            "At least 5 attacks",
            "Headache lasting 4-72 hours (untreated)",
            "At least 2 of: unilateral, pulsating, moderate/severe intensity, aggravation by routine physical activity",
            "At least 1 of: nausea/vomiting, photophobia and phonophobia",
        ],
        "acute_treatment": ["triptans (sumatriptan, rizatriptan, eletriptan)",
                            "NSAIDs (ibuprofen, naproxen)", "rimegepant",
                            "ubrogepant", "lasmiditan", "dihydroergotamine"],
        "preventive_treatment": ["topiramate", "propranolol", "amitriptyline",
                                  "valproate", "CGRP mAbs (erenumab, galcanezumab, fremanezumab)",
                                  "atogepant", "rimegepant"],
    },
    "migraine_with_aura": {
        "ichd3_code": "1.2",
        "criteria": [
            "At least 2 attacks",
            "One or more fully reversible aura symptoms: visual, sensory, speech/language, motor, brainstem, retinal",
            "At least 3 of: one aura symptom spreads gradually over >= 5 min, two or more symptoms occur in succession, "
            "each symptom lasts 5-60 min, at least one symptom is unilateral, aura accompanied/followed by headache within 60 min",
        ],
        "acute_treatment": ["triptans (after aura resolves)", "NSAIDs",
                            "gepants (rimegepant, ubrogepant)", "lasmiditan"],
        "preventive_treatment": ["same as migraine without aura",
                                  "avoid estrogen-containing contraceptives (stroke risk)"],
    },
    "chronic_migraine": {
        "ichd3_code": "1.3",
        "criteria": [
            "Headache on >= 15 days/month for > 3 months",
            "At least 8 days/month meeting migraine criteria or responsive to triptan/ergot",
            "Not better accounted for by another ICHD-3 diagnosis",
        ],
        "acute_treatment": ["limit acute medication to < 10-15 days/month to avoid MOH",
                            "triptans", "NSAIDs", "gepants"],
        "preventive_treatment": ["onabotulinumtoxinA (PREEMPT protocol)",
                                  "CGRP mAbs", "topiramate", "atogepant",
                                  "eptinezumab (quarterly IV)"],
    },
    "cluster": {
        "ichd3_code": "3.1",
        "criteria": [
            "At least 5 attacks",
            "Severe unilateral orbital/supraorbital/temporal pain lasting 15-180 min",
            "At least 1 ipsilateral: conjunctival injection, lacrimation, nasal congestion, rhinorrhea, "
            "forehead/facial sweating, miosis, ptosis, eyelid edema",
            "Frequency: 1 every other day to 8/day",
        ],
        "acute_treatment": ["high-flow oxygen (12-15 L/min via non-rebreather)",
                            "sumatriptan SC 6 mg", "zolmitriptan nasal spray"],
        "preventive_treatment": ["verapamil (first-line)", "galcanezumab",
                                  "lithium", "melatonin", "occipital nerve block",
                                  "short-course prednisone (bridge)"],
    },
    "tension_type": {
        "ichd3_code": "2.1/2.2/2.3",
        "criteria": [
            "Headache lasting 30 min to 7 days",
            "At least 2 of: bilateral, pressing/tightening (non-pulsating), mild-moderate intensity, "
            "not aggravated by routine physical activity",
            "No nausea/vomiting (mild nausea permitted in chronic form)",
            "No more than one of photophobia or phonophobia",
        ],
        "acute_treatment": ["NSAIDs (ibuprofen, naproxen, aspirin)",
                            "acetaminophen", "combination analgesics"],
        "preventive_treatment": ["amitriptyline", "nortriptyline", "mirtazapine",
                                  "physical therapy", "cognitive behavioral therapy",
                                  "biofeedback"],
    },
    "medication_overuse": {
        "ichd3_code": "8.2",
        "criteria": [
            "Headache on >= 15 days/month in a patient with pre-existing headache disorder",
            "Regular overuse of acute headache medication for > 3 months",
            "Triptans, ergots, opioids, combination analgesics: >= 10 days/month",
            "Simple analgesics: >= 15 days/month",
        ],
        "acute_treatment": ["withdrawal/detoxification of overused medication",
                            "bridge therapy (nerve block, short-course steroids, DHE protocol)"],
        "preventive_treatment": ["start preventive before or during withdrawal",
                                  "topiramate", "onabotulinumtoxinA", "CGRP mAbs",
                                  "atogepant"],
    },
    "trigeminal_neuralgia": {
        "ichd3_code": "13.1.1",
        "criteria": [
            "Recurrent paroxysms of unilateral facial pain in trigeminal distribution (V2/V3 > V1)",
            "Lasting a fraction of a second to 2 minutes",
            "Severe, electric shock-like, shooting, stabbing, or sharp quality",
            "Precipitated by innocuous stimuli (trigger zones)",
        ],
        "acute_treatment": ["carbamazepine (first-line)", "oxcarbazepine"],
        "preventive_treatment": ["carbamazepine", "oxcarbazepine", "baclofen",
                                  "lamotrigine", "microvascular decompression (MVD)",
                                  "percutaneous rhizotomy",
                                  "stereotactic radiosurgery (Gamma Knife)"],
    },
    "new_daily_persistent": {
        "ichd3_code": "4.10",
        "criteria": [
            "Daily and unremitting headache from onset (or < 3 days to becoming unremitting)",
            "Clearly remembered and unambiguous onset",
            "Headache present for > 3 months",
            "Not better accounted for by another ICHD-3 diagnosis",
        ],
        "acute_treatment": ["NSAIDs", "triptans (if migraine phenotype)",
                            "nerve blocks"],
        "preventive_treatment": ["topiramate", "gabapentin", "amitriptyline",
                                  "onabotulinumtoxinA", "CGRP mAbs",
                                  "doxycycline (anti-inflammatory)",
                                  "often treatment-refractory"],
    },
}


# ===================================================================
# MS DMT TIERS (3 tiers)
# ===================================================================

MS_DMT_TIERS: Dict[str, Dict[str, Any]] = {
    "platform": {
        "interferons": {
            "mechanism": "Immunomodulation via type I interferon receptor; shifts cytokine balance from Th1/Th17 to Th2",
            "efficacy_vs_platform": "Reference (30% ARR reduction vs placebo)",
            "safety_considerations": ["Flu-like symptoms", "injection site reactions",
                                       "hepatotoxicity", "depression", "lymphopenia"],
            "monitoring": ["CBC q3-6 months", "LFTs q3-6 months",
                           "thyroid function annually"],
        },
        "glatiramer": {
            "mechanism": "Myelin basic protein analogue; promotes Th2 shift and regulatory T cells",
            "efficacy_vs_platform": "Comparable to interferon beta (~30% ARR reduction)",
            "safety_considerations": ["Injection site reactions", "lipoatrophy",
                                       "immediate post-injection reaction (IPIR)"],
            "monitoring": ["No routine blood monitoring required"],
        },
        "teriflunomide": {
            "mechanism": "Dihydroorotate dehydrogenase inhibitor; blocks de novo pyrimidine synthesis in activated lymphocytes",
            "efficacy_vs_platform": "Comparable to interferons (~30-36% ARR reduction)",
            "safety_considerations": ["Hepatotoxicity", "teratogenicity (Pregnancy Category X)",
                                       "hair thinning", "GI symptoms", "peripheral neuropathy"],
            "monitoring": ["LFTs monthly for 6 months then periodically", "CBC",
                           "blood pressure", "TB screening",
                           "pregnancy test (requires accelerated elimination with cholestyramine if planning pregnancy)"],
        },
        "dimethyl_fumarate": {
            "mechanism": "Nrf2 pathway activation; anti-oxidant and anti-inflammatory effects",
            "efficacy_vs_platform": "Slightly superior to platform (~44-53% ARR reduction vs placebo)",
            "safety_considerations": ["Flushing", "GI symptoms (nausea, diarrhea, abdominal pain)",
                                       "lymphopenia (risk of PML if severe/prolonged)",
                                       "hepatotoxicity"],
            "monitoring": ["CBC q6 months (lymphocyte count critical)",
                           "LFTs periodically", "hold if lymphocytes < 500 for > 6 months"],
        },
    },
    "moderate": {
        "fingolimod": {
            "mechanism": "S1P receptor modulator; sequesters lymphocytes in lymph nodes",
            "efficacy_vs_platform": "Superior to interferon beta-1a IM (~54% ARR reduction vs placebo)",
            "safety_considerations": ["First-dose bradycardia/heart block", "macular edema",
                                       "infections (VZV, cryptococcus)", "PML (rare)",
                                       "hepatotoxicity", "skin cancer risk"],
            "monitoring": ["First-dose cardiac monitoring (6 hours)", "CBC q3-6 months",
                           "LFTs", "ophthalmology (baseline + 3-4 months)",
                           "VZV titer (vaccinate if negative before starting)",
                           "dermatology screening annually"],
        },
        "cladribine": {
            "mechanism": "Purine analogue; selective depletion of T and B lymphocytes",
            "efficacy_vs_platform": "Superior (~57% ARR reduction vs placebo)",
            "safety_considerations": ["Lymphopenia (expected, nadir at 2-3 months)",
                                       "herpes zoster reactivation", "malignancy concern (theoretical)",
                                       "teratogenicity"],
            "monitoring": ["CBC (lymphocyte count) before each course",
                           "do not start next course if lymphocytes < 800",
                           "cancer screening", "pregnancy test",
                           "VZV vaccination if seronegative"],
        },
    },
    "high_efficacy": {
        "natalizumab": {
            "mechanism": "Alpha-4 integrin antibody; blocks lymphocyte migration across BBB",
            "efficacy_vs_platform": "Highly superior (~68% ARR reduction vs placebo)",
            "safety_considerations": ["Progressive multifocal leukoencephalopathy (PML)",
                                       "risk stratified by JCV antibody index",
                                       "hepatotoxicity", "infusion reactions",
                                       "rebound disease activity on discontinuation"],
            "monitoring": ["JCV antibody index q6 months", "MRI q3-6 months (PML surveillance)",
                           "LFTs", "consider extended interval dosing (EID) if JCV+"],
        },
        "ocrelizumab": {
            "mechanism": "Anti-CD20 monoclonal antibody; depletes CD20+ B cells",
            "efficacy_vs_platform": "Highly superior (~46-47% ARR reduction vs IFN beta-1a; also effective in PPMS)",
            "safety_considerations": ["Infusion reactions", "infections (URI, herpes)",
                                       "hypogammaglobulinemia (long-term)", "hepatitis B reactivation",
                                       "possible malignancy risk (breast cancer signal in trials)"],
            "monitoring": ["CBC, immunoglobulin levels q6-12 months",
                           "hepatitis B/C screening", "vaccination status (complete before starting)",
                           "infusion q6 months", "cancer screening"],
        },
        "ofatumumab": {
            "mechanism": "Anti-CD20 monoclonal antibody (subcutaneous, targets distinct CD20 epitope)",
            "efficacy_vs_platform": "Highly superior (~50-58% ARR reduction vs teriflunomide)",
            "safety_considerations": ["Injection site reactions", "infections",
                                       "hypogammaglobulinemia", "hepatitis B reactivation"],
            "monitoring": ["CBC, immunoglobulin levels q6-12 months",
                           "hepatitis B/C screening", "vaccination status"],
        },
        "alemtuzumab": {
            "mechanism": "Anti-CD52 monoclonal antibody; profound lymphocyte depletion with immune reconstitution",
            "efficacy_vs_platform": "Highly superior (~49-55% ARR reduction vs IFN beta-1a SC)",
            "safety_considerations": ["Secondary autoimmunity (thyroid 30-40%, ITP 2-3%, anti-GBM nephropathy <1%)",
                                       "infusion reactions (severe)", "infections",
                                       "stroke/cervicocephalic arterial dissection (rare)",
                                       "malignancy risk", "REMS program"],
            "monitoring": ["CBC monthly for 48 months after last dose",
                           "creatinine and urinalysis monthly for 48 months",
                           "thyroid function q3 months for 48 months",
                           "skin exam annually", "HPV screening"],
        },
    },
}


# ===================================================================
# NEUROLOGY DRUGS
# ===================================================================

NEURO_DRUGS: List[Dict[str, Any]] = [
    # Stroke
    {"name": "alteplase", "class": "thrombolytic", "domain": "cerebrovascular", "mechanism": "plasminogen activator"},
    {"name": "tenecteplase", "class": "thrombolytic", "domain": "cerebrovascular", "mechanism": "plasminogen activator"},
    {"name": "aspirin", "class": "antiplatelet", "domain": "cerebrovascular", "mechanism": "COX inhibitor"},
    {"name": "clopidogrel", "class": "antiplatelet", "domain": "cerebrovascular", "mechanism": "P2Y12 inhibitor"},
    {"name": "apixaban", "class": "anticoagulant", "domain": "cerebrovascular", "mechanism": "factor Xa inhibitor"},
    {"name": "rivaroxaban", "class": "anticoagulant", "domain": "cerebrovascular", "mechanism": "factor Xa inhibitor"},
    # Dementia
    {"name": "lecanemab", "class": "anti-amyloid", "domain": "degenerative", "mechanism": "amyloid-beta antibody"},
    {"name": "donanemab", "class": "anti-amyloid", "domain": "degenerative", "mechanism": "amyloid-beta antibody"},
    {"name": "aducanumab", "class": "anti-amyloid", "domain": "degenerative", "mechanism": "amyloid-beta antibody"},
    {"name": "donepezil", "class": "cholinesterase_inhibitor", "domain": "degenerative", "mechanism": "AChE inhibitor"},
    {"name": "memantine", "class": "NMDA_antagonist", "domain": "degenerative", "mechanism": "NMDA receptor antagonist"},
    # Epilepsy
    {"name": "levetiracetam", "class": "antiseizure", "domain": "epilepsy", "mechanism": "SV2A modulator"},
    {"name": "valproate", "class": "antiseizure", "domain": "epilepsy", "mechanism": "multiple mechanisms"},
    {"name": "lamotrigine", "class": "antiseizure", "domain": "epilepsy", "mechanism": "sodium channel blocker"},
    {"name": "carbamazepine", "class": "antiseizure", "domain": "epilepsy", "mechanism": "sodium channel blocker"},
    {"name": "cenobamate", "class": "antiseizure", "domain": "epilepsy", "mechanism": "sodium channel + GABA"},
    {"name": "fenfluramine", "class": "antiseizure", "domain": "epilepsy", "mechanism": "serotonin modulator"},
    {"name": "cannabidiol", "class": "antiseizure", "domain": "epilepsy", "mechanism": "multiple mechanisms"},
    # MS
    {"name": "ocrelizumab", "class": "anti-CD20", "domain": "ms", "mechanism": "B-cell depletion"},
    {"name": "ofatumumab", "class": "anti-CD20", "domain": "ms", "mechanism": "B-cell depletion"},
    {"name": "natalizumab", "class": "integrin_inhibitor", "domain": "ms", "mechanism": "alpha-4 integrin blockade"},
    {"name": "siponimod", "class": "S1P_modulator", "domain": "ms", "mechanism": "S1P1/5 receptor modulator"},
    {"name": "dimethyl_fumarate", "class": "immunomodulator", "domain": "ms", "mechanism": "Nrf2 pathway activation"},
    {"name": "glatiramer_acetate", "class": "immunomodulator", "domain": "ms", "mechanism": "MBP analogue"},
    # Parkinson's
    {"name": "levodopa", "class": "dopamine_precursor", "domain": "movement", "mechanism": "DA precursor"},
    {"name": "pramipexole", "class": "dopamine_agonist", "domain": "movement", "mechanism": "D2/D3 agonist"},
    {"name": "ropinirole", "class": "dopamine_agonist", "domain": "movement", "mechanism": "D2/D3 agonist"},
    {"name": "safinamide", "class": "MAO-B_inhibitor", "domain": "movement", "mechanism": "MAO-B inhibitor"},
    {"name": "amantadine", "class": "NMDA_antagonist", "domain": "movement", "mechanism": "NMDA antagonist + DA release"},
    {"name": "valbenazine", "class": "VMAT2_inhibitor", "domain": "movement", "mechanism": "vesicular monoamine transporter 2 inhibitor"},
    # Headache
    {"name": "erenumab", "class": "CGRP_antibody", "domain": "headache", "mechanism": "CGRP receptor antibody"},
    {"name": "galcanezumab", "class": "CGRP_antibody", "domain": "headache", "mechanism": "CGRP ligand antibody"},
    {"name": "fremanezumab", "class": "CGRP_antibody", "domain": "headache", "mechanism": "CGRP ligand antibody"},
    {"name": "atogepant", "class": "CGRP_antagonist", "domain": "headache", "mechanism": "oral CGRP receptor antagonist"},
    {"name": "rimegepant", "class": "CGRP_antagonist", "domain": "headache", "mechanism": "oral CGRP receptor antagonist"},
    {"name": "topiramate", "class": "anticonvulsant", "domain": "headache", "mechanism": "multiple mechanisms"},
    # Neuromuscular
    {"name": "nusinersen", "class": "antisense", "domain": "neuromuscular", "mechanism": "SMN2 splicing modifier"},
    {"name": "risdiplam", "class": "SMN2_modifier", "domain": "neuromuscular", "mechanism": "oral SMN2 splicing modifier"},
    {"name": "onasemnogene", "class": "gene_therapy", "domain": "neuromuscular", "mechanism": "AAV9-SMN1 gene replacement"},
    {"name": "efgartigimod", "class": "FcRn_inhibitor", "domain": "neuromuscular", "mechanism": "neonatal Fc receptor inhibitor"},
    {"name": "tofersen", "class": "antisense", "domain": "neuromuscular", "mechanism": "SOD1 antisense oligonucleotide"},
    # Neuro-Oncology
    {"name": "temozolomide", "class": "alkylating", "domain": "neuro_oncology", "mechanism": "DNA alkylation"},
    {"name": "vorasidenib", "class": "IDH_inhibitor", "domain": "neuro_oncology", "mechanism": "dual IDH1/2 inhibitor"},
]


# ===================================================================
# NEUROLOGY GENES
# ===================================================================

NEURO_GENES: List[Dict[str, Any]] = [
    {"gene": "APP", "domain": "degenerative", "condition": "Early-onset Alzheimer's disease"},
    {"gene": "PSEN1", "domain": "degenerative", "condition": "Early-onset Alzheimer's disease"},
    {"gene": "PSEN2", "domain": "degenerative", "condition": "Early-onset Alzheimer's disease"},
    {"gene": "APOE", "domain": "degenerative", "condition": "Alzheimer's disease risk (e4 allele)"},
    {"gene": "MAPT", "domain": "degenerative", "condition": "Frontotemporal dementia, PSP"},
    {"gene": "GRN", "domain": "degenerative", "condition": "Frontotemporal dementia (progranulin)"},
    {"gene": "C9orf72", "domain": "degenerative", "condition": "FTD-ALS spectrum"},
    {"gene": "SNCA", "domain": "movement", "condition": "Parkinson's disease, Lewy body dementia"},
    {"gene": "LRRK2", "domain": "movement", "condition": "Parkinson's disease (autosomal dominant)"},
    {"gene": "GBA1", "domain": "movement", "condition": "Parkinson's disease risk, Gaucher carrier"},
    {"gene": "PARK2", "domain": "movement", "condition": "Early-onset Parkinson's disease"},
    {"gene": "PINK1", "domain": "movement", "condition": "Early-onset Parkinson's disease"},
    {"gene": "HTT", "domain": "movement", "condition": "Huntington disease (CAG repeat)"},
    {"gene": "ATP7B", "domain": "movement", "condition": "Wilson disease"},
    {"gene": "SCN1A", "domain": "epilepsy", "condition": "Dravet syndrome, GEFS+"},
    {"gene": "SCN2A", "domain": "epilepsy", "condition": "Epileptic encephalopathy"},
    {"gene": "KCNQ2", "domain": "epilepsy", "condition": "Benign familial neonatal seizures"},
    {"gene": "CDKL5", "domain": "epilepsy", "condition": "CDKL5 deficiency disorder"},
    {"gene": "TSC1", "domain": "epilepsy", "condition": "Tuberous sclerosis complex"},
    {"gene": "TSC2", "domain": "epilepsy", "condition": "Tuberous sclerosis complex"},
    {"gene": "DEPDC5", "domain": "epilepsy", "condition": "Familial focal epilepsy"},
    {"gene": "SLC6A1", "domain": "epilepsy", "condition": "Myoclonic-atonic epilepsy"},
    {"gene": "SMN1", "domain": "neuromuscular", "condition": "Spinal muscular atrophy"},
    {"gene": "DMD", "domain": "neuromuscular", "condition": "Duchenne/Becker muscular dystrophy"},
    {"gene": "SOD1", "domain": "neuromuscular", "condition": "Amyotrophic lateral sclerosis (familial)"},
    {"gene": "FUS", "domain": "neuromuscular", "condition": "ALS (familial)"},
    {"gene": "PMP22", "domain": "neuromuscular", "condition": "Charcot-Marie-Tooth type 1A"},
    {"gene": "NF1", "domain": "neuro_oncology", "condition": "Neurofibromatosis type 1"},
    {"gene": "NF2", "domain": "neuro_oncology", "condition": "Neurofibromatosis type 2"},
    {"gene": "IDH1", "domain": "neuro_oncology", "condition": "Glioma molecular subtyping"},
    {"gene": "IDH2", "domain": "neuro_oncology", "condition": "Glioma molecular subtyping"},
    {"gene": "MGMT", "domain": "neuro_oncology", "condition": "Glioblastoma (methylation predictor)"},
    {"gene": "TERT", "domain": "neuro_oncology", "condition": "Glioma prognostic marker"},
    {"gene": "BRAF", "domain": "neuro_oncology", "condition": "Pilocytic astrocytoma, pleomorphic xanthoastrocytoma"},
    {"gene": "AQP4", "domain": "ms", "condition": "Neuromyelitis optica spectrum disorder"},
    {"gene": "MOG", "domain": "ms", "condition": "MOG antibody-associated disease"},
    {"gene": "HLA-DRB1", "domain": "ms", "condition": "Multiple sclerosis susceptibility"},
    {"gene": "CACNA1A", "domain": "headache", "condition": "Familial hemiplegic migraine type 1"},
]


# ===================================================================
# CLINICAL SCALES
# ===================================================================

CLINICAL_SCALES: Dict[str, Dict[str, Any]] = {
    "nihss": {
        "full_name": "NIH Stroke Scale",
        "domain": "cerebrovascular",
        "range": "0-42",
        "max_score": 42,
        "interpretation": {
            "0": "No stroke symptoms",
            "1-4": "Minor stroke",
            "5-15": "Moderate stroke",
            "16-20": "Moderate to severe stroke",
            "21-42": "Severe stroke",
        },
    },
    "gcs": {
        "full_name": "Glasgow Coma Scale",
        "domain": "general",
        "range": "3-15",
        "max_score": 15,
        "interpretation": {
            "3-8": "Severe (coma)",
            "9-12": "Moderate",
            "13-15": "Mild",
        },
    },
    "moca": {
        "full_name": "Montreal Cognitive Assessment",
        "domain": "degenerative",
        "range": "0-30",
        "max_score": 30,
        "interpretation": {
            "26-30": "Normal cognition",
            "18-25": "Mild cognitive impairment",
            "10-17": "Moderate cognitive impairment",
            "0-9": "Severe cognitive impairment",
        },
    },
    "updrs": {
        "full_name": "MDS-UPDRS Part III (Motor Examination)",
        "domain": "movement",
        "range": "0-132",
        "max_score": 132,
        "interpretation": {
            "0-10": "Minimal motor findings",
            "11-22": "Mild motor impairment",
            "23-40": "Moderate motor impairment",
            "41-132": "Severe motor impairment",
        },
    },
    "edss": {
        "full_name": "Expanded Disability Status Scale",
        "domain": "ms",
        "range": "0-10",
        "max_score": 10.0,
        "interpretation": {
            "0": "Normal neurological exam",
            "1.0-3.5": "Ambulatory, mild disability",
            "4.0-5.5": "Walking limited, moderate disability",
            "6.0-6.5": "Requires walking aid",
            "7.0-9.5": "Wheelchair / bed-bound",
            "10.0": "Death due to MS",
        },
    },
    "mrs": {
        "full_name": "Modified Rankin Scale",
        "domain": "cerebrovascular",
        "range": "0-6",
        "max_score": 6,
        "interpretation": {
            "0": "No symptoms",
            "1": "No significant disability",
            "2": "Slight disability",
            "3": "Moderate disability (walks without assistance)",
            "4": "Moderately severe disability",
            "5": "Severe disability (bedridden)",
            "6": "Dead",
        },
    },
    "hit6": {
        "full_name": "Headache Impact Test-6",
        "domain": "headache",
        "range": "36-78",
        "max_score": 78,
        "interpretation": {
            "36-49": "Little/no impact",
            "50-55": "Some impact",
            "56-59": "Substantial impact",
            "60-78": "Severe impact",
        },
    },
    "alsfrs": {
        "full_name": "ALS Functional Rating Scale - Revised",
        "domain": "neuromuscular",
        "range": "0-48",
        "max_score": 48,
        "interpretation": {
            "39-48": "Mild functional impairment",
            "25-38": "Moderate functional impairment",
            "13-24": "Severe functional impairment",
            "0-12": "Very severe functional impairment",
        },
    },
    "aspects": {
        "full_name": "Alberta Stroke Program Early CT Score",
        "domain": "cerebrovascular",
        "range": "0-10",
        "max_score": 10,
        "interpretation": {
            "8-10": "Small infarct core (favorable for thrombectomy)",
            "6-7": "Moderate infarct core",
            "0-5": "Large infarct core (unfavorable for thrombectomy)",
        },
    },
    "hoehn_yahr": {
        "full_name": "Hoehn and Yahr Scale",
        "domain": "movement",
        "range": "0-5",
        "max_score": 5,
        "interpretation": {
            "0": "No signs of disease",
            "1": "Unilateral involvement only",
            "2": "Bilateral involvement without impairment of balance",
            "3": "Mild to moderate bilateral disease; physically independent",
            "4": "Severe disability; still able to walk or stand unassisted",
            "5": "Wheelchair bound or bedridden unless aided",
        },
    },
}


# ===================================================================
# PEDIATRIC NEURO-ONCOLOGY
# ===================================================================

PEDIATRIC_NEURO_ONCOLOGY: Dict[str, Dict[str, Any]] = {
    "methotrexate_leukoencephalopathy": {
        "description": "Methotrexate-induced leukoencephalopathy in pediatric ALL",
        "incidence": "3-10% of pediatric ALL patients on high-dose MTX",
        "risk_factors": [
            "Age <5 years",
            "MTHFR C677T polymorphism",
            "Cumulative intrathecal MTX >12 doses",
        ],
        "mri_findings": "Periventricular white matter T2 hyperintensities",
        "prevention": [
            "Leucovorin rescue",
            "Aggressive hydration",
            "Urinary alkalinization",
        ],
    },
    "vincristine_neuropathy": {
        "description": "Vincristine-induced peripheral neuropathy in children",
        "incidence": "30-40% of children; dose-limiting toxicity",
        "manifestations": [
            "Foot drop",
            "Constipation",
            "Jaw pain",
        ],
        "dose_cap": "2 mg per dose (absolute cap regardless of BSA)",
        "recovery": "Usually within 3 months of discontinuation",
    },
    "PRES": {
        "description": "Posterior reversible encephalopathy syndrome during induction chemotherapy",
        "incidence": "1-5% during induction",
        "presentation": "Seizures, headache, visual disturbance, altered mental status",
        "imaging": "Bilateral posterior white matter vasogenic edema on MRI",
        "management": "Blood pressure control, seizure management, hold offending agent",
    },
    "asparaginase_cerebral_thrombosis": {
        "description": "L-asparaginase-associated cerebral thrombosis",
        "incidence": "2-3% risk",
        "typical_location": "Sagittal sinus thrombosis",
        "mechanism": "Decreased antithrombin III, protein C, and protein S synthesis",
        "management": "Anticoagulation; resume asparaginase on case-by-case basis",
    },
    "cranial_radiation_effects": {
        "description": "Long-term neurotoxicity from cranial radiation in children",
        "effects": [
            "Neurocognitive decline (IQ loss ~2-4 points per Gy)",
            "Secondary meningiomas (latency 20+ years)",
            "Endocrinopathies (GH deficiency, hypothyroidism)",
            "White matter changes",
            "Cerebrovascular disease (moyamoya-like)",
        ],
    },
    "pediatric_brain_tumors": {
        "description": "Common pediatric brain tumors with neurological significance",
        "key_types": {
            "medulloblastoma": {
                "location": "Posterior fossa",
                "molecular_groups": ["WNT", "SHH", "Group 3", "Group 4"],
                "prognosis": "WNT best (>90% OS), Group 3 worst (~50% OS)",
            },
            "diffuse_midline_glioma": {
                "mutation": "H3K27M",
                "location": "Thalamus, brainstem (DIPG), spinal cord",
                "prognosis": "Median survival 9-11 months",
            },
        },
    },
}


# ===================================================================
# CONVENIENCE ACCESSORS
# ===================================================================


def get_domain_count() -> int:
    """Return the number of neurological disease domains."""
    return len(NEURO_DOMAINS)


def get_drug_count() -> int:
    """Return the number of curated neurology drugs."""
    return len(NEURO_DRUGS)


def get_gene_count() -> int:
    """Return the number of curated neurology genes."""
    return len(NEURO_GENES)


def get_scale_count() -> int:
    """Return the number of clinical scales."""
    return len(CLINICAL_SCALES)


def get_drugs_by_domain(domain: str) -> List[Dict[str, Any]]:
    """Return drugs filtered by neurological domain."""
    return [d for d in NEURO_DRUGS if d["domain"] == domain]


def get_genes_by_domain(domain: str) -> List[Dict[str, Any]]:
    """Return genes filtered by neurological domain."""
    return [g for g in NEURO_GENES if g["domain"] == domain]
