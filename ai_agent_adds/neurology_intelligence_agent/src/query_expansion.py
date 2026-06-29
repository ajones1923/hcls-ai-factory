"""Neurology Intelligence Agent -- Query Expansion Module.

Provides entity alias resolution, synonym mapping, and query expansion
for neurology-specific clinical queries.  Ensures that abbreviations,
brand names, and colloquial terms are normalized and expanded to improve
vector-search recall in the RAG pipeline.

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import logging
import re
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)


# ===================================================================
# ENTITY ALIASES (150+ neuro-specific aliases)
# ===================================================================

ENTITY_ALIASES: Dict[str, str] = {
    # -- General abbreviations --
    "TIA": "transient ischemic attack",
    "SAH": "subarachnoid hemorrhage",
    "ICH": "intracerebral hemorrhage",
    "ALS": "amyotrophic lateral sclerosis",
    "MS": "multiple sclerosis",
    "PD": "Parkinson's disease",
    "AD": "Alzheimer's disease",
    "GBM": "glioblastoma multiforme",
    "DBS": "deep brain stimulation",
    "tPA": "tissue plasminogen activator",
    "TPA": "tissue plasminogen activator",
    "ASM": "anti-seizure medication",
    "AED": "anti-epileptic drug",
    "ASMs": "anti-seizure medications",
    "AEDs": "anti-epileptic drugs",
    "DMT": "disease-modifying therapy",
    "DMTs": "disease-modifying therapies",
    "CGRP": "calcitonin gene-related peptide",
    "NfL": "neurofilament light chain",
    "NFL": "neurofilament light chain",
    "CSF": "cerebrospinal fluid",
    "EEG": "electroencephalogram",
    "EMG": "electromyography",
    "NCS": "nerve conduction studies",
    "LP": "lumbar puncture",
    "MRI": "magnetic resonance imaging",
    "CT": "computed tomography",
    "CTA": "CT angiography",
    "MRA": "MR angiography",
    "CTP": "CT perfusion",
    "DWI": "diffusion-weighted imaging",
    "FLAIR": "fluid-attenuated inversion recovery",
    "SWI": "susceptibility-weighted imaging",
    "DTI": "diffusion tensor imaging",
    "fMRI": "functional magnetic resonance imaging",
    "PET": "positron emission tomography",
    "SPECT": "single-photon emission computed tomography",
    "DAT": "dopamine transporter",
    "DAT scan": "dopamine transporter scan",
    "MIBG": "metaiodobenzylguanidine",

    # -- Clinical scales --
    "NIHSS": "NIH Stroke Scale",
    "GCS": "Glasgow Coma Scale",
    "MMSE": "Mini-Mental State Examination",
    "MoCA": "Montreal Cognitive Assessment",
    "CDR": "Clinical Dementia Rating",
    "EDSS": "Expanded Disability Status Scale",
    "UPDRS": "Unified Parkinson's Disease Rating Scale",
    "MDS-UPDRS": "Movement Disorder Society Unified Parkinson's Disease Rating Scale",
    "ALSFRS-R": "ALS Functional Rating Scale Revised",
    "ALSFRS": "ALS Functional Rating Scale Revised",
    "HIT-6": "Headache Impact Test 6",
    "MIDAS": "Migraine Disability Assessment",
    "mRS": "modified Rankin Scale",
    "ASPECTS": "Alberta Stroke Program Early CT Score",
    "KPS": "Karnofsky Performance Status",
    "RANO": "Response Assessment in Neuro-Oncology",
    "ABCD2": "ABCD2 stroke risk score",
    "H&Y": "Hoehn and Yahr scale",
    "HY": "Hoehn and Yahr scale",

    # -- Diseases / syndromes --
    "RRMS": "relapsing-remitting multiple sclerosis",
    "SPMS": "secondary progressive multiple sclerosis",
    "PPMS": "primary progressive multiple sclerosis",
    "CIS": "clinically isolated syndrome",
    "NMOSD": "neuromyelitis optica spectrum disorder",
    "NMO": "neuromyelitis optica",
    "MOGAD": "MOG antibody-associated disease",
    "ADEM": "acute disseminated encephalomyelitis",
    "bvFTD": "behavioral variant frontotemporal dementia",
    "FTD": "frontotemporal dementia",
    "PPA": "primary progressive aphasia",
    "svPPA": "semantic variant primary progressive aphasia",
    "nfvPPA": "nonfluent variant primary progressive aphasia",
    "lvPPA": "logopenic variant primary progressive aphasia",
    "DLB": "dementia with Lewy bodies",
    "LBD": "Lewy body dementia",
    "PSP": "progressive supranuclear palsy",
    "CBD": "corticobasal degeneration",
    "CBS": "corticobasal syndrome",
    "MSA": "multiple system atrophy",
    "MSA-C": "multiple system atrophy cerebellar type",
    "MSA-P": "multiple system atrophy parkinsonian type",
    "HD": "Huntington disease",
    "CJD": "Creutzfeldt-Jakob disease",
    "vCJD": "variant Creutzfeldt-Jakob disease",
    "NPH": "normal pressure hydrocephalus",
    "PCA": "posterior cortical atrophy",
    "MG": "myasthenia gravis",
    "GBS": "Guillain-Barre syndrome",
    "CIDP": "chronic inflammatory demyelinating polyneuropathy",
    "SMA": "spinal muscular atrophy",
    "DMD": "Duchenne muscular dystrophy",
    "BMD": "Becker muscular dystrophy",
    "CMT": "Charcot-Marie-Tooth disease",
    "IBM": "inclusion body myositis",
    "LEMS": "Lambert-Eaton myasthenic syndrome",
    "SE": "status epilepticus",
    "NCSE": "non-convulsive status epilepticus",
    "CSE": "convulsive status epilepticus",
    "JME": "juvenile myoclonic epilepsy",
    "CAE": "childhood absence epilepsy",
    "TLE": "temporal lobe epilepsy",
    "FLE": "frontal lobe epilepsy",
    "MTS": "mesial temporal sclerosis",
    "HS": "hippocampal sclerosis",
    "FCD": "focal cortical dysplasia",
    "TSC": "tuberous sclerosis complex",
    "LGS": "Lennox-Gastaut syndrome",
    "BECTS": "benign epilepsy with centrotemporal spikes",
    "PNES": "psychogenic non-epileptic seizures",
    "GTCS": "generalized tonic-clonic seizure",
    "GTC": "generalized tonic-clonic",
    "IIH": "idiopathic intracranial hypertension",
    "PTC": "pseudotumor cerebri",
    "MOH": "medication overuse headache",
    "NDPH": "new daily persistent headache",
    "TTH": "tension-type headache",
    "SUNCT": "short-lasting unilateral neuralgiform headache attacks with conjunctival injection and tearing",
    "SUNA": "short-lasting unilateral neuralgiform headache attacks with cranial autonomic features",
    "TAC": "trigeminal autonomic cephalalgia",
    "TN": "trigeminal neuralgia",
    "ON": "optic neuritis",
    "RBD": "REM sleep behavior disorder",
    "OSA": "obstructive sleep apnea",
    "CSA": "central sleep apnea",
    "RLS": "restless legs syndrome",
    "WED": "Willis-Ekbom disease",
    "PLMD": "periodic limb movement disorder",
    "PLMS": "periodic limb movements of sleep",
    "MSLT": "multiple sleep latency test",
    "PSG": "polysomnography",
    "ESS": "Epworth Sleepiness Scale",

    # -- Autoimmune / neuroimmunology --
    "NMDAR": "N-methyl-D-aspartate receptor",
    "NMDA": "N-methyl-D-aspartate",
    "LGI1": "leucine-rich glioma-inactivated 1",
    "CASPR2": "contactin-associated protein-like 2",
    "AQP4": "aquaporin-4",
    "MOG": "myelin oligodendrocyte glycoprotein",
    "GFAP": "glial fibrillary acidic protein",
    "AChR": "acetylcholine receptor",
    "MuSK": "muscle-specific kinase",
    "GAD65": "glutamic acid decarboxylase 65",
    "VGCC": "voltage-gated calcium channel",
    "VGKC": "voltage-gated potassium channel",
    "JCV": "JC virus",
    "PML": "progressive multifocal leukoencephalopathy",
    "SREAT": "steroid-responsive encephalopathy associated with autoimmune thyroiditis",
    "SPS": "stiff-person syndrome",

    # -- Neuro-oncology --
    "IDH": "isocitrate dehydrogenase",
    "MGMT": "O6-methylguanine-DNA methyltransferase",
    "1p19q": "1p/19q codeletion",
    "EGFR": "epidermal growth factor receptor",
    "PCNSL": "primary CNS lymphoma",
    "GBM": "glioblastoma",  # noqa: F601
    "WHO": "World Health Organization",
    "TMZ": "temozolomide",
    "SRS": "stereotactic radiosurgery",
    "WBRT": "whole brain radiation therapy",
    "TTFields": "tumor treating fields",
    "BBB": "blood-brain barrier",

    # -- Drug brand names -> generic --
    "Leqembi": "lecanemab",
    "Kisunla": "donanemab",
    "Aduhelm": "aducanumab",
    "Sinemet": "levodopa/carbidopa",
    "Madopar": "levodopa/benserazide",
    "Mirapex": "pramipexole",
    "Requip": "ropinirole",
    "Azilect": "rasagiline",
    "Xadago": "safinamide",
    "Symmetrel": "amantadine",
    "Gocovri": "amantadine extended-release",
    "Inbrija": "levodopa inhalation",
    "Duopa": "carbidopa/levodopa enteral suspension",
    "Aimovig": "erenumab",
    "Emgality": "galcanezumab",
    "Ajovy": "fremanezumab",
    "Qulipta": "atogepant",
    "Nurtec": "rimegepant",
    "Ubrelvy": "ubrogepant",
    "Reyvow": "lasmiditan",
    "Botox": "onabotulinumtoxinA",
    "Ocrevus": "ocrelizumab",
    "Kesimpta": "ofatumumab",
    "Tysabri": "natalizumab",
    "Gilenya": "fingolimod",
    "Mayzent": "siponimod",
    "Zeposia": "ozanimod",
    "Ponvory": "ponesimod",
    "Tecfidera": "dimethyl fumarate",
    "Vumerity": "diroximel fumarate",
    "Copaxone": "glatiramer acetate",
    "Aubagio": "teriflunomide",
    "Mavenclad": "cladribine",
    "Lemtrada": "alemtuzumab",
    "Avonex": "interferon beta-1a intramuscular",
    "Rebif": "interferon beta-1a subcutaneous",
    "Betaseron": "interferon beta-1b",
    "Plegridy": "peginterferon beta-1a",
    "Spinraza": "nusinersen",
    "Evrysdi": "risdiplam",
    "Zolgensma": "onasemnogene abeparvovec",
    "Vyvgart": "efgartigimod",
    "Soliris": "eculizumab",
    "Enspryng": "satralizumab",
    "Uplizna": "inebilizumab",
    "Qalsody": "tofersen",
    "Riluzole": "riluzole",
    "Rilutek": "riluzole",
    "Radicava": "edaravone",
    "Relyvrio": "sodium phenylbutyrate/taurursodiol",
    "Temodar": "temozolomide",
    "Voranigo": "vorasidenib",
    "Epidiolex": "cannabidiol",
    "Fintepla": "fenfluramine",
    "Diacomit": "stiripentol",
    "Banzel": "rufinamide",
    "Xcopri": "cenobamate",
    "Fycompa": "perampanel",
    "Briviact": "brivaracetam",
    "Keppra": "levetiracetam",
    "Lamictal": "lamotrigine",
    "Tegretol": "carbamazepine",
    "Trileptal": "oxcarbazepine",
    "Depakote": "valproate/divalproex",
    "Dilantin": "phenytoin",
    "Vimpat": "lacosamide",
    "Topamax": "topiramate",
    "Zonegran": "zonisamide",
    "Sabril": "vigabatrin",
    "Onfi": "clobazam",
    "Gabapentin": "gabapentin",
    "Neurontin": "gabapentin",
    "Lyrica": "pregabalin",
    "Ingrezza": "valbenazine",
    "Austedo": "deutetrabenazine",
    "Xenazine": "tetrabenazine",
    "Activase": "alteplase",
    "TNKase": "tenecteplase",
    "Plavix": "clopidogrel",
    "Eliquis": "apixaban",
    "Xarelto": "rivaroxaban",
    "Aricept": "donepezil",
    "Exelon": "rivastigmine",
    "Razadyne": "galantamine",
    "Namenda": "memantine",
    "Namzaric": "memantine/donepezil",
    "Zafgen": "ganaxolone",
    "Ztalmy": "ganaxolone",
    "Imitrex": "sumatriptan",
    "Maxalt": "rizatriptan",
    "Relpax": "eletriptan",
    "Zomig": "zolmitriptan",
    "Treximet": "sumatriptan/naproxen",
}


# ===================================================================
# SYNONYM MAPS (16 domain-specific synonym maps)
# ===================================================================

STROKE_MAP: Dict[str, List[str]] = {
    "stroke": ["cerebrovascular accident", "CVA", "brain attack", "cerebral infarction",
               "ischemic stroke", "hemorrhagic stroke", "brain infarct"],
    "tpa": ["alteplase", "tissue plasminogen activator", "thrombolytic", "clot buster",
            "IV tPA", "intravenous thrombolysis"],
    "thrombectomy": ["mechanical thrombectomy", "endovascular thrombectomy",
                     "clot retrieval", "stent retriever", "aspiration thrombectomy",
                     "neurointerventional procedure"],
    "lvo": ["large vessel occlusion", "proximal occlusion", "M1 occlusion",
            "ICA occlusion", "basilar occlusion"],
    "tia": ["transient ischemic attack", "mini-stroke", "warning stroke",
            "transient neurological deficit"],
    "hemorrhage": ["intracerebral hemorrhage", "ICH", "brain bleed",
                   "hypertensive hemorrhage", "cerebral hemorrhage"],
    "sah": ["subarachnoid hemorrhage", "aneurysmal SAH", "ruptured aneurysm",
            "berry aneurysm rupture"],
    "carotid": ["carotid stenosis", "carotid artery disease", "carotid endarterectomy",
                "CEA", "carotid stenting", "CAS"],
    "afib": ["atrial fibrillation", "AF", "cardioembolic", "cardiac embolism"],
}

DEMENTIA_MAP: Dict[str, List[str]] = {
    "alzheimers": ["Alzheimer's disease", "AD", "senile dementia", "amyloid plaques",
                   "neurofibrillary tangles", "tau pathology", "amyloid cascade"],
    "frontotemporal": ["frontotemporal dementia", "FTD", "Pick disease",
                       "frontotemporal lobar degeneration", "FTLD", "behavioral variant"],
    "lewy_body": ["Lewy body dementia", "DLB", "dementia with Lewy bodies",
                  "Lewy body disease", "alpha-synuclein dementia"],
    "vascular": ["vascular dementia", "vascular cognitive impairment", "VCI",
                 "multi-infarct dementia", "subcortical vascular dementia",
                 "Binswanger disease"],
    "mci": ["mild cognitive impairment", "MCI", "pre-dementia", "subjective cognitive decline",
            "SCD", "amnestic MCI"],
    "biomarkers": ["amyloid PET", "tau PET", "CSF biomarkers", "p-tau",
                   "amyloid-beta 42", "ATN framework", "NIA-AA criteria"],
    "anti_amyloid": ["lecanemab", "donanemab", "aducanumab", "amyloid immunotherapy",
                     "anti-amyloid antibody", "amyloid clearance"],
}

EPILEPSY_MAP: Dict[str, List[str]] = {
    "seizure": ["convulsion", "fit", "epileptic seizure", "ictus", "paroxysmal event",
                "epileptic event"],
    "focal": ["partial seizure", "focal onset", "focal aware", "focal impaired awareness",
              "complex partial", "simple partial"],
    "generalized": ["generalized seizure", "primary generalized", "absence seizure",
                    "myoclonic seizure", "tonic-clonic", "grand mal", "petit mal"],
    "status": ["status epilepticus", "SE", "prolonged seizure", "refractory SE",
               "super-refractory SE", "NCSE", "non-convulsive status"],
    "drug_resistant": ["refractory epilepsy", "drug-resistant epilepsy", "intractable epilepsy",
                       "medically refractory", "pharmacoresistant"],
    "surgery": ["epilepsy surgery", "temporal lobectomy", "amygdalohippocampectomy",
                "VNS", "vagus nerve stimulator", "RNS", "responsive neurostimulation",
                "corpus callosotomy", "hemispherectomy", "laser ablation", "LITT"],
    "dravet": ["Dravet syndrome", "SCN1A", "severe myoclonic epilepsy of infancy", "SMEI"],
}

MS_MAP: Dict[str, List[str]] = {
    "ms": ["multiple sclerosis", "demyelinating disease", "CNS demyelination",
           "white matter disease", "inflammatory demyelination"],
    "relapse": ["MS relapse", "exacerbation", "flare", "attack", "acute relapse",
                "pseudorelapse"],
    "lesion": ["demyelinating lesion", "MS plaque", "white matter lesion",
               "periventricular lesion", "juxtacortical lesion", "corpus callosum lesion",
               "Dawson fingers", "enhancing lesion", "black hole"],
    "nmo": ["neuromyelitis optica", "NMOSD", "Devic disease", "AQP4 antibody disease",
            "aquaporin-4"],
    "mog": ["MOG antibody disease", "MOGAD", "MOG-IgG", "myelin oligodendrocyte glycoprotein"],
    "ocb": ["oligoclonal bands", "OCB", "CSF oligoclonal bands", "intrathecal IgG synthesis"],
    "progression": ["disability progression", "PIRA", "progression independent of relapse activity",
                    "worsening", "secondary progression", "SPMS conversion"],
}

PARKINSONS_MAP: Dict[str, List[str]] = {
    "parkinsons": ["Parkinson's disease", "PD", "parkinsonism", "idiopathic PD",
                   "primary parkinsonism"],
    "tremor": ["resting tremor", "pill-rolling tremor", "action tremor",
               "postural tremor", "essential tremor", "ET"],
    "motor": ["bradykinesia", "rigidity", "akinesia", "dyskinesia", "levodopa-induced dyskinesia",
              "motor fluctuations", "wearing off", "on-off phenomenon"],
    "non_motor": ["REM sleep behavior disorder", "hyposmia", "anosmia", "constipation",
                  "depression", "anxiety", "cognitive impairment", "orthostatic hypotension",
                  "impulse control disorder"],
    "surgical": ["deep brain stimulation", "DBS", "subthalamic nucleus", "STN",
                 "globus pallidus interna", "GPi", "focused ultrasound", "FUS thalamotomy"],
    "genetics": ["LRRK2", "GBA1", "SNCA", "PARK2", "PINK1", "genetic Parkinson's"],
}

BRAIN_TUMOR_MAP: Dict[str, List[str]] = {
    "glioma": ["glioblastoma", "GBM", "astrocytoma", "oligodendroglioma",
               "diffuse glioma", "high-grade glioma", "low-grade glioma",
               "IDH-mutant", "IDH-wildtype"],
    "meningioma": ["meningeal tumor", "convexity meningioma", "skull base meningioma",
                   "atypical meningioma", "anaplastic meningioma"],
    "metastasis": ["brain metastasis", "brain met", "cerebral metastasis",
                   "leptomeningeal carcinomatosis", "leptomeningeal disease"],
    "molecular": ["IDH mutation", "MGMT methylation", "1p19q codeletion",
                  "TERT promoter mutation", "EGFR amplification", "H3K27M",
                  "BRAF V600E", "ATRX loss"],
    "treatment": ["temozolomide", "TMZ", "Stupp protocol", "radiation therapy",
                  "stereotactic radiosurgery", "Gamma Knife", "CyberKnife",
                  "tumor treating fields", "bevacizumab", "Avastin"],
    "pcnsl": ["primary CNS lymphoma", "PCNSL", "CNS lymphoma",
              "high-dose methotrexate", "intrathecal chemotherapy"],
}

HEADACHE_MAP: Dict[str, List[str]] = {
    "migraine": ["migraine headache", "hemicranial headache", "vascular headache",
                 "migraine with aura", "migraine without aura", "menstrual migraine",
                 "vestibular migraine", "chronic migraine", "episodic migraine"],
    "cluster": ["cluster headache", "suicide headache", "trigeminal autonomic cephalalgia",
                "TAC", "horton headache"],
    "tension": ["tension headache", "tension-type headache", "TTH", "muscle contraction headache",
                "stress headache"],
    "cgrp": ["CGRP", "calcitonin gene-related peptide", "CGRP inhibitor",
             "CGRP antibody", "CGRP receptor antagonist", "gepant", "anti-CGRP"],
    "overuse": ["medication overuse headache", "MOH", "rebound headache",
                "analgesic overuse", "triptan overuse"],
    "red_flags": ["thunderclap headache", "worst headache of life", "papilledema",
                  "new headache after 50", "progressive headache", "headache with fever",
                  "SNOOP criteria"],
    "trigeminal": ["trigeminal neuralgia", "tic douloureux", "facial pain",
                   "V2/V3 pain", "electric shock face pain"],
}

NEUROMUSCULAR_MAP: Dict[str, List[str]] = {
    "als": ["amyotrophic lateral sclerosis", "ALS", "Lou Gehrig's disease",
            "motor neuron disease", "MND", "progressive muscular atrophy"],
    "myasthenia": ["myasthenia gravis", "MG", "neuromuscular junction disorder",
                   "AChR antibody", "MuSK antibody", "myasthenic crisis",
                   "thymectomy", "thymus"],
    "gbs": ["Guillain-Barre syndrome", "GBS", "acute inflammatory demyelinating polyneuropathy",
            "AIDP", "Miller Fisher syndrome", "axonal GBS", "AMAN", "AMSAN"],
    "cidp": ["chronic inflammatory demyelinating polyneuropathy", "CIDP",
             "chronic demyelinating neuropathy", "IVIg", "plasma exchange"],
    "sma": ["spinal muscular atrophy", "SMA", "SMN1", "SMN2",
            "nusinersen", "risdiplam", "onasemnogene", "gene therapy SMA"],
    "neuropathy": ["peripheral neuropathy", "polyneuropathy", "mononeuropathy",
                   "diabetic neuropathy", "small fiber neuropathy",
                   "large fiber neuropathy", "autonomic neuropathy"],
    "myopathy": ["inflammatory myopathy", "dermatomyositis", "polymyositis",
                 "inclusion body myositis", "IBM", "necrotizing autoimmune myopathy",
                 "muscular dystrophy"],
}

EEG_MAP: Dict[str, List[str]] = {
    "normal": ["normal background", "posterior dominant rhythm", "PDR",
               "alpha rhythm", "symmetric background"],
    "epileptiform": ["spike", "sharp wave", "spike-and-wave", "polyspike",
                     "epileptiform discharge", "IED", "interictal epileptiform discharge"],
    "focal": ["focal slowing", "temporal slowing", "TIRDA", "focal spikes",
              "PLED", "lateralized periodic discharges", "LPD"],
    "generalized": ["generalized slowing", "diffuse slowing", "encephalopathy pattern",
                    "GRDA", "generalized rhythmic delta activity",
                    "generalized periodic discharges", "GPD"],
    "status": ["electrographic seizure", "electrographic status epilepticus",
               "NCSE pattern", "periodic discharges", "rhythmic delta activity"],
    "sleep": ["sleep architecture", "sleep spindles", "K-complex",
              "REM sleep", "NREM sleep", "sleep-wake cycle"],
    "pattern": ["burst suppression", "alpha coma", "theta coma",
                "triphasic waves", "FIRDA", "OIRDA", "hypsarrhythmia",
                "3 Hz spike-and-wave"],
}

NEUROIMAGING_MAP: Dict[str, List[str]] = {
    "mri_brain": ["brain MRI", "cranial MRI", "neuroimaging", "structural MRI",
                  "T1-weighted", "T2-weighted", "FLAIR", "contrast-enhanced MRI"],
    "mri_spine": ["spinal MRI", "cervical MRI", "thoracic MRI", "lumbar MRI",
                  "spinal cord imaging", "myelography"],
    "ct_head": ["CT head", "brain CT", "non-contrast CT head", "NCCT",
                "CT angiogram head and neck"],
    "perfusion": ["CT perfusion", "CTP", "MR perfusion", "PWI", "CBF", "CBV",
                  "MTT", "Tmax", "perfusion-diffusion mismatch"],
    "vascular": ["MRA", "CTA", "digital subtraction angiography", "DSA",
                 "catheter angiography", "cerebral angiogram",
                 "carotid duplex", "transcranial Doppler", "TCD"],
    "advanced": ["PET scan", "FDG-PET", "amyloid PET", "tau PET",
                 "DAT scan", "DaTscan", "SPECT", "MR spectroscopy",
                 "functional MRI", "tractography", "DTI"],
    "findings": ["infarct", "hemorrhage", "edema", "mass effect", "midline shift",
                 "hydrocephalus", "atrophy", "white matter hyperintensities",
                 "enhancement", "restricted diffusion"],
}

NEUROGENETICS_MAP: Dict[str, List[str]] = {
    "testing": ["genetic testing", "whole exome sequencing", "WES",
                "whole genome sequencing", "WGS", "gene panel",
                "targeted sequencing", "chromosomal microarray", "CMA"],
    "inheritance": ["autosomal dominant", "autosomal recessive", "X-linked",
                    "mitochondrial inheritance", "de novo mutation",
                    "germline variant", "somatic variant"],
    "variant": ["pathogenic variant", "likely pathogenic", "VUS",
                "variant of uncertain significance", "benign variant",
                "frameshift", "missense", "nonsense", "splice site"],
    "repeat": ["trinucleotide repeat", "CAG repeat", "CTG repeat",
               "CGG repeat", "repeat expansion", "anticipation"],
    "counseling": ["genetic counseling", "predictive testing", "carrier testing",
                   "prenatal testing", "preimplantation genetic diagnosis", "PGD"],
    "gene_therapy": ["gene therapy", "gene replacement", "antisense oligonucleotide",
                     "ASO", "RNA interference", "RNAi", "CRISPR",
                     "AAV vector", "adeno-associated virus"],
}

MOVEMENT_MAP: Dict[str, List[str]] = {
    "tremor": ["resting tremor", "action tremor", "postural tremor",
               "intention tremor", "essential tremor", "physiologic tremor",
               "Holmes tremor", "rubral tremor", "orthostatic tremor"],
    "dystonia": ["dystonia", "cervical dystonia", "torticollis", "blepharospasm",
                 "writer's cramp", "focal dystonia", "generalized dystonia",
                 "DYT1", "tardive dystonia", "oromandibular dystonia"],
    "chorea": ["chorea", "Huntington chorea", "Sydenham chorea",
               "chorea gravidarum", "hemiballismus", "ballism",
               "choreiform movements"],
    "ataxia": ["cerebellar ataxia", "spinocerebellar ataxia", "SCA",
               "Friedreich ataxia", "ataxia-telangiectasia",
               "sensory ataxia", "gait ataxia"],
    "tic": ["tic disorder", "Tourette syndrome", "motor tic", "vocal tic",
            "tic suppression", "CBIT", "habit reversal therapy"],
    "myoclonus": ["myoclonus", "cortical myoclonus", "subcortical myoclonus",
                  "action myoclonus", "palatal myoclonus", "opsoclonus-myoclonus"],
}

SLEEP_MAP: Dict[str, List[str]] = {
    "narcolepsy": ["narcolepsy", "excessive daytime sleepiness", "cataplexy",
                   "sleep attacks", "hypocretin deficiency", "orexin",
                   "narcolepsy type 1", "narcolepsy type 2"],
    "rbd": ["REM sleep behavior disorder", "RBD", "dream enactment",
            "REM without atonia", "prodromal synucleinopathy"],
    "apnea": ["sleep apnea", "obstructive sleep apnea", "OSA",
              "central sleep apnea", "CSA", "CPAP", "BiPAP",
              "apnea-hypopnea index", "AHI"],
    "insomnia": ["insomnia", "sleep-onset insomnia", "sleep maintenance insomnia",
                 "chronic insomnia", "CBT-I", "fatal familial insomnia"],
    "parasomnia": ["parasomnia", "sleepwalking", "somnambulism",
                   "sleep terrors", "confusional arousals",
                   "NREM parasomnia", "REM parasomnia"],
    "circadian": ["circadian rhythm disorder", "delayed sleep-wake phase",
                  "advanced sleep-wake phase", "irregular sleep-wake rhythm",
                  "non-24-hour sleep-wake disorder", "shift work disorder",
                  "jet lag"],
    "hypersomnia": ["idiopathic hypersomnia", "Kleine-Levin syndrome",
                    "recurrent hypersomnia", "excessive sleepiness",
                    "hypersomnolence"],
}

NEUROIMMUNOLOGY_MAP: Dict[str, List[str]] = {
    "autoimmune_encephalitis": ["autoimmune encephalitis", "AE",
                                "limbic encephalitis", "anti-NMDAR encephalitis",
                                "anti-LGI1 encephalitis", "anti-CASPR2 encephalitis",
                                "antibody-mediated encephalitis"],
    "nmda": ["anti-NMDA receptor encephalitis", "NMDAR encephalitis",
             "ovarian teratoma encephalitis", "dyskinetic encephalitis"],
    "lgi1": ["LGI1 antibody encephalitis", "anti-LGI1", "faciobrachial dystonic seizures",
             "FBDS", "hyponatremia encephalitis"],
    "caspr2": ["CASPR2 antibody disease", "anti-CASPR2", "Morvan syndrome",
               "neuromyotonia", "Isaac syndrome"],
    "paraneoplastic": ["paraneoplastic syndrome", "paraneoplastic cerebellar degeneration",
                       "anti-Hu", "anti-Yo", "anti-Ri", "anti-CV2",
                       "anti-amphiphysin", "onconeural antibody"],
    "stiff_person": ["stiff-person syndrome", "SPS", "GAD65 antibody",
                     "anti-GAD", "stiff-limb syndrome",
                     "progressive encephalomyelitis with rigidity and myoclonus"],
    "vasculitis": ["CNS vasculitis", "PACNS", "primary angiitis of the CNS",
                   "cerebral vasculitis", "neurosarcoidosis",
                   "Susac syndrome", "neuro-Behcet"],
}

NEUROREHAB_MAP: Dict[str, List[str]] = {
    "stroke_rehab": ["stroke rehabilitation", "post-stroke recovery",
                     "constraint-induced movement therapy", "CIMT",
                     "mirror therapy", "functional electrical stimulation"],
    "physical_therapy": ["physical therapy", "PT", "gait training",
                         "balance training", "vestibular rehabilitation",
                         "neurorehabilitation"],
    "occupational_therapy": ["occupational therapy", "OT", "ADL training",
                             "fine motor rehabilitation", "cognitive rehabilitation"],
    "speech_therapy": ["speech-language therapy", "SLP", "swallowing therapy",
                       "dysphagia management", "aphasia therapy",
                       "speech rehabilitation"],
    "neurostimulation": ["transcranial magnetic stimulation", "TMS",
                         "transcranial direct current stimulation", "tDCS",
                         "repetitive TMS", "rTMS", "theta burst stimulation"],
    "spasticity": ["spasticity management", "botulinum toxin", "baclofen pump",
                   "intrathecal baclofen", "ITB", "tone management"],
}

CSF_MAP: Dict[str, List[str]] = {
    "routine": ["CSF analysis", "lumbar puncture", "spinal tap",
                "opening pressure", "cell count", "protein", "glucose"],
    "infection": ["CSF culture", "meningitis", "encephalitis", "HSV PCR",
                  "bacterial meningitis", "viral meningitis",
                  "fungal meningitis", "cryptococcal antigen"],
    "autoimmune": ["oligoclonal bands", "OCB", "IgG index", "IgG synthesis rate",
                   "cell-based assay", "autoimmune panel"],
    "dementia": ["CSF amyloid-beta 42", "CSF phospho-tau", "CSF total tau",
                 "14-3-3 protein", "RT-QuIC", "ATN biomarkers"],
    "oncology": ["CSF cytology", "flow cytometry", "leptomeningeal disease",
                 "CSF protein elevation", "neoplastic meningitis"],
    "pressure": ["opening pressure", "intracranial pressure", "ICP",
                 "idiopathic intracranial hypertension", "IIH",
                 "normal pressure hydrocephalus", "NPH",
                 "CSF leak", "spontaneous intracranial hypotension"],
}

# Aggregate all synonym maps for unified lookup
NEURO_SYNONYMS: Dict[str, Dict[str, List[str]]] = {
    "stroke": STROKE_MAP,
    "dementia": DEMENTIA_MAP,
    "epilepsy": EPILEPSY_MAP,
    "ms": MS_MAP,
    "parkinsons": PARKINSONS_MAP,
    "brain_tumor": BRAIN_TUMOR_MAP,
    "headache": HEADACHE_MAP,
    "neuromuscular": NEUROMUSCULAR_MAP,
    "eeg": EEG_MAP,
    "neuroimaging": NEUROIMAGING_MAP,
    "neurogenetics": NEUROGENETICS_MAP,
    "movement": MOVEMENT_MAP,
    "sleep": SLEEP_MAP,
    "neuroimmunology": NEUROIMMUNOLOGY_MAP,
    "neurorehab": NEUROREHAB_MAP,
    "csf": CSF_MAP,
}


# ===================================================================
# WORKFLOW TERMS
# ===================================================================

# Lazy import to avoid circular dependency; uses string values from
# NeuroWorkflowType enum defined in src.models.

_WORKFLOW_TERMS: Dict[str, List[str]] = {
    "acute_stroke": [
        "stroke", "cerebrovascular", "tPA", "thrombectomy", "thrombolysis",
        "NIHSS", "ASPECTS", "large vessel occlusion", "LVO", "TIA",
        "hemorrhage", "SAH", "ICH", "carotid", "atrial fibrillation",
        "anticoagulation", "antiplatelet", "penumbra", "infarct core",
        "perfusion mismatch", "door-to-needle", "door-to-groin",
    ],
    "dementia_evaluation": [
        "dementia", "Alzheimer", "cognitive decline", "memory loss",
        "MoCA", "MMSE", "CDR", "amyloid", "tau", "biomarker",
        "frontotemporal", "Lewy body", "vascular dementia",
        "lecanemab", "donanemab", "cholinesterase inhibitor",
        "ATN", "NfL", "neurodegeneration",
    ],
    "epilepsy_focus": [
        "seizure", "epilepsy", "convulsion", "EEG", "epileptiform",
        "focal", "generalized", "status epilepticus", "Dravet",
        "Lennox-Gastaut", "anti-seizure", "ASM", "AED",
        "drug-resistant", "epilepsy surgery", "VNS", "RNS",
        "temporal lobe", "spike-and-wave",
    ],
    "brain_tumor": [
        "glioma", "glioblastoma", "GBM", "astrocytoma", "oligodendroglioma",
        "meningioma", "brain metastasis", "IDH", "MGMT", "1p19q",
        "temozolomide", "radiation", "KPS", "RANO",
        "tumor", "neoplasm", "mass lesion", "WHO grade",
    ],
    "ms_monitoring": [
        "multiple sclerosis", "MS", "demyelination", "relapse",
        "EDSS", "lesion", "gadolinium", "enhancing", "DMT",
        "ocrelizumab", "natalizumab", "fingolimod", "JCV",
        "PML", "oligoclonal bands", "NfL", "NMOSD", "MOG",
        "T2 lesion", "spinal cord lesion",
    ],
    "parkinsons_assessment": [
        "Parkinson", "parkinsonism", "tremor", "bradykinesia",
        "rigidity", "UPDRS", "Hoehn and Yahr", "levodopa",
        "dopamine", "DBS", "dyskinesia", "motor fluctuation",
        "wearing off", "LRRK2", "GBA1", "DAT scan",
        "REM sleep behavior disorder", "hyposmia",
    ],
    "headache_classification": [
        "headache", "migraine", "cluster headache", "tension-type",
        "CGRP", "triptan", "aura", "photophobia", "phonophobia",
        "HIT-6", "MIDAS", "chronic migraine", "medication overuse",
        "trigeminal neuralgia", "ICHD-3", "preventive therapy",
        "thunderclap", "new daily persistent",
    ],
    "neuromuscular_evaluation": [
        "ALS", "amyotrophic lateral sclerosis", "myasthenia gravis",
        "Guillain-Barre", "CIDP", "neuropathy", "EMG", "NCS",
        "nerve conduction", "ALSFRS", "weakness", "fasciculation",
        "SMA", "spinal muscular atrophy", "muscular dystrophy",
        "NMJ", "neuromuscular junction", "CK", "muscle biopsy",
    ],
    "general": [
        "neurology", "neurological", "brain", "spinal cord",
        "nervous system", "cranial nerve", "neuropathology",
        "neurophysiology", "neuroradiology",
    ],
}


# ===================================================================
# QUERY EXPANDER CLASS
# ===================================================================


class QueryExpander:
    """Expand neurology queries with aliases, synonyms, and workflow terms.

    The expander performs three operations:
    1. **Entity detection** -- identifies known abbreviations and brand names
       in the query text and resolves them to canonical forms.
    2. **Synonym expansion** -- augments the query with domain-specific
       synonyms drawn from the 16 synonym maps.
    3. **Workflow term injection** -- adds terms relevant to the detected
       or specified NeuroWorkflowType to improve recall.
    """

    def __init__(
        self,
        aliases: Optional[Dict[str, str]] = None,
        synonyms: Optional[Dict[str, Dict[str, List[str]]]] = None,
        workflow_terms: Optional[Dict[str, List[str]]] = None,
        max_expansion_terms: int = 30,
    ) -> None:
        self._aliases = aliases or ENTITY_ALIASES
        self._synonyms = synonyms or NEURO_SYNONYMS
        self._workflow_terms = workflow_terms or _WORKFLOW_TERMS
        self._max_expansion_terms = max_expansion_terms

        # Pre-compute a case-insensitive lookup for aliases
        self._alias_lower: Dict[str, str] = {
            k.lower(): v for k, v in self._aliases.items()
        }

        # Build a flat reverse index: synonym term -> category
        self._synonym_index: Dict[str, str] = {}
        for domain, category_map in self._synonyms.items():
            for category, terms in category_map.items():
                for term in terms:
                    self._synonym_index[term.lower()] = f"{domain}.{category}"

    # -----------------------------------------------------------------
    # Public API
    # -----------------------------------------------------------------

    def expand(self, query: str, workflow: Optional[str] = None) -> List[str]:
        """Return a list of expansion terms for the given query.

        Parameters
        ----------
        query : str
            Free-text clinical query.
        workflow : str, optional
            NeuroWorkflowType value (e.g. ``"acute_stroke"``).  If ``None``,
            the expander attempts to infer the workflow from query content.

        Returns
        -------
        list[str]
            Deduplicated list of expansion terms, capped at
            ``max_expansion_terms``.
        """
        terms: list[str] = []
        query_lower = query.lower()

        # 1. Resolve entities / aliases
        detected = self.detect_entities(query)
        for _alias, canonical in detected.items():
            terms.append(canonical)

        # 2. Synonym expansion
        for domain, category_map in self._synonyms.items():
            for category, syns in category_map.items():
                # Check if query mentions category key or any synonym
                if category.lower() in query_lower:
                    terms.extend(syns[:5])  # top 5 per match
                    continue
                for syn in syns:
                    if syn.lower() in query_lower:
                        terms.extend(syns[:5])
                        break

        # 3. Workflow-specific terms
        wf_key = workflow or self._infer_workflow(query_lower)
        wf_terms = self.get_workflow_terms(wf_key)
        terms.extend(wf_terms)

        # Deduplicate while preserving order, cap at limit
        seen: set[str] = set()
        unique: list[str] = []
        for t in terms:
            t_lower = t.lower()
            if t_lower not in seen and t_lower not in query_lower:
                seen.add(t_lower)
                unique.append(t)
        return unique[: self._max_expansion_terms]

    def detect_entities(self, query: str) -> Dict[str, str]:
        """Detect abbreviations and brand names in the query.

        Returns
        -------
        dict[str, str]
            Mapping of matched alias -> canonical expansion.
        """
        detected: Dict[str, str] = {}
        # Common English words that should not be treated as aliases even
        # if they happen to match (e.g. "on" -> ON, "or" could match, etc.)
        _STOP_WORDS = {
            "a", "an", "and", "are", "as", "at", "be", "by", "for", "from",
            "had", "has", "have", "he", "her", "his", "how", "i", "if", "in",
            "is", "it", "its", "my", "no", "not", "of", "on", "or", "our",
            "she", "so", "than", "that", "the", "then", "they", "this", "to",
            "up", "us", "was", "we", "what", "when", "who", "will", "with",
            "you",
        }
        # Tokenize into words and multi-word spans
        tokens = re.findall(r"[A-Za-z0-9\-']+", query)
        for token in tokens:
            token_lower = token.lower()
            if token_lower in _STOP_WORDS:
                continue
            # Exact-case match first (case-sensitive abbreviations like MS, PD)
            if token in self._aliases:
                detected[token] = self._aliases[token]
            elif token_lower in self._alias_lower:
                detected[token] = self._alias_lower[token_lower]
        return detected

    def get_workflow_terms(self, workflow: Optional[str] = None) -> List[str]:
        """Return search terms associated with the given workflow.

        Parameters
        ----------
        workflow : str, optional
            NeuroWorkflowType value.  Falls back to ``"general"``
            if not recognized.

        Returns
        -------
        list[str]
            Search terms for the workflow.
        """
        if workflow is None:
            return self._workflow_terms.get("general", [])
        return self._workflow_terms.get(workflow, self._workflow_terms.get("general", []))

    # -----------------------------------------------------------------
    # Internal helpers
    # -----------------------------------------------------------------

    def _infer_workflow(self, query_lower: str) -> str:
        """Best-effort workflow inference from query text.

        Scores each workflow by counting how many of its terms appear
        in the query and picks the highest-scoring one.
        """
        best_wf = "general"
        best_score = 0
        for wf_key, terms in self._workflow_terms.items():
            if wf_key == "general":
                continue
            score = sum(1 for t in terms if t.lower() in query_lower)
            if score > best_score:
                best_score = score
                best_wf = wf_key
        return best_wf if best_score > 0 else "general"
