"""EEG pattern and findings parser for the Neurology Intelligence Agent.

Parses EEG patterns and seeds 45 key EEG patterns and findings covering
normal variants, epileptiform discharges, seizure patterns, encephalopathy
patterns, coma patterns, and specialized monitoring findings.

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List

from .base import BaseIngestParser, IngestRecord

logger = logging.getLogger(__name__)


# ===================================================================
# SEED DATA: 45 EEG PATTERNS AND FINDINGS
# ===================================================================

EEG_PATTERNS: List[Dict[str, Any]] = [
    # --- Epileptiform Discharges ---
    {
        "pattern_id": "EEG-001",
        "name": "3 Hz Generalized Spike-and-Wave",
        "category": "epileptiform",
        "frequency": "3 Hz",
        "morphology": "Regular 3 Hz spike-and-wave complexes, bifrontally predominant, bilaterally synchronous",
        "clinical_correlation": "Childhood absence epilepsy; provoked by hyperventilation",
        "syndrome": "childhood_absence",
        "significance": "Hallmark of typical absence seizures; high yield with hyperventilation provocation",
    },
    {
        "pattern_id": "EEG-002",
        "name": "Polyspike-and-Wave Complexes",
        "category": "epileptiform",
        "frequency": "3-6 Hz",
        "morphology": "Multiple spikes (3-5) preceding slow wave; generalized, often with frontocentral maximum",
        "clinical_correlation": "Juvenile myoclonic epilepsy; photosensitivity common",
        "syndrome": "juvenile_myoclonic",
        "significance": "Strong association with JME; lifelong AED therapy typically required",
    },
    {
        "pattern_id": "EEG-003",
        "name": "Temporal Spike/Sharp Wave",
        "category": "epileptiform",
        "frequency": "intermittent",
        "morphology": "Unilateral anterior temporal sharp waves or spikes; F7/T3 or F8/T4 maximal",
        "clinical_correlation": "Temporal lobe epilepsy; mesial temporal sclerosis",
        "syndrome": "temporal_lobe",
        "significance": "Concordance with MRI and semiology predicts >80% seizure freedom post-surgery",
    },
    {
        "pattern_id": "EEG-004",
        "name": "Centrotemporal Spikes (Rolandic)",
        "category": "epileptiform",
        "frequency": "intermittent, activated by sleep",
        "morphology": "High-amplitude diphasic spikes at C3/C4-T3/T4 with horizontal dipole; tangential orientation",
        "clinical_correlation": "Benign epilepsy with centrotemporal spikes (BECTS); self-limited childhood epilepsy",
        "syndrome": "benign_rolandic",
        "significance": "Self-limited by adolescence; typically does not require treatment if seizures are infrequent",
    },
    {
        "pattern_id": "EEG-005",
        "name": "Hypsarrhythmia",
        "category": "epileptiform",
        "frequency": "chaotic",
        "morphology": "High-amplitude chaotic slow waves intermixed with multifocal spikes; no recognizable background organization",
        "clinical_correlation": "West syndrome (infantile spasms); age 3-12 months",
        "syndrome": "west",
        "significance": "Medical emergency; ACTH or vigabatrin first-line; prognosis depends on etiology",
    },
    {
        "pattern_id": "EEG-006",
        "name": "Slow Spike-and-Wave (<2.5 Hz)",
        "category": "epileptiform",
        "frequency": "1.5-2.5 Hz",
        "morphology": "Slow generalized spike-and-wave <2.5 Hz; diffuse, often with anterior predominance; bursts of fast rhythms",
        "clinical_correlation": "Lennox-Gastaut syndrome; refractory epilepsy",
        "syndrome": "lennox_gastaut",
        "significance": "Treatment-resistant; multiple seizure types; cognitive decline; consider VNS, corpus callosotomy",
    },
    {
        "pattern_id": "EEG-007",
        "name": "Generalized Paroxysmal Fast Activity (GPFA)",
        "category": "epileptiform",
        "frequency": "10-25 Hz",
        "morphology": "Generalized fast activity bursts during NREM sleep; low-medium amplitude",
        "clinical_correlation": "Lennox-Gastaut syndrome; tonic seizures during sleep",
        "syndrome": "lennox_gastaut",
        "significance": "Correlates with tonic seizures; highly specific for LGS; often subclinical",
    },
    {
        "pattern_id": "EEG-008",
        "name": "Focal Onset Seizure Pattern (Temporal)",
        "category": "ictal",
        "frequency": "5-9 Hz",
        "morphology": "Rhythmic temporal theta evolving in frequency, amplitude, and distribution; post-ictal slowing",
        "clinical_correlation": "Mesial temporal lobe seizure; may propagate to contralateral temporal lobe",
        "syndrome": "temporal_lobe",
        "significance": "Onset lateralization critical for surgical planning; propagation pattern informs network",
    },
    # --- Encephalopathy Patterns ---
    {
        "pattern_id": "EEG-009",
        "name": "Generalized Periodic Discharges (GPDs)",
        "category": "encephalopathy",
        "frequency": "0.5-2 Hz",
        "morphology": "Periodic sharp waves or complexes repeating at regular intervals; generalized distribution",
        "clinical_correlation": "CJD (triphasic morphology), anoxic brain injury, metabolic encephalopathy, status epilepticus",
        "syndrome": "none",
        "significance": "In CJD: bilateral synchronous triphasic waves at 1 Hz; must correlate with clinical context",
    },
    {
        "pattern_id": "EEG-010",
        "name": "Lateralized Periodic Discharges (LPDs/PLEDs)",
        "category": "encephalopathy",
        "frequency": "0.5-3 Hz",
        "morphology": "Periodic sharp or spike-wave complexes lateralized to one hemisphere; often temporal or frontotemporal",
        "clinical_correlation": "Acute structural lesion (stroke, herpes encephalitis, tumor); associated with seizures",
        "syndrome": "none",
        "significance": "High seizure risk (>60%); HSV encephalitis: temporal LPDs + fever; continuous EEG monitoring recommended",
    },
    {
        "pattern_id": "EEG-011",
        "name": "Triphasic Waves",
        "category": "encephalopathy",
        "frequency": "1.5-2.5 Hz",
        "morphology": "Triphasic morphology with anterior-to-posterior lag; positive-negative-positive waveform",
        "clinical_correlation": "Hepatic encephalopathy, uremic encephalopathy, metabolic encephalopathy",
        "syndrome": "none",
        "significance": "Classic for hepatic encephalopathy but non-specific; overlap with NCSE pattern; trial of benzodiazepine may help differentiate",
    },
    {
        "pattern_id": "EEG-012",
        "name": "Burst-Suppression Pattern",
        "category": "encephalopathy",
        "frequency": "variable bursts",
        "morphology": "Alternating high-voltage bursts (0.5-10s) and periods of suppression (<10 uV); bursts may contain spikes",
        "clinical_correlation": "Severe diffuse cerebral dysfunction; post-cardiac arrest; deep anesthesia; neonatal encephalopathy",
        "syndrome": "none",
        "significance": "In post-cardiac arrest: identical bursts suggest poor prognosis; highly reactive pattern may be favorable",
    },
    {
        "pattern_id": "EEG-013",
        "name": "Electrocerebral Inactivity (ECI)",
        "category": "encephalopathy",
        "frequency": "none",
        "morphology": "No cerebral electrical activity >2 uV at sensitivity 2 uV/mm; recording at least 30 minutes with full montage",
        "clinical_correlation": "Brain death determination; requires confirmation of no confounders (hypothermia, drugs)",
        "syndrome": "none",
        "significance": "Part of brain death protocol; must exclude hypothermia, CNS depressants, metabolic derangements",
    },
    # --- Normal Variants ---
    {
        "pattern_id": "EEG-014",
        "name": "Wicket Spikes",
        "category": "normal_variant",
        "frequency": "6-11 Hz",
        "morphology": "Arciform (wicket-shaped) spikes over temporal regions; monophasic; no aftergoing slow wave; occurs in drowsiness",
        "clinical_correlation": "Normal variant; adults >30 years; may be mistaken for epileptiform temporal spikes",
        "syndrome": "none",
        "significance": "Benign; no association with seizures; key distinguishing feature: absence of aftergoing slow wave",
    },
    {
        "pattern_id": "EEG-015",
        "name": "14 and 6 Hz Positive Bursts",
        "category": "normal_variant",
        "frequency": "14 Hz and 6 Hz",
        "morphology": "Comb-like positive sharp waves at 14 or 6 Hz; posterior temporal predominance; during drowsiness/light sleep",
        "clinical_correlation": "Normal variant; adolescents; no clinical significance",
        "syndrome": "none",
        "significance": "Benign; frequently overread as abnormal; most common in adolescence",
    },
    {
        "pattern_id": "EEG-016",
        "name": "Rhythmic Mid-Temporal Theta of Drowsiness (RMTD)",
        "category": "normal_variant",
        "frequency": "4-7 Hz",
        "morphology": "Rhythmic theta bursts over mid-temporal regions during drowsiness; notched or sinusoidal; previously called psychomotor variant",
        "clinical_correlation": "Normal variant; no association with epilepsy",
        "syndrome": "none",
        "significance": "Benign; frequently misinterpreted as temporal lobe seizure pattern",
    },
    {
        "pattern_id": "EEG-017",
        "name": "SREDA (Subclinical Rhythmic EEG Discharge of Adults)",
        "category": "normal_variant",
        "frequency": "5-7 Hz",
        "morphology": "Widespread rhythmic theta/delta; sudden onset; parietal/posterior temporal; may evolve in frequency",
        "clinical_correlation": "Normal variant; elderly; no association with seizures despite seizure-like morphology",
        "syndrome": "none",
        "significance": "Benign; highly mimics subclinical seizure; no treatment required",
    },
    {
        "pattern_id": "EEG-018",
        "name": "Breach Rhythm",
        "category": "normal_variant",
        "frequency": "mu range or beta",
        "morphology": "Increased amplitude and apparent sharpness over skull defect; due to reduced impedance through craniotomy site",
        "clinical_correlation": "Post-craniotomy; not epileptiform",
        "syndrome": "none",
        "significance": "Must not be misinterpreted as epileptiform; document skull defect history",
    },
    # --- ICU / Critical Care Patterns ---
    {
        "pattern_id": "EEG-019",
        "name": "Nonconvulsive Status Epilepticus (NCSE)",
        "category": "ictal",
        "frequency": ">2.5 Hz",
        "morphology": "Continuous or near-continuous epileptiform activity without prominent motor manifestations; may be focal or generalized",
        "clinical_correlation": "Altered consciousness without convulsions; post-convulsive SE; critically ill patients",
        "syndrome": "none",
        "significance": "Requires cEEG monitoring; found in 8-34% of comatose ICU patients; IV benzodiazepine trial diagnostic and therapeutic",
    },
    {
        "pattern_id": "EEG-020",
        "name": "Generalized Rhythmic Delta Activity (GRDA)",
        "category": "encephalopathy",
        "frequency": "1-4 Hz",
        "morphology": "Frontally predominant rhythmic delta activity; continuous or semi-continuous; no evolution",
        "clinical_correlation": "Diffuse encephalopathy; increased ICP; deep midline lesions; metabolic derangements",
        "syndrome": "none",
        "significance": "ACNS standardized terminology; formerly FIRDA; associated with increased seizure risk when + sharp waves",
    },
    {
        "pattern_id": "EEG-021",
        "name": "Lateralized Rhythmic Delta Activity (LRDA)",
        "category": "encephalopathy",
        "frequency": "1-4 Hz",
        "morphology": "Rhythmic delta activity lateralized to one hemisphere; focal structural correlate",
        "clinical_correlation": "Focal structural lesion (tumor, stroke, abscess); higher seizure risk than GRDA",
        "syndrome": "none",
        "significance": "ACNS terminology; carries 50-60% seizure risk; continuous EEG monitoring warranted",
    },
    {
        "pattern_id": "EEG-022",
        "name": "Stimulus-Induced Rhythmic, Periodic, or Ictal Discharges (SIRPIDs)",
        "category": "ictal",
        "frequency": "variable",
        "morphology": "Rhythmic or periodic discharges induced by alerting stimuli (sternal rub, voice, suctioning)",
        "clinical_correlation": "Critically ill patients; uncertain clinical significance; may represent cortical irritability",
        "syndrome": "none",
        "significance": "Document stimulus relationship; clinical significance debated; may warrant treatment if associated with clinical change",
    },
    # --- Sleep EEG Patterns ---
    {
        "pattern_id": "EEG-023",
        "name": "Sleep Spindles",
        "category": "normal_sleep",
        "frequency": "11-16 Hz",
        "morphology": "Waxing-waning sinusoidal bursts; 0.5-2 seconds; vertex/central maximum; NREM stage 2",
        "clinical_correlation": "Normal sleep architecture; absence may indicate thalamocortical dysfunction",
        "syndrome": "none",
        "significance": "Spindle asymmetry or absence can indicate structural lesion; preserved spindles favorable in encephalopathy",
    },
    {
        "pattern_id": "EEG-024",
        "name": "K-Complexes",
        "category": "normal_sleep",
        "frequency": "isolated",
        "morphology": "High-amplitude biphasic wave with initial sharp negative deflection; vertex maximum; NREM stage 2; may be spontaneous or evoked",
        "clinical_correlation": "Normal sleep; arousal response; staging marker for N2 sleep",
        "syndrome": "none",
        "significance": "Required for N2 staging; absence may indicate cortical dysfunction; may occur with spindles",
    },
    {
        "pattern_id": "EEG-025",
        "name": "Vertex Sharp Waves",
        "category": "normal_sleep",
        "frequency": "isolated",
        "morphology": "High-amplitude biphasic sharp transients at vertex (Cz); maximal in N1 sleep; bilateral and symmetric",
        "clinical_correlation": "Normal sleep onset; N1 staging marker",
        "syndrome": "none",
        "significance": "Must not be mistaken for epileptiform discharges; normal sleep architecture",
    },
    # --- Specific Syndrome Patterns ---
    {
        "pattern_id": "EEG-026",
        "name": "Dravet Syndrome EEG Evolution",
        "category": "epileptiform",
        "frequency": "variable",
        "morphology": "Initially normal; generalized spike-wave by age 2; photosensitivity; multifocal spikes; focal slowing",
        "clinical_correlation": "SCN1A-related Dravet syndrome; febrile status epilepticus onset",
        "syndrome": "dravet",
        "significance": "Early EEG may be normal; photosensitivity develops in 40%; avoid sodium channel blockers",
    },
    {
        "pattern_id": "EEG-027",
        "name": "ESES/CSWS Pattern",
        "category": "epileptiform",
        "frequency": "1.5-3 Hz",
        "morphology": "Continuous spike-and-wave during slow sleep occupying >85% of NREM sleep; marked activation from wake to sleep",
        "clinical_correlation": "Epileptic encephalopathy with ESES; Landau-Kleffner syndrome; cognitive regression",
        "syndrome": "unclassified",
        "significance": "Cognitive deterioration correlates with spike-wave index; treatment targets EEG normalization; steroids, benzodiazepines, or surgery",
    },
    {
        "pattern_id": "EEG-028",
        "name": "Extreme Delta Brush",
        "category": "encephalopathy",
        "frequency": "1-3 Hz delta with 20-30 Hz beta",
        "morphology": "Rhythmic delta activity with superimposed fast beta activity (brush-like); diffuse or posterior predominant",
        "clinical_correlation": "Anti-NMDA receptor encephalitis; highly specific marker",
        "syndrome": "none",
        "significance": "Present in ~30% of anti-NMDAR encephalitis cases; may indicate poor short-term outcome but not necessarily long-term",
    },
    # --- Neonatal EEG Patterns ---
    {
        "pattern_id": "EEG-029",
        "name": "Trace Discontinue Pattern (Neonatal)",
        "category": "neonatal",
        "frequency": "intermittent bursts",
        "morphology": "Discontinuous background with interburst intervals 6-12 seconds; higher amplitude bursts; age-dependent",
        "clinical_correlation": "Premature neonates; becomes continuous by 34-36 weeks; prolonged IBI concerning for encephalopathy",
        "syndrome": "none",
        "significance": "Normal for gestational age <34 weeks; IBI >30 seconds abnormal at any age; guides HIE prognostication",
    },
    {
        "pattern_id": "EEG-030",
        "name": "Neonatal Seizure Patterns",
        "category": "ictal",
        "frequency": "0.5-4 Hz",
        "morphology": "Focal rhythmic sharp activity evolving in frequency/amplitude; minimum 10 seconds duration; may be subtle/subclinical",
        "clinical_correlation": "Neonatal seizures; HIE, stroke, metabolic, infection; 80% may be electrographic-only",
        "syndrome": "none",
        "significance": "Continuous EEG monitoring recommended in at-risk neonates; clinical seizures may be decoupled from EEG",
    },
    # --- Additional Patterns ---
    {
        "pattern_id": "EEG-031",
        "name": "Photoparoxysmal Response (PPR)",
        "category": "epileptiform",
        "frequency": "3-6 Hz",
        "morphology": "Generalized spike-wave or polyspike-wave discharges during photic stimulation; grades I-IV (Waltz classification)",
        "clinical_correlation": "Photosensitive epilepsy; JME; IGE; genetic generalized epilepsy",
        "syndrome": "juvenile_myoclonic",
        "significance": "Grade IV (self-sustained after stimulus) is most clinically significant; lifestyle counseling for triggers",
    },
    {
        "pattern_id": "EEG-032",
        "name": "Periodic Lateralized Epileptiform Discharges with + Fast Activity (LPDs+F)",
        "category": "ictal",
        "frequency": "0.5-2 Hz with superimposed fast",
        "morphology": "LPDs with superimposed fast (beta/gamma) activity; lateralized; indicates high ictal potential",
        "clinical_correlation": "Acute brain injury with seizures; higher seizure risk than LPDs alone",
        "syndrome": "none",
        "significance": "ACNS terminology; +F modifier strongly predicts electrographic seizures; aggressive monitoring and treatment warranted",
    },
    {
        "pattern_id": "EEG-033",
        "name": "Alpha Coma Pattern",
        "category": "encephalopathy",
        "frequency": "8-13 Hz",
        "morphology": "Unreactive alpha-frequency activity diffusely; no anterior-posterior gradient; no reactivity to stimulation",
        "clinical_correlation": "Post-cardiac arrest; pontine lesion; drug intoxication (benzodiazepines, barbiturates)",
        "syndrome": "none",
        "significance": "In post-cardiac arrest: generally poor prognosis; in drug intoxication: potentially reversible; reactivity testing critical",
    },
    {
        "pattern_id": "EEG-034",
        "name": "Frontal Intermittent Rhythmic Delta Activity (FIRDA/GRDA)",
        "category": "encephalopathy",
        "frequency": "1-4 Hz",
        "morphology": "Rhythmic delta activity maximal frontally; intermittent; may be notched (GRDA+S)",
        "clinical_correlation": "Metabolic encephalopathy; increased ICP; deep midline lesions; toxic encephalopathy",
        "syndrome": "none",
        "significance": "Now classified as GRDA in ACNS terminology; non-specific for diffuse cerebral dysfunction",
    },
    {
        "pattern_id": "EEG-035",
        "name": "Mu Rhythm",
        "category": "normal_variant",
        "frequency": "8-13 Hz",
        "morphology": "Arch-shaped (comb-like) rhythm over central regions; attenuated by contralateral movement or motor imagery; unilateral or bilateral",
        "clinical_correlation": "Normal awake rhythm; may be asymmetric; attenuated by movement or imagery",
        "syndrome": "none",
        "significance": "Normal variant; must not be mistaken for focal slowing or epileptiform discharge; BCI applications",
    },
    # --- Expanded EEG Patterns (10 new entries) ---
    {
        "pattern_id": "EEG-036",
        "name": "GRDA (Generalized Rhythmic Delta Activity)",
        "category": "encephalopathy",
        "frequency": "1-4 Hz",
        "morphology": "Generalized rhythmic delta activity, frontally predominant; continuous or semi-continuous; no spatial or temporal evolution; may have notched morphology (GRDA+S)",
        "clinical_correlation": "Diffuse encephalopathy; metabolic derangements; toxic encephalopathy; increased intracranial pressure; deep midline structural lesions",
        "syndrome": "none",
        "significance": "ACNS 2021 standardized terminology; replaces FIRDA; when accompanied by sharps (+S) or fast activity (+F), seizure risk increases significantly; requires clinical correlation for treatment decisions",
    },
    {
        "pattern_id": "EEG-037",
        "name": "LRDA (Lateralized Rhythmic Delta Activity)",
        "category": "encephalopathy",
        "frequency": "1-4 Hz",
        "morphology": "Rhythmic delta activity lateralized to one hemisphere or focal region; semi-rhythmic with consistent morphology; no definite ictal evolution",
        "clinical_correlation": "Focal structural pathology (tumor, stroke, abscess); higher seizure risk than GRDA; often ipsilateral to structural lesion",
        "syndrome": "none",
        "significance": "ACNS 2021 terminology; carries 50-60% risk of seizures; continuous EEG monitoring recommended; when plus-modified (+F, +S, +R), seizure risk approaches that of frank electrographic seizures",
    },
    {
        "pattern_id": "EEG-038",
        "name": "GPDs (Generalized Periodic Discharges)",
        "category": "encephalopathy",
        "frequency": "0.5-2 Hz",
        "morphology": "Periodic sharp waves or complexes repeating at regular intervals; generalized distribution; may have triphasic morphology; inter-discharge interval relatively fixed",
        "clinical_correlation": "Metabolic encephalopathy (hepatic, uremic); toxic exposure; anoxic brain injury; CJD (1 Hz triphasic); status epilepticus on the ictal-interictal continuum",
        "syndrome": "none",
        "significance": "Context-dependent interpretation critical; in CJD: 1 Hz bilaterally synchronous with characteristic triphasic morphology; in metabolic: often triphasic with anterior-posterior lag; trial of benzodiazepine may differentiate from NCSE",
    },
    {
        "pattern_id": "EEG-039",
        "name": "LPDs (Lateralized Periodic Discharges)",
        "category": "encephalopathy",
        "frequency": "0.5-3 Hz",
        "morphology": "Periodic sharp or spike-wave complexes lateralized to one hemisphere; consistent inter-discharge interval; may have plus modifiers (+F, +R, +S)",
        "clinical_correlation": "Acute structural lesion (acute stroke, herpes simplex encephalitis, brain abscess, tumor); associated with high seizure risk (>60%)",
        "syndrome": "none",
        "significance": "Previously called PLEDs; HSV encephalitis: temporal LPDs with fever is classic; LPDs+F (with fast activity) have highest seizure risk; continuous EEG monitoring essential; may warrant empiric anticonvulsant treatment",
    },
    {
        "pattern_id": "EEG-040",
        "name": "BIRDs (Brief Potentially Ictal Rhythmic Discharges)",
        "category": "ictal",
        "frequency": "variable",
        "morphology": "Very brief (<10 seconds) runs of rhythmic activity that resemble seizure onset patterns but do not meet duration criteria for electrographic seizures; focal or generalized",
        "clinical_correlation": "Critically ill patients; acute brain injury; significance on ictal-interictal continuum debated; may represent brief subclinical seizures",
        "syndrome": "none",
        "significance": "ACNS 2021 terminology; clinical significance remains uncertain; some studies associate with worse outcomes in critically ill; treatment decisions should be individualized based on clinical context",
    },
    {
        "pattern_id": "EEG-041",
        "name": "Burst-Suppression Pattern",
        "category": "encephalopathy",
        "frequency": "variable bursts with suppression periods",
        "morphology": "Alternating periods of high-voltage mixed-frequency activity (bursts, 0.5-10s) and isoelectric or near-isoelectric suppression (<10 uV); burst content may include spikes or sharp waves",
        "clinical_correlation": "Deep sedation (barbiturate/propofol coma); post-cardiac arrest anoxic injury; severe diffuse cerebral dysfunction; therapeutic target in refractory status epilepticus",
        "syndrome": "none",
        "significance": "In therapeutic coma for refractory SE: burst-suppression is the treatment target; in post-cardiac arrest: identical/stereotyped bursts suggest poor prognosis; highly reactive pattern may indicate better outcome; burst-suppression ratio quantifiable",
    },
    {
        "pattern_id": "EEG-042",
        "name": "Alpha Coma Pattern",
        "category": "encephalopathy",
        "frequency": "8-13 Hz",
        "morphology": "Diffuse unreactive alpha-frequency activity without normal anterior-posterior gradient; no reactivity to noxious or auditory stimulation; monotonous and invariant",
        "clinical_correlation": "Post-cardiac arrest (most common); pontine/brainstem lesion; severe drug intoxication (benzodiazepines, barbiturates, opioids)",
        "syndrome": "none",
        "significance": "In post-cardiac arrest: generally indicates poor prognosis with diffuse cortical injury; in drug intoxication: potentially fully reversible; reactivity testing is crucial for prognostication; must not be confused with normal alpha rhythm",
    },
    {
        "pattern_id": "EEG-043",
        "name": "Spindle Coma Pattern",
        "category": "encephalopathy",
        "frequency": "11-16 Hz spindle-like activity",
        "morphology": "Sleep spindle-like activity in a comatose patient; vertex waves and K-complex-like waveforms may also be present; suggests preserved thalamocortical circuitry",
        "clinical_correlation": "Diffuse cortical injury with preserved brainstem and thalamic function; post-traumatic brain injury; hypoxic-ischemic encephalopathy",
        "syndrome": "none",
        "significance": "Generally more favorable prognosis than alpha or theta coma; indicates intact thalamocortical connections; may evolve to normal sleep patterns during recovery; serial EEG monitoring tracks evolution",
    },
    {
        "pattern_id": "EEG-044",
        "name": "ECI (Electrocerebral Inactivity)",
        "category": "encephalopathy",
        "frequency": "none (isoelectric)",
        "morphology": "No discernible cerebral electrical activity >2 uV; recording at sensitivity 2 uV/mm; minimum 30-minute recording; full 10-20 montage with interelectrode distances >10 cm",
        "clinical_correlation": "Brain death determination; must confirm absence of confounders: hypothermia (<36C), CNS depressant drugs, severe metabolic derangement",
        "syndrome": "none",
        "significance": "Confirmatory test in brain death protocol per AAN guidelines; technical requirements: interelectrode impedance 100-10,000 ohms; must test reactivity; two recordings 6+ hours apart recommended in some protocols; not sufficient alone for brain death determination",
    },
    {
        "pattern_id": "EEG-045",
        "name": "Wicket Spikes",
        "category": "normal_variant",
        "frequency": "6-11 Hz",
        "morphology": "Arciform (wicket-shaped) spikes or sharp transients over temporal regions; monophasic negative waves; no aftergoing slow wave; occur in trains during drowsiness or light sleep; unilateral or bilateral",
        "clinical_correlation": "Benign normal variant; adults >30 years; frequently misidentified as epileptiform temporal sharp waves leading to inappropriate epilepsy diagnosis",
        "syndrome": "none",
        "significance": "Key distinguishing features from true epileptiform discharges: absence of aftergoing slow wave, absence of disruption of background, occurrence in rhythmic trains; no association with seizure risk; no treatment required; one of the most commonly over-interpreted EEG patterns",
    },
]


# ===================================================================
# EEG PARSER IMPLEMENTATION
# ===================================================================


class EEGParser(BaseIngestParser):
    """Parse EEG patterns and findings for the Neurology Intelligence Agent.

    In seed mode, returns the curated EEG_PATTERNS list.

    Usage::

        parser = EEGParser()
        records, stats = parser.run()
    """

    def __init__(
        self,
        collection_manager: Any = None,
        embedder: Any = None,
    ) -> None:
        super().__init__(
            source_name="eeg",
            collection_manager=collection_manager,
            embedder=embedder,
        )

    def fetch(self, **kwargs: Any) -> List[Dict[str, Any]]:
        """Fetch EEG pattern data.

        Returns:
            List of raw EEG pattern dictionaries.
        """
        self.logger.info(
            "Using curated EEG pattern seed data (%d patterns)",
            len(EEG_PATTERNS),
        )
        return list(EEG_PATTERNS)

    def parse(self, raw_data: List[Dict[str, Any]]) -> List[IngestRecord]:
        """Parse raw EEG data into IngestRecord objects.

        Args:
            raw_data: List of EEG pattern dictionaries.

        Returns:
            List of IngestRecord objects.
        """
        records: List[IngestRecord] = []

        for entry in raw_data:
            pattern_id = entry.get("pattern_id", "")
            name = entry.get("name", "")
            category = entry.get("category", "")
            frequency = entry.get("frequency", "")
            morphology = entry.get("morphology", "")
            clinical_correlation = entry.get("clinical_correlation", "")
            syndrome = entry.get("syndrome", "none")
            significance = entry.get("significance", "")

            text = (
                f"EEG Pattern: {name}. "
                f"Category: {category}. "
                f"Frequency: {frequency}. "
                f"Morphology: {morphology}. "
                f"Clinical correlation: {clinical_correlation}. "
                f"Significance: {significance}."
            )

            record = IngestRecord(
                text=text,
                metadata={
                    "pattern_id": pattern_id,
                    "name": name,
                    "category": category,
                    "frequency": frequency,
                    "morphology": morphology,
                    "syndrome": syndrome,
                    "source_db": "eeg_seed",
                },
                collection_name="neuro_electrophysiology",
                record_id=pattern_id,
                source="eeg",
            )
            records.append(record)

        return records

    def validate_record(self, record: IngestRecord) -> bool:
        """Validate an EEG IngestRecord.

        Requirements:
            - text must be non-empty
            - must have pattern_id in metadata
            - must have name in metadata
            - must have category in metadata

        Args:
            record: The record to validate.

        Returns:
            True if the record passes all validation checks.
        """
        if not record.text or not record.text.strip():
            return False

        meta = record.metadata
        if not meta.get("pattern_id"):
            return False
        if not meta.get("name"):
            return False
        if not meta.get("category"):
            return False

        return True


def get_eeg_pattern_count() -> int:
    """Return the number of curated EEG patterns."""
    return len(EEG_PATTERNS)


def get_eeg_categories() -> List[str]:
    """Return a deduplicated sorted list of EEG categories."""
    categories = list({p["category"] for p in EEG_PATTERNS})
    categories.sort()
    return categories


def get_eeg_syndromes() -> List[str]:
    """Return a deduplicated sorted list of EEG-associated syndromes (excluding 'none')."""
    syndromes = list({p["syndrome"] for p in EEG_PATTERNS if p["syndrome"] != "none"})
    syndromes.sort()
    return syndromes
