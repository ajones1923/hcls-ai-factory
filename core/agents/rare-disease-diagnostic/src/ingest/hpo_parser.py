"""HPO (Human Phenotype Ontology) parser for the Rare Disease Diagnostic Agent.

Parses HPO ontology terms and seeds the top 150+ HPO terms with structured
data including HPO ID, name, definition, synonyms, and information content
(IC) scores for phenotype-driven diagnosis.

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List

from .base import BaseIngestParser, IngestRecord

logger = logging.getLogger(__name__)


# ===================================================================
# SEED DATA: TOP 150+ HPO TERMS
# ===================================================================

HPO_TERMS: List[Dict[str, Any]] = [
    {"hpo_id": "HP:0001250", "name": "Seizures", "definition": "Repeated occurrence of seizures", "synonyms": ["Epilepsy", "Convulsions", "Seizure disorder"], "ic_score": 2.1},
    {"hpo_id": "HP:0001252", "name": "Hypotonia", "definition": "Decreased muscle tone", "synonyms": ["Muscular hypotonia", "Low muscle tone", "Floppy infant"], "ic_score": 2.3},
    {"hpo_id": "HP:0001249", "name": "Intellectual disability", "definition": "Significant limitations in intellectual functioning", "synonyms": ["Mental retardation", "Cognitive impairment", "Learning disability"], "ic_score": 1.9},
    {"hpo_id": "HP:0000252", "name": "Microcephaly", "definition": "Head circumference more than 2 SD below mean", "synonyms": ["Small head", "Reduced head circumference"], "ic_score": 3.2},
    {"hpo_id": "HP:0001263", "name": "Global developmental delay", "definition": "Delay in achieving developmental milestones", "synonyms": ["Developmental delay", "Psychomotor retardation"], "ic_score": 1.8},
    {"hpo_id": "HP:0000256", "name": "Macrocephaly", "definition": "Head circumference more than 2 SD above mean", "synonyms": ["Large head", "Big head", "Megalocephaly"], "ic_score": 3.1},
    {"hpo_id": "HP:0002721", "name": "Immunodeficiency", "definition": "Impairment of the immune system", "synonyms": ["Immune deficiency", "Immunological deficiency"], "ic_score": 3.5},
    {"hpo_id": "HP:0001508", "name": "Failure to thrive", "definition": "Inadequate physical growth", "synonyms": ["Poor weight gain", "Growth failure", "FTT"], "ic_score": 2.5},
    {"hpo_id": "HP:0001510", "name": "Growth delay", "definition": "Growth below expected for age", "synonyms": ["Short stature", "Growth retardation", "Growth deficiency"], "ic_score": 2.4},
    {"hpo_id": "HP:0000365", "name": "Hearing impairment", "definition": "Reduced ability to detect sounds", "synonyms": ["Hearing loss", "Deafness", "Hypoacusis"], "ic_score": 2.8},
    {"hpo_id": "HP:0000639", "name": "Nystagmus", "definition": "Rhythmic involuntary oscillations of the eyes", "synonyms": ["Involuntary eye movements", "Eye oscillation"], "ic_score": 3.0},
    {"hpo_id": "HP:0000486", "name": "Strabismus", "definition": "Misalignment of the eyes", "synonyms": ["Squint", "Cross-eyed", "Crossed eyes"], "ic_score": 3.1},
    {"hpo_id": "HP:0000518", "name": "Cataract", "definition": "Opacity of the crystalline lens", "synonyms": ["Cataracts", "Lens opacity"], "ic_score": 3.3},
    {"hpo_id": "HP:0001156", "name": "Brachydactyly", "definition": "Short digits", "synonyms": ["Short fingers", "Short toes", "Stubby digits"], "ic_score": 3.8},
    {"hpo_id": "HP:0001166", "name": "Arachnodactyly", "definition": "Abnormally long and slender fingers", "synonyms": ["Spider fingers", "Long slender fingers"], "ic_score": 4.2},
    {"hpo_id": "HP:0002650", "name": "Scoliosis", "definition": "Lateral curvature of the spine", "synonyms": ["Spinal curvature", "Curved spine"], "ic_score": 2.7},
    {"hpo_id": "HP:0001371", "name": "Flexion contracture", "definition": "Permanent fixation of a joint in a bent position", "synonyms": ["Joint contracture", "Fixed flexion"], "ic_score": 3.0},
    {"hpo_id": "HP:0001382", "name": "Joint hypermobility", "definition": "Increased range of motion of joints", "synonyms": ["Hypermobile joints", "Lax joints", "Joint laxity"], "ic_score": 2.9},
    {"hpo_id": "HP:0002240", "name": "Hepatomegaly", "definition": "Enlargement of the liver", "synonyms": ["Enlarged liver", "Liver enlargement"], "ic_score": 2.6},
    {"hpo_id": "HP:0001744", "name": "Splenomegaly", "definition": "Enlargement of the spleen", "synonyms": ["Enlarged spleen", "Spleen enlargement"], "ic_score": 3.0},
    {"hpo_id": "HP:0001627", "name": "Abnormal heart morphology", "definition": "Structural abnormality of the heart", "synonyms": ["Cardiac malformation", "Heart defect", "Congenital heart defect"], "ic_score": 2.2},
    {"hpo_id": "HP:0001631", "name": "Atrial septal defect", "definition": "Defect in the interatrial septum", "synonyms": ["ASD", "Hole in the heart"], "ic_score": 3.5},
    {"hpo_id": "HP:0001629", "name": "Ventricular septal defect", "definition": "Defect in the interventricular septum", "synonyms": ["VSD"], "ic_score": 3.4},
    {"hpo_id": "HP:0001638", "name": "Cardiomyopathy", "definition": "Disease of the heart muscle", "synonyms": ["Myocardial disease", "Heart muscle disease"], "ic_score": 2.8},
    {"hpo_id": "HP:0000098", "name": "Tall stature", "definition": "Height above 97th percentile", "synonyms": ["Overgrowth", "Tall", "Large stature"], "ic_score": 3.2},
    {"hpo_id": "HP:0000347", "name": "Micrognathia", "definition": "Abnormally small mandible", "synonyms": ["Small jaw", "Receding chin", "Mandibular hypoplasia"], "ic_score": 3.4},
    {"hpo_id": "HP:0000175", "name": "Cleft palate", "definition": "Incomplete closure of the palate", "synonyms": ["Palatal cleft", "Cleft of palate"], "ic_score": 3.3},
    {"hpo_id": "HP:0000028", "name": "Cryptorchidism", "definition": "Undescended testes", "synonyms": ["Undescended testis", "Testicular maldescent"], "ic_score": 3.1},
    {"hpo_id": "HP:0000951", "name": "Abnormality of the skin", "definition": "Skin abnormality", "synonyms": ["Skin defect", "Dermatological abnormality"], "ic_score": 1.5},
    {"hpo_id": "HP:0000957", "name": "Cafe-au-lait spot", "definition": "Light brown flat skin lesion", "synonyms": ["CAL spot", "Coffee-colored macule"], "ic_score": 3.6},
    {"hpo_id": "HP:0001263", "name": "Global developmental delay", "definition": "A delay in achieving developmental milestones across multiple domains", "synonyms": ["GDD", "Psychomotor delay"], "ic_score": 1.8},
    {"hpo_id": "HP:0001290", "name": "Generalized hypotonia", "definition": "Reduced muscle tone involving all muscle groups", "synonyms": ["Global hypotonia", "Diffuse hypotonia"], "ic_score": 2.5},
    {"hpo_id": "HP:0002020", "name": "Gastroesophageal reflux", "definition": "Retrograde flow of gastric contents into esophagus", "synonyms": ["GERD", "Acid reflux", "GER"], "ic_score": 2.7},
    {"hpo_id": "HP:0002093", "name": "Respiratory insufficiency", "definition": "Inadequate ventilation to maintain normal blood gases", "synonyms": ["Respiratory failure", "Ventilatory insufficiency"], "ic_score": 2.6},
    {"hpo_id": "HP:0001270", "name": "Motor delay", "definition": "Delay in motor development", "synonyms": ["Delayed motor milestones", "Motor retardation"], "ic_score": 2.2},
    {"hpo_id": "HP:0001344", "name": "Absent speech", "definition": "Complete absence of speech development", "synonyms": ["Nonverbal", "No speech"], "ic_score": 3.5},
    {"hpo_id": "HP:0000729", "name": "Autistic behavior", "definition": "Behavioral features characteristic of autism spectrum disorder", "synonyms": ["Autism", "ASD features", "Autistic traits"], "ic_score": 2.9},
    {"hpo_id": "HP:0001256", "name": "Intellectual disability, mild", "definition": "IQ between 50 and 69", "synonyms": ["Mild mental retardation", "Mild cognitive impairment"], "ic_score": 2.8},
    {"hpo_id": "HP:0002342", "name": "Intellectual disability, moderate", "definition": "IQ between 35 and 49", "synonyms": ["Moderate mental retardation"], "ic_score": 3.2},
    {"hpo_id": "HP:0010864", "name": "Intellectual disability, severe", "definition": "IQ between 20 and 34", "synonyms": ["Severe mental retardation"], "ic_score": 3.6},
    {"hpo_id": "HP:0000707", "name": "Abnormality of the nervous system", "definition": "Neurological abnormality", "synonyms": ["CNS abnormality", "Neurological defect"], "ic_score": 1.3},
    {"hpo_id": "HP:0001939", "name": "Abnormality of metabolism/homeostasis", "definition": "Metabolic abnormality", "synonyms": ["Metabolic disorder", "Inborn error of metabolism"], "ic_score": 1.4},
    {"hpo_id": "HP:0003077", "name": "Hyperlipidemia", "definition": "Elevated blood lipid levels", "synonyms": ["High cholesterol", "Dyslipidemia"], "ic_score": 3.0},
    {"hpo_id": "HP:0001943", "name": "Hypoglycemia", "definition": "Abnormally low blood glucose", "synonyms": ["Low blood sugar", "Hypoglycaemia"], "ic_score": 2.8},
    {"hpo_id": "HP:0003155", "name": "Elevated alkaline phosphatase", "definition": "Increased serum alkaline phosphatase", "synonyms": ["High ALP", "Elevated ALP"], "ic_score": 3.3},
    {"hpo_id": "HP:0003236", "name": "Elevated creatine kinase", "definition": "Increased serum creatine kinase", "synonyms": ["High CK", "Elevated CK", "Hyperckemia"], "ic_score": 3.1},
    {"hpo_id": "HP:0001903", "name": "Anemia", "definition": "Reduction in hemoglobin or red blood cell count", "synonyms": ["Anaemia", "Low hemoglobin"], "ic_score": 2.4},
    {"hpo_id": "HP:0001873", "name": "Thrombocytopenia", "definition": "Abnormally low platelet count", "synonyms": ["Low platelets", "Platelet deficiency"], "ic_score": 2.9},
    {"hpo_id": "HP:0100543", "name": "Cognitive impairment", "definition": "Deficit in cognitive function", "synonyms": ["Cognitive decline", "Cognitive deficit"], "ic_score": 2.0},
    {"hpo_id": "HP:0002360", "name": "Sleep disturbance", "definition": "Abnormality of sleep", "synonyms": ["Sleep disorder", "Insomnia", "Sleep problems"], "ic_score": 2.5},
    {"hpo_id": "HP:0001251", "name": "Ataxia", "definition": "Lack of voluntary coordination of muscle movements", "synonyms": ["Cerebellar ataxia", "Gait ataxia", "Uncoordinated movement"], "ic_score": 2.6},
    {"hpo_id": "HP:0001257", "name": "Spasticity", "definition": "Velocity-dependent increase in muscle tone", "synonyms": ["Hypertonia", "Increased muscle tone", "Spastic"], "ic_score": 2.5},
    {"hpo_id": "HP:0000238", "name": "Hydrocephalus", "definition": "Accumulation of cerebrospinal fluid within the brain", "synonyms": ["Water on the brain", "CSF accumulation"], "ic_score": 3.4},
    {"hpo_id": "HP:0002119", "name": "Ventriculomegaly", "definition": "Enlarged cerebral ventricles", "synonyms": ["Dilated ventricles", "Ventricular enlargement"], "ic_score": 3.3},
    {"hpo_id": "HP:0002126", "name": "Polymicrogyria", "definition": "Excessive number of small convolutions on the brain surface", "synonyms": ["PMG", "Too many small gyri"], "ic_score": 5.2},
    {"hpo_id": "HP:0007370", "name": "Aplasia/Hypoplasia of the corpus callosum", "definition": "Absent or underdeveloped corpus callosum", "synonyms": ["Agenesis of corpus callosum", "ACC", "Thin corpus callosum"], "ic_score": 3.8},
    {"hpo_id": "HP:0000750", "name": "Delayed speech and language development", "definition": "Slower than normal development of speech and language", "synonyms": ["Speech delay", "Language delay", "Late talker"], "ic_score": 2.1},
    {"hpo_id": "HP:0002079", "name": "Hypoplasia of the corpus callosum", "definition": "Underdevelopment of the corpus callosum", "synonyms": ["Thin corpus callosum"], "ic_score": 4.0},
    {"hpo_id": "HP:0010804", "name": "Tented upper lip vermilion", "definition": "Upper lip with a tented appearance", "synonyms": ["Tented upper lip"], "ic_score": 5.5},
    {"hpo_id": "HP:0000431", "name": "Wide nasal bridge", "definition": "Increased breadth of the nasal bridge", "synonyms": ["Broad nasal bridge", "Flat nasal bridge"], "ic_score": 3.5},
    {"hpo_id": "HP:0000316", "name": "Hypertelorism", "definition": "Increased distance between the eyes", "synonyms": ["Wide-set eyes", "Widely spaced eyes"], "ic_score": 3.4},
    {"hpo_id": "HP:0000508", "name": "Ptosis", "definition": "Drooping of the upper eyelid", "synonyms": ["Droopy eyelid", "Blepharoptosis"], "ic_score": 3.2},
    {"hpo_id": "HP:0000369", "name": "Low-set ears", "definition": "Position of ears below normal", "synonyms": ["Low-set pinnae", "Lowly set ears"], "ic_score": 3.5},
    {"hpo_id": "HP:0001999", "name": "Abnormal facial shape", "definition": "Unusual shape or structure of the face", "synonyms": ["Dysmorphic facies", "Coarse facies", "Facial dysmorphism"], "ic_score": 2.0},
    {"hpo_id": "HP:0000272", "name": "Malar flattening", "definition": "Reduced prominence of the malar region", "synonyms": ["Flat cheekbones", "Malar hypoplasia"], "ic_score": 3.7},
    {"hpo_id": "HP:0000445", "name": "Wide nose", "definition": "Increased width of the nose", "synonyms": ["Broad nose"], "ic_score": 3.4},
    {"hpo_id": "HP:0000194", "name": "Open mouth", "definition": "Inability to keep the mouth closed", "synonyms": ["Mouth gaping", "Open-mouthed appearance"], "ic_score": 3.8},
    {"hpo_id": "HP:0011220", "name": "Prominent forehead", "definition": "Forward protrusion of the forehead", "synonyms": ["Frontal bossing", "Bossed forehead"], "ic_score": 3.3},
    {"hpo_id": "HP:0000400", "name": "Macrotia", "definition": "Abnormally large ears", "synonyms": ["Large ears", "Big ears"], "ic_score": 3.9},
    {"hpo_id": "HP:0010880", "name": "Increased nuchal translucency", "definition": "Abnormally increased NT on prenatal ultrasound", "synonyms": ["Increased NT", "Thickened nuchal fold"], "ic_score": 4.5},
    {"hpo_id": "HP:0001562", "name": "Oligohydramnios", "definition": "Decreased amniotic fluid volume", "synonyms": ["Low amniotic fluid"], "ic_score": 3.9},
    {"hpo_id": "HP:0001561", "name": "Polyhydramnios", "definition": "Excess amniotic fluid volume", "synonyms": ["Too much amniotic fluid", "Hydramnios"], "ic_score": 3.6},
    {"hpo_id": "HP:0001622", "name": "Premature birth", "definition": "Birth before 37 weeks of gestation", "synonyms": ["Preterm birth", "Prematurity"], "ic_score": 2.5},
    {"hpo_id": "HP:0001518", "name": "Small for gestational age", "definition": "Birth weight below 10th percentile for gestational age", "synonyms": ["SGA", "Low birth weight", "IUGR"], "ic_score": 2.8},
    {"hpo_id": "HP:0003270", "name": "Abdominal distension", "definition": "Enlargement of the abdomen", "synonyms": ["Distended abdomen", "Bloating"], "ic_score": 2.7},
    {"hpo_id": "HP:0002910", "name": "Elevated hepatic transaminases", "definition": "Increased liver enzymes AST and ALT", "synonyms": ["Elevated liver enzymes", "High transaminases"], "ic_score": 2.8},
    {"hpo_id": "HP:0000119", "name": "Abnormality of the genitourinary system", "definition": "Genitourinary abnormality", "synonyms": ["GU malformation"], "ic_score": 1.6},
    {"hpo_id": "HP:0000107", "name": "Renal cyst", "definition": "Fluid-filled sac in the kidney", "synonyms": ["Kidney cyst", "Cystic kidney"], "ic_score": 3.4},
    {"hpo_id": "HP:0000083", "name": "Renal insufficiency", "definition": "Decreased kidney function", "synonyms": ["Kidney failure", "CKD", "Chronic kidney disease"], "ic_score": 2.7},
    {"hpo_id": "HP:0001395", "name": "Hepatic fibrosis", "definition": "Scarring of the liver", "synonyms": ["Liver fibrosis", "Hepatic scarring"], "ic_score": 3.5},
    {"hpo_id": "HP:0002910", "name": "Elevated hepatic transaminases", "definition": "Increased serum levels of hepatic transaminases", "synonyms": ["High liver enzymes"], "ic_score": 2.8},
    {"hpo_id": "HP:0000978", "name": "Bruising susceptibility", "definition": "Tendency to bruise easily", "synonyms": ["Easy bruising", "Ecchymoses"], "ic_score": 2.9},
    {"hpo_id": "HP:0001928", "name": "Abnormality of coagulation", "definition": "Abnormality of the blood clotting system", "synonyms": ["Coagulopathy", "Bleeding disorder", "Clotting abnormality"], "ic_score": 2.3},
    {"hpo_id": "HP:0004322", "name": "Short stature", "definition": "Height below 3rd percentile for age and sex", "synonyms": ["Dwarfism", "Growth deficiency", "Small stature"], "ic_score": 2.0},
    {"hpo_id": "HP:0000821", "name": "Hypothyroidism", "definition": "Deficiency of thyroid hormone", "synonyms": ["Underactive thyroid", "Low thyroid"], "ic_score": 2.9},
    {"hpo_id": "HP:0000836", "name": "Hyperthyroidism", "definition": "Excess thyroid hormone", "synonyms": ["Overactive thyroid", "Thyrotoxicosis"], "ic_score": 3.2},
    {"hpo_id": "HP:0000938", "name": "Osteopenia", "definition": "Decreased bone mineral density", "synonyms": ["Low bone density", "Thin bones"], "ic_score": 3.0},
    {"hpo_id": "HP:0002757", "name": "Recurrent fractures", "definition": "Multiple bone fractures with minimal trauma", "synonyms": ["Frequent fractures", "Pathological fractures", "Bone fragility"], "ic_score": 3.2},
    {"hpo_id": "HP:0001290", "name": "Generalized hypotonia", "definition": "Reduced tone affecting all muscle groups", "synonyms": ["Global hypotonia"], "ic_score": 2.5},
    {"hpo_id": "HP:0003202", "name": "Skeletal muscle atrophy", "definition": "Loss of skeletal muscle mass", "synonyms": ["Muscle wasting", "Muscle atrophy"], "ic_score": 2.8},
    {"hpo_id": "HP:0003701", "name": "Proximal muscle weakness", "definition": "Weakness affecting proximal muscles", "synonyms": ["Limb-girdle weakness", "Proximal myopathy"], "ic_score": 2.9},
    {"hpo_id": "HP:0002515", "name": "Waddling gait", "definition": "Gait characterized by lateral trunk movements", "synonyms": ["Duck-like walk", "Myopathic gait"], "ic_score": 3.5},
    {"hpo_id": "HP:0003198", "name": "Myopathy", "definition": "Muscle disease", "synonyms": ["Muscle disorder", "Muscular disease"], "ic_score": 2.3},
    {"hpo_id": "HP:0002878", "name": "Respiratory failure", "definition": "Inability of the respiratory system to maintain adequate gas exchange", "synonyms": ["Lung failure", "Ventilatory failure"], "ic_score": 2.4},
    {"hpo_id": "HP:0030680", "name": "Abnormality of cardiovascular system morphology", "definition": "Structural abnormality of the cardiovascular system", "synonyms": ["CV malformation", "Cardiac structural defect"], "ic_score": 1.6},
    {"hpo_id": "HP:0001537", "name": "Umbilical hernia", "definition": "Protrusion through the umbilical ring", "synonyms": ["Navel hernia"], "ic_score": 3.3},
    {"hpo_id": "HP:0000768", "name": "Pectus carinatum", "definition": "Protrusion of the anterior chest wall", "synonyms": ["Pigeon chest", "Protruding chest"], "ic_score": 3.7},
    {"hpo_id": "HP:0000767", "name": "Pectus excavatum", "definition": "Depression of the anterior chest wall", "synonyms": ["Funnel chest", "Sunken chest"], "ic_score": 3.5},
    {"hpo_id": "HP:0001548", "name": "Overgrowth", "definition": "Excessive growth compared to expected", "synonyms": ["Somatic overgrowth", "Gigantism"], "ic_score": 3.0},
    {"hpo_id": "HP:0001007", "name": "Hirsutism", "definition": "Excessive terminal hair growth", "synonyms": ["Excess body hair", "Hypertrichosis"], "ic_score": 3.4},
    {"hpo_id": "HP:0000023", "name": "Inguinal hernia", "definition": "Protrusion of abdominal contents through the inguinal canal", "synonyms": ["Groin hernia"], "ic_score": 3.2},
    # --- Additional clinically important HPO terms ---
    {"hpo_id": "HP:0000545", "name": "Myopia", "definition": "Nearsightedness due to excessive axial length of the eye", "synonyms": ["Nearsightedness", "Short-sightedness"], "ic_score": 2.8},
    {"hpo_id": "HP:0003073", "name": "Hypoalbuminemia", "definition": "Abnormally low serum albumin level", "synonyms": ["Low albumin", "Reduced serum albumin"], "ic_score": 3.1},
    {"hpo_id": "HP:0001324", "name": "Muscle weakness", "definition": "Reduced ability to generate force with skeletal muscles", "synonyms": ["Weakness", "Muscular weakness", "Motor weakness"], "ic_score": 2.0},
    {"hpo_id": "HP:0002015", "name": "Dysphagia", "definition": "Difficulty in swallowing", "synonyms": ["Swallowing difficulty", "Difficulty swallowing"], "ic_score": 2.6},
    {"hpo_id": "HP:0000076", "name": "Vesicoureteral reflux", "definition": "Retrograde flow of urine from the bladder into the ureter", "synonyms": ["VUR", "Ureteral reflux"], "ic_score": 3.5},
    {"hpo_id": "HP:0001635", "name": "Congestive heart failure", "definition": "Inability of the heart to pump blood adequately", "synonyms": ["Heart failure", "CHF", "Cardiac failure"], "ic_score": 2.7},
    {"hpo_id": "HP:0001250", "name": "Seizures (tonic-clonic)", "definition": "Generalized tonic-clonic seizures with loss of consciousness", "synonyms": ["Grand mal seizures", "Convulsions"], "ic_score": 2.3},
    {"hpo_id": "HP:0002069", "name": "Bilateral tonic-clonic seizure", "definition": "Seizure with bilateral tonic and clonic phases", "synonyms": ["Generalized seizure", "Grand mal"], "ic_score": 2.5},
    {"hpo_id": "HP:0011153", "name": "Focal seizure", "definition": "Seizure originating from one hemisphere", "synonyms": ["Partial seizure", "Focal onset seizure"], "ic_score": 2.7},
    {"hpo_id": "HP:0002123", "name": "Generalized myoclonic seizure", "definition": "Brief shock-like involuntary muscle jerks", "synonyms": ["Myoclonic jerks", "Myoclonus epilepsy"], "ic_score": 3.0},
    {"hpo_id": "HP:0002197", "name": "Generalized-onset seizure", "definition": "Seizure with initial activation of both hemispheres", "synonyms": ["Generalized epilepsy"], "ic_score": 2.2},
    {"hpo_id": "HP:0000716", "name": "Depression", "definition": "Persistent feeling of sadness or loss of interest", "synonyms": ["Depressive disorder", "Major depression", "Low mood"], "ic_score": 2.3},
    {"hpo_id": "HP:0000739", "name": "Anxiety", "definition": "Excessive worry or apprehension", "synonyms": ["Anxiety disorder", "Anxiousness"], "ic_score": 2.4},
    {"hpo_id": "HP:0001328", "name": "Specific learning disability", "definition": "Impairment in a specific area of learning", "synonyms": ["Learning disorder", "Dyslexia"], "ic_score": 2.8},
    {"hpo_id": "HP:0002360", "name": "Sleep abnormality", "definition": "Disrupted or abnormal sleep pattern", "synonyms": ["Sleep disorder", "Dysomnia"], "ic_score": 2.5},
    {"hpo_id": "HP:0001635", "name": "Heart failure", "definition": "Inability of the heart to maintain adequate circulation", "synonyms": ["Cardiac insufficiency", "CHF"], "ic_score": 2.7},
    {"hpo_id": "HP:0001250", "name": "Epileptic spasms", "definition": "Clusters of brief tonic contractions in infancy", "synonyms": ["Infantile spasms", "West syndrome"], "ic_score": 3.2},
    {"hpo_id": "HP:0002133", "name": "Status epilepticus", "definition": "Prolonged seizure lasting more than 5 minutes", "synonyms": ["Continuous seizure", "Prolonged seizure"], "ic_score": 3.0},
    {"hpo_id": "HP:0001276", "name": "Hypertonia", "definition": "Increased resistance to passive movement", "synonyms": ["Increased tone", "Rigidity"], "ic_score": 2.4},
    {"hpo_id": "HP:0002059", "name": "Cerebral atrophy", "definition": "Loss of brain tissue volume", "synonyms": ["Brain atrophy", "Cerebral volume loss"], "ic_score": 3.3},
    {"hpo_id": "HP:0002500", "name": "Abnormal cerebral white matter morphology", "definition": "Structural abnormality of cerebral white matter", "synonyms": ["White matter abnormality", "Leukoencephalopathy"], "ic_score": 3.0},
    {"hpo_id": "HP:0000602", "name": "Ophthalmoplegia", "definition": "Paralysis of extraocular muscles", "synonyms": ["Eye movement paralysis", "Gaze palsy"], "ic_score": 3.5},
    {"hpo_id": "HP:0000648", "name": "Optic atrophy", "definition": "Degeneration of the optic nerve", "synonyms": ["Optic nerve atrophy", "Optic nerve pallor"], "ic_score": 3.4},
    {"hpo_id": "HP:0000589", "name": "Coloboma", "definition": "Developmental defect of the eye", "synonyms": ["Ocular coloboma", "Iris coloboma"], "ic_score": 4.0},
    {"hpo_id": "HP:0000568", "name": "Microphthalmos", "definition": "Abnormally small eye", "synonyms": ["Small eye", "Microphthalmia"], "ic_score": 4.2},
    {"hpo_id": "HP:0001263", "name": "Psychomotor retardation", "definition": "Slowing of thought and physical movement", "synonyms": ["Psychomotor slowing", "Mental and motor delay"], "ic_score": 2.0},
    {"hpo_id": "HP:0000112", "name": "Nephropathy", "definition": "Disease of the kidney", "synonyms": ["Kidney disease", "Renal disease"], "ic_score": 2.5},
    {"hpo_id": "HP:0000093", "name": "Proteinuria", "definition": "Excess protein in urine", "synonyms": ["Protein in urine", "Albuminuria"], "ic_score": 2.8},
    {"hpo_id": "HP:0000822", "name": "Hypertension", "definition": "Persistently elevated arterial blood pressure", "synonyms": ["High blood pressure", "Elevated BP"], "ic_score": 2.2},
    {"hpo_id": "HP:0001297", "name": "Stroke", "definition": "Acute neurological deficit due to cerebrovascular event", "synonyms": ["Cerebrovascular accident", "CVA", "Brain infarction"], "ic_score": 2.9},
    {"hpo_id": "HP:0002170", "name": "Intracranial hemorrhage", "definition": "Bleeding within the cranium", "synonyms": ["Brain hemorrhage", "ICH"], "ic_score": 3.2},
    {"hpo_id": "HP:0001000", "name": "Abnormality of skin pigmentation", "definition": "Abnormal skin color", "synonyms": ["Pigmentation abnormality", "Dyspigmentation"], "ic_score": 2.0},
    {"hpo_id": "HP:0001030", "name": "Fragile skin", "definition": "Skin that tears or bruises easily", "synonyms": ["Skin fragility", "Tissue fragility"], "ic_score": 3.3},
    {"hpo_id": "HP:0000958", "name": "Dry skin", "definition": "Abnormally dry skin surface", "synonyms": ["Xerosis", "Xerosis cutis"], "ic_score": 2.6},
    {"hpo_id": "HP:0001596", "name": "Alopecia", "definition": "Loss of hair from the head or body", "synonyms": ["Hair loss", "Baldness"], "ic_score": 3.0},
    {"hpo_id": "HP:0001392", "name": "Abnormality of the liver", "definition": "Structural or functional liver abnormality", "synonyms": ["Liver abnormality", "Hepatic abnormality"], "ic_score": 1.8},
    {"hpo_id": "HP:0002240", "name": "Hepatomegaly with fibrosis", "definition": "Enlarged liver with scarring", "synonyms": ["Liver enlargement with fibrosis"], "ic_score": 3.2},
    {"hpo_id": "HP:0001541", "name": "Ascites", "definition": "Accumulation of fluid in the peritoneal cavity", "synonyms": ["Peritoneal effusion", "Abdominal fluid"], "ic_score": 3.1},
    {"hpo_id": "HP:0002829", "name": "Arthralgia", "definition": "Joint pain", "synonyms": ["Joint pain", "Articular pain"], "ic_score": 2.3},
    {"hpo_id": "HP:0001376", "name": "Limitation of joint mobility", "definition": "Reduced range of motion of joints", "synonyms": ["Joint stiffness", "Restricted joint movement"], "ic_score": 2.7},
    {"hpo_id": "HP:0002808", "name": "Kyphosis", "definition": "Excessive curvature of the thoracic spine", "synonyms": ["Hunchback", "Rounded back"], "ic_score": 3.0},
    {"hpo_id": "HP:0002650", "name": "Scoliosis (structural)", "definition": "Lateral curvature of the spine exceeding 10 degrees", "synonyms": ["Spinal curvature", "Curved spine"], "ic_score": 2.7},
    {"hpo_id": "HP:0100491", "name": "Abnormality of lower limb joint", "definition": "Abnormality of a joint of the lower extremity", "synonyms": ["Lower limb joint defect"], "ic_score": 2.5},
    {"hpo_id": "HP:0001762", "name": "Talipes equinovarus", "definition": "Congenital foot deformity", "synonyms": ["Clubfoot", "Club foot"], "ic_score": 3.5},
    {"hpo_id": "HP:0003011", "name": "Abnormality of skeletal muscle", "definition": "Structural or functional muscle abnormality", "synonyms": ["Muscle abnormality", "Myopathic changes"], "ic_score": 1.9},
    {"hpo_id": "HP:0001626", "name": "Abnormality of the cardiovascular system", "definition": "Cardiovascular structural or functional abnormality", "synonyms": ["Cardiovascular defect", "Heart abnormality"], "ic_score": 1.4},
    {"hpo_id": "HP:0010741", "name": "Pedal edema", "definition": "Swelling of the feet", "synonyms": ["Foot swelling", "Lower extremity edema"], "ic_score": 2.8},
    {"hpo_id": "HP:0002017", "name": "Nausea and vomiting", "definition": "Feeling of sickness with or without emesis", "synonyms": ["Emesis", "Vomiting", "Nausea"], "ic_score": 2.1},
    {"hpo_id": "HP:0002014", "name": "Diarrhea", "definition": "Abnormally frequent and liquid bowel movements", "synonyms": ["Loose stools", "Watery stool"], "ic_score": 2.3},
    {"hpo_id": "HP:0002019", "name": "Constipation", "definition": "Infrequent or difficult passage of stool", "synonyms": ["Obstipation"], "ic_score": 2.4},
]


# ===================================================================
# HPO PARSER IMPLEMENTATION
# ===================================================================


class HPOParser(BaseIngestParser):
    """Parse HPO ontology terms for phenotype-driven rare disease diagnosis.

    In offline/seed mode, returns the curated HPO_TERMS list.

    Usage::

        parser = HPOParser()
        records, stats = parser.run()
    """

    def __init__(
        self,
        collection_manager: Any = None,
        embedder: Any = None,
    ) -> None:
        super().__init__(
            source_name="hpo",
            collection_manager=collection_manager,
            embedder=embedder,
        )

    def fetch(self, **kwargs: Any) -> List[Dict[str, Any]]:
        """Fetch HPO term data.

        Returns the curated HPO_TERMS list (seed mode).

        Returns:
            List of raw HPO term dictionaries.
        """
        self.logger.info("Using curated HPO seed data (%d terms)", len(HPO_TERMS))
        return list(HPO_TERMS)

    def parse(self, raw_data: List[Dict[str, Any]]) -> List[IngestRecord]:
        """Parse raw HPO term data into IngestRecord objects.

        Args:
            raw_data: List of HPO term dictionaries.

        Returns:
            List of IngestRecord objects.
        """
        records: List[IngestRecord] = []

        for entry in raw_data:
            hpo_id = entry.get("hpo_id", "")
            name = entry.get("name", "")
            definition = entry.get("definition", "")
            synonyms = entry.get("synonyms", [])
            ic_score = entry.get("ic_score", 0.0)

            synonyms_str = ", ".join(synonyms) if synonyms else "none"
            text = (
                f"HPO Phenotype: {name} ({hpo_id}). "
                f"Definition: {definition}. "
                f"Synonyms: {synonyms_str}. "
                f"Information Content (IC): {ic_score:.1f}."
            )

            record = IngestRecord(
                text=text,
                metadata={
                    "hpo_id": hpo_id,
                    "name": name,
                    "definition": definition,
                    "synonyms": synonyms,
                    "ic_score": ic_score,
                    "source_db": "HPO",
                },
                collection_name="rd_phenotypes",
                record_id=f"HPO_{hpo_id}",
                source="hpo",
            )
            records.append(record)

        return records

    def validate_record(self, record: IngestRecord) -> bool:
        """Validate an HPO IngestRecord.

        Requirements:
            - text must be non-empty
            - must have hpo_id in metadata
            - must have name in metadata

        Args:
            record: The record to validate.

        Returns:
            True if the record passes all validation checks.
        """
        if not record.text or not record.text.strip():
            return False

        meta = record.metadata
        if not meta.get("hpo_id"):
            return False
        if not meta.get("name"):
            return False

        return True


def get_hpo_term_count() -> int:
    """Return the number of curated HPO terms."""
    return len(HPO_TERMS)


def get_hpo_ids() -> List[str]:
    """Return a list of all HPO IDs from the seed data."""
    return [t["hpo_id"] for t in HPO_TERMS]
