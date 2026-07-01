#!/usr/bin/env python3
"""Generate biomarker_critical_values.json and biomarker_discordance_rules.json
reference files for the precision biomarker agent."""

import json
import os

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data', 'reference')

# ---------------------------------------------------------------------------
# 1. Critical Values (21 entries)
# ---------------------------------------------------------------------------

CRITICAL_VALUES = [
    {
        "id": "crit_platelet_very_low",
        "biomarker": "Platelet Count",
        "loinc_code": "777-3",
        "critical_low": None,
        "critical_high": 50.0,
        "severity": "critical",
        "escalation_target": "hematology_oncall",
        "clinical_action": "Hold all anticoagulants and antiplatelets. Order stat peripheral smear. Evaluate for DIC, TTP/HUS, or marrow failure. Consider platelet transfusion if active bleeding or planned procedure.",
        "cross_checks": ["INR", "Fibrinogen", "D-Dimer", "LDH", "Haptoglobin"],
        "text_chunk": "A platelet count below 50 x10^9/L represents a critical thrombocytopenia requiring immediate evaluation. Spontaneous bleeding risk increases significantly below this threshold, and surgical hemostasis cannot be assured. Common etiologies include DIC, TTP/HUS, heparin-induced thrombocytopenia, bone marrow failure, and massive transfusion dilution. Cross-check with coagulation studies and hemolysis markers to differentiate consumptive from production-failure causes."
    },
    {
        "id": "crit_platelet_low",
        "biomarker": "Platelet Count",
        "loinc_code": "777-3",
        "critical_low": None,
        "critical_high": 100.0,
        "severity": "warning",
        "escalation_target": "ordering_provider",
        "clinical_action": "Review medication list for offending agents. Repeat CBC to confirm. Evaluate trend. Consider hematology consult if declining trend.",
        "cross_checks": ["INR", "Fibrinogen", "Mean Platelet Volume"],
        "text_chunk": "A platelet count below 100 x10^9/L indicates mild-to-moderate thrombocytopenia warranting monitoring and medication review. While spontaneous bleeding is uncommon at this level, a declining trend may signal evolving pathology such as immune thrombocytopenia, drug-induced marrow suppression, or early DIC. Mean platelet volume can help distinguish destruction (high MPV) from underproduction (low MPV). Serial monitoring and trend analysis guide urgency of intervention."
    },
    {
        "id": "crit_glucose_high",
        "biomarker": "Glucose",
        "loinc_code": "2345-7",
        "critical_low": None,
        "critical_high": 250.0,
        "severity": "urgent",
        "escalation_target": "endocrinology_oncall",
        "clinical_action": "Check serum ketones and arterial blood gas. Evaluate for DKA or HHS. Initiate IV insulin protocol if confirmed. Monitor electrolytes every 2 hours.",
        "cross_checks": ["HbA1c", "Serum Ketones", "Bicarbonate", "Anion Gap", "Serum Osmolality"],
        "text_chunk": "A blood glucose exceeding 250 mg/dL is an urgent finding that may indicate diabetic ketoacidosis (DKA) or hyperosmolar hyperglycemic state (HHS). Immediate evaluation should include serum ketones, arterial blood gas, and electrolyte panel to classify the hyperglycemic emergency. Concurrent HbA1c helps distinguish new-onset diabetes from chronic poor control. Prompt insulin therapy and aggressive fluid resuscitation reduce morbidity and mortality."
    },
    {
        "id": "crit_glucose_low",
        "biomarker": "Glucose",
        "loinc_code": "2345-7",
        "critical_low": 50.0,
        "critical_high": None,
        "severity": "critical",
        "escalation_target": "rapid_response_team",
        "clinical_action": "Administer IV dextrose (D50) immediately. Check insulin and C-peptide levels before correction if etiology unclear. Monitor glucose every 15 minutes until stable above 100 mg/dL.",
        "cross_checks": ["Insulin", "C-Peptide", "Cortisol", "Beta-Hydroxybutyrate"],
        "text_chunk": "A blood glucose below 50 mg/dL is a critical hypoglycemic emergency with imminent risk of seizures, loss of consciousness, and permanent neurological injury. Whipple's triad (symptoms, documented low glucose, resolution with treatment) confirms clinical significance. Etiologies include insulin overdose, sulfonylurea use, insulinoma, adrenal insufficiency, and sepsis. Paired insulin and C-peptide drawn before glucose correction help determine the underlying cause."
    },
    {
        "id": "crit_potassium_high",
        "biomarker": "Potassium",
        "loinc_code": "2823-3",
        "critical_low": None,
        "critical_high": 6.0,
        "severity": "critical",
        "escalation_target": "rapid_response_team",
        "clinical_action": "Obtain stat 12-lead ECG. Administer IV calcium gluconate for cardiac membrane stabilization. Initiate insulin-glucose infusion and nebulized albuterol for potassium shifting. Evaluate for hemolysis artifact.",
        "cross_checks": ["Creatinine", "eGFR", "Bicarbonate", "ECG"],
        "text_chunk": "A serum potassium above 6.0 mEq/L is a life-threatening hyperkalemia with risk of fatal cardiac arrhythmia including ventricular fibrillation and asystole. ECG changes (peaked T waves, widened QRS, sine wave pattern) indicate myocardial toxicity requiring emergent treatment. Common causes include acute kidney injury, potassium-sparing diuretics, ACE inhibitors, rhabdomyolysis, and tumor lysis syndrome. Pseudohyperkalemia from hemolyzed specimens must be excluded before initiating aggressive therapy."
    },
    {
        "id": "crit_potassium_low",
        "biomarker": "Potassium",
        "loinc_code": "2823-3",
        "critical_low": 3.0,
        "critical_high": None,
        "severity": "urgent",
        "escalation_target": "ordering_provider",
        "clinical_action": "Obtain stat ECG. Begin IV potassium replacement with cardiac monitoring. Check magnesium level concurrently as hypomagnesemia impairs potassium repletion.",
        "cross_checks": ["Magnesium", "Bicarbonate", "Digoxin Level", "ECG"],
        "text_chunk": "A serum potassium below 3.0 mEq/L is a significant hypokalemia that predisposes to life-threatening cardiac arrhythmias, particularly in patients on digoxin. ECG findings may include U waves, ST depression, and T wave flattening. Concurrent hypomagnesemia must be corrected as it prevents effective potassium repletion. Common causes include diuretic therapy, vomiting, diarrhea, and renal tubular acidosis. IV replacement with cardiac monitoring is required at this severity."
    },
    {
        "id": "crit_inr_high",
        "biomarker": "INR",
        "loinc_code": "6301-6",
        "critical_low": None,
        "critical_high": 4.0,
        "severity": "critical",
        "escalation_target": "hematology_oncall",
        "clinical_action": "Hold warfarin. Assess for active bleeding. Consider vitamin K (IV or oral) based on bleeding risk. For life-threatening bleeding, administer 4-factor PCC or FFP. Recheck INR in 6 hours.",
        "cross_checks": ["Platelet Count", "Hemoglobin", "Fibrinogen", "Liver Function Tests"],
        "text_chunk": "An INR above 4.0 represents a critical coagulopathy with markedly elevated bleeding risk, particularly intracranial and gastrointestinal hemorrhage. In warfarin-treated patients, this indicates supratherapeutic anticoagulation requiring immediate dose adjustment or reversal. Hepatic dysfunction, drug interactions (azole antifungals, amiodarone, fluoroquinolones), and vitamin K deficiency are common contributing factors. Cross-referencing platelet count and fibrinogen helps assess overall hemostatic capacity and guides transfusion decisions."
    },
    {
        "id": "crit_sodium_high",
        "biomarker": "Sodium",
        "loinc_code": "2951-2",
        "critical_low": None,
        "critical_high": 150.0,
        "severity": "urgent",
        "escalation_target": "nephrology_oncall",
        "clinical_action": "Assess volume status and urine osmolality. Calculate free water deficit. Correct sodium no faster than 10 mEq/L per 24 hours to avoid cerebral edema. Evaluate for diabetes insipidus.",
        "cross_checks": ["Urine Osmolality", "Urine Sodium", "Serum Osmolality", "ADH Level"],
        "text_chunk": "A serum sodium above 150 mEq/L indicates significant hypernatremia, most commonly from free water deficit due to inadequate intake, diabetes insipidus, or excessive insensible losses. Rapid correction risks cerebral edema, so sodium should be lowered no faster than 10 mEq/L per 24 hours. Urine osmolality distinguishes central diabetes insipidus (dilute urine) from dehydration (concentrated urine). Elderly and critically ill patients with impaired thirst mechanisms are at highest risk."
    },
    {
        "id": "crit_sodium_low",
        "biomarker": "Sodium",
        "loinc_code": "2951-2",
        "critical_low": 125.0,
        "critical_high": None,
        "severity": "critical",
        "escalation_target": "rapid_response_team",
        "clinical_action": "Assess neurological status. For symptomatic hyponatremia, administer 3% hypertonic saline (100 mL bolus over 10 min, may repeat twice). Correct no faster than 8 mEq/L in 24 hours to avoid osmotic demyelination syndrome.",
        "cross_checks": ["Urine Osmolality", "Urine Sodium", "Serum Osmolality", "TSH", "Cortisol"],
        "text_chunk": "A serum sodium below 125 mEq/L is a critical hyponatremia with risk of cerebral edema, seizures, and respiratory arrest. Acute symptomatic hyponatremia requires emergent treatment with hypertonic saline, but correction must not exceed 8 mEq/L in 24 hours to prevent osmotic demyelination syndrome (central pontine myelinolysis). Etiologies include SIADH, heart failure, cirrhosis, thiazide diuretics, and adrenal insufficiency. Urine studies (osmolality and sodium) are essential for determining the underlying mechanism."
    },
    {
        "id": "crit_hemoglobin_low",
        "biomarker": "Hemoglobin",
        "loinc_code": "718-7",
        "critical_low": 7.0,
        "critical_high": None,
        "severity": "critical",
        "escalation_target": "rapid_response_team",
        "clinical_action": "Type and crossmatch for 2 units pRBCs. Assess hemodynamic stability. Identify bleeding source. Activate massive transfusion protocol if ongoing hemorrhage with hemodynamic instability.",
        "cross_checks": ["Reticulocyte Count", "Haptoglobin", "LDH", "MCV", "Iron Studies"],
        "text_chunk": "A hemoglobin below 7 g/dL is a critical anemia requiring transfusion consideration per restrictive transfusion guidelines, and mandatory transfusion in the setting of active hemorrhage or hemodynamic instability. Reticulocyte count distinguishes hypoproliferative (marrow failure, iron/B12 deficiency) from hyperproliferative (hemolysis, acute blood loss) etiologies. Hemolysis workup (haptoglobin, LDH, indirect bilirubin) should be performed if reticulocyte count is elevated. Rapid clinical assessment of perfusion and oxygenation guides transfusion urgency."
    },
    {
        "id": "crit_calcium_high",
        "biomarker": "Calcium (Total)",
        "loinc_code": "17861-6",
        "critical_low": None,
        "critical_high": 12.0,
        "severity": "urgent",
        "escalation_target": "endocrinology_oncall",
        "clinical_action": "Initiate aggressive IV normal saline hydration (200-300 mL/hr). Measure PTH, PTHrP, and 1,25-dihydroxyvitamin D. Consider IV bisphosphonate (zoledronic acid) if malignancy-related. Monitor ECG for shortened QT.",
        "cross_checks": ["PTH", "PTHrP", "Phosphorus", "Albumin", "Ionized Calcium", "Vitamin D"],
        "text_chunk": "A total serum calcium above 12 mg/dL constitutes significant hypercalcemia with risk of cardiac arrhythmias, renal injury, and altered mental status. The two most common causes are primary hyperparathyroidism (elevated PTH) and malignancy (suppressed PTH with elevated PTHrP or 1,25-vitamin D). Ionized calcium should be checked to confirm true hypercalcemia, as total calcium is affected by albumin levels. Aggressive IV hydration is the first-line treatment regardless of etiology."
    },
    {
        "id": "crit_calcium_low",
        "biomarker": "Calcium (Total)",
        "loinc_code": "17861-6",
        "critical_low": 7.0,
        "critical_high": None,
        "severity": "critical",
        "escalation_target": "rapid_response_team",
        "clinical_action": "Obtain stat ionized calcium. If confirmed, administer IV calcium gluconate with cardiac monitoring. Check magnesium and replete if low. Evaluate PTH and vitamin D levels.",
        "cross_checks": ["Ionized Calcium", "Magnesium", "Phosphorus", "PTH", "Vitamin D", "Albumin"],
        "text_chunk": "A total serum calcium below 7 mg/dL is a critical hypocalcemia with risk of tetany, laryngospasm, seizures, and cardiac arrest from prolonged QT interval. Ionized calcium confirms true hypocalcemia versus pseudohypocalcemia from hypoalbuminemia. Common causes include post-thyroidectomy hypoparathyroidism, severe vitamin D deficiency, magnesium depletion, and acute pancreatitis. IV calcium gluconate with cardiac monitoring is required for symptomatic or severe cases."
    },
    {
        "id": "crit_troponin_i_high",
        "biomarker": "Troponin I (High-Sensitivity)",
        "loinc_code": "89579-7",
        "critical_low": None,
        "critical_high": 26.0,
        "severity": "critical",
        "escalation_target": "cardiology_oncall",
        "clinical_action": "Activate STEMI/NSTEMI protocol. Obtain stat 12-lead ECG. Administer aspirin, heparin, and consider P2Y12 inhibitor. Arrange emergent cardiac catheterization if ECG shows ST changes.",
        "cross_checks": ["BNP", "CK-MB", "ECG", "D-Dimer", "Creatinine"],
        "text_chunk": "A high-sensitivity troponin I above 26 ng/L (above the 99th percentile upper reference limit) indicates myocardial injury and is the cornerstone biomarker for diagnosing acute myocardial infarction when combined with a rise-and-fall pattern. Serial measurements at 0 and 3 hours establish the kinetic pattern. Non-ACS causes of troponin elevation include myocarditis, pulmonary embolism, sepsis, heart failure exacerbation, and renal failure. ECG correlation and clinical context are essential for distinguishing type 1 MI from type 2 MI and non-ischemic myocardial injury."
    },
    {
        "id": "crit_wbc_high",
        "biomarker": "WBC Count",
        "loinc_code": "6690-2",
        "critical_low": None,
        "critical_high": 30.0,
        "severity": "urgent",
        "escalation_target": "hematology_oncall",
        "clinical_action": "Order stat peripheral blood smear with manual differential. Evaluate for leukemia, leukemoid reaction, or severe infection. If blasts present, arrange emergent hematology consultation.",
        "cross_checks": ["Differential Count", "Peripheral Smear", "LDH", "Uric Acid", "CRP"],
        "text_chunk": "A WBC count above 30 x10^9/L is an urgent finding requiring peripheral smear review to distinguish reactive leukocytosis from hematologic malignancy. Leukemoid reactions (WBC >50) can occur with severe infections, burns, and corticosteroid use. Presence of blasts or immature cells on smear mandates emergent hematology evaluation for acute leukemia. LDH and uric acid elevation in the context of high WBC raises concern for tumor lysis syndrome."
    },
    {
        "id": "crit_wbc_low",
        "biomarker": "WBC Count",
        "loinc_code": "6690-2",
        "critical_low": 2.0,
        "critical_high": None,
        "severity": "critical",
        "escalation_target": "hematology_oncall",
        "clinical_action": "Obtain absolute neutrophil count (ANC). If ANC <500, initiate neutropenic fever precautions. Reverse isolation if indicated. Evaluate for drug-induced agranulocytosis, aplastic anemia, or marrow infiltration.",
        "cross_checks": ["ANC", "Differential Count", "Reticulocyte Count", "Peripheral Smear"],
        "text_chunk": "A WBC count below 2 x10^9/L indicates severe leukopenia with high risk of life-threatening infections, particularly if the absolute neutrophil count (ANC) falls below 500 cells/uL (severe neutropenia). Common causes include chemotherapy-induced myelosuppression, drug-induced agranulocytosis (clozapine, methimazole, carbimazole), aplastic anemia, and viral bone marrow suppression. Neutropenic precautions and empiric broad-spectrum antibiotics at the first sign of fever are essential to prevent sepsis-related mortality."
    },
    {
        "id": "crit_creatinine_high",
        "biomarker": "Creatinine",
        "loinc_code": "2160-0",
        "critical_low": None,
        "critical_high": 4.0,
        "severity": "urgent",
        "escalation_target": "nephrology_oncall",
        "clinical_action": "Assess volume status and urine output. Obtain renal ultrasound to exclude obstruction. Review nephrotoxic medications. Check potassium, bicarbonate, and phosphorus. Evaluate for emergent dialysis indications.",
        "cross_checks": ["eGFR", "BUN", "Potassium", "Bicarbonate", "Phosphorus", "Urinalysis"],
        "text_chunk": "A serum creatinine above 4.0 mg/dL indicates severe renal impairment (estimated GFR typically <15 mL/min) and requires urgent nephrology evaluation. Emergent dialysis indications include refractory hyperkalemia, severe metabolic acidosis, volume overload unresponsive to diuretics, and uremic encephalopathy or pericarditis. Renal ultrasound is essential to exclude obstructive uropathy. All nephrotoxic medications (NSAIDs, aminoglycosides, contrast agents) must be immediately discontinued."
    },
    {
        "id": "crit_bilirubin_total_high",
        "biomarker": "Total Bilirubin",
        "loinc_code": "1975-2",
        "critical_low": None,
        "critical_high": 10.0,
        "severity": "urgent",
        "escalation_target": "gastroenterology_oncall",
        "clinical_action": "Fractionate to direct and indirect bilirubin. Obtain liver function panel (AST, ALT, ALP, GGT, albumin). Order right upper quadrant ultrasound. Evaluate for biliary obstruction, acute hepatitis, or hemolytic crisis.",
        "cross_checks": ["Direct Bilirubin", "AST", "ALT", "ALP", "GGT", "Albumin", "INR"],
        "text_chunk": "A total bilirubin above 10 mg/dL indicates severe hyperbilirubinemia requiring urgent evaluation for biliary obstruction, acute hepatic failure, or massive hemolysis. Fractionation into direct (conjugated) and indirect (unconjugated) bilirubin guides the differential diagnosis. Predominantly direct hyperbilirubinemia suggests biliary obstruction or hepatocellular injury, while indirect predominance indicates hemolysis or Gilbert syndrome exacerbation. Concurrent evaluation of liver synthetic function (INR, albumin) assesses severity of hepatic impairment."
    },
    {
        "id": "crit_egfr_low",
        "biomarker": "eGFR (CKD-EPI)",
        "loinc_code": "98979-8",
        "critical_low": 15.0,
        "critical_high": None,
        "severity": "critical",
        "escalation_target": "nephrology_oncall",
        "clinical_action": "Urgent nephrology referral for dialysis planning or transplant evaluation. Review and adjust all renally cleared medications. Assess for uremic symptoms. Evaluate vascular access options.",
        "cross_checks": ["Creatinine", "BUN", "Potassium", "Phosphorus", "PTH", "Hemoglobin"],
        "text_chunk": "An eGFR below 15 mL/min/1.73m2 represents Stage 5 chronic kidney disease (kidney failure) and necessitates urgent nephrology involvement for renal replacement therapy planning. Patients at this stage are at high risk for uremic complications including encephalopathy, pericarditis, coagulopathy, and severe electrolyte derangements. All renally cleared medications require dose adjustment or discontinuation. Concurrent anemia (erythropoietin deficiency) and secondary hyperparathyroidism are expected complications requiring management."
    },
    {
        "id": "crit_free_t4_high",
        "biomarker": "Free T4",
        "loinc_code": "3024-7",
        "critical_low": None,
        "critical_high": 5.0,
        "severity": "critical",
        "escalation_target": "endocrinology_oncall",
        "clinical_action": "Evaluate for thyroid storm (Burch-Wartofsky score). Initiate PTU or methimazole, beta-blocker, and stress-dose corticosteroids if thyroid storm suspected. Check TSH, Free T3, and thyroid antibodies.",
        "cross_checks": ["TSH", "Free T3", "Anti-TPO", "TSI", "Thyroglobulin"],
        "text_chunk": "A free T4 above 5 ng/dL represents severe thyrotoxicosis with risk of thyroid storm, a life-threatening endocrine emergency with mortality exceeding 20% if untreated. Clinical features include hyperthermia, tachycardia, altered mental status, and heart failure. The Burch-Wartofsky Point Scale quantifies thyroid storm severity. Treatment requires multi-modal therapy: thionamides to block new hormone synthesis, iodine to block release, beta-blockers for adrenergic symptoms, and corticosteroids to inhibit T4-to-T3 conversion."
    },
    {
        "id": "crit_tsh_high",
        "biomarker": "TSH",
        "loinc_code": "11580-8",
        "critical_low": None,
        "critical_high": 50.0,
        "severity": "urgent",
        "escalation_target": "endocrinology_oncall",
        "clinical_action": "Check Free T4 and Free T3. Evaluate for myxedema coma if altered mental status, hypothermia, or bradycardia present. Initiate IV levothyroxine if myxedema coma suspected. Start oral levothyroxine replacement.",
        "cross_checks": ["Free T4", "Free T3", "Anti-TPO", "Cortisol"],
        "text_chunk": "A TSH above 50 mIU/L indicates severe primary hypothyroidism with profound thyroid hormone deficiency. If accompanied by altered mental status, hypothermia, and bradycardia, myxedema coma must be considered -- a medical emergency with mortality rates of 30-60%. Cortisol should be checked before initiating thyroid replacement, as coexistent adrenal insufficiency (Schmidt syndrome) requires concurrent glucocorticoid therapy. Most cases result from autoimmune thyroiditis (Hashimoto disease) confirmed by anti-TPO antibodies."
    },
    {
        "id": "crit_lactate_high",
        "biomarker": "Lactate",
        "loinc_code": "2524-7",
        "critical_low": None,
        "critical_high": 4.0,
        "severity": "critical",
        "escalation_target": "critical_care_oncall",
        "clinical_action": "Assess tissue perfusion and hemodynamic status. Initiate aggressive fluid resuscitation if sepsis suspected. Identify and treat underlying cause of tissue hypoxia. Serial lactate measurements every 2-4 hours to guide resuscitation adequacy.",
        "cross_checks": ["Bicarbonate", "Anion Gap", "Base Excess", "Procalcitonin", "Blood Cultures"],
        "text_chunk": "A serum lactate above 4.0 mmol/L indicates severe tissue hypoperfusion or impaired lactate clearance and is associated with significantly increased mortality, particularly in sepsis where it defines septic shock per the Sepsis-3 criteria. Type A lactic acidosis results from tissue hypoxia (shock, severe anemia, cardiac arrest), while Type B occurs without overt hypoxia (hepatic failure, malignancy, medication effects such as metformin or linezolid). Serial lactate clearance (>10% reduction per 2 hours) is a validated resuscitation target and predictor of survival. Concurrent metabolic acidosis with elevated anion gap supports the clinical significance of the lactate elevation."
    }
]

# ---------------------------------------------------------------------------
# 2. Discordance Rules (12 entries)
# ---------------------------------------------------------------------------

DISCORDANCE_RULES = [
    {
        "id": "disc_ferritin_low_tsat_high",
        "name": "Contradictory Iron Studies",
        "biomarker_a": "Ferritin",
        "biomarker_b": "Transferrin Saturation (TSAT)",
        "condition": "Ferritin LOW and TSAT HIGH",
        "expected_relationship": "Ferritin and TSAT typically move in the same direction; low ferritin with high TSAT is physiologically contradictory.",
        "differential_diagnosis": [
            "Specimen mislabeling or pre-analytical error",
            "Acute hepatocellular necrosis releasing ferritin stores with concurrent iron release",
            "Recent parenteral iron infusion with incomplete incorporation",
            "Hereditary hemochromatosis with concurrent inflammation suppressing ferritin"
        ],
        "agent_handoff": ["lab_quality_agent", "iron_metabolism_agent"],
        "priority": "high",
        "text_chunk": "Ferritin and transferrin saturation normally trend concordantly: both low in iron deficiency, both elevated in iron overload. When ferritin is low but TSAT is elevated, this discordance is physiologically contradictory and should raise suspicion for specimen error or a complex pathophysiological state. Acute hepatocellular injury can release stored iron, temporarily elevating TSAT while inflammatory cytokines suppress ferritin. Recent parenteral iron administration may also produce this pattern before reticuloendothelial incorporation is complete. This combination warrants repeat specimen collection and clinical correlation before diagnostic conclusions."
    },
    {
        "id": "disc_ft3_high_tsh_normal",
        "name": "Thyroid Hormone-TSH Discordance",
        "biomarker_a": "Free T3",
        "biomarker_b": "TSH",
        "condition": "Free T3 HIGH and TSH NORMAL",
        "expected_relationship": "Elevated Free T3 should suppress TSH via negative feedback on the hypothalamic-pituitary-thyroid axis.",
        "differential_diagnosis": [
            "TSH-secreting pituitary adenoma (thyrotropinoma)",
            "Thyroid hormone resistance syndrome (THRB mutation)",
            "Assay interference (heterophilic antibodies, biotin)",
            "T3-predominant early Graves disease with delayed TSH suppression",
            "Familial dysalbuminemic hyperthyroxinemia (FDH)"
        ],
        "agent_handoff": ["endocrine_agent", "pituitary_agent"],
        "priority": "high",
        "text_chunk": "Elevated Free T3 with a non-suppressed (normal) TSH violates the expected negative feedback relationship of the hypothalamic-pituitary-thyroid axis. In normal physiology, elevated thyroid hormones suppress TSH to below the reference range. This discordance pattern is characteristic of TSH-secreting pituitary adenomas (central hyperthyroidism) or thyroid hormone resistance syndrome caused by THRB gene mutations. Assay interference from heterophilic antibodies or biotin supplementation must be excluded by running the sample on an alternative platform. A pituitary MRI and alpha-subunit measurement are indicated if TSH-secreting adenoma is suspected. Thyroid hormone resistance is confirmed by genetic testing for THRB mutations."
    },
    {
        "id": "disc_ggt_isolated_high",
        "name": "Isolated GGT Elevation",
        "biomarker_a": "GGT",
        "biomarker_b": "AST/ALT/ALP",
        "condition": "GGT HIGH and AST, ALT, ALP all NORMAL",
        "expected_relationship": "GGT elevation from hepatobiliary disease is typically accompanied by elevation of at least one other liver enzyme.",
        "differential_diagnosis": [
            "Alcohol use (GGT is the most sensitive marker of alcohol consumption)",
            "Enzyme-inducing medications (phenytoin, barbiturates, carbamazepine)",
            "Non-alcoholic fatty liver disease (early stage)",
            "Pancreatic disease",
            "Metabolic syndrome and obesity"
        ],
        "agent_handoff": ["hepatology_agent", "substance_use_agent"],
        "priority": "medium",
        "text_chunk": "Isolated GGT elevation without concomitant AST, ALT, or ALP abnormalities most commonly reflects alcohol consumption, as GGT is the most sensitive (though not specific) biomarker of alcohol use with a half-life of 14-26 days. Enzyme-inducing medications including phenytoin, barbiturates, and carbamazepine upregulate GGT without causing hepatocellular injury. Early non-alcoholic fatty liver disease may present with isolated GGT before transaminase elevation becomes apparent. Pancreatic pathology (chronic pancreatitis, pancreatic carcinoma) can elevate GGT independently of liver enzymes. Clinical correlation with medication history, alcohol screening questionnaires (AUDIT-C), and hepatic ultrasound guides further evaluation."
    },
    {
        "id": "disc_fibrinogen_high_crp_normal",
        "name": "Acute Phase Reactant Discordance",
        "biomarker_a": "Fibrinogen",
        "biomarker_b": "hs-CRP",
        "condition": "Fibrinogen HIGH and hs-CRP NORMAL",
        "expected_relationship": "Both are acute phase reactants that typically elevate concurrently during inflammatory states.",
        "differential_diagnosis": [
            "Estrogen effect (oral contraceptives, pregnancy, hormone replacement therapy)",
            "Nephrotic syndrome (hepatic overproduction of fibrinogen)",
            "Smoking (isolated fibrinogen elevation)",
            "Residual fibrinogen elevation after resolved inflammation (slower clearance than CRP)",
            "Genetic polymorphisms in FGB gene causing constitutively elevated fibrinogen"
        ],
        "agent_handoff": ["coagulation_agent", "inflammation_agent"],
        "priority": "medium",
        "text_chunk": "Fibrinogen and hs-CRP are both hepatic acute phase reactants that typically rise in parallel during systemic inflammation. When fibrinogen is elevated but hs-CRP remains normal, non-inflammatory causes of fibrinogen elevation should be considered. Estrogen states (pregnancy, oral contraceptives, HRT) potently stimulate hepatic fibrinogen synthesis without triggering CRP production. Nephrotic syndrome drives compensatory hepatic protein synthesis including fibrinogen. Fibrinogen has a longer half-life (4-5 days) than CRP (19 hours), so discordance may represent resolving inflammation where CRP has normalized but fibrinogen remains elevated. FGB gene polymorphisms can cause constitutively elevated fibrinogen levels without inflammatory significance."
    },
    {
        "id": "disc_cystatin_c_high_egfr_normal",
        "name": "Renal Filtration Marker Discordance",
        "biomarker_a": "Cystatin C",
        "biomarker_b": "eGFR (Creatinine-based)",
        "condition": "Cystatin C HIGH and eGFR (creatinine-based) NORMAL",
        "expected_relationship": "Both markers should reflect the same glomerular filtration rate; discordance suggests one marker is confounded.",
        "differential_diagnosis": [
            "Creatinine underestimating renal impairment due to low muscle mass (sarcopenia, cachexia, amputation)",
            "Early diabetic nephropathy (cystatin C rises before creatinine)",
            "Thyroid dysfunction (hyperthyroidism increases cystatin C independently)",
            "High-dose corticosteroid therapy (elevates cystatin C non-renally)",
            "Obesity (adipose tissue produces cystatin C)"
        ],
        "agent_handoff": ["nephrology_agent", "geriatric_agent"],
        "priority": "high",
        "text_chunk": "Elevated cystatin C with a normal creatinine-based eGFR suggests occult renal impairment masked by low creatinine generation, a scenario common in elderly patients with sarcopenia, amputees, and cachectic patients with chronic illness. Cystatin C is produced at a constant rate by all nucleated cells and is not influenced by muscle mass, diet, or tubular secretion, making it a more reliable filtration marker in these populations. Early diabetic nephropathy may elevate cystatin C before creatinine rises, providing an opportunity for earlier intervention. However, non-renal confounders of cystatin C must be considered: hyperthyroidism, high-dose corticosteroids, and obesity can elevate cystatin C independently of GFR. The CKD-EPI cystatin C or combined creatinine-cystatin C equation should be used for more accurate GFR estimation."
    },
    {
        "id": "disc_adiponectin_high_hdlp_low",
        "name": "Metabolic Discordance",
        "biomarker_a": "Adiponectin",
        "biomarker_b": "HDL Particle Number (HDL-P)",
        "condition": "Adiponectin HIGH and HDL-P LOW",
        "expected_relationship": "Adiponectin normally promotes HDL biogenesis via ABCA1 upregulation; high adiponectin should correlate with high HDL-P.",
        "differential_diagnosis": [
            "Adiponectin resistance (receptor downregulation despite high ligand)",
            "ABCA1 loss-of-function variants impairing HDL biogenesis",
            "Advanced chronic kidney disease (impaired adiponectin clearance with independent HDL suppression)",
            "Anorexia nervosa or severe caloric restriction (elevated adiponectin, depleted lipoproteins)",
            "Heart failure (adiponectin paradox -- elevated levels associated with worse prognosis)"
        ],
        "agent_handoff": ["metabolic_agent", "cardiovascular_risk_agent"],
        "priority": "medium",
        "text_chunk": "Adiponectin normally stimulates hepatic HDL biogenesis through ABCA1 pathway upregulation, so elevated adiponectin with low HDL particle number is metabolically discordant. This pattern may indicate adiponectin resistance, analogous to insulin resistance, where elevated ligand concentrations fail to produce the expected downstream effects. ABCA1 loss-of-function genetic variants (as seen in Tangier disease spectrum) impair HDL formation regardless of adiponectin signaling. In advanced CKD, impaired renal clearance raises adiponectin while uremia independently suppresses HDL production. The adiponectin paradox in heart failure, where higher levels paradoxically predict worse outcomes, may also present with this discordance pattern. Genetic testing for ABCA1 variants and assessment of renal function help clarify the underlying mechanism."
    },
    {
        "id": "disc_hba1c_normal_glucose_critical",
        "name": "Glycemic Discordance",
        "biomarker_a": "HbA1c",
        "biomarker_b": "Fasting Glucose",
        "condition": "HbA1c NORMAL and Fasting Glucose CRITICAL HIGH (>250 mg/dL)",
        "expected_relationship": "HbA1c reflects 2-3 month average glycemia; a critically elevated fasting glucose should be accompanied by elevated HbA1c unless the hyperglycemia is acute.",
        "differential_diagnosis": [
            "New-onset diabetes mellitus (hyperglycemia too recent to affect HbA1c)",
            "Stress hyperglycemia (acute illness, surgery, corticosteroid-induced)",
            "Hemoglobin variants (HbS, HbC, HbE) causing falsely low HbA1c",
            "Hemolytic anemia or recent blood loss reducing HbA1c (shortened RBC lifespan)",
            "Recent blood transfusion diluting glycated hemoglobin"
        ],
        "agent_handoff": ["endocrine_agent", "glycemic_management_agent"],
        "priority": "high",
        "text_chunk": "A normal HbA1c in the setting of critically elevated fasting glucose (>250 mg/dL) represents a clinically important discordance requiring careful interpretation. The most benign explanation is new-onset hyperglycemia too recent to impact the 2-3 month HbA1c average, as seen in stress hyperglycemia from acute illness, surgery, or corticosteroid administration. However, conditions that falsely lower HbA1c must be excluded: hemoglobin variants (HbS, HbC, HbE) interfere with many HbA1c assays, and hemolytic anemias shorten red blood cell lifespan reducing glycation time. Recent blood transfusion can dilute glycated hemoglobin with donor cells. Fructosamine (reflecting 2-3 week glycemic average) can serve as an alternative glycation marker when HbA1c is unreliable. This discordance should never lead to dismissal of the acute hyperglycemia."
    },
    {
        "id": "disc_ldlc_apob_discordant",
        "name": "Lipid Particle Discordance",
        "biomarker_a": "LDL-C",
        "biomarker_b": "ApoB",
        "condition": "LDL-C DISCORDANT with ApoB (low LDL-C with high ApoB or vice versa)",
        "expected_relationship": "LDL-C and ApoB typically correlate, as each LDL particle carries one ApoB molecule. Discordance indicates abnormal particle composition.",
        "differential_diagnosis": [
            "Small dense LDL phenotype (high ApoB, low/normal LDL-C -- each particle carries less cholesterol)",
            "Metabolic syndrome and insulin resistance (drives small dense LDL formation)",
            "Hypertriglyceridemia causing Friedewald LDL-C underestimation",
            "Familial hypobetalipoproteinemia (low ApoB, variable LDL-C)",
            "PCSK9 gain-of-function variants (high ApoB relative to LDL-C)"
        ],
        "agent_handoff": ["cardiovascular_risk_agent", "lipid_agent"],
        "priority": "high",
        "text_chunk": "Discordance between LDL-C and ApoB has significant implications for cardiovascular risk stratification. When ApoB is elevated but LDL-C is normal or low, the patient has a high number of cholesterol-depleted small dense LDL particles, each carrying one ApoB molecule but less cholesterol per particle. This pattern, common in metabolic syndrome and insulin resistance, means LDL-C underestimates true atherogenic particle burden. Mendelian randomization studies confirm that ApoB (particle number) is the primary causal driver of atherosclerosis, not cholesterol content per se. Conversely, high LDL-C with low ApoB indicates fewer but cholesterol-enriched large buoyant particles, which may carry lower risk. When Friedewald-calculated LDL-C is used with triglycerides >150 mg/dL, LDL-C is systematically underestimated, making ApoB the preferred risk metric. Treatment decisions should be guided by ApoB when discordant with LDL-C."
    },
    {
        "id": "disc_platelet_inr_compound",
        "name": "Compound Hemorrhage Risk",
        "biomarker_a": "Platelet Count",
        "biomarker_b": "INR",
        "condition": "Platelet Count CRITICAL LOW (<50) and INR CRITICAL HIGH (>4.0)",
        "expected_relationship": "Independent coagulation defects; their co-occurrence compounds hemorrhage risk multiplicatively.",
        "differential_diagnosis": [
            "Disseminated intravascular coagulation (DIC)",
            "Severe liver failure (impaired thrombopoietin and clotting factor synthesis)",
            "Massive transfusion coagulopathy (dilutional)",
            "Heparin-induced thrombocytopenia with concurrent warfarin overlap",
            "Thrombotic microangiopathy with superimposed anticoagulation"
        ],
        "agent_handoff": ["hematology_agent", "critical_care_agent", "transfusion_agent"],
        "priority": "high",
        "text_chunk": "The simultaneous occurrence of critical thrombocytopenia (platelet count <50 x10^9/L) and critical coagulopathy (INR >4.0) represents a compound hemorrhage emergency where bleeding risk is multiplicatively increased. This dual defect impairs both primary hemostasis (platelet plug formation) and secondary hemostasis (fibrin clot stabilization). Disseminated intravascular coagulation (DIC) is the most dangerous etiology, characterized by consumptive coagulopathy with simultaneous thrombosis and hemorrhage. Severe hepatic failure produces this pattern through combined thrombopoietin deficiency and impaired clotting factor synthesis. Fibrinogen, D-dimer, and peripheral smear for schistocytes are essential for distinguishing DIC from other causes. Management requires simultaneous attention to both defects: platelet transfusion, FFP or PCC for coagulopathy, and treatment of the underlying trigger."
    },
    {
        "id": "disc_tsh_ft4_discordant",
        "name": "TSH-Free T4 Axis Discordance",
        "biomarker_a": "TSH",
        "biomarker_b": "Free T4",
        "condition": "TSH and Free T4 moving in the same direction (both elevated or both suppressed) instead of the expected inverse relationship",
        "expected_relationship": "TSH and Free T4 have an inverse log-linear relationship via hypothalamic-pituitary-thyroid negative feedback.",
        "differential_diagnosis": [
            "TSH-secreting pituitary adenoma (both elevated)",
            "Thyroid hormone resistance syndrome (both elevated, THRB mutation)",
            "Central hypothyroidism (both low -- pituitary or hypothalamic failure)",
            "Non-thyroidal illness / euthyroid sick syndrome (both low in severe illness)",
            "Assay interference (biotin, heterophilic antibodies)",
            "Recent thyroid hormone dose adjustment (TSH lags Free T4 by 6-8 weeks)"
        ],
        "agent_handoff": ["endocrine_agent", "pituitary_agent"],
        "priority": "high",
        "text_chunk": "The TSH-Free T4 relationship is the cornerstone of thyroid function assessment, with a well-characterized inverse log-linear feedback loop. When TSH and Free T4 move in the same direction (both elevated or both suppressed), this discordance indicates either a central axis defect or assay artifact. Both elevated: consider TSH-secreting pituitary adenoma (thyrotropinoma) or thyroid hormone resistance syndrome, both rare but clinically significant diagnoses requiring pituitary MRI and genetic testing respectively. Both suppressed: suggests central (secondary) hypothyroidism from pituitary or hypothalamic disease, or non-thyroidal illness syndrome in acutely ill patients. Biotin supplementation causes predictable assay interference (falsely low TSH, falsely high Free T4 on streptavidin-biotin platforms). A 48-hour biotin washout and repeat testing should precede extensive workup. Recent levothyroxine dose changes cause transient discordance as TSH equilibration lags 6-8 weeks behind Free T4 changes."
    },
    {
        "id": "disc_b12_normal_mma_high",
        "name": "Functional B12 Deficiency",
        "biomarker_a": "Vitamin B12",
        "biomarker_b": "Methylmalonic Acid (MMA)",
        "condition": "Vitamin B12 NORMAL and MMA ELEVATED",
        "expected_relationship": "Normal B12 should provide adequate cofactor for methylmalonyl-CoA mutase, keeping MMA in normal range.",
        "differential_diagnosis": [
            "Functional (cellular) B12 deficiency despite normal serum levels",
            "B12 bound to inactive carriers (haptocorrin) rather than active transcobalamin (holotranscobalamin)",
            "Renal insufficiency (impaired MMA excretion)",
            "Intestinal bacterial overgrowth (propionate metabolism increasing MMA)",
            "Methylmalonic acidemia (inborn error of metabolism, rare in adults)"
        ],
        "agent_handoff": ["hematology_agent", "nutrition_agent"],
        "priority": "medium",
        "text_chunk": "Elevated methylmalonic acid (MMA) with a normal serum B12 level indicates functional (tissue-level) B12 deficiency, a condition increasingly recognized as clinically significant. Serum B12 measures total cobalamin bound to both active transcobalamin II (holoTC, ~20%) and inactive haptocorrin (~80%). A patient may have normal total B12 but insufficient active holoTC, resulting in cellular B12 deficiency detectable only by metabolic markers (MMA, homocysteine). HoloTC (active B12) is a more specific marker and should be measured when this discordance is identified. Renal insufficiency must be excluded as a confounder, since impaired renal excretion elevates MMA independent of B12 status. A therapeutic trial of B12 supplementation with MMA normalization confirms functional deficiency. This discordance is common in elderly patients with subclinical B12 depletion and in patients on metformin or proton pump inhibitors."
    },
    {
        "id": "disc_vitd_normal_pth_high",
        "name": "Secondary Hyperparathyroidism",
        "biomarker_a": "25-Hydroxyvitamin D",
        "biomarker_b": "PTH (Intact)",
        "condition": "25-Hydroxyvitamin D NORMAL and PTH ELEVATED",
        "expected_relationship": "Adequate vitamin D should suppress PTH via calcium-mediated negative feedback; elevated PTH with normal vitamin D suggests autonomous parathyroid activity or calcium sensing defect.",
        "differential_diagnosis": [
            "Primary hyperparathyroidism (autonomous parathyroid adenoma or hyperplasia)",
            "Chronic kidney disease (impaired 1,25-dihydroxyvitamin D production despite normal 25-OH-D)",
            "Calcium malabsorption (celiac disease, bariatric surgery, chronic pancreatitis)",
            "Lithium-induced hyperparathyroidism (altered calcium set point)",
            "Familial hypocalciuric hypercalcemia (CaSR mutation)"
        ],
        "agent_handoff": ["endocrine_agent", "nephrology_agent", "bone_metabolism_agent"],
        "priority": "medium",
        "text_chunk": "Elevated PTH with normal 25-hydroxyvitamin D levels creates a diagnostic challenge that requires systematic evaluation. In primary hyperparathyroidism, one or more parathyroid glands autonomously secrete PTH independent of calcium and vitamin D status, making this the most common cause in outpatient settings. Serum calcium (usually elevated in primary, normal or low in secondary) is the key differentiating test. In chronic kidney disease, although 25-OH vitamin D (the storage form) may be normal, the kidney fails to convert it to active 1,25-dihydroxyvitamin D, leading to calcium malabsorption and compensatory PTH elevation. Post-bariatric surgery patients and those with celiac disease may have calcium malabsorption despite normal vitamin D, stimulating secondary PTH elevation. Lithium therapy alters the calcium-sensing receptor set point, requiring higher calcium levels to suppress PTH. A comprehensive workup includes serum calcium, phosphorus, 1,25-dihydroxyvitamin D, 24-hour urine calcium, and renal function."
    }
]


def main():
    os.makedirs(DATA_DIR, exist_ok=True)

    # Write critical values
    cv_path = os.path.join(DATA_DIR, 'biomarker_critical_values.json')
    with open(cv_path, 'w') as f:
        json.dump(CRITICAL_VALUES, f, indent=2)
    print(f"Wrote {len(CRITICAL_VALUES)} critical values to {cv_path}")

    # Write discordance rules
    dr_path = os.path.join(DATA_DIR, 'biomarker_discordance_rules.json')
    with open(dr_path, 'w') as f:
        json.dump(DISCORDANCE_RULES, f, indent=2)
    print(f"Wrote {len(DISCORDANCE_RULES)} discordance rules to {dr_path}")


if __name__ == '__main__':
    main()
