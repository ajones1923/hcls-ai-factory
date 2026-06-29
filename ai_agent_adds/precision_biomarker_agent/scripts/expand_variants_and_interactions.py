#!/usr/bin/env python3
"""
Expand biomarker_genetic_variants.json and biomarker_drug_interactions.json
with new entries for AJ carrier panel, inflammatory, neurological, bone health,
cardiovascular, and drug-biomarker interaction variants.

Reads existing files, checks for duplicate IDs, appends new entries, and writes back.
"""

import json
import os

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VARIANTS_PATH = os.path.join(BASE_DIR, "data", "reference", "biomarker_genetic_variants.json")
INTERACTIONS_PATH = os.path.join(BASE_DIR, "data", "reference", "biomarker_drug_interactions.json")


def load_json(path):
    with open(path, "r") as f:
        return json.load(f)


def save_json(path, data):
    with open(path, "w") as f:
        json.dump(data, f, indent=2)
    print(f"  Wrote {len(data)} total entries to {path}")


# ---------------------------------------------------------------------------
# Part 1: New Genetic Variants
# ---------------------------------------------------------------------------

NEW_VARIANTS = [
    # --- AJ Carrier Panel Variants (6) ---
    {
        "id": "gv_gba_rs76763715",
        "gene": "GBA",
        "rs_id": "rs76763715",
        "risk_allele": "A",
        "protective_allele": "G",
        "effect_size": "N370S mutation; Gaucher Disease carrier frequency 1/15 in Ashkenazi Jewish population; compound risk with APOE E4 for Parkinson's/Lewy body dementia",
        "mechanism": "GBA encodes glucocerebrosidase, a lysosomal enzyme that degrades glucosylceramide. The N370S (rs76763715) variant reduces enzyme activity, causing glucosylceramide accumulation. Heterozygous carriers typically do not develop Gaucher disease but have a 5-8x increased risk of Parkinson's disease and Lewy body dementia through impaired alpha-synuclein degradation via the autophagy-lysosomal pathway. The risk is compounded in APOE E4 carriers, creating a synergistic neurodegenerative risk profile.",
        "disease_associations": "Gaucher Disease (homozygous),Parkinson's Disease,Lewy Body Dementia,Lysosomal Storage Disorder",
        "text_chunk": "GBA rs76763715 (N370S) is the most common Gaucher Disease mutation, with carrier frequency of 1/15 in the Ashkenazi Jewish population. Heterozygous carriers have 5-8x increased risk of Parkinson's disease and Lewy body dementia via impaired lysosomal alpha-synuclein clearance. This risk is compounded in APOE E4 carriers. Homozygous N370S causes Type 1 (non-neuronopathic) Gaucher disease treatable with enzyme replacement therapy. Carrier screening is recommended for AJ individuals planning families."
    },
    {
        "id": "gv_brca1_rs80357914",
        "gene": "BRCA1",
        "rs_id": "rs80357914",
        "risk_allele": "delAG",
        "protective_allele": "wt",
        "effect_size": "185delAG founder mutation; 60-80% lifetime breast cancer risk; 20-40% ovarian cancer risk; carrier frequency 1/40 in Ashkenazi Jewish population",
        "mechanism": "BRCA1 encodes a tumor suppressor protein essential for homologous recombination DNA double-strand break repair. The 185delAG (rs80357914) frameshift mutation produces a truncated, non-functional protein, abolishing BRCA1-mediated DNA repair. This leads to genomic instability and accumulation of oncogenic mutations. Carriers have markedly elevated lifetime risks of breast (60-80%) and ovarian (20-40%) cancer. BRCA1-deficient tumors are sensitive to platinum chemotherapy and PARP inhibitors (olaparib, niraparib) due to synthetic lethality.",
        "disease_associations": "Hereditary Breast Cancer,Ovarian Cancer,Fallopian Tube Cancer,Peritoneal Cancer",
        "text_chunk": "BRCA1 185delAG (rs80357914) is an Ashkenazi Jewish founder mutation with carrier frequency of 1/40. It abolishes BRCA1 DNA repair function, conferring 60-80% lifetime breast cancer and 20-40% ovarian cancer risk. BRCA1-deficient tumors respond to PARP inhibitors (olaparib, niraparib) and platinum chemotherapy via synthetic lethality. Risk-reducing mastectomy and salpingo-oophorectomy are evidence-based prevention options. Enhanced screening with breast MRI starting at age 25 is recommended for carriers."
    },
    {
        "id": "gv_brca2_rs80359550",
        "gene": "BRCA2",
        "rs_id": "rs80359550",
        "risk_allele": "delT",
        "protective_allele": "wt",
        "effect_size": "6174delT founder mutation; 45-70% lifetime breast cancer risk; 10-20% ovarian cancer risk; elevated prostate cancer risk; carrier frequency 1/80 in Ashkenazi Jewish population",
        "mechanism": "BRCA2 encodes a mediator of RAD51-dependent homologous recombination DNA repair. The 6174delT (rs80359550) frameshift mutation produces a truncated protein unable to recruit RAD51 to double-strand breaks. This impairs high-fidelity DNA repair, promoting error-prone repair pathways and genomic instability. Male carriers have significantly elevated prostate cancer risk (up to 5-8x by age 65) and breast cancer risk. Like BRCA1, BRCA2-deficient tumors exhibit BRCAness and respond to PARP inhibitors and platinum agents.",
        "disease_associations": "Hereditary Breast Cancer,Ovarian Cancer,Prostate Cancer,Pancreatic Cancer",
        "text_chunk": "BRCA2 6174delT (rs80359550) is an Ashkenazi Jewish founder mutation with carrier frequency of 1/80. It disrupts RAD51-mediated DNA repair, conferring 45-70% lifetime breast cancer risk, 10-20% ovarian cancer risk, and significantly elevated prostate and pancreatic cancer risk in male carriers. PARP inhibitors are effective in BRCA2-deficient tumors. Male carriers should undergo enhanced prostate cancer screening starting at age 40. Cascade genetic testing of family members is strongly recommended."
    },
    {
        "id": "gv_hexa_1278instatc",
        "gene": "HEXA",
        "rs_id": "1278insTATC",
        "risk_allele": "insTATC",
        "protective_allele": "wt",
        "effect_size": "Tay-Sachs carrier frequency 1/30 in Ashkenazi Jewish population; complete hexosaminidase A deficiency in homozygotes",
        "mechanism": "HEXA encodes the alpha subunit of hexosaminidase A, a lysosomal enzyme that degrades GM2 ganglioside in neurons. The 1278insTATC frameshift insertion produces a non-functional alpha subunit. Homozygous individuals accumulate GM2 ganglioside in neurons, causing progressive neurodegeneration. Infantile Tay-Sachs presents with developmental regression at 6 months, cherry-red macula, and death by age 4-5. Carrier detection is performed via enzyme activity assay (reduced but functional in carriers) combined with molecular testing.",
        "disease_associations": "Tay-Sachs Disease (homozygous),GM2 Gangliosidosis,Neurodegenerative Disease",
        "text_chunk": "HEXA 1278insTATC is the most common Tay-Sachs mutation in the Ashkenazi Jewish population (carrier frequency 1/30). It causes complete loss of hexosaminidase A activity in homozygotes, leading to fatal GM2 ganglioside accumulation in neurons. Carrier screening via enzyme activity assay and molecular testing is standard of care for AJ individuals. Population-based carrier screening programs have reduced Tay-Sachs incidence by over 90% in the AJ population. Late-onset forms exist with residual enzyme activity."
    },
    {
        "id": "gv_fancc_ivs4",
        "gene": "FANCC",
        "rs_id": "IVS4+4A>T",
        "risk_allele": "T",
        "protective_allele": "A",
        "effect_size": "Fanconi Anemia Type C carrier frequency 1/89 in Ashkenazi Jewish population; bone marrow failure and cancer predisposition in homozygotes",
        "mechanism": "FANCC encodes Fanconi anemia complementation group C protein, a key component of the Fanconi anemia DNA interstrand crosslink repair pathway. The IVS4+4A>T splice site mutation causes exon 4 skipping, producing a non-functional protein. Homozygous individuals develop Fanconi anemia Type C, characterized by progressive bone marrow failure (pancytopenia), congenital anomalies (skeletal, renal, cardiac), and markedly elevated cancer risk including acute myeloid leukemia and head/neck squamous cell carcinoma.",
        "disease_associations": "Fanconi Anemia Type C (homozygous),Bone Marrow Failure,Acute Myeloid Leukemia,Squamous Cell Carcinoma",
        "text_chunk": "FANCC IVS4+4A>T is the predominant Fanconi Anemia Type C mutation in the Ashkenazi Jewish population (carrier frequency 1/89). Homozygotes develop progressive bone marrow failure, congenital anomalies, and extreme cancer susceptibility. The FANCC protein is essential for DNA interstrand crosslink repair. Carrier screening is included in expanded AJ carrier panels. Affected individuals require hematologic monitoring, often bone marrow transplantation, and cancer surveillance. Chromosomal breakage testing confirms diagnosis."
    },
    {
        "id": "gv_aspa_e285a",
        "gene": "ASPA",
        "rs_id": "E285A",
        "risk_allele": "C",
        "protective_allele": "A",
        "effect_size": "Canavan Disease carrier frequency 1/40 in Ashkenazi Jewish population; N-acetylaspartoacylase deficiency in homozygotes",
        "mechanism": "ASPA encodes aspartoacylase (N-acetylaspartoacylase), which hydrolyzes N-acetylaspartic acid (NAA) in the brain. The E285A missense mutation abolishes enzyme activity, causing NAA accumulation in the brain. NAA is an osmolyte that damages myelin through osmotic swelling and impaired oligodendrocyte lipid synthesis. Homozygous individuals develop Canavan disease, a fatal leukodystrophy with spongy white matter degeneration, macrocephaly, progressive motor and cognitive deterioration, and death typically in the first decade.",
        "disease_associations": "Canavan Disease (homozygous),Leukodystrophy,N-acetylaspartic Aciduria",
        "text_chunk": "ASPA E285A is the most common Canavan Disease mutation in the Ashkenazi Jewish population (carrier frequency 1/40). It abolishes N-acetylaspartoacylase activity, causing toxic N-acetylaspartic acid accumulation and spongy white matter degeneration. Canavan disease is a fatal leukodystrophy presenting in infancy with macrocephaly and progressive neurological decline. Carrier screening is included in standard AJ panels. Elevated urine NAA is diagnostic. Gene therapy trials are ongoing as potential treatment."
    },

    # --- Inflammatory Variants (4) ---
    {
        "id": "gv_il6_rs1800795",
        "gene": "IL6",
        "rs_id": "rs1800795",
        "risk_allele": "G",
        "protective_allele": "C",
        "effect_size": "G allele associated with higher IL-6 production, elevated baseline CRP, and enhanced inflammatory response",
        "mechanism": "IL6 encodes interleukin-6, a pleiotropic cytokine central to acute-phase inflammatory response, immune regulation, and hematopoiesis. The rs1800795 G>C promoter variant (-174G/C) modulates IL-6 transcription. The G allele is associated with higher constitutive and inducible IL-6 expression, leading to elevated baseline CRP levels, enhanced acute-phase response, and increased susceptibility to systemic inflammation. This variant modifies cardiovascular risk through inflammation-dependent atherogenesis.",
        "disease_associations": "Systemic Inflammation,Cardiovascular Disease,Type 2 Diabetes,Elevated CRP",
        "text_chunk": "IL6 rs1800795 (-174G/C) is a promoter variant that modulates IL-6 cytokine production. The G allele drives higher IL-6 expression, resulting in elevated baseline CRP and enhanced inflammatory response. GG homozygotes have significantly higher CRP levels independent of BMI or infection. This variant is critical for interpreting CRP and inflammatory biomarkers, as genetically elevated CRP may not indicate active disease. It also modifies cardiovascular risk assessment using hs-CRP as a risk stratifier."
    },
    {
        "id": "gv_tnf_rs1800629",
        "gene": "TNF",
        "rs_id": "rs1800629",
        "risk_allele": "A",
        "protective_allele": "G",
        "effect_size": "A allele increases TNF-alpha production by 2-fold; associated with increased sepsis susceptibility and severity",
        "mechanism": "TNF encodes tumor necrosis factor alpha, a master pro-inflammatory cytokine that initiates the inflammatory cascade, activates NF-kB signaling, and mediates septic shock pathophysiology. The rs1800629 G>A promoter variant (-308G/A) increases TNF-alpha transcription. The A allele is associated with 2-fold higher TNF-alpha production, increased susceptibility to severe sepsis, higher sepsis mortality, and more aggressive inflammatory diseases including rheumatoid arthritis and inflammatory bowel disease.",
        "disease_associations": "Sepsis Susceptibility,Rheumatoid Arthritis,Inflammatory Bowel Disease,Cerebral Malaria",
        "text_chunk": "TNF rs1800629 (-308G/A) is a promoter variant affecting TNF-alpha production. The A allele doubles TNF-alpha expression, increasing susceptibility to severe sepsis and inflammatory autoimmune diseases. AA carriers have higher sepsis mortality and more aggressive rheumatoid arthritis. This variant affects interpretation of inflammatory biomarkers (ESR, CRP, ferritin) and may predict response to TNF-alpha inhibitor therapy (infliximab, adalimumab). It is relevant for critical care risk stratification."
    },
    {
        "id": "gv_crp_rs1130864",
        "gene": "CRP",
        "rs_id": "rs1130864",
        "risk_allele": "T",
        "protective_allele": "C",
        "effect_size": "T allele raises baseline CRP by 0.3-0.6 mg/L independent of inflammatory status",
        "mechanism": "CRP encodes C-reactive protein, the primary acute-phase reactant synthesized by hepatocytes in response to IL-6 signaling. The rs1130864 variant in the CRP gene 3'-UTR affects mRNA stability and translation efficiency. The T allele is associated with constitutively higher CRP levels independent of inflammation, infection, or metabolic status. This genetic elevation of CRP confounds the use of hs-CRP for cardiovascular risk stratification and infection monitoring.",
        "disease_associations": "Elevated Baseline CRP,Modified Cardiovascular Risk Assessment,Confounded Inflammatory Markers",
        "text_chunk": "CRP rs1130864 is a variant in the CRP gene itself that raises baseline CRP levels by 0.3-0.6 mg/L independent of inflammation. The T allele increases CRP mRNA stability, producing constitutively higher protein levels. This is critical for interpreting hs-CRP in cardiovascular risk assessment, as genetically elevated CRP may misclassify patients into higher risk categories. Mendelian randomization studies using this variant suggest CRP is a marker, not a causal factor, in cardiovascular disease."
    },
    {
        "id": "gv_il10_rs1800896",
        "gene": "IL10",
        "rs_id": "rs1800896",
        "risk_allele": "A",
        "protective_allele": "G",
        "effect_size": "G allele increases IL-10 production; A allele associated with reduced anti-inflammatory capacity and increased autoimmune risk",
        "mechanism": "IL10 encodes interleukin-10, the primary anti-inflammatory cytokine that suppresses macrophage activation, reduces pro-inflammatory cytokine production, and promotes regulatory T-cell function. The rs1800896 promoter variant (-1082G/A) modulates IL-10 transcription. The A allele reduces IL-10 production, impairing anti-inflammatory feedback and increasing susceptibility to chronic inflammatory conditions, autoimmune diseases, and transplant rejection. The G allele supports robust anti-inflammatory responses.",
        "disease_associations": "Autoimmune Disease,Transplant Rejection,Chronic Inflammation,Inflammatory Bowel Disease",
        "text_chunk": "IL10 rs1800896 (-1082G/A) is a promoter variant affecting anti-inflammatory IL-10 cytokine production. The A allele reduces IL-10 expression, impairing anti-inflammatory feedback loops and increasing susceptibility to autoimmune diseases and chronic inflammation. AA homozygotes have lower IL-10 levels with reduced capacity to resolve inflammation. This variant helps interpret persistent inflammatory biomarker elevation and may predict response to immunosuppressive therapy and transplant outcomes."
    },

    # --- Neurological Variants (3) ---
    {
        "id": "gv_trem2_rs75932628",
        "gene": "TREM2",
        "rs_id": "rs75932628",
        "risk_allele": "T",
        "protective_allele": "C",
        "effect_size": "R47H variant; 2-4x increased Alzheimer's disease risk; comparable to single APOE E4 allele",
        "mechanism": "TREM2 encodes triggering receptor expressed on myeloid cells 2, a microglial surface receptor that senses lipid debris, apoptotic neurons, and amyloid-beta. The R47H variant (rs75932628 T allele) impairs TREM2 ligand binding, reducing microglial phagocytic clearance of amyloid plaques and damaged neurons. This leads to impaired microglial clustering around plaques, defective debris clearance, and sustained neuroinflammation. TREM2 R47H confers Alzheimer's risk comparable to a single APOE E4 allele but through an independent microglial mechanism.",
        "disease_associations": "Alzheimer's Disease,Frontotemporal Dementia,Nasu-Hakola Disease (homozygous),Neuroinflammation",
        "text_chunk": "TREM2 rs75932628 (R47H) is a rare microglial function variant conferring 2-4x increased Alzheimer's disease risk, comparable to a single APOE E4 allele. The variant impairs microglial phagocytosis of amyloid-beta plaques and neuronal debris, promoting neuroinflammation. TREM2 operates through a pathway independent of APOE, meaning carriers of both have compounded risk. This variant is relevant for neurological biomarker interpretation and emerging anti-amyloid and microglial-targeted therapies."
    },
    {
        "id": "gv_tomm40_rs10524523",
        "gene": "TOMM40",
        "rs_id": "rs10524523",
        "risk_allele": "VL",
        "protective_allele": "S",
        "effect_size": "Very long (VL) poly-T repeat associated with earlier Alzheimer's onset; modifies APOE E3 risk",
        "mechanism": "TOMM40 encodes translocase of outer mitochondrial membrane 40, essential for mitochondrial protein import. The rs10524523 poly-T length polymorphism in intron 6 is in strong linkage disequilibrium with the APOE locus. Very long (VL) poly-T repeats are associated with earlier age of Alzheimer's onset, particularly in APOE E3/E3 individuals who would otherwise be at average risk. The mechanism may involve impaired mitochondrial function and reduced energy metabolism in neurons, contributing to neurodegeneration.",
        "disease_associations": "Alzheimer's Disease (modified onset age),Mitochondrial Dysfunction,Cognitive Decline",
        "text_chunk": "TOMM40 rs10524523 is a poly-T length polymorphism near the APOE locus that modifies Alzheimer's disease age of onset. Very long (VL) poly-T repeats are associated with earlier onset, particularly in APOE E3/E3 carriers who would otherwise be at average risk. The variant affects mitochondrial protein import and neuronal energy metabolism. TOMM40 genotyping may refine Alzheimer's risk prediction beyond APOE genotype alone, especially for E3 homozygotes seeking personalized risk assessment."
    },
    {
        "id": "gv_bdnf_rs6265",
        "gene": "BDNF",
        "rs_id": "rs6265",
        "risk_allele": "A",
        "protective_allele": "G",
        "effect_size": "Val66Met substitution; A (Met) allele reduces activity-dependent BDNF secretion by 25-30%; affects hippocampal function and memory",
        "mechanism": "BDNF encodes brain-derived neurotrophic factor, a key neurotrophin supporting neuronal survival, synaptic plasticity, and hippocampal long-term potentiation (memory formation). The Val66Met variant (rs6265 A allele) impairs intracellular trafficking and activity-dependent secretion of mature BDNF from neurons. Met carriers have reduced hippocampal volume, impaired episodic memory, and altered stress response. The variant also affects exercise-induced BDNF release, neuroplasticity after stroke, and antidepressant response.",
        "disease_associations": "Memory Impairment,Depression,Anxiety Disorders,Reduced Neuroplasticity,Altered Antidepressant Response",
        "text_chunk": "BDNF rs6265 (Val66Met) affects brain-derived neurotrophic factor secretion. The Met (A) allele reduces activity-dependent BDNF release by 25-30%, leading to smaller hippocampal volume and impaired episodic memory. Met carriers show reduced neuroplasticity and may have altered response to antidepressants and cognitive rehabilitation. Exercise-induced BDNF elevation is also blunted in Met carriers. This variant is relevant for interpreting cognitive biomarkers and personalizing neurological and psychiatric treatment approaches."
    },

    # --- Bone Health Variants (2) ---
    {
        "id": "gv_col1a1_rs1800012",
        "gene": "COL1A1",
        "rs_id": "rs1800012",
        "risk_allele": "T",
        "protective_allele": "G",
        "effect_size": "Sp1 binding site variant; T allele associated with reduced bone mineral density, increased fracture risk (OR 1.3-1.5)",
        "mechanism": "COL1A1 encodes the alpha-1 chain of type I collagen, the primary structural protein of bone, skin, and tendons. The rs1800012 variant affects an Sp1 transcription factor binding site in intron 1, altering the ratio of alpha-1 to alpha-2 collagen chains. The T allele increases alpha-1 chain production, producing collagen with abnormal alpha-1/alpha-2 stoichiometry that has reduced mechanical strength. This leads to lower bone mineral density, increased osteoporotic fracture risk, and potentially altered bone quality independent of BMD.",
        "disease_associations": "Osteoporosis,Reduced Bone Mineral Density,Vertebral Fractures,Osteopenia",
        "text_chunk": "COL1A1 rs1800012 (Sp1) is a type I collagen variant affecting bone structure. The T allele alters collagen alpha-1/alpha-2 chain ratios, reducing bone mechanical strength and lowering bone mineral density (OR 1.3-1.5 for fracture). TT homozygotes have significantly increased vertebral fracture risk. This variant is important for interpreting DEXA scan results and bone turnover markers (CTX, P1NP, osteocalcin). It may guide earlier osteoporosis intervention in carriers, particularly postmenopausal women."
    },
    {
        "id": "gv_esr1_rs2234693",
        "gene": "ESR1",
        "rs_id": "rs2234693",
        "risk_allele": "C",
        "protective_allele": "T",
        "effect_size": "PvuII variant; C allele associated with lower bone mineral density in postmenopausal women; modifies estrogen response",
        "mechanism": "ESR1 encodes estrogen receptor alpha, the primary mediator of estrogen effects on bone metabolism including osteoblast survival, osteoclast apoptosis, and calcium homeostasis. The rs2234693 (PvuII) variant in intron 1 affects ESR1 transcription through altered binding of the B-myb transcription factor. The C allele is associated with reduced estrogen receptor expression in bone, diminished bone-protective estrogen signaling, and accelerated postmenopausal bone loss. The effect is most pronounced after menopause when circulating estrogen declines.",
        "disease_associations": "Postmenopausal Osteoporosis,Reduced Bone Mineral Density,Fracture Risk,Modified HRT Response",
        "text_chunk": "ESR1 rs2234693 (PvuII) affects estrogen receptor alpha expression in bone. The C allele reduces estrogen-mediated bone protection, leading to lower bone mineral density particularly in postmenopausal women. CC homozygotes have accelerated postmenopausal bone loss and increased fracture risk. This variant is relevant for interpreting bone density (DEXA) and bone turnover markers, and may influence decisions about hormone replacement therapy and osteoporosis prevention strategies in postmenopausal women."
    },

    # --- Cardiovascular Extended (3) ---
    # NOTE: gv_pcsk9_rs11591147 already exists; skip it
    {
        "id": "gv_apoa5_rs662799",
        "gene": "APOA5",
        "rs_id": "rs662799",
        "risk_allele": "C",
        "protective_allele": "T",
        "effect_size": "C allele raises triglycerides by 15-30%; associated with hypertriglyceridemia and metabolic syndrome",
        "mechanism": "APOA5 encodes apolipoprotein A-V, a potent modulator of plasma triglyceride levels that activates lipoprotein lipase and inhibits hepatic VLDL production. The rs662799 promoter variant (-1131T>C) reduces APOA5 transcription, lowering circulating ApoA-V levels. Reduced ApoA-V impairs triglyceride-rich lipoprotein catabolism, raising plasma triglycerides by 15-30%. This variant is one of the strongest common genetic determinants of triglyceride levels and interacts with dietary fat intake.",
        "disease_associations": "Hypertriglyceridemia,Metabolic Syndrome,Cardiovascular Disease,Pancreatitis Risk",
        "text_chunk": "APOA5 rs662799 (-1131T>C) is a promoter variant affecting apolipoprotein A-V expression. The C allele reduces ApoA-V production, raising triglycerides by 15-30% through impaired lipoprotein lipase activation. CC homozygotes have significantly elevated triglyceride levels and increased pancreatitis risk. This variant is important for interpreting triglyceride biomarker levels and personalizing dietary fat intake recommendations. It interacts with omega-3 supplementation response and statin-fibrate combination therapy decisions."
    },
    {
        "id": "gv_pon1_rs662",
        "gene": "PON1",
        "rs_id": "rs662",
        "risk_allele": "G",
        "protective_allele": "A",
        "effect_size": "Q192R substitution; R (G) allele alters paraoxonase activity and HDL antioxidant function; modifies cardiovascular risk",
        "mechanism": "PON1 encodes paraoxonase 1, an HDL-associated enzyme that prevents LDL oxidation by hydrolyzing lipid peroxides and homocysteine thiolactone. The Q192R variant (rs662 G allele, Arg192) alters substrate-specific enzyme activity: R192 has higher paraoxonase activity but lower lactonase and arylesterase activity. The net effect is reduced HDL antioxidant capacity against lipid peroxidation, potentially diminishing the cardioprotective properties of HDL. This variant affects HDL quality beyond simple HDL-C measurement.",
        "disease_associations": "Modified HDL Function,Cardiovascular Disease,Oxidative Stress,Organophosphate Sensitivity",
        "text_chunk": "PON1 rs662 (Q192R) affects paraoxonase 1, an HDL-associated antioxidant enzyme. The R (G) allele alters enzyme substrate specificity, reducing HDL antioxidant capacity against lipid peroxidation. This variant affects HDL functional quality beyond HDL-C levels, explaining why some individuals with normal HDL-C still have elevated cardiovascular risk. It is relevant for interpreting oxidized LDL, HDL function assays, and may influence the clinical significance of measured HDL-C levels."
    },
]

# ---------------------------------------------------------------------------
# Part 2: New Drug-Biomarker Interactions
# ---------------------------------------------------------------------------

NEW_INTERACTIONS = [
    # --- Medication-Biomarker Effects (male patient context) ---
    {
        "id": "dxi_atorvastatin_biomarker_effects",
        "drug": "Atorvastatin",
        "gene": "SLCO1B1/HMGCR",
        "interaction_type": "biomarker_effect",
        "severity": "moderate",
        "alternative": "Monitor LDL, CK, CoQ10, GGT; supplement CoQ10 if depleted; consider SLCO1B1 genotyping for myopathy risk",
        "text_chunk": "Atorvastatin affects multiple biomarkers: it lowers LDL and Total Cholesterol (primary therapeutic effect) but depletes CoQ10 by inhibiting the mevalonate pathway (shared with cholesterol synthesis). CoQ10 depletion may contribute to statin myopathy, detectable by CK elevation. GGT may be mildly elevated due to hepatic enzyme induction. SLCO1B1 genotype (rs4149056) significantly affects statin myopathy risk, with *5 carriers having higher systemic exposure and 17x increased myopathy risk for simvastatin. CoQ10 supplementation (100-200 mg/day) may mitigate muscle symptoms."
    },
    {
        "id": "dxi_lisinopril_biomarker_effects",
        "drug": "Lisinopril",
        "gene": "ACE/REN",
        "interaction_type": "biomarker_effect",
        "severity": "moderate",
        "alternative": "Monitor potassium, creatinine, BUN regularly; ARB (losartan) if ACE inhibitor cough develops",
        "text_chunk": "Lisinopril (ACE inhibitor) affects multiple renal and electrolyte biomarkers. It raises serum potassium by reducing aldosterone-mediated potassium excretion, creating hyperkalemia risk especially with renal impairment, potassium supplements, or potassium-sparing diuretics. Creatinine may increase initially (up to 30%) due to reduced efferent arteriolar tone; this is expected and usually stabilizes. BUN may rise proportionally. ACE genotype (I/D polymorphism) may influence blood pressure response magnitude. Regular monitoring of potassium, creatinine, and BUN is essential, particularly in patients with CKD."
    },
    {
        "id": "dxi_fish_oil_biomarker_effects",
        "drug": "Fish Oil (Omega-3)",
        "gene": "FADS1/FADS2",
        "interaction_type": "biomarker_effect",
        "severity": "minor",
        "alternative": "Monitor INR if on anticoagulants; adjust dose based on omega-3 index and triglyceride response",
        "text_chunk": "Fish oil (omega-3 fatty acids, EPA/DHA) affects several biomarkers. It lowers triglycerides by 15-30% at therapeutic doses (2-4 g/day) by reducing hepatic VLDL production. LDL may slightly increase (especially with DHA-predominant formulations). Fish oil inhibits platelet aggregation and may prolong bleeding time, increasing INR when combined with anticoagulants (warfarin, aspirin). It raises the omega-3 index toward cardioprotective levels (>8%). FADS1 genotype affects endogenous omega-3 synthesis and supplementation requirements. Platelet function should be monitored in patients on dual antiplatelet or anticoagulant therapy."
    },
    {
        "id": "dxi_l_methylfolate_biomarker_effects",
        "drug": "L-Methylfolate",
        "gene": "MTHFR",
        "interaction_type": "biomarker_effect",
        "severity": "minor",
        "alternative": "Standard folic acid for MTHFR wild-type; monitor homocysteine response",
        "text_chunk": "L-Methylfolate (5-MTHF) is the bioactive form of folate that bypasses the MTHFR enzyme. It directly lowers homocysteine by serving as the methyl donor for methionine synthase-mediated homocysteine remethylation. Serum folate levels increase with supplementation. L-Methylfolate is particularly important for MTHFR C677T TT homozygotes (70% reduced enzyme activity) who cannot efficiently convert folic acid to active 5-MTHF. Homocysteine reduction of 20-30% is expected with adequate methylfolate supplementation. It also supports SAMe-dependent methylation reactions relevant to DNA repair, neurotransmitter synthesis, and epigenetic regulation."
    },
    {
        "id": "dxi_vitamin_d3_biomarker_effects",
        "drug": "Vitamin D3 Supplementation",
        "gene": "VDR/GC",
        "interaction_type": "biomarker_effect",
        "severity": "minor",
        "alternative": "Adjust dose based on 25-OH Vitamin D levels and VDR/GC genotype; monitor calcium",
        "text_chunk": "Vitamin D3 (cholecalciferol) supplementation raises 25-OH Vitamin D levels, the primary biomarker for vitamin D status. It can elevate serum calcium through enhanced intestinal calcium absorption, particularly at high doses (>4000 IU/day) or in patients with granulomatous diseases. Vitamin D3 suppresses PTH by normalizing calcium-vitamin D axis feedback. VDR genotype (rs2228570) affects biological response to circulating vitamin D. GC genotype (rs2282679) affects vitamin D binding protein levels and measured 25-OH Vitamin D. Monitoring calcium, 25-OH Vitamin D, and PTH is recommended during supplementation."
    },
    {
        "id": "dxi_coq10_biomarker_effects",
        "drug": "CoQ10 Supplementation",
        "gene": "COQ2/HMGCR",
        "interaction_type": "biomarker_effect",
        "severity": "minor",
        "alternative": "Ubiquinol form for better absorption; standard dose 100-200 mg/day with statin therapy",
        "text_chunk": "CoQ10 (ubiquinone/ubiquinol) supplementation restores CoQ10 levels depleted by statin therapy. Statins inhibit HMG-CoA reductase, reducing mevalonate pathway flux needed for both cholesterol and CoQ10 synthesis. CoQ10 supplementation at 100-200 mg/day may reduce statin-associated myalgia and CK elevation. It may provide modest blood pressure reduction (3-5 mmHg systolic) through improved endothelial function. CoQ10 supports mitochondrial electron transport chain function. Ubiquinol form has superior bioavailability compared to ubiquinone, particularly in older adults."
    },

    # --- Medication-Biomarker Effects (female patient context) ---
    {
        "id": "dxi_oral_contraceptives_biomarker_effects",
        "drug": "Oral Contraceptives (Ethinyl Estradiol/Norgestimate)",
        "gene": "F5/F2/SERPINE1",
        "interaction_type": "biomarker_effect",
        "severity": "major",
        "alternative": "Progestin-only methods or non-hormonal contraception for high VTE risk; IUD for Factor V Leiden carriers",
        "text_chunk": "Oral contraceptives containing ethinyl estradiol have widespread biomarker effects. They increase SHBG (2-3x), raising total hormone levels while decreasing free testosterone. TBG rises, increasing total T4 (but not free T4, which remains normal). Triglycerides increase by 30-50% through hepatic VLDL stimulation. Fibrinogen and clotting factors (VII, VIII, X) increase, elevating VTE risk 3-5x (higher with Factor V Leiden or Prothrombin G20210A). CRP may rise due to hepatic acute-phase protein induction. Free testosterone decreases, which is therapeutic for PCOS/acne but may cause low libido."
    },
    {
        "id": "dxi_iron_bisglycinate_biomarker_effects",
        "drug": "Iron Bisglycinate",
        "gene": "HFE/TMPRSS6",
        "interaction_type": "biomarker_effect",
        "severity": "minor",
        "alternative": "IV iron (ferric carboxymaltose) for malabsorption or intolerance; monitor ferritin and TSAT",
        "text_chunk": "Iron bisglycinate supplementation raises ferritin (the primary iron storage biomarker) and transferrin saturation (TSAT). It supports hemoglobin synthesis in iron-deficiency anemia. Bisglycinate chelate form has superior absorption and GI tolerability compared to ferrous sulfate. HFE genotype should be checked before supplementation in patients with unexplained elevated ferritin, as C282Y carriers risk iron overload. Target ferritin is typically 50-100 ng/mL for iron-replete status. Over-supplementation can cause ferritin >200 ng/mL with potential oxidative stress. TMPRSS6 variants affect hepcidin regulation and iron absorption efficiency."
    },
    {
        "id": "dxi_prenatal_dha_biomarker_effects",
        "drug": "Prenatal DHA",
        "gene": "FADS1/FADS2",
        "interaction_type": "biomarker_effect",
        "severity": "minor",
        "alternative": "Algal DHA for vegetarian/vegan patients; adjust dose based on omega-3 index",
        "text_chunk": "Prenatal DHA supplementation (200-600 mg/day) increases the omega-3 index and DHA levels specifically. DHA is the predominant omega-3 fatty acid in fetal brain and retinal tissue, essential for neural development. Maternal DHA status directly affects fetal neural development through placental transfer. FADS1/FADS2 genotype affects endogenous DHA synthesis from shorter-chain precursors; carriers of low-activity alleles require more preformed DHA. Prenatal DHA may slightly increase LDL-C (a known DHA-specific effect). Adequate maternal DHA (omega-3 index >8%) is associated with reduced preterm birth risk."
    },

    # --- Additional Drug-Gene Interactions (6) ---
    {
        "id": "dxi_tacrolimus_cyp3a5",
        "drug": "Tacrolimus",
        "gene": "CYP3A5",
        "interaction_type": "substrate",
        "severity": "major",
        "alternative": "Genotype-guided dosing: CYP3A5 expressors (*1/*1, *1/*3) need 1.5-2x standard dose; non-expressors (*3/*3) use standard dose",
        "text_chunk": "Tacrolimus is a calcineurin inhibitor immunosuppressant critical for transplant rejection prevention, primarily metabolized by CYP3A5 (and CYP3A4). CYP3A5 expressors (*1/*1, *1/*3) have significantly higher tacrolimus clearance, requiring 1.5-2x higher doses to achieve therapeutic trough levels (5-15 ng/mL). Non-expressors (*3/*3) need standard or lower doses. CYP3A5*1 is more prevalent in African Americans (~60%) versus Europeans (~15%), contributing to racial disparities in transplant dosing. CPIC provides genotype-guided dosing recommendations. Therapeutic drug monitoring (TDM) is essential regardless of genotype."
    },
    {
        "id": "dxi_phenytoin_cyp2c9",
        "drug": "Phenytoin",
        "gene": "CYP2C9",
        "interaction_type": "substrate",
        "severity": "major",
        "alternative": "Levetiracetam or valproic acid (non-CYP2C9 anticonvulsants); CYP2C9 PM requires 25-50% dose reduction",
        "text_chunk": "Phenytoin is a narrow therapeutic index anticonvulsant (10-20 mcg/mL) primarily metabolized by CYP2C9 (with minor CYP2C19 contribution). CYP2C9 poor metabolizers (*2/*2, *2/*3, *3/*3) have dramatically reduced phenytoin clearance with nonlinear accumulation kinetics, leading to toxicity: ataxia, nystagmus, diplopia, cognitive impairment, and potentially fatal cardiovascular collapse. CPIC recommends 25-50% dose reduction for CYP2C9 PMs with intensive TDM. HLA-B*15:02 testing is also required before initiation in Southeast Asian populations to prevent Stevens-Johnson syndrome."
    },
    {
        "id": "dxi_aripiprazole_cyp2d6",
        "drug": "Aripiprazole",
        "gene": "CYP2D6",
        "interaction_type": "substrate",
        "severity": "moderate",
        "alternative": "Reduce dose by 50-67% for CYP2D6 poor metabolizers; consider quetiapine or ziprasidone as alternatives",
        "text_chunk": "Aripiprazole is an atypical antipsychotic metabolized by CYP2D6 (and CYP3A4) to its active metabolite dehydro-aripiprazole. CYP2D6 poor metabolizers have approximately 80% higher aripiprazole exposure, increasing risk of akathisia, extrapyramidal symptoms, somnolence, and QTc prolongation. FDA labeling recommends 50-67% dose reduction for known CYP2D6 PMs (e.g., 5-10 mg instead of 15-30 mg). Concurrent use of strong CYP2D6 inhibitors (paroxetine, fluoxetine, quinidine) phenocopies poor metabolizer status and requires similar dose reduction."
    },
    {
        "id": "dxi_carvedilol_cyp2d6",
        "drug": "Carvedilol",
        "gene": "CYP2D6",
        "interaction_type": "substrate",
        "severity": "moderate",
        "alternative": "Start at lower dose for CYP2D6 PMs; bisoprolol or nebivolol as CYP2D6-independent alternatives",
        "text_chunk": "Carvedilol is a non-selective beta-blocker with alpha-1 blocking activity, metabolized primarily by CYP2D6 (and CYP2C9). CYP2D6 poor metabolizers have 2-3x higher plasma concentrations, increasing risk of hypotension, bradycardia, dizziness, and fatigue. DPWG recommends starting at the lowest available dose in CYP2D6 PMs with careful uptitration and monitoring of heart rate and blood pressure. Bisoprolol (renally cleared) is a CYP2D6-independent beta-blocker alternative. Carvedilol is widely used in heart failure, making genotype-guided dosing particularly important for this vulnerable population."
    },
    {
        "id": "dxi_losartan_cyp2c9",
        "drug": "Losartan",
        "gene": "CYP2C9",
        "interaction_type": "substrate",
        "severity": "moderate",
        "alternative": "Valsartan, irbesartan, or candesartan (ARBs less dependent on CYP2C9 activation)",
        "text_chunk": "Losartan is a prodrug requiring CYP2C9-mediated oxidation to its active metabolite E-3174 (EXP 3174), which has 10-40x greater AT1 receptor affinity. CYP2C9 poor metabolizers (*2/*3, *3/*3) have reduced conversion to the active metabolite, potentially resulting in suboptimal blood pressure control and reduced renal protective effects. CYP2C9 *3 carriers show approximately 50% lower active metabolite levels. Alternative ARBs (valsartan, candesartan, irbesartan) are active drugs that do not require CYP2C9 bioactivation and may be preferred for CYP2C9 PMs."
    },
    {
        "id": "dxi_lamotrigine_ugt1a4",
        "drug": "Lamotrigine",
        "gene": "UGT1A4",
        "interaction_type": "substrate",
        "severity": "moderate",
        "alternative": "Dose adjustment with OCP co-administration; TDM recommended; levetiracetam as alternative anticonvulsant",
        "text_chunk": "Lamotrigine is an anticonvulsant and mood stabilizer primarily cleared by UGT1A4 glucuronidation. UGT1A4 variants (*2, *3) affect glucuronidation efficiency and lamotrigine clearance. Critically, oral contraceptives (OCPs) containing ethinyl estradiol induce UGT1A4, reducing lamotrigine levels by approximately 50%, potentially causing seizure breakthrough. Conversely, OCP discontinuation can double lamotrigine levels, risking toxicity (dizziness, ataxia, Stevens-Johnson syndrome). Women on lamotrigine require dose adjustment when starting or stopping OCPs, with therapeutic drug monitoring to maintain levels of 3-14 mcg/mL."
    },
]


def main():
    # --- Load existing data ---
    print("Loading existing data files...")
    variants = load_json(VARIANTS_PATH)
    interactions = load_json(INTERACTIONS_PATH)

    existing_variant_ids = {v["id"] for v in variants}
    existing_interaction_ids = {i["id"] for i in interactions}

    print(f"  Existing genetic variants: {len(variants)}")
    print(f"  Existing drug interactions: {len(interactions)}")

    # --- Add new variants (skip duplicates) ---
    print("\nAdding new genetic variants...")
    added_v = 0
    skipped_v = 0
    for v in NEW_VARIANTS:
        if v["id"] in existing_variant_ids:
            print(f"  SKIP (duplicate): {v['id']}")
            skipped_v += 1
        else:
            variants.append(v)
            existing_variant_ids.add(v["id"])
            added_v += 1
            print(f"  ADDED: {v['id']}")

    print(f"\n  Variants added: {added_v}, skipped: {skipped_v}")
    save_json(VARIANTS_PATH, variants)

    # --- Add new interactions (skip duplicates) ---
    print("\nAdding new drug interactions...")
    added_i = 0
    skipped_i = 0
    for i in NEW_INTERACTIONS:
        if i["id"] in existing_interaction_ids:
            print(f"  SKIP (duplicate): {i['id']}")
            skipped_i += 1
        else:
            interactions.append(i)
            existing_interaction_ids.add(i["id"])
            added_i += 1
            print(f"  ADDED: {i['id']}")

    print(f"\n  Interactions added: {added_i}, skipped: {skipped_i}")
    save_json(INTERACTIONS_PATH, interactions)

    # --- Summary ---
    print("\n" + "=" * 60)
    print("EXPANSION COMPLETE")
    print(f"  Genetic variants: {len(variants)} total ({added_v} new)")
    print(f"  Drug interactions: {len(interactions)} total ({added_i} new)")
    print("=" * 60)


if __name__ == "__main__":
    main()
