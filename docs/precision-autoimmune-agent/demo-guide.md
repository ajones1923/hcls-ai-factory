# Precision Autoimmune Intelligence Agent -- Demo Guide

**Author:** Adam Jones
**Date:** March 2026

> Step-by-step walkthrough of all 10 Streamlit UI tabs using three demo patients. Use this guide for live demos.

---

## Prerequisites

Before starting the demo, ensure:

1. **Milvus** is running on `localhost:19530` with 14 collections created
2. **Streamlit UI** is running on `http://localhost:8531`
3. **Anthropic API key** is configured (for RAG queries in Tab 1)
4. **Demo patient data** has been ingested (clinical PDFs for 9 patients)

Verify everything is ready:

```bash
cd ai_agent_adds/precision_autoimmune_agent

# Create collections
python3 scripts/setup_collections.py

# Start API (seeds demo data on first run)
uvicorn api.main:app --host 0.0.0.0 --port 8532 &

# Ingest demo patient PDFs
curl -X POST http://localhost:8532/ingest/demo-data

# Start UI
streamlit run app/autoimmune_ui.py --server.port 8531
```

---

## Demo Flow Overview

| Tab | Time | Key Talking Points |
|---|---|---|
| 1. Clinical Query | 3 min | RAG-powered Q&A with evidence from 14 collections, citation scoring |
| 2. Patient Analysis | 3 min | Full pipeline: antibodies + HLA + activity + flare + biologics |
| 3. Document Ingest | 2 min | Upload clinical PDFs, automatic chunking and embedding |
| 4. Diagnostic Odyssey | 3 min | Timeline visualization, diagnostic delay analysis, missed opportunities |
| 5. Autoantibody Panel | 3 min | Interactive antibody interpretation, disease associations with sensitivity/specificity |
| 6. HLA Analysis | 2 min | HLA-disease lookup, odds ratios, shared epitope explanation |
| 7. Disease Activity | 2 min | DAS28-CRP, SLEDAI-2K, BASDAI scoring with threshold visualization |
| 8. Flare Prediction | 2 min | Biomarker-based flare risk, contributing vs protective factors |
| 9. Therapy Advisor | 3 min | Biologic therapy selection with PGx, contraindications, monitoring |
| 10. Knowledge Base | 2 min | Collection stats, evidence explorer, knowledge version info |
| **Total** | **~25 min** | |

---

## Patient 1: Sarah Mitchell (SLE) -- HCLS-AUTO-2026-SM001

### Key Narrative

> A 28-year-old female presenting with a classic SLE serological profile. Her diagnostic odyssey spans multiple years -- from initial complaints of fatigue and joint pain dismissed as "stress" through progressive multi-organ involvement culminating in lupus nephritis. She has 27+ clinical documents spanning PCP visits, rheumatology consultations, nephrology referral, renal biopsy, and longitudinal lab monitoring.

### Tab 5: Autoantibody Panel

1. **Load Sarah Mitchell's autoantibody results:**
   - ANA: 1:640, homogeneous pattern -- **POSITIVE**
   - anti-dsDNA: elevated -- **POSITIVE**
   - anti-Smith: positive -- **POSITIVE**
   - anti-SSA/Ro: positive -- **POSITIVE**
2. **Point out disease associations:**
   - ANA positive in 95% of SLE (sensitivity 0.95, specificity 0.65)
   - anti-dsDNA: specificity 0.95 for SLE -- this is highly diagnostic
   - anti-Smith: specificity 0.99 for SLE -- pathognomonic
3. **Discuss the pattern:** Homogeneous ANA pattern correlates with anti-dsDNA positivity. The combination of ANA + anti-dsDNA + anti-Smith is essentially diagnostic for SLE.

**Key talking point:** *"ANA alone is sensitive but not specific -- many healthy individuals are ANA positive. But when you see anti-dsDNA at specificity 0.95 and anti-Smith at 0.99, you have a serological profile that is virtually pathognomonic for SLE."*

### Tab 7: Disease Activity

1. Navigate to Disease Activity tab
2. **Select SLEDAI-2K scoring for SLE**
3. **Enter biomarkers:** CRP elevated, anti-dsDNA rising, complement C3/C4 dropping, proteinuria present
4. **Show threshold visualization:** Remission (<0), Low (<4), Moderate (<8), High (>=12)
5. **Point out multi-organ involvement:** Renal (proteinuria), hematologic (lymphopenia), immunologic (low complement, rising anti-dsDNA)

**Key talking point:** *"SLEDAI-2K captures the multi-organ nature of SLE. This patient scores high not just because of joint involvement but because of renal, hematologic, and immunologic criteria -- each contributing independently to the overall activity score."*

### Tab 8: Flare Prediction

1. Navigate to Flare Prediction tab
2. **Show SLE flare pattern biomarkers:** anti-dsDNA titer, complement C3, complement C4, lymphocyte count, proteinuria
3. **Highlight contributing factors:**
   - Rising anti-dsDNA titer (>25% increase over 60 days)
   - Falling complement C3 and C4 (>20% decrease over 30 days)
4. **Protective factors:** Sun avoidance, medication adherence
5. **Recommended monitoring:** Repeat anti-dsDNA, complement, and urinalysis in 2-4 weeks

**Key talking point:** *"The combination of rising anti-dsDNA with falling complement is the classic SLE flare warning pattern. This serological shift often precedes clinical flare by 4-8 weeks -- giving us a window for preemptive intervention."*

### Tab 9: Therapy Advisor

1. Navigate to Therapy Advisor tab
2. **Show Belimumab recommendation:**
   - Drug class: BLyS inhibitor
   - Mechanism: Anti-BLyS (BAFF) monoclonal antibody
   - PGx: TNFSF13B genotype affects baseline BLyS levels and response
   - Monitoring: Depression screening, immunoglobulin levels, SLEDAI-2K tracking
3. **Show Rituximab as second-line option:**
   - Drug class: Anti-CD20 B-cell depleter
   - PGx: FCGR3A V158F affects ADCC efficacy
   - Monitoring: Immunoglobulin levels, CD19/CD20 counts, PML surveillance

**Key talking point:** *"Belimumab is the first biologic specifically approved for SLE. It targets BLyS/BAFF, the B-cell survival factor that drives autoantibody production. The PGx consideration here is that TNFSF13B genotype affects baseline BLyS levels, which can predict treatment response."*

### Tab 4: Diagnostic Odyssey

1. Navigate to Diagnostic Odyssey tab
2. **Load Sarah Mitchell's timeline** -- 27+ clinical documents
3. **Show the timeline visualization:**
   - First symptom: fatigue and joint pain
   - Multiple PCP visits over months before rheumatology referral
   - Key turning point: autoimmune panel revealing ANA + anti-dsDNA
   - Progressive multi-organ involvement leading to nephrology referral
   - Renal biopsy confirming lupus nephritis
4. **Point out diagnostic delay metrics:**
   - Number of specialists seen
   - Misdiagnosis events (if any)
   - Time from first symptom to diagnosis

**Key talking point:** *"The average diagnostic delay for SLE is 4.6 years. This agent reconstructs the entire diagnostic odyssey from fragmented clinical records, identifying missed opportunities and turning points -- the kind of pattern recognition that could accelerate diagnosis for future patients."*

---

## Patient 2: Rachel Thompson (Rheumatoid Arthritis) -- HCLS-AUTO-2026-RT002

### Key Narrative

> A 52-year-old female with seropositive rheumatoid arthritis. RF and anti-CCP positive, HLA-DRB1*04:01 carrier (shared epitope). Her profile demonstrates the convergence of serology, genetics, and disease activity scoring in guiding treatment decisions.

### Tab 5: Autoantibody Panel

1. **Load Rachel Thompson's antibody results:**
   - RF (Rheumatoid Factor): elevated -- **POSITIVE** (sensitivity 0.70, specificity 0.85 for RA)
   - anti-CCP: elevated -- **POSITIVE** (sensitivity 0.67, specificity 0.95 for RA)
2. **Discuss dual positivity:** RF + anti-CCP together has higher positive predictive value than either alone. anti-CCP specificity of 0.95 makes this highly diagnostic.
3. **Show serology scoring:** Per 2010 ACR/EULAR criteria, high-positive RF or anti-CCP = 3 points toward classification threshold of 6

**Key talking point:** *"Anti-CCP is the game-changer in RA diagnosis. While RF has been used for decades, it is not specific to RA -- it appears in Sjogren's, hepatitis C, and even healthy elderly. Anti-CCP at 95% specificity targets citrullinated peptides, which are the actual pathogenic target in RA."*

### Tab 6: HLA Analysis

1. Navigate to HLA Analysis tab
2. **Enter HLA-DRB1*04:01:**
   - Disease: Rheumatoid Arthritis
   - Odds Ratio: 4.2
   - PMID: 20301572
   - Note: Shared epitope hypothesis -- citrullinated peptide binding
3. **Explain the shared epitope:** HLA-DRB1*04:01 encodes a specific amino acid sequence in the peptide-binding groove that preferentially presents citrullinated peptides to T cells, driving the anti-CCP response.

**Key talking point:** *"Shared epitope HLA-DRB1*04:01 means higher risk of erosive disease. This is not just a risk allele -- it is the mechanistic link between genetics and pathology. The shared epitope binds citrullinated peptides, which is why anti-CCP positive patients with this HLA allele have the most aggressive disease course."*

### Tab 7: Disease Activity

1. Navigate to Disease Activity tab
2. **Show DAS28-CRP scoring for RA:**
   - Components: tender joints (28-count), swollen joints (28-count), CRP, patient global VAS
   - Thresholds: Remission (<2.6), Low (2.6-3.2), Moderate (3.2-5.1), High (>5.1)
3. **Also show CDAI and SDAI as alternatives** -- note these do not require inflammatory markers

**Key talking point:** *"DAS28-CRP is the workhorse of RA disease activity monitoring. The treat-to-target strategy aims for remission (DAS28 <2.6) or at minimum low disease activity (<3.2). Every 3-6 month assessment drives treatment escalation decisions."*

### Tab 9: Therapy Advisor

1. Navigate to Therapy Advisor tab
2. **Show multiple TNF inhibitor options:**
   - Adalimumab: PGx -- HLA-DRB1*03:01 associated with anti-drug antibody formation
   - Etanercept: PGx -- HLA-DRB1 shared epitope correlates with better TNF-i response
3. **Show Tocilizumab as alternative:**
   - IL-6 receptor inhibitor
   - PGx: IL6R Asp358Ala (rs2228145) affects receptor binding
   - Key caveat: CRP unreliable for infection monitoring on IL-6R blockade
4. **Show Tofacitinib (JAK inhibitor):**
   - PGx: CYP3A4 and CYP2C19 metabolism -- dose adjustment with strong inhibitors
   - VTE risk in patients >50 with CV risk factors

**Key talking point:** *"For seropositive RA with shared epitope, TNF inhibitors are first-line biologics. But the choice between adalimumab, tocilizumab, and tofacitinib depends on PGx factors. Rachel's HLA-DRB1*04:01 actually predicts better response to TNF inhibitors, while IL6R polymorphisms might predict tocilizumab response."*

---

## Patient 3: James Cooper (T1D + Celiac Overlap) -- HCLS-AUTO-2026-JC003

### Key Narrative

> A 19-year-old male presenting with autoimmune overlap syndrome: Type 1 Diabetes and Celiac Disease. He presented with DKA and was subsequently found to have celiac disease on screening. HLA-DQ2/DQ8 positive, linking both conditions through shared genetic susceptibility. His clinical records span his DKA presentation, endocrinology management, GI consultation, duodenal biopsy, and dietary management.

### Tab 5: Autoantibody Panel

1. **Load James Cooper's dual antibody panel:**
   - T1D antibodies: GAD65 antibody positive, IA-2 antibody positive
   - Celiac antibodies: anti-tTG IgA elevated (sensitivity 0.93, specificity 0.97 for celiac)
2. **Show dual-disease serological profile:** GAD65 confirms autoimmune beta-cell destruction while anti-tTG confirms celiac. The co-occurrence is not coincidence -- it is driven by shared HLA susceptibility.

**Key talking point:** *"James has two autoimmune diseases driven by the same genetic susceptibility. GAD65 antibodies confirm T1D, while anti-tTG at 93% sensitivity and 97% specificity confirms celiac. The overlap is not random -- 5-10% of T1D patients have celiac disease."*

### Tab 6: HLA Analysis

1. Navigate to HLA Analysis tab
2. **Enter HLA-DQB1*02:01 (DQ2):**
   - Celiac disease: OR=7.0 (PMID:17554300) -- forms the HLA-DQ2 heterodimer with DQA1*05:01
   - Type 1 Diabetes: OR=3.0 (PMID:17554300)
3. **Enter HLA-DQB1*03:02 (DQ8):**
   - Type 1 Diabetes: OR=6.5 (PMID:17554300) -- HLA-DQ8 heterodimer
4. **Explain the shared pathway:** Both celiac and T1D are driven by HLA class II molecules that present specific antigens to CD4+ T cells -- deamidated gliadin peptides in celiac, beta-cell antigens in T1D.

**Key talking point:** *"HLA-DQ2 and DQ8 are the genetic bridge between Type 1 Diabetes and Celiac Disease. DQ2 confers an odds ratio of 7.0 for celiac -- the strongest HLA association outside of HLA-B27 and ankylosing spondylitis. When you see both conditions in one patient, the HLA typing explains why."*

### Tab 4: Diagnostic Odyssey

1. Navigate to Diagnostic Odyssey tab
2. **Load James Cooper's timeline:**
   - DKA presentation (acute)
   - Autoantibody testing confirming T1D
   - Routine celiac screening (recommended for all new T1D diagnoses)
   - GI consultation and duodenal biopsy
   - HLA typing confirming DQ2/DQ8
   - Dietitian consultation for combined GFD + T1D management
3. **Show overlap syndrome detection:**
   - t1d_celiac_overlap syndrome detected
   - Shared markers: anti-tTG IgA, HLA-DQ2, HLA-DQ8
   - Prevalence: 5-10% of T1D patients have celiac disease

**Key talking point:** *"This is the overlap syndrome demo. The agent automatically detects the T1D-celiac overlap based on shared antibody markers and HLA typing, and flags the 5-10% co-occurrence rate. In clinical practice, this should trigger automatic celiac screening for every newly diagnosed T1D patient -- and this agent enforces that screening protocol."*

### Tab 8: Flare Prediction

1. Navigate to Flare Prediction tab
2. **Show T1D flare patterns:**
   - Early warning: HbA1c rising, C-peptide falling, GAD65 antibody titer
   - Protective: stable HbA1c, preserved C-peptide, insulin dose stable, time-in-range >70%
3. **Show celiac flare patterns:**
   - Early warning: anti-tTG rising (>20 U/mL), ferritin falling (<15 ng/mL)
   - Protective: normal tTG, strict gluten-free diet, normal hemoglobin
4. **Cross-disease management complexity:** Managing both conditions simultaneously requires coordinating GFD compliance (celiac) with carbohydrate counting (T1D)

**Key talking point:** *"The challenge with autoimmune overlap syndromes is cross-disease management. A celiac flare with villous atrophy causes malabsorption that destabilizes glucose control in T1D. The agent monitors both disease patterns simultaneously and alerts when biomarker shifts in one disease threaten stability in the other."*

---

## Troubleshooting

| Issue | Resolution |
|---|---|
| Streamlit not starting | Check port 8531 is free: `lsof -i :8531` |
| "No collections found" | Run `python3 scripts/setup_collections.py` |
| Empty evidence results | Check Milvus is running: `curl localhost:19530/healthz` |
| RAG returns empty | Verify Anthropic API key: `echo $AUTO_ANTHROPIC_API_KEY` |
| PDF ingestion fails | Ensure PyPDF2 is installed: `pip install PyPDF2` |
| FHIR validation errors | Verify export.py Patient resource in bundle |
| HLA lookup returns nothing | Use format `HLA-DRB1*04:01` (with HLA- prefix) |
| Flare prediction shows base risk only | Ensure biomarker values are provided (CRP, ESR, complement) |
| No therapy recommendations | Ensure at least one diagnosed_condition is set |

---

## Quick Reset

If you need to reset the demo environment:

```bash
# Re-create all collections (drops and recreates)
python3 scripts/setup_collections.py

# Restart API
uvicorn api.main:app --host 0.0.0.0 --port 8532 &

# Re-ingest demo data
curl -X POST http://localhost:8532/ingest/demo-data

# Restart UI
# Ctrl+C the streamlit process, then:
streamlit run app/autoimmune_ui.py --server.port 8531
```
