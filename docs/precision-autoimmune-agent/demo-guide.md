# Precision Autoimmune Intelligence Agent -- Demo Guide

**Author:** Adam Jones
**Date:** March 2026

> Step-by-step walkthrough of all 10 Streamlit UI tabs using the 9 demo patients. Use this guide for live demos, recorded walkthroughs, and self-guided exploration.

---

## Prerequisites

Before starting the demo, ensure:

1. **Milvus** is running on `localhost:19530` with 14 autoimmune collections created
2. **Streamlit UI** is running on `http://localhost:8531`
3. **FastAPI server** is running on `http://localhost:8532`
4. **Anthropic API key** is configured (for RAG queries in Tab 1 and patient analysis in Tab 2)
5. **Demo patient data** has been generated and ingested

Verify everything is ready:

```bash
cd ai_agent_adds/precision_autoimmune_agent

# Create collections and seed knowledge
python3 scripts/setup_collections.py --seed

# Generate and ingest 9 demo patients
python3 scripts/generate_demo_patients.py

# Run tests to verify
python3 -m pytest tests/ -v
# -> 455 tests passed
```

### Quick Health Check

```bash
# Check Milvus
curl -s http://localhost:19530/healthz
# Expected: {"status":"OK"}

# Check FastAPI
curl -s http://localhost:8532/healthz
# Expected: {"status":"healthy","collections":14,"milvus":"connected"}

# Check Streamlit
curl -s http://localhost:8531/_stcore/health
# Expected: "ok"
```

If any service is not responding, see the Troubleshooting section at the end of this guide.

---

## Demo Flow Overview

| Tab | Time | Key Talking Points |
|---|---|---|
| 1. Clinical Query | 3 min | RAG search across 14 collections, streaming responses, citation scoring |
| 2. Patient Analysis | 4 min | Full pipeline: autoantibodies + HLA + activity + flare + therapy |
| 3. Document Ingest | 2 min | PDF upload, chunking, embedding, collection storage |
| 4. Diagnostic Odyssey | 3 min | Timeline visualization, delays, pattern recognition |
| 5. Autoantibody Panel | 3 min | Panel interpretation, disease associations, titer analysis |
| 6. HLA Analysis | 2 min | Allele-disease mapping, odds ratios, population context |
| 7. Disease Activity | 3 min | 20 scoring systems, threshold-based classification |
| 8. Flare Prediction | 3 min | Biomarker patterns, risk stratification, protective factors |
| 9. Therapy Advisor | 3 min | 22 biologics, PGx filtering, contraindications |
| 10. Knowledge Base | 1 min | Collection stats, knowledge version, coverage summary |
| **Total** | **~27 min** | |

### Recommended Demo Paths

For a **15-minute executive demo**, focus on:
- Tab 1 (Clinical Query) with Sarah Mitchell -- 3 min
- Tab 2 (Patient Analysis) with Sarah Mitchell -- 4 min
- Tab 4 (Diagnostic Odyssey) with Sarah Mitchell -- 3 min
- Tab 9 (Therapy Advisor) with Rachel Thompson -- 3 min
- Tab 10 (Knowledge Base) -- 2 min

For a **10-minute technical demo**, focus on:
- Tab 1 (Clinical Query) -- show streaming RAG -- 2 min
- Tab 3 (Document Ingest) -- show PDF pipeline -- 2 min
- Tab 5 (Autoantibody Panel) with James Cooper -- 2 min
- Tab 6 (HLA Analysis) with David Park -- 2 min
- Tab 10 (Knowledge Base) -- collection stats -- 2 min

---

## Patient 1: Sarah Mitchell -- SLE (Systemic Lupus Erythematosus)

### Key Narrative

> A 34-year-old woman with a 3-year diagnostic odyssey from initial fatigue and joint pain through multiple specialists before receiving a definitive SLE diagnosis. Positive ANA 1:640 homogeneous, anti-dsDNA positive, low C3/C4 complement. HLA-DRB1*03:01 carrier. Currently on hydroxychloroquine and mycophenolate with history of lupus nephritis (Class IV). 27 clinical documents spanning 2022-2025.

### Tab 1: Clinical Query

**Setup:** Navigate to Tab 1 (Clinical Query).

**Demo steps:**

1. Type the query: `What is Sarah Mitchell's current disease status and flare risk?`
2. Click Submit and watch the streaming response appear token by token
3. **Point out key elements in the response:**
   - Evidence sources from multiple collections (clinical documents, lab results, autoantibody panels)
   - Citation confidence levels displayed as badges (High, Medium, Low)
   - The response references specific lab values from ingested documents
   - The answer is grounded in domain-specific knowledge, not general LLM output
4. Type a follow-up query: `Compare rituximab vs belimumab for Sarah's refractory lupus nephritis`
5. **Show collection sources** in the evidence panel:
   - `autoimmune_biologic_therapies` -- drug mechanism and indications
   - `autoimmune_clinical_documents` -- Sarah's actual clinical history
   - `autoimmune_literature` -- published evidence for rituximab vs belimumab
   - `autoimmune_pgx_rules` -- pharmacogenomic considerations

**Expected output highlights:**
- Streaming response with 3-8 evidence citations
- Citations labeled by collection (e.g., "[Biologic Therapy: rituximab]", "[Clinical Document: sarah_mitchell_progress_note_2024-03]")
- Relevance scores for each citation (High >= 0.80, Medium >= 0.60)

**Key talking point:** *"The RAG engine searches 14 specialized collections in parallel -- from autoantibody panels to biologic therapy databases to published literature -- and synthesizes a grounded answer with citations. This is not a general-purpose chatbot; every response is anchored in domain-specific knowledge."*

### Tab 2: Patient Analysis

**Setup:** Navigate to Tab 2 (Patient Analysis).

**Demo steps:**

1. Select **Sarah Mitchell** from the patient dropdown
2. Click **Run Full Analysis**
3. Wait for the 5-step pipeline to complete (typically 10-30 seconds)
4. Walk through each section of the output:

   **Section 1: Autoantibody Interpretation**
   - ANA: 1:640, homogeneous pattern -- positive
   - Anti-dsDNA: positive -- highly specific for SLE (specificity ~95%)
   - Anti-Smith: positive -- pathognomonic for SLE
   - Show how each antibody maps to disease associations with sensitivity/specificity data

   **Section 2: HLA Associations**
   - HLA-DRB1*03:01 -- OR 2.4 for SLE, OR 3.1 for Sjogren's
   - Explain the concept of odds ratios and multi-disease risk

   **Section 3: Disease Activity**
   - SLEDAI-2K score calculation
   - Activity level classification (remission/low/moderate/high/very high)
   - Score components breakdown

   **Section 4: Flare Prediction**
   - Risk score (0.0-1.0) with classification (imminent/high/moderate/low)
   - Contributing factors (rising anti-dsDNA, falling complement C3)
   - Protective factors (stable hydroxychloroquine adherence)

   **Section 5: Biologic Therapy Recommendations**
   - Ranked list of indicated biologics
   - PGx considerations for each drug
   - Contraindications and required monitoring

5. **Point out cross-agent findings** -- How biomarker results from the Biomarker Agent can feed into autoimmune analysis

**Expected output highlights:**
- Complete analysis result with all 5 sections populated
- Critical alerts highlighted (if disease activity is high or flare risk is imminent)
- Therapy recommendations with PGx context

**Key talking point:** *"In under 30 seconds, we go from raw clinical data to a complete autoimmune intelligence report -- autoantibody interpretation, HLA risk analysis, disease activity scoring, flare prediction, and therapy recommendations. A rheumatologist would need hours to synthesize this from scattered records."*

### Tab 3: Document Ingest

**Setup:** Navigate to Tab 3 (Document Ingest).

**Demo steps:**

1. Show the existing document list for Sarah Mitchell -- 27 clinical PDFs already ingested
2. **Explain the ingestion pipeline** by walking through the stages:
   - PDF upload -> page-by-page text extraction (PyPDF2)
   - Text chunking: 2,500 character chunks with 200 character overlap
   - Document type classification (7 types: lab report, progress note, imaging, pathology, genetic, referral, medication list)
   - Medical specialty detection (11 specialties)
   - Entity extraction: autoantibody names (29 patterns), lab test values (45 patterns)
   - BGE-small-en-v1.5 embedding (384 dimensions)
   - Milvus vector insert into `autoimmune_clinical_documents` collection
3. Point out the document type breakdown for Sarah:
   - Progress notes (rheumatology, nephrology, dermatology)
   - Lab reports (autoantibody panels, CBC, CMP, complement levels)
   - Imaging reports (chest X-ray, renal ultrasound)
   - Pathology (renal biopsy -- Class IV lupus nephritis)
   - Referral letters
4. Optionally upload a new test PDF to demonstrate the live pipeline

**Expected output highlights:**
- Document list with type classification, specialty, provider, and date
- Chunk count per document (typical: 3-8 chunks per PDF page)
- Successful embedding and storage confirmation

**Key talking point:** *"Every clinical document is chunked, embedded, and stored as searchable vectors. When you ask a question, the RAG engine can pull relevant passages from any document in the patient's history -- not just structured data, but the actual clinical narrative."*

### Tab 4: Diagnostic Odyssey

**Setup:** Navigate to Tab 4 (Diagnostic Odyssey).

**Demo steps:**

1. Select **Sarah Mitchell** from the patient dropdown
2. **Show the timeline visualization** -- a chronological view of the diagnostic journey from first symptoms (March 2022) to SLE diagnosis (May 2023)
3. Walk through the key events on the timeline:
   - **2022-03:** Annual physical -- fatigue dismissed as stress (symptom onset)
   - **2022-06:** Follow-up -- joint pain reported, CBC/CMP ordered (first labs)
   - **2022-10:** ANA ordered after persistent symptoms (first autoimmune test)
   - **2023-01:** Dermatology referral for malar rash (referral)
   - **2023-04:** Rheumatology referral -- full autoimmune panel ordered (turning point)
   - **2023-05:** Definitive SLE diagnosis -- ANA 1:640, anti-dsDNA+, anti-Smith+ (diagnosis)
   - **2023-09:** Nephrology consult -- proteinuria detected (complication)
   - **2023-10:** Renal biopsy -- Class IV lupus nephritis (procedure)
4. **Highlight the diagnostic delay** -- 14 months from first symptom to correct diagnosis
5. **Show event classification** -- each event is classified by type (symptom onset, diagnosis, misdiagnosis, lab result, imaging, biopsy, referral, ER visit, treatment start, treatment change, flare)
6. **Show the specialists seen** -- primary care, dermatology, rheumatology, nephrology
7. Point out the pattern: initial dismissal -> delayed testing -> specialist referral -> diagnosis

**Expected output highlights:**
- Chronological timeline with color-coded event types
- Days-from-first-symptom calculation for each event
- Total diagnostic delay prominently displayed

**Key talking point:** *"The diagnostic odyssey analyzer reveals the full timeline of Sarah's journey -- 14 months from first symptoms to diagnosis, 4 different specialists, multiple dismissals. The average lupus patient sees 4 physicians over 4 years before diagnosis. This tool makes that delay visible and actionable."*

### Tab 5: Autoantibody Panel

**Setup:** Navigate to Tab 5 (Autoantibody Panel).

**Demo steps:**

1. Select **Sarah Mitchell** from the patient dropdown
2. **Show Sarah's autoantibody panel:**
   - ANA: 1:640, homogeneous pattern -- positive
   - Anti-dsDNA: positive -- highly specific for SLE (specificity ~95%)
   - Anti-Smith: positive -- pathognomonic for SLE (specificity ~99%)
   - C3: low (48 mg/dL, normal 90-180) -- complement consumption
   - C4: low (6 mg/dL, normal 16-47) -- complement consumption
3. **Show disease association scoring:**
   - Each positive antibody maps to diseases with sensitivity and specificity data
   - The combination score shows high confidence for SLE
   - Individual antibody associations are displayed with evidence quality
4. **Contrast with ANA alone:**
   - ANA 1:640 homogeneous is found in 5-15% of healthy individuals
   - The specificity comes from the combination: ANA + anti-dsDNA + anti-Smith + low complements

**Expected output highlights:**
- Antibody results table with positive/negative status
- Disease association matrix showing which diseases each antibody is associated with
- Combined confidence score for the most likely diagnosis

**Key talking point:** *"ANA alone is found in 5-15% of healthy individuals. But ANA 1:640 homogeneous plus anti-dsDNA plus anti-Smith plus low complements -- that pattern is pathognomonic for SLE. The autoantibody interpreter understands these combinations, not just individual results."*

### Tab 6: HLA Analysis

**Setup:** Navigate to Tab 6 (HLA Analysis).

**Demo steps:**

1. Select **Sarah Mitchell** from the patient dropdown
2. **Show Sarah's HLA profile:**
   - HLA-DRB1*03:01 -- OR 2.4 for SLE, OR 3.1 for Sjogren's
3. **Point out multi-disease risk:**
   - This single allele confers risk for SLE, Sjogren's, T1D, Graves, and Celiac
   - Explain that this is why autoimmune patients often develop multiple conditions
4. **Show the odds ratio context:**
   - Compare Sarah's HLA-DRB1*03:01 (OR 2.4 for SLE) to HLA-B*27:05 (OR 87.4 for AS)
   - Explain the wide range of HLA effect sizes across diseases
5. **Discuss population context:**
   - HLA allele frequencies vary by ethnic background
   - The odds ratios in the knowledge base are derived from published studies

**Expected output highlights:**
- HLA allele table with associated diseases and odds ratios
- Multi-disease risk profile visualization
- Population frequency context

**Key talking point:** *"HLA-DRB1*03:01 is a shared risk allele for SLE, Sjogren's, Type 1 Diabetes, Graves, and Celiac. This is why autoimmune patients often develop multiple conditions -- they carry genetic risk across several diseases."*

---

## Patient 2: Rachel Thompson -- Mixed Connective Tissue Disease (MCTD)

### Key Narrative

> A 38-year-old woman with Mixed Connective Tissue Disease (MCTD), presenting with overlapping features of SLE, SSc, and myositis. Anti-U1 RNP strongly positive, Raynaud's phenomenon, swollen fingers, and inflammatory myopathy. Being managed with immunosuppressive therapy.

### Tab 2: Patient Analysis (Rachel Thompson)

**Setup:** Navigate to Tab 2 (Patient Analysis).

**Demo steps:**

1. Select **Rachel Thompson** from the patient dropdown
2. Click **Run Full Analysis**
3. **Highlight the overlap syndrome detection:**
   - The agent identifies MCTD as an overlap of SLE, SSc, and myositis features
   - Anti-U1 RNP is the hallmark antibody for MCTD
   - Disease activity must be tracked across multiple organ systems
4. **Compare with Sarah's SLE analysis:**
   - Sarah has a single disease (SLE) with clear classification criteria
   - Rachel's MCTD requires evaluating criteria from multiple diseases simultaneously

**Expected output highlights:**
- Overlap syndrome identification (MCTD = SLE + SSc + myositis features)
- Multi-system disease activity assessment
- Therapy recommendations addressing the predominant features

**Key talking point:** *"Rachel's case demonstrates why overlap syndromes are the hardest challenge in autoimmune medicine. The agent evaluates features across SLE, SSc, and myositis simultaneously and recommends therapy based on which features predominate."*

### Tab 7: Disease Activity

**Setup:** Navigate to Tab 7 (Disease Activity).

**Demo steps:**

1. Select **Rachel Thompson** from the patient dropdown
2. **Show disease activity scores:**
   - MCTD activity assessment with component scores for each overlap domain
   - Joint involvement score (from SLE/RA criteria)
   - Myositis markers (CK, aldolase) contributing to myopathy assessment
   - Skin thickening score (from SSc modified Rodnan)
   - Serositis assessment
3. **Compare scoring systems:**
   - Show how SLEDAI-2K captures the lupus features
   - Show how the modified Rodnan skin score captures SSc features
   - Show how CK/aldolase track myositis activity
4. **Point out the treatment monitoring complexity:**
   - MCTD requires simultaneous monitoring of multiple organ systems
   - Different scoring systems are needed for different disease components

**Expected output highlights:**
- Multi-system disease activity breakdown
- Activity level classification per domain
- Composite assessment with treatment implications

**Key talking point:** *"Rachel's MCTD demonstrates why overlap syndromes are challenging -- she has features spanning SLE, SSc, and myositis. The disease activity scorer evaluates each component and tracks the overall disease trajectory."*

### Tab 8: Flare Prediction

**Setup:** Navigate to Tab 8 (Flare Prediction).

**Demo steps:**

1. Select **Rachel Thompson** from the patient dropdown
2. **Show the flare risk assessment:**
   - Overall risk score (0.0-1.0)
   - Risk classification (imminent >= 0.8, high >= 0.6, moderate >= 0.4, low < 0.4)
   - 90-day prediction window
3. **Walk through contributing factors:**
   - Rising CRP with changing anti-U1 RNP titers
   - New organ involvement (e.g., worsening Raynaud's)
   - Recent immunosuppression changes
4. **Show protective factors:**
   - Stable medication adherence
   - Normal inflammatory markers in other domains
5. **Discuss monitoring recommendations:**
   - Frequency of lab draws based on risk level
   - Clinical assessments to schedule
   - Warning signs for patients

**Expected output highlights:**
- Risk score with classification badge
- Contributing factors with individual risk contributions
- Protective factors with stabilization evidence
- Monitoring frequency recommendations

**Key talking point:** *"Flare prediction is not about a single lab value -- it is about patterns. Rising CRP with changing anti-U1 RNP titers and new organ involvement suggests an imminent flare. The predictor identifies these patterns across 13 biomarker signals per disease."*

### Tab 9: Therapy Advisor

**Setup:** Navigate to Tab 9 (Therapy Advisor).

**Demo steps:**

1. Select **Rachel Thompson** from the patient dropdown
2. **Show biologic recommendations for MCTD:**
   - Immunosuppressive therapies targeting overlap features
   - Organ-specific treatment considerations
   - Biologic options ranked by evidence and predominant features
3. **Walk through each recommendation:**
   - Drug name, class, and mechanism of action
   - Indicated diseases and evidence level
   - PGx considerations (e.g., CYP3A4 for JAK inhibitors, FCGR3A for rituximab)
   - Contraindications (active infections, hepatitis B, pregnancy)
   - Required monitoring (TB screening, hepatitis B serology, CBC, LFTs)
4. **Show the PGx filtering:**
   - If Rachel has relevant pharmacogenomic variants, show how they modify recommendations
   - Explain that PGx data is integrated from the Biomarker Agent

**Expected output highlights:**
- Ranked therapy list with evidence levels
- PGx considerations per drug
- Contraindication flags
- Monitoring requirements table

**Key talking point:** *"Rachel's MCTD requires a nuanced therapy approach -- the advisor considers which overlap features predominate and recommends treatments accordingly. Each recommendation includes PGx considerations, contraindications, and required monitoring."*

---

## Patient 3: James Cooper -- T1D + Celiac Overlap

### Key Narrative

> A 19-year-old male with Type 1 Diabetes, now presenting with GI symptoms and positive anti-tTG antibodies. HLA-DQB1*02:01 carrier (shared risk for both T1D and Celiac). Demonstrates autoimmune polyendocrine overlap syndrome.

### Tab 5: Autoantibody Panel (James Cooper)

**Setup:** Navigate to Tab 5 (Autoantibody Panel).

**Demo steps:**

1. Select **James Cooper** from the patient dropdown
2. **Show James's autoantibody panel:**
   - Anti-GAD65: positive -- highly specific for T1D (specificity ~98%)
   - Anti-IA2: positive -- T1D-specific islet cell antibody
   - Anti-ZnT8: positive -- zinc transporter antibody, T1D marker
   - Anti-tTG IgA: elevated -- tissue transglutaminase, celiac marker
   - Anti-EMA: positive -- endomysial antibody, celiac confirmation
3. **Show disease association scoring:**
   - Anti-GAD65 + anti-IA2 + anti-ZnT8 map strongly to T1D
   - Anti-tTG IgA + anti-EMA map strongly to Celiac Disease
   - The agent recognizes this as a T1D-Celiac overlap pattern
4. **Contrast with Sarah's single-disease panel:**
   - Sarah's antibodies all point to one disease (SLE)
   - James's antibodies span two distinct autoimmune diseases

**Expected output highlights:**
- Dual disease association mapping
- Overlap syndrome identification
- Individual antibody sensitivity/specificity data

**Key talking point:** *"James has autoantibodies spanning two diseases -- anti-GAD65 and anti-IA2 for Type 1 Diabetes, plus anti-tTG and anti-EMA for Celiac. This is not coincidence. HLA-DQB1*02:01 confers shared risk for both conditions, and up to 10% of T1D patients develop Celiac. The autoimmune agent recognizes these overlap patterns."*

### Tab 6: HLA Analysis (James Cooper)

**Setup:** Navigate to Tab 6 (HLA Analysis).

**Demo steps:**

1. Select **James Cooper** from the patient dropdown
2. **Show James's HLA profile:**
   - HLA-DQB1*02:01 -- OR 7.0 for Celiac, OR 3.0 for T1D
   - HLA-DRB1*03:01 -- shared risk for T1D, SLE, Graves, Celiac
3. **Explain the HLA-DQ2 heterodimer:**
   - HLA-DQB1*02:01 forms part of the HLA-DQ2 molecule
   - HLA-DQ2 is the strongest known genetic risk factor for celiac disease
   - Nearly all celiac patients carry HLA-DQ2 or HLA-DQ8
4. **Cross-disease risk profile:**
   - Show the polyautoimmune risk assessment
   - Discuss monitoring recommendations for additional autoimmune conditions

**Expected output highlights:**
- HLA allele table with dual disease associations
- Odds ratios for both T1D and Celiac
- Cross-disease risk visualization

**Key talking point:** *"HLA-DQB1*02:01 forms part of the HLA-DQ2 heterodimer -- the strongest genetic risk factor for celiac disease. The same allele also confers significant T1D risk. This is why autoimmune overlap syndromes are not random -- they are genetically programmed."*

---

## Patient 4: David Park -- Ankylosing Spondylitis

### Key Narrative

> A 45-year-old male with ankylosing spondylitis (AS), diagnosed after a 3-year diagnostic odyssey. HLA-B*27:05 carrier (OR 87.4 for AS). Presents with chronic lower back pain, sacroiliitis on imaging, and recurrent anterior uveitis. Currently on a TNF inhibitor.

### Tab 6: HLA Analysis (David Park)

**Setup:** Navigate to Tab 6 (HLA Analysis).

**Demo steps:**

1. Select **David Park** from the patient dropdown
2. **Show David's HLA profile:**
   - HLA-B*27:05 -- OR 87.4 for Ankylosing Spondylitis
   - This is the strongest HLA-disease association in the entire knowledge base
3. **Compare odds ratios across diseases:**
   - HLA-B*27:05 for AS: OR 87.4 (extremely strong)
   - HLA-DRB1*03:01 for SLE: OR 2.4 (moderate)
   - HLA-DQB1*02:01 for Celiac: OR 7.0 (strong)
4. **Discuss clinical implications:**
   - Over 90% of AS patients are HLA-B*27 positive
   - But only 5-10% of HLA-B*27 carriers develop AS
   - HLA-B*27 testing is a key part of the AS diagnostic workup

**Expected output highlights:**
- HLA-B*27:05 prominently displayed with highest odds ratio
- Disease association confidence metrics
- Broad allele matching (B*27:05 maps to B*27 group)

**Key talking point:** *"HLA-B*27:05 has the strongest disease association in our knowledge base -- an odds ratio of 87.4 for ankylosing spondylitis. Over 90% of AS patients carry this allele. Yet only 5-10% of carriers develop AS, showing that HLA is necessary but not sufficient."*

### Tab 7: Disease Activity (David Park)

**Setup:** Navigate to Tab 7 (Disease Activity).

**Demo steps:**

1. Select **David Park** from the patient dropdown
2. **Show BASDAI scoring:**
   - Bath Ankylosing Spondylitis Disease Activity Index
   - Components: fatigue, spinal pain, peripheral pain, enthesitis, morning stiffness severity, morning stiffness duration
   - Score interpretation: < 4 = low activity, >= 4 = high activity
3. **Discuss treatment decisions based on BASDAI:**
   - BASDAI >= 4 on two occasions 12 weeks apart -> eligible for biologic therapy
   - This is the ASAS/EULAR recommendation for biologic initiation

**Expected output highlights:**
- BASDAI score with component breakdown
- Activity level classification
- Treatment eligibility assessment

### Tab 4: Diagnostic Odyssey (David Park)

**Setup:** Navigate to Tab 4 (Diagnostic Odyssey).

**Demo steps:**

1. Select **David Park** from the patient dropdown
2. **Show the 3-year diagnostic journey:**
   - Initial presentation with lower back pain to primary care
   - Misattributed to mechanical back pain / disc disease
   - Physical therapy without improvement
   - Uveitis episode triggers immunology referral
   - HLA-B*27 testing and sacroiliac joint imaging
   - Definitive AS diagnosis by rheumatology
3. **Highlight the pattern:** inflammatory back pain misdiagnosed as mechanical

**Key talking point:** *"David's 3-year odyssey from first symptoms to AS diagnosis is typical -- inflammatory back pain in young men is frequently misattributed to mechanical causes. The diagnostic odyssey tool makes this pattern visible."*

---

## Patient 5: Linda Chen -- Sjogren's Syndrome

### Key Narrative

> A 45-year-old woman with Sjogren's Syndrome presenting with sicca syndrome (dry eyes, dry mouth), positive anti-SSA/Ro and anti-SSB/La antibodies, abnormal Schirmer test, and parotid gland involvement. HLA-DRB1*03:01 carrier.

### Tab 5: Autoantibody Panel (Linda Chen)

**Setup:** Navigate to Tab 5 (Autoantibody Panel).

**Demo steps:**

1. Select **Linda Chen** from the patient dropdown
2. **Show Linda's panel:**
   - Anti-SSA/Ro: positive -- the most characteristic antibody for Sjogren's (sensitivity ~70%)
   - Anti-SSB/La: positive -- highly specific for Sjogren's (specificity ~95%)
   - ANA: positive (1:320, speckled pattern)
   - RF (Rheumatoid Factor): positive
3. **Explain the SSA/SSB pattern:**
   - Anti-SSA alone is found in SLE, Sjogren's, and neonatal lupus
   - Anti-SSA + anti-SSB together is the hallmark of primary Sjogren's
   - The speckled ANA pattern is consistent with anti-SSA/SSB

**Key talking point:** *"Anti-SSA/Ro is not specific to Sjogren's -- it is also found in SLE. But anti-SSA plus anti-SSB together, with sicca symptoms and a speckled ANA pattern, is the hallmark of primary Sjogren's. The agent understands these patterns."*

---

## Patient 6: Maya Rodriguez -- POTS / hEDS / MCAS

### Key Narrative

> A 28-year-old woman with the POTS/hEDS/MCAS triad (Postural Orthostatic Tachycardia Syndrome, hypermobile Ehlers-Danlos Syndrome, Mast Cell Activation Syndrome). This case demonstrates a complex diagnostic odyssey through multiple specialties before the triad is recognized.

### Tab 4: Diagnostic Odyssey (Maya Rodriguez)

**Setup:** Navigate to Tab 4 (Diagnostic Odyssey).

**Demo steps:**

1. Select **Maya Rodriguez** from the patient dropdown
2. **Show the extended diagnostic odyssey:**
   - Years of unexplained symptoms: dizziness, fatigue, GI distress, joint hypermobility
   - Multiple specialists: cardiology (POTS evaluation), gastroenterology, rheumatology, allergy/immunology
   - Tilt table testing for POTS diagnosis
   - Beighton score for hEDS assessment
   - Mast cell mediator testing for MCAS
3. **Highlight the triad recognition:**
   - POTS, hEDS, and MCAS frequently co-occur
   - The agent's overlap syndrome detection identifies this as the POTS/hEDS/MCAS triad
   - Each component was diagnosed separately over months to years

**Key talking point:** *"Maya's case represents the POTS/hEDS/MCAS triad -- three conditions that frequently co-occur but are often diagnosed separately over years. The overlap syndrome detector recognizes this pattern and alerts clinicians to evaluate all three."*

---

## Patient 7: Emma Williams -- Multiple Sclerosis (RRMS)

### Key Narrative

> A 34-year-old woman with Relapsing-Remitting Multiple Sclerosis (RRMS). Presented with optic neuritis, MRI showing periventricular lesions, oligoclonal bands in CSF. HLA-DRB1*15:01 carrier (the strongest genetic risk factor for MS).

### Tab 7: Disease Activity (Emma Williams)

**Setup:** Navigate to Tab 7 (Disease Activity).

**Demo steps:**

1. Select **Emma Williams** from the patient dropdown
2. **Show EDSS scoring:**
   - Expanded Disability Status Scale (0-10)
   - Functional system scores (visual, brainstem, pyramidal, cerebellar, sensory, bowel/bladder, cerebral)
3. **Discuss DMT (disease-modifying therapy) selection:**
   - Based on disease activity level and EDSS trajectory
   - High-efficacy vs moderate-efficacy DMTs

**Key talking point:** *"Emma's RRMS requires continuous monitoring with the EDSS scale. The disease activity scorer tracks her functional system scores and helps guide disease-modifying therapy escalation decisions."*

---

## Patient 8: Karen Foster -- Systemic Sclerosis (dcSSc)

### Key Narrative

> A 48-year-old woman with diffuse cutaneous systemic sclerosis (dcSSc). Anti-Scl-70 (anti-topoisomerase I) positive, Raynaud's phenomenon, skin thickening, and interstitial lung disease (ILD).

### Tab 9: Therapy Advisor (Karen Foster)

**Setup:** Navigate to Tab 9 (Therapy Advisor).

**Demo steps:**

1. Select **Karen Foster** from the patient dropdown
2. **Show therapy recommendations for dcSSc:**
   - Immunosuppressive therapies for skin and lung involvement
   - Mycophenolate or cyclophosphamide for ILD
   - Nintedanib for progressive ILD
   - Rituximab as an alternative
3. **Discuss organ-specific monitoring:**
   - Pulmonary function tests (FVC, DLCO) for ILD tracking
   - Modified Rodnan skin score for skin progression
   - Echocardiography for pulmonary hypertension screening

**Key talking point:** *"Karen's dcSSc with ILD requires therapy that addresses both skin and lung disease. The advisor ranks treatments by evidence level and includes organ-specific monitoring protocols."*

---

## Patient 9: Michael Torres -- Graves' Disease

### Key Narrative

> A 41-year-old male with Graves' Disease. Positive thyroid-stimulating immunoglobulin (TSI), elevated free T4, suppressed TSH. HLA-DRB1*03:01 carrier.

### Tab 2: Patient Analysis (Michael Torres)

**Setup:** Navigate to Tab 2 (Patient Analysis).

**Demo steps:**

1. Select **Michael Torres** from the patient dropdown
2. Click **Run Full Analysis**
3. **Show Burch-Wartofsky scoring:**
   - Point system for thyroid storm assessment
   - Components: temperature, CNS effects, GI dysfunction, heart rate, heart failure, atrial fibrillation, precipitant history
   - Score interpretation: < 25 unlikely, 25-44 impending, >= 45 thyroid storm
4. **Show therapy recommendations:**
   - Antithyroid drugs (methimazole, PTU) vs radioactive iodine (RAI) vs thyroidectomy
   - PGx considerations for drug metabolism

**Key talking point:** *"Michael's Graves' Disease case demonstrates the Burch-Wartofsky scoring system for thyroid storm assessment -- a clinical emergency. The agent calculates the score and provides treatment recommendations including the antithyroid drug vs RAI decision framework."*

---

## Tab 10: Knowledge Base

**Setup:** Navigate to Tab 10 (Knowledge Base).

**Demo steps (use with any patient or standalone):**

1. **Show knowledge base version and statistics:**
   - Version: 2.0.0
   - 22 HLA alleles mapped with disease associations and odds ratios
   - 24 autoantibodies characterized with sensitivity/specificity data
   - 22 biologic therapies catalogued with PGx considerations
   - 20 disease activity scoring systems across 13 diseases
   - 13 flare prediction biomarker patterns
   - 10 ACR/EULAR classification criteria sets
   - 9 overlap syndrome patterns
2. **Show collection statistics:**
   - Vector counts per collection
   - Total vectors across all 14 collections
   - Index type and similarity metric
3. **Show source citations:**
   - ACR/EULAR classification criteria (with publication years)
   - CPIC pharmacogenomic guidelines
   - HLA disease association databases
   - FDA biologic drug approvals

**Key talking point:** *"The knowledge base is curated from validated clinical sources -- ACR/EULAR classification criteria, HLA disease association databases, CPIC pharmacogenomic guidelines, and FDA biologic approvals. Every recommendation traces back to published evidence."*

---

## Sample API Interactions

For technical audiences, demonstrate the REST API alongside the UI:

### Query via API

```bash
curl -X POST http://localhost:8532/v1/query \
  -H "Content-Type: application/json" \
  -d '{
    "question": "What biologics are recommended for refractory lupus nephritis?",
    "patient_id": "sarah_mitchell"
  }'
```

**Expected response:** JSON with `answer`, `evidence` (list of citations with scores), and `knowledge_used` fields.

### Patient Analysis via API

```bash
curl -X POST http://localhost:8532/v1/analyze \
  -H "Content-Type: application/json" \
  -d '{
    "patient_id": "sarah_mitchell",
    "include_autoantibodies": true,
    "include_hla": true,
    "include_activity": true,
    "include_flare": true,
    "include_therapy": true
  }'
```

**Expected response:** JSON with sections for autoantibody_interpretation, hla_associations, disease_activity, flare_prediction, and therapy_recommendations.

### Export via API

```bash
# FHIR R4 export
curl -X POST http://localhost:8532/v1/export/fhir \
  -H "Content-Type: application/json" \
  -d '{"patient_id": "sarah_mitchell"}' \
  -o sarah_mitchell_fhir.json

# PDF report
curl -X GET http://localhost:8532/v1/report/sarah_mitchell/pdf \
  -o sarah_mitchell_report.pdf
```

---

## Troubleshooting

### Common Issues

| Issue | Symptoms | Resolution |
|---|---|---|
| Streamlit not starting | Port 8531 not responding | Check port is free: `lsof -i :8531`. Kill conflicting process if needed. |
| "No collections found" | Empty results in all tabs | Run `python3 scripts/setup_collections.py --seed` to create and seed collections. |
| RAG returns empty responses | Query returns but with no citations | Check Milvus is running: `curl localhost:19530/healthz`. Verify collections are seeded. |
| Document ingest fails | Error on PDF upload | Check PDF file size (< 50 MB) and format. Ensure PyPDF2 can read the file. Corrupted PDFs will fail. |
| No demo patients available | Patient dropdown is empty | Run `python3 scripts/generate_demo_patients.py` to generate and ingest demo data. |
| LLM features unavailable | Queries return without synthesis | Verify `AUTO_ANTHROPIC_API_KEY` or `ANTHROPIC_API_KEY` environment variable is set with a valid key. |
| HLA analysis empty | No HLA associations shown | Ensure the selected patient has HLA profile data loaded. Not all patients have HLA genotyping. |
| Slow query responses | Queries take > 60 seconds | Check network latency to Anthropic API. Verify Milvus is not under memory pressure (`docker stats`). |
| Collection count mismatch | Fewer than 14 collections | Run `python3 scripts/setup_collections.py --seed` with `--drop-existing` flag to recreate all collections. |

### Milvus Connection Issues

```bash
# Check if Milvus container is running
docker ps | grep milvus

# Check Milvus logs for errors
docker logs milvus-standalone --tail 50

# Restart Milvus if needed
docker restart milvus-standalone

# Wait for Milvus to be ready (30-60 seconds)
until curl -s http://localhost:19530/healthz | grep -q OK; do
  echo "Waiting for Milvus..."
  sleep 5
done
echo "Milvus is ready"
```

### Anthropic API Issues

```bash
# Verify API key is set
echo $ANTHROPIC_API_KEY | head -c 10
# Should show: sk-ant-api

# Test API connectivity
curl -s https://api.anthropic.com/v1/messages \
  -H "x-api-key: $ANTHROPIC_API_KEY" \
  -H "content-type: application/json" \
  -H "anthropic-version: 2023-06-01" \
  -d '{"model":"claude-sonnet-4-6","max_tokens":10,"messages":[{"role":"user","content":"hi"}]}' \
  | python3 -m json.tool
```

### Performance Optimization

If demo performance is slow:

1. **Reduce collections searched:** Set `AUTO_TOP_K_PER_COLLECTION=3` (default 5)
2. **Lower evidence threshold:** Set `AUTO_SCORE_THRESHOLD=0.50` (default 0.40) to reduce low-relevance results
3. **Check Docker resource allocation:** Ensure Docker has at least 8 GB RAM allocated
4. **Monitor Milvus memory:** `docker stats milvus-standalone` -- Milvus should use < 4 GB for 14 collections

---

## Quick Reset

If you need to reset the demo environment completely:

```bash
cd ai_agent_adds/precision_autoimmune_agent

# Re-create collections and re-seed knowledge (drops and recreates)
python3 scripts/setup_collections.py --seed --drop-existing

# Re-generate demo patient data
python3 scripts/generate_demo_patients.py

# Verify all tests pass
python3 -m pytest tests/ -v
# -> 455 tests passed

# Restart Streamlit UI
# Ctrl+C the streamlit process, then:
streamlit run app/autoimmune_ui.py --server.port 8531

# Restart FastAPI server
# Ctrl+C the uvicorn process, then:
uvicorn api.main:app --host 0.0.0.0 --port 8532 --workers 2
```

### Docker Reset

```bash
# Stop all autoimmune agent containers
docker compose -f docker-compose.yml down

# Rebuild and restart
docker compose -f docker-compose.yml up -d --build

# Watch setup progress
docker compose logs -f autoimmune-setup

# Verify services are healthy
docker compose ps
```

---

!!! warning "Clinical Decision Support Disclaimer"
    The Precision Autoimmune Agent is a clinical decision support research tool for autoimmune disease assessment. It is not FDA-cleared and is not intended as a standalone diagnostic device. All recommendations should be reviewed by qualified healthcare professionals. Apache 2.0 License.
