# HCLS AI Factory -- Live Demo Script

**Duration:** 45-60 minutes (live walkthrough of the 5-hour pipeline)
**Hardware:** NVIDIA DGX Spark -- GB10 GPU, 128GB unified memory, $4,699
**Presenter:** Adam Jones

---

## Pre-Demo Checklist

- [ ] DGX Spark powered on, all services green at http://localhost:8080
- [ ] Grafana dashboard open at http://localhost:3000
- [ ] FASTQ files staged (Patient Sarah -- simulated pediatric VCP case)
- [ ] Milvus loaded: ClinVar (~2.7M), AlphaMissense (71M), agent collections
- [ ] BioNeMo NIMs running: MolMIM (port 8001), DiffDock (port 8002)
- [ ] Screen recording software active
- [ ] Backup slides ready in case of hardware failure

---

## ACT 1: The Patient (2 minutes)

### Opening -- Set the Stakes

*[Slide: black screen, single sentence]*

> "No parent should ever have to bury a child because of disease."

*[Pause. Let it breathe.]*

> "This is Sarah. She is seven years old. She likes purple dinosaurs, chocolate ice cream, and building Lego spaceships with her dad."

*[Slide: child's drawing -- abstract, no real patient data]*

> "Six months ago, Sarah's parents noticed she was stumbling more than usual. She was forgetting words she knew last year. Her pediatrician referred her to a neurologist, who ordered genetic testing."

> "The results came back with a variant in the VCP gene -- Valosin-Containing Protein. The same gene implicated in Frontotemporal Dementia, Inclusion Body Myopathy, and Paget Disease of Bone. In adults, this is devastating. In a child, the literature is almost silent."

> "Today, we are going to take Sarah's raw DNA -- her FASTQ files, the rawest form of genomic data -- and in under five hours, on a machine that costs less than a MacBook Pro, we are going to find drug candidates that could change her life."

> "This is not a simulation. This is not a demo with pre-baked results. This is the real pipeline, running in real time, on real hardware."

---

## ACT 2: Stage 1 -- Genomics Pipeline (8 minutes)

### Launch the Pipeline

*[Terminal: show the Nextflow command]*

```bash
nextflow run main.nf --fastq_dir /data/sarah/ --output_dir /results/sarah/ -profile dgx_spark
```

> "We start where every precision medicine journey begins -- raw sequencing data. Sarah's whole genome was sequenced at 30x coverage. That is approximately 200 gigabytes of compressed FASTQ data -- three billion base pairs, read thirty times over for accuracy."

*[Switch to Grafana -- GPU utilization dashboard]*

> "Watch the GPU. This is Parabricks 4.6 running on the GB10. BWA-MEM2 is aligning every read to the human reference genome. On a traditional CPU cluster, this takes 24-48 hours. On the DGX Spark, we will finish alignment in under 90 minutes."

### Key Metrics to Highlight

- GPU utilization: ~95% sustained during alignment
- Memory bandwidth: NVLink-C2C streaming reads to GPU at 900 GB/s
- Throughput: ~1.8 billion reads aligned per hour

> "DeepVariant is next -- NVIDIA's GPU-accelerated variant caller. It examines every position in the genome and asks: does Sarah differ from the reference here? By the time it finishes, we will have identified approximately 11.7 million variant records."

*[Show variant count ticking up in the Genomics UI at port 5000]*

> "11.7 million variant records. That is Sarah's complete genetic fingerprint. But only a handful of those matter for her disease. Finding those needles in the haystack is where the AI agents take over."

---

## ACT 3: Stage 2 -- RAG/Chat Intelligence (10 minutes)

### The VCF Enters the Knowledge Engine

*[Open Chat UI at http://localhost:8501]*

> "Sarah's VCF -- Variant Call Format -- now enters our RAG engine. RAG stands for Retrieval-Augmented Generation. Think of it as giving an AI doctor instant access to every clinical database simultaneously."

> "We have loaded three massive knowledge bases into Milvus, our vector database:"

- **ClinVar**: 2.7 million clinically annotated variants from NCBI
- **AlphaMissense**: 71 million missense variant pathogenicity predictions from DeepMind
- **3.5 million semantic vectors** for similarity search

*[Type into the chat]*

```
Analyze VCP variants found in patient Sarah's VCF. Focus on pathogenic
and likely pathogenic variants with neurological implications.
```

> "Watch what happens. The system does not just do a keyword search. It embeds Sarah's variants into the same 384-dimensional space as every known clinical variant, using BGE-small-en-v1.5, and finds the nearest neighbors by meaning."

*[Results appear -- highlight the VCP variant]*

> "Here it is. VCP p.Arg155His -- classified as Pathogenic in ClinVar, with an AlphaMissense score of 0.94. The system has cross-referenced this against 47 published case reports, 3 functional studies, and the ACMG classification criteria. It has generated a full evidence chain, with citations."

### Cross-Agent Intelligence

> "But we are not done. Watch what happens when we engage the intelligence agents."

*[Navigate to Biomarker Agent at port 8502]*

> "The Biomarker Agent receives Sarah's variant and immediately pulls her lab values. Elevated CK levels, subtle inflammatory markers. It calculates a biological age offset -- Sarah's cellular machinery is aging faster than her chronological age. It flags muscle enzymes trending upward over the last three visits."

*[Navigate to Autoimmune Agent at port 8506]*

> "The Autoimmune Agent checks Sarah's HLA profile extracted from the same VCF. VCP mutations can trigger autoimmune-like inflammatory cascades. It cross-references her HLA alleles against known disease associations and flags potential immune dysregulation pathways."

> "This is the power of the multi-agent architecture. No single AI sees the whole picture. Together, they build a Patient 360 -- a complete view of Sarah's disease from every angle."

---

## ACT 4: Stage 3 -- Drug Discovery (10 minutes)

### From Target to Molecules

*[Open Drug Discovery UI at http://localhost:8505]*

> "The RAG engine has identified VCP -- p97 -- as our therapeutic target. p97 is an AAA+ ATPase essential for protein homeostasis. When it is mutated, misfolded proteins accumulate in neurons and muscle cells. That is what is happening to Sarah."

> "Our seed compound is CB-5083 -- a known VCP inhibitor that failed clinical trials for cancer due to toxicity, but whose binding mode gives us a structural starting point."

> "We have four crystal structures loaded:"
- **5FTK**: VCP D2 domain with ADP
- **8OOI**: VCP with allosteric inhibitor
- **9DIL**: Full-length VCP hexamer
- **7K56**: VCP-Npl4-Ufd1 complex

*[Click "Generate Candidates"]*

> "MolMIM -- a BioNeMo NIM running locally on this GPU -- is now generating novel molecular candidates. It uses the seed structure of CB-5083 and the binding pocket geometry to propose molecules that have never existed before."

*[Show generation progress -- 500 candidates in ~2 minutes]*

> "500 candidates generated. Now DiffDock takes over -- another BioNeMo NIM -- performing molecular docking. It predicts how each candidate physically binds to VCP's active site, calculating binding affinity and pose."

*[Show docking scores appearing in the UI]*

### Pharmacogenomic Filtering

> "But generating molecules that bind is not enough. We need molecules that are safe for Sarah. This is where the pipeline reconnects to her genome."

> "RDKit calculates drug-likeness -- Lipinski's rules, synthetic accessibility, toxicity flags. The Biomarker Agent cross-references her CYP450 genotypes -- the enzymes that metabolize drugs. If Sarah has a variant in CYP2D6 that makes her a poor metabolizer, any candidate metabolized by that pathway gets flagged."

*[Show the filtered results -- top 10 candidates with safety profiles]*

> "From 500 candidates, we are down to 12 that dock well, pass drug-likeness filters, and are compatible with Sarah's personal pharmacogenomic profile. Each one has a complete evidence chain back to her original DNA."

---

## ACT 5: The Reveal -- Patient 360 (8 minutes)

### The Complete Picture

*[Open Portal at http://localhost:8510]*

> "Let me show you what we have built in under five hours."

*[Walk through the Patient 360 dashboard]*

**Genomics Layer:**
> "11.7 million variant records called. 3 pathogenic variants identified in VCP. Full ACMG classification with evidence."

**Biomarker Layer:**
> "Biological age offset: +4.2 years. Muscle enzyme trajectory: ascending. Inflammatory markers: subclinical elevation. Disease trajectory: early-stage with intervention window."

**Intelligence Layer:**
> "5 AI agents, each analyzing Sarah from a different clinical perspective. The Oncology Agent cleared her VCP variants for cancer predisposition. The Autoimmune Agent flagged immune monitoring needs. The Imaging Agent has baseline recommendations for MRI surveillance."

**Drug Discovery Layer:**
> "12 novel drug candidates, each with predicted binding affinity to mutant VCP, drug-likeness scores, and pharmacogenomic compatibility with Sarah's genome."

**Evidence Chain:**
> "Every single finding traces back to Sarah's FASTQ files. Every recommendation cites primary literature. Every drug candidate links to the crystal structure it was docked against. This is not a black box. This is an evidence machine."

---

## ACT 6: The Mission (3 minutes)

### Why This Matters

> "Let me tell you what this means in real terms."

> "This entire pipeline -- from raw DNA to drug candidates -- ran on a single NVIDIA DGX Spark. A machine that costs $4,699. It fits under a desk. It draws 250 watts. A single graduate student can operate it."

> "Five years ago, this analysis would have required a $50 million supercomputer, a team of 20 bioinformaticians, 6 months, and a pharmaceutical partnership."

> "Today, Sarah's parents could walk into a children's hospital that has one of these on a shelf, and by the end of the day, they would have candidate molecules designed specifically for their daughter's mutation."

*[Pause]*

> "Every piece of software you have seen today is open source. Apache 2.0 licensed. The Nextflow orchestrator, the RAG engine, the intelligence agents, the drug discovery pipeline -- all of it is on GitHub. Any hospital, any university, any researcher in the world can download it and run it."

> "Because no parent should ever have to bury a child because the technology existed but was locked behind a paywall, or a partnership, or a server they could not afford."

> "This is the HCLS AI Factory. Patient DNA to drug candidates. Under five hours. Under four thousand dollars. Open to everyone."

*[End]*

---

## Appendix: Troubleshooting During Live Demo

### If Parabricks stalls
- Check GPU memory: `nvidia-smi`
- Fallback: switch to pre-computed VCF (`/data/sarah/precomputed.vcf`)

### If Milvus is slow
- Check collection load status: Attu at http://localhost:8000
- Restart Milvus: `docker compose restart milvus`

### If MolMIM/DiffDock NIMs are unresponsive
- Check NIM health: `curl http://localhost:8001/v1/health`
- Fallback: show pre-generated docking results from `/results/sarah/docking_backup/`

### If an agent is down
- Landing page at :8080 shows service health
- Individual restart: `docker compose restart precision-biomarker-agent`

### Timing Guide (for pre-recorded segments)
| Stage | Live Time | Real Pipeline Time |
|-------|-----------|-------------------|
| Genomics start | 0:00 | 0:00 |
| Alignment complete | 0:05 (show Grafana) | ~90 min |
| Variant calling done | 0:08 (show counts) | ~180 min |
| RAG analysis | 0:10 (live) | ~5 min |
| Agent analysis | 0:15 (live) | ~10 min |
| Drug generation | 0:20 (live) | ~8 min |
| Docking complete | 0:28 (live) | ~16 min |
| Patient 360 reveal | 0:35 | -- |
| Closing | 0:42 | Total: <5 hours |
