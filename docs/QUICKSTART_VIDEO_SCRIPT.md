---
tags:
  - Video
  - Demo
  - Quickstart
---

# HCLS AI Factory -- 5-Minute Quickstart Video Script

Shot-by-shot storyboard for a screen recording demo of the HCLS AI Factory: an end-to-end precision medicine platform that transforms Patient DNA into Drug Candidates in under 5 hours on a single NVIDIA DGX Spark ($3,999). Apache 2.0 licensed, open-source. Created by Adam Jones (14+ years genomic research).

---

## Overview

| Field | Detail |
|---|---|
| **Duration** | 5:00 (five minutes) |
| **Format** | Screen recording with voice-over narration; picture-in-picture camera optional |
| **Audience** | Bioinformaticians, computational biologists, translational researchers, AI-in-healthcare engineers |
| **Tools needed** | OBS Studio or ScreenFlow, terminal (font size 18+), browser (Chrome/Firefox, 1920x1080), microphone |
| **Recording resolution** | 1920x1080 @ 30 fps, exported as MP4 (H.264) |
| **Tone** | Conversational but technical -- assume the viewer understands VCFs, variant calling, and molecular docking |

---

## Pre-Recording Checklist

- [ ] All services running (`./demo.sh --status` shows green across the board)
- [ ] Browser tabs pre-loaded: Landing (:8080), Genomics (:5000), Chat (:8501), Drug Discovery (:8505), CAR-T (:8521), Imaging (:8525), Oncology (:8526)
- [ ] Terminal font size 18+ with dark background for readability
- [ ] `.env.example` open in an editor tab (VS Code or nano) for the Setup section
- [ ] Chat UI history cleared so the VCP query is typed fresh on camera
- [ ] Drug Discovery UI reset to show a clean run for VCP/CB-5083
- [ ] Close all notifications, Slack, email -- clean desktop

---

## Section 1: Opening (0:00 -- 0:30)

| Timecode | Shot | Screen Content | Narration | Notes |
|---|---|---|---|---|
| 0:00 | 1 | Terminal, clean prompt. | "This is the HCLS AI Factory -- an open-source platform that takes patient DNA and produces ranked drug candidates in under five hours, on a single NVIDIA DGX Spark that costs three thousand nine hundred ninety-nine dollars." | Speak steadily. Let the price point land. |
| 0:07 | 2 | Type: `git clone https://github.com/ajones1923/hcls-ai-factory.git` -- show clone output scrolling. | "Everything is Apache 2.0. One git clone and you have the entire platform -- genomics, evidence reasoning, drug discovery, and three intelligence agents." | Let the clone finish or cut to completed output. |
| 0:18 | 3 | Type: `cd hcls-ai-factory && ls -la` -- show top-level directory listing with pipeline folders visible. | "Three pipeline stages, a Nextflow orchestrator, a landing page, monitoring, and full test suites. Let me show you how to get it running." | Hold on directory listing for 3 seconds so viewers can read folder names. |

**B-roll suggestion:** Quick cut to a photo of the DGX Spark hardware sitting on a desk, overlaid with the text "GB10 GPU / 128 GB unified LPDDR5x / 20 ARM cores / $3,999". Hold for 3 seconds.

**Transition:** Direct cut to editor/terminal.

---

## Section 2: Setup (0:30 -- 1:00)

| Timecode | Shot | Screen Content | Narration | Notes |
|---|---|---|---|---|
| 0:30 | 4 | Show `.env.example` in editor. Highlight the `NGC_API_KEY`, `ANTHROPIC_API_KEY`, and `NVIDIA_API_KEY` lines. | "Configuration is a single .env file. You need three API keys: NGC for Parabricks, Anthropic for Claude, and NVIDIA for the BioNeMo NIM endpoints -- MolMIM and DiffDock." | Zoom in on the key lines. Do not show real key values. |
| 0:40 | 5 | Terminal: type `cp .env.example .env` then `nano .env` (briefly show editing, then exit). | "Copy the example, drop in your keys, and you are configured. Everything else has sensible defaults -- Milvus on 19530, cloud NIM mode for ARM64 compatibility, Claude as the LLM provider." | Keep the nano view brief -- 3 seconds max. |
| 0:50 | 6 | Terminal: type `./setup-data.sh --status`. Show the status dashboard output with download progress or completion checkmarks. | "The setup script handles all data acquisition -- reference genomes, ClinVar, AlphaMissense, PDB structures. Run it with --status to see where you stand." | If data is already downloaded, the dashboard will show green checkmarks. That is ideal. |

**B-roll suggestion:** None needed. Fast-paced terminal work keeps attention.

**Transition:** Direct cut to `./demo.sh` execution.

---

## Section 3: Launch (1:00 -- 1:30)

| Timecode | Shot | Screen Content | Narration | Notes |
|---|---|---|---|---|
| 1:00 | 7 | Terminal: type `./demo.sh`. Show the NVIDIA-green ASCII banner printing, prerequisites check (Docker, Python, Ollama), and service startup sequence. | "One command launches everything. The demo script checks prerequisites, starts Milvus, spins up each pipeline UI, launches the intelligence agents, and opens the landing page." | Let the banner print fully -- it is visually distinctive. |
| 1:12 | 8 | Terminal continues: show services coming online -- "Waiting for Milvus... READY", "Waiting for RAG Streamlit... READY", etc. Green checkmarks appearing one by one. | "Each service gets a health check. When you see green across the board, you are ready to go. Full cold start takes about two to three minutes." | Hold on the checkmarks. Viewers will want to see every service come up. |
| 1:20 | 9 | Browser: Landing page at `localhost:8080`. Show the health grid with all services green, the three-stage pipeline diagram, and service links. | "The landing page gives you the full picture -- every service, its port, its health status. This is your control panel. All green. Let us walk through each stage." | Mouse over a few service tiles to show they are clickable. |

**B-roll suggestion:** Optional lower-third overlay listing all ports: "Landing=8080 / Genomics=5000 / RAG API=5001 / Chat=8501 / Drug Discovery=8505 / Portal=8510 / CAR-T=8521 / Imaging=8525 / Oncology=8526".

**Transition:** Click the Genomics tile on the landing page, or direct-cut to `:5000`.

---

## Section 4: Stage 1 -- GPU Genomics (1:30 -- 2:00)

| Timecode | Shot | Screen Content | Narration | Notes |
|---|---|---|---|---|
| 1:30 | 10 | Browser: Genomics Portal at `localhost:5000`. Show the HG002 sample info panel -- sample ID, input FASTQ size (~200 GB), reference genome. | "Stage one is GPU-accelerated genomics. We are using NVIDIA Parabricks 4.6 with BWA-MEM2 for alignment and DeepVariant for variant calling. The demo sample is HG002 from the Genome in a Bottle consortium -- about 200 gigabytes of paired-end whole-genome sequencing data." | Emphasize "GPU-accelerated" -- this is the core differentiator over CPU pipelines. |
| 1:40 | 11 | Show the pipeline output section: variant count (11.7M), runtime (120-240 min on DGX Spark), accuracy (>99% concordance). | "On the DGX Spark, alignment through variant calling completes in two to four hours. The output is 11.7 million variant calls at greater than 99 percent concordance. On a CPU cluster, this same work takes 24 to 48 hours." | Pause briefly after "24 to 48 hours" to let the comparison land. |
| 1:50 | 12 | Scroll to show VCF output summary, or show a terminal snippet of a VCF header with the DeepVariant version tag. | "The VCF flows directly into Stage 2 -- no manual handoff, no file conversion, no waiting. This is the point where most traditional workflows stop and hand off to a separate bioinformatics team." | This sets up the seamless pipeline narrative. |

**B-roll suggestion:** Side-by-side comparison graphic: "CPU Cluster: 24-48 hrs, $50K-500K+" vs. "DGX Spark: 2-4 hrs, $3,999". Hold for 4 seconds.

**Transition:** Click through to Chat UI or direct-cut to `:8501`.

---

## Section 5: Stage 2 -- RAG/Chat Evidence Engine (2:00 -- 3:00)

| Timecode | Shot | Screen Content | Narration | Notes |
|---|---|---|---|---|
| 2:00 | 13 | Browser: Chat UI at `localhost:8501`. Clean interface, empty chat history, input box visible. | "Stage two is where we turn variants into actionable targets. This is a retrieval-augmented generation system backed by Milvus with 3.56 million searchable vectors -- ClinVar variants, AlphaMissense pathogenicity predictions, and a curated knowledge base covering 201 genes across 13 therapeutic areas." | Speak at a measured pace. These numbers matter to the audience. |
| 2:15 | 14 | Type into the chat box: "What is known about VCP mutations in frontotemporal dementia? What variants are pathogenic and what makes VCP a druggable target?" Press enter. | "Let me query for VCP -- Valosin-Containing Protein -- in the context of frontotemporal dementia. This is a real research question with a known druggable target." | Type at a readable speed. Viewers will want to see the query. |
| 2:25 | 15 | Show Claude's response streaming in. The response should include: ClinVar variant matches with clinical significance, AlphaMissense pathogenicity scores, structural and functional context for VCP, and a druggability assessment. | "Claude synthesizes evidence from all three collections in under five seconds. You can see ClinVar hits with clinical significance ratings, AlphaMissense pathogenicity scores for specific missense variants, and a druggability assessment grounded in the structural data." | Let the response stream for several seconds. Do not rush past it. The quality of the synthesis is a key differentiator. |
| 2:42 | 16 | Scroll through the response to show specific variant details -- rsIDs, amino acid changes, significance classifications. Highlight any Milvus hit counts if displayed in the UI. | "Every claim is backed by vector-retrieved evidence. This is not a hallucinating chatbot -- it is a grounded reasoning engine. 85 percent of the 201 genes in the knowledge base have confirmed druggable targets, and VCP is one of them." | Point the cursor at specific evidence citations as you speak. |
| 2:55 | 17 | Optionally show a second query or show the sidebar with collection statistics (vector counts, embedding model info). | "You can drill deeper -- ask about specific variants, compare across diseases, explore structural implications. The system handles it conversationally, but every answer traces back to indexed evidence." | Keep this brief. The point is to show depth without burning time. |

**B-roll suggestion:** Lower-third overlay: "3.56M vectors / ClinVar ~2.7M / AlphaMissense 71M / 201 genes / 13 therapeutic areas / <5 sec query latency".

**Transition:** Narrate the handoff: "VCP is our target. Now let us generate drug candidates." Direct-cut to `:8505`.

---

## Section 6: Stage 3 -- Drug Discovery (3:00 -- 4:00)

| Timecode | Shot | Screen Content | Narration | Notes |
|---|---|---|---|---|
| 3:00 | 18 | Browser: Drug Discovery UI at `localhost:8505`. Show the target input panel with VCP selected, seed compound CB-5083 entered, PDB structures listed (5FTK, 8OOI, 9DIL, 7K56). | "Stage three takes the validated target and generates novel drug candidates. We are targeting VCP with CB-5083 as the seed compound -- a known VCP inhibitor. The system pulls crystal structures from RCSB PDB automatically." | Point to each PDB ID as you mention it. |
| 3:12 | 19 | Show the molecule generation step: MolMIM producing SMILES strings, generation progress indicator. | "MolMIM, one of NVIDIA's BioNeMo NIMs, generates structurally novel molecules seeded from CB-5083. It runs in cloud mode against the NVIDIA health API -- no local GPU container needed, which is critical for ARM64 compatibility on the DGX Spark." | Mention cloud NIM explicitly -- this is an architecture decision viewers will care about. |
| 3:25 | 20 | Show the docking step: DiffDock progress, confidence scores appearing for each molecule-protein pair. | "Each generated molecule gets docked against the VCP crystal structures using DiffDock. We run 10 poses per molecule and score them for binding affinity and geometric confidence." | Let a few docking scores appear on screen. |
| 3:40 | 21 | Show the ranked candidates table: columns for molecule ID, SMILES, docking score, QED, Lipinski pass/fail, composite score. Scroll through the top 10-20 candidates. | "The output is a ranked table of drug candidates scored on docking affinity, drug-likeness via QED, and Lipinski Rule of Five compliance. Our top VCP candidate shows a 39 percent composite improvement over the CB-5083 seed compound." | Slow down on "39 percent composite improvement" -- this is the headline result. |
| 3:52 | 22 | Optionally click into a top candidate to show its 2D structure visualization or 3D docking pose if the UI supports it. | "The entire drug discovery stage -- from target input to ranked candidates -- runs in 8 to 16 minutes. That is structure retrieval, molecule generation, 3D conformer creation, molecular docking, scoring, and report generation." | Summarize the sub-steps to reinforce the automation. |

**B-roll suggestion:** Split-screen showing the MolMIM SMILES output on the left and the ranked table on the right. Overlay text: "100 candidates / 10 poses each / 8-16 min total".

**Transition:** "Beyond the core pipeline, we have built three intelligence agents. Let me show you." Direct-cut to agent UIs.

---

## Section 7: Intelligence Agents (4:00 -- 4:30)

| Timecode | Shot | Screen Content | Narration | Notes |
|---|---|---|---|---|
| 4:00 | 23 | Browser: CAR-T Intelligence Agent at `localhost:8521`. Show the main interface with collection selector, evidence search, or a sample query result. | "The CAR-T Intelligence Agent is a specialized evidence engine for chimeric antigen receptor T-cell therapy. Ten dedicated collections, 6,266 vectors, and it handles everything from target antigen analysis to manufacturing protocol evidence. 241 tests, all passing." | Quick tour -- do not linger. 10 seconds max per agent. |
| 4:10 | 24 | Browser: Imaging Intelligence Agent at `localhost:8525`. Show the interface with NIM service indicators (VISTA-3D, MAISI, VILA-M3, Llama-3). | "The Imaging Intelligence Agent integrates four NVIDIA NIM microservices -- VISTA-3D for segmentation, MAISI for synthetic imaging, VILA-M3 for visual question answering, and Llama-3 for report generation. Ten collections, 539 tests. It produces FHIR R4 DiagnosticReports." | Mention FHIR R4 -- it signals clinical interoperability to the audience. |
| 4:20 | 25 | Browser: Precision Oncology Agent at `localhost:8526`. Show the MTB packet generation interface or a case summary view. | "The Precision Oncology Agent handles molecular tumor board workflows -- case creation, therapy ranking, clinical trial matching, and MTB packet generation. Eleven collections, 516 tests, FHIR R4 bundle export." | This is the newest agent. Keep it crisp. |

**B-roll suggestion:** Lower-third overlay table: "CAR-T: 10 collections, 241 tests / Imaging: 10 collections, 4 NIMs, 539 tests / Oncology: 11 collections, 516 tests / Total: 1,296 tests in 3.78 sec".

**Transition:** "Let me bring it all together." Direct-cut back to landing page.

---

## Section 8: Closing (4:30 -- 5:00)

| Timecode | Shot | Screen Content | Narration | Notes |
|---|---|---|---|---|
| 4:30 | 26 | Browser: Return to Landing page at `localhost:8080`. All services green. | "Here is the full picture. Every service healthy. Three pipeline stages, three intelligence agents, vector database, monitoring -- all running on a single workstation." | Let the green health grid fill the screen. |
| 4:38 | 27 | Hold on landing page. Overlay or narrate the key numbers. | "Let me put the numbers in context. End-to-end, this platform goes from raw FASTQ to ranked drug candidates in under five hours. The traditional approach takes 6 to 18 months -- that is a 99 percent reduction in time. The hardware is a $3,999 DGX Spark versus the $50,000 to $500,000 you would typically spend on cluster infrastructure and software licenses." | Speak deliberately. These are the numbers the viewer will remember. |
| 4:50 | 28 | Show the GitHub URL in the browser address bar or navigate to the GitHub repo page. | "171 druggable targets across 13 therapeutic areas. 3.56 million searchable vectors. 1,296 agent tests running in under 4 seconds. And every line of it is Apache 2.0 on GitHub." | Rattle off the stats with confidence. |
| 4:55 | 29 | Browser on GitHub repo page, or terminal with the repo URL displayed. Clean ending frame. | "Clone it, extend it, build on it. The link is in the description. I am Adam Jones -- thanks for watching." | End with a clean pause. Hold the frame for 3 seconds of silence before the video ends. |

**B-roll suggestion:** End card with GitHub URL (`github.com/ajones1923/hcls-ai-factory`), Apache 2.0 badge, and DGX Spark photo. Hold for 5 seconds.

**Transition:** Fade to black.

---

## Summary of Key Numbers (Reference for Narration)

Keep these numbers handy during recording. They appear throughout the script but are consolidated here for quick reference during retakes.

| Metric | Value |
|---|---|
| End-to-end time | < 5 hours |
| Traditional approach | 6-18 months |
| Time reduction | ~99% |
| Hardware cost | $3,999 (DGX Spark) |
| Traditional infrastructure | $50K-500K+ |
| GPU | NVIDIA GB10 Grace Blackwell |
| Memory | 128 GB unified LPDDR5x |
| CPU | 20 ARM cores (Grace) |
| Variant calls (Stage 1) | ~11.7 million |
| Stage 1 runtime | 120-240 min |
| Stage 1 accuracy | >99% concordance |
| CPU baseline for Stage 1 | 24-48 hours |
| Searchable vectors (Milvus) | 3.56 million |
| ClinVar records | ~2.7 million |
| AlphaMissense records | 71 million |
| Knowledge base genes | 201 across 13 therapeutic areas |
| Druggable targets | 171 (85% of 201) |
| Query latency (Stage 2) | < 5 seconds |
| Stage 3 runtime | 8-16 min |
| Top VCP candidate improvement | +39% composite over CB-5083 |
| CAR-T Agent tests | 241 |
| Imaging Agent tests | 539 |
| Precision Oncology Agent tests | 516 |
| Total agent tests | 1,296 in 3.78 sec |
| License | Apache 2.0 |

---

## Post-Production Notes

- **Intro/Outro cards:** Use the project banner from `docs/diagrams/hcls-ai-factory-diagram.png` as the opening card. GitHub repo URL and Apache 2.0 badge on the end card.
- **Lower thirds:** Use a consistent lower-third style for all metric overlays. White text on a semi-transparent dark bar. Keep font size large enough to read at 720p.
- **Music:** Optional low-energy ambient track under the narration. Keep it subtle -- the audience is technical and will tune out anything distracting. No music during the terminal and chat UI sections where the viewer is reading code.
- **Captions:** Generate closed captions. Many viewers in academic and clinical settings watch without sound.
- **Thumbnail:** DGX Spark hardware photo with overlay text: "Patient DNA to Drug Candidates / 5 Hours / $3,999".
- **Video description:** Include the GitHub URL, the three-stage pipeline summary, hardware specs, and a link to the full documentation.
- **Upload targets:** YouTube (primary), LinkedIn (native upload for better reach), project docs site (embedded).
