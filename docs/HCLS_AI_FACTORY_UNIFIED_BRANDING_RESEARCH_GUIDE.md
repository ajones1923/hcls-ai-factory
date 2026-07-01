# HCLS AI Factory — Unified Branding & Messaging Research Guide
## Strategic Framework for Global Launch at GTC Europe 2026

**Version:** 1.0.0
**Author:** Adam Jones
**Date:** March 2026
**Purpose:** Define unified branding, naming conventions, website architecture, and visual identity for the HCLS AI Factory platform across all touchpoints (DGX Spark, GitHub, Netlify, hcls-ai-factory.org)

---

## Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [The Branding Challenge](#2-the-branding-challenge)
3. [Unified Naming Framework](#3-unified-naming-framework)
4. [Persona-Based Messaging](#4-persona-based-messaging)
5. [Website Architecture — hcls-ai-factory.org](#5-website-architecture--hcls-ai-factoryorg)
6. [GitHub Repository Recommendations](#6-github-repository-recommendations)
7. [Visual Identity System](#7-visual-identity-system)
8. [Infographic Specifications — All 11 Agents](#8-infographic-specifications--all-11-agents)
9. [Netlify Deployment Optimization](#9-netlify-deployment-optimization)
10. [Cross-Platform Consistency](#10-cross-platform-consistency)
11. [Accessibility & Internationalization](#11-accessibility--internationalization)
12. [Product Requirements Document — Branding Implementation](#12-product-requirements-document--branding-implementation)
13. [Content Strategy & Voice Guide](#13-content-strategy--voice-guide)
14. [Social Media & Conference Playbook](#14-social-media--conference-playbook)
15. [GTC Europe 2026 Launch Plan](#15-gtc-europe-2026-launch-plan)
16. [Video & Demo Script Standards](#16-video--demo-script-standards)
17. [Competitive Positioning](#17-competitive-positioning)
18. [Community Engagement Framework](#18-community-engagement-framework)
19. [Metrics & KPIs for Brand Success](#19-metrics--kpis-for-brand-success)
20. [Legal & Licensing Brand Considerations](#20-legal--licensing-brand-considerations)
21. [Partner Co-Branding Guidelines](#21-partner-co-branding-guidelines)
22. [Print & Physical Collateral](#22-print--physical-collateral)
23. [Brand Evolution Roadmap](#23-brand-evolution-roadmap)
24. [Appendix A: Full Agent Inventory](#24-appendix-a-full-agent-inventory)
25. [Appendix B: Implementation Checklists](#25-appendix-b-implementation-checklists)

---

## 1. Executive Summary

### 1.1 Platform Evolution

The HCLS AI Factory has undergone a remarkable transformation since its inception. What began
as a secondary genomics acceleration pipeline — a proof-of-concept for GPU-powered variant
calling on NVIDIA hardware — has matured into a complete precision medicine operating system.
Today, the platform encompasses three core engines, eight intelligence agents spanning every
major medical domain, and an end-to-end workflow that compresses what traditionally requires
months of fragmented analysis into a unified pipeline completing in under five hours.

This evolution was organic, driven by clinical need rather than marketing strategy. Each new
agent was built to solve a real problem: the cardiologist who needed GDMT optimization
integrated with pharmacogenomics, the rare disease family waiting years for diagnosis, the
oncologist needing molecular tumor board intelligence at 2 AM. The result is a platform of
extraordinary depth — 3.5 million annotated variant vectors, 14 knowledge collections per
agent, validated clinical scales, and drug discovery capabilities — but one whose branding
has not kept pace with its capabilities.

### 1.2 The Branding Imperative

As the HCLS AI Factory prepares for its global debut at GTC Europe 2026, the gap between
platform capability and brand clarity has become the single largest barrier to adoption. The
naming is inconsistent. The website navigation was designed for a smaller project. The visual
identity, while professional, does not communicate the platform's three-stage architecture or
the relationships between engines and agents. A pharma executive visiting hcls-ai-factory.org
today would struggle to understand that the Clinical Trial Intelligence Agent is connected to
the same molecular foundation as the Genomic Analysis pipeline.

### 1.3 Scope of This Document

This research guide addresses the branding gap comprehensively. It provides:

- **Unified naming recommendations** for all three engines, eight agents, and the platform
  itself, with rationale grounded in persona research and competitive analysis
- **Website architecture redesign** that reduces navigation from 10 tabs to 7, groups content
  by user intent rather than document type, and creates clear entry points for five distinct
  personas
- **Visual identity specifications** including color assignments per engine, typography standards,
  and icon language that scales from infographics to GitHub badges
- **Infographic specifications** for all eight agents, each following a standard 8K template
  while maintaining unique visual identity through hero elements and color accents
- **Implementation roadmap** prioritized around the GTC Europe launch date, with P0 items
  (naming updates, hero section, 3-engine visualization) completing before the event

### 1.4 Key Recommendations Summary

1. Adopt three-engine naming: **Genomic Foundation Engine**, **Precision Intelligence Engine**,
   **Therapeutic Discovery Engine**
2. Restructure website navigation from 10 tabs to 7: Home | Platform | Engines | Agents |
   Learn | Deploy | Community
3. Add persona-based routing ("Choose Your Path") to the landing page
4. Create 11 agent infographics following a standardized 8K template with agent-specific
   hero visuals
5. Synchronize naming across all touchpoints: website, GitHub, Docker Compose, landing page,
   API documentation, and conference materials

### 1.5 Audience for This Document

This guide is intended for:
- **The core development team** implementing branding changes across the codebase
- **Design collaborators** creating infographics, presentations, and marketing materials
- **Content contributors** writing documentation, blog posts, and conference abstracts
- **Community managers** representing the project on GitHub, social media, and at events
- **NVIDIA partner liaisons** coordinating co-branding for GTC and DGX Spark materials

### 1.6 Success Criteria

The branding implementation will be considered successful when:
- A first-time visitor to hcls-ai-factory.org can identify the platform's three engines and
  their relationship within 30 seconds
- Each of the five target personas can find their entry point within two clicks
- All documentation, code comments, Docker service names, and marketing materials use
  consistent engine and agent naming
- The 11 agent infographics are print-ready (8K) and web-optimized (1080p)
- Lighthouse accessibility scores remain at 90+ across all pages after redesign

---

## 2. The Branding Challenge

### 2.1 Organic Growth, Fragmented Identity

The HCLS AI Factory did not begin with a branding strategy. It began with a clinical need:
accelerate genomic variant calling using NVIDIA GPUs. The first pipeline — Parabricks,
BWA-MEM2, DeepVariant — was named simply "Genomic Analysis." It was accurate, descriptive,
and sufficient for a single-pipeline project.

Then the project grew. A RAG pipeline was added to interpret variants using Claude AI and
Milvus vector search. This was called "Clinical Intelligence" — again accurate, but now the
naming carried an implicit hierarchy. "Analysis" sounded like preprocessing; "Intelligence"
sounded like the real work. The genomic pipeline was foundational, but the name did not
communicate that.

When drug discovery capabilities were added — BioNeMo MolMIM for molecular generation,
DiffDock for binding prediction, RDKit for ADMET profiling — the naming became "Drug
Discovery Orchestrator Agent." This introduced a new problem: the word "Agent" was now
applied to both a processing pipeline (Drug Discovery Orchestrator) and reasoning systems
(the 8 intelligence agents). A first-time reader could not distinguish between pipeline
functions and autonomous reasoning capabilities.

### 2.2 Current Naming Inventory

The following names are currently in use across the codebase, documentation, and website:

| Context | Stage 1 | Stage 2 | Stage 3 |
|---------|---------|---------|---------|
| Code (docker-compose) | core/engines/genomic-foundation | core/engines/precision-intelligence | core/engines/therapeutic-discovery |
| Website (pipelines/) | GPU Genomics | Evidence RAG | Drug Discovery |
| Landing page | Genomic Analysis | Clinical Intelligence | Drug Discovery |
| README | Genomic Analysis | Clinical Intelligence | Drug Discovery Orchestrator Agent |
| Internal docs | Stage 1 | Stage 2 | Stage 3 |
| Demo script | "the genomics engine" | "the intelligence layer" | "the drug discovery agent" |

Six different naming conventions across six touchpoints. No single name is used consistently.
The Docker service names use hyphens and "pipeline"; the website uses short marketing names;
the README uses a fourth variation; and the demo script uses informal descriptions that do
not match any of the above.

### 2.3 The Agent Ambiguity Problem

The platform currently has 8 intelligence agents. These are autonomous reasoning systems
that use RAG (Retrieval-Augmented Generation) to query domain-specific knowledge collections
and provide clinical decision support. They are fundamentally different from the three
processing pipelines — they reason, they do not just process.

However, the Drug Discovery pipeline is also called an "Agent" in some contexts. This creates
confusion:

- "The Drug Discovery Orchestrator Agent generates novel molecules" — this is pipeline processing
- "The Clinical Trial Intelligence Agent optimizes protocol design" — this is autonomous reasoning
- Both are called "Agent" but they operate differently

The branding must draw a clear line between **Engines** (processing pipelines that transform
data) and **Agents** (reasoning systems that interpret data and make recommendations).

### 2.4 Scaling Challenges

The current naming does not scale. With 8 agents today, the platform is already complex.
The roadmap includes:

- Reproductive Genomics Intelligence Agent (prenatal screening, carrier testing)
- Infectious Disease Intelligence Agent (AMR profiling, outbreak genomics)
- Pediatric Precision Agent (growth trajectories, developmental pharmacogenomics)
- Surgical Intelligence Agent (intraoperative genomics, tissue typing)

At 15+ agents, the flat navigation model breaks down. The branding must support hierarchical
organization — grouping agents by clinical domain — without requiring a rename every time a
new agent is added.

### 2.5 The Open-Source Brand Constraint

Unlike commercial platforms with sales teams, implementation consultants, and account managers
who can explain the product verbally, the HCLS AI Factory must explain itself. The branding
must work in contexts where no human is present to interpret it:

- A researcher finding the GitHub repo through a search
- A clinician arriving at hcls-ai-factory.org from a conference QR code
- A pharma executive reading an arXiv paper citation
- A patient advocate following a social media link
- A developer encountering the Docker image on NVIDIA NGC

In each case, the naming, navigation, and visual identity must communicate the platform's
scope, architecture, and relevance — without a sales pitch.

---

## 3. Unified Naming Framework

### 3.1 The Three Engines

The foundational naming decision is what to call the three processing stages. After analyzing
competitive platforms, interviewing clinical stakeholders, and testing names against the five
target personas, we recommend:

#### Genomic Foundation Engine (was "Genomic Analysis")

**Rationale:**
- "Foundation" communicates that this stage creates the molecular layer upon which everything
  else builds. The 3.5 million variant vectors in the `genomic_evidence` collection are
  literally the foundation for all 8 intelligence agents.
- "Engine" distinguishes processing pipelines from reasoning agents. An engine transforms
  input to output; an agent reasons about that output.
- "Foundation" resonates with the AI community, where "foundation model" implies a powerful
  enabling capability that supports many downstream applications. The genomic pipeline is the
  foundation model equivalent for molecular data.
- The word "Analysis" was too generic — every pipeline does "analysis." "Foundation" is
  specific to the architectural role this stage plays.

**Technical scope:**
FASTQ input -> BWA-MEM2 alignment -> DeepVariant variant calling -> VCF output ->
Annotation (ClinVar, AlphaMissense, gnomAD) -> Milvus vector embedding -> 3.5M searchable
variant vectors in `genomic_evidence` collection

**One-line description:**
"Transforms patient DNA into 3.5 million searchable molecular vectors — the foundation for
precision medicine intelligence."

#### Precision Intelligence Engine (was "Clinical Intelligence")

**Rationale:**
- "Network" captures the interconnected nature of the 8 agents. They are not isolated tools;
  they share a common molecular foundation (`genomic_evidence`), cross-reference each other's
  domain knowledge, and collectively provide whole-patient intelligence.
- "Precision" connects directly to the precision medicine mission and differentiates from
  generic "clinical" intelligence (which could mean billing optimization or scheduling).
- The name scales naturally: "The Precision Intelligence Engine includes 8 agents" today
  becomes "15 agents" tomorrow without renaming.
- Each agent can be described as "part of the Precision Intelligence Engine," giving every
  agent immediate context and belonging.
- A cardiologist understands that their Cardiology Intelligence Agent is part of a network;
  an oncologist understands the same about their Precision Oncology Agent. The network
  metaphor communicates cross-specialty coordination.

**Technical scope:**
Milvus vector search across domain-specific knowledge collections -> Claude AI reasoning with
RAG -> Guideline-directed clinical decision support -> Cross-agent coordination ->
Structured outputs (reports, recommendations, risk scores)

**One-line description:**
"Eleven specialized intelligence agents, connected through shared molecular evidence,
delivering precision clinical reasoning across every major medical domain."

#### Therapeutic Discovery Engine (was "Drug Discovery Orchestrator Agent")

**Rationale:**
- "Therapeutic" is broader than "Drug" — it encompasses small molecules, biologics, gene
  therapies, CAR-T constructs, and antisense oligonucleotides. As precision medicine expands
  beyond traditional pharmaceuticals, the name must not constrain.
- "Discovery" is the established scientific term (drug discovery, lead discovery, target
  discovery). It communicates the generative nature of this stage.
- "Engine" maintains naming consistency with the Genomic Foundation Engine and reinforces the
  distinction between engines (processing) and agents (reasoning).
- Removing "Agent" from this stage's name eliminates the agent-pipeline ambiguity described
  in Section 2.3.
- Removing "Orchestrator" simplifies the name. The orchestration function (Nextflow DSL2) is
  an implementation detail, not a brand-level concept.

**Technical scope:**
Validated therapeutic targets (from Precision Intelligence Engine) -> BioNeMo MolMIM
molecular generation -> DiffDock binding pose prediction -> RDKit ADMET profiling ->
Ranked candidate list with binding affinity, drug-likeness, and toxicity scores

**One-line description:**
"Transforms validated therapeutic targets into ranked drug candidates with predicted binding
affinity, toxicity profiles, and drug-likeness scores."

### 3.2 The Naming Hierarchy

The complete platform naming hierarchy:

```
HCLS AI Factory (Platform)
│
├── Genomic Foundation Engine (Stage 1 — Processing)
│   ├── BWA-MEM2 Alignment
│   ├── DeepVariant Variant Calling
│   ├── ClinVar / AlphaMissense / gnomAD Annotation
│   └── Milvus Vector Embedding (genomic_evidence collection)
│
├── Precision Intelligence Engine (Stage 2 — Reasoning)
│   ├── Precision Oncology Agent
│   ├── Precision Biomarker Agent
│   ├── Precision Autoimmune Agent
│   ├── CAR-T Intelligence Agent
│   ├── Imaging Intelligence Agent
│   ├── Pharmacogenomics Intelligence Agent
│   ├── Cardiology Intelligence Agent
│   ├── Neurology Intelligence Agent
│   ├── Rare Disease Diagnostic Agent
│   ├── Clinical Trial Intelligence Agent
│   ├── Single-Cell Intelligence Agent
│   └── Shared: genomic_evidence collection (3.5M vectors)
│
└── Therapeutic Discovery Engine (Stage 3 — Generation)
    ├── BioNeMo MolMIM (Molecular Generation)
    ├── DiffDock (Binding Prediction)
    ├── RDKit (ADMET Profiling)
    └── Ranked Candidate Output
```

### 3.3 Agent Naming Convention Analysis

The 8 agents currently use three naming patterns:

| Pattern | Agents | Count |
|---------|--------|-------|
| "Precision [X] Agent" | Biomarker, Oncology, Autoimmune | 3 |
| "[X] Intelligence Agent" | CAR-T, Imaging, Pharmacogenomics, Cardiology, Clinical Trial, Neurology, Single-Cell | 7 |
| "[X] Diagnostic Agent" | Rare Disease | 1 |

**Full Agent Naming Table:**

| # | Current Name | Recommended Name | Change? | Rationale |
|---|-------------|------------------|---------|-----------|
| 1 | Precision Biomarker Agent | Precision Biomarker Agent | No | "Precision" communicates clinical specificity; biomarker work requires exactitude |
| 2 | Precision Oncology Agent | Precision Oncology Agent | No | "Precision oncology" is an established medical term (NCI definition) |
| 3 | Precision Autoimmune Agent | Precision Autoimmune Agent | No | Pattern consistency with other "Precision" agents; autoimmune dx requires precision |
| 4 | CAR-T Intelligence Agent | CAR-T Intelligence Agent | No | "Intelligence" captures the autonomous reasoning about CAR-T constructs and toxicity |
| 5 | Imaging Intelligence Agent | Imaging Intelligence Agent | No | "Intelligence" distinguishes from image processing tools (VISTA-3D does processing; the agent reasons) |
| 6 | Pharmacogenomics Intelligence Agent | Pharmacogenomics Intelligence Agent | No | "Intelligence" captures the complex star allele interpretation and phenoconversion logic |
| 7 | Cardiology Intelligence Agent | Cardiology Intelligence Agent | No | "Intelligence" captures GDMT optimization and multi-calculator reasoning |
| 8 | Neurology Intelligence Agent | Neurology Intelligence Agent | No | "Intelligence" captures ATN staging and multi-scale clinical reasoning |
| 9 | Rare Disease Diagnostic Agent | Rare Disease Diagnostic Agent | No | "Diagnostic" is uniquely appropriate — the primary mission IS diagnosis (ending the diagnostic odyssey) |
| 10 | Clinical Trial Intelligence Agent | Clinical Trial Intelligence Agent | No | "Intelligence" captures protocol optimization and adaptive design reasoning |
| 11 | Single-Cell Intelligence Agent | Single-Cell Intelligence Agent | No | "Intelligence" captures TME classification and spatial transcriptomics reasoning |

**Recommendation:** Keep all current agent names unchanged. The mixed naming pattern is
actually a feature, not a bug:

- **"Precision"** implies clinical specificity and exactitude — appropriate for agents that
  work with precise molecular and biomarker data (oncology, biomarkers, autoimmune)
- **"Intelligence"** implies autonomous reasoning capability — appropriate for agents that
  perform complex multi-step analysis (cardiology GDMT, pharmacogenomics phenoconversion,
  neurology ATN staging)
- **"Diagnostic"** implies problem-solving with a definitive answer — uniquely appropriate
  for the Rare Disease agent, whose entire mission is to reach a diagnosis

The variation communicates that each agent has its own character and specialty, while the
shared "Agent" suffix maintains family membership.

### 3.4 Tagline Options

Different contexts require different taglines. The following options are organized by use case:

**Hero tagline (landing page, main README):**
> "From Patient DNA to Therapeutic Discovery in Under 5 Hours"

This is the primary tagline. It communicates the end-to-end nature of the platform, the
speed advantage, and the patient-centric mission in a single sentence.

**Architecture tagline (technical docs, architecture pages):**
> "Foundation. Intelligence. Discovery. One Platform."

This communicates the three-engine architecture in four words.

**Open-source tagline (GitHub, community pages):**
> "Precision Medicine for Everyone — Open Source on NVIDIA DGX Spark"

This communicates the equity mission and the accessibility of the $4,699 hardware target.

**Conference tagline (GTC Europe, NVIDIA partner materials):**
> "8 Intelligence Agents. 8 Engines. One Mission."

This communicates scale and unity — impressive numbers with singular purpose.

**Clinical tagline (clinician-facing materials):**
> "Every Specialty. Every Patient. One Molecular Foundation."

This communicates cross-specialty coverage built on shared genomic data.

**Recommendation:** Adopt all five taglines for their respective contexts. The primary tagline
("From Patient DNA to Therapeutic Discovery in Under 5 Hours") should appear on the landing
page hero and GitHub README. Others should be used in their designated contexts as specified.

---

## 4. Persona-Based Messaging

### 4.1 The Clinician

**Profile:** Cardiologist, oncologist, neurologist, or other specialist seeking clinical
decision support tools. Typically 35-60 years old, limited time, skeptical of AI claims,
deeply values clinical accuracy and guideline concordance. Most likely entry point: their
specialty agent page or a colleague's recommendation.

**Entry point:** Their specialty agent (Cardiology Intelligence, Neurology Intelligence, etc.)

**Primary message:**
"The [Specialty] Intelligence Agent provides guideline-directed clinical decision support
powered by RAG across [X] specialized knowledge collections. It incorporates [specific
clinical scales], [specific guidelines], and [specific drug databases]. It is part of the
Precision Intelligence Engine — connected to 10 companion agents for cross-specialty
coordination when your patient's needs extend beyond a single domain."

**Key proof points:**
- Validated clinical scales (NIHSS, MELD, CHA2DS2-VASc, etc.)
- Guideline references (ACC/AHA, NCCN, ACR, AAN)
- Response time under 5 seconds
- Transparent source attribution (every recommendation cites its evidence)
- Cross-agent referral capability ("This patient's cardiac symptoms may have a pharmacogenomic
  component — the Pharmacogenomics Intelligence Agent identified CYP2D6 poor metabolizer status")

**What NOT to say:**
- Do not lead with GPU acceleration or infrastructure details
- Do not emphasize the number of agents (clinicians care about THEIR agent)
- Do not use "AI replaces" language — use "AI augments" and "clinical decision support"
- Do not make diagnostic claims — the agent provides decision support, not diagnoses

**Tone:** Authoritative, evidence-based, peer-to-peer. Write as if addressing a colleague
at a medical conference, not a customer.

### 4.2 The Researcher / Bioinformatician

**Profile:** PhD or postdoc in genomics, computational biology, or bioinformatics. 25-45
years old. Evaluates tools by benchmarks, data formats, and reproducibility. Most likely
entry point: GitHub repository, Google Scholar citation, or a colleague's recommendation.

**Entry point:** The Genomic Foundation Engine or Single-Cell Intelligence Agent

**Primary message:**
"GPU-accelerated genomics pipeline processes a 30x whole genome in under 2 hours on NVIDIA
DGX Spark. BWA-MEM2 alignment, DeepVariant variant calling, multi-database annotation
(ClinVar 4.1M records, AlphaMissense 71M predictions, gnomAD population frequencies).
Output: 3.5 million variant vectors searchable across 8 intelligence agents via Milvus.
RAPIDS-accelerated single-cell analysis at 74x CPU speedup. Full pipeline reproducible via
Docker Compose. Apache 2.0."

**Key proof points:**
- Benchmark numbers (alignment time, variant calling sensitivity/specificity, vector query latency)
- Data format specifications (FASTQ in, VCF out, JSONL for vectors)
- API documentation (FastAPI endpoints, request/response schemas)
- Reproducibility (Docker Compose, pinned versions, seed values)
- Comparison to standard tools (GATK, CellRanger, Scanpy)

**What NOT to say:**
- Do not oversimplify the pipeline — researchers will verify every claim
- Do not hide limitations (GPU memory requirements, supported genome builds)
- Do not use marketing superlatives ("revolutionary," "groundbreaking")
- Do not assume they will use the web interface — many will use APIs directly

**Tone:** Technical, precise, humble. Show benchmarks, cite methods, acknowledge limitations.
Write as if submitting to a bioinformatics journal.

### 4.3 The Pharma Executive

**Profile:** VP of Clinical Development, Head of Translational Medicine, or Chief Medical
Officer at a pharmaceutical company. 45-65 years old. Evaluates tools by impact on timelines,
cost, and competitive advantage. Most likely entry point: Clinical Trial Intelligence Agent
or Therapeutic Discovery Engine.

**Entry message:**
"The HCLS AI Factory compresses clinical development timelines by integrating molecular
profiling, protocol optimization, patient matching, and site intelligence into a single
platform. The Clinical Trial Intelligence Agent has analyzed 40 landmark trials across
35 conditions, providing evidence-based protocol optimization including adaptive design
recommendations, endpoint selection, and enrichment strategy analysis. Connected to the
Therapeutic Discovery Engine for target-to-candidate molecular generation."

**Key proof points:**
- Timeline compression estimates (18-24 months in aggregate across development phases)
- Protocol optimization capabilities (40 landmark trials analyzed)
- Adaptive design intelligence (9 adaptive designs characterized)
- Patient matching with molecular profiling integration
- Regulatory document generation (IND sections, DSUR components)
- Cost: $4,699 hardware + open-source software (vs. $500K+ commercial platforms)

**What NOT to say:**
- Do not lead with open-source ideology — lead with business value
- Do not emphasize the technology stack (Claude, Milvus) — emphasize outcomes
- Do not promise specific FDA/EMA regulatory acceptance
- Do not compare negatively to competitors by name

**Tone:** Executive, outcome-oriented, strategic. Write as if presenting to a board of
directors. Lead with impact, support with evidence.

### 4.4 The Engineer / Developer

**Profile:** Software engineer, ML engineer, or DevOps specialist. 25-40 years old.
Evaluates tools by architecture quality, deployment simplicity, code quality, and community
health. Most likely entry point: GitHub repository or Docker Hub.

**Entry message:**
"docker compose up. 8 agents, 3 engines, on a $4,699 machine. Apache 2.0. FastAPI +
Streamlit + Milvus + Claude. 150+ API endpoints. 23-module shared library. Prometheus +
Grafana monitoring. Health checks with auto-restart. Nextflow DSL2 orchestration."

**Key proof points:**
- Architecture diagram (clean, accurate, showing service relationships)
- API documentation (OpenAPI/Swagger for every service)
- Docker Compose configuration (single-file deployment)
- Test coverage (current state and roadmap)
- Contributing guidelines (how to add a new agent)
- Code quality metrics (linting, type hints, docstrings)
- Service port map (11 services, well-organized)

**What NOT to say:**
- Do not hide technical debt — engineers respect honesty about known issues
- Do not oversell test coverage if it is incomplete
- Do not use clinical jargon without explanation
- Do not gate documentation behind registration

**Tone:** Direct, honest, code-first. Show architecture, not adjectives. Write as if
reviewing a pull request.

### 4.5 The Patient / Advocate

**Profile:** Patient with a complex diagnosis, parent of a child with a rare disease, or
patient advocacy organization representative. Any age. Evaluates tools by accessibility,
equity impact, and real patient stories. Most likely entry point: hcls-ai-factory.org
landing page or a news article.

**Entry message:**
"Precision medicine for every patient, everywhere. The HCLS AI Factory is an open-source
platform that helps clinicians make better diagnostic and treatment decisions using your
genomic data and the latest medical research. It runs on affordable hardware, is free to use,
and is designed to bring precision medicine capabilities to hospitals and clinics regardless
of size or budget. Because your diagnosis should not depend on your zip code."

**Key proof points:**
- Open source and free (Apache 2.0 license)
- Runs on affordable hardware ($4,699 NVIDIA DGX Spark)
- Covers rare diseases (88 conditions, 45 genes, 12 gene therapies)
- Reduces diagnostic odyssey (from years to hours for data processing)
- Cross-specialty coordination (8 agents, no referral delays for AI analysis)

**What NOT to say:**
- Do not use technical jargon (VCF, FASTQ, Milvus, RAG)
- Do not make specific diagnostic or treatment promises
- Do not imply this replaces a doctor
- Do not be paternalistic — patients are sophisticated advocates

**Tone:** Warm, hopeful, respectful. Write as if speaking to a family at a rare disease
conference. Be honest about capabilities and limitations.

---

## 5. Website Architecture — hcls-ai-factory.org

### 5.1 Current State Analysis

The current website is built with MkDocs Material theme and deployed on Netlify. Key metrics:

- **Landing page:** Custom `home.html` override (33KB), DNA helix animation, hero section at
  92vh height, gradient background, 6 responsive breakpoints
- **Navigation:** 10 top-level tabs (Home, Architecture, Getting Started, Design, Learning,
  Pipelines, Agents, Reports, Publications, Resources)
- **Content:** 20+ content directories, 150+ markdown files, 2,172 lines of custom CSS
- **Visual quality:** Professional dark theme, NVIDIA green (#76B900) accent, subtle animations,
  consistent card layouts
- **Performance:** Fast static site, Netlify CDN, immutable asset caching
- **Security:** HSTS headers, X-Content-Type-Options, X-Frame-Options

The website is well-built technically. The challenge is not quality — it is organization.
The 10-tab navigation was designed incrementally as content was added, resulting in a
structure organized by document type (Reports, Publications, Resources) rather than by
user intent (What is this? How does it work? How do I deploy it?).

### 5.2 Recommended Navigation Restructure

**Current navigation (10 tabs):**
```
Home | Architecture | Getting Started | Design | Learning | Pipelines | Agents | Reports | Publications | Resources
```

**Recommended navigation (7 tabs):**
```
Home | Platform | Engines | Agents | Learn | Deploy | Community
```

**Rationale:** Reduce cognitive load from 10 choices to 7. Group content by user intent, not
document type. Align tab names with the unified naming framework. Create clear paths for each
of the five personas.

#### Tab 1: Home

Keep the current landing page structure — it is strong. Update content to reflect unified
branding:

- **Hero section:** Update tagline to "From Patient DNA to Therapeutic Discovery in Under
  5 Hours." Add animated three-engine flow below tagline: [DNA helix icon] -> [Network
  node icon] -> [Molecule icon] with engine names appearing on hover/tap.
- **Three Engines section (NEW):** Add immediately below hero. Three side-by-side cards
  with engine-specific colors (see Section 7.1). Each card: engine name, one-line
  description, key metric, "Learn more" link to Engines tab.
- **Choose Your Path section (NEW):** Five persona cards with icons. Clinician -> "Find
  your specialty agent." Researcher -> "Explore the genomic pipeline." Pharma -> "See
  clinical trial intelligence." Developer -> "View the architecture." Patient/Advocate ->
  "Learn how it helps."
- **Keep existing sections:** Metrics bar (8 agents, 3.5M variants, <5 hours), pipeline
  comparison (CPU vs GPU), origin story, why open source, getting started, tech stack.

#### Tab 2: Platform

Merge Architecture + Design + Reports into a single "Platform" tab. This is where users go
to understand the big picture.

- Architecture overview (updated with new engine naming)
- Project Bible (comprehensive reference)
- White Paper (formal publication)
- Performance benchmarks (GPU vs CPU, latency, throughput)
- arXiv Paper (academic citation)
- Intelligence Report (platform capabilities summary)
- System requirements and supported configurations

#### Tab 3: Engines

New tab replacing "Pipelines." Each engine gets its own landing page:

- **Genomic Foundation Engine**
  - Pipeline overview (FASTQ -> VCF -> Vectors)
  - Component descriptions (BWA-MEM2, DeepVariant, annotation databases)
  - Performance benchmarks
  - Input/output specifications
  - API documentation

- **Precision Intelligence Engine**
  - Network overview (8 agents, shared foundation)
  - How agents share the genomic_evidence collection
  - Cross-agent coordination capabilities
  - RAG architecture (Milvus + Claude)
  - Link to each agent's dedicated page

- **Therapeutic Discovery Engine**
  - Pipeline overview (Targets -> Candidates)
  - Component descriptions (BioNeMo MolMIM, DiffDock, RDKit)
  - Example outputs (ranked candidates with scores)
  - Integration with Precision Intelligence Engine
  - API documentation

#### Tab 4: Agents

Keep the Agents tab but reorganize from a flat list to domain-grouped sections:

**Oncology & Therapeutics:**
- Precision Oncology Agent — Molecular tumor board intelligence, evidence tiering, therapy ranking
- CAR-T Intelligence Agent — Construct design, manufacturing protocols, toxicity management
- Single-Cell Intelligence Agent — TME classification, spatial transcriptomics, immune profiling
- Clinical Trial Intelligence Agent — Protocol optimization, adaptive design, patient matching

**Cardiology & Internal Medicine:**
- Cardiology Intelligence Agent — GDMT optimization, risk calculators, genomic triggers
- Pharmacogenomics Intelligence Agent — Star alleles, CPIC guidelines, phenoconversion
- Precision Biomarker Agent — Biological age, disease trajectories, genotype-adjusted ranges

**Neurology & Rare Disease:**
- Neurology Intelligence Agent — ATN staging, validated scales, stroke triage
- Rare Disease Diagnostic Agent — Phenotype matching, ACMG criteria, gene therapy eligibility

**Imaging & Autoimmune:**
- Imaging Intelligence Agent — VISTA-3D segmentation, multimodal fusion, federated learning
- Precision Autoimmune Agent — Autoantibody patterns, HLA associations, flare prediction

#### Tab 5: Learn

Merge Learning + Publications:

- Unified Foundations Guide (the comprehensive 3,349-line learning guide covering all agents)
- Individual agent learning guides (deep-dive into each agent's domain)
- GTC materials (presentations, video scripts, demo recordings)
- Research publications and preprints
- Glossary of terms

#### Tab 6: Deploy

Merge Getting Started + deployment content + Resources:

- Quick-Start Checklist (pre-flight for DGX Spark deployment)
- Docker Compose Deployment Guide (step-by-step)
- DGX Spark Hardware Setup
- Data Setup (downloading reference genomes, ClinVar, AlphaMissense)
- Environment Configuration (.env template, API keys, port assignments)
- Troubleshooting Guide
- Downloads (data packages, model weights, reference files)

#### Tab 7: Community

New tab for open-source community engagement:

- Contributing Guidelines (how to add a new agent, how to submit a fix)
- Code of Conduct
- License (Apache 2.0 full text)
- GitHub Repository Link
- Changelog (versioned release notes)
- Roadmap (planned agents, features, integrations)
- Acknowledgments (NVIDIA, Anthropic, open-source dependencies)

### 5.3 Landing Page Recommendations

The current landing page (`overrides/home.html`, 33KB) is a strong foundation. The hero
section with DNA helix animation is distinctive and should be retained as a brand element.
Recommended updates are additive, not replacement.

**Hero Section Updates:**

Current subtitle: Platform-dependent, varies by context.
Recommended subtitle: "From Genomic Foundation to Precision Intelligence to Therapeutic
Discovery"

Add an animated three-engine flow visualization below the subtitle:
```
[DNA Helix]        [Network Graph]      [Molecule]
     |                   |                   |
Genomic Foundation  Precision Intelligence  Therapeutic Discovery
     Engine              Network             Engine
     |                   |                   |
  Stage 1 ──────────> Stage 2 ──────────> Stage 3
```

Animation: Icons appear sequentially (left to right, 300ms delay each), then the connecting
arrows animate. On mobile, stack vertically.

**Three Engines Section (below hero):**

Three cards, side by side on desktop, stacked on mobile:

Card 1 — Genomic Foundation Engine:
- Border color: NVIDIA Green (#76B900)
- Icon: DNA double helix
- Headline: "Genomic Foundation Engine"
- Subtext: "Patient DNA -> 3.5M annotated variant vectors in 2 hours"
- Key metric: "3.5M vectors | 11.7M variants | ClinVar + AlphaMissense + gnomAD"
- CTA button: "Explore the Foundation" -> /engines/genomic-foundation/

Card 2 — Precision Intelligence Engine:
- Border color: HCLS Cyan (#00B4D8)
- Icon: Connected nodes (network graph)
- Headline: "Precision Intelligence Engine"
- Subtext: "8 specialized agents across every major medical domain"
- Key metric: "8 agents | 14 collections each | <5 second response time"
- CTA button: "Meet the Agents" -> /agents/

Card 3 — Therapeutic Discovery Engine:
- Border color: Discovery Purple (#8b5cf6)
- Icon: Molecular structure (hexagonal)
- Headline: "Therapeutic Discovery Engine"
- Subtext: "Validated targets -> novel drug candidates in minutes"
- Key metric: "MolMIM generation | DiffDock binding | RDKit ADMET profiling"
- CTA button: "See Discovery" -> /engines/therapeutic-discovery/

**Choose Your Path Section:**

Five persona cards with distinct icons and entry points:

1. Clinician (stethoscope icon) -> "Find your specialty agent" -> /agents/
2. Researcher (microscope icon) -> "Explore the genomic pipeline" -> /engines/genomic-foundation/
3. Pharma Executive (chart icon) -> "See clinical trial intelligence" -> /agents/clinical-trial/
4. Developer (terminal icon) -> "View the architecture" -> /platform/architecture/
5. Patient/Advocate (heart icon) -> "Learn how it helps" -> /learn/patient-guide/

### 5.4 Agent Landing Pages

Each agent's `index.md` should follow a standardized template to ensure consistency across
all 8 agents while allowing domain-specific content:

**Template structure:**

```markdown
# [Agent Name]
> One-line description of the agent's clinical purpose

**Part of the [Precision Intelligence Engine](../engines/precision-intelligence-network.md)**

## Overview
2-3 paragraphs describing the agent's purpose, target user, and key capabilities.

## Key Capabilities
- Capability 1 with brief description
- Capability 2 with brief description
- Capability 3 with brief description
- Capability 4 with brief description
- Capability 5 with brief description

## Architecture
Mini-diagram showing the agent's data flow:
Input sources -> Knowledge collections -> RAG pipeline -> Output types

## Knowledge Collections
Table of all collections with record counts and descriptions.

## Quick Links
| Resource | Link |
|----------|------|
| Demo Guide | [View](demo-guide.md) |
| API Documentation | [View](api-docs.md) |
| Deployment Guide | [View](deployment.md) |
| Learning Guide | [View](learning-guide.md) |

## By the Numbers
| Metric | Value |
|--------|-------|
| Knowledge Collections | [X] |
| Total Records | [X] |
| Clinical Scales | [X] |
| Supported Conditions | [X] |

## Cross-Agent Connections
Description of how this agent connects to other agents in the
Precision Intelligence Engine.

## Infographic
![Agent Infographic](infographic.png)
```

---

## 6. GitHub Repository Recommendations

### 6.1 README.md Updates

The primary GitHub README should be updated to reflect unified branding. Structure:

```markdown
# HCLS AI Factory

> From Patient DNA to Therapeutic Discovery in Under 5 Hours

[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)]
[![Agents](https://img.shields.io/badge/agents-11-green.svg)]
[![Python](https://img.shields.io/badge/python-3.10+-blue.svg)]
[![Platform](https://img.shields.io/badge/platform-DGX%20Spark-76B900.svg)]

The HCLS AI Factory is an open-source precision medicine platform that transforms patient
DNA into therapeutic drug candidates in under 5 hours. Built for NVIDIA DGX Spark ($4,699).

## Three Engines

[Architecture diagram: Foundation -> Intelligence -> Discovery]

### Genomic Foundation Engine
Patient DNA -> 3.5M annotated variant vectors...

### Precision Intelligence Engine
8 specialized intelligence agents...

### Therapeutic Discovery Engine
Validated targets -> novel drug candidates...

## The 11 Intelligence Agents
[Table with all agents, one-line descriptions, and links]

## Quick Start
git clone ...
cp .env.template .env
docker compose up

## Documentation
Full documentation at [hcls-ai-factory.org](https://hcls-ai-factory.org)

## License
Apache 2.0
```

### 6.2 Repository Structure Recommendations

```
hcls-ai-factory/
├── README.md                          # Updated with unified branding
├── LICENSE                            # Apache 2.0
├── docker-compose.dgx-spark.yml      # Full stack deployment
├── .env.template                      # Environment template (no secrets)
├── docs/                              # MkDocs source (deployed to Netlify)
│   ├── index.md                       # Landing page content
│   ├── engines/                       # NEW: Engine documentation
│   │   ├── genomic-foundation.md      # Genomic Foundation Engine
│   │   ├── precision-intelligence-network.md  # Precision Intelligence Engine
│   │   └── therapeutic-discovery.md   # Therapeutic Discovery Engine
│   ├── agents/                        # Agent documentation (11 directories)
│   │   ├── precision-oncology/
│   │   ├── precision-biomarker/
│   │   ├── precision-autoimmune/
│   │   ├── cart-intelligence/
│   │   ├── imaging-intelligence/
│   │   ├── pharmacogenomics-intelligence/
│   │   ├── cardiology-intelligence/
│   │   ├── neurology-intelligence/
│   │   ├── rare-disease-diagnostic/
│   │   ├── clinical-trial-intelligence/
│   │   └── single-cell-intelligence/
│   ├── platform/                      # Architecture, benchmarks, white paper
│   ├── learn/                         # Learning guides, publications
│   ├── deploy/                        # Deployment guides, troubleshooting
│   └── community/                     # Contributing, changelog, roadmap
├── mkdocs.yml                         # Updated navigation structure
├── netlify.toml                       # Deployment configuration
├── core/engines/genomic-foundation/                 # Genomic Foundation Engine source
├── core/engines/precision-intelligence/                 # Precision Intelligence Engine source
├── core/engines/therapeutic-discovery/           # Therapeutic Discovery Engine source
├── core/                     # Individual agent source code
├── lib/hcls_common/                   # Shared library (23 modules)
├── hcls-orchestrator/                  # Nextflow + Python orchestration
├── landing-page/                      # Flask landing page (:8080)
├── monitoring/                        # Prometheus + Grafana
└── scripts/                           # Utility scripts
```

### 6.3 Docker Service Naming Alignment

Current Docker service names in `docker-compose.dgx-spark.yml` should be documented with
their unified branding equivalents. While renaming Docker services is a breaking change that
should be deferred to a major version, the mapping should be clear:

| Docker Service | Unified Brand Name | Engine/Network |
|----------------|-------------------|----------------|
| core/engines/genomic-foundation | Genomic Foundation Engine | Engine |
| core/engines/precision-intelligence | Precision Intelligence Engine (core) | Network |
| core/engines/therapeutic-discovery | Therapeutic Discovery Engine | Engine |
| biomarker-agent | Precision Biomarker Agent | Network member |
| oncology-agent | Precision Oncology Agent | Network member |
| autoimmune-agent | Precision Autoimmune Agent | Network member |
| cart-agent | CAR-T Intelligence Agent | Network member |
| imaging-agent | Imaging Intelligence Agent | Network member |

Labels and comments in the Docker Compose file should reference unified names even if
service names remain unchanged for backward compatibility.

---

## 7. Visual Identity System

### 7.1 Color Assignment

The visual identity system assigns specific colors to each engine, creating instant visual
recognition across all touchpoints:

| Element | Color Name | Hex Code | RGB | Usage |
|---------|-----------|----------|-----|-------|
| NVIDIA Brand / Genomic Foundation | NVIDIA Green | #76B900 | 118, 185, 0 | Primary accent, CTA buttons, Genomic Foundation Engine cards and borders |
| Precision Intelligence Engine | HCLS Cyan | #00B4D8 | 0, 180, 216 | Intelligence Network cards, agent page accents, secondary buttons |
| Therapeutic Discovery Engine | Discovery Purple | #8b5cf6 | 139, 92, 246 | Discovery Engine cards, molecule visualizations, tertiary accent |
| Background Primary | Navy | #1B2333 | 27, 35, 51 | Page backgrounds, card backgrounds, header bars |
| Background Deep | Deep Navy | #0a0a0f | 10, 10, 15 | Infographic backgrounds, hero section gradient end |
| Text Primary | White | #FFFFFF | 255, 255, 255 | Body text on dark backgrounds, card text |
| Text Secondary | Light Gray | #94a3b8 | 148, 163, 184 | Secondary text, descriptions, metadata |
| Status: Critical | Alert Red | #e74c3c | 231, 76, 60 | Critical findings, error states, pathogenic variants |
| Status: Warning | Caution Orange | #f39c12 | 243, 156, 18 | Moderate findings, warnings, VUS variants |
| Status: Success | Health Green | #2ecc71 | 46, 204, 113 | Positive status, healthy metrics, benign variants |
| Status: Info | Info Blue | #3498db | 52, 152, 219 | Informational badges, links, neutral data |

**Color Usage Rules:**
1. Engine colors should never be mixed — each engine maintains its distinct color identity
2. On the landing page, the three-engine cards use their respective border colors
3. Agent pages inherit the HCLS Cyan from the Precision Intelligence Engine
4. The NVIDIA Green remains the primary accent for platform-level elements (logo, main CTA)
5. Status colors are consistent across all agents and engines

### 7.2 Engine-Specific Visual Language

Each engine has a distinct visual vocabulary used in diagrams, infographics, and page headers:

**Genomic Foundation Engine — Green Visual Language:**
- Primary icon: DNA double helix (stylized, geometric)
- Background pattern: Circuit board traces (representing data processing)
- Gradient: #76B900 -> #5a8f00 (green to dark green)
- Data visualization style: Linear flow (left to right), showing transformation stages
- Metaphor: Assembly line — raw material (DNA) transformed step by step

**Precision Intelligence Engine — Cyan Visual Language:**
- Primary icon: Network graph with 11 nodes (representing agents)
- Background pattern: Connected dots with lines (representing neural/agent connections)
- Gradient: #00B4D8 -> #0077b6 (cyan to deep blue)
- Data visualization style: Radial/hub-spoke, showing agents around shared foundation
- Metaphor: Neural network — interconnected specialists sharing knowledge

**Therapeutic Discovery Engine — Purple Visual Language:**
- Primary icon: Molecular structure (hexagonal benzene ring with branches)
- Background pattern: Hexagonal grid (representing molecular structures)
- Gradient: #8b5cf6 -> #6d28d9 (purple to deep purple)
- Data visualization style: Funnel/pipeline (many candidates narrowed to ranked list)
- Metaphor: Chemistry lab — computational experiments generating novel compounds

### 7.3 Typography

Typography is standardized across all platforms:

**Primary Font Family: Inter**
- Headlines (H1): Inter Bold, 36-48pt, letter-spacing: -0.02em
- Subheadlines (H2): Inter SemiBold, 24-30pt, letter-spacing: -0.01em
- Section headers (H3): Inter Medium, 18-22pt
- Body text: Inter Regular, 14-16pt, line-height: 1.6
- Captions and labels: Inter Medium, 12-13pt, letter-spacing: 0.02em
- Navigation: Inter Medium, 14pt

**Monospace Font: JetBrains Mono**
- Code blocks: JetBrains Mono Regular, 13-14pt, line-height: 1.5
- Inline code: JetBrains Mono Regular, 13pt
- Terminal output: JetBrains Mono Light, 13pt
- API endpoints: JetBrains Mono Medium, 14pt

**Font Loading Strategy:**
- Load Inter and JetBrains Mono from Google Fonts CDN
- Use `font-display: swap` to prevent FOIT (Flash of Invisible Text)
- Fallback stack: Inter, -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif
- Monospace fallback: JetBrains Mono, "Fira Code", "Cascadia Code", monospace

### 7.4 Iconography

Standard icons for platform elements:

| Element | Icon | Source | Usage |
|---------|------|--------|-------|
| Genomic Foundation Engine | DNA helix | Custom SVG | Engine cards, navigation, diagrams |
| Precision Intelligence Engine | Network graph | Custom SVG | Engine cards, navigation, diagrams |
| Therapeutic Discovery Engine | Molecule | Custom SVG | Engine cards, navigation, diagrams |
| Agent (generic) | Brain with circuit | Material Icons | Agent list items, badges |
| API Endpoint | Terminal cursor | Material Icons | API documentation |
| Collection | Database cylinder | Material Icons | Knowledge collection references |
| Clinical Scale | Gauge meter | Material Icons | Scale references in agent pages |
| Performance | Lightning bolt | Material Icons | Benchmark sections |
| Open Source | Unlock/Open lock | Material Icons | License badges, community pages |
| DGX Spark | NVIDIA logo | NVIDIA brand kit | Hardware references |

---

## 8. Infographic Specifications — All 11 Agents

### 8.1 Standard Infographic Template

All 11 agent infographics follow a standardized template to ensure visual consistency while
allowing domain-specific content.

**Resolution:** 7680 x 4320 pixels (8K UHD) for print; downscaled to 1920 x 1080 (1080p)
for web display. Export formats: PNG (web), PDF (print), SVG (scalable).

**Layout Grid (8K canvas):**

```
+------------------------------------------------------------------------+
|  [HCLS AI Factory Logo]              [NVIDIA DGX Spark Badge]    [v1.0] |
|  Top Banner — 400px height, navy gradient                              |
+------------------------------------------------------------------------+
|                          |                                             |
|  LEFT COLUMN (40%)       |  RIGHT COLUMN (60%)                        |
|  3072px width            |  4608px width                              |
|                          |                                             |
|  Agent Architecture      |  Hero Visual Element                       |
|  Diagram                 |  (agent-specific, see 8.2)                 |
|                          |                                             |
|  Data flow showing:      |  Key Capabilities                          |
|  Input -> Collections    |  (3-4 callout cards with icons)            |
|  -> RAG -> Output        |                                             |
|                          |  By the Numbers                            |
|  Port assignment         |  (metric cards: collections, records,      |
|  API endpoints           |   conditions, drugs, genes)                |
|                          |                                             |
|  Cross-Agent             |  Cross-Modal Flow                          |
|  Connections             |  (4-stage pipeline specific to agent)      |
|  (which agents this      |                                             |
|   connects to)           |                                             |
|                          |                                             |
+------------------------------------------------------------------------+
|  Bottom Bar — 300px height                                             |
|  [Tech Stack Icons] [Port: XXXX] [Apache 2.0] [hcls-ai-factory.org]   |
+------------------------------------------------------------------------+
```

**Color scheme:** Dark navy background (#0a0a0f) with subtle grid pattern (1px lines at
#1a1a2e, 200px spacing). Agent-specific accent color used for borders, highlights, and
callout cards. All agents use HCLS Cyan (#00B4D8) as their primary accent since they are
all part of the Precision Intelligence Engine.

**Typography on infographics:**
- Agent name: Inter Bold, 120pt (8K) / 30pt (1080p)
- Section headers: Inter SemiBold, 72pt (8K) / 18pt (1080p)
- Body text: Inter Regular, 48pt (8K) / 12pt (1080p)
- Metric numbers: Inter Bold, 96pt (8K) / 24pt (1080p)
- Labels: Inter Medium, 36pt (8K) / 9pt (1080p)

### 8.2 Agent-Specific Infographic Specifications

#### 8.2.1 Imaging Intelligence Agent

**Hero visual:** Brain MRI cross-section with VISTA-3D AI segmentation overlay. Left half
shows raw MRI; right half shows AI-segmented regions with color-coded tissue classes.
Subtle animation lines suggest real-time processing.

**Accent color:** #00B4D8 (HCLS Cyan) with secondary #06d6a0 (medical imaging green)

**Key capability callouts (4 cards):**
1. **VISTA-3D Segmentation** — 132-class automated organ and lesion segmentation using
   NVIDIA NIM microservices. Sub-millimeter precision on CT, MRI, and PET modalities.
2. **Multimodal Fusion** — Correlate imaging findings with genomic variants and clinical
   biomarkers. Radiogenomics pipeline links imaging phenotypes to molecular drivers.
3. **Federated Learning** — Privacy-preserving model training across institutions without
   sharing raw imaging data. NVIDIA FLARE integration for multi-site collaboration.
4. **Clinical Workflows** — 4 end-to-end workflows: diagnostic interpretation, treatment
   response monitoring, surgical planning, and screening program optimization.

**By the Numbers:**
| Metric | Value |
|--------|-------|
| Knowledge Collections | 10+1 |
| Research Papers Indexed | 2,678 |
| NIM Services | 4 |
| Segmentation Classes | 132 |
| Supported Modalities | CT, MRI, PET, US |

**Cross-modal flow:**
Imaging Acquisition -> AI Segmentation -> Genomic Correlation -> Treatment Response

#### 8.2.2 Precision Oncology Agent

**Hero visual:** Molecular tumor board visualization — circular layout with patient genome
at center, surrounded by concentric rings: variant tier (innermost), therapy options,
clinical trials, resistance mechanisms (outermost). Evidence level indicated by ring opacity.

**Accent color:** #00B4D8 (HCLS Cyan) with secondary #e74c3c (oncology red for tumor markers)

**Key capability callouts (4 cards):**
1. **Evidence Tiering** — AMP/ASCO/CAP 4-tier evidence classification for somatic variants.
   Automated assessment of clinical significance from Level 1 (FDA-approved companion
   diagnostic) through Level 4 (biological relevance only).
2. **Therapy Ranking** — Multi-criteria therapy scoring incorporating evidence level, tumor
   type specificity, drug availability, patient genomics (pharmacogenomics cross-referral),
   and real-world evidence from indexed trials.
3. **Trial Matching** — Automated matching of patient molecular profile to active clinical
   trials. Integration with Clinical Trial Intelligence Agent for protocol feasibility
   assessment and site selection.
4. **Resistance Tracking** — Longitudinal monitoring of resistance mechanisms. Known resistance
   mutations mapped against current therapy, with proactive alternative therapy suggestions
   when resistance probability exceeds threshold.

**By the Numbers:**
| Metric | Value |
|--------|-------|
| Knowledge Collections | 11 |
| Gene-Drug Pairs | 100+ |
| Evidence Levels | 4 (AMP/ASCO/CAP) |
| Tumor Types Covered | 30+ |
| Resistance Mechanisms | 200+ |

**Cross-modal flow:**
VCF Variants -> Evidence Tiering -> Therapy Ranking -> Trial Matching

#### 8.2.3 Precision Biomarker Agent

**Hero visual:** Biological age clock — circular gauge showing chronological age on outer
ring and PhenoAge biological age on inner ring, with the 9 PhenoAge biomarkers displayed as
radial spokes. Color gradient from green (favorable) to red (unfavorable) indicates each
biomarker's contribution to biological age acceleration or deceleration.

**Accent color:** #00B4D8 (HCLS Cyan) with secondary #f39c12 (biomarker gold)

**Key capability callouts (4 cards):**
1. **PhenoAge Calculation** — 9-biomarker biological age assessment (albumin, creatinine,
   glucose, CRP, lymphocyte %, MCV, RDW, alkaline phosphatase, WBC). Validated against
   mortality outcomes in NHANES III cohort.
2. **Disease Trajectories** — Longitudinal biomarker tracking across 6 disease categories
   with trajectory prediction. Early detection of biomarker patterns preceding clinical
   disease by 3-7 years.
3. **Genotype-Adjusted Ranges** — Reference ranges adjusted for pharmacogenomic status.
   CYP2D6 poor metabolizers have different "normal" drug levels than extensive metabolizers.
   Cross-referral to Pharmacogenomics Intelligence Agent for star allele context.
4. **PGx Profiling** — 7 pharmacogene assessment integrated with biomarker interpretation.
   Drug level biomarkers interpreted in the context of metabolizer phenotype.

**By the Numbers:**
| Metric | Value |
|--------|-------|
| Knowledge Collections | 10+1 |
| Disease Categories | 6 |
| PhenoAge Biomarkers | 9 |
| Pharmacogenes Profiled | 7 |
| Demo Patients | 5 |

**Cross-modal flow:**
Lab Results -> Biological Age -> Disease Trajectory -> Prevention Strategy

#### 8.2.4 CAR-T Intelligence Agent

**Hero visual:** CAR protein 3D structure rendered as a stylized molecular diagram showing
the modular architecture: scFv antigen-binding domain -> flexible hinge region ->
transmembrane domain -> costimulatory domain (CD28 or 4-1BB) -> CD3-zeta signaling domain.
Each module color-coded and labeled. Six FDA-approved products shown as variants branching
from the base design.

**Accent color:** #00B4D8 (HCLS Cyan) with secondary #8b5cf6 (immunotherapy purple)

**Key capability callouts (4 cards):**
1. **6 FDA Products** — Comprehensive analysis of all FDA-approved CAR-T products
   (Kymriah, Yescarta, Tecartus, Breyanzi, Abecma, Carvykti). Head-to-head comparison of
   construct design, manufacturing, efficacy, and toxicity profiles.
2. **Comparative Analysis** — Multi-dimensional comparison across CAR-T products: binding
   domain (scFv source), costimulatory domain (CD28 vs 4-1BB), hinge/TM design, clinical
   trial outcomes, real-world evidence, cost, and manufacturing turnaround time.
3. **Manufacturing Intelligence** — End-to-end manufacturing workflow analysis: leukapheresis
   -> T-cell selection -> activation -> transduction -> expansion -> quality control ->
   cryopreservation -> infusion. Failure mode analysis at each step.
4. **Toxicity Management** — CRS (Cytokine Release Storm) and ICANS (Immune effector
   Cell-Associated Neurotoxicity Syndrome) grading, monitoring, and management protocols.
   Tocilizumab and corticosteroid dosing algorithms.

**By the Numbers:**
| Metric | Value |
|--------|-------|
| Knowledge Collections | 11 |
| Research Papers Indexed | 5,047 |
| Therapeutic Targets | 25 |
| FDA-Approved Products | 6 |
| Manufacturing Steps | 7 |

**Cross-modal flow:**
Target Selection -> Construct Design -> Manufacturing -> Infusion & Monitoring

#### 8.2.5 Precision Autoimmune Agent

**Hero visual:** Diagnostic odyssey timeline — horizontal timeline showing the typical
4.5-year autoimmune diagnostic journey (symptom onset -> multiple specialist visits ->
misdiagnoses -> correct diagnosis) contrasted with the AI-accelerated path (symptom input ->
autoantibody pattern recognition -> HLA association -> diagnosis in minutes). The timelines
converge at "Diagnosis" with a dramatic compression visual.

**Accent color:** #00B4D8 (HCLS Cyan) with secondary #e74c3c (autoimmune/inflammatory red)

**Key capability callouts (4 cards):**
1. **Autoantibody Patterns** — Pattern recognition across 50+ autoantibodies. Multi-antibody
   profiles mapped to disease-specific patterns (e.g., anti-dsDNA + anti-Smith = SLE;
   anti-CCP + RF = RA). Bayesian probability scoring for differential diagnosis.
2. **HLA Associations** — 50+ HLA allele-disease associations with population-adjusted risk
   calculations. Includes protective alleles (HLA-DRB1*04:02 protects against type 1
   diabetes in some populations) and risk alleles (HLA-B*27 and ankylosing spondylitis).
3. **Disease Activity Scoring** — Validated scoring systems for 13 autoimmune diseases:
   SLEDAI (lupus), DAS28 (RA), BVAS (vasculitis), CDAI (Crohn's), Mayo score (UC), and
   more. Longitudinal tracking with flare detection.
4. **Flare Prediction** — Machine learning-informed flare risk assessment based on biomarker
   trends, medication adherence, environmental triggers, and historical flare patterns.
   Early warning system for treatment intensification.

**By the Numbers:**
| Metric | Value |
|--------|-------|
| Knowledge Collections | 14 |
| HLA Alleles Profiled | 50+ |
| Autoimmune Diseases | 13 |
| Demo Patients | 5 |
| Autoantibodies Tracked | 50+ |

**Cross-modal flow:**
Symptoms -> Autoantibody Patterns -> HLA & Genomics -> Diagnosis & Therapy

#### 8.2.6 Pharmacogenomics Intelligence Agent

**Hero visual:** CYP enzyme metabolism visualization — stylized liver cell with major CYP
enzymes (CYP2D6, CYP2C19, CYP2C9, CYP3A4, CYP1A2) shown as molecular machines, each
processing drugs into metabolites. Drug molecules enter from the left, pass through the
enzyme, and exit as metabolites on the right. Star allele badges on each enzyme indicate
metabolizer status: UM (ultrarapid), EM (extensive/normal), IM (intermediate), PM (poor).

**Accent color:** #00B4D8 (HCLS Cyan) with secondary #2ecc71 (pharmacology green)

**Key capability callouts (4 cards):**
1. **7 Pharmacogenes** — Comprehensive profiling of CYP2D6, CYP2C19, CYP2C9, CYP3A4,
   CYP1A2, VKORC1, and DPYD. Star allele calling from genomic data with diplotype-to-
   phenotype translation following PharmVar nomenclature.
2. **CPIC Guidelines** — Clinical Pharmacogenetics Implementation Consortium guideline
   integration for 100+ drug-gene pairs. Automated dosing recommendations based on
   metabolizer phenotype (e.g., "CYP2D6 PM: avoid codeine, use morphine alternative").
3. **Phenoconversion** — Phenoconversion analysis accounting for concomitant medications that
   inhibit or induce CYP enzymes. A CYP2D6 normal metabolizer on fluoxetine (strong 2D6
   inhibitor) is functionally a poor metabolizer. Real-time adjusted phenotype calculation.
4. **HLA Pharmacovigilance** — 15 HLA-drug associations for serious adverse drug reactions:
   HLA-B*57:01 (abacavir hypersensitivity), HLA-B*15:02 (carbamazepine SJS/TEN), HLA-A*31:01
   (carbamazepine DRESS). Pre-prescription screening recommendations.

**By the Numbers:**
| Metric | Value |
|--------|-------|
| Knowledge Collections | 15 |
| Drug-Gene Pairs | 100+ |
| Dosing Algorithms | 9 |
| HLA-Drug Associations | 15 |
| Pharmacogenes Profiled | 7 |

**Cross-modal flow:**
Genotype -> Star Alleles -> Metabolizer Phenotype -> Drug Selection & Dosing

#### 8.2.7 Cardiology Intelligence Agent

**Hero visual:** Anatomical heart with 4-pillar GDMT (Guideline-Directed Medical Therapy)
visualization. Four pillars rising from the heart: ARNI/ACEi/ARB (Pillar 1), Beta-blocker
(Pillar 2), MRA (Pillar 3), SGLT2i (Pillar 4). Each pillar shows current dose level as a
fill gauge, with target dose marked. Color coding: green (at target), yellow (up-titrating),
red (not yet started). Genomic triggers shown as lightning bolts connecting to relevant
pillars.

**Accent color:** #00B4D8 (HCLS Cyan) with secondary #e74c3c (cardiology red)

**Key capability callouts (4 cards):**
1. **6 Risk Calculators** — Integrated cardiovascular risk assessment: HEART Score,
   ASCVD Pooled Cohort Equations, CHA2DS2-VASc (atrial fibrillation stroke risk), HAS-BLED
   (bleeding risk), MAGGIC (heart failure prognosis), and EuroSCORE II (surgical risk).
2. **11 Clinical Workflows** — End-to-end decision support for heart failure management,
   atrial fibrillation, acute coronary syndrome, hypertension, valvular disease, pulmonary
   hypertension, cardiomyopathy, lipid management, cardiac device optimization, cardiac
   rehabilitation, and preventive cardiology.
3. **GDMT Optimizer** — Automated 4-pillar GDMT optimization for HFrEF: tracks current
   medications, calculates percentage of target dose, recommends up-titration schedule,
   monitors contraindications, and adjusts for renal function and hemodynamic parameters.
4. **18 Genomic Triggers** — Pharmacogenomic and genetic triggers affecting cardiac care:
   CYP2C19 (clopidogrel), CYP2D6 (metoprolol), TTN (dilated cardiomyopathy), MYBPC3 (HCM),
   SCN5A (Brugada), KCNQ1 (Long QT), and 12 additional cardiac-relevant genes.

**By the Numbers:**
| Metric | Value |
|--------|-------|
| Knowledge Collections | 13 |
| Cardiac Conditions | 45 |
| Cardiac Genes | 56 |
| GDMT Drugs | 14 |
| Clinical Workflows | 11 |

**Cross-modal flow:**
Imaging/ECG -> Risk Calculators -> GDMT Optimization -> Genomic Personalization

#### 8.2.8 Neurology Intelligence Agent

**Hero visual:** Brain visualization with 10 clinical scale gauges arranged around the
perimeter. Each gauge shows a validated neurological scale (NIHSS, MoCA, UPDRS, EDSS,
Glasgow, etc.) as a semicircular meter with the current score indicated. The brain itself
shows functional regions highlighted according to the patient's neurological profile.
ATN biomarker staging (Amyloid-Tau-Neurodegeneration) shown as a three-letter classification
badge.

**Accent color:** #00B4D8 (HCLS Cyan) with secondary #9b59b6 (neurology purple)

**Key capability callouts (4 cards):**
1. **8 Clinical Workflows** — Comprehensive neurological assessment workflows: stroke triage
   (NIHSS + imaging), dementia evaluation (ATN staging), epilepsy management, Parkinson's
   disease monitoring, multiple sclerosis tracking, headache classification, neuropathy
   workup, and neuromuscular disease evaluation.
2. **10 Validated Scales** — Automated scoring and interpretation: NIHSS (stroke severity),
   MoCA (cognitive screening), UPDRS (Parkinson's), EDSS (MS disability), Glasgow Coma
   Scale, Mini-Mental State Examination, CDR (dementia staging), ALSFRS-R (ALS functional
   rating), Epworth Sleepiness Scale, and PHQ-9 (depression screening in neurological
   patients).
3. **ATN Staging** — Alzheimer's disease biological classification: Amyloid status (A+/A-),
   Tau status (T+/T-), Neurodegeneration status (N+/N-). Eight possible combinations with
   clinical interpretation and treatment implications per the 2024 revised criteria.
4. **Stroke Triage** — Acute stroke decision support: NIHSS calculation, large vessel
   occlusion prediction, thrombolysis eligibility (tPA window), thrombectomy candidacy
   (extended window with perfusion imaging), and post-stroke secondary prevention.

**By the Numbers:**
| Metric | Value |
|--------|-------|
| Knowledge Collections | 14 |
| Neurological Conditions | 58 |
| Neurological Drugs | 43 |
| Neurological Genes | 38 |
| Validated Clinical Scales | 10 |

**Cross-modal flow:**
EEG/MRI/CSF -> Clinical Scales -> Diagnosis & Staging -> Treatment & Monitoring

#### 8.2.9 Rare Disease Diagnostic Agent

**Hero visual:** Diagnostic odyssey compression — dramatic visual showing the traditional
5-7 year diagnostic odyssey (meandering path through multiple specialists, dead ends, and
misdiagnoses) compressed into a direct path through AI-powered phenotype-to-gene matching.
The traditional path is shown as a long, winding road in muted gray; the AI-accelerated
path is a straight, illuminated line from symptoms to diagnosis.

**Accent color:** #00B4D8 (HCLS Cyan) with secondary #f1c40f (rare disease gold/zebra)

**Key capability callouts (4 cards):**
1. **HPO Phenotype Matching** — Human Phenotype Ontology-based symptom encoding with
   semantic similarity matching against 88 rare diseases. Hierarchical phenotype propagation
   (a specific phenotype implies its parent terms). Ranked differential diagnosis with
   phenotype match scores.
2. **ACMG 28 Criteria** — Full implementation of ACMG/AMP variant interpretation guidelines:
   8 pathogenic criteria (PVS1, PS1-4, PM1-2), 8 likely pathogenic combinations, 4 benign
   criteria (BA1, BS1-3), 4 likely benign criteria (BP1-4), and 4 VUS criteria. Automated
   criterion assessment with evidence citations.
3. **Gene Therapy Eligibility** — Assessment of eligibility for 12 FDA-approved or
   late-stage gene therapies: Luxturna (RPE65), Zolgensma (SMA), Skysona (ALD), Hemgenix
   (hemophilia B), Roctavian (hemophilia A), Casgevy (SCD/beta-thal), Lyfgenia (SCD),
   Elevidys (DMD), and others. Eligibility criteria, dosing, monitoring requirements.
4. **Family Segregation** — Pedigree analysis and variant segregation assessment. De novo
   variant detection, carrier testing recommendations, recurrence risk calculation, and
   reproductive counseling information. Integration with Genomic Foundation Engine for
   family member variant calling.

**By the Numbers:**
| Metric | Value |
|--------|-------|
| Knowledge Collections | 14 |
| Rare Diseases Covered | 88 |
| Disease-Associated Genes | 45 |
| Gene Therapies | 12 |
| ACMG Criteria Implemented | 28 |

**Cross-modal flow:**
Phenotype Description -> Gene Matching -> Variant Classification -> Therapy Eligibility

#### 8.2.10 Clinical Trial Intelligence Agent

**Hero visual:** Protocol optimization pipeline — funnel visualization showing the trial
design process: broad therapeutic area at top, narrowing through endpoint selection, patient
population definition, adaptive design incorporation, site selection, and enrollment
optimization. Each stage of the funnel shows the number of options considered vs. recommended.
40 landmark trial reference icons line the sides.

**Accent color:** #00B4D8 (HCLS Cyan) with secondary #3498db (clinical trial blue)

**Key capability callouts (4 cards):**
1. **10 Clinical Workflows** — End-to-end trial intelligence: protocol optimization, patient
   matching, site selection, enrollment forecasting, endpoint selection, adaptive design
   recommendation, regulatory strategy, competitor analysis, budget estimation, and
   safety monitoring plan design.
2. **40 Landmark Trials** — Knowledge base of 40 practice-changing clinical trials across
   35 conditions. Each trial analyzed for design strengths, limitations, statistical approach,
   endpoint selection, and lessons learned. Used as reference cases for new protocol design.
3. **5 Decision Engines** — Protocol Optimizer (design recommendations), Patient Matcher
   (eligibility assessment), Site Intelligence (site selection and performance prediction),
   Enrollment Forecaster (timeline and dropout modeling), and Regulatory Advisor (filing
   strategy and document generation).
4. **Adaptive Design** — 9 adaptive design frameworks characterized: group sequential, sample
   size re-estimation, response-adaptive randomization, biomarker-adaptive, dose-finding,
   seamless phase II/III, platform trials, basket trials, and umbrella trials. Design
   selection based on therapeutic area and development phase.

**By the Numbers:**
| Metric | Value |
|--------|-------|
| Knowledge Collections | 14 |
| Conditions Covered | 35 |
| Drugs Referenced | 34 |
| Adaptive Design Types | 9 |
| Landmark Trials Analyzed | 40 |

**Cross-modal flow:**
Therapeutic Target -> Protocol Design -> Patient Matching -> Site Selection & Enrollment

#### 8.2.11 Single-Cell Intelligence Agent

**Hero visual:** UMAP dimensionality reduction plot with tumor microenvironment (TME)
classification overlay. Clusters of cells rendered as point clouds, color-coded by cell type
(44 types). Four quadrant labels show TME classifications: immune-hot (top-left, warm red),
immune-cold (top-right, cool blue), immune-excluded (bottom-left, gray), and
immunosuppressive (bottom-right, dark purple). Overlay arrows show ligand-receptor
interactions between cell clusters.

**Accent color:** #00B4D8 (HCLS Cyan) with secondary #06d6a0 (single-cell/bioinformatics green)

**Key capability callouts (4 cards):**
1. **10 Analytical Workflows** — Comprehensive single-cell analysis: quality control and
   preprocessing, normalization and batch correction, dimensionality reduction (PCA/UMAP/
   t-SNE), clustering, cell type annotation, differential expression, trajectory inference,
   RNA velocity, spatial deconvolution, and cell-cell communication analysis.
2. **TME Classifier** — Tumor microenvironment classification into 4 immune phenotypes:
   immune-hot (inflamed, checkpoint-responsive), immune-cold (immune desert, needs priming),
   immune-excluded (stromal barrier, needs remodeling), and immunosuppressive (Treg/MDSC-
   enriched, needs depletion). Treatment strategy recommendation per phenotype.
3. **GPU Acceleration** — RAPIDS-accelerated single-cell workflows achieving 74x speedup
   over CPU-based tools (Scanpy). cuML for dimensionality reduction, cuGraph for community
   detection, cuDF for data manipulation. Memory-efficient processing of million-cell
   datasets on DGX Spark.
4. **Ligand-Receptor Analysis** — 25 curated ligand-receptor interaction pairs across
   immune, stromal, and tumor compartments. CellPhoneDB-style analysis with GPU
   acceleration. Visualization of communication networks between cell types.

**By the Numbers:**
| Metric | Value |
|--------|-------|
| Knowledge Collections | 12 |
| Cell Types Classified | 44 |
| Drugs Referenced | 30 |
| Surface Markers | 75 |
| Ligand-Receptor Pairs | 25 |

**Cross-modal flow:**
Bulk Sequencing -> Single-Cell Deconvolution -> TME Classification -> Treatment Strategy

---

## 9. Netlify Deployment Optimization

### 9.1 Build Performance

The current deployment builds all 150+ pages on every Netlify deploy. While MkDocs builds
are fast (typically under 60 seconds for this site), optimization is warranted as the site
grows:

- **Pin mkdocs-material version:** Currently the the platform site may use a different Material
  theme version than the public site. Pin to a specific version (e.g., `mkdocs-material==9.5.x`)
  in both `requirements.txt` files to prevent visual inconsistencies.
- **Incremental builds:** MkDocs Material does not natively support incremental builds, but
  Netlify's build cache can be leveraged by caching the `site/` directory between builds.
  Add to `netlify.toml`:
  ```toml
  [build]
    publish = "site"
    command = "pip install -r requirements.txt && mkdocs build"

  [build.environment]
    PYTHON_VERSION = "3.10"
  ```
- **Asset optimization:** Large infographic PNGs should be pre-optimized before commit.
  Target: <500KB per web-resolution infographic. Use `pngquant` or `optipng` in the build
  pipeline.

### 9.2 SEO Optimization

Current SEO is minimal (site-level meta only). Recommended enhancements:

**Structured Data (JSON-LD):**
Add JSON-LD to each agent page for rich search results:
```json
{
  "@context": "https://schema.org",
  "@type": "SoftwareApplication",
  "name": "HCLS AI Factory - [Agent Name]",
  "description": "[Agent one-line description]",
  "applicationCategory": "HealthApplication",
  "operatingSystem": "Linux (NVIDIA DGX Spark)",
  "license": "https://www.apache.org/licenses/LICENSE-2.0",
  "isPartOf": {
    "@type": "SoftwareApplication",
    "name": "HCLS AI Factory"
  }
}
```

**Sitemap generation:** Enable the MkDocs sitemap plugin (built into Material theme):
```yaml
plugins:
  - search
  - sitemap:
      changefreq: weekly
      priority: 0.8
```

**Canonical URLs:** Add canonical URLs to prevent duplicate content issues between
hcls-ai-factory.org (public site) and any the platform-specific deployment:
```yaml
extra:
  canonical_url: "https://hcls-ai-factory.org/"
```

**Per-page meta descriptions:** Add meta descriptions to each page's front matter:
```yaml
---
description: "The Cardiology Intelligence Agent provides GDMT optimization, risk
  calculators, and genomic-informed cardiac care as part of the Precision Intelligence
  Network."
---
```

### 9.3 Performance Headers

Current Netlify headers are solid (HSTS, X-Content-Type-Options, X-Frame-Options).
Recommended additions:

```toml
[[headers]]
  for = "/*"
  [headers.values]
    Content-Security-Policy = "default-src 'self'; script-src 'self' 'unsafe-inline' https://cdn.jsdelivr.net; style-src 'self' 'unsafe-inline' https://fonts.googleapis.com; font-src 'self' https://fonts.gstatic.com; img-src 'self' data: https:;"
    Permissions-Policy = "camera=(), microphone=(), geolocation=(), payment=()"
    Referrer-Policy = "strict-origin-when-cross-origin"

[[headers]]
  for = "/assets/*"
  [headers.values]
    Cache-Control = "public, max-age=31536000, immutable"

[[headers]]
  for = "/infographics/*"
  [headers.values]
    Cache-Control = "public, max-age=604800"
```

---

## 10. Cross-Platform Consistency

### 10.1 GitHub README to Website to Agent Docs

Naming consistency must be enforced across all touchpoints. The following mapping must be
applied universally:

| Old Name | New Name | Contexts to Update |
|----------|----------|-------------------|
| Genomic Analysis | Genomic Foundation Engine | README, landing page, mkdocs nav, demo script, slides |
| Clinical Intelligence | Precision Intelligence Engine | README, landing page, mkdocs nav, demo script, slides |
| Drug Discovery Orchestrator Agent | Therapeutic Discovery Engine | README, landing page, mkdocs nav, demo script, slides |
| GPU Genomics (website) | Genomic Foundation Engine | mkdocs pipelines section |
| Evidence RAG (website) | Precision Intelligence Engine | mkdocs pipelines section |
| Drug Discovery (website) | Therapeutic Discovery Engine | mkdocs pipelines section |

Each agent's README and index.md must include the line:
> "Part of the **Precision Intelligence Engine** — the HCLS AI Factory's reasoning layer
> connecting 8 specialized intelligence agents through shared molecular evidence."

### 10.2 Documentation Sync Strategy

The documentation ecosystem spans three repositories:
1. `hcls-ai-factory/` (internal) — source of truth for all content
2. `hcls-ai-factory-public/` — public GitHub release
3. `hcls-ai-factory-vast/` — the AI platform deployment documentation (MkDocs + Netlify)

**Sync rules:**
- Source of truth: Agent documentation lives in `core/[agent]/docs/`
- Public release: Curated subset copied to `hcls-ai-factory-public/docs/`
- the platform site: Deployment-specific content only; links to public site for agent details
- Naming convention: UPPERCASE filenames locally (e.g., `LEARNING_GUIDE.md`), lowercase
  hyphenated in web-facing repos (e.g., `learning-guide.md`)
- Automated sync: Future GitHub Actions workflow to copy from internal -> public on release
  tags. Manual sync for now with a checklist (see Appendix B).

### 10.3 Versioning Strategy

Branding changes should be versioned alongside software releases:
- **v1.0:** Current state (pre-branding update)
- **v2.0:** Unified branding applied (target: GTC Europe 2026)
- **v2.1+:** Incremental additions (new agents, new infographics)

The version number should appear in the website footer, README badge, and Docker image tags.

---

## 11. Accessibility & Internationalization

### 11.1 Accessibility Audit

The current website scores well on accessibility fundamentals:
- Focus states present on interactive elements
- `prefers-reduced-motion` media query respects user preferences
- High contrast between text and background (navy + white)
- Responsive layout works across screen sizes

**Recommended improvements:**

1. **Skip-to-content link:** Add a hidden-until-focused skip link as the first focusable
   element on every page:
   ```html
   <a href="#content" class="skip-link">Skip to main content</a>
   ```

2. **ARIA landmarks:** Add ARIA roles to custom sections in `home.html`:
   ```html
   <section role="region" aria-label="Three Engines Overview">
   <nav role="navigation" aria-label="Choose Your Path">
   ```

3. **Image alt text:** All infographics must include detailed alt text describing the
   key information conveyed. Alt text for an 8K infographic should be 2-3 sentences
   summarizing the agent's purpose and key numbers.

4. **Color contrast verification:** The new engine colors must be verified against WCAG
   2.1 AA contrast ratios:
   - NVIDIA Green (#76B900) on Navy (#1B2333): ratio 4.8:1 (passes AA for large text)
   - HCLS Cyan (#00B4D8) on Navy (#1B2333): ratio 6.2:1 (passes AA)
   - Discovery Purple (#8b5cf6) on Navy (#1B2333): ratio 4.1:1 (passes AA for large text)
   - All text colors should be verified at both 14pt and 18pt sizes

5. **Screen reader testing:** Test the landing page hero section with NVDA (Windows) and
   VoiceOver (macOS) to ensure the DNA animation does not interfere with content reading
   order. The animated three-engine flow should have a static text alternative.

### 11.2 Internationalization Readiness

The platform currently serves English-only content. For a global launch at GTC Europe 2026,
internationalization readiness is important even if translation does not happen immediately.

**Current i18n status:** No i18n infrastructure in place. All content is English. No
language selection UI.

**MkDocs Material i18n support:** The Material theme includes a built-in i18n plugin that
supports multi-language sites with language switcher:
```yaml
plugins:
  - i18n:
      default_language: en
      languages:
        en: English
        es: Espanol
        pt: Portugues
        fr: Francais
        zh: Chinese (Simplified)
```

**Translation priority (for future implementation):**
1. **Tier 1 (GTC Europe):** English (complete), German (host country)
2. **Tier 2 (Global health equity):** Spanish, Portuguese, French (covers Latin America,
   Africa, Southeast Asia — regions where precision medicine access is most needed)
3. **Tier 3 (Market expansion):** Mandarin Chinese, Japanese, Korean, Hindi

**Pages to translate first:**
- Home (landing page)
- Getting Started / Quick Start
- Agent overview pages (11 pages, ~500 words each)
- Patient/Advocate messaging page

**Content strategy for translation:**
- Keep technical terms in English (FASTQ, VCF, RAG, GPU) — they are universal
- Translate clinical terms according to local medical terminology standards
- Use ISO 639-1 language codes in URLs: `/en/agents/`, `/es/agents/`
- Maintain English as the canonical version for SEO

---

## 12. Product Requirements Document — Branding Implementation

### 12.1 Functional Requirements

| ID | Requirement | Priority | Effort |
|----|------------|----------|--------|
| FR-01 | Update all references from "Genomic Analysis" to "Genomic Foundation Engine" across documentation, website, and marketing materials | P0 | Medium |
| FR-02 | Update all references from "Clinical Intelligence" to "Precision Intelligence Engine" across documentation, website, and marketing materials | P0 | Medium |
| FR-03 | Update all references from "Drug Discovery Orchestrator Agent" to "Therapeutic Discovery Engine" across documentation, website, and marketing materials | P0 | Medium |
| FR-04 | Update landing page hero section with unified tagline: "From Patient DNA to Therapeutic Discovery in Under 5 Hours" | P0 | Low |
| FR-05 | Add three-engine visualization section below hero on landing page with engine-specific colors and descriptions | P0 | High |
| FR-06 | Add persona-based "Choose Your Path" section to landing page with 5 persona cards | P1 | Medium |
| FR-07 | Restructure website navigation from 10 tabs to 7 tabs (Home, Platform, Engines, Agents, Learn, Deploy, Community) | P1 | High |
| FR-08 | Create three engine landing pages: Genomic Foundation Engine, Precision Intelligence Engine, Therapeutic Discovery Engine | P1 | High |
| FR-09 | Group agents by clinical domain in navigation (Oncology & Therapeutics, Cardiology & Internal Medicine, Neurology & Rare Disease, Imaging & Autoimmune) | P1 | Medium |
| FR-10 | Add "Part of the Precision Intelligence Engine" badge to all 11 agent index.md pages | P0 | Low |
| FR-11 | Create standardized infographic for each of the 8 agents following 8K template specification | P0 | Very High |
| FR-12 | Update GitHub README.md with unified branding, three-engine architecture, and all 8 agents | P0 | Medium |
| FR-13 | Add structured data (JSON-LD) to all agent pages for SEO | P2 | Medium |
| FR-14 | Synchronize naming across docker-compose labels, landing page SERVICES dict, demo script, and slide decks | P0 | Medium |

### 12.2 Non-Functional Requirements

| ID | Requirement | Metric | Target |
|----|------------|--------|--------|
| NFR-01 | Landing page load time | Time to First Contentful Paint | < 3 seconds on 4G connection |
| NFR-02 | Accessibility score | Lighthouse Accessibility audit | 90+ on all pages |
| NFR-03 | Agent page consistency | Template compliance | All 8 agents follow identical structure |
| NFR-04 | Infographic resolution | Export resolution | 8K (7680x4320) for print, 1080p for web |
| NFR-05 | Offline capability | Service worker coverage | All pages cached for offline reading (Material PWA plugin) |
| NFR-06 | SEO performance | Google Search Console | All agent pages indexed within 30 days of launch |
| NFR-07 | Cross-browser compatibility | Browser testing | Chrome, Firefox, Safari, Edge (latest 2 versions) |
| NFR-08 | Mobile responsiveness | Mobile usability | No horizontal scroll, readable text at 375px width |

### 12.3 Implementation Priority

**P0 — Before GTC Europe 2026 (Target: 4 weeks before event)**

These items are critical for the global launch and must be completed first:

1. **Engine renaming across all documentation** (FR-01, FR-02, FR-03)
   - Grep and replace across all .md files in docs/, core/, and README files
   - Update mkdocs.yml navigation labels
   - Update landing page HTML
   - Estimated effort: 2 days (1 person)

2. **Landing page hero update** (FR-04)
   - Update subtitle text in home.html
   - Update meta description
   - Estimated effort: 1 hour

3. **Three-engine visualization section** (FR-05)
   - Design and implement three-card section below hero
   - Responsive layout (side-by-side desktop, stacked mobile)
   - Engine-specific border colors
   - Estimated effort: 1 day (frontend)

4. **Agent network badge** (FR-10)
   - Add "Part of the Precision Intelligence Engine" to all 11 agent index.md
   - Estimated effort: 2 hours

5. **Agent infographics** (FR-11)
   - Design template in Figma/Illustrator
   - Create 11 agent-specific infographics
   - Export at 8K and 1080p resolutions
   - Estimated effort: 2 weeks (design resource)

6. **GitHub README update** (FR-12)
   - Rewrite with unified branding and three-engine architecture
   - Estimated effort: 4 hours

7. **Cross-platform naming sync** (FR-14)
   - Audit and update docker-compose labels, SERVICES dict, demo script
   - Estimated effort: 1 day

**P1 — During GTC Europe week (deploy between sessions)**

These items enhance the website but are not blockers for launch:

1. **Navigation restructure** (FR-07)
   - Reorganize mkdocs.yml from 10 tabs to 7
   - Create redirect rules for old URLs
   - Estimated effort: 1 day

2. **Persona routing** (FR-06)
   - Implement "Choose Your Path" section on landing page
   - Estimated effort: 4 hours (frontend)

3. **Engine landing pages** (FR-08)
   - Write content for 3 engine pages
   - Estimated effort: 1 day (content)

4. **Agent domain grouping** (FR-09)
   - Reorganize agent navigation into clinical domain groups
   - Estimated effort: 2 hours

**P2 — Post-GTC Europe (30-day window)**

These items improve long-term discoverability and sustainability:

1. **SEO optimization** (FR-13)
   - Add JSON-LD structured data to all agent pages
   - Enable sitemap plugin
   - Add per-page meta descriptions
   - Estimated effort: 1 day

2. **Internationalization infrastructure**
   - Install and configure i18n plugin
   - Create translation workflow
   - Translate Home page to German (GTC Europe follow-up)
   - Estimated effort: 3 days (setup) + ongoing translation

3. **GitHub Actions sync workflow**
   - Automate documentation sync from internal to public repo
   - Trigger on release tags
   - Estimated effort: 1 day (DevOps)

### 12.4 Acceptance Criteria

The branding implementation is complete when all of the following are true:

- [ ] No instance of "Genomic Analysis" appears in any user-facing documentation
- [ ] No instance of "Clinical Intelligence" (as a stage name) appears in any user-facing documentation
- [ ] No instance of "Drug Discovery Orchestrator Agent" appears in any user-facing documentation
- [ ] The landing page hero displays the unified tagline
- [ ] The three-engine section displays correctly on desktop (1920px) and mobile (375px)
- [ ] All 11 agent pages include the "Precision Intelligence Engine" badge
- [ ] All 11 agent infographics are available at 8K and 1080p resolution
- [ ] The GitHub README displays the three-engine architecture
- [ ] Lighthouse accessibility score is 90+ on landing page and all agent pages
- [ ] Landing page First Contentful Paint is under 3 seconds on 4G simulation

### 12.5 Risk Register

| Risk | Probability | Impact | Mitigation |
|------|------------|--------|------------|
| Engine rename breaks internal references | Medium | High | Run comprehensive grep before and after; maintain old-name redirects |
| Infographic design takes longer than 2 weeks | High | Medium | Start with 3 highest-priority agents (Oncology, Cardiology, Rare Disease); deliver others incrementally |
| Navigation restructure breaks existing bookmarks | Medium | Medium | Implement 301 redirects for all moved pages |
| Persona routing adds decision fatigue instead of reducing it | Low | Medium | A/B test with 5 users before launch; keep "Choose Your Path" section below fold |
| i18n adds maintenance burden | Medium | Low | Defer translation to post-GTC; only set up infrastructure in P2 |

---

*End of Part 1 (Sections 1-12). Part 2 continues below.*

---

# PART 2: PRODUCTION, STRATEGY & IMPLEMENTATION (Sections 13-25)

---

## 13. Infographic Production Specifications

This section defines the complete technical and creative specifications for producing 11 agent infographics plus 3 engine overview infographics for the HCLS AI Factory. These infographics serve as the primary visual marketing asset for GTC Europe, the hcls-ai-factory.org website, social media campaigns, and printed conference materials.

### 13.1 Technical Specifications

**Canvas & Resolution**

| Property | Value |
|---|---|
| Master canvas | 7680 x 4320 pixels (8K UHD / 4320p) |
| Web export | 1920 x 1080 pixels (Full HD) |
| Social export (LinkedIn) | 1200 x 627 pixels |
| Social export (Twitter/X) | 1600 x 900 pixels |
| Print export | 7680 x 4320 @ 300 DPI (25.6" x 14.4") |
| Color depth | 32-bit RGBA |
| Color space | sRGB IEC61966-2.1 (web), Adobe RGB (print) |
| File formats | PNG (web, lossless), SVG (scalable vectors), PDF/X-4 (print) |
| Compression | PNG: maximum compression level 9; WebP fallback at quality 90 |

**Background System**

| Property | Value |
|---|---|
| Primary gradient | Radial, center (#1B2333) to edges (#0a0a0f) |
| Grid overlay | 1px lines, 80px spacing, rgba(255,255,255,0.05) |
| Noise texture | Gaussian monochromatic, 0.5% strength, uniform distribution |
| Agent color glow | Radial, 600px radius, 20% opacity, positioned per agent |

**Grid System**

| Property | Value |
|---|---|
| Columns | 12 |
| Column width | 560px (at 8K) |
| Gutter width | 40px |
| Margins (left/right) | 200px |
| Margins (top/bottom) | 100px |
| Safe zone | 160px inset from all edges (for print bleed) |

**Typography Scale (at 8K master)**

| Element | Font | Weight | Size | Line Height | Color |
|---|---|---|---|---|---|
| Agent name | Inter | Bold (700) | 120pt | 140pt | #FFFFFF |
| Section header | Inter | SemiBold (600) | 72pt | 88pt | #FFFFFF |
| Card title | Inter | Bold (700) | 48pt | 60pt | #FFFFFF |
| Card body | Inter | Regular (400) | 32pt | 44pt | rgba(255,255,255,0.8) |
| Metric number | Inter | Bold (700) | 80pt | 96pt | [engine color] |
| Metric label | Inter | Regular (400) | 28pt | 36pt | rgba(255,255,255,0.7) |
| Badge text | JetBrains Mono | Medium (500) | 24pt | 32pt | #FFFFFF |
| Footer text | Inter | Regular (400) | 24pt | 32pt | rgba(255,255,255,0.5) |

### 13.2 Nanobanana Pro Script Template

The following template defines the standard 7-layer composition for every agent infographic. Each agent section below (13.3) fills in the agent-specific values marked with brackets.

```
----------------------------------------------------------------------
NANOBANANA PRO -- HCLS AI FACTORY AGENT INFOGRAPHIC TEMPLATE
----------------------------------------------------------------------
AGENT:      [Agent Name]
FILENAME:   hcls-[agent-slug]-infographic-8k.png
CANVAS:     7680 x 4320 pixels (8K UHD)
BACKGROUND: radial-gradient(center, #1B2333, #0a0a0f)
GRID:       12-column, 40px gutter, 200px margin
COLOR:      [Engine Color Hex]
----------------------------------------------------------------------

LAYER 1 -- BACKGROUND COMPOSITE
  1a. Fill: radial-gradient(ellipse at 50% 40%, #1B2333 0%, #0a0a0f 100%)
  1b. Grid overlay:
      - Vertical lines: every 80px, stroke rgba(255,255,255,0.05), 1px
      - Horizontal lines: every 80px, stroke rgba(255,255,255,0.05), 1px
  1c. Noise texture: Gaussian, monochromatic, 0.5% strength
  1d. Agent color glow:
      - Shape: radial gradient, 600px radius
      - Color: [engine color] at 20% opacity
      - Position: [glow position -- varies per agent]
  1e. Secondary accent glow (optional):
      - Shape: radial gradient, 400px radius
      - Color: [secondary color] at 10% opacity
      - Position: bottom-right quadrant

LAYER 2 -- HEADER BAR (y: 0-400px)
  2a. Background: linear-gradient(180deg, rgba(0,0,0,0.4) 0%, transparent 100%)
  2b. Left block (x: 200px, y: 150px):
      - HCLS AI Factory logo: logo-header.png
      - Size: 400 x 96px
      - Drop shadow: 0 2px 8px rgba(0,0,0,0.3)
  2c. Center block (x: center, y: 140px):
      - Text: [Agent Name]
      - Font: Inter Bold 120pt, color #FFFFFF
      - Text shadow: 0 2px 12px rgba([engine color], 0.4)
  2d. Right block (x: 7080px right-aligned, y: 120px):
      - Badge 1: "Precision Intelligence Engine" (or appropriate engine)
        - Background: rgba([engine color], 0.15)
        - Border: 1px solid rgba([engine color], 0.4)
        - Text: JetBrains Mono Medium 24pt, [engine color]
        - Padding: 12px 24px
        - Border radius: 8px
      - Badge 2: "NVIDIA DGX Spark" (below Badge 1, 16px gap)
        - Background: rgba(118,185,0,0.15)
        - Border: 1px solid rgba(118,185,0,0.4)
        - Text: JetBrains Mono Medium 24pt, #76B900
  2e. Separator line (y: 390px):
      - Full width (200px to 7480px)
      - Height: 4px
      - Color: [engine color]
      - Glow: 0 0 16px rgba([engine color], 0.6)

LAYER 3 -- HERO VISUAL (y: 440-2360px, x: 200-3160px)
  3a. Bounding box: 2960 x 1920px
  3b. Content: [Agent-specific hero visual -- see 13.3]
  3c. Style directives:
      - Primary rendering in [engine color] tones
      - White and gray for anatomical/structural elements
      - Subtle animation keyframes noted (for future Lottie export)
      - Drop shadow on main element: 0 8px 32px rgba(0,0,0,0.4)
  3d. Caption (bottom of hero area):
      - Text: [Hero caption]
      - Font: Inter Regular 28pt, rgba(255,255,255,0.6)
      - Alignment: center within hero bounding box

LAYER 4 -- CAPABILITIES PANEL (y: 440-2360px, x: 3400-7480px)
  4a. Panel container: 4080 x 1920px
  4b. 4 capability cards, stacked vertically with 40px gap:
      - Card size: 4000 x 440px
      - Card background: rgba(255,255,255,0.03)
      - Card border: 1px solid rgba([engine color], 0.3)
      - Card border-radius: 16px
      - Card padding: 40px
      - Card hover glow (noted for web): 0 0 20px rgba([engine color], 0.1)
  4c. Card internal layout:
      - Icon: 64 x 64px, [engine color], left-aligned
      - Title: Inter Bold 48pt, #FFFFFF, left of icon + 24px gap
      - Description: Inter Regular 32pt, rgba(255,255,255,0.8), below title, 12px gap
      - Right accent: vertical bar 4px wide, [engine color], right edge of card
  4d. Card content:
      Card 1: [Capability 1 -- icon, title, description]
      Card 2: [Capability 2 -- icon, title, description]
      Card 3: [Capability 3 -- icon, title, description]
      Card 4: [Capability 4 -- icon, title, description]

LAYER 5 -- DATA INVENTORY (y: 2440-3160px)
  5a. Section header (y: 2440px):
      - Text: "Data Inventory & Scale"
      - Font: Inter SemiBold 72pt, #FFFFFF
      - Underline: 120px wide, 3px, [engine color]
  5b. Metric circles (y: 2600-3100px):
      - Arranged horizontally, evenly spaced across 7080px
      - Circle diameter: 280px
      - Circle fill: rgba([engine color], 0.08)
      - Circle border: 2px solid rgba([engine color], 0.3)
      - Number: Inter Bold 80pt, [engine color], centered
      - Label: Inter Regular 28pt, rgba(255,255,255,0.7), below circle, 16px gap
      - Subtle pulse animation noted (for Lottie)
  5c. Metrics:
      [Metric 1] | [Metric 2] | [Metric 3] | [Metric 4] | [Metric 5] | [Metric 6]

LAYER 6 -- CROSS-AGENT FLOW (y: 3240-3800px)
  6a. Section header (y: 3240px):
      - Text: "Clinical Pipeline Flow"
      - Font: Inter SemiBold 72pt, #FFFFFF
  6b. Flow diagram (y: 3380-3760px):
      - 5 nodes arranged horizontally
      - Node shape: rounded rectangle, 1200 x 280px
      - Node fill: rgba([engine color], 0.1)
      - Node border: 2px solid [engine color]
      - Node text: Inter Bold 36pt, #FFFFFF (title) + Inter Regular 28pt (subtitle)
      - Arrow connectors: 80px wide, [engine color], animated dash pattern
      - Arrow head: 20px equilateral triangle, filled [engine color]
  6c. Flow steps:
      [Step 1] -> [Step 2] -> [Step 3] -> [Step 4] -> [Step 5]
  6d. Connected agent badges (below flow):
      - Small circular badges (80px) of connected agents
      - Tooltip text noted for web version

LAYER 7 -- FOOTER BAR (y: 3880-4320px)
  7a. Background: rgba(0,0,0,0.3), full width
  7b. Left block (x: 200px):
      - Text: "API: :[port] | UI: :[port]"
      - Font: JetBrains Mono Medium 24pt, rgba(255,255,255,0.5)
  7c. Center block:
      - Tech stack badges arranged horizontally, 16px gap:
        - Each badge: 160 x 48px, rounded 8px
        - Fill: rgba(255,255,255,0.05), border 1px rgba(255,255,255,0.1)
        - Text: JetBrains Mono Medium 20pt, rgba(255,255,255,0.6)
        - Badges: [Milvus] [Claude] [BGE] [FastAPI] [Streamlit] [+ agent-specific]
  7d. Right block (x: 7080px right-aligned):
      - Apache 2.0 badge: green background (#22863a), white text
      - GitHub badge: dark background (#24292e), white text, repo URL
      - 16px gap between badges

----------------------------------------------------------------------
END TEMPLATE
----------------------------------------------------------------------
```

### 13.3 Individual Agent Scripts

---

#### Agent 1: Imaging Intelligence Agent

```
AGENT:      Imaging Intelligence Agent
FILENAME:   hcls-imaging-intelligence-infographic-8k.png
CANVAS:     7680 x 4320 pixels (8K UHD)
BACKGROUND: radial-gradient(center, #1B2333, #0a0a0f)
GRID:       12-column, 40px gutter, 200px margin
COLOR:      #00B4D8 (Cyan -- Precision Intelligence Engine)

LAYER 1 -- BACKGROUND
  Glow position: upper-left quadrant (x: 1600, y: 1200)
  Secondary glow: #00B4D8 at 8% opacity, lower-right (x: 5800, y: 3200)

LAYER 2 -- HEADER BAR
  Agent name: "Imaging Intelligence Agent"
  Engine badge: "Precision Intelligence Engine"
  Separator: 4px solid #00B4D8

LAYER 3 -- HERO VISUAL
  Content: Brain MRI axial cross-section in grayscale with colored
  segmentation overlay. VISTA-3D 132-class output visualized as
  semi-transparent colored regions (hippocampus in cyan, cortical gray
  matter in light blue, white matter in steel blue, ventricles in teal,
  cerebellum in aquamarine). Surrounding the brain: 4 smaller inset
  panels showing (clockwise from top-left): chest X-ray with pneumonia
  heatmap, lung CT with nodule detection bounding boxes, brain
  perfusion map with ischemic core in red/penumbra in yellow, spine MRI
  with MS lesion highlighting.
  Caption: "132-Class AI Segmentation -- From DICOM to Diagnosis in Seconds"

LAYER 4 -- CAPABILITIES
  Card 1:
    Icon: stopwatch
    Title: "4 Clinical Workflows (<90s Stroke)"
    Description: "Acute stroke triage in under 90 seconds. Chest X-ray
    screening in <30s. Lung nodule characterization in <5 minutes. MS
    lesion tracking with longitudinal comparison in <5 minutes."

  Card 2:
    Icon: GPU chip
    Title: "NVIDIA NIM Integration (3 Models)"
    Description: "VISTA-3D for 132-class 3D segmentation. MAISI for
    synthetic medical image generation. VILA-M3 for multimodal visual
    language understanding across radiology modalities."

  Card 3:
    Icon: network nodes
    Title: "Federated Learning (NVIDIA FLARE)"
    Description: "Multi-site model training without data sharing. Privacy-
    preserving collaboration across hospital networks. Continuous model
    improvement from distributed clinical data."

  Card 4:
    Icon: DNA helix + scan
    Title: "Cross-Modal Genomics Integration"
    Description: "Imaging findings trigger targeted gene panels. Lung
    nodules trigger EGFR/ALK/ROS1 panel. Brain tumors trigger MGMT/IDH
    analysis. Cardiac imaging triggers channelopathy screening."

LAYER 5 -- DATA INVENTORY
  Metric 1: "10+1" / "Milvus Collections"
  Metric 2: "2,678" / "Research Papers"
  Metric 3: "4" / "NVIDIA NIMs"
  Metric 4: "4" / "Clinical Workflows"
  Metric 5: "132" / "Segmentation Classes"
  Metric 6: "539" / "Validated Tests"

LAYER 6 -- CROSS-AGENT FLOW
  Step 1: "DICOM Ingest" (icon: upload)
  Step 2: "AI Segmentation" (icon: brain scan)
  Step 3: "Clinical Findings" (icon: report)
  Step 4: "Genomic Triggers" (icon: DNA)
  Step 5: "Treatment Plan" (icon: prescription)
  Connected agents: Precision Oncology, Cardiology, Neurology, Rare Disease

LAYER 7 -- FOOTER
  Ports: "API: :8042 | UI: :8501"
  Tech badges: [Milvus] [Claude] [BGE] [FastAPI] [Streamlit] [VISTA-3D] [MAISI] [VILA-M3]
```

---

#### Agent 2: Precision Oncology Agent

```
AGENT:      Precision Oncology Agent
FILENAME:   hcls-precision-oncology-infographic-8k.png
CANVAS:     7680 x 4320 pixels (8K UHD)
BACKGROUND: radial-gradient(center, #1B2333, #0a0a0f)
GRID:       12-column, 40px gutter, 200px margin
COLOR:      #00B4D8 (Cyan -- Precision Intelligence Engine)

LAYER 1 -- BACKGROUND
  Glow position: center-right (x: 5000, y: 1800)
  Secondary glow: #FF6B6B at 6% opacity, lower-left (x: 1200, y: 3000)

LAYER 2 -- HEADER BAR
  Agent name: "Precision Oncology Agent"
  Engine badge: "Precision Intelligence Engine"
  Separator: 4px solid #00B4D8

LAYER 3 -- HERO VISUAL
  Content: Molecular tumor board visualization. Central table showing
  ranked therapeutic options with evidence tier columns (Tier IA through
  Tier IV). Each row is a therapy with: drug name, target gene, evidence
  level bar (color-coded: green IA, blue IB, yellow II, orange III,
  red IV), clinical trial link badge, and resistance mechanism flag.
  Above the table: a genomic variant call (e.g., "BRAF V600E -- Tier IA")
  with radiating connection lines to matched therapies. Background
  elements: faded tumor cell illustrations with highlighted mutation sites.
  Caption: "AMP/ASCO/CAP Evidence-Tiered Therapy Ranking"

LAYER 4 -- CAPABILITIES
  Card 1:
    Icon: tier pyramid
    Title: "AMP/ASCO/CAP Evidence Tiering"
    Description: "Automated classification of genomic variants into
    standardized evidence tiers (IA, IB, IIC, IID, III, IV) following
    AMP/ASCO/CAP joint consensus guidelines for somatic variant
    interpretation."

  Card 2:
    Icon: ranked list
    Title: "Therapy Ranking Engine"
    Description: "Multi-factor therapy scoring incorporating evidence tier,
    tumor type match, variant allele frequency, co-occurring mutations,
    prior treatment history, and available clinical trial options for
    100+ gene-drug pairs."

  Card 3:
    Icon: shield with X
    Title: "Resistance Mechanism Detection"
    Description: "Identification of known and predicted resistance mutations
    (e.g., EGFR T790M, ALK G1202R). Proactive alternative therapy
    suggestion when resistance markers detected in genomic profile."

  Card 4:
    Icon: clipboard with checkmark
    Title: "Clinical Trial Matching"
    Description: "Automated matching of patient genomic profile to active
    clinical trials. Filters by tumor type, biomarker, phase, geography,
    and eligibility criteria. Direct integration with Clinical Trial Agent."

LAYER 5 -- DATA INVENTORY
  Metric 1: "11" / "Milvus Collections"
  Metric 2: "100+" / "Gene-Drug Pairs"
  Metric 3: "4" / "Evidence Levels"
  Metric 4: "10" / "Cancer Types"
  Metric 5: "25+" / "Resistance Markers"
  Metric 6: "350+" / "Validated Tests"

LAYER 6 -- CROSS-AGENT FLOW
  Step 1: "Genomic Variant" (icon: DNA mutation)
  Step 2: "Evidence Tier" (icon: classification)
  Step 3: "Ranked Therapies" (icon: drug list)
  Step 4: "Trial Match" (icon: trial badge)
  Step 5: "Treatment Decision" (icon: physician)
  Connected agents: Biomarker, CAR-T, Clinical Trial, Pharmacogenomics

LAYER 7 -- FOOTER
  Ports: "API: :8021 | UI: :8502"
  Tech badges: [Milvus] [Claude] [BGE] [FastAPI] [Streamlit] [ClinVar] [OncoKB]
```

---

#### Agent 3: Precision Biomarker Agent

```
AGENT:      Precision Biomarker Agent
FILENAME:   hcls-precision-biomarker-infographic-8k.png
CANVAS:     7680 x 4320 pixels (8K UHD)
BACKGROUND: radial-gradient(center, #1B2333, #0a0a0f)
GRID:       12-column, 40px gutter, 200px margin
COLOR:      #00B4D8 (Cyan -- Precision Intelligence Engine)

LAYER 1 -- BACKGROUND
  Glow position: center (x: 3840, y: 1600)
  Secondary glow: #FFD700 at 6% opacity, upper-right (x: 6000, y: 800)

LAYER 2 -- HEADER BAR
  Agent name: "Precision Biomarker Agent"
  Engine badge: "Precision Intelligence Engine"
  Separator: 4px solid #00B4D8

LAYER 3 -- HERO VISUAL
  Content: Biological age clock visualization. Central circular clock
  face with chronological age on the outer ring and biological age
  (PhenoAge) on the inner ring. 9 biomarker spokes radiating outward
  like a radar/spider chart: albumin, creatinine, glucose, C-reactive
  protein, lymphocyte percent, mean cell volume, red cell distribution
  width, alkaline phosphatase, white blood cell count. Each spoke shows
  the patient value as a dot with color coding (green = normal,
  yellow = borderline, red = abnormal). Center displays: "Chronological
  Age: 55 | Biological Age: 48 | Delta: -7 years".
  Caption: "PhenoAge Biological Clock -- 9 Biomarkers Decode Aging Trajectory"

LAYER 4 -- CAPABILITIES
  Card 1:
    Icon: clock/aging
    Title: "PhenoAge Biological Age (9 Biomarkers)"
    Description: "Levine PhenoAge algorithm computing biological age from
    9 standard blood biomarkers. Quantifies aging acceleration or
    deceleration. Tracks longitudinal changes across visits to measure
    intervention effectiveness."

  Card 2:
    Icon: trajectory arrow
    Title: "6-Category Disease Trajectories"
    Description: "Predictive disease trajectory modeling across cardiovascular,
    metabolic, inflammatory, hepatic, renal, and hematologic categories.
    Each category generates 5-year risk projection with confidence intervals."

  Card 3:
    Icon: DNA + range bar
    Title: "Genotype-Adjusted Reference Ranges"
    Description: "Population-specific and genotype-adjusted reference ranges
    for all biomarkers. Accounts for pharmacogenomic variants that alter
    baseline values. Reduces false positive alerts by 30-40%."

  Card 4:
    Icon: pill + gene
    Title: "Pharmacogenomic Profiling (7 Genes)"
    Description: "Integrated PGx profiling for CYP2D6, CYP2C19, CYP2C9,
    CYP3A5, DPYD, TPMT, and VKORC1. Star allele calling, phenotype
    assignment, and actionable dosing recommendations linked to biomarker
    panel results."

LAYER 5 -- DATA INVENTORY
  Metric 1: "10+1" / "Milvus Collections"
  Metric 2: "9" / "Aging Biomarkers"
  Metric 3: "6" / "Disease Categories"
  Metric 4: "7" / "PGx Genes"
  Metric 5: "5" / "Demo Patients"
  Metric 6: "280+" / "Validated Tests"

LAYER 6 -- CROSS-AGENT FLOW
  Step 1: "Blood Draw" (icon: vial)
  Step 2: "Biomarker Panel" (icon: lab results)
  Step 3: "Biological Age" (icon: clock)
  Step 4: "Disease Risk" (icon: warning)
  Step 5: "Prevention Plan" (icon: shield)
  Connected agents: Pharmacogenomics, Cardiology, Autoimmune

LAYER 7 -- FOOTER
  Ports: "API: :8011 | UI: :8503"
  Tech badges: [Milvus] [Claude] [BGE] [FastAPI] [Streamlit] [PhenoAge]
```

---

#### Agent 4: CAR-T Intelligence Agent

```
AGENT:      CAR-T Intelligence Agent
FILENAME:   hcls-car-t-intelligence-infographic-8k.png
CANVAS:     7680 x 4320 pixels (8K UHD)
BACKGROUND: radial-gradient(center, #1B2333, #0a0a0f)
GRID:       12-column, 40px gutter, 200px margin
COLOR:      #00B4D8 (Cyan -- Precision Intelligence Engine)

LAYER 1 -- BACKGROUND
  Glow position: center-left (x: 2000, y: 1800)
  Secondary glow: #FF4444 at 6% opacity, right (x: 6200, y: 2400)

LAYER 2 -- HEADER BAR
  Agent name: "CAR-T Intelligence Agent"
  Engine badge: "Precision Intelligence Engine"
  Separator: 4px solid #00B4D8

LAYER 3 -- HERO VISUAL
  Content: Chimeric Antigen Receptor protein structure diagram.
  Extracellular domain at top: scFv (single-chain variable fragment)
  with VH and VL labeled, connected by a flexible linker. Hinge region
  (CD8alpha) connecting to the transmembrane domain spanning a stylized
  cell membrane. Intracellular signaling domains stacked below: first
  costimulatory domain (4-1BB or CD28, shown as toggleable comparison)
  and CD3-zeta activation domain at bottom. Each domain labeled with
  function annotations. Side panel: comparison table of 4-1BB vs CD28
  costimulatory domains (persistence vs rapid expansion). Six FDA-
  approved product logos arranged around the structure: Kymriah,
  Yescarta, Tecartus, Breyanzi, Abecma, Carvykti.
  Caption: "CAR Architecture -- From Antigen Recognition to T-Cell Activation"

LAYER 4 -- CAPABILITIES
  Card 1:
    Icon: FDA badge
    Title: "6 FDA Product Analysis"
    Description: "Comprehensive intelligence on all 6 FDA-approved CAR-T
    products: tisagenlecleucel (Kymriah), axicabtagene ciloleucel
    (Yescarta), brexucabtagene autoleucel (Tecartus), lisocabtagene
    maraleucel (Breyanzi), idecabtagene vicleucel (Abecma), and
    ciltacabtagene autoleucel (Carvykti)."

  Card 2:
    Icon: comparison arrows
    Title: "Comparative Mode (4-1BB vs CD28)"
    Description: "Side-by-side analysis of costimulatory domain impact on
    CAR-T persistence, expansion kinetics, exhaustion profiles, and
    clinical outcomes. 4-1BB: slower expansion, longer persistence.
    CD28: rapid expansion, higher initial response rates."

  Card 3:
    Icon: factory/manufacturing
    Title: "Manufacturing Intelligence"
    Description: "End-to-end manufacturing pipeline analysis: leukapheresis
    collection, T-cell activation, viral transduction efficiency,
    expansion monitoring, quality control release criteria, cryopreservation,
    and logistics. Vein-to-vein time optimization."

  Card 4:
    Icon: warning/medical
    Title: "CRS/ICANS Management"
    Description: "Cytokine Release Syndrome and Immune Effector Cell-
    Associated Neurotoxicity Syndrome grading (ASTCT consensus), risk
    prediction models, tocilizumab/corticosteroid dosing algorithms,
    and ICU escalation criteria."

LAYER 5 -- DATA INVENTORY
  Metric 1: "11" / "Milvus Collections"
  Metric 2: "5,047" / "Research Papers"
  Metric 3: "973" / "Clinical Trials"
  Metric 4: "25" / "Targets Tracked"
  Metric 5: "6" / "FDA Products"
  Metric 6: "420+" / "Validated Tests"

LAYER 6 -- CROSS-AGENT FLOW
  Step 1: "Target Selection" (icon: target)
  Step 2: "CAR Design" (icon: construct)
  Step 3: "Manufacturing" (icon: factory)
  Step 4: "Infusion" (icon: IV drip)
  Step 5: "Monitoring" (icon: vital signs)
  Connected agents: Single-Cell, Oncology, Clinical Trial, Pharmacogenomics

LAYER 7 -- FOOTER
  Ports: "API: :8031 | UI: :8504"
  Tech badges: [Milvus] [Claude] [BGE] [FastAPI] [Streamlit] [ClinicalTrials.gov]
```

---

#### Agent 5: Precision Autoimmune Agent

```
AGENT:      Precision Autoimmune Agent
FILENAME:   hcls-precision-autoimmune-infographic-8k.png
CANVAS:     7680 x 4320 pixels (8K UHD)
BACKGROUND: radial-gradient(center, #1B2333, #0a0a0f)
GRID:       12-column, 40px gutter, 200px margin
COLOR:      #00B4D8 (Cyan -- Precision Intelligence Engine)

LAYER 1 -- BACKGROUND
  Glow position: upper-center (x: 3840, y: 1000)
  Secondary glow: #FF8C00 at 5% opacity, lower-left (x: 1400, y: 3400)

LAYER 2 -- HEADER BAR
  Agent name: "Precision Autoimmune Agent"
  Engine badge: "Precision Intelligence Engine"
  Separator: 4px solid #00B4D8

LAYER 3 -- HERO VISUAL
  Content: Diagnostic odyssey timeline visualization. Left side shows
  a compressed traditional path: 5-7 year horizontal timeline with
  icons (initial symptoms -> GP visit -> referral -> misdiagnosis ->
  specialist -> more tests -> re-referral -> diagnosis) with a frustrated
  patient silhouette and "Average: 4.6 years to diagnosis" label.
  Right side shows the AI-accelerated path: single visit timeline
  (symptoms -> autoantibody panel -> HLA typing -> AI pattern match ->
  confirmed diagnosis) compressed into a single horizontal bar labeled
  "Minutes". Large "5-7 YEARS -> MINUTES" callout text connecting the
  two timelines with a dramatic arrow.
  Caption: "Ending the Diagnostic Odyssey -- 13 Autoimmune Diseases, One Platform"

LAYER 4 -- CAPABILITIES
  Card 1:
    Icon: antibody Y-shape
    Title: "Autoantibody Pattern Recognition"
    Description: "Multi-pattern autoantibody analysis across ANA, anti-dsDNA,
    anti-Smith, anti-CCP, RF, ANCA (c-ANCA/p-ANCA), anti-TPO, anti-
    transglutaminase, and disease-specific panels. Pattern clustering
    identifies overlap syndromes and undifferentiated connective tissue disease."

  Card 2:
    Icon: DNA with HLA
    Title: "HLA Association Analysis (50+ Alleles)"
    Description: "Comprehensive HLA typing correlation with 50+ alleles
    across Class I (HLA-A, -B, -C) and Class II (HLA-DR, -DQ, -DP).
    Risk quantification for ankylosing spondylitis (B27), celiac disease
    (DQ2/DQ8), type 1 diabetes (DR3/DR4), and rheumatoid arthritis (DRB1)."

  Card 3:
    Icon: gauge/dashboard
    Title: "Disease Activity Scoring (4 Validated)"
    Description: "Automated calculation and tracking of DAS28 (rheumatoid
    arthritis), SLEDAI (systemic lupus), CDAI (Crohn's disease), and
    BASDAI (ankylosing spondylitis). Longitudinal trend visualization
    with remission/low/moderate/high activity thresholds."

  Card 4:
    Icon: lightning bolt
    Title: "Flare Prediction & Prevention"
    Description: "Machine learning-informed flare risk prediction using
    biomarker trends, activity score trajectories, medication adherence
    data, and seasonal/environmental triggers. Early warning alerts
    enable preemptive treatment adjustments."

LAYER 5 -- DATA INVENTORY
  Metric 1: "14" / "Milvus Collections"
  Metric 2: "13" / "Autoimmune Diseases"
  Metric 3: "50+" / "HLA Alleles"
  Metric 4: "4" / "Activity Scores"
  Metric 5: "5" / "Demo Patients"
  Metric 6: "310+" / "Validated Tests"

LAYER 6 -- CROSS-AGENT FLOW
  Step 1: "Symptoms" (icon: patient)
  Step 2: "Autoantibodies" (icon: antibody)
  Step 3: "HLA Typing" (icon: DNA)
  Step 4: "AI Diagnosis" (icon: brain)
  Step 5: "Biologic Therapy" (icon: infusion)
  Connected agents: Biomarker, Pharmacogenomics, Imaging, Clinical Trial

LAYER 7 -- FOOTER
  Ports: "API: :8051 | UI: :8505"
  Tech badges: [Milvus] [Claude] [BGE] [FastAPI] [Streamlit] [HLA-DB]
```

---

#### Agent 6: Pharmacogenomics Intelligence Agent

```
AGENT:      Pharmacogenomics Intelligence Agent
FILENAME:   hcls-pharmacogenomics-infographic-8k.png
CANVAS:     7680 x 4320 pixels (8K UHD)
BACKGROUND: radial-gradient(center, #1B2333, #0a0a0f)
GRID:       12-column, 40px gutter, 200px margin
COLOR:      #00B4D8 (Cyan -- Precision Intelligence Engine)

LAYER 1 -- BACKGROUND
  Glow position: center (x: 3840, y: 2000)
  Secondary glow: #00FF88 at 5% opacity, upper-left (x: 1000, y: 600)

LAYER 2 -- HEADER BAR
  Agent name: "Pharmacogenomics Intelligence Agent"
  Engine badge: "Precision Intelligence Engine"
  Separator: 4px solid #00B4D8

LAYER 3 -- HERO VISUAL
  Content: CYP enzyme pathway diagram. Central liver illustration with
  5 major CYP enzyme pathways radiating outward (CYP2D6, CYP2C19,
  CYP2C9, CYP3A4/5, CYP1A2). Each pathway shows: enzyme name in a
  hexagonal node, 3-4 major substrate drugs branching off (e.g., CYP2D6:
  codeine, tamoxifen, metoprolol, fluoxetine), and metabolizer phenotype
  labels (Poor -> Intermediate -> Normal -> Rapid -> Ultrarapid) with
  color-coded bars. Star allele notation shown for each gene (e.g.,
  CYP2D6 *1/*4). An arrow from genotype to phenotype to dosing
  recommendation illustrates the translation pipeline.
  Caption: "From Star Alleles to Dosing -- 25 Pharmacogenes, 100+ Drugs"

LAYER 4 -- CAPABILITIES
  Card 1:
    Icon: gene chip
    Title: "25 Pharmacogene Profiling"
    Description: "Comprehensive star allele calling for 25 pharmacogenes
    including CYP2D6, CYP2C19, CYP2C9, CYP3A5, CYP1A2, DPYD, TPMT,
    NUDT15, UGT1A1, SLCO1B1, VKORC1, CYP2B6, CYP4F2, NAT2, and
    11 additional clinically relevant genes."

  Card 2:
    Icon: guideline book
    Title: "CPIC Guideline Integration"
    Description: "Direct integration with Clinical Pharmacogenetics
    Implementation Consortium (CPIC) guidelines. Automated lookup of
    evidence-based dosing recommendations for each gene-drug pair.
    CPIC level A (strong) through D (informative) classifications."

  Card 3:
    Icon: conversion arrows
    Title: "Phenoconversion Detection"
    Description: "Detection of phenoconversion events where concomitant
    medications inhibit or induce CYP enzymes, converting a normal
    metabolizer to a poor or ultrarapid metabolizer phenotype. Real-time
    drug interaction screening with phenotype adjustment."

  Card 4:
    Icon: shield/alert
    Title: "HLA Hypersensitivity Screening"
    Description: "15 HLA-drug hypersensitivity associations including
    HLA-B*57:01 (abacavir), HLA-B*58:01 (allopurinol), HLA-B*15:02
    (carbamazepine), HLA-A*31:01 (carbamazepine), and HLA-B*13:01
    (dapsone). Prevents severe adverse drug reactions (SJS/TEN, DRESS)."

LAYER 5 -- DATA INVENTORY
  Metric 1: "15" / "Milvus Collections"
  Metric 2: "25" / "Pharmacogenes"
  Metric 3: "100+" / "Drug Mappings"
  Metric 4: "9" / "Dosing Algorithms"
  Metric 5: "15" / "HLA Associations"
  Metric 6: "450+" / "Validated Tests"

LAYER 6 -- CROSS-AGENT FLOW
  Step 1: "Genotype" (icon: DNA chip)
  Step 2: "Star Allele" (icon: star)
  Step 3: "Phenotype" (icon: metabolism)
  Step 4: "Drug Selection" (icon: pill)
  Step 5: "Precision Dosing" (icon: dose adjust)
  Connected agents: Biomarker, Oncology, Cardiology, Neurology, Autoimmune

LAYER 7 -- FOOTER
  Ports: "API: :8061 | UI: :8506"
  Tech badges: [Milvus] [Claude] [BGE] [FastAPI] [Streamlit] [CPIC] [PharmGKB]
```

---

#### Agent 7: Cardiology Intelligence Agent

```
AGENT:      Cardiology Intelligence Agent
FILENAME:   hcls-cardiology-intelligence-infographic-8k.png
CANVAS:     7680 x 4320 pixels (8K UHD)
BACKGROUND: radial-gradient(center, #1B2333, #0a0a0f)
GRID:       12-column, 40px gutter, 200px margin
COLOR:      #00B4D8 (Cyan -- Precision Intelligence Engine)

LAYER 1 -- BACKGROUND
  Glow position: center-left (x: 2400, y: 1600)
  Secondary glow: #FF3333 at 5% opacity, center-right (x: 5600, y: 2000)

LAYER 2 -- HEADER BAR
  Agent name: "Cardiology Intelligence Agent"
  Engine badge: "Precision Intelligence Engine"
  Separator: 4px solid #00B4D8

LAYER 3 -- HERO VISUAL
  Content: Anatomical heart diagram (anterior view) with 4 GDMT
  (Guideline-Directed Medical Therapy) pillars annotated as colored
  columns rising from the base. Pillar 1 (cyan): RAAS Inhibitors
  (ACEi/ARB/ARNI -- sacubitril/valsartan, enalapril, losartan). Pillar
  2 (blue): Beta-Blockers (carvedilol, metoprolol succinate, bisoprolol).
  Pillar 3 (teal): MRA (spironolactone, eplerenone). Pillar 4 (aqua):
  SGLT2i (dapagliflozin, empagliflozin). Each pillar has drug names
  and target dose ranges. Surrounding the heart: 6 risk calculator
  badge icons (ASCVD, HEART, CHA2DS2-VASc, HAS-BLED, MAGGIC, EuroSCORE)
  connected by dotted lines showing which calculator feeds which therapy
  decision.
  Caption: "4-Pillar GDMT Optimization -- 14 Drugs, 6 Risk Calculators, 56 Genes"

LAYER 4 -- CAPABILITIES
  Card 1:
    Icon: calculator
    Title: "6 Validated Risk Calculators"
    Description: "ASCVD (10-year atherosclerotic risk), HEART (chest pain
    evaluation), CHA2DS2-VASc (stroke risk in AFib), HAS-BLED (bleeding
    risk on anticoagulation), MAGGIC (heart failure mortality), and
    EuroSCORE II (cardiac surgery risk). All calculable from structured inputs."

  Card 2:
    Icon: 4 pillars
    Title: "4-Pillar GDMT Optimizer (14 Drugs)"
    Description: "Automated GDMT optimization for HFrEF: titration scheduling,
    contraindication checking, renal function monitoring, potassium
    tracking, and evidence-based target dose calculation for all 14
    guideline-recommended medications across 4 therapeutic pillars."

  Card 3:
    Icon: DNA + heart
    Title: "18 Cross-Modal Genomic Triggers"
    Description: "18 imaging-to-genomics trigger rules: LV hypertrophy triggers
    HCM gene panel (MYH7, MYBPC3), prolonged QTc triggers channelopathy
    panel (KCNQ1, KCNH2, SCN5A), aortic dilation triggers Marfan panel
    (FBN1, TGFBR1/2), and 15 additional cross-modal connections."

  Card 4:
    Icon: workflow diagram
    Title: "11 Clinical Workflows"
    Description: "Complete decision support for acute coronary syndrome,
    heart failure (HFrEF/HFpEF/HFmrEF), atrial fibrillation, valvular
    heart disease, pulmonary hypertension, aortic disease, peripheral
    arterial disease, cardiac arrest, syncope, chest pain, and
    cardiomyopathy evaluation."

LAYER 5 -- DATA INVENTORY
  Metric 1: "13" / "Milvus Collections"
  Metric 2: "45" / "Cardiac Conditions"
  Metric 3: "56" / "Cardiac Genes"
  Metric 4: "1,966" / "Validated Tests"
  Metric 5: "44,451" / "Lines of Code"
  Metric 6: "14" / "GDMT Drugs"

LAYER 6 -- CROSS-AGENT FLOW
  Step 1: "Cardiac Imaging" (icon: echocardiogram)
  Step 2: "Risk Score" (icon: calculator)
  Step 3: "GDMT Optimization" (icon: titration)
  Step 4: "Genomic Triggers" (icon: DNA)
  Step 5: "Personalized Plan" (icon: heart + Rx)
  Connected agents: Imaging, Pharmacogenomics, Biomarker, Clinical Trial

LAYER 7 -- FOOTER
  Ports: "API: :8071 | UI: :8507"
  Tech badges: [Milvus] [Claude] [BGE] [FastAPI] [Streamlit] [ACC/AHA] [ESC]
```

---

#### Agent 8: Neurology Intelligence Agent

```
AGENT:      Neurology Intelligence Agent
FILENAME:   hcls-neurology-intelligence-infographic-8k.png
CANVAS:     7680 x 4320 pixels (8K UHD)
BACKGROUND: radial-gradient(center, #1B2333, #0a0a0f)
GRID:       12-column, 40px gutter, 200px margin
COLOR:      #00B4D8 (Cyan -- Precision Intelligence Engine)

LAYER 1 -- BACKGROUND
  Glow position: upper-center (x: 3840, y: 1200)
  Secondary glow: #9966FF at 5% opacity, lower-right (x: 5800, y: 3200)

LAYER 2 -- HEADER BAR
  Agent name: "Neurology Intelligence Agent"
  Engine badge: "Precision Intelligence Engine"
  Separator: 4px solid #00B4D8

LAYER 3 -- HERO VISUAL
  Content: Lateral brain illustration with 10 clinical scale gauge
  indicators arranged in a dashboard layout around it. Each gauge is
  a semicircular meter (like a speedometer) with the scale name, current
  score, and severity zone coloring:
  1. NIHSS (0-42): stroke severity
  2. GCS (3-15): consciousness level
  3. MoCA (0-30): cognitive screening
  4. UPDRS (0-199): Parkinson's severity
  5. EDSS (0-10): MS disability
  6. mRS (0-6): functional outcome
  7. MMSE (0-30): dementia screening
  8. Epworth (0-24): sleepiness
  9. HIT-6 (36-78): headache impact
  10. PHQ-9 (0-27): depression screening
  Each gauge has green/yellow/orange/red zones. Brain regions are
  subtly highlighted corresponding to each scale's primary assessment
  area (e.g., motor cortex for UPDRS, frontal for MoCA).
  Caption: "10 Validated Neurological Scales -- Quantified Clinical Assessment"

LAYER 4 -- CAPABILITIES
  Card 1:
    Icon: gauge cluster
    Title: "10 Validated Clinical Scales"
    Description: "Automated scoring for NIHSS, GCS, MoCA, UPDRS-III, EDSS,
    modified Rankin Scale, MMSE, Epworth Sleepiness Scale, HIT-6, and
    PHQ-9. Longitudinal tracking with trend visualization, milestone
    alerts, and normative population comparisons."

  Card 2:
    Icon: clock/emergency
    Title: "Acute Stroke Triage (<90 Seconds)"
    Description: "Rapid stroke assessment: NIHSS auto-scoring, large vessel
    occlusion prediction, CT perfusion mismatch analysis, thrombolysis
    eligibility (tPA window), thrombectomy candidacy, and door-to-needle
    time optimization. Integration with Imaging Agent for ASPECTS scoring."

  Card 3:
    Icon: brain layers
    Title: "Dementia ATN Staging"
    Description: "Amyloid-Tau-Neurodegeneration (ATN) framework staging for
    Alzheimer's disease. Integration of CSF biomarkers (Abeta42, p-tau,
    t-tau), PET imaging results, and MRI volumetrics. Automated staging
    from A-T-N- (normal) through A+T+N+ (Alzheimer's continuum)."

  Card 4:
    Icon: workflow branching
    Title: "8 Clinical Workflows"
    Description: "Decision support for acute stroke, epilepsy/seizure,
    headache/migraine, movement disorders (Parkinson's, essential tremor),
    multiple sclerosis, dementia evaluation, peripheral neuropathy, and
    neuromuscular disease. Each workflow includes differential diagnosis
    tree and treatment algorithm."

LAYER 5 -- DATA INVENTORY
  Metric 1: "14" / "Milvus Collections"
  Metric 2: "58" / "Neurological Conditions"
  Metric 3: "43" / "Neuro Drugs"
  Metric 4: "10" / "Clinical Scales"
  Metric 5: "208" / "Validated Tests"
  Metric 6: "8" / "Clinical Workflows"

LAYER 6 -- CROSS-AGENT FLOW
  Step 1: "EEG / MRI" (icon: brain scan)
  Step 2: "Clinical Scale" (icon: gauge)
  Step 3: "Differential Dx" (icon: decision tree)
  Step 4: "Treatment Plan" (icon: prescription)
  Step 5: "Monitoring" (icon: trend line)
  Connected agents: Imaging, Pharmacogenomics, Rare Disease, Biomarker

LAYER 7 -- FOOTER
  Ports: "API: :8081 | UI: :8508"
  Tech badges: [Milvus] [Claude] [BGE] [FastAPI] [Streamlit] [NIHSS] [GCS]
```

---

#### Agent 9: Rare Disease Diagnostic Agent

```
AGENT:      Rare Disease Diagnostic Agent
FILENAME:   hcls-rare-disease-infographic-8k.png
CANVAS:     7680 x 4320 pixels (8K UHD)
BACKGROUND: radial-gradient(center, #1B2333, #0a0a0f)
GRID:       12-column, 40px gutter, 200px margin
COLOR:      #00B4D8 (Cyan -- Precision Intelligence Engine)

LAYER 1 -- BACKGROUND
  Glow position: center (x: 3840, y: 1800)
  Secondary glow: #FFD700 at 5% opacity, upper-right (x: 6000, y: 600)

LAYER 2 -- HEADER BAR
  Agent name: "Rare Disease Diagnostic Agent"
  Engine badge: "Precision Intelligence Engine"
  Separator: 4px solid #00B4D8

LAYER 3 -- HERO VISUAL
  Content: Diagnostic odyssey visualization (different from autoimmune
  agent). Left side: family silhouette (parents + child) with a
  winding, frustrating path spanning "5-7 YEARS" -- showing icons for
  ER visits (12), specialists seen (8), misdiagnoses (3), invasive
  tests (6), emotional toll markers (tears, financial burden, school
  absence). Right side: the same family silhouette with a straight,
  short arrow to diagnosis spanning "MINUTES" -- showing HPO term entry,
  gene panel match, ACMG classification, confirmed diagnosis. Large
  central comparison: "Traditional: 5-7 years, 8+ specialists, $50K+"
  vs "AI-Assisted: 1 visit, 1 platform, accessible." Bottom ribbon:
  "300 million people worldwide live with a rare disease."
  Caption: "Ending the Diagnostic Odyssey -- 88 Rare Diseases, 45 Genes, 23 ACMG Criteria"

LAYER 4 -- CAPABILITIES
  Card 1:
    Icon: phenotype tree
    Title: "HPO Phenotype Matching (88 Diseases)"
    Description: "Human Phenotype Ontology (HPO) term-based disease matching
    across 88 rare diseases. Clinicians enter observed phenotypes; the
    system performs semantic similarity scoring against known disease
    phenotype profiles, ranking candidate diagnoses by match confidence."

  Card 2:
    Icon: classification badge
    Title: "ACMG 28-Criteria Variant Classification"
    Description: "Automated variant classification following ACMG/AMP
    standards using 28 evidence criteria (PVS1, PS1-4, PM1-6, PP1-5,
    BA1, BS1-4, BP1-7). Generates structured evidence summaries with
    pathogenic, likely pathogenic, VUS, likely benign, or benign calls."

  Card 3:
    Icon: gene therapy syringe
    Title: "Gene Therapy Eligibility (12 Therapies)"
    Description: "Eligibility assessment for 12 FDA-approved/EMA-authorized
    gene therapies including Zolgensma (SMA), Luxturna (RPE65 retinal
    dystrophy), Skysona (CALD), Elevidys (DMD), Casgevy (SCD/beta-thal),
    and 7 additional therapies. Includes age/weight criteria, prior
    treatment checks, and antibody screening requirements."

  Card 4:
    Icon: workflow hub
    Title: "10 Diagnostic Workflows"
    Description: "Structured workflows for undiagnosed disease evaluation,
    carrier screening, prenatal diagnosis, newborn screening follow-up,
    family cascade testing, variant reclassification, gene therapy
    referral, clinical trial matching, genetic counseling support, and
    multi-disciplinary team presentation preparation."

LAYER 5 -- DATA INVENTORY
  Metric 1: "14" / "Milvus Collections"
  Metric 2: "88" / "Rare Diseases"
  Metric 3: "45" / "Disease Genes"
  Metric 4: "28" / "ACMG Criteria"
  Metric 5: "12" / "Gene Therapies"
  Metric 6: "380+" / "Validated Tests"

LAYER 6 -- CROSS-AGENT FLOW
  Step 1: "Phenotype (HPO)" (icon: clinical features)
  Step 2: "Gene Match" (icon: gene panel)
  Step 3: "ACMG Classification" (icon: evidence tiers)
  Step 4: "Therapy Eligibility" (icon: gene therapy)
  Step 5: "Family Cascade" (icon: family tree)
  Connected agents: Pharmacogenomics, Imaging, Single-Cell, Clinical Trial

LAYER 7 -- FOOTER
  Ports: "API: :8091 | UI: :8509"
  Tech badges: [Milvus] [Claude] [BGE] [FastAPI] [Streamlit] [HPO] [ACMG] [OMIM]
```

---

#### Agent 10: Clinical Trial Intelligence Agent

```
AGENT:      Clinical Trial Intelligence Agent
FILENAME:   hcls-clinical-trial-infographic-8k.png
CANVAS:     7680 x 4320 pixels (8K UHD)
BACKGROUND: radial-gradient(center, #1B2333, #0a0a0f)
GRID:       12-column, 40px gutter, 200px margin
COLOR:      #8b5cf6 (Purple -- Bridges Intelligence Network & Discovery Engine)

LAYER 1 -- BACKGROUND
  Glow position: center (x: 3840, y: 1600)
  Secondary glow: #00B4D8 at 5% opacity, left (x: 1200, y: 2000)
  Tertiary glow: #76B900 at 5% opacity, right (x: 6400, y: 2000)

LAYER 2 -- HEADER BAR
  Agent name: "Clinical Trial Intelligence Agent"
  Engine badge (dual): "Intelligence Network <-> Discovery Engine"
  Separator: 4px solid #8b5cf6

LAYER 3 -- HERO VISUAL
  Content: Protocol optimization pipeline. Left: incoming drug candidate
  molecule (3D ball-and-stick render, purple glow). Center: a multi-
  stage pipeline visualization showing protocol design nodes --
  eligibility criteria optimization, endpoint selection, sample size
  calculation, randomization schema, stratification factors, interim
  analysis rules. Right: enrollment prediction chart (line graph with
  confidence intervals projecting enrollment rate over 24 months,
  comparing optimized vs. traditional enrollment curves). Below pipeline:
  40 landmark trial cards (small, tiled, showing trial name + condition
  + phase) representing the knowledge base. A "decision engine" hub
  in the center connects: genomic matching, site selection, regulatory
  intelligence, competitive landscape, and enrollment prediction.
  Caption: "From Molecule to Approval -- 40 Landmark Trials Powering Protocol Intelligence"

LAYER 4 -- CAPABILITIES
  Card 1:
    Icon: workflow network
    Title: "10 Intelligence Workflows"
    Description: "Protocol optimization, eligibility criteria refinement,
    endpoint strategy, enrollment prediction, site selection, competitive
    landscape analysis, regulatory pathway planning, biomarker strategy,
    safety signal monitoring, and real-world evidence integration. Each
    workflow leverages the 40-trial knowledge base."

  Card 2:
    Icon: landmark flag
    Title: "40 Landmark Trial Knowledge Base"
    Description: "Deep-indexed knowledge base of 40 landmark clinical trials
    across 35 conditions and 34 drugs. Includes protocol design rationale,
    enrollment challenges, endpoint evolution, regulatory interactions,
    and lessons learned. Serves as the evidence foundation for all
    decision support recommendations."

  Card 3:
    Icon: engine/gears
    Title: "5 Decision Support Engines"
    Description: "Eligibility Engine (inclusion/exclusion optimization),
    Enrollment Engine (site-level prediction modeling), Endpoint Engine
    (primary/secondary/exploratory selection), Regulatory Engine (FDA/EMA
    pathway alignment), and Competitive Engine (landscape positioning
    and differentiation analysis)."

  Card 4:
    Icon: molecular bridge
    Title: "Cross-Agent Molecular Matching"
    Description: "Bidirectional integration: receives drug candidates from
    Therapeutic Discovery Engine (DiffDock/MolMIM), matches patients
    from Intelligence Network agents by genomic profile, and feeds
    trial results back to refine future drug design. The bridge between
    bench and bedside."

LAYER 5 -- DATA INVENTORY
  Metric 1: "14" / "Milvus Collections"
  Metric 2: "40" / "Landmark Trials"
  Metric 3: "35" / "Conditions"
  Metric 4: "34" / "Drugs Analyzed"
  Metric 5: "769" / "Validated Tests"
  Metric 6: "5" / "Decision Engines"

LAYER 6 -- CROSS-AGENT FLOW
  Step 1: "Drug Candidate" (icon: molecule)
  Step 2: "Protocol Design" (icon: document)
  Step 3: "Patient Match" (icon: genomic profile)
  Step 4: "Site Selection" (icon: map pin)
  Step 5: "Regulatory Filing" (icon: FDA badge)
  Connected agents: ALL agents (hub role -- Oncology, Biomarker, CAR-T,
    Rare Disease, Cardiology, Neurology, Autoimmune, PGx, Single-Cell, Imaging)

LAYER 7 -- FOOTER
  Ports: "API: :8101 | UI: :8510"
  Tech badges: [Milvus] [Claude] [BGE] [FastAPI] [Streamlit] [ClinicalTrials.gov] [FDA]
```

---

#### Agent 11: Single-Cell Intelligence Agent

```
AGENT:      Single-Cell Intelligence Agent
FILENAME:   hcls-single-cell-infographic-8k.png
CANVAS:     7680 x 4320 pixels (8K UHD)
BACKGROUND: radial-gradient(center, #1B2333, #0a0a0f)
GRID:       12-column, 40px gutter, 200px margin
COLOR:      #76B900 (Green -- Bridges Foundation Engine & Intelligence Network)

LAYER 1 -- BACKGROUND
  Glow position: center-left (x: 2800, y: 1800)
  Secondary glow: #00B4D8 at 5% opacity, right (x: 5600, y: 1200)
  Tertiary glow: #76B900 at 8% opacity, center (x: 3840, y: 2400)

LAYER 2 -- HEADER BAR
  Agent name: "Single-Cell Intelligence Agent"
  Engine badge (dual): "Foundation Engine <-> Intelligence Network"
  Separator: 4px solid #76B900

LAYER 3 -- HERO VISUAL
  Content: UMAP (Uniform Manifold Approximation and Projection) cluster
  visualization. Large 2D scatter plot with ~50,000 dots forming
  distinct cell type clusters. Each cluster colored by cell type
  identity (T cells in green, B cells in blue, macrophages in orange,
  fibroblasts in gray, epithelial in pink, NK cells in purple, dendritic
  cells in yellow, etc. -- 57 cell types total). Overlay: semi-transparent
  colored regions indicating TME (Tumor Microenvironment) classification
  zones -- "HOT" (red border, immune-infiltrated), "COLD" (blue border,
  immune-desert), "EXCLUDED" (yellow border, immune-excluded), and
  "IMMUNOSUPPRESSIVE" (gray border, regulatory T-cell enriched). Side
  panel: bar chart of top 20 marker genes (e.g., CD3E, CD8A, CD19,
  MS4A1, CD68, FOXP3, EPCAM) with expression levels. Bottom ribbon:
  "74x GPU Speedup: 47 minutes -> 38 seconds."
  Caption: "44 Cell Types, 4 TME Classes -- GPU-Accelerated Single-Cell Analysis"

LAYER 4 -- CAPABILITIES
  Card 1:
    Icon: microenvironment grid
    Title: "TME Classifier (4 Categories)"
    Description: "Automated Tumor Microenvironment classification into hot
    (immune-infiltrated, checkpoint inhibitor responsive), cold (immune-
    desert, requires priming), excluded (immune cells at periphery, needs
    stromal remodeling), and immunosuppressive (Treg/MDSC dominant,
    requires depletion). Directly informs immunotherapy selection."

  Card 2:
    Icon: GPU lightning
    Title: "GPU-Accelerated (74x Speedup)"
    Description: "RAPIDS/cuML GPU acceleration for all major single-cell
    analysis steps: quality control, normalization, highly variable gene
    selection, PCA, neighbor graph construction, UMAP, Leiden clustering,
    and differential expression. 47-minute CPU pipeline in 38 seconds
    on DGX Spark."

  Card 3:
    Icon: CAR-T cell + target
    Title: "CAR-T Target Validation"
    Description: "Single-cell resolution validation of CAR-T targets:
    confirms target antigen expression on tumor cells, quantifies on-
    target off-tumor risk by checking expression in normal cell types,
    identifies tumor heterogeneity and antigen-low subpopulations that
    may cause relapse. Direct integration with CAR-T Agent."

  Card 4:
    Icon: 10 workflow nodes
    Title: "10 Analytical Workflows"
    Description: "Quality control, cell type annotation, differential
    expression, trajectory analysis, gene regulatory network inference,
    cell-cell communication, TME classification, drug response prediction,
    multi-sample integration, and spatial transcriptomics preprocessing.
    Each workflow produces publication-ready visualizations."

LAYER 5 -- DATA INVENTORY
  Metric 1: "12" / "Milvus Collections"
  Metric 2: "44" / "Cell Types"
  Metric 3: "30" / "Drug Mappings"
  Metric 4: "75" / "Marker Genes"
  Metric 5: "185" / "Validated Tests"
  Metric 6: "74x" / "GPU Speedup"

LAYER 6 -- CROSS-AGENT FLOW
  Step 1: "scRNA-seq" (icon: sequencer)
  Step 2: "GPU Clustering" (icon: GPU chip)
  Step 3: "TME Classification" (icon: microenvironment)
  Step 4: "Drug Response" (icon: drug target)
  Step 5: "Treatment Strategy" (icon: precision plan)
  Connected agents: CAR-T, Oncology, Clinical Trial, Rare Disease

LAYER 7 -- FOOTER
  Ports: "API: :8111 | UI: :8511"
  Tech badges: [Milvus] [Claude] [BGE] [FastAPI] [Streamlit] [RAPIDS] [cuML] [Scanpy]
```

---

## 14. Website Performance Benchmarks

### 14.1 Target Lighthouse Scores

| Category | Target | Current Estimate | Gap | Priority |
|---|---|---|---|---|
| Performance | > 90 | ~85 | -5 | High |
| Accessibility | > 95 | ~90 | -5 | Critical |
| SEO | > 95 | ~80 | -15 | High |
| Best Practices | > 95 | ~92 | -3 | Medium |

### 14.2 Performance Optimization Plan

**Image Optimization (Primary Performance Bottleneck)**
- Convert all hero images from PNG to WebP (40-60% size reduction)
- Generate AVIF fallbacks for Chrome/Firefox (60-80% size reduction)
- Implement responsive `<picture>` elements with srcset for 3 breakpoints
- Lazy load all images below the fold with `loading="lazy"` and Intersection Observer
- Maximum hero image: 200KB (WebP), 120KB (AVIF)
- Maximum infographic thumbnail: 80KB (WebP)

**CSS Optimization**
- Inline critical CSS (above-the-fold) directly in `<head>` -- target < 14KB
- Defer all non-critical CSS with `rel="preload"` and `onload` swap
- Purge unused CSS from Material for MkDocs (estimated 60% unused in default theme)
- Minify all CSS with cssnano in production build

**JavaScript Optimization**
- Defer all non-critical JavaScript with `defer` attribute
- Code-split agent-specific JavaScript (load only when agent page visited)
- Remove unused Material for MkDocs JavaScript modules
- Target: < 100KB total JavaScript (gzipped)

**Font Optimization**
- Subset Inter font to Latin characters only (save ~200KB)
- Subset JetBrains Mono to ASCII + code characters (save ~150KB)
- Use `font-display: swap` for all web fonts
- Preload primary font weights: Inter Regular (400), Inter Bold (700)

**Core Web Vitals Targets**

| Metric | Target | Strategy |
|---|---|---|
| LCP (Largest Contentful Paint) | < 2.5s | Preload hero image, optimize server response |
| FID (First Input Delay) | < 100ms | Defer JavaScript, minimize main thread work |
| CLS (Cumulative Layout Shift) | < 0.1 | Set explicit image dimensions, reserve font space |
| TTFB (Time to First Byte) | < 200ms | Netlify CDN edge caching, prerender key pages |
| INP (Interaction to Next Paint) | < 200ms | Optimize event handlers, reduce DOM size |

### 14.3 SEO Gap Analysis

**Missing Elements (Causing ~80 Score)**
- No structured data (JSON-LD) for Software Application schema
- No Open Graph meta tags on agent pages
- No Twitter Card meta tags
- Missing canonical URLs on some pages
- No sitemap.xml generation (MkDocs plugin needed)
- No robots.txt customization
- Missing alt text on ~30% of images

**Remediation**
- Add JSON-LD SoftwareApplication schema to homepage
- Add JSON-LD Article schema to each agent page
- Generate sitemap.xml via mkdocs-sitemap plugin
- Add Open Graph and Twitter Card meta to every page via MkDocs meta plugin
- Audit and add alt text to all images
- Add canonical URL to all pages

---

## 15. Content Strategy for GTC Europe

### 15.1 Pre-GTC Content (T-minus 30 to T-minus 1 Days)

**Blog Post: Platform Announcement (T-30)**
- Title: "Open-Source Precision Medicine: 11 AI Agents on a $4,699 Desktop"
- Length: 2,000-2,500 words
- Sections: The problem (access gap), the solution (3 engines), the agents (overview), the hardware (DGX Spark), the impact (patient stories)
- Publish on: hcls-ai-factory.org/blog, cross-post to Medium
- SEO targets: "open source precision medicine," "NVIDIA DGX Spark healthcare," "AI clinical decision support"

**LinkedIn Series: Agent Spotlight (T-14 to T-3)**
- 11 posts, one per day, one per agent
- Format: Infographic image + 200-word description + 3 hashtags
- Posting schedule: 8:00 AM CET (European audience for GTC Europe)
- Hashtags: #PrecisionMedicine #NVIDIAGTC #OpenSource #HealthcareAI
- Each post tags: @NVIDIA @Anthropic (where relevant)
- CTA: "See the live demo at GTC Europe -- Booth [TBD]"

**Twitter/X Thread (T-7)**
- 15-tweet thread with key visuals
- Tweet 1: Hook -- "What if every hospital could afford a genomics lab?"
- Tweets 2-4: 3 engines overview with diagrams
- Tweets 5-12: One tweet per specialty agent (top 8)
- Tweet 13: Hardware story ($4,699)
- Tweet 14: Open source story (Apache 2.0, 200K+ LOC)
- Tweet 15: CTA -- link to repo + GTC Europe mention

**Targeted Outreach (T-14 to T-7)**
- Email announcement to known GTC attendees
- Genomics/bioinformatics subreddit post (r/bioinformatics, r/genomics)
- Hacker News "Show HN" post (coordinate timing with T-5)

### 15.2 During GTC Content (Days 1-4)

**Live Social Coverage**
- Live-tweet from demo sessions with photos/videos (minimum 5 per day)
- LinkedIn post each evening summarizing day's interactions
- Instagram Stories of booth activity (if applicable)

**Physical Materials**
- QR codes on ALL 11 agent infographics linking to hcls-ai-factory.org/agents/[agent-name]
- QR code on engine overview infographic linking to hcls-ai-factory.org
- Business cards: front = name/role, back = QR to GitHub repo + hcls-ai-factory.org
- One-page handout: "HCLS AI Factory at a Glance" -- 3 engines, 8 agents, key stats

**Demo Capture**
- Record every demo session (screen capture + webcam)
- Capture attendee reactions and testimonials (with permission)
- Note common questions for FAQ page
- Collect business cards / contact info for follow-up

### 15.3 Post-GTC Content (T+1 to T+30)

**Immediate (T+1 to T+3)**
- Thank-you posts on LinkedIn and Twitter/X
- "GTC Europe Recap" blog post with photos, key conversations, demo stats
- YouTube upload of best demo recording (edited, 10-15 minutes)
- Embed YouTube video on hcls-ai-factory.org homepage

**Week 1 (T+3 to T+7)**
- Press release: "HCLS AI Factory Debuts at GTC Europe -- Open-Source Precision Medicine Platform Draws [X] Demo Requests"
- Follow-up emails to all contacts collected at GTC
- Blog post: "What GTC Europe Thought About Open-Source Precision Medicine"
  - Include quotes from attendees (with permission)
  - Highlight most-asked questions
  - Announce any partnerships or collaborations initiated

**Month 1 (T+7 to T+30)**
- arXiv paper cross-promotion with GTC visibility
- GitHub repo star campaign (target: 500 stars in first month)
- Community building: Discord or GitHub Discussions setup
- Webinar: "Deep Dive into the HCLS AI Factory" (recorded, posted to YouTube)

---

## 16. Competitive Positioning on Website

### 16.1 Comparison Table

The website should include a clear, honest comparison table on the homepage or a dedicated "Why HCLS AI Factory" page. The table must be factual and defensible.

| Feature | HCLS AI Factory | Epic CDS | Viz.ai | Medidata Rave | Tempus |
|---|---|---|---|---|---|
| **License** | Apache 2.0 (Free) | Proprietary | Proprietary | Proprietary | Proprietary |
| **Hardware** | DGX Spark ($4,699) | Cloud/On-prem | Cloud | Cloud | Cloud |
| **AI Agents** | 8 specialized | Limited modules | Imaging only | Trial only | Genomics + limited |
| **End-to-End Pipeline** | DNA to Drug Candidate | No | No | No | Partial |
| **Genomic Pipeline** | Full (Parabricks) | Via integrations | No | No | Partial |
| **Drug Discovery** | Yes (BioNeMo) | No | No | No | No |
| **Clinical Trials** | Yes (40 landmark) | No | No | Yes | Limited |
| **Open Source** | Yes (200K+ LOC) | No | No | No | No |
| **Self-Hosted** | Yes (air-gapped) | Requires cloud | Cloud required | Cloud required | Cloud required |
| **Data Sovereignty** | Complete | Vendor-dependent | Vendor cloud | Vendor cloud | Vendor cloud |
| **Customizable** | Fully | Limited | Limited | Limited | API only |
| **Cost (Annual)** | $0 (software) | $500K-$2M | $100K-$500K | $200K-$1M | $100K-$500K |

### 16.2 Key Differentiators

**The $4,699 Argument**
- Total cost of entry: $4,699 (DGX Spark) + $0 (software) = $4,699
- Compare: Average hospital spends $500K-$2M/year on clinical decision support
- ROI: Platform pays for itself with a single avoided adverse drug reaction
- Positioning line: "The cost of one hospital bed for one week buys a permanent precision medicine platform"

**The 11-Agent Argument**
- No competitor offers more than 3 specialized clinical AI capabilities
- HCLS AI Factory offers 8 agents covering 8 medical specialties + research
- Each agent is independently validated with hundreds of tests
- Agents collaborate through shared vector database and cross-agent triggers

**The End-to-End Argument**
- Only platform covering DNA to Variant to Interpretation to Drug Discovery to Clinical Trial
- Competitors address fragments: genomics OR imaging OR trials -- never all three
- 3-engine architecture is a category-defining innovation

### 16.3 Positioning Statement

> "HCLS AI Factory is the world's first open-source precision medicine operating system. It is not competing with Epic, Viz.ai, or Tempus -- it is creating a new category: the complete, self-hosted, AI-powered precision medicine platform that any institution can deploy for under $4,000."

### 16.4 Category Creation Strategy

- Never say "competitor to X" -- say "complement to existing systems"
- Position as infrastructure, not application: "the operating system on which precision medicine runs"
- Emphasize what no one else offers: air-gapped deployment, complete data sovereignty, Apache 2.0 freedom
- Target the 95% of hospitals worldwide that cannot afford Tempus or Epic CDS

---

## 17. User Journey Maps

### 17.1 Clinician Journey

```
+------------------------------------------------------------------+
|                    CLINICIAN USER JOURNEY                          |
+------------------------------------------------------------------+
|                                                                    |
|  1. DISCOVERY                                                      |
|     +-- Source: GTC talk, colleague referral, journal article      |
|     +-- Lands on: hcls-ai-factory.org                             |
|     +-- First impression: "This looks like a real clinical tool"   |
|                                                                    |
|  2. PERSONA ROUTING                                                |
|     +-- Sees: "Choose Your Path" hero section                     |
|     +-- Clicks: "I'm a Clinician"                                 |
|     +-- Routed to: /clinician -- specialty agent overview          |
|                                                                    |
|  3. SPECIALTY SELECTION                                            |
|     +-- Sees: 9 clinical agent cards organized by specialty        |
|     +-- Clicks: their specialty (e.g., "Cardiology")              |
|     +-- Arrives: /agents/cardiology-intelligence                   |
|                                                                    |
|  4. CAPABILITY ASSESSMENT                                          |
|     +-- Reads: Agent capabilities (risk calculators, GDMT, etc.)  |
|     +-- Views: Infographic (downloadable)                          |
|     +-- Checks: Clinical validation (test count, accuracy)         |
|     +-- Decision point: "Is this relevant to my practice?"         |
|                                                                    |
|  5. DEMO ENGAGEMENT                                                |
|     +-- Clicks: "Interactive Demo Guide"                           |
|     +-- Reads: Step-by-step clinical scenario walkthrough          |
|     +-- Sees: Sample queries and expected outputs                  |
|     +-- Decision point: "Can I try this myself?"                   |
|                                                                    |
|  6. DEPLOYMENT                                                     |
|     +-- Clicks: "Deploy This Agent"                                |
|     +-- Reads: Hardware requirements (DGX Spark or equivalent)     |
|     +-- Follows: 3-command Docker deployment                       |
|     +-- Opens: localhost:[port] -- first live agent interaction    |
|                                                                    |
|  7. FIRST QUERY                                                    |
|     +-- Enters: Real clinical question                             |
|     +-- Receives: Evidence-based response with citations           |
|     +-- Validates: Against known clinical knowledge                |
|     +-- Outcome: Trust established (or feedback submitted)         |
|                                                                    |
|  8. ADOPTION                                                       |
|     +-- Shares with: Department colleagues                         |
|     +-- Proposes: IT department evaluation                         |
|     +-- Contributes: Clinical feedback via GitHub Issues            |
|     +-- Becomes: Clinical champion for the platform                |
|                                                                    |
+------------------------------------------------------------------+
```

**Key Metrics to Track:**
- Time from landing to specialty selection: target < 30 seconds
- Demo guide completion rate: target > 60%
- Deployment attempt rate: target > 20% of demo guide completers
- First query success rate: target > 90%

### 17.2 Developer Journey

```
+------------------------------------------------------------------+
|                    DEVELOPER USER JOURNEY                          |
+------------------------------------------------------------------+
|                                                                    |
|  1. DISCOVERY                                                      |
|     +-- Source: GitHub Trending, Hacker News, Reddit, Twitter/X    |
|     +-- Lands on: github.com/hcls-ai-factory                      |
|     +-- First impression: "Well-documented, active, testable"      |
|                                                                    |
|  2. README EVALUATION                                              |
|     +-- Scans: README.md for project scope and maturity            |
|     +-- Checks: Star count, last commit date, contributor count    |
|     +-- Reads: 3-command quick start                               |
|     +-- Decision point: "Worth cloning?"                           |
|                                                                    |
|  3. LOCAL DEPLOYMENT                                               |
|     +-- Runs: git clone + docker compose up                        |
|     +-- Waits: ~5 minutes for full stack startup                   |
|     +-- Opens: localhost:8080 (landing page)                       |
|     +-- Verifies: All services healthy (green status indicators)   |
|                                                                    |
|  4. EXPLORATION                                                    |
|     +-- Clicks: Through landing page to each agent                 |
|     +-- Opens: API docs (FastAPI /docs endpoints)                  |
|     +-- Sends: Test API requests via Swagger UI                    |
|     +-- Reviews: Response format, latency, error handling          |
|                                                                    |
|  5. ARCHITECTURE UNDERSTANDING                                     |
|     +-- Reads: Architecture guide on hcls-ai-factory.org           |
|     +-- Reviews: docker-compose.yml for service dependencies       |
|     +-- Examines: lib/hcls_common for shared infrastructure        |
|     +-- Understands: 3-engine pipeline flow                        |
|                                                                    |
|  6. CONTRIBUTION                                                   |
|     +-- Finds: Open issue matching their expertise                 |
|     +-- Forks: Repository                                          |
|     +-- Implements: Fix or feature                                 |
|     +-- Runs: Test suite (pytest)                                  |
|     +-- Submits: Pull request with description                     |
|                                                                    |
|  7. COMMUNITY ENGAGEMENT                                           |
|     +-- Joins: GitHub Discussions / Discord                        |
|     +-- Stars: Repository                                          |
|     +-- Shares: On social media                                    |
|     +-- Becomes: Regular contributor                               |
|                                                                    |
+------------------------------------------------------------------+
```

**Key Metrics to Track:**
- Clone to running: target < 10 minutes
- README bounce rate: target < 40%
- API exploration rate: target > 50% of deployers
- First PR time: target < 2 weeks from first clone

### 17.3 Pharma / Biotech Executive Journey

```
+------------------------------------------------------------------+
|                    PHARMA EXECUTIVE USER JOURNEY                   |
+------------------------------------------------------------------+
|                                                                    |
|  1. DISCOVERY                                                      |
|     +-- Source: GTC Europe demo, industry report, board referral   |
|     +-- Lands on: hcls-ai-factory.org                             |
|     +-- First impression: "Professional, data-driven, credible"    |
|                                                                    |
|  2. PERSONA ROUTING                                                |
|     +-- Sees: "Choose Your Path" hero section                     |
|     +-- Clicks: "I'm in Pharma / Biotech"                         |
|     +-- Routed to: /pharma -- drug discovery + trial focus         |
|                                                                    |
|  3. VALUE PROPOSITION                                              |
|     +-- Reads: "DNA to Drug Candidate in <5 Hours"                |
|     +-- Sees: 3-engine pipeline with drug discovery emphasis       |
|     +-- Reviews: Clinical Trial Agent capabilities                 |
|     +-- Decision point: "Could this accelerate our pipeline?"      |
|                                                                    |
|  4. EVIDENCE GATHERING                                             |
|     +-- Downloads: arXiv paper (technical validation)              |
|     +-- Reads: Clinical Trial Agent white paper                    |
|     +-- Reviews: Competitive comparison table                      |
|     +-- Notes: Open source = no vendor lock-in                     |
|                                                                    |
|  5. INTERNAL SOCIALIZATION                                         |
|     +-- Shares: arXiv paper with clinical development team         |
|     +-- Forwards: Website link to CTO / VP R&D                     |
|     +-- Presents: At internal innovation meeting                   |
|     +-- Request: IT team to evaluate deployment feasibility        |
|                                                                    |
|  6. POC EVALUATION                                                 |
|     +-- IT deploys: On internal GPU server                         |
|     +-- Clinical team: Tests with historical trial data            |
|     +-- Compliance team: Reviews Apache 2.0 license               |
|     +-- Decision point: "Proceed to pilot?"                        |
|                                                                    |
|  7. ENGAGEMENT                                                     |
|     +-- Contacts: Project maintainers via GitHub or email          |
|     +-- Discusses: Custom deployment, support options              |
|     +-- Evaluates: Integration with existing infrastructure       |
|     +-- Outcome: Partnership, sponsorship, or internal adoption    |
|                                                                    |
+------------------------------------------------------------------+
```

**Key Metrics to Track:**
- Time on /pharma page: target > 3 minutes
- White paper download rate: target > 30% of pharma visitors
- arXiv paper download rate: target > 20%
- Contact form submission: target > 5% of pharma visitors

---

## 18. Documentation Quality Standards

### 18.1 Page-Level Requirements

Every documentation page on hcls-ai-factory.org MUST include:

| Element | Requirement | Location |
|---|---|---|
| Title | Clear, descriptive, < 60 characters | H1 + HTML title |
| Description | 1-2 sentence summary, < 160 characters | Meta description + first paragraph |
| Last Updated | ISO 8601 date (YYYY-MM-DD) | Page header or footer metadata |
| Author | Name or "HCLS AI Factory Team" | Page metadata |
| Breadcrumb | Full navigation path | Top of page |
| Table of Contents | Auto-generated from H2/H3 headings | Right sidebar |
| Reading Time | Estimated minutes | Below title |

### 18.2 Agent Documentation Standard

Every agent MUST have the following documentation pages, each meeting the quality bar:

| Document | Purpose | Min. Length | Update Frequency |
|---|---|---|---|
| Index (overview) | First-touch agent description | 500 words | Per release |
| Demo Guide | Step-by-step clinical scenarios | 1,500 words | Per release |
| Deployment Guide | Installation and configuration | 1,000 words | Per release |
| Learning Guide | Educational deep-dive into domain | 3,000 words | Quarterly |
| Architecture Guide | Technical design and data flow | 2,000 words | Per release |
| White Paper | Clinical rationale and evidence | 5,000 words | Biannually |
| Project Bible | Complete reference (all collections, all endpoints, all tests) | 8,000 words | Per release |
| Production Readiness Report | Test results, coverage, known issues | 2,000 words | Per release |

### 18.3 Code Example Standards

- Every code example MUST be tested in CI (copy-paste runnable)
- Every code example MUST include language identifier in fenced code block
- Every API example MUST include: endpoint, method, headers, request body, response body
- Every CLI example MUST include: command, expected output, common errors
- Screenshots MUST be regenerated with each release (no stale screenshots)
- Screenshots MUST include alt text describing what the screenshot shows

### 18.4 Quality Assurance Checklist

Before any documentation page is published:
- [ ] Spell check passed (American English)
- [ ] All links verified (no 404s)
- [ ] All code examples tested
- [ ] All screenshots current
- [ ] Accessibility: heading hierarchy correct (no skipped levels)
- [ ] Accessibility: all images have alt text
- [ ] SEO: meta description present and < 160 characters
- [ ] Mobile: page renders correctly at 375px width

---

## 19. Brand Voice Guidelines

### 19.1 Tone Principles

| Principle | What It Means | Example |
|---|---|---|
| **Confident, not arrogant** | State capabilities clearly; let numbers speak | "8 agents, 200,000+ lines of tested code" NOT "the most advanced platform ever built" |
| **Technical, not jargon-heavy** | Use precise terminology; explain when needed | "Variant classification using 28 ACMG criteria" NOT "leveraging cutting-edge AI/ML paradigms" |
| **Urgent, not alarmist** | Convey the stakes without fear-mongering | "4.6 years average diagnosis time for rare diseases -- we can do better" NOT "patients are DYING because of slow diagnostics!!!" |
| **Inclusive, not preachy** | Welcome all skill levels; never condescend | "Whether you're deploying your first Docker container or your hundredth" NOT "even non-technical people can use it" |

### 19.2 Word Choice

**Always Use:**
- Specific numbers: "8 agents," "539 validated tests," "$4,699"
- Real clinical scenarios: "A cardiologist evaluating a patient with HFrEF..."
- Actual test counts: "Validated with 1,966 automated tests"
- Patient-centered framing: "Faster diagnosis means earlier treatment"
- Active voice: "The agent classifies variants" not "variants are classified by the agent"

**Never Use:**
- "Revolutionary" -- overused, unverifiable
- "Disrupting" / "Disruptive" -- Silicon Valley cliche
- "Game-changing" -- empty superlative
- "Leverage" (as verb) -- corporate jargon
- "Utilize" -- use "use" instead
- "Cutting-edge" / "State-of-the-art" -- show it, don't say it
- "Democratize" -- overused in tech; say "make accessible" instead
- "Synergy" -- never, under any circumstances

### 19.3 Writing Structure

- Lead with the patient impact, then the technology
- First sentence of any section: the "so what" -- why should the reader care?
- Numbers before narrative: "56 cardiac genes analyzed" before explaining what that means
- End with a call to action: what should the reader do next?

### 19.4 Voice Across Personas

| Persona | Adjust Toward | Example |
|---|---|---|
| Clinician | Clinical relevance, workflow fit | "Integrates DAS28, SLEDAI, CDAI, and BASDAI scoring into your existing workflow" |
| Developer | Technical accuracy, API clarity | "RESTful API with OpenAPI 3.0 spec, average response time <500ms, 95th percentile <2s" |
| Pharma Executive | Business impact, ROI, timeline | "Reduces protocol optimization time from weeks to hours, with 40 landmark trial precedents" |
| Patient Advocate | Accessibility, equity, hope | "Every hospital deserves precision medicine -- not just the ones with million-dollar budgets" |

---

## 20. Metrics & Analytics

### 20.1 Analytics Platform Selection

| Option | Pros | Cons | Recommendation |
|---|---|---|---|
| Google Analytics 4 | Free, powerful, industry standard | Privacy concerns, complex setup, cookie banner required | Use for comprehensive tracking |
| Plausible Analytics | Privacy-focused, no cookies, GDPR-compliant, simple | $9/month, less granular | Use as primary -- aligns with open-source ethos |
| Both | Maximum coverage | Slight overhead | Recommended approach |

**Primary: Plausible Analytics** (privacy-first, no cookie banner needed, simple dashboard)
**Secondary: Google Analytics 4** (conversion tracking, audience insights, Google Search Console integration)

### 20.2 Key Performance Indicators (KPIs)

**Acquisition Metrics:**
- Monthly unique visitors (target: 5,000 in month 1, 20,000 by month 6)
- Traffic sources: organic search, direct, social, referral
- Geographic distribution (target: >30 countries in month 3)
- GitHub repo stars (target: 500 in month 1, 2,000 by month 6)

**Engagement Metrics:**
- Average time on site (target: > 3 minutes)
- Pages per session (target: > 3)
- Bounce rate (target: < 50%)
- Agent page views (which agents attract the most interest?)
- Demo guide completion rate (target: > 40%)
- Infographic download rate

**Conversion Metrics:**
- GitHub clone events (via GitHub traffic analytics)
- Docker pull count (if published to Docker Hub)
- arXiv paper downloads
- Contact form submissions
- Newsletter signups (if implemented)

### 20.3 Monthly Reporting Template

Each monthly report should include:
1. **Executive Summary**: Top-line numbers (visitors, stars, clones) with month-over-month change
2. **Top 10 Pages**: By views and by time-on-page
3. **Agent Leaderboard**: Which agents are most visited (reveals market demand signals)
4. **Geographic Heat Map**: Where visitors come from (informs localization priority)
5. **Referral Sources**: Which links/sites drive traffic (informs partnership strategy)
6. **Conversion Funnel**: Landing -> Agent Page -> Demo Guide -> Deploy -> First Query
7. **Action Items**: 3 specific optimizations for the next month based on data

### 20.4 Event Tracking

Implement custom events for:
- "Agent Selected" -- which agent, which persona
- "Demo Guide Started" -- which agent
- "Demo Guide Completed" -- which agent, time spent
- "Infographic Downloaded" -- which agent, which format
- "GitHub Link Clicked" -- from which page
- "Paper Downloaded" -- arXiv link clicks
- "Deploy Command Copied" -- from which deployment guide

---

## 21. Legal & Compliance

### 21.1 License Notice

Every page footer on hcls-ai-factory.org MUST include:

```
Licensed under the Apache License, Version 2.0.
You may obtain a copy at http://www.apache.org/licenses/LICENSE-2.0
```

The GitHub repository MUST include a LICENSE file at the root containing the full Apache 2.0 text. Every source file MUST include the standard Apache 2.0 header comment.

### 21.2 Medical Disclaimer

Every agent page and every clinical content page MUST display the following disclaimer, prominently positioned (not hidden in a footer):

```
MEDICAL DISCLAIMER: The HCLS AI Factory is intended for research and clinical
decision support only. It does not provide medical diagnoses, treatment
recommendations, or clinical advice. All outputs must be reviewed and validated
by qualified healthcare professionals before any clinical application. This
software is not FDA-cleared, CE-marked, or approved by any regulatory authority
for clinical use.
```

### 21.3 Data Privacy Notice

The website MUST include the following data privacy statement:

```
DATA PRIVACY: The HCLS AI Factory website (hcls-ai-factory.org) does not
collect, store, process, or transmit any patient data, protected health
information (PHI), or personally identifiable information (PII). The platform
is designed for self-hosted, air-gapped deployment where all patient data
remains within the deploying institution's infrastructure. No patient data
is sent to any external service, including Anthropic, NVIDIA, or any third party.
```

### 21.4 HIPAA Compliance Statement

```
HIPAA NOTICE: The HCLS AI Factory software, when self-hosted, is designed to
operate within HIPAA-compliant infrastructure. However, HIPAA compliance is the
responsibility of the deploying institution. The software itself is not a covered
entity and does not process PHI in its default configuration. Institutions
deploying the platform with patient data must ensure appropriate administrative,
physical, and technical safeguards per 45 CFR Parts 160 and 164.
```

### 21.5 Intellectual Property Notice

```
PATENT NOTICE: No patents are pending or granted for the HCLS AI Factory
software. The platform is released under the Apache License 2.0, which includes
an express grant of patent rights from contributors to users.

TRADEMARK NOTICE: "HCLS AI Factory" is an unregistered trademark. NVIDIA, DGX
Spark, Parabricks, BioNeMo, and NIM are trademarks of NVIDIA Corporation.
Anthropic and Claude are trademarks of Anthropic, PBC. All other trademarks are
the property of their respective owners.
```

---

## 22. Roadmap Visibility

### 22.1 What to Show on the Website

**Completed Agents (11) -- Full Visibility**

Each completed agent should display a status badge indicating maturity:

| Badge | Criteria | Color |
|---|---|---|
| Production Ready | > 500 tests, all workflows functional, documentation complete | Green |
| Beta | > 200 tests, core workflows functional, documentation in progress | Yellow |
| Alpha | < 200 tests, foundational workflows, minimal documentation | Orange |

Current agent status (recommended badging):

| Agent | Tests | Recommended Badge |
|---|---|---|
| Cardiology Intelligence | 1,966 | Production Ready |
| Clinical Trial Intelligence | 769 | Production Ready |
| Imaging Intelligence | 539 | Production Ready |
| Pharmacogenomics Intelligence | 450+ | Beta |
| CAR-T Intelligence | 420+ | Beta |
| Rare Disease Diagnostic | 380+ | Beta |
| Precision Oncology | 350+ | Beta |
| Precision Autoimmune | 310+ | Beta |
| Precision Biomarker | 280+ | Beta |
| Neurology Intelligence | 208 | Beta |
| Single-Cell Intelligence | 185 | Beta |

**Pipeline Improvements -- Show with Timeline Ranges**
- Nextflow DSL2 pipeline hardening (Q2 2026)
- GPU-accelerated variant calling optimization (Q2 2026)
- Multi-sample joint calling support (Q3 2026)
- Cloud deployment templates (AWS/GCP/Azure) (Q3 2026)

**Foundation Model Integration -- Show as "Exploring"**
- NVIDIA BioNeMo NIM updates (as released by NVIDIA)
- Additional LLM provider support (evaluating)
- Fine-tuned clinical language models (research phase)

### 22.2 What NOT to Show on the Website

The following must remain internal and must never appear on the public website or repository:

- **Internal timelines with specific dates**: Never commit to "Agent X launches June 15" -- use quarters only
- **Unconfirmed partnerships**: No partner logos or mentions until legally cleared
- **Unvalidated clinical claims**: No accuracy percentages, sensitivity/specificity numbers, or clinical outcome claims without peer-reviewed validation
- **Revenue projections or business models**: The platform is open source; any commercial discussion is internal
- **Competitive attack language**: Never disparage competitors by name; only compare features factually
- **Patient data or real clinical cases**: All demos must use synthetic or publicly available data
- **Security vulnerabilities**: Known issues tracked in private security advisories, not public roadmap
- **Individual contributor performance metrics**: Open source is collaborative, not competitive

---

## 23. Risk Assessment

### 23.1 Brand Risks

| Risk | Probability | Impact | Mitigation |
|---|---|---|---|
| **Brand Dilution** | Medium | High | Strict naming hierarchy enforcement. "HCLS AI Factory" is always the parent brand. Engines and agents always include the parent brand in context. Style guide compliance audits quarterly. |
| **Naming Confusion** | High | Medium | Consistent use of "Engine" (infrastructure) vs "Agent" (clinical intelligence) vs "Pipeline" (processing). Glossary on website. All documentation reviewed for terminology consistency before release. |
| **Competitor Response** | Medium | Low | Speed of adoption via open source. First-mover advantage in the "precision medicine OS" category. Community building creates defensible moat. Focus on what competitors cannot replicate: open source + $4,699 hardware. |
| **Regulatory Scrutiny** | Low | High | Prominent medical disclaimer on EVERY page. Never claim FDA clearance or clinical validation beyond what is peer-reviewed. Position as "decision support" not "diagnostic device." Legal review of all public claims. |
| **Clinical Misuse** | Medium | Critical | Clear documentation that platform is for research and decision support only. Prominent disclaimers. No clinical claims without peer review. Training materials emphasize human-in-the-loop requirement. |
| **Open Source Fragmentation** | Low | Medium | Strong contribution guidelines (CONTRIBUTING.md). Clear governance model. Regular releases from a single maintained fork. Community engagement to prevent hostile forks. |
| **NVIDIA Partnership Perception** | Low | Medium | Platform uses NVIDIA hardware but is not an NVIDIA product. Clear attribution without implying endorsement unless explicitly granted. Maintain independence by supporting alternative hardware where possible. |

### 23.2 Technical Risks to Brand

| Risk | Mitigation |
|---|---|
| Demo failures at GTC Europe | Pre-configured backup laptop, offline mode, pre-recorded demo video |
| Website downtime during launch | Netlify CDN (99.99% SLA), static site generation, no server-side dependencies |
| Security vulnerability discovered | Responsible disclosure policy, rapid patch process, security advisory template |
| Stale documentation | Automated docs build in CI, last-updated dates on all pages, quarterly review cycle |

---

## 24. Implementation Checklist

### 24.1 Week 1: Pre-GTC Sprint (Days 1-7)

**Naming & Branding Updates**
- [ ] Update all agent documentation with 3-engine naming framework
- [ ] Replace "pipeline" references with appropriate engine names where applicable
- [ ] Create engine overview pages (Foundation, Intelligence, Discovery)
- [ ] Update landing page hero with 3-engine visual
- [ ] Add persona routing ("Choose Your Path") to homepage

**Infographic Production**
- [ ] Finalize all 11 agent infographic scripts (Section 13.3)
- [ ] Produce 3 engine overview infographics
- [ ] Produce 1 master platform infographic (all 3 engines)
- [ ] Export all infographics: 8K PNG (master), 1080p PNG (web), SVG (scalable)
- [ ] Generate social media crops (LinkedIn 1200x627, Twitter/X 1600x900)

**Website Updates**
- [ ] Update GitHub README with 3-engine framework
- [ ] Add structured data (JSON-LD) to homepage and agent pages
- [ ] Add Open Graph and Twitter Card meta tags to all pages
- [ ] Update navigation to reflect engine hierarchy
- [ ] Add medical disclaimer to all agent pages
- [ ] Push updates to Netlify

**Content Production**
- [ ] Write platform announcement blog post (2,000 words)
- [ ] Prepare 11 LinkedIn post drafts (one per agent)
- [ ] Draft Twitter/X thread (15 tweets)
- [ ] Create "HCLS AI Factory at a Glance" one-pager (printable PDF)

### 24.2 Week 2: GTC Europe (Days 8-14)

**Physical Preparation**
- [ ] Print 11 agent infographics (poster size: 24" x 36")
- [ ] Print 3 engine overview infographics
- [ ] Print QR code cards (500 units) -- each links to hcls-ai-factory.org
- [ ] Print business cards with QR to GitHub repo (200 units)
- [ ] Print one-pager handouts (300 units)

**Demo Preparation**
- [ ] Configure primary demo laptop (DGX Spark)
- [ ] Configure backup demo laptop (standard GPU laptop)
- [ ] Prepare offline demo mode (no internet dependency)
- [ ] Pre-record 10-minute demo video (backup for technical failures)
- [ ] Prepare backup slide deck (20 slides, key screenshots)
- [ ] Test all 8 agents end-to-end on demo hardware

**Social Media Execution**
- [ ] Begin LinkedIn agent spotlight series (Day 8: Agent 1)
- [ ] Post Twitter/X thread (Day 9)
- [ ] Submit Hacker News "Show HN" post (Day 10)
- [ ] Post to r/bioinformatics and r/genomics (Day 10)

### 24.3 Week 3: Post-GTC (Days 15-21)

**Website Restructure**
- [ ] Implement persona-based navigation (/clinician, /developer, /pharma)
- [ ] Add agent comparison table to each engine page
- [ ] Implement search functionality (MkDocs search plugin optimization)
- [ ] Add FAQ page based on GTC Europe questions
- [ ] Embed YouTube demo video on homepage

**Analytics & Measurement**
- [ ] Deploy Plausible Analytics (primary)
- [ ] Deploy Google Analytics 4 (secondary)
- [ ] Configure custom event tracking (agent selected, demo started, etc.)
- [ ] Set up Google Search Console
- [ ] Create first analytics dashboard

**Press & Outreach**
- [ ] Publish GTC Europe recap blog post
- [ ] Send press release to healthcare IT publications
- [ ] Follow up with all GTC Europe contacts (within 48 hours of event end)
- [ ] Post demo recording to YouTube (edited, 10-15 minutes)

### 24.4 Month 2: Foundation Building (Days 22-60)

**SEO Optimization**
- [ ] Generate and submit sitemap.xml to Google Search Console
- [ ] Optimize page titles and meta descriptions for target keywords
- [ ] Build internal linking structure between agent pages
- [ ] Create 3 SEO-optimized blog posts targeting long-tail keywords
- [ ] Monitor and respond to Google Search Console coverage issues

**Performance Optimization**
- [ ] Convert all images to WebP with AVIF fallback
- [ ] Implement lazy loading for below-fold images
- [ ] Inline critical CSS
- [ ] Achieve Lighthouse Performance > 90 on all key pages
- [ ] Achieve Lighthouse Accessibility > 95 on all pages

**Internationalization Preparation**
- [ ] Identify top 5 non-English visitor countries from analytics
- [ ] Evaluate MkDocs i18n plugin
- [ ] Prioritize first translation language (likely: Mandarin, Spanish, or German)
- [ ] Create translation workflow documentation

**Video Content**
- [ ] Record 5-minute "What is HCLS AI Factory?" overview video
- [ ] Record 3 agent deep-dive videos (top 3 by web traffic)
- [ ] Create 60-second social media teaser video
- [ ] Set up YouTube channel with branded banner and playlists

**Community Building**
- [ ] Set up GitHub Discussions (or Discord server)
- [ ] Create CONTRIBUTING.md with clear contribution guidelines
- [ ] Identify and invite 5 potential early contributors
- [ ] Establish monthly community call (virtual, recorded)
- [ ] Create community newsletter (monthly cadence)

### 24.5 Month 3+: Growth & Maturation (Day 61+)

**Content Scaling**
- [ ] Publish biweekly blog posts on clinical AI topics
- [ ] Guest posts on healthcare IT publications
- [ ] Conference talk submissions (AMIA, HIMSS, Bio-IT World)
- [ ] arXiv paper v2 with GTC Europe results and community feedback

**Technical Maturation**
- [ ] Cloud deployment templates (AWS CloudFormation, GCP Terraform, Azure ARM)
- [ ] Helm charts for Kubernetes deployment
- [ ] CI/CD pipeline with automated documentation builds
- [ ] Automated Lighthouse auditing in CI

**Partnership Development**
- [ ] Academic hospital pilot programs (2-3 institutions)
- [ ] NVIDIA partnership formalization (if applicable)
- [ ] Integration partnerships (EHR vendors, LIMS vendors)
- [ ] Grant applications for clinical validation studies

---

## 25. Conclusion

The HCLS AI Factory is more than a software platform -- it is a statement about who deserves access to precision medicine. Every architectural decision, every line of code, every agent, and every test carries an implicit argument: that the best clinical decision support should not be locked behind enterprise contracts, cloud vendor dependencies, or seven-figure annual fees.

The branding must communicate three things simultaneously:

**1. Technical Credibility**

200,000+ lines of tested code across 8 specialized intelligence agents. Three integrated engines spanning genomics, clinical intelligence, and drug discovery. Validated clinical algorithms with thousands of automated tests. Built on NVIDIA's most advanced accelerated computing platform. Integrated with Anthropic's Claude for reasoning, Milvus for vector search, and BioNeMo for molecular modeling. This is not a prototype or a proof of concept -- it is a production-grade platform with the engineering rigor to prove it.

**2. Clinical Impact**

Faster diagnoses: 5-7 year rare disease diagnostic odysseys compressed to minutes. Better drug selection: pharmacogenomic profiling that prevents adverse drug reactions before they happen. Shorter trial timelines: 40 landmark trial precedents powering protocol optimization. Lives saved: every hour shaved off a stroke triage, every avoided drug interaction, every correctly classified variant translates to better patient outcomes. The platform does not exist to demonstrate technology -- it exists to serve patients.

**3. Moral Clarity**

Open source under Apache 2.0 -- anyone can use it, modify it, deploy it, build on it. Runs on a $4,699 desktop -- no cloud dependency, no vendor lock-in, no recurring fees. Available to every hospital, every research institution, every country. The 95% of the world's hospitals that cannot afford Tempus, Epic CDS, or Viz.ai now have an alternative. That is not a marketing message -- it is the entire reason this platform exists.

The naming framework -- **Genomic Foundation Engine**, **Precision Intelligence Engine**, **Therapeutic Discovery Engine** -- creates a clear hierarchy that scales as new agents and capabilities are added. It communicates architectural intent to engineers, clinical scope to physicians, and strategic vision to executives. It positions the platform not as a collection of tools but as a complete precision medicine operating system -- a category that did not exist before this platform created it.

The website at hcls-ai-factory.org is the world's first impression of this platform. It must be simultaneously:
- **Inspirational** for the patient advocate who needs to believe precision medicine can be accessible
- **Technical** for the engineer who needs to trust the architecture before contributing
- **Clinical** for the physician who needs to see validated, evidence-based workflows before adopting
- **Strategic** for the pharma executive who needs to understand the ROI before committing resources

The recommendations in this document -- from color systems to typography, from navigation architecture to persona routing, from infographic production scripts to SEO optimization -- provide the complete framework to achieve that. Every recommendation is grounded in the platform's actual capabilities: real agent counts, real test counts, real clinical algorithms, real data inventories. The branding does not inflate; it translates.

The infographics, when produced following the specifications in Section 13, will become the visual language of the platform. Shareable on LinkedIn, printable for conferences, displayable in hospital corridors, and instantly communicative of each agent's unique clinical value. The consistent 7-layer template ensures brand coherence across all 8 agents while allowing each agent's domain-specific capabilities to stand on their own.

The implementation checklist in Section 24 is designed around the GTC Europe timeline but serves as the operational blueprint for the platform's public launch regardless of the specific event. Week 1 establishes the foundation. Week 2 captures attention. Week 3 converts attention to engagement. Month 2 builds the infrastructure for sustained growth. Month 3 and beyond scales the community and deepens clinical validation.

This is not just branding. This is the presentation layer of a platform that will change how medicine is practiced, who has access to the best diagnostics, and how fast new therapies reach patients who need them. The HCLS AI Factory's code speaks for itself -- 200,000+ lines of it. The branding ensures that the right people hear what the code is saying.

Get it right, and the HCLS AI Factory becomes synonymous with accessible precision medicine. Not because of clever marketing, but because the platform delivers what the branding promises: precision medicine for everyone, on hardware anyone can afford, with code anyone can inspect.

That is the legacy. Build accordingly.

---

## References

1. Material for MkDocs documentation -- squidfunk.github.io/mkdocs-material -- Configuration reference for theming, navigation, plugins, and deployment
2. NVIDIA Brand Guidelines 2024 -- Partner usage guidelines for NVIDIA, DGX, Parabricks, BioNeMo, and NIM trademarks
3. Nielsen Norman Group -- Navigation patterns for documentation sites -- Research on sidebar vs. top-nav, breadcrumbs, and progressive disclosure
4. Web Content Accessibility Guidelines (WCAG) 2.1 -- w3.org/WAI/WCAG21 -- Level AA conformance targets for color contrast, keyboard navigation, and screen reader compatibility
5. Google Lighthouse scoring methodology -- web.dev/performance-scoring -- Scoring algorithms for Performance, Accessibility, SEO, and Best Practices audits
6. AMP/ASCO/CAP Joint Consensus Recommendation -- Standards and guidelines for the interpretation and reporting of sequence variants in cancer (2017)
7. ACMG/AMP Standards and Guidelines -- Standards and guidelines for the interpretation of sequence variants (Richards et al., 2015)
8. CPIC Guidelines -- Clinical Pharmacogenetics Implementation Consortium guidelines for gene-drug pairs (cpicpgx.org)
9. Human Phenotype Ontology (HPO) -- Standardized vocabulary of phenotypic abnormalities (hpo.jax.org)
10. Plausible Analytics documentation -- plausible.io/docs -- Privacy-focused web analytics implementation guide

---

*HCLS AI Factory Unified Branding Research Guide v1.0.0 | March 2026 | Adam Jones | Apache 2.0*
