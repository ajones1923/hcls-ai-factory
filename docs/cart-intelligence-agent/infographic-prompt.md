# Nano Banana Pro — HCLS AI Factory CAR-T Intelligence Agent on NVIDIA DGX Spark Infographic

## IMPORTANT: Read this entire prompt before generating. This describes a single technical architecture infographic — NOT a slide deck. Every element below appears on ONE canvas. Keep text SHORT inside every box — prefer 2-3 word labels and single-line descriptions. Avoid paragraphs inside boxes. Use icons and arrows to carry meaning instead of text.

---

## OVERALL LAYOUT AND STYLE

Create a professional technical architecture infographic in landscape orientation (16:9 aspect ratio). Clean, structured, authoritative — matching enterprise white paper aesthetics.

**Canvas:** White background (#FFFFFF). Organized with clear section boundaries and generous spacing between components.

**Typography:**
- Title: Large, bold, sans-serif (Inter or Helvetica), deep navy (#1B2333)
- Subtitle: Smaller, medium gray (#666666)
- Section headers: Bold, navy (#1B2333) with NVIDIA green (#76B900) left-border accent
- Component labels: Medium bold, 2-5 words maximum per label
- Body text: 10pt minimum, dark gray (#333333), maximum 2 short lines per box
- Metric callouts: Bold inside rounded pill badges

**Color Palette:**
- NVIDIA Green: #76B900 — NVIDIA components, infrastructure, metric badges
- Deep Navy: #1B2333 — titles, dark bars
- Teal: #1AAFCC — data flow, retrieval components
- Purple: #8B5CF6 — Claude / LLM components
- Amber: #F5A623 — knowledge graph, query expansion
- Emerald: #059669 — comparative mode, outputs
- Light Gray: #F5F5F5 — card backgrounds
- White: #FFFFFF — canvas, text on dark backgrounds

**Visual Elements:**
- Rounded-corner rectangles (8px radius)
- Thin-line monochrome icons (not emoji)
- Color-coded arrows: gray (data flow), teal (retrieval), emerald (comparative)
- Metric badges: small rounded pills, white text on colored background
- NVIDIA logo in title bar and infrastructure footer
- No VAST branding anywhere

---

## CANVAS STRUCTURE (Top to Bottom, 7 bands)

### ━━━ BAND 1: TITLE BAR ━━━

**Left badge:** "CAR-T Intelligence Agent" in green (#76B900) pill

**Center:**
- **Title:** "CAR-T Intelligence Agent"
- **Subtitle:** "Cross-Functional CAR-T Intelligence on NVIDIA DGX Spark"
- **Tagline:** "GB10 Superchip | 128 GB Unified Memory | $3,999"

**Right — Legend** (compact):
```
● Literature (5,047)  ● Trials (973)  ● Constructs (6)
● Assays (45)  ● Manufacturing (30)  +6 more collections
→ Data Flow   ⇒ Comparative Mode
```

---

### ━━━ BAND 2: DATA SOURCES (left column) ━━━

Vertical stack of 5 small cards with icon + bold label + one-line detail:

1. **PubMed** [journal] — 5,047 abstracts via NCBI E-utilities
2. **ClinicalTrials.gov** [clipboard] — 973 trials via API v2
3. **FDA Products** [shield] — 6 approved CAR-T therapies
4. **Landmark Papers** [star] — 45 assay records (ELIANA, ZUMA-1, KarMMa)
5. **Manufacturing** [factory] — 30 CMC/process records

Arrows flow right labeled: "fetch → parse → embed → store"

---

### ━━━ BAND 3: DGX SPARK PLATFORM (center, largest section) ━━━

Enclosed in a container with NVIDIA green border.

**Header bar (green #76B900, white text):** "NVIDIA DGX Spark — GB10 | 128 GB Unified Memory"

#### Layer 3A: DATA LAYER

**Left label bar:** "① Data Layer"

**Box: Embedding** [teal border]
- "BGE-small-en-v1.5"
- "384-dim vectors"
- Badge: "33M params"

**Box: Milvus 2.4** [green border, largest box in this layer]
- "Vector Database — 11 Collections"
- 10 color-coded mini-cards in two rows, each showing only collection name and count:
  - Row 1: Literature: 5,047 | Trials: 973 | Constructs: 6 | Assays: 45 | Manufacturing: 30
  - Row 2: Safety: 40 | Biomarkers: 43 | Regulatory: 25 | Sequences: 27 | Real-World: 30
- Small note: "+1 read-only: Genomic Evidence (3.56M)"
- Below: "IVF_FLAT | COSINE | nlist=1024"
- Badge: "6,266 owned vectors"

#### Layer 3B: EXECUTION LAYER

**Left label bar:** "② Execution Layer"

Flow runs left-to-right with arrows connecting each box:

**Box: User Query** [green border]
- One example: "Why do CD19 CAR-T therapies fail?"

**Box: Comparative Detection** [emerald border]
- "Auto-detect X vs Y"
- Two arrows: "YES → Dual Path" / "NO → Standard"

**Box: Query Expansion** [amber border]
- "169 keywords → 1,496 terms"
- Badge: "12 maps"

**Box: Knowledge Graph** [amber border]
- "25 targets | 15 biomarkers | 8 toxicities"
- Badge: "70+ entities"

**Box: Parallel Search** [green border, wide]
- "11 collections searched simultaneously"
- 11 small colored squares in a row (one per collection)
- Badge: "12-16 ms"

**Box: Comparative Mode** [emerald border, below standard path]
- "Dual retrieval: Entity A vs Entity B"
- Badge: "~365 ms"

#### Layer 3C: SYNTHESIS LAYER

**Left label bar:** "③ Synthesis Layer"

**Box: Claude Sonnet 4.6** [purple border, largest]
- "Anthropic API"
- "Streaming RAG + citations"
- Badges: "~24 sec" | "~30 sec comparative"

**Box: Report Generation** [teal border]
- "Markdown | JSON | PDF"
- "Clickable PMID + NCT links"

---

### ━━━ BAND 4: OUTPUT MODES (right-center column) ━━━

4 stacked cards [emerald border], each with bold label + one line:

1. **RAG Response** — Grounded narrative with citations
2. **Comparative Analysis** — Side-by-side tables + pros/cons
3. **Knowledge Graph** — Interactive network visualization
4. **Image Verification** — Claude Vision claim checking

---

### ━━━ BAND 5: INTERFACES (right edge) ━━━

3 destination cards [purple border]:

1. **Streamlit UI** — NVIDIA dark theme, port 8521
2. **FastAPI REST** — 7 endpoints, port 8522
3. **Export** — Markdown, JSON, PDF

---

### ━━━ BAND 6: HCLS AI FACTORY INTEGRATION (bottom strip) ━━━

Light indigo background (#E0E7FF).

**Header:** "HCLS AI Factory — Cross-Pipeline Integration"

5 horizontal boxes, each with bold label + one line:

1. **Genomics** — Shared Milvus, Parabricks 4.6
2. **RAG/Chat** — Shared embeddings + collections
3. **Drug Discovery** — BioNeMo NIMs
4. **Imaging Agent** — Cross-modal enrichment (future)
5. **Biomarker Agent** — CRS prediction (future)

Dashed arrows connect upward to DGX Spark container.

---

### ━━━ BAND 7: NVIDIA INFRASTRUCTURE BAR (bottom) ━━━

Full-width bar, NVIDIA green (#76B900) background, white text. NVIDIA logo on left.

5 columns, each with bold header + 2 short lines:

| DGX Spark | Milvus 2.4 | BGE-small | Claude Sonnet 4.6 | Docker Compose |
|---|---|---|---|---|
| GB10 Superchip | IVF_FLAT / COSINE | 384-dim vectors | Streaming RAG | 6 services |
| 128 GB, $3,999 | 11 collections | Asymmetric encoding | ~25 sec e2e | Ports 8521, 8522 |

---

## ANNOTATIONS (keep to 3 maximum)

**Badge 1** (top-left of DGX section):
- "5 CAR-T Development Stages in Every Query"

**Badge 2** (center-right floating):
- "6 FDA-Approved Products: Kymriah, Yescarta, Tecartus, Breyanzi, Abecma, Carvykti"

**Badge 3** (bottom-right):
- "Open Source — Apache 2.0 | github.com/ajones1923/cart-intelligence-agent"

---

## TEXT DENSITY RULES

- **Maximum 2 lines of body text per box** — if you need more, use a badge instead
- **Labels are 2-5 words** — never full sentences inside component boxes
- **Tables have maximum 2 data rows** (header + 1-2 rows)
- **Metric badges replace inline numbers** — put numbers in pills, not in body text
- **No code blocks, no file paths, no directory trees** on the canvas
- **Prefer icons + arrows over text** to show relationships

---

## WHAT THIS DIAGRAM MUST COMMUNICATE AT A GLANCE

1. Single NVIDIA DGX Spark ($3,999) runs everything
2. 5 data sources feed 11 Milvus collections (6,266+ vectors, 10 owned + 1 read-only)
3. Three-layer architecture: Data → Execution → Synthesis
4. Knowledge graph (70+ entities) + query expansion (12 maps, 1,496 terms) enrich every search
5. Comparative mode auto-detects "X vs Y" queries
6. Claude Sonnet 4.6 generates grounded answers with citations (~25 sec)
7. Part of the broader HCLS AI Factory ecosystem

The overall impression: a complete CAR-T intelligence platform on a single desktop supercomputer — clean enough to read every label at a glance.

---

*End of Prompt*
