# hcls-ai-factory.org — Site Upgrade Plan

*A detailed, honesty-first approach to make the public site amazing over the next ~2 months, in
lockstep with building the factory up to what the site claims (8 engines · 8 agents · TSC).*

Status: **planning** · Owner: Adam · Horizon: ~8 weeks · Source of truth for capability status:
`lib/hcls_common/capabilities.json` (the registry) — the site must never claim ahead of it.

---

## 0. The one idea that resolves everything

You're about to point a lot of interested people at this site **while the factory is still maturing
toward the roster the site advertises.** The naive move — make it look finished — is exactly the trap
the honesty discipline exists to prevent, and it's what experts and skeptics punish hardest.

So invert it. **Make the site's integrity its wow.** The centerpiece is a *live capability-maturity
matrix* rendered from the registry, where every engine/agent/program carries an honest badge —
`live · verified · preclinical · research-use · roadmap · gated` — and as you build over the two
months, **the site visibly goes green in public.** "Amazing" and "honest" stop fighting and become
the same axis. Nobody else in this space leads with machine-verified honesty. That is the hook.

Three consequences follow from this idea, and they organize the whole plan:
1. **The site renders *from* the registry**, so it is self-honest by construction and cannot drift.
2. **"Roadmap" is not a weakness to hide — it's the story of the two months**, shown with dignity.
3. **Every "wow" element must survive an expert's scrutiny**, because the same polish that thrills the
   public is what a clinician or skeptic inspects first (see `broad-general-persona`).

---

## 1. What the site is *for* (the single job)

The site does not fulfill the mission — no page can. Its job is to **convert a wave of attention into
the four things that actually compound:**

- **Forks & clones** (builders) — "you can run this yourself" is the one claim a competitor can't dismiss.
- **Expert trust** (clinicians, geneticists, oncologists, pharma R&D) — converts on *honesty alone*;
  their public pushback is expensive to earn back, so never give them a gap between claim and reality.
- **Collaborators & institutions** (the MONAI model: credibility = who uses it + citations + governance).
- **Dignified true hope** (patients/families) — the highest-stakes audience; restraint is the design.

Everything below is judged by whether it moves one of these — not by pageviews.

### The five personas → what each needs on the site
| Persona | Lands wanting | Give them (section) | Never |
|---|---|---|---|
| AI-curious public | one legible "wow" + "it's yours" | Hero + mission cut + short honest demo | infer "cures disease now" |
| Domain experts | precision, tools running, limits stated | Engine/Agent pages, citations, "run it yourself", registry | spectacle outrunning evidence |
| Builders / OSS | a forkable, documented repo | Quickstart, hardware, reproducibility, Apache-2.0 | polish over a non-reproducible repo |
| Patients / families | concrete, true, understated hope | Mission page, TSC story (dignified) | hope the project can't deliver |
| Skeptics / competitors | something to attack | Openness, maturity matrix, "fork it and run it" | any claim/reality gap |

**Two-cut principle, applied as IA (not one compromised page):** depth for experts *and* clarity for
the public coexist because the information architecture routes them to different pages — the mission
cut and the technical cut already model this; extend it site-wide.

---

## 2. Where the site is today (baseline, July 2026)

**Good bones, honest voice, but thin and generic.**
- **Stack:** GitHub → MkDocs (Material) → Netlify → hcls-ai-factory.org; `mkdocs build --strict`,
  deploy-on-merge, Python 3.12, pinned `requirements-docs.txt`. Solid; keep it.
- **Published:** exactly **2 pages** — `Home` (mission line, 8/8/1 cards, "What it is — honestly")
  and `Capability Brief` (technical + mission cuts, an `architecture.svg`). Everything else in
  `docs/` is intentionally excluded from the build.
- **Design:** *stock* Material — default teal / deep-orange palette, `material/dna` logo icon,
  dark/light toggle. No custom CSS, no wordmark, **no OG/social share image**, no favicon set,
  no custom typography.
- **Missing:** per-engine pages, per-agent pages, a TSC deep-dive, the D1–D7 demo portfolio,
  a "run it yourself" quickstart, any published diagrams/infographics, a live status/maturity
  matrix, a citations/evidence page, a governance page, `llms.txt`, analytics, a roadmap/changelog.
- **Latent goldmine (unpublished in `docs/`):** white paper (DGX Spark), learning guides
  (foundations/advanced), demo guide, deployment guide, project bible, architecture-diagram research,
  mindmap, **30+ Nano-Banana infographic prompts**, pediatric demo workflows. Raw material for ~80%
  of the new pages already exists — it needs curating and honesty-labeling, not writing from scratch.
- **Where the depth actually lives today:** two hand-authored, **monolithic HTML files** —
  `brief/capability-brief.html` (521 lines, technical cut) and `brief/mission.html` (206 lines,
  mission cut) — polished and print-ready, but **un-themed, not MkDocs-searchable, not in the nav,
  with ports and honesty statuses hard-coded inline** and only manually reconciled to the registry.
  They're a real asset (great content) *and* a liability (a second, drift-prone source of truth).

### 2.1 Known defects to fix first (some are honesty issues, not polish)
- **Homepage roster bug:** the "8 Engines" card on `index.md` lists only **seven** domains and
  mislabels the count in prose. A stale/short count on the front page is exactly the drift the canon
  guards against — **fix in Phase 0.**
- **Stale roster in excluded files:** several `infographic_prompts/*` say
  "Three_Engine_Architecture_**11_Agents**" (old 3·11 count) — inconsistent with 8·8·1. Purge or
  re-label before any of it gets published; also delete the stray `*.md.bak`.
- **No registry↔site drift guard:** nothing checks that the site's `live`/`preclinical`/`planned`
  labels match `capabilities.json`; ports/statuses in the HTML must be hand-synced. Automate this.
- **CI doesn't build the docs:** `ci.yml` runs lint/tests/registry only — the site is gated solely
  by Netlify's `--strict`. Add a docs-build job so a broken site fails a PR *before* Netlify.
- **OG is one plugin away:** mkdocs-material ships a **`social` plugin** that auto-generates OG cards
  — enabling it (+ `pillow`, `cairosvg` in `requirements-docs.txt`) closes the biggest sharing gap fast.

---

## 3. Target information architecture (the sitemap)

Grow from 2 pages to a real, drill-down structure. Every page carries: a one-line "what it is," an
**honest status badge**, real inputs → outputs, "how to run / how to verify," and citations.

```
Home                          upgraded hero: thesis + live maturity summary + "run it yourself"
The Factory/
  Overview                    the Platform → Engine → Agent → Program model; architecture diagram
  8 Engines/                  index (matrix) + one page each (E1…E8): what · status · I/O · run · cite
  8 Intelligence Agents/      index (matrix) + one page each (A1…A8): domain · status · reasoning · cite
  One Platform                registry · MCP tool-surface · composer · MLOps · governance gates
Disease Programs/
  TSC (flagship)              the weight→compression→hope vertical; how the whole factory converges
Demos/                        D1–D7 portfolio (the proving ground) — each with honesty flags
Run It Yourself/              quickstart · hardware (DGX Spark + elastic burst) · reproducibility · datasets
Honesty & Governance/
  Maturity Matrix             THE centerpiece — live, from the registry
  Honesty Ledger              preclinical/roadmap/gated/research-use, stated plainly
  Decision-support posture    never autonomous diagnosis; regulatory stance
  Citations & Evidence        every clinical claim → primary source; verified-on-real-data proofs
Capability Brief              keep (technical + mission cuts) — the deep single-doc read
About / Mission               the "why," with dignity (patients/families land here)
Roadmap / Changelog           the two-month build, in public
```

Machine layer (see §6): `llms.txt`, `sitemap.xml`, per-page OG + schema.

---

## 4. Design system — beyond stock Material (make it beautiful *and* trustworthy)

The goal is "considered clinical," not "AI-generated startup." Move deliberately off the defaults.

- **Palette:** pick a real, chosen palette — a deep clinical navy/ink as ground, one warm accent for
  hope (not the default deep-orange), a restrained teal for structure, and a semantic status set that
  is *separate from the accent*: `live=green · verified=green+check · preclinical=amber ·
  research-use=blue · roadmap=slate · gated=violet`. Neutral greys biased slightly toward the accent.
- **Typography:** pair a characterful display face (headlines/hero) with a clean body face and a mono
  for code/data. Inline as `@font-face` data URIs or self-host — **do not** rely on a font CDN
  (Netlify is fine, but keep the build reproducible and offline-safe). Set a type scale and hold it.
- **Identity:** a real wordmark/logo (beyond the stock DNA glyph), a favicon set, consistent iconography.
- **Custom CSS + overrides:** add `extra_css:` and a Material `overrides/` theme dir; style the hero,
  the card grid, the status badges, and the matrix. Keep it tasteful — most pages are documents, not
  billboards.
- **Real visuals, published:** turn the best of the 30+ infographic prompts into actual assets via the
  `nanobanana` tool (already installed), and publish the architecture diagram. Tag every rendered
  image `real data · derived · illustrative` per the honesty rule — a gorgeous figure that implies a
  specific patient's truth is exactly what an expert catches.
- **The maturity matrix as a visual language:** the status badges are the site's signature UI. Design
  them once, use them everywhere (engine cards, agent cards, demo flags, the matrix).
- **Motion with restraint:** a subtle hero reveal and the matrix filling in; respect
  `prefers-reduced-motion`. No gratuitous animation — it reads as AI-generated and undercuts trust.

---

## 5. The trust/credibility layer (what actually converts experts & skeptics)

This is where the site earns the room. Build these as first-class, not footnotes:

- **"Run it yourself," front and center.** A real quickstart (clone → the gate → a real run). The
  single most disarming thing on the site. Show it *running*, not just its results.
- **Verified-on-real-data badges.** Where a capability is proven against real input (e.g. E1's
  VAF/mosaicism + QC verified on real HG002), say so and link the recorded output. This is the
  difference between `live` and `verified` on the matrix.
- **Citations for every clinical claim.** CPIC, ClinVar, AlphaMissense, ACMG SF v3.3, the EXIST
  trial (everolimus), CAC/ASCVD guidelines — each claim traces to a primary source. A dedicated
  Citations page + inline footnotes. "If a claim can't be cited, it doesn't go on screen."
- **License clarity.** Apache-2.0, commercial use explicitly welcome (open ≠ non-commercial).
- **The registry, published.** Show that the site's statuses come from `capabilities.json` — the
  honesty is mechanical, not marketing.
- **Governance & lineage.** The governed-app gates + 21 CFR Part 11 lineage, shown as design intent.
- **Reproducibility.** Pinned versions, the **public-datasets manifest with NF-9 version pins**
  (the central store you just built) — "here's exactly what data, at exactly which version."

---

## 6. AI-discoverability (a real 2026 differentiator, done honestly)

AI agents increasingly *are* the first visitor. In a regulated health vertical, `llms.txt` is how you
route them to your **honest, canonical** surfaces instead of hype — which is perfectly on-brand here.

- **`llms.txt` at the site root:** H1 (project name) + a 1–3 sentence blockquote that becomes the
  model's mental model ("open-source precision-medicine platform; decision support, not diagnosis;
  8 engines, 8 agents, TSC; honest maturity labels"), then a curated link map to the canonical pages
  (Overview, each Engine/Agent, TSC, Run It Yourself, Honesty & Governance, the registry). Optionally
  an `llms-full.txt` with the expanded brief. *(Adopted by OpenAI, Anthropic, Stripe, Cloudflare,
  Vercel; Chrome Lighthouse 13.3 now audits for it.)*
- **Clean Markdown + MCP:** the factory already exposes an **MCP tool-surface** — advertise it, so an
  agent can not only read the site but discover the governed tools. This is a genuinely unusual,
  credible flex.
- **Structured metadata:** per-page `<title>`/description, canonical URLs, JSON-LD
  (SoftwareApplication / Dataset / MedicalWebPage as appropriate), OG/Twitter cards, `sitemap.xml`.

---

## 7. The 2-month execution plan (phased, in lockstep with the capability build)

The organizing rule: **the site grows green as the factory matures.** Each phase both ships site work
*and* flips badges as real capabilities come online. Nothing here needs RunPod/GPUs to start.

### Phase 0 — Foundation & the honesty spine (Week 0–1)
The fastest path to "this looks intentional," and it sets the honest baseline (mostly roadmap — fine).
- [ ] **Fix the homepage roster bug** (7→8 engine domains + correct count) and purge the stale
      `3·11` infographic prompts + the `*.md.bak`. Honesty defects go first.
- [ ] Design system v1: palette, typography, `extra_css` + `overrides/`, status-badge components.
- [ ] **Enable mkdocs-material `social` plugin** (+ `pillow`, `cairosvg`) → auto OG cards; add a
      favicon set + per-page SEO metadata (closes the biggest sharing gap immediately).
- [ ] **`llms.txt`** + `sitemap.xml` + JSON-LD.
- [ ] **The live maturity matrix**, rendered from `capabilities.json` (a small `mkdocs-gen-files` /
      macro script → a Markdown table/grid). This is the centerpiece; ship it early.
- [ ] IA skeleton: create the nav + stub pages for all sections (honest "roadmap" placeholders).
- [ ] **Registry-driven status + ports:** make pages read ports/statuses from `capabilities.json`
      (via the gen-files script) so the two HTML cuts stop being a second source of truth.
- [ ] A CI drift-guard: site status labels must match the registry (extend `validate_registry.py`
      or add a docs check) **and add a `mkdocs build --strict` job to `ci.yml`** so a broken site
      fails the PR before Netlify.

### Phase 1 — Depth: engines, agents, platform (Weeks 2–4)
- [ ] Per-engine pages (E1–E8) and per-agent pages (A1–A8), each curated from the latent `docs/`
      content + the roster canon, with honest badges. Flip badges to `live`/`verified` as each is proven.
- [ ] One Platform page (registry · MCP · composer · MLOps · governance).
- [ ] Publish the architecture diagram + the first batch of infographics (nanobanana), provenance-tagged.
- [ ] TSC flagship page — the weight→compression→hope vertical, dignified, decision-support framed.

### Phase 2 — Proof: run it yourself + evidence (Weeks 5–7)
- [ ] "Run It Yourself" quickstart (clone → gate → a real run), hardware + elastic-burst honesty.
- [ ] Recorded runs / screenshots / asciinema of real capabilities; verified-on-real-data proofs.
- [ ] Citations & Evidence page (every clinical claim → primary source).
- [ ] Governance page (gates + Part 11 lineage); reproducibility + public-datasets manifest (NF-9).
- [ ] Demos section (D1–D7) with per-demo honesty flags.

### Phase 3 — Launch polish (Week 8)
- [ ] Performance + accessibility pass (Lighthouse, keyboard focus, reduced-motion, alt text).
- [ ] Analytics (privacy-respecting) to measure the *right* conversions (forks/clones/inbound).
- [ ] A launch moment: the mission cut + a short, honest demo video (two-cut: technical + mission),
      share cards, and a coordinated post. Hold the honesty line hardest exactly here.
- [ ] Final README ↔ site ↔ registry three-way reconciliation.

---

## 8. Keeping it honest & in-sync (site governance)

- **Three-way sync** on every framing/port/status change: `README.md` ↔ `docs/` (site) ↔ registry.
  The merge republishes; the CI drift-guard fails the build on mismatch (Pillars 15 & 16).
- **Deploy-on-merge, `--strict`, fail on broken links** (already configured — keep).
- **Neutral repo** (pre-commit guard): no vendor/alternate-edition branding in any tracked site file.
- **Honesty gate before merge:** every new page checked against the honesty ledger; anything not
  plainly `live` is labeled; every clinical claim cited. When drama and accuracy conflict, accuracy wins.

---

## 9. How you'll know it's amazing (metrics that aren't vanity)

- GitHub **stars / forks / clones** trend (builders acted).
- **Expert & collaborator inbound** (the MONAI signal — who reaches out, who cites, who forks).
- **Depth-page engagement** (experts reading engine/agent/TSC pages, not just the hero).
- **`llms.txt` / agent hits** and how AI assistants summarize the project back (test it).
- Qualitatively: **can a skeptical expert find nothing that overclaims?** That's the bar.

---

## Appendix A — Latent `docs/` content → target pages (reuse, don't rewrite)
| Existing (unpublished) | Feeds |
|---|---|
| `HCLS_AI_FACTORY_WHITE_PAPER_DGX_SPARK.md` | The Factory / Overview + Run It Yourself (hardware) |
| `HCLS_AI_FACTORY_LEARNING_GUIDE_*` | Per-engine/agent explainers; "learn the foundations" track |
| `HCLS_AI_FACTORY_DEMO_GUIDE.md`, pediatric demo workflows | Demos (D1–D7) + TSC flagship |
| `HCLS_AI_FACTORY_ARCHITECTURE_DIAGRAM_RESEARCH.md`, mindmap | Overview architecture diagram |
| `infographic_prompts/*` (30+) | Published figures via `nanobanana` (provenance-tagged) |
| `HCLS_AI_FACTORY_DGX_SPARK_DEPLOYMENT_GUIDE.md`, `RUNBOOK.md` | Run It Yourself / quickstart |
| `HCLS_AI_FACTORY_PROJECT_BIBLE.md` | About / Mission + canonical framing |
| `brief/` (capability-brief.html, mission.html, architecture.svg) | Keep as the Capability Brief |

## Appendix B — Immediately buildable now (no GPUs/RunPod, no new credentials)
Design system · OG card + favicon · `llms.txt` + sitemap + JSON-LD · the maturity matrix from the
registry · the IA skeleton + stub pages · engine/agent/TSC/platform pages from latent content ·
citations page · the registry-sync CI guard. → **This is all of Phase 0 and most of Phase 1.**
