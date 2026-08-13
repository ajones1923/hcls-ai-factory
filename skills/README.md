# HCLS AI Factory — Skills

This folder holds the **skills** that keep the factory's build aligned with its foundation and
strategy: durable, versioned guidance a developer (or an AI assistant) references while building
the next generation of the platform. Two kinds of thing live here:

- **Invocable skills** — a directory with a `SKILL.md` file (YAML frontmatter + instructions)
  that can be loaded and referenced by name.
- **Reference notes** — plain text/markdown that captures strategy or standing context.

## Layout

```
skills/
├── README.md                          # this file
├── hcls-core-vision-mission/          # ⭐ NORTH STAR (WHY) — the mission; consult first
│   └── SKILL.md                       #   the founding mission, vision, and principles
├── engines-agents-disease-programs/        # ⭐ NORTH STAR (WHAT) — the canonical capability roster
│   └── SKILL.md                       #   8 engines · 8 agents · 1 disease program (TSC), with ports
├── AI Factory Core Pillars.txt        # the 16 foundational pillars (strategic reference)
├── build-housekeeping-standards/      # invocable skill — HOW we build
│   └── SKILL.md                       #   best-practices checklist for every build
├── demo-foundation-alignment/         # invocable skill — WHAT we build (to the demo portfolio)
│   └── SKILL.md                       #   anchors development to the D1–D7 demos + coverage matrix
├── development/                       # building & shipping capabilities
│   ├── capability-delivery-playbook/SKILL.md   # end-to-end "ship a capability" checklist
│   ├── disease-program-authoring/SKILL.md      # add a disease-program vertical (TSC template)
│   └── model-integration/SKILL.md              # add a model / NIM / frontier model cleanly
├── governance/                        # honesty & compliance
│   ├── clinical-claim-honesty/SKILL.md         # verify + label every clinical claim before it ships
│   └── regulatory-compliance-posture/SKILL.md  # not-a-device, Part 11, HIPAA/GDPR, synthetic-only
├── clinical/                          # clinical-domain grounding (keep claims correct)
│   ├── clinical-genomics-standards/SKILL.md     # ACMG/ClinVar/AlphaMissense, mosaicism, ACMG SF
│   ├── pharmacogenomics-cpic/SKILL.md           # star alleles → CPIC dosing / safety interlocks
│   └── oncology-mtb/SKILL.md                    # MTB packets, fusion-first pediatric, trial match
├── communication/                     # telling the story
│   ├── demo-script-authoring/SKILL.md           # write a demo to weight→compression→hope + two-cut
│   └── release-and-site-publishing/SKILL.md     # cut a release → Netlify → hcls-ai-factory.org
├── branding/                          # naming + messaging conventions
│   └── naming-and-messaging/SKILL.md            # how to name engines/agents/programs + on-brand copy
├── personas/                          # audience-awareness / messaging
│   ├── broad-general-persona/SKILL.md           # the 5 personas + guardrails for any outreach
│   ├── persona-ai-curious-public/SKILL.md
│   ├── persona-domain-experts/SKILL.md
│   ├── persona-builders-oss/SKILL.md
│   ├── persona-patients-families/SKILL.md
│   └── persona-skeptics/SKILL.md
├── assistant/                          # tooling the assistant drives (audio, film, media)
│   └── eleven-labs/
│       ├── elevenlabs-video-production/SKILL.md  # PLAN: the 18-film engine/agent series
│       ├── video-demos-alignment/SKILL.md        # bind films to the demos doc + honesty ledger
│       └── hcls-film-series-pipeline/SKILL.md    # ⭐ BUILT: the shipped 7-film series — read
│                                                 #   before changing anything in tmp/videos/
└── architecture/                      # one invocable skill per Core Pillar (16)
    ├── 01-compute-dgx-spark-remote-gpus/SKILL.md
    ├── 02-storage-and-data-layer/SKILL.md
    ├── 03-networking-and-ingress/SKILL.md
    ├── 04-containers-and-orchestration-runtime/SKILL.md
    ├── 05-structured-databases/SKILL.md
    ├── 06-vector-databases-and-embeddings/SKILL.md
    ├── 07-message-bus-and-async/SKILL.md
    ├── 08-inference-serving/SKILL.md
    ├── 09-ai-orchestration/SKILL.md
    ├── 10-observability/SKILL.md
    ├── 11-security-and-secrets/SKILL.md
    ├── 12-cost-and-economics/SKILL.md
    ├── 13-reliability-and-operations/SKILL.md
    ├── 14-ease-of-deployment/SKILL.md
    ├── 15-github-structure-and-presentation/SKILL.md
    └── 16-github-netlify-site/SKILL.md
```

### `hcls-core-vision-mission/SKILL.md` ⭐ (north star)
The founding **mission, vision, and principles** — *"No one should die of a disease we could have
understood in time."* This is the guiding context for **all** work: it governs what to build, how
to frame it, and what to promise, and every other skill exists to serve it. **Consult it first, on
every task**; when a technical choice and the mission seem to conflict, the mission wins.

### `engines-agents-disease-programs/SKILL.md` ⭐ (north star — what the factory *is*)
The canonical capability roster — the authoritative, current picture of the **8 engines, 8
intelligence agents, and the 1 disease program (TSC)** with ports and what each does, plus the
platform layer, frontier models, and the honesty ledger. Consult it whenever naming or counting
engines/agents or creating any artifact (doc, figure, demo, video, README, deck) that must reflect
the canonical framing — so nothing reverts to a stale count or invents a capability. Distilled from
the "8 Engines - 8 Agents Focus" paper; pairs with `hcls-core-vision-mission` (why) and
`demo-foundation-alignment` (how it's proven).

### `AI Factory Core Pillars.txt`
The 16 pillars the platform is built on — compute, storage, networking, orchestration,
databases, embeddings, messaging, inference serving, observability, security, cost, reliability,
ease of deployment, and GitHub/site presentation. This is the strategic "what the factory is made
of" reference; skills operationalize the housekeeping side of these pillars.

### `build-housekeeping-standards/SKILL.md`
The standing checklist for **how** we build: the pre-commit gate, golden paths for adding an
engine or agent, capability registration (the registry drift-guard), governance wiring, the
dependency baseline, neutrality/secret/path/data hygiene, docs-stay-true, and verify-then-commit
discipline. Reference it whenever adding or modifying an engine, agent, capability, dependency,
test, or doc. It maps back to the Core Pillars (esp. Security, Reliability, Ease of Deployment,
and GitHub Structure).

### `demo-foundation-alignment/SKILL.md`
The standing rule for **what** we build: every new engine, agent, model, or capability must be
anchored to the **demonstration foundation** — the D1–D7 demo portfolio and its capability-coverage
matrix. It carries the portfolio structure, the shock-and-hope principle, the honesty ledger, and
a checklist so each capability earns a demo home and keeps the coverage matrix green. Consult it
when planning or reviewing any capability. It is the bridge between the north star (*why*) and
`build-housekeeping-standards` (*how*).

### `development/`
Building & shipping capabilities. `capability-delivery-playbook` is the end-to-end "ship a
capability" checklist that ties the whole framework together (mission → demo home → build →
register → verify → honesty → docs → gate). `disease-program-authoring` covers adding a new
disease-program vertical (TSC as template). `model-integration` covers adding a model / NIM /
frontier model cleanly (serving mode, license gate, honest status, verification).

### `governance/`
Honesty & compliance. `clinical-claim-honesty` is the active discipline for verifying and labeling
every clinical claim before it ships anywhere (the honesty ledger). `regulatory-compliance-posture`
holds the conservative posture: not a medical device, decision-support-not-diagnosis, 21 CFR Part
11 reproducibility, HIPAA/GDPR, synthetic-and-public-data-only.

### `clinical/`
Clinical-domain grounding that keeps the factory's claims correct: `clinical-genomics-standards`
(ACMG/ClinVar/AlphaMissense, mosaicism, ACMG SF), `pharmacogenomics-cpic` (star alleles → CPIC
dosing + safety interlocks), and `oncology-mtb` (MTB packets, fusion-first pediatric, trial match).
Consult the relevant one when building or reviewing that clinical area.

### `communication/`
Telling the story. `demo-script-authoring` (write a demo to the weight→compression→hope arc + the
two-cut principle) and `release-and-site-publishing` (cut a release → Netlify → hcls-ai-factory.org).

### `branding/`
Naming and messaging conventions. `naming-and-messaging` defines how the factory names engines (by
function + "Engine"), agents (by domain + "Intelligence Agent"), and disease programs, why
"Intelligence" is reserved for the reasoning agents, how to resolve an engine/agent name collision
(the Analysis → Compute → drop-qualifier order, e.g. Single-Cell Analysis Engine vs. Single-Cell
Intelligence Agent), and the on-brand messaging rules. Consult when naming anything or writing
factory-facing copy.

### `personas/`
Audience-awareness and messaging skills for outreach. `broad-general-persona` carries all five
personas + guardrails — consult it for **any** artifact people outside the codebase will read or
watch. Five per-persona drill-downs go deeper: `persona-ai-curious-public`, `persona-domain-experts`,
`persona-builders-oss`, `persona-patients-families`, `persona-skeptics`.

### `architecture/`
One invocable skill per **Core Pillar** — 16 in total, numbered `01`–`16` to match
`AI Factory Core Pillars.txt`. Each is a best-practice standard for that pillar (compute, storage,
networking, containers, structured databases, vector databases + embeddings, messaging/async,
inference serving, AI orchestration, observability, security, cost, reliability, ease of
deployment, GitHub structure, and the Netlify site), grounded in the real factory stack. Reference
the relevant pillar skill when designing, building, operating, or reviewing that layer;
`build-housekeeping-standards` is the cross-cutting "how we build" skill that ties them together.

## Skill convention

Each invocable skill is its **own directory** containing a `SKILL.md`:

```markdown
---
name: my-skill-name          # kebab-case, matches the directory name
description: >-              # one clear paragraph on WHEN to use it (drives discovery)
  What this skill is for and the situations that should trigger it.
---

# Title

The instructions / checklist body.
```

Keep skills **neutral** (no proprietary vendor or alternate-edition branding — the same rule the
repo's pre-commit guard enforces) so they commit cleanly and stay publishable. Reference real
files by path; keep them actionable and scannable.

## Where skills live and how they're discovered

- **In the repo (`skills/`)** — versioned, reviewable, and shareable. This is the canonical copy
  and what ships with the open-source project.
- **In `.claude/skills/`** — how Claude Code auto-discovers invocable skills in a session.
  `.claude/` is gitignored (local only), so hooking a skill up for auto-load never publishes.

To make the repo skills auto-loadable, sync them into `.claude/skills/` as real copies (chosen
over symlinks for maximum discovery compatibility across Claude Code builds):

```bash
./scripts/sync-skills.sh      # rebuilds .claude/skills/ from skills/ (run after editing a skill)
```

The repo `skills/` tree is the single source of truth; `sync-skills.sh` re-copies it, so re-run
it whenever you add or edit a skill. Start a fresh Claude Code session for newly synced skills to
be picked up. Every skill in this tree is synced this way — the two north stars, the cross-cutting
`build-housekeeping-standards` and `demo-foundation-alignment`, the development / governance /
clinical / communication / branding / persona skills, and the 16 pillar skills.

## Adding a new skill

1. Create `skills/<skill-name>/SKILL.md` with the frontmatter + body above.
2. Verify it's neutral and passes the repo's checks (`./scripts/install-hooks.sh` then a trial
   `git add` + commit runs the guard).
3. Run `./scripts/sync-skills.sh` to copy it into `.claude/skills/` for in-session auto-load.
4. Add a one-line entry to the Layout section of this README.
