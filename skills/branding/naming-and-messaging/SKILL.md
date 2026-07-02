---
name: naming-and-messaging
description: >-
  The HCLS AI Factory's naming and messaging conventions — how to name any engine, agent, or disease
  program, and how to talk about the factory so names stay consistent and copy stays honest and
  on-brand. Consult when naming a new capability, resolving a name collision (an engine and an agent
  in the same domain), or writing any factory-facing copy (README, deck, figure, post, docs). Keeps
  "Intelligence" meaningful, the counts canonical (8·8·1), and the claims honest.
---

# Naming & Messaging — the factory's brand rules

Headline (fixed): **Eight Engines · Eight Intelligence Agents · One Platform**, plus the disease
programs it powers. Keep every name and every line true to that and to the mission. Fuller
reference: the "HCLS AI Factory - Branding and Messaging" document.

## The naming model — three tiers, one rule each
| Tier | Named by | Pattern | Example |
|---|---|---|---|
| **Engine** | its **function** | `<Function> Engine` | Genomic Foundation Engine · Therapeutic Discovery Engine · Single-Cell Analysis Engine |
| **Intelligence Agent** | its **domain** | `<Domain> Intelligence Agent` (or `<Domain> Agent`) | Neurology Intelligence Agent · Pharmacogenomics Agent |
| **Disease Program** | its **condition** | a vertical | Tuberous Sclerosis Complex |

- **Engines are the machinery** — they *compute* and produce real artifacts. Name them for **what
  they do** (Foundation, Discovery, Imaging, Analysis…), **not** for their layer. Engines don't say
  "Compute" in their names.
- **Agents are the reasoning** — **"Intelligence" is the agent word.** It marks the
  reasoning/decision-support layer.
- **Disease programs** are named for the condition.

## Why "Intelligence" is reserved for agents
The two words are the two layers of one pipeline: `raw data → compute (Engine) → intelligence (Agent)
→ decision`. An engine computes; an Intelligence agent reasons over the result. Reserving
"Intelligence" for agents is what makes it mean something — so **don't sprinkle "Intelligence" onto
engines** for polish. (Same shape as a database engine vs. the analytics layer, or a compiler vs. an
IDE.)

## Resolving a name collision (an engine AND an agent in one domain)
When a domain has both — the clearest case is **single-cell** — a flat list makes the shared name
look like a duplicate. **Fix it by presentation, not architecture:** give the *engine* a functional
qualifier and let the *agent* keep "Intelligence." Order of preference for the engine's qualifier:

1. **Functional descriptor** (e.g. "Analysis") — *brand-friendliest*; reads like the other
   function-named engines, contrasts with "Intelligence" without a cold layer-word.
2. **"Compute"** — *crispest/most literal*, but a layer-word that sticks out among function-named engines.
3. **Drop the qualifier** — rely on "Engine" vs. "Intelligence Agent" with a light
   `(computes)`/`(reasons)` tag in shared lists.

**In diagrams:** render the pair as **one stacked capability** (Intelligence Agent cap on an Engine
base, "computes → interprets" link) — never two identical side-by-side nodes.

**Worked example — single-cell (the standard):**
- **Single-Cell Analysis Engine** (`singlecell-compute`, :8573) — scanpy analysis → cell-type clusters.
- **Single-Cell Intelligence Agent** (`single-cell-intelligence-agent`, :8130) — reasons over them.
- Why *Analysis*: it's literally what the engine does, reads like the other engines, and contrasts
  with *Intelligence* without the cold "Compute."

## Messaging principles
- **Honesty is load-bearing** — decision support, not diagnosis; label preclinical/roadmap/research-use;
  a `live` capability is never mock-served. Overstatement betrays the mission.
- **Openness is the hook** — not "AI did this," but *"AI did this — and it's yours"* (Apache-2.0,
  commercial-OK).
- **Lead with the canonical framing** — 8 engines · 8 agents · 1 platform + the disease programs;
  never a stale count (e.g. "six engines").
- **Match the audience; accuracy wins over spectacle** (see `broad-general-persona`).
- **Use the exact canonical names** — verify against the capability roster / registry; consistent
  names, palette, and voice read as a product.

## Do / Don't
**Do:** name engines by function + "Engine"; give agents "Intelligence Agent"; qualify an engine
(prefer "Analysis") only when it collides with a same-domain agent; stack the pair in diagrams; keep
counts at 8·8·1; keep the honesty labels in copy.
**Don't:** put "Intelligence" on an engine for polish; use a stale count; invent a name not in the
canon; let one artifact's copy exceed the honesty ledger; render a colliding engine/agent as two
identical nodes.

## Related
- `engines-agents-disease-programs` — the canonical roster these rules name.
- `hcls-core-vision-mission` — the mission the messaging serves.
- `broad-general-persona` · `clinical-claim-honesty` — audience + honesty discipline in copy.
- `figurelabs-illustrations` — applies the "stack single-cell / current-framing" rules in figures.
