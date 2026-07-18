# hcls-ai-factory.org — Visual Design Plan

*How to make the site as visually stunning as possible — without spending a single unit of that
polish on hype. For this project, stunning and credible are the same axis.*

Status: **plan** · Companion to `SITE_UPGRADE_PLAN.md` (IA/content) and the Website Excellence
Research Report. Scope: **visual craft** — the layer the current site has barely touched.

---

## 0. The thesis (what "stunning" means here)

The current site is honest, fast, and complete — but it *looks like documentation*: stock Material
chrome (teal bar, sidebar, right-hand ToC), the mission line, a paragraph, an inline diagram. It
reads. It does not **land**.

The best technical/scientific product sites — **Linear, Vercel, Stripe, DeepMind AlphaFold** — are
stunning through **restraint and craft**, never decoration: one considered ground, one bold
typographic statement, ruthless removal of clutter, and a single signature visual. Crucially, for a
**precision-medicine, honesty-first** project, spectacle that reads as hype *destroys* the exact
trust the project sells (the domain-expert and skeptic personas punish it publicly). So the design
rule is unusual and freeing:

> **Every gram of visual polish must make the site feel *more* trustworthy, not more exciting.**
> The signature "wow" is the honesty itself — the live **Capability Maturity Matrix going green** —
> rendered as a beautiful, data-dashboard-grade centerpiece. No competitor in medical AI leads with
> machine-verified honesty *as the hero*. That is stunning **and** unattackable.

## 1. The one structural move: two experiences, one site

MkDocs Material is a *docs* theme; forcing a landing page through the docs template is why the home
looks like a page. The fix is the standard best-in-class pattern (Material's own site, FastAPI, etc.):

- **Home = a designed landing page** that *breaks out* of the docs chrome via a full `overrides/home.html`
  template — no sidebar, no ToC, full-bleed hero, its own sections. This is where all the visual
  ambition goes.
- **Everything else = Material docs**, refined but familiar — the depth pages stay scannable and
  fast. Experts drill; they don't want a landing page on the E1 spec.

This is the two-cut principle expressed in *layout*: a landing cut and a reference cut.

## 2. Design system v2 — off the defaults

The current palette is **stock Material teal + deep-orange** — one of the recognizable
"AI-generated / un-chosen" tells. We already have the hard part (self-hosted **Fraunces** serif +
**IBM Plex Sans/Mono**, a real wordmark, the status-badge component). v2 adds a *chosen* palette,
a hero, motion, and depth.

### 2.1 Palette — "clinical nocturne"
A deep, quiet, luminous dark ground for the landing (the Linear/Vercel move, but clinical, not neon),
carried into an airy light mode for the docs. Named tokens (evolves the existing infographic family
so the whole identity stays coherent):

| Token | Hex | Use |
|---|---|---|
| `--ink` | `#0B1722` | landing ground (near-black, blue bias) |
| `--ink-2` | `#132B3E` | panels, the header, card grounds |
| `--paper` | `#F7F7F3` | docs/light ground (warm off-white, not pure white) |
| `--jade` | `#3E8E8A` | structural accent (evolved from teal — a touch desaturated, less "default") |
| `--signal` | `#E0A65E` | the one warm "hope" accent — used *rarely* (a single CTA glow, the "Hope" beat) |
| status set | green/green-check/amber/blue/slate/violet | the maturity badges (already built) — the site's signature color language |

The accent discipline: **one warm accent, used once per view.** Everything else is ink, jade, and
the semantic status colors. A pure-grey neutral reads as unconsidered — bias the greys slightly
toward the ink.

### 2.2 Typography — make the serif do the work
Fraunces is already loaded but under-used. Push it: an **oversized serif hero headline** (the mission
line) with tight tracking and `text-wrap: balance`; IBM Plex Sans for body at a comfortable 65-char
measure; IBM Plex Mono for every port/status/`code` token (data should *look* like data). A real type
scale, held everywhere. This pairing (characterful serif + technical sans) is deliberately **not** the
Inter/Space-Grotesk default the whole category uses.

### 2.3 The signature visual — the matrix as a living dashboard
The maturity matrix is currently an honest table. Make it the **hero centerpiece**: a compact,
animated-on-load status board (badges settling into place; a subtle "9 verified" counter) with
data-dashboard elegance — filterable, `tabular-nums`, a legend. It's the thing a skeptic screenshots
*in the project's favour*. This is where "amazing" and "honest" become one object.

### 2.4 Motion — with restraint (this is where hype creeps in)
- A single orchestrated **page-load reveal** on the hero (headline fades up; the matrix badges settle;
  the "going green" counter ticks once). Nothing loops.
- Scroll-reveal on section entry (subtle, 1 per section), hover micro-states on cards/CTAs.
- **Always** respect `prefers-reduced-motion` (kill all of it). No parallax, no ambient particles —
  those read as AI-startup, the opposite of clinical.

### 2.5 Components
Elevate the repeating pieces to a component language: status badge (built) · **capability card**
(icon, name, one-line, badge, hover-lift) · **hero** · section eyebrow labels (mono, letter-spaced) ·
**stat block** (e.g. "9 verified · 25 live · in hours") · citation footnote · provenance tag on every
figure. Glassmorphic panels on the dark hero (used once), sharp cards in the docs.

## 3. The homepage, section by section (the landing)

A single scroll, each section a deliberate beat (weight → compression → hope):

1. **Hero** (dark, full-bleed): wordmark nav · the **mission line as an oversized serif headline** ·
   one calm subhead ("Patient DNA → therapeutic candidates, in hours — open-source, one box") · two
   CTAs (**Explore the factory** / **Run it yourself**, the repo link above the fold) · the **live
   maturity board** as the right-hand signature. One bounding honesty clause, visible.
2. **The pipeline** — the "in hours, not months" DNA→drug arc (the existing infographic, re-cut clean),
   with a stat block.
3. **The factory** — 8 Engines · 8 Agents · 1 Program · One Platform, as four elegant capability-card
   rows (registry-sourced), each linking into the docs.
4. **Honest by construction** — the matrix + the "9 verified on real data" proof, the differentiator
   stated plainly with the badge legend.
5. **The flagship** — a dignified TSC beat (weight→compression→hope), understated.
6. **Run it / fork it** — the builder CTA (one-command bring-up, Apache-2.0 commercial-OK).
7. **Footer** — mission line repeat, license, GitHub, the honest one-liner.

Everything in it is registry-sourced or cited; the honesty invariants (INV-1…9) apply at full force —
the landing is exactly where an expert looks for the crack.

## 4. Depth-page polish (keep Material, refine it)
- A slightly custom Material palette (the v2 tokens), refined card CSS, the mono-for-data rule,
  section eyebrows, provenance tags on figures, the badge component everywhere.
- The **matrix page** becomes a proper dashboard (filters, counts, legend) — not just on the home.
- Engine/agent pages: a small header band with the badge + ports as a clean spec, the six-block spine.

## 5. Technical implementation (honest + reproducible — non-negotiable)
- `overrides/home.html` (extends `main.html`, `{% block content %}`), `overrides/partials/*` where
  needed; disable `navigation.instant` interplay on home or scope styles to a `.home` body class.
- All visual work in `extra_css` (v2 tokens + hero + components) + tiny vanilla JS for the load reveal
  and the matrix animation (no framework, no CDN — **self-hosted, offline-safe, CON-2**).
- **Performance budget** (a conversion + credibility concern — a 1s delay ≈ −7% conversions): inline
  critical CSS, lazy-load below-the-fold images, keep the woff2 subset, ship the hero as static + a
  progressive-enhancement reveal. Lighthouse perf/a11y/SEO all green.
- **Accessibility**: WCAG-minded contrast on the dark hero, color-plus-text status (already true),
  visible focus, reduced-motion, alt text with provenance.
- Keep the strict build + drift-guard gating every change; nothing in the redesign hard-codes a
  status (all badges stay registry-sourced).

## 6. Phased delivery (each phase shippable, gated, reversible)
- **V1 — the hero (highest impact):** `overrides/home.html` + the dark hero + v2 palette tokens + the
  serif headline + the live matrix board + CTAs. This alone transforms the first impression.
- **V2 — the landing scroll:** sections 2–7, the capability-card rows, the stat blocks, scroll-reveal.
- **V3 — the dashboard matrix + depth polish:** the filterable matrix, refined docs cards, eyebrows,
  provenance tags.
- **V4 — motion + performance + launch pass:** the orchestrated load reveal, Lighthouse pass,
  OG-card refresh to match, the two-cut launch.

## 7. Reference exemplars (steal the *restraint*, not the theme)
**Linear** (typographic dark hero, zero clutter) · **Vercel** (dark gradient field, binary CTAs) ·
**Stripe** (craft + precision, credibility) · **DeepMind AlphaFold** (scientific gravitas, a database
as a public good) · **Material for MkDocs' own site** (the docs-with-a-landing pattern to copy
mechanically). The move each shares: **one ground, one headline, one signature visual, nothing else.**

## The bottom line
The site is already *honest and complete* — the rarest, hardest part. What it lacks is a **first
impression worthy of the work**. The plan spends all of its visual ambition on one thing: making the
honesty *beautiful* — a designed landing whose signature is a live, going-green integrity board.
Stunning, and impossible for a skeptic to attack, because it hides nothing.
