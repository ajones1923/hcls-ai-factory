---
title: Getting hcls-ai-factory.org to 10/10
description: An honest assessment of the public site and the specific gaps between where it is now and best-in-class, with the work ordered by what a visiting clinician actually notices.
---

# Getting hcls-ai-factory.org to 10/10

**Date:** 2026-08-16 · Assessed against the live site, not the source.

The site is the main focal point for all the material, so it is scored as a *destination for a
visiting clinician or researcher* — not as a docs build.

---

## Current: **8 / 10**

| Dimension | Score | Basis |
|---|---:|---|
| Factual accuracy | **9.8** | 0 port mismatches · 0 status contradictions · 44 CPIC pairs verified · 0 gene-locus errors |
| Honesty framing | **9.5** | Public mission wording correct on all 4 pages; 8–12 decision-support mentions per clinical page |
| Depth of material | **9** | 69 subject pages · 11 audit documents · maturity matrix generated from the registry |
| Media | **9** | 22 videos, **all now in Adam's clone**; 17 new infographics |
| Visual consistency | 8 | one illustration family; a few older figures predate it |
| **Reproducibility** | **4** | **the weakest link — see below** |
| Discoverability | 5 | repo homepage unset, no topics, no sitemap submission |
| Scholarly credibility | 5 | no citable paper |

---

## The one that costs the most: reproducibility

"Run It Yourself" is **146 lines across three pages** and tells a visitor to `git clone`. It never
mentions:

- the **`.venv`** the platform needs (`hcls_common` will not import without it)
- **`HCLS_API_KEY`** — the security gate
- **`scripts/run_all_tests.py`** — the one command that proves the platform works
- **`scripts/run_demo.py`** — the 17 demos

A researcher following that page today hits exactly the wall documented in `BUILD_GUIDE.md` §1:
a stale editable install and `No module named 'hcls_common'`. They conclude the project does not
work and leave.

**This is the single highest-value fix on the site.** The content already exists — `BUILD_GUIDE.md`
is a tested, cold-clone-to-running path. The run section should be rewritten from it.

> Claiming "fork it and run it" while the documented path fails is an honesty problem, not just a
> docs problem. It is the same class as a capability marked `live` that nothing serves.

---

## What 10/10 requires

### Tier 1 — credibility with the audience you actually have

| # | Gap | Why it matters |
|---|---|---|
| **1** | **Rewrite `/run/` from `BUILD_GUIDE.md`** | The reproducibility claim is the project's differentiator. Today it does not survive contact. |
| **2** | **Publish a citable paper** (batch ask 11) | Clinicians and researchers cite. Without a DOI or arXiv ID the work cannot enter the literature. Add `CITATION.cff` linkage. |
| **3** | **Make one demo runnable end to end by a visitor** | E8 single-cell already runs on CPU with no gated software — 2,700 cells → 9 clusters. Publish its transcript and the exact command. |
| **4** | **Set repo homepage + topics** | Both empty on a public Apache-2.0 repo. Free discoverability. |

### Tier 2 — completeness

| # | Gap |
|---|---|
| 5 | Add the decision-support frame to **genomic-foundation** and **precision-intelligence** — 5 of 17 subjects never state it, and both reach clinical interpretation |
| ~~6~~ | ~~Label CYP3A4 pairs as substrate, not CPIC~~ — **done**, and the real error was different: the KB claimed CPIC references CYP3A4 for tacrolimus dosing, when that guideline is CYP3A5-based |
| 7 | Surface the **honesty ledger** more prominently — it is the project's strongest differentiator and currently sits three clicks deep |
| 8 | Add a **"what is not live yet"** page. The maturity matrix carries this, but a plain-language version is more persuasive than a table |

### Tier 3 — polish

| # | Gap |
|---|---|
| 9 | Refresh `docs/infographic_prompts/` to the 8·8·1 roster (says "6 agents"/"7 agents"). Unpublished, but in a public repo |
| 10 | Consistent poster frames for all 22 videos |
| 11 | Submit the sitemap to Google/Bing Search Console |
| 12 | Add `og:` / Twitter card metadata for link previews |

---

## What is already best-in-class — do not disturb

- **The honesty framework.** LIVE / REPRESENTATIVE / BURST, the maturity matrix generated from the
  registry rather than hand-written, and the decision-support frame stated inline rather than in a
  footer. Very few projects in this space do this at all.
- **Pharmacogenomics content.** 44 CPIC pairs across 13 genes, every one verified genuine, with
  correct variant identifiers (VKORC1 rs9923231, SLCO1B1 rs4149056).
- **The public mission wording.** *"No one should wait years for a disease we could have understood
  in a day."* Softer than the internal line and correct for a family-facing audience. It is
  consistent across all 4 pages that use it.
- **The generated subject guides.** Ports and test counts come from `scripts/subject_facts.py`, so
  they cannot silently drift from the code.

---

## Ordered plan

```
1. Rewrite /run/ from BUILD_GUIDE.md      ← biggest credibility gain, content already exists
2. Publish E8's transcript as a live proof point
3. Repo homepage + topics
4. Decision-support frame on the 2 engines that lack it
5. CYP3A4 labelling
6. The paper (ask 11)
```

Items 1–5 are a day's work and move the site from 8 to ~9.3. Item 6 is what makes it a 10, and it is
the one that cannot be rushed.
