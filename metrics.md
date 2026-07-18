# Site Metrics — governance record

A tracked, weekly-updated record of the metrics that actually compound for hcls-ai-factory.org — **not
a growth dashboard.** Reach metrics (pageviews, session duration, bounce) are deliberately excluded:
they measure attention, which is not the job. Per the Website Upgrade PRD (§A.3), this file holds the
CSV trends and the month's AI-summary transcripts, updated at the weekly reconciliation.

## What we measure (and why)

| ID | Signal | How it's captured | Target |
|---|---|---|---|
| **SM-1** | **Builders acted** — stars, forks, unique clones; fork:star ratio | `gh api` + GitHub Insights → Traffic, weekly to the CSV below | trend up vs. baseline (absolute TBD — OQ-1) |
| **SM-2** | **Experts converted** — inbound from clinicians/geneticists/oncologists/pharma; citations | tagged intake + CITATION.cff / Scholar backlinks | ≥1 qualified inbound/active week (TBD — OQ-2) |
| **SM-3** | **Depth-page engagement** — engine/agent/TSC/governance pages read, not just hero | cookieless analytics funnel (tool TBD — OQ-5) | depth-share above hero-only bounce (TBD — OQ-3) |
| **SM-4** | **Agent discoverability** — `llms.txt` hits; how assistants summarize the project back | server-log grep + monthly AI-summary audit | served + fetched; monthly summaries pass with zero overclaim |
| **SM-5** | **The "no overclaim" bar** (primary, qualitative) | pre-launch red-team + open "found an overclaim?" issue channel | **zero** open overclaim findings at each release gate |

## Weekly trend (SM-1) — fill at each reconciliation

| Week (ISO) | Stars | Forks | Unique clones (14d) | fork:star | Notes |
|---|---|---|---|---|---|
| _baseline (DEP-6)_ | TBD | TBD | TBD | TBD | capture before launch |

## Monthly AI-summary audit (SM-4)

Prompt ChatGPT / Claude / Gemini / Perplexity: *"What is the HCLS AI Factory?"* and score the answer
against the [honesty ledger](docs/honesty/ledger.md). Record the transcript + any overclaim found.

| Month | Assistant | Overclaim? | Transcript link / note |
|---|---|---|---|
| _pending first run_ | — | — | — |

## Overclaim findings (SM-5)

| Date | Source (red-team / issue) | Finding | Status / time-to-fix |
|---|---|---|---|
| _none yet_ | — | — | — |

*This is a governance record. It is updated by hand at the weekly reconciliation, not by a SaaS cockpit.*
