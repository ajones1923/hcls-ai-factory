# What to Do Before the Gated Software Arrives

**Date:** 2026-08-16 · Written after the full audit, with the explicit constraint that
**Parabricks, the BioNeMo NIMs, Chai-1, CUDA PyTorch and ESMFold are not yet installed.**

Everything below is achievable **without any gated component.**

---

## The strategic point

The platform's weakest measured dimension is not code quality or accuracy — it is that
**none of it is running.** 8,402 tests pass; zero engine or agent containers exist on the box.

That matters most *right now*, for one reason:

> When the gated software lands, you want to slot it into a system that already works. Otherwise
> you debug the NIM integration and the platform bring-up at the same time, and you cannot tell
> which layer is failing.

Bring it up first. Everything else on this list is secondary to that.

---

## Tier 1 — do these before touching gated software

### 1. Bring the platform up and keep it up

Nothing here is blocked. Milvus is already running with 12 collections, 673 GB is free, and all 17
subjects are declared in compose.

```bash
cp .env.example .env && ${EDITOR:-nano} .env    # 4 vars; compose will not parse without them
docker compose -f docker-compose.dgx-spark.yml build
docker compose -f docker-compose.dgx-spark.yml up -d
.venv/bin/python scripts/run_demo.py --check-all
```

**Disable the cron supervisor first** (`crontab -l | grep -v health-monitor.sh | crontab -`) or it
will restart services underneath you every 5 minutes.

**What this unlocks immediately:** 11 of 17 demos become genuinely LIVE — every agent, plus
precision-intelligence, oncology and cardiology. None of them need a gated component.

### 2. Seed the agent corpora

Milvus currently holds **12 collections, all `imaging_*`**. Every other agent would answer from an
empty index and look broken while working correctly.

**12 subjects ship seed/ingest scripts** — this is a supported path:

```
cart (13)  ·  precision-oncology (15)  ·  clinical-imaging (14)  ·  clinical-trial (2)
neurology (2)  ·  pharmacogenomics (2)  ·  precision-biomarker (2)  ·  single-cell (2)
rare-disease (2)  ·  cardiology (2)  ·  precision-autoimmune (1)  ·  precision-intelligence (1)
```

An agent demo without a seeded corpus is not a demo.

### 3. Turn the security gate on — one environment variable

Authentication is now implemented for all 12 endpoints but **ships off by default**:

```bash
echo "HCLS_API_KEY=$(openssl rand -hex 24)" >> .env
```

Verified behaviour: 401 without a key, 401 with a wrong one, `/health` still public. This is the
single highest value-to-effort change available. Security 7.5 → ~9.

### 4. Verify every `live` capability actually answers

Two capabilities were found registered `live` while nothing bound their ports. Once the platform is
up, close the loop permanently — `BUILD_GUIDE.md` §7 has the probe script. Consider adding it to CI
as a nightly job against a running stack.

---

## Tier 2 — remove the latent hazards

### 5. Rename the 12 stdlib-shadowing modules

`src/collections.py` in 11 subjects (plus `precision-autoimmune/config/logging.py`). Imports are
package-qualified today, so it works — but any config that puts `src/` on `sys.path` kills the
interpreter *before collection*:

```
ImportError: cannot import name 'namedtuple' from partially initialized module 'collections'
```

This produced two false readings during the audit and it will bite a Dockerfile eventually. Rename
to `vector_collections.py`. Mechanical, and there are 8,402 tests to verify against.

### 6. Choose one runtime model

`docker-compose.dgx-spark.yml` versus `health-monitor.sh` under cron. Two half-maintained models is
exactly why the port map drifted on 8 of 13 subjects. Per-service venvs exist for only 4 of ~20
services. **Recommend standardising on compose** and reducing `health-monitor.sh` to container
health checks.

### 7. Make the governance gates actually run

They are installed on all 12 services but **handlers still do not call them**. `install_governance`
is middleware; `require_valid_input()` and `honesty_flags()` must be invoked inside the endpoint.
Until they are, `X-HCLS-Governed` will be honestly absent on nearly every request.

Start with the three that emit the most consequential output: pharmacogenomics, precision-oncology,
rare-disease.

---

## Tier 3 — credibility polish (cheap, visible)

| # | Item | Effort |
|---|---|---|
| 8 | Add the decision-support frame to **genomic-foundation** and **precision-intelligence** — 5 of 17 subjects never state it, and these two reach clinical interpretation | low |
| 9 | Label **CYP3A4** → midazolam/cyclosporine as substrate relationships, not CPIC guideline pairs. Accuracy 9.5 → 10 | trivial |
| 10 | Set the repo **homepage** (`https://hcls-ai-factory.org`) and **topics** — both empty on a public Apache-2.0 repo | trivial |
| 11 | Refresh `docs/infographic_prompts/` to the 8·8·1 roster (says "6 agents"/"7 agents"). Not published, but it is in a public repo | low |
| 12 | Get the `subjects` CI job green and drop `continue-on-error` | medium |
| 13 | Agree a **test-depth floor** for clinical-output subjects. Range is 4 → 1,966 tests | medium |

---

## What NOT to do yet

- **Do not promote any capability to `live` in anticipation.** Status follows evidence. That rule
  is what caught the two unreachable capabilities.
- **Do not relabel demos E1/E3/E7** from REPRESENTATIVE. They stay that way until the NIMs answer.
- **Do not adopt Git LFS.** Evaluated and rejected — it trades repo size for per-build bandwidth on
  a git-linked Netlify site, and adds a pointer-file failure mode.
- **Do not install docs dependencies into the platform venv.** `requirements-docs.txt` drags in
  `zarr` and `fast-array-utils`, whose pytest plugins abort collection for three subjects. Use
  `.venv-docs`.

---

## Suggested order

```
.env + HCLS_API_KEY  ──►  compose up  ──►  seed corpora  ──►  verify live capabilities
                                                │
                                                ├──►  11 demos genuinely LIVE
                                                └──►  THEN start gated software
```

By the time Parabricks and the NIMs arrive, you want a running, authenticated, seeded platform with
11 demos already proven. Then each gated component is a single, isolated change you can verify on
its own.
