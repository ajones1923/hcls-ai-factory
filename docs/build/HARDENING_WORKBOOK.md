# Workbook — Hardening the HCLS AI Factory

Companion to [HARDENING_PRD.md](HARDENING_PRD.md). Work it in order; each step carries the command
and the check that settles it. Where a step has already cost time on this machine, the trap is
written next to it.

```bash
cd ~/projects/hcls-ai-factory
export PY=.venv/bin/python
```

---

## Phase 0 — Close the door  (~1 hour)

### 0.1 See the problem yourself first

```bash
ip -4 addr show | grep -oP 'inet \K[\d.]+' | grep -v 127.0.0.1   # your LAN address
LAN=$(ip -4 addr show | grep -oP 'inet \K[\d.]+' | grep -v '^127\.' | head -1)

curl -s -m 6 -o /dev/null -w "health: %{http_code}\n" http://$LAN:8522/health
curl -s -m 20 -o /dev/null -w "clinical POST: %{http_code}\n" \
  -X POST -H 'Content-Type: application/json' -d '{"question":"test"}' http://$LAN:8522/search
```

Both return **200** today. Anyone on the network can ask a clinical agent a question, get generated
prose with dosing and management guidance, and spend your API key doing it.

### 0.2 Answer H-D1 before changing anything

Three options; pick one and write it in the PRD §7:

| Option | Do this | Cost |
|---|---|---|
| **Localhost only** | services bind `127.0.0.1` | portal unreachable from other machines |
| **Firewall the range** | deny 8500–8600 except from trusted hosts | demos still work on-box |
| **Keep exposure, require the key** | turn on `HCLS_API_KEY` | every client must send a header |

Exposed **and** unauthenticated is the one combination nobody would choose deliberately. Anything
else is fine if it is chosen.

### 0.3 Turn on the API key

```bash
sed -i 's/^# HCLS_API_KEY=/HCLS_API_KEY=/' .env
grep -c '^HCLS_API_KEY=' .env        # expect 1
./health-monitor.sh restart all      # services read it at startup
```

Verify fail-closed on a real endpoint:

```bash
curl -s -o /dev/null -w "no key:    %{http_code}\n" -X POST -H 'Content-Type: application/json' \
  -d '{"question":"x"}' http://localhost:8522/search
curl -s -o /dev/null -w "wrong key: %{http_code}\n" -H "X-API-Key: wrong" -X POST \
  -H 'Content-Type: application/json' -d '{"question":"x"}' http://localhost:8522/search
curl -s -o /dev/null -w "right key: %{http_code}\n" -H "X-API-Key: $(grep '^HCLS_API_KEY=' .env | cut -d= -f2)" \
  -X POST -H 'Content-Type: application/json' -d '{"question":"x"}' http://localhost:8522/search
```

Expect **401 / 401 / 200**.

> **Trap — this will break the demos and the eval.** `run_demo.py` and `run_clinical_eval.py` send
> no auth header. Turning the key on without updating them turns 17/17 into 0/17 and 6/6 into six
> ERRORs. Update both to read `HCLS_API_KEY` from the environment **in the same change**, then:
>
> ```bash
> for k in E4 A1 A6; do $PY scripts/run_demo.py $k >/dev/null && echo "$k ok"; done
> $PY scripts/run_clinical_eval.py | tail -2
> ```

### 0.4 Firewall the service range (if that was the decision)

```bash
sudo ufw status verbose                      # ufw reports active but is not blocking these ports
sudo ufw deny proto tcp from any to any port 8500:8600
sudo ufw reload
curl -s -m 6 -o /dev/null -w "from LAN: %{http_code}\n" http://$LAN:8522/health   # expect 000
curl -s -m 6 -o /dev/null -w "on box:  %{http_code}\n" http://localhost:8522/health # expect 200
```

Leave 80/443 alone — Caddy and the site ride on those.

### 0.5 Fix the license so GitHub can see it

```bash
gh repo view ajones1923/hcls-ai-factory --json licenseInfo -q '.licenseInfo.spdxId // "none"'
# -> none

curl -sL https://www.apache.org/licenses/LICENSE-2.0.txt -o /tmp/apache.txt
diff <(sed -n '1,190p' LICENSE) <(sed -n '1,190p' /tmp/apache.txt) | head -20
```

The body matches; the file is missing the **APPENDIX** and the copyright line placement that
`licensee` keys on. Replace with the canonical text, then re-apply the attribution:

```bash
cp LICENSE /tmp/LICENSE.bak
cp /tmp/apache.txt LICENSE
# then put "Copyright 2026 Adam Jones" into the APPENDIX placeholder at the end
git diff --stat LICENSE
```

Push, wait a minute, then:

```bash
gh repo view ajones1923/hcls-ai-factory --json licenseInfo -q '.licenseInfo.spdxId'   # Apache-2.0
```

### 🚦 Phase 0 gate

- LAN clinical POST refused (401 or 000)
- 17/17 demos and 6/6 eval still pass **with** the key on
- GitHub reports `Apache-2.0`

---

## Phase 1 — Supply chain  (~half a day)

### 1.1 Least-privilege CI tokens

Add to the top of `.github/workflows/ci.yml`, above `jobs:`:

```yaml
permissions:
  contents: read
```

No job in this workflow writes to the repo. Verify the next run still passes all six.

### 1.2 Pin actions by SHA

```bash
grep -ohE 'uses: [^ ]+' .github/workflows/*.yml | sort -u
gh api repos/actions/checkout/git/ref/tags/v4 -q .object.sha
gh api repos/actions/setup-python/git/ref/tags/v5 -q .object.sha
```

Replace `@v4` / `@v5` with the SHA, keeping the tag in a trailing comment.

### 1.3 Patch the platform venv

```bash
$PY -m pip install --upgrade pip accelerate
.venv/bin/pip-audit --local 2>/dev/null | tail -5      # expect clean
```

> Only **2** packages in the platform venv are affected. A bare `pip-audit` reports ~40 because it
> audits the OS Python too (`ufw`, `python-apt`, `louis`) — those are Ubuntu's, not the project's.
> Do not let that number drive a panic.

### 1.4 Triage the 14 dependabot PRs

```bash
gh pr list --state open --limit 20
```

Suggested handling:

| Group | Action |
|---|---|
| `github-actions` bumps | merge after 1.2 |
| minor/patch python group (#132) | merge, then run the gate |
| `anthropic >= 0.116` | merge — the SDK in use is already newer than the pin |
| `streamlit >= 1.58` | merge, then confirm the 4 UIs still serve |
| **`pymilvus 2.x → 3.0.0`** (#45, #49, #52, #54) | **hold — H-D3** |

> **Trap — pymilvus 3.0 is a major bump.** 2.6.8 is installed and working, and the agents'
> retrieval was repaired only days ago by fixing a client-API mismatch (`_CollectionManager`
> vs `MilvusClient`, `param=` vs `search_params=`, Hit-objects vs dicts). A second client-API
> change would look identical and be very hard to attribute. If you merge these, run
> `$PY scripts/run_clinical_eval.py` immediately after — it is the only detector.

After each merge:

```bash
ruff check --select E9,F82,F811,F706,F707 core lib scripts
( cd lib/hcls_common && ../../$PY -m pytest -q )
$PY scripts/run_all_tests.py | tail -2
$PY scripts/run_clinical_eval.py | tail -2
```

### 1.5 Make CI actually gate

```bash
gh api -X PUT repos/ajones1923/hcls-ai-factory/branches/main/protection \
  -f 'required_status_checks[strict]=true' \
  -f 'required_status_checks[contexts][]=Lint (real-bug rules)' \
  -f 'required_status_checks[contexts][]=Subject suites (8 engines, 8 agents, 1 program)' \
  -f 'required_status_checks[contexts][]=Capability registry validation' \
  -f 'required_status_checks[contexts][]=Docs site (strict build)' \
  -f 'required_pull_request_reviews[required_approving_review_count]=0' \
  -F 'enforce_admins=true' -F 'restrictions=null'
```

Verify by pushing a deliberately failing branch and confirming the merge button is blocked.

### 🚦 Phase 1 gate

`pip-audit --local` clean · fewer than 5 open PRs · a red run cannot reach `main`.

---

## Phase 2 — De-duplicate  ✅ done 2026-09-16

**Gate respected:** the clinical eval was expanded from 6 to 23 cases *before* any refactor, then
to 25. It is the only thing that can tell you a refactor changed an answer.

### 2.1 Extract the RAG path — done, and **smaller than this PRD claimed**

The PRD asserted "twelve `rag_engine.py` files, 13,821 LOC, 15–67% similar" and proposed
extracting the shared RAG path. Measured rather than assumed, that premise does not hold. Hashing
every function body with names, strings and comments normalised away:

| function | in files | distinct bodies |
|---|---|---|
| `query` | 12 | **12** |
| `retrieve` | 7 | 7 |
| `_rerank_results` | 5 | 5 |
| `_build_context` | 5 | 5 |
| `_search_collection` | 5 | **1** |
| `_save_conversation` / `_load_conversation` / `_cleanup_expired_conversations` | 5 | **1** |

The engines have genuinely **diverged**, they are not copies: twelve distinct `query` bodies in
twelve files. Forcing them behind one abstraction would invent coupling that the code does not
have, and would put the domain-specific weighting and prompts — the part that is supposed to
differ per agent — behind a shared seam.

What *was* duplicated is real and now shared:

- `hcls_common.vector_search.search_collection` — 88 lines × 5, identical in **every line of
  code**, differing only in one docstring example and one comment. It is the hottest path in the
  platform: every clinical answer passes through it once per collection searched.
- `hcls_common.conversation_store.ConversationStore` — the session-memory trio, byte-identical in
  five subjects and **tested in none of them**.

~520 lines removed; 24 tests added where there were none.

**The extraction also exposed a latent fault in the platform library.** Adding a `hcls_common`
import to an engine that had none made a cardiology test fail with:

```
ValueError: Duplicated timeseries in CollectorRegistry: {'hcls_milvus_search_seconds', …}
```

Every Prometheus collector in `hcls_common` was created at MODULE level inside a
`try: … except ImportError:` block, so a duplicate registration — a `ValueError`, not an
`ImportError` — killed the import of the module and every caller with it. Reaching it needs
nothing exotic: `mock.patch("src.rag_engine.X")` resolves its target by importing
`src.rag_engine`, and if the test already imported a bare `rag_engine`, Python holds two module
objects for the same file and runs its imports twice.

Seven library modules carried that fragility (19 registrations). `hcls_common.metrics.metric()`
now reuses a collector already registered under the same name, and re-raises any `ValueError`
that is *not* a re-registration so genuine misuse still fails. 8 tests.

A library module must not explode because it was imported twice — and this one would have, for
any subject that later imported `hcls_common` from a module reachable under two names.

**The honest conclusion is the finding.** "13,821 duplicated lines" was an estimate from file
sizes and surface similarity. The duplication that survives measurement is ~4% of that. Recording
the refutation is worth more than delivering the number the plan first promised.

### 2.2 Extract the LLM client — done (PR #140)

Eight per-service `_LLMClient` classes, byte-for-byte alike apart from which `settings.LLM_MODEL`
they read. That is why the `temperature` removal had to be applied in nine places and was missed
in eight of them — every agent's synthesis 400'd, each route caught it and fell back to a stub, so
the fleet reported healthy while answering "Search completed. See evidence passages below."

Now `hcls_common.service_llm`: −361 lines, 22 tests, one edit for the next model change. The only
genuine per-service difference — the default system prompt, which five of the eight baked in —
stays a parameter.

```bash
git grep -l 'class _LLMClient' -- '*/api/main.py' | wc -l    # 0
```

### 2.3 Stop shadowing the stdlib — done (PR #141)

Eleven subjects shipped `src/collections.py`; **nine files import both it and the real
`collections`**, so putting a subject's `src/` on `PYTHONPATH` killed the interpreter before
collection. The harness worked around it by withholding `src/` from nine of seventeen subjects —
meaning those nine suites ran against less code than CI claimed.

Renamed to `src/vector_collections.py`; imports, `mock.patch` string targets (these resolve by
string, so a missed one fails *silently*), and 108 prose references updated; the withholding guard
removed.

| | subjects | passed | failed | errors | `src/` withheld |
|---|---|---|---|---|---|
| before | 17 | 8397 | 0 | 0 | **9** |
| after | 17 | 8397 | 0 | 0 | **0** |

```bash
( cd core/agents/cart && PYTHONPATH="$PWD:$PWD/src" ../../../$PY -m pytest -q )   # 415 passed
```

### 🚦 Phase 2 gate — ✅ met

One LLM client · one search path · one conversation store · bare `pytest` works in every subject ·
eval unchanged (25/25, no answer altered by any of the three refactors).

The gate originally read "one RAG path". It is recorded as met on the measured finding above
rather than the assumed one: the twelve `query` implementations are not duplicates and were
deliberately left alone.

---

## Phase 3 — Polish

### 3.1 Get the video out of git

```bash
git ls-files -z 'docs/assets/videos/*.mp4' | xargs -0 du -ch | tail -1   # 192M of a 213M repo
```

Options: Git LFS, or host the files and reference them. Either way, history still carries them —
a `--depth 1` clone is the metric that matters:

```bash
git clone --depth 1 https://github.com/ajones1923/hcls-ai-factory /tmp/clonetest 2>&1 | tail -2
du -sh /tmp/clonetest && rm -rf /tmp/clonetest
```

### 3.2 Doc links and commands — ✅ done 2026-09-16

The raw link scan reports 40 hits; **29 resolve to build-time generated pages**
(`honesty/maturity-matrix.md`, `factory/engines/<id>.md`, `brief/README.md`) and are not broken.
`mkdocs build --strict` passing is the authority on those. Eleven were genuinely broken, pointing
at `licensing.md`, `DATA_SETUP.md`, `demo-guide.md` and friends — files that were never written.
They are repointed at the real documents, or de-linked where the prose already carries the fact.

**The larger finding was not the links.** Chasing `DATA_SETUP.md` turned up 26 references across
five guides instructing the reader to run a repository-root `setup-data.sh` with `--all` / `--stage1` /
`--stage2` / `--stage3` / `--status` — a script documented with a complete flag interface, called *"the
recommended approach"*, and **never written**. `start-services.sh` and `run_pipeline.py` were the
same. The documented first step of a public open-source repository failed immediately, and the
manual procedure that does work was labelled "for reference" underneath it.

Corrected to what exists: `core/engines/genomic-foundation/run.sh {check,login,download,reference,test,full}`
for Stage 1, the explicit manual procedure for Stages 2 and 3 (no automated downloader ships for
ClinVar or AlphaMissense — that is now stated rather than implied), and `./start-factory.sh` for
services.

`scripts/check_docs.py` now guards both, and runs in CI beside `mkdocs build --strict`:

```bash
.venv/bin/python scripts/check_docs.py          # docs/ + the genomics & intelligence engines
.venv/bin/python scripts/check_docs.py --all    # every tracked markdown file
```

It knows about the generated pages, follows `cd` inside a fenced block, and treats
NVFlare's provisioning-created `start.sh` / `fl_admin.sh` as legitimately absent — a check that
cries wolf is one nobody runs. It deliberately stops at *does a script by this name exist in the
repo at all*: proving the reader's working directory needs a shell interpreter, and the weaker
rule catches the whole defect class with no false alarms.

> A broken link is a nuisance. A confident instruction to run a command that does not exist is
> the honesty failure this platform exists to avoid — the documentation equivalent of a `live`
> capability that is mock-served.

### 3.3 Stale ports — ✅ done 2026-09-16

`docs/HCLS_AI_FACTORY_DEMO_GUIDE.md` and `docs/HCLS_AI_FACTORY_DGX_SPARK_DEPLOYMENT_GUIDE.md`
both printed pre-convention ports. The deployment guide's agent table had **eleven of eleven rows
wrong**, including the retired 81xx API block and rare-disease with its two ports transposed.

Corrected against `docs/build/PORT_MAP.md` and **verified against the running fleet** rather than
the registry alone — all eleven documented API ports answer `/health` with 200.

The failure mode here is worse than a stale number. The retired ports are mostly dead:

```
8103 → refused   8107 → refused   8128 → refused   8130 → refused   8134 → refused
8529 → 200
```

`8529` answers. The old guide listed it as **neurology**; it is **precision-biomarker's** API since
the 2026-08-15 sweep. A reader following that guide would not see a connection error — they would
get fluent, confident clinical answers from the wrong agent. A dead port fails loudly; a
reassigned one fails silently, which is the same class of fault as a `live` capability that is
mock-served.

Both guides now carry a banner naming `PORT_MAP.md` and the registry as the authority, so the
provenance is kept rather than quietly overwritten.

**Still stale, deliberately left:** the retired ports also appear in point-in-time documents —
`HCLS_AI_FACTORY_v1.3.0_RELEASE_REPORT.md`, `build/GAP_ANALYSIS.md`, `build/PRD.md`, the
architecture-research and infographic-prompt files (8107 in 17 files, 8134 in 15, 8128 and 8130 in
11 each). Those record what was true when written and rewriting them would destroy that. A
registry-backed port checker was prototyped and **rejected**: legitimate infrastructure ports
(8501 Streamlit, 8510 portal) and non-port four-digit numbers (a test count of 8397) make it cry
wolf, and a check that cries wolf is one nobody runs.

### 3.4 Test-depth floor — ✅ decided 2026-09-16 (H-D5)

**A raw test count is the wrong instrument.** By absolute count the alarming subject is the
single-cell *engine* with 4 tests — but it is 146 lines of deterministic scanpy glue and 4 tests is
proportionate. clinical-imaging looks well covered at 1,365 tests and is 65,000 lines, which is
thinner than it looks. Normalising by source size inverts the ranking:

| subject | tests | src LOC | per 100 |
|---|---|---|---|
| single-cell *(agent)* | 185 | 19,332 | **1.0** |
| rare-disease-diagnostic | 206 | 20,713 | **1.0** |
| neurology | 208 | 20,849 | **1.0** |
| tuberous-sclerosis | 92 | 6,370 | **1.4** |
| clinical-imaging | 1,365 | 65,227 | 2.1 |
| … | | | |
| cardiology | 1,966 | 33,853 | 5.8 |
| single-cell *(engine)* | 4 | 146 | 2.7 — *fine, and exempt anyway* |

**The metric:** tests per 100 non-test source lines.
**The scope:** subjects that emit **clinical output** — prose or a recommendation a clinician could
act on. Listed explicitly in `scripts/check_test_depth.py` with the reason for each, because
inferring it from "has a clinical eval case" is circular: the two worst faults found on 2026-09-16
(precision-biomarker, clinical-imaging) were both in subjects that had **no** eval case. A thin
suite on a data-loading engine costs a broken build; a thin suite on a subject that generates
clinical prose costs a wrong answer that reads exactly like a right one.

**The floor: 2.0** — half the median of 2.8. Below that is not a number anyone has to argue about.

**Enforced as a ratchet, not a cliff.** Four clinical subjects sit below 2.0. Making the floor
blocking today would either turn CI red on merge or force ~600 tests to be written in a hurry, and
tests written to satisfy a number are worth nothing. So:

| rule | status |
|---|---|
| **RATCHET** — no clinical subject falls below its recorded baseline (±0.15 noise) | **blocking**, in CI |
| **FLOOR 2.0** — a *new* clinical subject starts at or above it | **blocking**, in CI |
| **FLOOR 2.0** for the four already below | **advisory**, recorded here |

Backsliding is what actually decays, and the ratchet catches it the day it happens. Verified by
simulation: deleting 600 pharmacogenomics tests trips `REGRESSION` and exits 1.

```bash
$PY scripts/run_all_tests.py --json /tmp/t.json
$PY scripts/check_test_depth.py /tmp/t.json            # report
$PY scripts/check_test_depth.py /tmp/t.json --enforce  # the CI gate
$PY scripts/check_test_depth.py /tmp/t.json --update-baseline   # deliberate, reviewed in a PR
```

Baseline: `docs/build/test_depth_baseline.json`. Raising it is a PR like any other; lowering it
should be argued for in the PR description, not done quietly.

**The four below the floor, in priority order** — all three 1.0 agents are ~20,000 lines emitting
clinical prose, which is the exact profile of the faults found this week:
`neurology` · `rare-disease-diagnostic` · `single-cell (agent)` · then `tuberous-sclerosis` (1.4,
and it has its own construct-validity eval, so it is the least urgent).

---

### 3.5 Data licences — ✅ decided 2026-09-16 (H-D2)

H-D2 read: *"AlphaMissense is CC BY-NC-SA (non-commercial) across 132 tracked files, while the
platform is Apache-2.0 and the site welcomes commercial use."*

**Measured, the conflict is not in the repository.** No third-party dataset is redistributed here
— `git ls-files | grep -iE 'alphamissense|clinvar|oncokb|cosmic'` returns nothing. The 134 files
are code that *reads* the data and docs that explain how to obtain it. Apache-2.0 code reading
CC BY-NC-SA data on a user's own machine is the user's obligation, not a defect in the licence
grant on this code.

**The real defect was that the repository contradicted itself.**
`core/agents/rare-disease-diagnostic/docs/RARE_DISEASE_DIAGNOSTIC_AGENT_RESEARCH_PAPER.md` listed
AlphaMissense as **CC BY 4.0**, while `docs/build/ACQUISITION_MANIFEST.md` correctly said
**CC BY-NC-SA 4.0, non-commercial**. The error ran in the permissive direction — it told a
commercial reader the data was free to use. Corrected.

**Resolution: `DATA_LICENSES.md` at the repository root**, linked from the README's licence
section. It states the governing fact (nothing is redistributed, with the command to verify it),
names the three datasets that are gated in practice (AlphaMissense non-commercial; OncoKB and
COSMIC licence-required; OMIM registration), and lists **every** external dataset the platform
reads with a link to its own terms page — deliberately linking to the authority rather than
asserting terms in a table that will go stale.

It also records that the platform degrades cleanly without AlphaMissense *by construction*:
`ingest_vcf.py` checks `.exists()`, logs "AlphaMissense annotation will be skipped", and continues
on ClinVar significance alone. A first draft of that sentence claimed an `HCLS_SKIP_ALPHAMISSENSE`
env flag; no such flag exists, and it was replaced with the mechanism that does. A licence file
asserting a feature that is not there would be the same fault this platform keeps finding.

---

## Traps already paid for on this machine

1. **A count is not a cause.** `run_all_tests.py` reported "errors 36" with no traceback and cost a
   full CI cycle. It now prints the failing test ids and the pytest tail.
2. **`MilvusClient` connects eagerly** in `__init__` — construct it inside the lifespan try, never
   at module scope, or every import fails where Milvus is absent.
3. **`exec 9>&-` is load-bearing in three places** in `health-monitor.sh`. A started service
   inherits the flock fd and holds the lock for its own lifetime.
4. **`str(DataType.FLOAT)` is `"10"`, not `"FLOAT"`.** Use `.name`.
5. **Fetch is not persist.** Three ingest scripts logged "N validated" and wrote nothing.
6. **Absence of a name is not absence of behaviour.** `pharmacogenomics` was reported broken
   because it lacked a `_persist` function; it had always persisted by another route.
