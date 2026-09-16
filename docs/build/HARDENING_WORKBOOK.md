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

⚠️ **The baseline is recorded from CI, not from this box, and the first run proved why.** A
locally-captured baseline failed instantly in CI: **CI collects fewer tests than a developer
machine** — 300 cart tests there against 415 here, 665 biomarker against 709, 1,325 imaging
against 1,365 — because optional and GPU-gated dependencies are absent, so those suites skip or
fail to collect. A local run now measures at or above the CI baseline and passes, which is the
right asymmetry: the gate cannot be satisfied by a machine with more installed than the gate has.

That gap (~9% of the suite) is worth knowing on its own: **roughly 700 tests that run here never
run in CI.** It is not a defect this decision fixes, but it is the reason the number in CI is the
one that counts.

LOC is counted over `git ls-files` rather than a filesystem walk, so a scratch file in the working
tree cannot move the ratio.

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

### 3.6 pymilvus 3.0 — ✅ decided 2026-09-16 (H-D3)

Four dependabot PRs (#45 #49 #52 #54) proposing pymilvus 3.0 were held on the reasoning that
2.6.8 works and a client-API mismatch had already broken retrieval in five services once, so a
second would present identically. **Tested rather than assumed, that fear does not hold for 3.0 —
and the hold was hiding a different problem.**

Measured in an isolated venv against the live Milvus server:

| | |
|---|---|
| every symbol the code imports (`connections`, `utility`, `Collection`, `MilvusClient`, …) | **present in 3.0.0** |
| ORM connect + `list_collections` + `Collection.load()` against the running server | **works** — 113 collections, entities returned |
| `MilvusClient` against the same server | **works** |
| `hcls_common` Milvus + vector-search tests under 3.0.0 | **37 passed** |

**But every ORM call now warns:**

```
PyMilvusDeprecationWarning: `connections.connect` is an ORM-style PyMilvus API
and will be removed in PyMilvus 3.1. Use `MilvusClient` instead.
```

The codebase uses that ORM surface at roughly **70 call sites** (`connections` 25, `utility` 23,
`Collection` 22). So 3.0 is safe and **3.1 is the cliff** — and holding the PRs does not avoid that
work, it only delays discovering it.

**The hold was also masking drift.** The requirement files claimed `==2.4.1` (×4), `>=2.4.0,<2.6`,
`>=2.4.0` and `>=2.4` across 14 files while **2.6.8 was actually installed** — two of those pins
the running fleet already violates. The declared dependency and the deployed one had nothing to do
with each other.

**Decision: `pymilvus>=2.6.8,<3.1`, applied consistently to all 13 live requirement files.**

- it matches what is actually installed, so the files stop lying;
- it permits 3.0, which is tested and working, without forcing an upgrade on a clinical fleet;
- **`<3.1` is the point.** It converts the ORM removal from something that arrives silently in a
  fresh install into a deliberate decision someone has to make. The four dependabot PRs are closed
  in favour of this, since their intent — allow 3.0 — is satisfied.

**Follow-on work, now visible instead of deferred:** migrate ~70 ORM call sites to `MilvusClient`
before anyone lifts `<3.1`. `hcls_common.vector_search.search_collection` already handles both
client shapes, so the search path is done; what remains is connection management and
`utility.*` calls. Run `scripts/run_clinical_eval.py` after — it is the only detector that notices
retrieval quietly returning nothing.

---

### 3.7 Grounding and corpus — ✅ measured 2026-09-16

**The eval graded the answer and never asked whether retrieval contributed.** That gap hid two
things at once.

**44 of 113 Milvus collections are empty.** Not broken — *empty*. The service is up, the
collection is loaded, the query succeeds, and the result set has nothing in it. There is no error
anywhere in that chain, so the model answers from its own knowledge and the reply looks exactly
like a sourced one.

| subject | collections | empty | vectors |
|---|---|---|---|
| precision-autoimmune | 13 | **8** | 116 |
| neurology | 13 | 3 | 168 |
| clinical-trial | 13 | **11** | 178 |
| pharmacogenomics | 14 | 0 | 240 |
| single-cell | 11 | **8** | 279 |
| rare-disease-diagnostic | 13 | **10** | 440 |
| clinical-imaging | 12 | 2 | 440 |
| cart | 10 | 2 | 879 |
| precision-biomarker | 13 | 0 | 1,244 |

**Eight of ten subjects hold fewer than 500 vectors.** Every one of them was passing the clinical
eval.

**The single-cell agent answers with no evidence at all:**

```
/v1/sc/query -> answer 2,658 chars · evidence: [] · guidelines_cited: [] · confidence: 0.3
```

It passed `sc-tcell-marker` on model knowledge alone. The agent is honest enough to report
confidence 0.3; nothing surfaced it.

**Two checks now make this visible:**

- `scripts/run_clinical_eval.py` reads the evidence count out of whatever key a service uses and
  returns a new **`UNGROUNDED`** verdict — a right-looking answer with zero retrieved passages is
  not a pass for a platform whose claim is RAG over a curated corpus, and it counts as a failure.
  Verified: precision-biomarker `PASS 30 evidence`, single-cell `UNGROUNDED 0 evidence`.
- `scripts/check_corpus.py` reports per-subject collections, empties and vector counts against a
  floor, and runs inside `scripts/reboot_check.py`.

Neither fixes the corpus. Seeding it is content work — `docs/build/CORPUS_SEEDING.md` — but it is
now a number someone can see rather than an absence nobody can.

### 3.8 One pre-commit guard, not two — ✅ done 2026-09-16

The repo shipped `.pre-commit-config.yaml` (gitleaks, `check-added-large-files --maxkb=5120`,
yaml/json checks) **and** a hand-written `scripts/pre-commit-hook.sh`. Git calls the hand-written
one, and it never delegated to the framework — so the documented config was dead weight and
gitleaks never ran locally.

The hand-written hook also **explicitly exempted `docs/assets/videos/*.mp4`** from its 5 MB limit.
That exemption is exactly how 899 MB of video reached the history: 123 blobs, 84 of them
superseded re-encodes of the same 22 files, stripped in the H-D4 rewrite that recovered 233 MB per
clone.

Now: the hook runs `pre-commit run` when the framework is installed, and re-committing a video
requires an explicit `HCLS_ALLOW_VIDEO_COMMIT=1`. Adding a video is still allowed; doing it by
accident is not. Verified both paths — blocked without the flag, permitted with it.

---

### 3.9 The honesty gate was withholding correct answers — ✅ fixed 2026-09-16

The 27-case eval came back 24/27, and one of the misses was the **flagship** question:

> *Which genes cause tuberous sclerosis complex and which pathway is dysregulated?*

Run five times against the same service, the same question returned:

```
run 1: WITHHELD by honesty gate    585 chars
run 2: correct (names TSC1/TSC2) 4836 chars
run 3: WITHHELD by honesty gate    585 chars
run 4: correct (names TSC1/TSC2) 4888 chars
```

**The gate was blocking a textbook genetics answer about half the time**, on
*"Diagnostic-certainty overclaim"* — because a correct answer naturally says how the diagnosis is
established. A gate that withholds the right answer half the time teaches people to route around
it, and it made the eval non-deterministic, which is why "27/27" was not reproducible.

**Fix: subject-scope the diagnostic-certainty rules**, the same treatment the regulatory rules got
when enforcement was turned on. They block on a **self-reference or a patient reference** and
degrade to `warn` for statements about how diagnosis works in general:

| | |
|---|---|
| "The result confirms the diagnosis of the disorder." | **block** |
| "The patient's diagnosis is confirmed by the variant." | **block** |
| "Our analysis provides a definitive diagnosis." | **block** |
| "Genetic testing confirms the diagnosis in 85% of cases." | warn |
| "A definitive diagnosis requires molecular confirmation of TSC1 or TSC2." | warn |

A `warn` is published with the disclaimer attached; only a `block` or a refuted claim is withheld.
Verified live afterwards: **5 of 5 runs correct**, none withheld.

**Writing that test found two pre-existing bugs in the safety rules themselves.**

1. **`\b100%\b` never matched anything.** `\b` after `%` requires a following word character, so
   *"100% of cases"*, *"100% response rate"* and *"a 100%-effective drug"* all passed the
   absolute-certainty rule. The existing test only appeared to cover it because the sentence also
   said "cures", which tripped a different rule. The single commonest overclaim token had never
   fired.
2. **Only the active voice was matched.** *"The patient's diagnosis is confirmed"* — the more
   natural way to say the dangerous thing — was not caught at all.

Both fixed and tested. `_SELF_REF` was also widened to the platform's own output nouns
("the result", "the report", "the finding"), without which the sentence the rule exists for
degraded to a warning.

---

### 3.10 The flagship's corpus — ✅ measured and repaired 2026-09-16

`scripts/check_corpus.py` reported **zero collections** for the Tuberous Sclerosis program, the
flagship disease vertical. Four separate things were true at once:

**1. It runs on the in-memory store.** `TSC_USE_MILVUS: "1"` is set in the program's
`docker-compose.yml`, but the supervisor launches it with uvicorn, so that never applies. Its
corpus is rebuilt from `SEED_CORPUS` at startup — never stale, but invisible to every corpus
check, dashboard and operator, and unable to grow beyond what is hard-coded. Documented in
`.env.example`; the runtime default is deliberately left alone, because with a static corpus
in-memory is defensible and switching the flagship's store is not a change to make silently.

**2. That corpus is six chunks.** Six. For the flagship.

**3. The Milvus path never flushed.** `client.insert()` with no flush means `num_entities` keeps
reporting **0** while the rows sit in a growing segment. The loader prints *"Ingested 6 chunks"*
and every observer sees an empty collection — success reported, nothing visible. Exactly the
failure mode this platform keeps finding, in the ingest path this time. Fixed.

**4. `upsert()` was an append.** The collection was created with `auto_id=True`, which
`lib/hcls_common/ingest_persist.py` explicitly warns against — *"PK differs per collection, never
auto_id; a content-derived id makes re-ingest idempotent."* Running `scripts/load_rag.py` three
times produced **18 rows of a 6-chunk corpus**, and duplicate passages skew retrieval while
looking like a healthy corpus. Now a deterministic `blake2b` content key + `client.upsert()`.
Verified: three consecutive loader runs → **6 live rows, 6 distinct keys**.

A collection created before this change is detected and left alone rather than dropped — by then
it may hold real ingested literature — and the store logs how to migrate it.

**And the measurement tool was wrong too.** `check_corpus.py` used `num_entities`, which counts
soft-deleted rows until compaction: after the upsert fix, `tsc_literature` read **18** while
holding **6**. It now uses `count(*)`. The "44 empty" figure was unaffected (soft deletes only
inflate, never deflate) and was re-verified by flushing every reportedly-empty collection: all 44
are genuinely empty.

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
