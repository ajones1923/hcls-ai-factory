# Workbook — Finishing the HCLS AI Factory

**Date:** 2026-08-17 · **Companion:** [FACTORY_COMPLETION_PRD.md](FACTORY_COMPLETION_PRD.md)
**Repo-only.** Not published to hcls-ai-factory.org.

Work this in order. Every step has a **✅ Done when** you can check without judgement. Do not move to
the next phase until its gate passes — §8 of the PRD explains why, and it is the single biggest
time-saver in this document.

**Set this once per shell:**

```bash
cd /home/adam/projects/hcls-ai-factory
export PY=.venv/bin/python
```

---

## Phase 0 · Foundation

*No accounts. No downloads over 1 GB. Everything here is unblocked today.*

### 0.1 Create `.env` — there is currently no file at all

```bash
cp .env.example .env
chmod 600 .env
```

Fill the four required values. Leave gated ones blank for now:

```bash
echo "HCLS_API_KEY=$(openssl rand -hex 24)"        >> .env
echo "HCLS_MINIO_USER=hcls"                        >> .env
echo "HCLS_MINIO_PASSWORD=$(openssl rand -hex 16)" >> .env
echo "GRAFANA_PASSWORD=$(openssl rand -hex 12)"    >> .env
```

> **`minioadmin/minioadmin` is the default in some compose examples. Never keep it.** It is a known
> remediation item and an obvious way in.

**✅ Done when:** `docker compose -f docker-compose.dgx-spark.yml config >/dev/null` prints no
warnings about unset variables.

### 0.2 Prove the security gate fails closed

The gate ships **off** deliberately (trusted-network posture). You just turned it on — verify it.

```bash
$PY -c "
from fastapi import FastAPI; from fastapi.testclient import TestClient
import os; os.environ['HCLS_API_KEY']='testkey'
from hcls_common.api_auth import install_api_key_auth
app=FastAPI(); install_api_key_auth(app,'demo')
@app.get('/q')
def q(): return {'ok':True}
c=TestClient(app)
print('no key   ->', c.get('/q').status_code)
print('wrong    ->', c.get('/q', headers={'X-API-Key':'nope'}).status_code)
print('right    ->', c.get('/q', headers={'X-API-Key':'testkey'}).status_code)
print('health   ->', c.get('/health').status_code)"
```

**✅ Done when:** `401 / 401 / 200 / 200`. Health must stay open or probes fail closed and the
platform looks dead.

### 0.3 Free the GPU — it cannot allocate right now

99 GiB of the 119 GiB unified memory is host page cache. On unified memory the cache and the GPU
compete for the same pages, so CUDA fails at context creation.

```bash
free -g                                    # look at the buff/cache column
sync && sudo sh -c 'echo 3 > /proc/sys/vm/drop_caches'
free -g                                    # cache should collapse

core/engines/therapeutic-discovery/venv/bin/python -c "
import torch
print('device:', torch.cuda.get_device_name(0))
print('free/total GiB:', [round(x/2**30,1) for x in torch.cuda.mem_get_info()])
a=torch.randn(4096,4096,device='cuda'); b=a@a; torch.cuda.synchronize()
print('matmul ok, peak GiB', round(torch.cuda.max_memory_allocated()/2**30,2))"
```

**✅ Done when:** `mem_get_info` returns numbers instead of raising, and the matmul completes.

> **Make this permanent (R5).** Add the drop to whatever brings the platform up, before any GPU
> service starts. A manual step you have to remember is a step that fails at 2 a.m.

### 0.4 Install the free packages and binaries

```bash
$PY -m pip install \
  pysam cyvcf2 rdkit py3dmol stmol \
  monai monai-deploy-app-sdk SimpleITK pyradiomics scikit-image \
  e3nn einops wandb openai langgraph nemoguardrails \
  minio redis psycopg2-binary nvidia-ml-py \
  opentelemetry-api opentelemetry-sdk \
  plotly pyvis PyPDF2 streamlit-chat imageio-ffmpeg decorator fhir.resources

sudo apt-get update && sudo apt-get install -y \
  samtools bcftools tabix ffmpeg openjdk-17-jre-headless

curl -s https://get.nextflow.io | bash && sudo mv nextflow /usr/local/bin/
```

GATK: download the release zip from <https://github.com/broadinstitute/gatk/releases>, unzip, put
`gatk` on `PATH`.

> **Never install `mkdocs-material` into `.venv`.** It pulls `zarr` and `fast-array-utils`, whose
> pytest plugins abort collection for three subjects. That cost 373 tests once already. Docs
> toolchain lives in `.venv-docs`.

> **`nvidia-ml-py` provides the `pynvml` import.** Package and module names differ.

**✅ Done when:**

```bash
$PY -c "import pysam,cyvcf2,rdkit,monai,SimpleITK,e3nn,openai,langgraph,minio,redis,pynvml,plotly,fhir.resources; print('packages ok')"
for b in samtools bcftools tabix bgzip nextflow gatk ffmpeg; do command -v $b >/dev/null || echo "MISSING $b"; done; echo "binaries checked"
```

### 0.5 Decide D3, then stand up HCLS's Milvus

**Port `:19530` is already taken** by `imaging-milvus-standalone`, and it holds only the twelve
`imaging_*` collections. The other seven agents have nowhere to write.

Answer **D3** (PRD §7) first. The lower-risk option is a second instance on its own port:

```bash
docker ps --format '{{.Names}} {{.Ports}}' | grep -i milvus     # confirm the conflict
# then either: point HCLS at a NEW port in .env and compose,
# or: accept the shared instance and namespace every HCLS collection hcls_*
```

**✅ Done when:** HCLS connects to a Milvus it owns and `utility.list_collections()` returns
without the imaging stack's names being at risk from an HCLS reset.

### 🚦 Phase 0 gate

- [ ] `.env` exists, `compose config` clean
- [ ] auth returns 401/401/200/200
- [ ] CUDA allocates and completes a matmul
- [ ] all packages import; 7 binaries on `PATH`
- [ ] D3 answered; HCLS Milvus reachable

---

## Phase 1 · Bring-up

### 1.1 Answer D1 — compose or host processes

You currently maintain both. The compose header calls itself the portable target; `health-monitor.sh`
runs the live system; the per-service `venv/` directories that path needs exist for **3 of ~20**
services. Pick one and demote the other in writing.

**Recommendation: compose.** 26 services are already declared, it is reproducible for a stranger
following the quickstart, and it removes the venv-per-service problem entirely.

### 1.2 Build the images

There are **zero** HCLS engine/agent images on the box today.

```bash
docker compose -f docker-compose.dgx-spark.yml build 2>&1 | tee /tmp/build.log
grep -iE "error|failed" /tmp/build.log | head
```

> **The portal Dockerfile's build context is the repo root.** A root `.dockerignore` is
> **required** — without it the build tars hundreds of GB.

**✅ Done when:** `docker images | grep -c hcls` ≥ 20 and the log has no errors.

### 1.3 Start, and watch what actually comes up

```bash
docker compose -f docker-compose.dgx-spark.yml up -d
sleep 60
docker compose -f docker-compose.dgx-spark.yml ps --format 'table {{.Name}}\t{{.Status}}'
docker compose -f docker-compose.dgx-spark.yml ps | grep -v healthy   # the ones to chase
```

**✅ Done when:** every service reports healthy. Chase them one at a time — `docker compose logs
<service> --tail 50`.

### 1.4 Seed the agent corpora — this is what unlocks the demos

Eight agents are RAG services. Without a populated vector DB they return an honest 503 and **eight
demos stay unavailable**. This is the highest-leverage step in the whole workbook.

```bash
ls scripts/*seed* scripts/*ingest* core/agents/*/scripts/*seed* 2>/dev/null
# then per agent, e.g.:
$PY scripts/seed_corpus.py --agent pharmacogenomics
```

**✅ Done when:** each agent's collections exist and a query returns hits:

```bash
$PY -c "
from pymilvus import connections, utility
connections.connect(host='localhost', port='19530')
print(sorted(utility.list_collections()))"
```

### 1.5 Port convention — do not fight it

**Registry advertises the UI; the API is UI + 1.** CI cross-checks the registry against the process
supervisor, so a port changed in one place and not the other fails the build. That guard is
deliberate; if it fires, fix the drift rather than the guard.

### 1.6 Survive a reboot

```bash
grep -c "restart: unless-stopped" docker-compose.dgx-spark.yml   # should cover every service
sudo reboot
# after it returns, wait 5 minutes:
$PY scripts/run_demo.py --check-all
```

**✅ Done when:** readiness after reboot equals readiness before it. Anything that dies on reboot is
a service you will be restarting by hand forever — the cardiac imaging demo already has this
problem, per its own operating doc.

### 🚦 Phase 1 gate

- [ ] D1 answered and recorded
- [ ] all compose services healthy
- [ ] corpora seeded for all 8 agents
- [ ] `--check-all` ≥ 13/17
- [ ] readiness survives a reboot unattended

---

## Phase 2 · Free accounts

| Credential | Where | Unlocks |
|---|---|---|
| `ANTHROPIC_API_KEY` | console.anthropic.com | reasoning in every agent — **11 of 17 demos** |
| HuggingFace token | huggingface.co | BGE embeddings (242 refs), ESMFold weights |
| `NCBI_API_KEY` | ncbiinsights.ncbi.nlm.nih.gov | PubMed ingest 3/s → 10/s |

```bash
$PY -m pip install sentence-transformers
$PY -c "from sentence_transformers import SentenceTransformer as S; S('BAAI/bge-small-en-v1.5'); print('BGE ok')"
$PY -c "
import anthropic
print(anthropic.Anthropic().messages.create(
  model='claude-sonnet-4-5-20250929', max_tokens=20,
  messages=[{'role':'user','content':'reply OK'}]).content[0].text)"
```

> **The code pins `claude-sonnet-4-20250514` in 36 places and `claude-sonnet-4-6` in 23.** Both are
> old. Decide a model policy and set it in one place rather than 59.

**✅ Done when:** an agent returns cited output end-to-end:

```bash
curl -s -X POST localhost:8508/query -H "X-API-Key: $HCLS_API_KEY" \
     -H 'Content-Type: application/json' -d @demo/requests/pharmacogenomics_query.json | head -40
```

### 🚦 Phase 2 gate
- [ ] BGE loads · Claude answers · one agent returns cited output with sources

---

## Phase 3 · Gated software — the long pole

**Everything before and after this phase is independent of it.** If NGC stalls, keep moving.

### 3.1 NGC

```bash
# ngc.nvidia.com -> Setup -> API Key -> Generate
echo "NGC_API_KEY=<paste>" >> .env
echo "$NGC_API_KEY" | docker login nvcr.io --username '$oauthtoken' --password-stdin
```

**✅ Done when:** *Login Succeeded*. Nothing below works otherwise.

### 3.2 Answer D2 first — before building anything on it

```bash
docker manifest inspect nvcr.io/nim/nvidia/molmim:1.0.0 | grep -i architecture
docker manifest inspect nvcr.io/nim/mit/diffdock:2.2.0  | grep -i architecture
```

If `amd64` only: they can **never** run on this aarch64 box. Then MolMIM/DiffDock demos are
permanently **BURST**, the narration says "elastic burst", and you build the remote path now rather
than discovering it later.

### 3.3 Pull in dependency order

```
CUDA PyTorch (G1) ──► Parabricks (G2) ──► real FASTQ→VCF
        └──► MolMIM (G3) + DiffDock (G4) ──► the flagship discovery demo
```

### 3.4 Real genomics input

The FASTQ files under `core/engines/genomic-foundation/data/input/` are **0-byte placeholders**, and
**GRCh38 is not present**. Answer **D5** — a downsampled chromosome may prove the pipeline for a
fraction of 200 GB.

```bash
# GRCh38 (~3 GB)
# https://ftp.ncbi.nlm.nih.gov/genomes/  -> primary assembly + index
# HG002 (~200 GB, or a subset)
# https://ftp-trace.ncbi.nlm.nih.gov/giab/
```

**✅ Done when:** Parabricks writes a real VCF and the store reports **Ts/Tv ≈ 2.0** — the same
number the site already publishes, now produced here.

### 🚦 Phase 3 gate
- [ ] NGC login succeeds
- [ ] **D2 answered and written down**
- [ ] CUDA PyTorch behind the platform venv
- [ ] Parabricks produces a real VCF with Ts/Tv ≈ 2.0

---

## Phase 4 · Demonstrations

### 4.1 Close the gap to 17/17

```bash
$PY scripts/run_demo.py --check-all
```

The runner **refuses to run a LIVE demo whose service is unreachable** rather than returning a canned
result. When it says not-ready, that is the honest state — fix the service or change the label.

### 4.2 Run each and record a transcript

```bash
for d in E1 E2 E3 E4 E5 E6 E7 E8 A1 A2 A3 A4 A5 A6 A7 A8 P1; do
  echo "=== $d ==="; $PY scripts/run_demo.py $d 2>&1 | tail -20
done
ls demo/transcripts/
```

**✅ Done when:** 17 transcripts exist and re-running produces the same output twice.

### 4.3 Make the labels true

Every demo carries exactly one label — **LIVE / REPRESENTATIVE / BURST**. After the runs, reconcile:

- `LIVE` only where the service actually answered on this box
- `BURST` where it ran on a remote GPU — say "elastic burst", never "all on one box"
- `REPRESENTATIVE` where the result was pre-computed — **say so out loud in the demo**

> The catalogue's headline counts were hand-written and disagreed with its own body once already.
> **Generate them from `--check-all`** (R20) so they cannot drift again.

### 🚦 Phase 4 gate
- [ ] 17/17 prerequisites met
- [ ] 17 transcripts recorded and reproducible
- [ ] every label matches what actually ran
- [ ] catalogue counts generated, not typed

---

## Phase 5 · Quality debt

| # | Task | Why |
|---|---|---|
| R21 | Remove `continue-on-error` from the subject-suite job and make it green | CI currently reports **success while that job fails**. For a project whose argument is honesty-by-construction, that is the wrong irony to ship. |
| R11 | `validate_registry.py --probe` — fail CI on a `live` capability that does not answer | Turns the honesty rule from a periodic audit into an enforced invariant. |
| R22 | Call `require_valid_input()` / `honesty_flags()` in handlers | Coverage is **1/12** for both. The governance page states this honestly; closing it makes the claim strong instead of merely honest. |
| R23 | Rename the 12 stdlib-shadowing modules | `src/collections.py` shadows the stdlib in 11 subjects; it has produced two false audit results already. |
| R24 | structural-biology aggregator on `:8581`, or leave `planned` | Five services, no aggregator. Status follows evidence. |
| R25 | Agree a test-depth floor for clinical-output subjects | Range is 4 → 1,966. Pass rate alone is uninformative. |

---

## Traps already paid for on this machine

Each of these cost real time. They are listed so they cost nothing the second time.

1. **PEP 668** — system Python is externally managed. Everything goes in a venv.
2. **Corrupt editable install** — if `import hcls_common` fails, check
   `~/.local/lib/python3.12/site-packages/__editable___hcls_common_*_finder.py`. An empty `MAPPING`
   dict means it is broken; delete both `__editable__*` files and reinstall. This blocked the merge
   gate entirely.
3. **Never run bare `pytest` across subjects** — `src/collections.py` shadows the stdlib. Use
   `scripts/run_all_tests.py`, which withholds `src/` where needed.
4. **`zarr` / `fast-array-utils` register pytest plugins that abort collection** — the harness passes
   `-p no:zarr -p no:fast_array_utils`. `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1` "fixes" it and breaks 38
   async tests instead.
5. **A uniform result across every subject means the harness is wrong, not the code.** All 17 once
   looked broken because `--timeout` was passed without `pytest-timeout` installed.
6. **Measure through the right interpreter.** The platform venv is CPU-only; CUDA lives in
   `core/engines/therapeutic-discovery/venv`. A measurement through one interpreter is a claim about
   that interpreter, not the machine — this produced a wrong hardware score in a published document.
7. **The page cache eats unified memory.** Drop it before GPU work.
8. **`:19530` belongs to the imaging stack.** Check before assuming it is yours.
9. **Health checks prove liveness, not progress.** Check that work is *advancing*, not just that a
   port answers.
10. **Check CI before polling a deploy.** A `.gitignore` rule once swallowed a directory, CI failed,
    and the deploy never ran — time was spent watching a deploy that was never coming.

---

## One-command health check

Keep this to hand; it answers "where am I" in ten seconds.

```bash
#!/usr/bin/env bash
cd /home/adam/projects/hcls-ai-factory
echo "== services ==";  docker compose -f docker-compose.dgx-spark.yml ps 2>/dev/null | grep -c healthy
echo "== registry ==";  .venv/bin/python scripts/validate_registry.py | tail -2
echo "== tests ==";     .venv/bin/python scripts/run_all_tests.py 2>&1 | tail -1
echo "== demos ==";     .venv/bin/python scripts/run_demo.py --check-all 2>&1 | tail -1
echo "== gpu ==";       nvidia-smi --query-gpu=utilization.gpu,power.draw --format=csv,noheader
echo "== memory ==";    free -g | awk 'NR==2{print "  cache",$6"G of",$2"G"}'
```

---

## Progress tracker

| Phase | Gate | Done |
|---|---|---|
| 0 · Foundation | `.env`, auth closed, CUDA allocates, HCLS Milvus | ☐ |
| 1 · Bring-up | services healthy, corpora seeded, survives reboot | ☐ |
| 2 · Free accounts | BGE + Claude + one cited agent answer | ☐ |
| 3 · Gated | NGC, D2 answered, Parabricks VCF at Ts/Tv ≈ 2.0 | ☐ |
| 4 · Demos | 17/17, transcripts, honest labels | ☐ |
| 5 · Debt | CI blocking and green, probe in CI | ☐ |

**Acceptance (PRD §3):** A1 every `live` capability probes · A2 17/17 demos run · A3 survives reboot ·
A4 clean clone reproduces · A5 site claims nothing the box does not deliver.
