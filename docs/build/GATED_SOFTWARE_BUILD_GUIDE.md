# Gated Software — Build Guide

Step-by-step install for everything the HCLS AI Factory needs that requires credentials or an
accepted licence. Companion to `GATED_SOFTWARE_PRD.md`.

**Target:** DGX Spark — GB10, 128 GB unified LPDDR5x, 20 Grace cores, **aarch64**.
Architecture matters throughout: several NVIDIA artefacts are x86-only, and that changes the plan.

> **Rule for this entire guide:** a capability's registry status changes **only after** a real input
> has produced a real recorded output. Two capabilities were already found marked `live` while
> unreachable. Do not repeat that.

---

## Step 0 — the ungated installs (do these first, no credentials)

```bash
cd /home/adam/projects/hcls-ai-factory
.venv/bin/pip install scanpy anndata pysam cyvcf2 rdkit

# Nextflow — the orchestrator named throughout the architecture, currently absent
curl -s https://get.nextflow.io | bash && sudo mv nextflow /usr/local/bin/
nextflow -version
```

Verify nothing regressed:

```bash
.venv/bin/python scripts/run_all_tests.py     # expect 17 subjects, 8,402 passed
```

---

## Step 1 — NGC account and entitlements  *(GR-1)*

1. Create/sign in at <https://ngc.nvidia.com>, then **Setup → API Key → Generate**.
2. Store it — never commit it:

```bash
printf 'NGC_API_KEY=<paste>\n' >> .env      # .env is gitignored (.gitignore:76-77,102)
```

3. Log in to the registry:

```bash
echo "$NGC_API_KEY" | docker login nvcr.io --username '$oauthtoken' --password-stdin
```

**Checkpoint:** `docker login` must report *Login Succeeded*. If it does not, nothing below will
work — stop and resolve entitlements first.

---

## Step 2 — CUDA PyTorch for GB10 / aarch64  *(G1, GR-2)*

This is the keystone. **Prefer the NGC container over pip wheels** — the container ships a runtime
matched to the GPU, and aarch64 + Blackwell (sm_121) wheel coverage on the public index is
unreliable.

```bash
docker pull nvcr.io/nvidia/pytorch:25.01-py3          # check NGC for the current aarch64 tag
docker run --rm --gpus all nvcr.io/nvidia/pytorch:25.01-py3 \
  python -c "import torch; print(torch.__version__, torch.cuda.is_available(), torch.cuda.get_device_name(0))"
```

**Expected:** `True` and `NVIDIA GB10`.

If you must use pip instead, install into the repo venv from NVIDIA's index and verify the compute
capability is actually supported:

```bash
.venv/bin/pip install --index-url https://download.pytorch.org/whl/cu126 torch torchvision
.venv/bin/python -c "import torch; print(torch.cuda.get_arch_list())"   # must include sm_121
```

**Do not accept `torch.cuda.is_available() == True` as sufficient.** Confirm real utilisation:

```bash
nvidia-smi --query-gpu=utilization.gpu,power.draw --format=csv -l 1 &
.venv/bin/python -c "
import torch; a=torch.randn(8000,8000,device='cuda')
for _ in range(50): a=a@a.T/1e4
torch.cuda.synchronize(); print('ok')"
```

Utilisation must leave 0% and power must rise above the ~13–15 W idle floor measured on this box.

---

## Step 3 — Parabricks  *(G2, GR-3)*

```bash
docker pull nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1
```

Run a real sample end to end (HG002 is publicly consented — the correct sample for a public demo):

```bash
docker run --rm --gpus all \
  -v /home/adam/projects/hcls-ai-factory/hcls-ai-factory-core-data:/data \
  nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1 \
  pbrun fq2bam --ref /data/reference/GRCh38.fa \
    --in-fq /data/fastq/HG002_R1.fastq.gz /data/fastq/HG002_R2.fastq.gz \
    --out-bam /data/out/HG002.bam
```

**Record the wall-clock.** The platform's public claim rests on it. If Parabricks is unavailable,
the genomics engine serves pre-computed results and every demo must say **REPRESENTATIVE**.

---

## Step 4 — BioNeMo NIMs: MolMIM and DiffDock  *(G3/G4, GR-4)*

```bash
docker pull nvcr.io/nim/nvidia/molmim:latest
docker pull nvcr.io/nim/mit/diffdock:latest

docker run -d --gpus all --name molmim  -p 8001:8000 -e NGC_API_KEY nvcr.io/nim/nvidia/molmim:latest
docker run -d --gpus all --name diffdock -p 8002:8000 -e NGC_API_KEY nvcr.io/nim/mit/diffdock:latest

curl -f http://localhost:8001/v1/health/ready && curl -f http://localhost:8002/v1/health/ready
```

**⚠️ Architecture check.** If either image is x86-only, it will not run on this box. Then:

- run them on remote GPUs over the private mesh, and
- describe it as **elastic burst** in every demo and script — never "all on one box".

Prove they work before promoting status:

```bash
# one real generation and one real pose, saved as evidence
curl -s -X POST http://localhost:8001/generate -H 'Content-Type: application/json' \
  -d '{"smi":"CC(=O)Oc1ccccc1C(=O)O","num_molecules":5}' | tee /tmp/molmim-proof.json
```

---

## Step 5 — ESMFold weights  *(G5, GR-5)*

Accept the model licence on HuggingFace (`facebook/esmfold_v1`), then:

```bash
export HF_TOKEN=<token>          # keep in .env, never committed
cd core/engines/structural-biology
PYTHONPATH=.:./src ../../../.venv/bin/uvicorn esmfold_service:_app_factory --factory --port 8570
curl -f http://localhost:8570/health
```

> `structural-biology-engine` is `planned` because **nothing binds :8581** — this engine is five
> services (8570 esmfold, 8571 protein-search, 8576 developability, 8577 immunogenicity,
> 8578 proteinmpnn) with no aggregator. Either build an aggregator on 8581 or leave the engine-level
> capability `planned` and promote only the model-level ones.

---

## Step 6 — RFdiffusion + DGL  *(G7 — needs G1 first)*

```bash
cd core/engines/structural-biology/vendor_rfdiffusion
../../../../.venv/bin/pip install -e .
../../../../.venv/bin/pip install dgl -f https://data.dgl.ai/wheels/cu126/repo.html
```

Weights are gated separately; follow the RFdiffusion repository instructions.

**Do not build the vendored tree into the engine image.** It is excluded from both the Dockerfile
and `scripts/run_all_tests.py` on purpose.

---

## Step 7 — reconcile the registry with reality  *(GR-6)*

```bash
.venv/bin/python - <<'PY'
import json, urllib.request
caps=json.load(open("lib/hcls_common/capabilities.json"))["capabilities"]
bad=[]
for c in caps:
    if c.get("status")!="live" or not c.get("endpoint"): continue
    try: code=urllib.request.urlopen(f"http://{c['endpoint']}/health", timeout=3).status
    except Exception as e: code=type(e).__name__; bad.append(c["id"])
    print(f"{c['id']:42s} {c['endpoint']:22s} {code}")
print("\nlive but unreachable:", bad or "none")
PY
```

Anything listed as unreachable must be fixed **or** demoted to `planned`. That check is how
`singlecell-compute` and `structural-biology-engine` were caught.

Then run the full gate:

```bash
.venv/bin/python scripts/validate_registry.py     # enforces UI/UI+1 and cross-checks health-monitor.sh
.venv/bin/python scripts/run_all_tests.py
```

---

## Step 8 — label the 17 demos honestly  *(GR-7)*

Each demo is marked exactly one of:

| Label | Meaning |
|---|---|
| **LIVE** | Ran now, in front of the audience, on real input, on this box |
| **REPRESENTATIVE** | Pre-computed or curated result standing in for a long step — say so |
| **BURST** | Ran live, but on remote GPUs over the private mesh — say "elastic burst" |

A capability that failed Step 7 cannot be demoed as LIVE.

---

## Rollback

Every step is additive. To undo:

```bash
docker rm -f molmim diffdock
docker rmi nvcr.io/nim/nvidia/molmim:latest nvcr.io/nim/mit/diffdock:latest
.venv/bin/pip uninstall -y torch torchvision dgl
```

The repo venv is disposable — `rm -rf .venv` and rebuild per `BUILD_GUIDE.md` §1.
