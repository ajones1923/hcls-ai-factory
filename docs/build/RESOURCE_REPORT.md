# batch-01 — Resource Consumption Report

**Run:** `batch-01-hcls-aif-core-updates.txt` · **Host:** NVIDIA DGX Spark (GB10, 20 Grace cores,
128 GB unified LPDDR5x, aarch64) · **Date:** 2026-08-15

Sampling ran from the top of the batch to its close, because these figures cannot be
reconstructed afterwards.
`tmp/batch-01-run/sampler.sh` writes one CSV row every 30 s; raw data is in
`tmp/batch-01-run/resources.csv`.

---

## 1. Headline

| Metric | Value |
|---|---|
| Sampled window | **203 minutes**, 406 samples |
| CPU | mean **4.4 %**, peak **34.2 %** of 20 cores |
| Memory | peak **19 GB of 119 GB** (16.3 %) |
| **GPU utilisation** | **0 % — never engaged** |
| GPU power | 12.3 – 16.0 W (idle floor) |
| GPU energy | **45.9 Wh** — entirely idle draw |
| Disk delta | **+13 GB** |
| Disk free | 770 GB |
| Network | 3.4 GB in · 4.4 GB out |

## 2. The finding buried in these numbers

**The GPU was never used.** 0 % utilisation across 203 minutes on a machine whose value proposition
is the GB10. Power never left the 12–16 W idle band. The 45.9 Wh of GPU energy was spent doing
nothing.

Cause: `torch 2.10.0+cpu`, `torch.cuda.is_available() = False`. Every ML path on this box is
CPU-bound or non-functional. Unblocking it is gated behind NGC credentials
(`GATED_SOFTWARE_PRD.md`, G1).

Peak memory of **19 GB of 119** tells the same story from the other side: a full 17-subject,
8,402-test sweep used 16 % of the machine. The hardware is not the constraint — provisioning is.

## 3. CPU profile

Mean 4.4 % against a 34.2 % peak. The peaks align with the three full test sweeps and the
`pip install` of duckdb / statsmodels / biopython / peft. Nothing in this run was compute-bound;
the work was I/O and analysis.

Grace cores available: 20. Effective parallelism used: low — pytest ran single-process per subject.
A parallel harness (`pytest -n`) would cut sweep time materially and is worth adding.

## 4. Disk

+13 GB over the run, from:

- the repo `.venv` (created this session, `--system-site-packages`)
- duckdb, statsmodels, biopython, peft and their transitive dependencies
- pytest caches across 17 subjects

770 GB free. Not a constraint.

## 5. Network

3.4 GB in / 4.4 GB out. Inbound is PyPI wheels. Outbound is **not** attributable to this run alone —
the box concurrently runs the Marketing AI Factory (`intel-capture-*`), a cardiac-imaging Milvus
stack, and Prometheus/Grafana/DCGM. Counters are host-wide, so treat these two figures as an upper
bound, not as batch-01's own traffic.

## 6. AI model usage

| Model | Use |
|---|---|
| Claude Opus 4.6 (this session) | all analysis, code changes and documents |
| — | **no** ElevenLabs calls — no narration was generated |
| — | **no** Nano Banana / Gemini calls — no images were generated |
| — | **no** local model inference — GPU never engaged |

Token consumption is not measurable from inside this environment; the harness does not expose a
counter. What can be stated precisely is that **no paid generation API was called** during batch-01:
no TTS credits, no image credits, no NIM inference.

## 7. Honest caveats

1. **The sampled window is not wall-clock session time.** Sampling began a few minutes after the
   batch started, and the 203 minutes include intervals where the machine was idle between
   commands. It measures the run's footprint, not continuous compute.
2. **The host is shared.** CPU, memory, disk and network figures include the Marketing AI Factory,
   the cardiac demo stack and the monitoring containers. Only the GPU figure is unambiguous — it was
   0 % throughout, so nothing on the box used it either.
3. **GPU energy is idle draw**, not work. Reporting 45.9 Wh should not be read as compute cost.
4. **Token counts are unavailable**, not zero. See §6.

## 8. Reproducing

```bash
# start sampling before any work
tmp/batch-01-run/sampler.sh &

# summarise at any point
.venv/bin/python - <<'PY'
import csv
r=list(csv.DictReader(open("tmp/batch-01-run/resources.csv")))
f=lambda k:[float(x[k]) for x in r if x[k] not in ("NA","")]
cpu,gp=f("cpu_pct"),f("gpu_pwr_w")
print(f"{len(r)} samples · CPU mean {sum(cpu)/len(cpu):.1f}% · GPU {sum(gp)*30/3600:.1f} Wh")
PY
```

Columns: `ts, cpu_pct, mem_used_gb, mem_pct, gpu_util_pct, gpu_pwr_w, disk_used_gb, disk_free_gb,
net_rx_mb, net_tx_mb, load1`.
