# Reboot check (PRD A3)

**A3 is the one acceptance criterion still unmet:** *the platform survives a reboot unattended.*
Everything else has been verified by running it. This has not — no reboot has happened since the
bring-up, so the boot path is the only part of the system still taken on trust.

Run this after the next reboot. It takes about five minutes.

## What was hardened in advance

Three things were known to break across a restart and are now handled in `health-monitor.sh`,
which cron runs every 5 minutes:

| Failure | Why it is invisible | Fix |
|---|---|---|
| Services start with **no `ANTHROPIC_API_KEY`** | They degrade to retrieval-only and still report healthy | the supervisor sources `.env` |
| **Milvus collections are not loaded** | Search returns `collection not loaded`; the service is up | `scripts/load_collections.py` runs each tick |
| **GPU cannot allocate** (page cache holds unified memory) | Small allocations succeed, large models fail | cache drop during bring-up (needs root) |

Each of these fails *silently* — the status table still reads 32/32. That is why they are checked
explicitly below rather than inferred from service health.

## The check

```bash
cd ~/projects/hcls-ai-factory

# wait ~2 min after boot for the first cron tick, then:
./health-monitor.sh status                         # expect: 32/32, All systems operational
.venv/bin/python scripts/load_collections.py       # expect: failed 0, ~3,400 vectors
.venv/bin/python scripts/validate_registry.py --probe   # expect: "Every live endpoint answered"
.venv/bin/python scripts/run_demo.py --check-all   # expect: 17/17
```

Then prove the three silent failures did not happen:

```bash
# 1. the key survived  (expect "Anthropic LLM client initialized", NOT "LLM features disabled")
grep -i 'LLM client init\|LLM features disabled' logs/cart.log | tail -1

# 2. retrieval survived  (expect hits > 0, not "collection not loaded")
curl -s -m 90 -X POST -H 'Content-Type: application/json' \
  -d '{"question":"CD19 CAR-T toxicity"}' http://localhost:8522/search | head -c 200

# 3. the answers are still clinically correct  (expect 6/6)
.venv/bin/python scripts/run_clinical_eval.py
```

## If something is wrong

| Symptom | Cause | Fix |
|---|---|---|
| Services healthy, answers have no citations | key not loaded | check `.env` exists and is readable; `./health-monitor.sh restart <svc>` |
| `collection not loaded` | Milvus restarted before the loader ran | `.venv/bin/python scripts/load_collections.py` |
| A service will not start | its `venv` symlink or the platform `.venv` is missing | `docs/build/VENV_AND_RUNTIME.md` |
| `Skipped tick: previous run still in progress` forever | a started service inherited the lock fd | should be impossible now (`exec 9>&-` in both launch branches) — if it recurs, check `ls -l /proc/<svc-pid>/fd/9` |
| GPU warning in the log | page cache holds unified memory | `sudo sysctl -w vm.drop_caches=3` (the supervisor does this automatically only with root) |

## Recording the result

A3 is met when `status` reads 32/32 and `--check-all` reads 17/17 **without anyone touching
anything** after the boot. If it takes a manual step, A3 is not met — write down which step, because
that step is the actual defect.
