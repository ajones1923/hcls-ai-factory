# Deploying the platform to the DGX Spark

This is about putting **repo code onto the running box**. Publishing the *website* is a different
job with a different runbook — `DEPLOY_RUNBOOK.md` (Netlify).

Written 2026-11-04, when `v2.2.0` had been tagged and released for a day and the box was still
serving September code. "Released" and "deployed" had drifted into two different claims, which is
the same defect class this review program exists to close — one level up, in the operations.

---

## What a deploy actually is here

There is **no systemd unit, no container, and no image** for any of the 32 platform services.
Checked 2026-11-04: `systemctl list-units` knows about `docker`, `pmie` and `tailscaled` and
nothing else of ours. Milvus, Grafana, Prometheus and the exporters are Docker Compose; everything
else is a `nohup` child of `health-monitor.sh`, started with `cd <service dir> && ./venv/bin/...`.

Three consequences, and they are the whole reason this document exists:

1. **The deploy unit is a `git checkout` in the live tree.** There is nothing to build or push.
2. **A checkout does not restart anything.** A running service holds its code in memory and keeps
   serving the old version until something restarts it. So the checkout changes what *future*
   starts will run, immediately, and what is being served, not at all.
3. **`cron` runs `health-monitor.sh fix` every five minutes**, and `fix` starts anything DOWN while
   leaving anything HEALTHY alone (`start_service` returns early on `check_service`). So a service
   that is *down* deploys itself on the next tick. One that is *up* does not.

Point 2 and point 3 together are the trap. After a checkout the box is in a mixed state: new code
on disk, old code in memory, and the next crash or reboot silently promotes any service to the new
version at an unpredictable time. **The checkout is the commit point for everything.** Staging only
controls when each service notices. Do not do the checkout and walk away.

## Why there is no blue/green, and what to have instead

Tempting, because `health-monitor.sh` derives every path from `SCRIPT_DIR` — a second complete
checkout would run the fleet from itself with no edits. The venvs do not block it either: 20 of 22
service venvs are **symlinks** to the shared `/home/adam/projects/hcls-ai-factory/.venv`, so a
second tree needs the symlinks recreated, not 8.5 GB of installs. Only
`therapeutic-discovery` (6.6 GB) and `precision-intelligence` (1.8 GB) hold real venvs, and those
can be pointed at the existing ones.

What blocks it is **ports**. Blue and green would both want `:8527`, `:8541`, `:19530`. On one box
with a registry that advertises fixed ports there is no flip that keeps both alive, so you cannot
get zero downtime this way.

You can still get **fast rollback**, which is the half worth having: keep two checkouts and change
which path the supervisor runs from. That change is a line in `crontab`, so it is Adam's to make,
and it is not needed for the deploy below. For a research-use demo platform, rollback measured in
a minute is sufficient and downtime measured in seconds is not a requirement.

---

## Blocker: the live tree is not on `main`

Measured 2026-11-04:

```
$ git -C /home/adam/projects/hcls-ai-factory branch --show-current
provenance/zenodo-metadata
$ git rev-list --left-right --count HEAD...origin/main
1       121
```

The live tree is parked on a feature branch, **one commit ahead and 121 behind**. That commit is
`84f69b4` (`.zenodo.json`), it exists only on this branch and in PR #154, and it is why the box is
on September code: a checkout for that PR was never returned to `main`.

So `git checkout main` today would *discard a commit that is not merged anywhere*. **Merge PR #154
first.** Then the live commit is contained in `main`, the checkout is an ordinary fast-forward, and
nothing is lost. Also merge #162 and #163, which carry the ratchets and the class-B fixes.

Also present: an untracked `explainer-engine/` directory. It survives a checkout untouched. Leave
it alone; do not `git add -A` or `git clean` on the live tree.

---

## Deploy, in waves by blast radius

### Before anything — record the rollback target and the pre-state

```bash
cd /home/adam/projects/hcls-ai-factory
git rev-parse HEAD | tee /home/adam/deploy-rollback-sha.txt   # where to go back to
./health-monitor.sh status --json > /home/adam/deploy-pre-state.json
python scripts/probe_live_endpoints.py                        # registry-side truth
```

The two measurements answer different questions and you want both. `status` probes the
**supervisor's** roster; `probe_live_endpoints.py` probes what the **capability registry**
advertises. On 2026-11-04 they disagreed: `status` said 31/32 healthy while the registry had nine
`live` capabilities nothing served. A deploy that fixes one and not the other has not finished.

### The fast-forward

```bash
git checkout main && git pull --ff-only
grep -c '^    "' health-monitor.sh    # expect 42; September's copy has 33
```

### Wave 0 — the nine agent UIs (no restart, no risk to anything working)

These nine registry entries are marked `live` and nothing listens on them:

```
8126 cardiology  8507 pharmacogenomics  8521 cart       8523 imaging  8526 oncology
8531 autoimmune  8538 clinical-trial    8540 single-cell 8544 rare-disease
```

All nine ship an `app/*_ui.py`, all nine were added to the supervisor roster on 2026-11-03 after
being hand-started on a scratch port and answering `/healthz` with 200, and all nine have a venv
with `streamlit` present — **no installs required**. Nothing is running on those ports, so starting
them cannot break anything that currently works. This is the whole of the deploy's upside with none
of its risk: nine false `live` claims become true.

Do nothing. The next `*/5` tick starts them. Then:

```bash
./health-monitor.sh status | grep -i ui
python scripts/probe_live_endpoints.py && python scripts/check_endpoints.py
```

`check_endpoints.py` does not fail on recovery, so the nine coming up is a silent pass. Update the
baseline deliberately afterward, with a reason, rather than letting the snapshot drift.

### Wave 1 — leaf APIs, one at a time

The agent APIs have no dependents inside the platform. Restart individually and stop at the first
failure, rather than `restart all`:

```bash
for svc in oncology cart imaging autoimmune pharmacogenomics cardiology \
           clinical-trial rare-disease neurology single-cell biomarker; do
  ./health-monitor.sh restart "$svc" || { echo "STOP: $svc"; break; }
done
```

### Wave 2 — shared and upstream, last

`milvus`, `rag-api`, `rag-chat`, `genomics`, `landing`, and `drug-portal` (`:8510`, already
DEGRADED before the deploy — do not read its failure afterward as caused by the deploy). Milvus
does not reload collections after a restart; `health-monitor.sh` has a load step for exactly that,
but verify the vector count against the pre-state before calling it done.

### Rollback, per service

```bash
git checkout "$(cat /home/adam/deploy-rollback-sha.txt)"
./health-monitor.sh restart <svc>
```

Per service, not globally. The failure you are rolling back is almost certainly in one of them, and
reverting the fleet to September to fix one agent trades a small problem for the big one.

---

## What to check before calling it deployed

| Check | Command | Expected |
|---|---|---|
| Supervisor roster is the new one | `grep -c '^    "' health-monitor.sh` | 42 |
| Supervisor view | `./health-monitor.sh status` | 40/41 (`drug-portal` known) |
| Registry view | `python scripts/probe_live_endpoints.py` | 24/24 |
| Ratchet agrees | `python scripts/check_endpoints.py` | exit 0 |
| Corpus did not regress | `python scripts/check_corpus.py` | exit 0 |
| Tree is on a merged commit | `git branch -r --contains HEAD` | includes `origin/main` |

The last row is the one that would have prevented this whole situation.
