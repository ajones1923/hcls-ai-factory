# Deploy Runbook — hcls-ai-factory.org

**Target:** `ephemeral-biscuit-62e53b.netlify.app` → `https://hcls-ai-factory.org`

The site is **built and verified**; only publication is outstanding.

```
mkdocs build --strict     PASS
126 HTML pages            includes 69 subject pages + 11 build docs
22 videos embedded        15 re-narrated in Adam's clone (2026-08-16)
```

Build with the isolated docs toolchain (**not** the platform venv — see the warning at the end):

```bash
.venv-docs/bin/python -m mkdocs build --strict
```

---

## Why this is blocked

There are **no Netlify credentials on this machine**. `~/.config/netlify/config.json` contains only
`telemetryDisabled` and a `cliId` — no auth token, no user, no linked site. `netlify status` reports
*Not logged in*, and `netlify login` is an interactive browser flow.

Knowing the site URL does not grant access. One of the two paths below is required.

---

## Path A — direct deploy (recommended)

Publishes the built artifact. **Adds nothing to git history.**

1. Create a personal access token: Netlify → **User settings → Applications → New access token**
2. Then:

```bash
cd /home/adam/projects/hcls-ai-factory
export NETLIFY_AUTH_TOKEN=<token>
.venv-docs/bin/python -m mkdocs build --strict
netlify deploy --prod --dir=site --site=ephemeral-biscuit-62e53b
```

Or interactively, once:

```bash
netlify login && netlify link --name ephemeral-biscuit-62e53b
netlify deploy --prod --dir=site
```

## Path B — git push (the historical path)

`netlify.toml` builds from GitHub on merge to `main`. This is how the site has been updated before.

**Cost:** the change carries **~121 MB**, almost all of it the 15 re-narrated videos. Those files are
already tracked (44 video objects in git), so replacing them writes ~121 MB of **new blobs that stay
in history permanently**. Every future re-narration repeats that. Path A avoids it entirely.

```bash
git checkout main && git merge cardiac-demo-and-milvus-fixes   # currently 6 commits ahead
git add -A && git commit -m "docs: 85 subject assets, build audit, re-narrated site videos"
git push origin main          # Netlify builds automatically
```

⚠️ **Before pushing, decide the sync scope.** The working tree contains `.env`-adjacent files and
1.4 GB of genomic data under `hcls-ai-factory-core-data/`. `.gitignore` covers `.env`, `.env.*` and
`*.env`, but a blind `git add -A` should still be reviewed:

```bash
git status --porcelain | wc -l          # 47 files
git add -An | grep -iE "env|key|secret|token|\.pem|core-data"   # dry run — must be empty
```

---

## What is still not re-narrated

7 of 22 videos remain in the George placeholder voice:

| Video | Why |
|---|---|
| `cardiology-intelligence-agent` | narration data exists, **beat images missing** |
| `genomics-engine`, `structural-biology-engine` | built by a different script |
| `frontier-bionemo-chai`, `frontier-chai-bionemo` | different builder |
| `tsc-two-journeys`, `tsc-usecase` | different builder |

Re-narrate the 15 that are covered with:

```bash
cd tmp/videos
set -a; . /home/adam/projects/services/demo-pipeline/eleven-labs/.env; set +a
export ELEVEN_VOICE_ID=2IhYZeguWsdW83BKiLE4      # asserted; without it you get George
python rebuild_site_videos_adam.py --list
```

Originals (George) are backed up with checksums at
`/home/adam/backups/hcls-ai-factory/2026-08-16-site-videos-george/`.
**One exception:** `singlecell-compute.mp4` in that backup is the *new* Adam version — its George
original was overwritten before the backup ran.

---

## ⚠️ Build with `.venv-docs`, never `.venv`

Installing `requirements-docs.txt` into the platform venv drags in `zarr` and `fast-array-utils`,
whose pytest plugins import an application Settings model at startup and **kill test collection for
three subjects** (8,402 → 8,029 passing). `scripts/run_all_tests.py` now disables those two plugins
specifically — not all autoloading, because the suites genuinely need `pytest-asyncio`.

Keep the toolchains separate:

| venv | Purpose |
|---|---|
| `.venv` | platform: `hcls_common`, tests, demos, registry |
| `.venv-docs` | mkdocs and its plugins only |
