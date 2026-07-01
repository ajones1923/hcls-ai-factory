---
name: release-and-site-publishing
description: >-
  How to cut a versioned release and publish the public docs site for the HCLS AI Factory. Use when
  tagging a release, writing release notes, updating CHANGELOG/CITATION, or shipping/reconciling the
  docs site at https://hcls-ai-factory.org. Encodes the semver + annotated-tag + GitHub Release
  pattern, the green-CI precondition, CITATION sync, docs-as-code deploy-on-merge to Netlify, and the
  hard separation from the runtime landing dashboard.
---

# Release & Site Publishing — cut the tag, let CI publish the site

Two related but distinct acts: **cutting a release** (a semver tag + GitHub Release, gated on green
CI) and **publishing the public docs site** (docs-as-code → Netlify → the custom domain). Both are
automated by construction — you don't hand-deploy either.

## Cutting a release

Follow the v2.0.0 pattern. A release is a promise the repo can keep, so it ships only from a clean,
gated tree.

1. **CI must be green first.** The merge gate is the root `.github/workflows/ci.yml` — three jobs
   run for real, never `|| true`:
   - **Lint (real-bug rules):** `ruff check --select E9,F82,F811,F706,F707 core lib scripts`
   - **Platform library:** `hcls_common` test matrix (Python 3.11/3.12) via `pytest -q`
   - **Capability registry validation:** `python scripts/validate_registry.py`
   Do not tag a red tree. If a job can't run reliably yet, say so in a comment — never fake a pass.
2. **Pick a semantic version.** MAJOR for breaking changes, MINOR for additive capabilities, PATCH
   for fixes. Current release is v2.0.0.
3. **Create an annotated git tag** (`git tag -a vX.Y.Z -m "..."`) — annotated, not lightweight, so
   the tag carries author, date, and message. Push the tag.
4. **Publish a GitHub Release** off that tag with human release notes: highlights, new
   engines/agents/capabilities, breaking changes + migration, honest labels for anything
   preclinical/roadmap. This is the v2.0.0 pattern.
5. **Update the CHANGELOG** in the same change — a curated, human-readable history, newest first,
   grouped (Added / Changed / Fixed / Removed). The release notes and the CHANGELOG entry agree.
6. **Keep `CITATION.cff` in sync.** Bump its `version:` to match the tag in the *same* change so the
   GitHub "Cite this repository" widget resolves correctly. CITATION currently reads `2.0.0`.

## Publishing the public site (docs-as-code)

The public documentation site is deploy-on-merge, static, and the world's first meeting with the
factory — it must be automated, in-sync, and shareable.

- **Pipeline:** GitHub → static docs build → Netlify → `https://hcls-ai-factory.org` (custom domain
  + TLS). The v1.3.0 release published its Markdown pages this way.
- **Source of truth is the repo.** Docs are Markdown under `docs/` and per-component READMEs; the
  site is *generated*, never hand-authored in a CMS.
- **Deploy on merge to `main` only.** No hand-deploys, no dragging a build folder. Let CI build and
  Netlify publish. Preview locally before merge; keep the build reproducible (pinned generator +
  in-repo config for nav/theme/redirects).
- **Ship OG/share cards + SEO metadata:** an OG image, per-page title + description, canonical URLs,
  and a sitemap, so links unfurl well in Slack/social and pages are indexed.
- **Custom domain + TLS.** Own the apex/`www` DNS and automatic cert renewal; never serve the
  canonical site off a raw platform subdomain.
- **Fail the build on broken internal links / missing nav pages** rather than publishing a dead
  page.

## Keep the site and README in sync (the drift trap)

The most common failure is **drift**: the site and the `README.md` disagree on the
**Eight Engines · Eight Intelligence Agents · One Platform** framing, the counts, or the port map.
Reconcile them in the *same* change: update `README.md` → update `docs/` → the merge republishes the
site. The site and the repo front door must tell one story.

## Don't conflate the docs site with the landing dashboard

Two different "front pages", two different repos:

- **Public docs site** — static, on Netlify, at `https://hcls-ai-factory.org`; this skill.
- **Runtime landing page** — the Flask health dashboard on `:8080` (`landing-page/`), which ships
  from its **own separate repo, `hcls-landing-page-dev`**, intentionally *not* in this monorepo.

They serve different audiences and deploy differently. Don't put the dashboard in this monorepo, and
don't couple the stateless docs site to any factory service or port.

## Do / Don't

**Do:** confirm CI is green before tagging; use semver + an annotated tag + a GitHub Release with
real notes; keep CHANGELOG and `CITATION.cff` in sync in the same change; edit Markdown and let CI
publish the site; ship OG/SEO metadata; verify the custom domain + TLS; reconcile the site with the
README on any framing/port change.

**Don't:** tag a red tree or add `|| true` to a CI step; cut a lightweight tag or skip the release
notes; let CITATION lag the tag (breaks "Cite this repository"); hand-deploy or edit the live site
out-of-band; commit build output; let the site drift from the README; put the `:8080` landing
dashboard in this monorepo.

## Release + publish checklist

1. `ruff check --select E9,F82,F811,F706,F707 core lib scripts` → clean.
2. `( cd lib/hcls_common && pip install -e ".[dev]" -q && pytest -q )` → green.
3. `python scripts/validate_registry.py` → `OK`. (CI enforces all three; confirm green on `main`.)
4. Choose the semver bump; write CHANGELOG entry + GitHub Release notes (honest labels included).
5. Bump `CITATION.cff` `version:` to match — same change.
6. `git tag -a vX.Y.Z -m "..."`; push the tag; publish the GitHub Release.
7. For docs: PR the Markdown, reconcile README ↔ `docs/`, merge to `main`, let CI build + Netlify
   deploy; verify the live site, share-card unfurl, and TLS.

## Related
- `15-github-structure-and-presentation` — front-door docs, real CI, tagged releases, CITATION sync.
- `16-github-netlify-site` — the docs-as-code → Netlify → custom-domain pipeline and its pitfalls.
- `build-housekeeping-standards` — the non-negotiable gate the release rides on; docs-match-code.
