---
name: 16-github-netlify-site
description: >-
  Best-practice standards for Pillar 16 (GitHub → Netlify → hcls-ai-factory.org) of the HCLS AI Factory. Use
  when designing, building, operating, or reviewing the public docs site and its deploy pipeline. Concrete
  triggers: editing docs that publish to the site, changing the docs build/nav, updating OG/SEO metadata,
  wiring or debugging the Netlify deploy, or reconciling the site with the README.
---

# Pillar 16 — GitHub → Netlify → hcls-ai-factory.org

The public documentation site is deployed docs-as-code: Markdown in GitHub builds to a static site and
publishes to Netlify at the custom domain https://hcls-ai-factory.org on merge to `main`. The site is how the
world first meets the factory, so it must be automated, in-sync, and shareable.

## In the HCLS AI Factory
- **Pipeline:** GitHub → MkDocs static build → Netlify → https://hcls-ai-factory.org (custom domain + TLS).
  The v1.3.0 release published 142 Markdown pages this way (MkDocs + Netlify).
- **Source of truth is the repo:** the docs live as Markdown under `docs/` and per-component READMEs; the
  site is generated, never hand-authored in a CMS.
- **Deploy on merge:** publish is CI-driven on merge to `main` — no hand-deploys, no dragging a build folder.
- **Discoverability:** an OG/share image (`og-image`) and SEO metadata (title, description, canonical URL)
  so links unfurl well in Slack/social and the site is indexable.
- **Separation of concerns:** the runtime landing page (the Flask health dashboard on `:8080`,
  `landing-page/`) is a *different* thing from the public docs site, and lives in its own repo
  (`hcls-landing-page-dev`) — intentionally **not** in this monorepo. Don't conflate the two.

## Best-practice standards
- Treat docs as code: everything that publishes is Markdown in the repo, reviewed via PR, built by the site
  generator — no out-of-band edits to the live site.
- Deploy on merge to `main` only; keep the build reproducible (pinned generator + config) so a merge always
  produces the same site.
- Own the domain and TLS: custom domain https://hcls-ai-factory.org with automatic certificate renewal;
  never serve the canonical site off a raw platform subdomain.
- Ship social/SEO metadata with the site: an OG/share card, per-page title + description, canonical URLs, and
  a sitemap so links preview well and pages are indexed.
- Keep the site in sync with the README's *Eight Engines · Eight Intelligence Agents · One Platform* framing
  and port map — the site and the repo front door must tell the same story.
- Fail the build on broken internal links / missing nav pages rather than publishing a dead page.
- Keep the docs build reasonably fast and its config in-repo (nav, theme, redirects) so contributors can
  preview locally before merge.

## Do / Don't
**Do:** edit Markdown in the repo and let CI publish; keep OG/SEO metadata current; verify the custom domain
+ TLS; preview the docs build locally before merge; reconcile the site with the README whenever framing or
ports change.
**Don't:** hand-deploy or edit the live site out-of-band; let the site drift from the README (stale engine/
agent counts or ports); publish broken links; put the Flask landing-page dashboard in this monorepo (it's a
separate repo); commit build output — publish it, don't track it.

## Wiring it in
```bash
# docs-as-code publish flow (conceptual)
#   edit docs/*.md  →  PR  →  merge to main  →  CI builds MkDocs  →  Netlify deploys
#                                                                     → https://hcls-ai-factory.org

# preview locally before merge (pin the generator/config in-repo)
mkdocs serve            # local preview; build must be reproducible in CI

# what to keep in the repo, versioned:
#   docs/            markdown source
#   mkdocs config    nav, theme, redirects
#   og-image + SEO   share card + per-page metadata
```
Sync checklist on any framing/port change: update `README.md` → update `docs/` → the merge republishes the
site. The landing dashboard (`:8080`) is out of scope here — it ships from `hcls-landing-page-dev`.

## Pitfalls (single-box DGX Spark / ARM / this factory)
- The docs site is stateless static hosting on Netlify — it does **not** run on the DGX Spark and has none of
  the box's reliability concerns; don't couple it to factory services or ports.
- Easiest failure is drift: the site and the README fall out of sync on engine/agent counts, ports, or the
  hero framing — reconcile them in the same change (ties directly to Pillar 15).
- Don't confuse the two "front pages": the public **docs site** (Netlify, this pillar) vs. the runtime
  **landing dashboard** (Flask `:8080`, separate repo). They serve different audiences and deploy
  differently.
- Custom-domain DNS + TLS is easy to leave half-configured — verify the apex/`www` resolve and the cert
  auto-renews, so the canonical URL never lands on an untrusted or platform-subdomain page.

## Related
- Pillars: 15-github-structure-and-presentation, 14-ease-of-deployment, 03-networking-and-ingress
- build-housekeeping-standards
