# GitHub to Netlify Deployment Flow

This document explains how content moves from your local development environment through GitHub to the live Netlify-hosted website.

---

## Overview

```
┌─────────────┐      ┌─────────────┐      ┌─────────────┐      ┌─────────────┐
│   Local     │      │   GitHub    │      │   Netlify   │      │    Live     │
│   Machine   │ ───▶ │    Repo     │ ───▶ │   Build     │ ───▶ │   Website   │
│             │ git  │             │ hook │             │ CDN  │             │
│ docs/*.md   │ push │ main branch │      │ mkdocs      │      │ hcls-ai-    │
│ mkdocs.yml  │      │             │      │ build       │      │ factory.org │
└─────────────┘      └─────────────┘      └─────────────┘      └─────────────┘
```

---

## Step-by-Step Flow

### 1. Local Development

You edit files on your local machine:

| File/Folder | Purpose |
|-------------|---------|
| `docs/*.md` | Markdown content pages |
| `docs/overrides/` | Custom HTML templates (home.html, 404.html) |
| `docs/stylesheets/extra.css` | Custom styling |
| `docs/javascripts/extra.js` | Custom JavaScript |
| `docs/assets/` | Images, logos, diagrams |
| `mkdocs.yml` | Site configuration (nav, theme, plugins) |
| `netlify.toml` | Netlify build and hosting configuration |

**Local preview:**
```bash
mkdocs serve
# Opens at http://127.0.0.1:8000
```

---

### 2. Git Push to GitHub

When you push changes:

```bash
git add .
git commit -m "Update content"
git push origin main
```

Your changes go to:
- **Repository:** `github.com/NVIDIA/hcls-ai-factory`
- **Branch:** `main`

---

### 3. Netlify Webhook Trigger

GitHub and Netlify are connected via a **webhook**:

1. GitHub sends a POST request to Netlify when `main` branch updates
2. Netlify receives notification of new commit `464237c...`
3. Netlify queues a new build

**How this was set up:**
- In Netlify dashboard → Site settings → Build & deploy → Continuous deployment
- Linked to GitHub repository with OAuth authorization
- Netlify automatically created the webhook in your GitHub repo

You can see the webhook at:
`GitHub repo → Settings → Webhooks → netlify.com webhook`

---

### 4. Netlify Build Process

Netlify spins up a build container and executes the commands in `netlify.toml`:

```toml
[build]
  command = "pip install 'mkdocs-material>=9.5' && mkdocs build"
  publish = "site"

[build.environment]
  PYTHON_VERSION = "3.12"
```

**What happens:**

1. **Clone** — Netlify clones your repo at the new commit
2. **Environment** — Sets up Python 3.12
3. **Install** — Runs `pip install 'mkdocs-material>=9.5'`
4. **Build** — Runs `mkdocs build` which:
   - Reads `mkdocs.yml` configuration
   - Processes all `docs/*.md` files
   - Applies the Material theme
   - Applies custom templates from `docs/overrides/`
   - Copies static assets
   - Generates static HTML in the `site/` directory
5. **Publish** — Takes the `site/` folder as the deployment artifact

**Build output structure:**
```
site/
├── index.html          ← From docs/index.md + home.html template
├── architecture/
│   └── index.html      ← From docs/architecture.md
├── quickstart/
│   └── index.html      ← From docs/quickstart.md
├── stage-1-genomics/
│   └── index.html
├── stage-2-rag/
│   └── index.html
├── stage-3-drug-discovery/
│   └── index.html
├── assets/
│   ├── stylesheets/
│   │   └── extra.css
│   ├── javascripts/
│   │   └── extra.js
│   └── logo-header.png
├── search/
│   └── search_index.json
└── 404.html
```

---

### 5. Netlify CDN Deployment

After build succeeds:

1. **Atomic deploy** — New version staged without affecting live site
2. **Asset hashing** — Static assets get content-based hashes for cache busting
3. **CDN distribution** — Files pushed to Netlify's global edge network
4. **DNS cutover** — Traffic routes to new version instantly

**Your CDN edges:**
- North America, Europe, Asia, Australia, South America
- Users get content from nearest edge location
- Typical latency: 20-50ms globally

---

### 6. HTTP Headers Applied

Netlify applies headers from `netlify.toml`:

```toml
# Static assets cached for 1 year (immutable)
[[headers]]
  for = "/assets/*"
  [headers.values]
    Cache-Control = "public, max-age=31536000, immutable"

# Search index cached for 1 day
[[headers]]
  for = "/search/*"
  [headers.values]
    Cache-Control = "public, max-age=86400"

# HTML pages always revalidated
[[headers]]
  for = "/*.html"
  [headers.values]
    Cache-Control = "public, max-age=0, must-revalidate"

# Security headers on all responses
[[headers]]
  for = "/*"
  [headers.values]
    X-Frame-Options = "SAMEORIGIN"
    X-Content-Type-Options = "nosniff"
    X-XSS-Protection = "1; mode=block"
    Referrer-Policy = "strict-origin-when-cross-origin"
    Permissions-Policy = "camera=(), microphone=(), geolocation=()"
    Strict-Transport-Security = "max-age=31536000; includeSubDomains"
```

**What these do:**

| Header | Purpose |
|--------|---------|
| `Cache-Control: immutable` | Browser never re-requests assets (fast loads) |
| `Cache-Control: must-revalidate` | HTML always checked for updates |
| `X-Frame-Options` | Prevents clickjacking attacks |
| `X-Content-Type-Options` | Prevents MIME sniffing |
| `Strict-Transport-Security` | Forces HTTPS for 1 year |
| `Permissions-Policy` | Disables camera/mic/location APIs |

---

### 7. Custom Domain Routing

Your domain `hcls-ai-factory.org` is configured:

1. **DNS records** point to Netlify:
   - `A` record → Netlify load balancer IP
   - `CNAME` for `www` → Netlify subdomain
2. **SSL certificate** — Netlify auto-provisions via Let's Encrypt
3. **HTTPS redirect** — HTTP requests automatically redirect to HTTPS

---

## Complete Request Flow

When someone visits `https://hcls-ai-factory.org/quickstart/`:

```
User Browser
     │
     ▼
┌─────────────────────────────────────────────────────────────────┐
│ 1. DNS Lookup: hcls-ai-factory.org → Netlify edge IP            │
└─────────────────────────────────────────────────────────────────┘
     │
     ▼
┌─────────────────────────────────────────────────────────────────┐
│ 2. TLS Handshake: Verify SSL certificate                        │
└─────────────────────────────────────────────────────────────────┘
     │
     ▼
┌─────────────────────────────────────────────────────────────────┐
│ 3. Edge Cache Check: Is /quickstart/index.html cached?          │
│    • HIT → Return cached response (< 50ms)                      │
│    • MISS → Fetch from origin, cache, return                    │
└─────────────────────────────────────────────────────────────────┘
     │
     ▼
┌─────────────────────────────────────────────────────────────────┐
│ 4. Response with headers:                                       │
│    • Content-Type: text/html                                    │
│    • Cache-Control: public, max-age=0, must-revalidate          │
│    • Strict-Transport-Security: max-age=31536000                │
└─────────────────────────────────────────────────────────────────┘
     │
     ▼
┌─────────────────────────────────────────────────────────────────┐
│ 5. Browser renders HTML, requests assets:                       │
│    • /assets/stylesheets/extra.css (cached 1 year)              │
│    • /assets/logo-header.png (cached 1 year)                    │
│    • /assets/javascripts/extra.js (cached 1 year)               │
└─────────────────────────────────────────────────────────────────┘
```

---

## Build Logs

You can view build logs at:
- **Netlify Dashboard** → Deploys → Click any deploy → View build log

Example successful build log:
```
12:34:56 PM: Build ready to start
12:34:57 PM: Cloning repository...
12:34:58 PM: Preparing Git Reference refs/heads/main
12:35:00 PM: Installing pip dependencies
12:35:15 PM: Running build command: mkdocs build
12:35:16 PM: INFO - Building documentation...
12:35:17 PM: INFO - Cleaning site directory
12:35:20 PM: INFO - Documentation built in 3.42 seconds
12:35:21 PM: Uploading artifacts
12:35:25 PM: Deploy live at https://hcls-ai-factory.org
```

---

## Triggering Rebuilds

### Automatic (default)
Every `git push` to `main` triggers a build.

### Manual
- Netlify Dashboard → Deploys → "Trigger deploy" → "Deploy site"

### Clear cache and rebuild
- Netlify Dashboard → Deploys → "Trigger deploy" → "Clear cache and deploy site"

Use this if you suspect caching issues.

---

## Rollback

If a deploy breaks something:

1. Go to Netlify Dashboard → Deploys
2. Find the last working deploy
3. Click "..." menu → "Publish deploy"

The previous version goes live instantly (no rebuild needed).

---

## Summary

| Stage | What Happens | Time |
|-------|--------------|------|
| Git push | Code uploaded to GitHub | 1-2 sec |
| Webhook | GitHub notifies Netlify | < 1 sec |
| Build queue | Netlify schedules build | 0-30 sec |
| Build | mkdocs generates static site | 15-30 sec |
| Deploy | Files pushed to CDN edges | 5-10 sec |
| **Total** | Push to live | **~1 minute** |

Your content flows from local Markdown files → GitHub repository → Netlify build server → global CDN → user's browser, all automatically on every push.
