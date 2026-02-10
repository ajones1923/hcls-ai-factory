# Deployment Guide

## Default: Local Development

All services run on localhost without TLS. Access the landing page at `http://localhost:8080`.

## Production: HTTPS with Caddy

[Caddy](https://caddyserver.com/) provides automatic HTTPS via Let's Encrypt.

### Quick Start

```bash
cd deploy/caddy

# Localhost (self-signed certificate)
docker compose up -d

# Production (auto-provisions Let's Encrypt cert)
DOMAIN=yourdomain.com docker compose up -d
```

### Requirements

- Docker and Docker Compose installed
- Port 80 and 443 available
- For real TLS: a domain name pointing to your server's public IP

### Architecture

```
Client (HTTPS:443) → Caddy → localhost:8080 (Landing)
                            → localhost:5000 (Genomics)
                            → localhost:5001 (RAG/Chat)
                            → localhost:8501 (Chat UI)
                            → localhost:8505 (Drug Discovery)
                            → localhost:3000 (Grafana)
```

Caddy handles TLS termination, certificate renewal, and HTTP→HTTPS redirect automatically.
