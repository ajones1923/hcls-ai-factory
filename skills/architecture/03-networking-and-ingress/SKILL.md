---
name: 03-networking-and-ingress
description: >-
  Best-practice standards for Pillar 03 (Networking & Ingress) of the HCLS AI Factory. Use when
  designing, building, operating, or reviewing how traffic reaches a service — the TLS edge, the port
  map, the internal compose network, and CORS. Concrete triggers: adding a service port, editing the
  Caddyfile, setting CORS origins, or wiring one service to call another.
---

# Pillar 03 — Networking & Ingress

How requests reach the factory and how services find each other: a single TLS edge (Caddy) in front of
the landing page, a documented port map, an internal docker-compose network for service-to-service
calls, and env-driven CORS — with no hardcoded LAN IPs in shipping code.

## In the HCLS AI Factory
- **TLS edge — Caddy (`Caddyfile`).** HTTPS terminates at Caddy and reverse-proxies to the landing page
  (`reverse_proxy localhost:8080`), which fronts the whole factory. `tls internal` uses Caddy's managed
  internal CA; `:80` redirects to HTTPS. The disease-program (TSC) engine gets its own dedicated HTTPS
  port block (`:8543` → API `:8560`) rather than a subpath, so the FastAPI API proxies cleanly.
- **The port map (canonical, from `README.md`):** 8080 landing · 5000 genomics · 5001 Precision
  Intelligence (RAG) · 8505 Therapeutic Discovery · 8510 portal · 19530 vector DB (Milvus) · agents in
  the 85xx range (8521–8544) · protein/structural services **8570–8578** (8570 ESMFold, 8571 sequence
  search, 8572 ADMET, 8573 single-cell, 8574 molecule gen, 8575 variant store, 8576 developability) ·
  3000 Grafana · 9099 Prometheus.
- **Internal network:** compose services address each other by service DNS name, not IP — the landing
  page reaches agents via `http://precision-oncology-agent:8000`, `http://cardiology-...:8126`, etc.
  (see the `landing-page` env in `docker-compose.dgx-spark.yml`). Only the edge and the mapped ports
  are externally reachable.
- **CORS / origins come from env** (`lib/hcls_common/config.py`), not literals; the LLM/embedding
  provider allow-lists live there too.

## Best-practice standards
- **One TLS edge, one front door.** Terminate TLS at Caddy and proxy to the landing page; don't expose
  raw service ports to the outside when a proxied route will do.
- **Every new service gets a port in the documented map**, in its range (agents 85xx, protein 857x), and
  the map in `README.md` is updated in the same change — no undocumented ports.
- **Talk to peers by service name on the internal network**, never by host IP; only the edge knows the box.
- **Never hardcode a LAN IP in shipping code.** The one place a concrete IP is legitimate is the Caddy
  edge config (it must issue a cert for the address you reach it by); application code derives origins
  from env (`$HCLS_ROOT`-style config), not `192.168.x.x`.
- **CORS origins are an explicit env allow-list** — no `*` on anything that returns clinical output.
- **Bind services to the compose network, publish only what needs publishing** (UI + API ports); keep
  etcd/MinIO internal-only unless a console is deliberately exposed.
- **`:80` always redirects to HTTPS** so plain-http hits still land.

## Do / Don't
**Do:** front everything with Caddy; give TSC/heavy APIs a dedicated HTTPS port (not a fragile subpath);
reference peers as `http://<service-name>:<port>`; pull CORS origins and provider allow-lists from env;
update the port map when you add a port.
**Don't:** hardcode `192.168.68.107` in Python/JS; set `Access-Control-Allow-Origin: *` on agent APIs;
publish internal infra ports (etcd/MinIO) to the LAN by default; invent an off-map port; assume
service-to-service calls should go through the public edge.

## Wiring it in
```caddyfile
# Caddyfile — TLS edge fronting the landing page (the landing page fans out to every service)
192.168.68.107, localhost, 127.0.0.1 {
    tls internal
    reverse_proxy localhost:8080
}
:80 { redir https://{host}{uri} }
```
```yaml
# docker-compose.dgx-spark.yml — service-to-service by DNS name, not IP
landing-page:
  environment:
    ONCOLOGY_URL: http://precision-oncology-agent:8000
    CARDIOLOGY_URL: http://cardiology-intelligence-agent:8126
```

## Pitfalls (single-box DGX Spark / ARM / this factory)
- **`tls internal` = self-signed internal CA** — clients see an untrusted cert until they trust Caddy's
  root; that's expected on a single LAN box, not a bug.
- **The repo cert (`certs/cert.pem`) is unused by default** — its `key.pem` is mode-600/owned by the
  user and unreadable by the caddy service; copy both to a caddy-readable path to switch off `tls internal`.
- **Subpath proxying breaks FastAPI apps** that assume root — give real APIs a dedicated port block
  (as TSC does on `:8543`).
- **Port collisions on one box are easy** — the 85xx agent range and 857x protein range are tight; check
  the map before claiming a port.
- **Live agents run as native host processes** (not the compose containers) — service-name DNS only
  resolves inside the compose network, so host-mode deployments must use localhost + real ports.

## Related
- Pillars: 04-containers-and-orchestration-runtime, 11-security-and-secrets, 10-observability
- build-housekeeping-standards
