---
name: 11-security-and-secrets
description: >-
  Best-practice standards for Pillar 11 (Security & Secrets) of the HCLS AI Factory. Use when designing,
  building, operating, or reviewing secret handling, input/output governance, disclosure policy, and data
  hygiene. Concrete triggers: adding an env var or credential, wiring governance onto a service, touching
  docker-compose creds, handling patient data, or reviewing anything before it is committed or published.
---

# Pillar 11 — Security & Secrets

Security here is layered and boring on purpose: keep secrets out of git, keep the repo neutral and
free of PHI, gate every service's inputs and outputs, and ship an honest disclosure policy — for an
Apache-2.0, clinical-adjacent platform that is explicitly **not a medical device**.

## In the HCLS AI Factory
- **Secrets live in `.env` (gitignored); only `.env.example` placeholders ship.** Compose credentials
  are env-substituted — required secrets use the fail-fast form `${VAR:?set VAR in .env}` (e.g.
  `HCLS_MINIO_USER`, `HCLS_MINIO_PASSWORD`, `GRAFANA_PASSWORD`, `ANTHROPIC_API_KEY`); safe fallbacks
  use `${VAR:-default}`. No default credentials are baked into `docker-compose.dgx-spark.yml`.
- **Commit-time leak guard** — `scripts/pre-commit-hook.sh` (dependency-free; installed by
  `scripts/install-hooks.sh`) blocks oversized files (>5 MB), secret shapes (`sk-ant-…`, `AKIA…`,
  `-----BEGIN … PRIVATE KEY`, `gh[pousr]_…`), and non-neutral vendor branding. Belt-and-suspenders:
  `.pre-commit-config.yaml` runs **gitleaks**, `detect-private-key`, `check-added-large-files`
  (`--maxkb=5120`), and JSON/YAML checks.
- **Governance gates on every service** — the input-validation gate + output-honesty gate wire in one
  line via `hcls_common.api_gate` (`create_governed_app()` / `install_governance()`).
  `require_valid_input()` rejects bad enums / missing-required with HTTP 422 and clamps out-of-range
  numerics; `honesty_flags()` / `assert_publishable()` run the honesty register in `verify_gate.py`
  (blocks FDA-clearance / diagnostic-certainty / cure / absolute-safety overclaims; warns on
  missing research-use disclaimers).
- **Input sanitization** — `lib/hcls_common/security.py` validates/sanitizes free-text search, gene
  names (HGNC), chromosomes, patient/sample IDs, SMILES, PDB IDs, filenames, and rate-limits;
  `validate_milvus_filter` token-whitelists boolean filter expressions against injection.
- **Dependency + license gate** — `lib/hcls_common/license_gate.py` requires exact-pinned deps
  (21 CFR Part 11 reproducibility) and blocks non-commercial licenses / flags copyleft.
- **Disclosure & data policy** — `SECURITY.md`: private disclosure to adam.m.jones@mac.com,
  5-/15-business-day windows, "**not a medical device**", and the rule that the repo ships **only
  synthetic and public reference data**; real patient data is handled under HIPAA/GDPR off-repo.

## Best-practice standards
- **Env-based secrets only.** Keys/tokens come from env vars or a local `.env`; ship a matching
  `.env.example` placeholder with every new secret. Never hardcode a key, and never a `/home/<user>`
  path or LAN IP in shipping code.
- **No default creds in composes.** Use `${VAR:?…}` for anything that must be set, `${VAR:-default}`
  only for genuinely safe defaults.
- **Govern every service by construction.** Prefer `create_governed_app(...)`; call
  `require_valid_input(capability_id, payload)` in POST handlers and `honesty_flags(text)` /
  `assert_publishable(text, llm=...)` on generated clinical text. Governance is inherited, not optional.
- **Least privilege + sanitize at the edge.** Validate every free-text/gene/variant/SMILES/patient-id
  at the boundary with `hcls_common.security`; whitelist, never blacklist, vector-database filter tokens.
- **Never commit PHI.** Synthetic/public reference data only; the honesty gate labels mock outputs and
  blocks clinical overclaims, but that is a backstop, not permission to load real patient data.
- **Install the local guard once per clone** (`./scripts/install-hooks.sh`); bypass (`--no-verify`)
  only with a stated real reason.
- **Keep deps pinned and permissively licensed** — run `license_gate.audit(...)` in CI intent; a
  non-commercial dependency is a blocker for an Apache-2.0 platform.
- **Report privately.** Route vulnerabilities through `SECURITY.md`, never a public issue.

## Do / Don't
**Do:** put secrets in `.env`, ship `.env.example`; use `${VAR:?…}` in compose; `create_governed_app` +
`require_valid_input` + `honesty_flags`; sanitize inputs via `hcls_common.security`; run the pre-commit
guard; keep only synthetic/public data.
**Don't:** commit real keys or PHI; bake default credentials into a compose; expose an ungoverned route;
blacklist injection tokens instead of whitelisting; advertise the platform as FDA-cleared or a diagnosis;
hardcode machine-specific paths/IPs.

## Wiring it in
```python
# governed by construction (preferred for new services)
from hcls_common.api_gate import create_governed_app, require_valid_input, honesty_flags
app = create_governed_app("cart", capability_id="cart-intelligence-agent")

@app.post("/analyze")
def analyze(payload: dict):
    payload = require_valid_input("cart-intelligence-agent", payload)   # 422 on bad input
    answer = run(payload)
    return {"answer": answer, "honesty": honesty_flags(answer)}          # deterministic overclaim scan
```
```yaml
# docker-compose.dgx-spark.yml — env-substituted creds, no defaults for secrets
environment:
  GF_SECURITY_ADMIN_PASSWORD: ${GRAFANA_PASSWORD:?GRAFANA_PASSWORD must be set in .env}
```
Guards: `./scripts/install-hooks.sh` (native hook) and `pre-commit install` (gitleaks/private-key/large-file).

## Pitfalls (single-box DGX Spark / ARM / this factory)
- **One box = one blast radius.** Every service, the vector database, and the metrics stack share the
  host; a leaked credential or an ungoverned route exposes the whole factory. Bind internal-only
  services to the local network, not `0.0.0.0` on a public interface.
- **The neutrality guard also blocks legitimate prose.** The pre-commit hook rejects the exact vendor
  tokens — write "database"/"vector database" (lowercase), and reach for `--no-verify` only on a real
  false positive (and say why).
- **RunPod burst crosses a network boundary.** Secrets that reach a remote GPU travel the private
  Tailscale mesh; keep them env-injected on the RunPod side too and never in a committed manifest.
- **`nvidia-smi`/host tooling isn't a trust boundary.** Out-of-scope per `SECURITY.md` are issues
  needing pre-existing local/root access and vendored third-party tools (`vendor_proteinmpnn/`,
  `vendor_rfdiffusion/`) — route those upstream.
- **Data lives outside git for a reason.** The 711 GB working tree (VCFs, weights, PDBs) is gitignored;
  a stray `git add -f` of patient-adjacent data is both a security and a size violation.

## Related
- Pillars: 09-ai-orchestration, 10-observability, 12-cost-and-economics
- build-housekeeping-standards
