# Genomic Foundation — Advanced Learning Guide

For engineers extending or operating this subject.

**Source:** `core/engines/genomic-foundation` · 15 Python files · 3,264 LOC · 8 test files

## Registered capabilities

| Capability | Type | Status | Endpoint |
|---|---|---|---|
| `genomics-engine` | engine | **live** | localhost:5000 |
| `variant-store` | service | **live** | localhost:8575 |
| `mosaicism-vaf` | stage | **planned** | localhost:8575 |
| `acmg-secondary-findings` | stage | **live** | — |
| `gwas-association` | stage | **live** | — |

**:5000** — a single-service portal, so there is no separate API port

## Principal modules

### `web-portal/app/server.py`

`rate_limit`, `require_api_key`, `save_pipeline_state`, `load_pipeline_state`, `load_config`, `save_config`

- **`rate_limit`** — Apply rate limiting if available, otherwise pass through.
- **`require_api_key`** — Require X-API-Key header for dangerous endpoints.
- **`save_pipeline_state`** — Save pipeline state to disk for persistence across restarts
- **`load_pipeline_state`** — Load pipeline state from disk

### `src/variant_store.py`

`VariantStore`


### `src/acmg_sf.py`

`is_on_panel`, `is_reportable`, `secondary_findings`, `panel_summary`

- **`is_reportable`** — Reportable secondary finding: gene on the ACMG SF panel AND a (likely) pathogenic call.
- **`secondary_findings`** — Filter annotated variants to reportable ACMG SF secondary findings.

### `web-portal/app/security.py`

`add_security_headers`, `generate_csrf_token`, `verify_csrf_token`, `SimpleAuthenticator`, `require_local_access`, `RateLimiter`

- **`add_security_headers`** — Add security headers to response.
- **`generate_csrf_token`** — Generate a CSRF token.
- **`verify_csrf_token`** — Verify a CSRF token using constant-time comparison.
- **`SimpleAuthenticator`** — Simple API key authentication for the portal.

### `web-portal/app/validation.py`

`validate_step_name`, `validate_log_type`, `validate_config_key`, `validate_config_value`, `sanitize_path`, `validate_patient_id`

- **`validate_step_name`** — Validate pipeline step name.
- **`validate_log_type`** — Validate log type name.
- **`validate_config_key`** — Validate configuration key name.
- **`validate_config_value`** — Validate configuration value.


## Dependencies

`Flask-CORS==6.0.2`, `Flask==3.1.2`, `loguru==0.7.3`, `psutil==7.2.1`, `pynvml==13.0.1`, `python-dotenv==1.2.1`

## Running the tests

```bash
.venv/bin/python scripts/run_all_tests.py genomic-foundation
```

Two traps the shared harness handles, which a hand-rolled `pytest` invocation will not:

1. Several subjects ship `src/collections.py`, which **shadows the Python standard library**. Putting
   their `src/` on `PYTHONPATH` kills the interpreter before collection.
2. `structural-biology/vendor_rfdiffusion/` is vendored third-party code needing gated GPU packages
   and is excluded.

## Operational notes

Alignment and variant calling require NVIDIA Parabricks, which is not installed on this box. Results shown today are pre-computed.

Before changing a port, read [`../../build/PORT_MAP.md`](https://github.com/ajones1923/hcls-ai-factory/blob/main/docs/build/PORT_MAP.md). The convention is
enforced by `scripts/validate_registry.py`, which also cross-checks `health-monitor.sh` — a port
change in one place and not the other fails the build.

## Extending it

1. Add or change code under `core/engines/genomic-foundation`.
2. Keep the capability entry in `lib/hcls_common/capabilities.json` truthful — **a `live` status must
   answer a health probe**. Two capabilities were found registered `live` with nothing bound to their
   ports; do not add a third.
3. Run the gate: `ruff`, `pytest lib/hcls_common`, `validate_registry.py`, `run_all_tests.py`.
