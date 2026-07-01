# TSC Intelligence Engine — PHI-Readiness Gate

**Status (2026-06-29): SYNTHETIC-ONLY. No PHI. Not FDA-cleared. SaMD posture undetermined.**

The engine today runs entirely on a deterministic synthetic 50-patient cohort, every output is
watermarked `SYNTHETIC`, and HIPAA/FDA are explicitly out of scope. This document is the **hard gate**:
*every* item below must be satisfied before any real patient data (PHI) touches this engine. It is a
checklist, not an authorization — clearing it is necessary, not sufficient; institutional sign-off
(IRB/Privacy/Security) is the final authority.

## Control gate

| # | Control | Status | Detail / what remains |
|---|---|---|---|
| 1 | **Request authentication** on API + surfaces | ⚙️ built, OFF by default | `X-API-Key` middleware (P1-7). Fail-closed when `TSC_API_KEY` is set; **must be set** before PHI. Surfaces (Streamlit) need their own auth/proxy too. |
| 2 | **No default credentials** | ✅ done | `tsc:tsc`/`minioadmin` removed → env-driven fail-closed `${VAR:?}` (P0-3); real secrets in gitignored `.env` (mode 600). |
| 3 | **Network exposure** — data tier off `0.0.0.0` | ✅ done | postgres/redis/minio/milvus/ollama bound `127.0.0.1` (P0-6). The API/surfaces front doors still need a TLS reverse proxy (P1-7 TLS, gated) before any non-LAN exposure. |
| 4 | **Encryption at rest** | ❌ required | NVMe is unencrypted. LUKS on the data disk (and backup drive) is a prerequisite for PHI — needs a reformat/reboot (operator action, P1-2). |
| 5 | **PHI-adjacent inference stays local** | ⚙️ configurable | Route any PHI-bearing text through the local Ollama tier, **never** the frontier API. The model router supports local fallback; enforce it (no `TSC_ANTHROPIC_API_KEY` use on PHI prompts) by policy + config. |
| 6 | **Audit logging** of access + decisions | ❌ required | Add an access/audit log (who read which patient, when) before multi-user PHI use. |
| 7 | **Durable, backed-up event store** | ✅ done | Event-sourced SQLite/Postgres persists across restarts (P3-1); named-volume backup + DR drill (P0-2/P1-6). For PHI use the Postgres backend (`TSC_POSTGRES_DSN`) on encrypted storage. |
| 8 | **Honesty register intact** | ✅ verified | fail-to-pending, calibrated 50/90% intervals, `SYNTHETIC` watermark, offline-mode self-labeling, computed-from-data ACMG (validator-authoritative). Audited 2026-06-29 (P3-6) — no regression. |
| 9 | **SaMD / HIPAA / IRB scope resolved** | ⛔ institutional | Determine SaMD classification, execute BAA/DUA, obtain IRB approval. **Adam + institution only.** |
| 10 | **Clinician review of generated content** (AC-8) | ⛔ Adam | Every draft is decision-support requiring clinician sign-off (the variant draft already has a `sign-off` gate). |

Legend: ✅ done · ⚙️ built/configurable, must be activated · ❌ required, not yet built · ⛔ institutional/operator.

## The one-line summary
**Synthetic demo today; before PHI you must (at minimum):** set `TSC_API_KEY`, enable disk encryption,
pin inference to local for PHI prompts, add audit logging, move to the Postgres backend on encrypted
storage, and clear SaMD/HIPAA/IRB. Items 2, 3, 7, 8 are already done; 1 and 5 are built and need
activating; 4, 6 are real build/ops work; 9, 10 are institutional.

*Companion: see `../../../ai-assistant/AI_FACTORY_UPGRADE_PLAN.md` §P3-5 / SEC-17.*
