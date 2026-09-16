# PRD — Hardening the HCLS AI Factory

**Date:** 2026-09-16 · **Owner:** Adam Jones · **Companion:** [HARDENING_WORKBOOK.md](HARDENING_WORKBOOK.md)
**Status of this document:** repo-only. Not published to hcls-ai-factory.org.

---

## 1. The problem in one paragraph

The factory now **runs**: 32/32 services, 17/17 demonstrations on real input, ~4,000 corpus vectors,
a green merge gate, and a registry whose every `live` claim answers a probe. What has not had a
pass is everything *around* the code — how it is secured, how it is structured, how it is licensed
to the public, and how much of it is duplicated. This document covers that surface. It is not a
rewrite and it is not new capability; it is the difference between a platform that works and one
that can be handed to someone else.

**One finding is urgent and the rest are not.** Clinical decision-support endpoints currently
answer unauthenticated requests from the local network. Everything else here can be scheduled.

---

## 2. Measured state (2026-09-16)

Every figure was produced by running a check.

| Dimension | State |
|---|---|
| Services listening on `0.0.0.0` | **45** |
| Clinical POST from the LAN IP, no credential | **200 OK** — verified against `192.168.68.107:8522/search` |
| `HCLS_API_KEY` | generated, **commented out** in `.env`; gate is fail-closed only once set |
| Entrypoints *able* to enforce auth | 12 of 12 |
| GitHub license detection | **`license=none`** — `LICENSE` is 190 lines vs the canonical 201, APPENDIX absent |
| Branch protection on `main` | required checks **off**, required reviews **off**, `enforce_admins` **false** |
| Open dependabot PRs | **14**, oldest 2026-07-02 (~2.5 months) |
| CI workflow `permissions:` block | **absent** — `GITHUB_TOKEN` takes the repo default |
| Actions pinned | by tag (`@v4`), not SHA |
| Vulnerable packages *in the platform venv* | **2** (`pip` 24.0, `accelerate`) — the other 38 flagged are OS Python |
| Modules shadowing the stdlib | **11** (`src/collections.py`) |
| `rag_engine.py` implementations | **12**, **13,821 LOC**, 15–67% similar to each other |
| Per-service `_LLMClient` copies | **8** (+1 shared) — the `temperature` bug lived in all nine |
| Broken internal doc links | **16** genuinely missing (a further 24 resolve to build-time generated pages) |
| Tracked video in git | **192 MB across 22 files** — ~90% of a 213 MB repo, no LFS |
| Test depth spread | **4 → 1,966** per subject |

---

## 3. What "hardened" means

Five acceptance criteria. Each is a command whose output settles it.

| # | Done when | Command that proves it |
|---|---|---|
| **H1** | No clinical endpoint answers an unauthenticated request from off-box | `curl -X POST http://<lan-ip>:8522/search` returns **401** |
| **H2** | GitHub shows the repo as Apache-2.0 | `gh repo view --json licenseInfo` returns `Apache-2.0` |
| **H3** | A red CI run cannot be merged to `main` | required status checks listed in branch protection |
| **H4** | One bug fixed in the RAG path is fixed everywhere | the 3-fault chain has one home, not twelve |
| **H5** | A clean clone is under 50 MB | `git clone --depth 1` transfers < 50 MB |

---

## 4. Requirements

**P0** blocks release · **P1** blocks the phase · **P2** quality.

### 4.1 Security

| # | Requirement | Pri | Done when |
|---|---|---|---|
| S1 | Uncomment `HCLS_API_KEY`; verify fail-closed on every one of the 12 entrypoints | **P0** | 401 without a key, 401 with a wrong key, 200 with the right one |
| S2 | Bind services to `127.0.0.1`, or block 8500–8600 at the firewall, or both | **P0** | LAN probe refused |
| S3 | Add `permissions:` to every workflow (least privilege, `contents: read`) | P1 | no job holds write scope it does not use |
| S4 | Pin actions by commit SHA | P2 | `uses:` lines carry 40-char SHAs |
| S5 | Upgrade `pip` and `accelerate` in the platform venv | P1 | `pip-audit --local` clean |
| S6 | Decide whether a LAN-exposed demo posture is intended at all | **P0** | recorded in §7 |

> **Why S1/S2 are P0 and nothing else is.** The header already reports only the gates that ran,
> and the endpoints *can* enforce auth — the gate is simply off. On a trusted single box that was
> a defensible demo posture. It stopped being defensible when the box acquired an API key, a
> seeded clinical corpus, and a LAN address: an unauthenticated caller can now obtain generated
> clinical prose with dosing and management guidance, and spend the key doing it.

### 4.2 Licensing and GitHub

| # | Requirement | Pri | Done when |
|---|---|---|---|
| G1 | Replace `LICENSE` with the canonical Apache-2.0 text including the APPENDIX | **P0** | GitHub shows Apache-2.0 |
| G2 | Triage the 14 dependabot PRs; merge the safe, close the rest with a reason | P1 | open count < 5 |
| G3 | Turn on branch protection: require the 6 CI checks, require a PR | P1 | a red run cannot merge |
| G4 | Add `CODEOWNERS` and a PR template | P2 | both present |
| G5 | Move tracked video to LFS or out of the repo | P2 | H5 passes |

> **G1 is P0 because the whole positioning depends on it.** The project's pitch is "Apache-2.0 on
> purpose — take it, build a business on it." GitHub currently displays no license, and corporate
> scanners read an unlicensed public repo as *do not use*. The text is present but non-canonical,
> so `licensee` cannot classify it. This is a five-minute fix guarding the central claim.
> **The `pymilvus 2.x → 3.0.0` PRs are a major version bump** — 2.6.8 is what is installed and
> working, and the agents' retrieval was only just repaired. Do not merge those without running
> the clinical eval after.

### 4.3 Structure

| # | Requirement | Pri | Done when |
|---|---|---|---|
| C1 | Rename the 11 `src/collections.py` → `vector_collections.py` | P1 | bare `pytest` works without the harness workaround |
| C2 | Extract the shared RAG path into `hcls_common` | P1 | the search/flatten/route logic has one home |
| C3 | Extract the per-service `_LLMClient` into `hcls_common` | P1 | one place to fix a model or parameter change |
| C4 | Raise test depth where output is clinical (R25) | P2 | floor agreed and met |

> **C2 and C3 are the ones that pay.** This is not tidiness: in the last week a single
> `_CollectionManager`/`param=`/result-shape chain broke retrieval in **five** services
> identically, and a `temperature` change broke synthesis in **eight**, because each has its own
> copy. Every future model or client change costs 8–12 edits and will be applied to some and
> missed on others — which is exactly how both bugs shipped. `ingest_persist.py` is the pattern
> to follow: one implementation, the traps documented in it.

### 4.4 Documentation

| # | Requirement | Pri | Done when |
|---|---|---|---|
| D1 | Fix or remove the 16 broken internal links | P2 | link check clean |
| D2 | Mark superseded guides as superseded | P2 | `HCLS_AI_FACTORY_DEMO_GUIDE.md` (pre-convention ports) carries a banner |
| D3 | Keep the 509-file doc surface, but publish only what is current | P2 | site nav unchanged; stale files labelled |

---

## 5. Phases

```
Phase 0  Close the door          S1 S2 S6 G1        ~1 hour      ← do this first
   │     gate: LAN probe refused; GitHub shows Apache-2.0
   ▼
Phase 1  Supply chain            S3 S4 S5 G2 G3     ~half a day
   │     gate: pip-audit --local clean; red CI cannot merge; <5 open PRs
   ▼
Phase 2  De-duplicate            C2 C3 C1           ~2-3 days
   │     gate: one RAG path, one LLM client; bare pytest works
   ▼
Phase 3  Polish                  G4 G5 D1 D2 D3 C4  ongoing
         gate: clean clone <50 MB; link check clean
```

Phases 0 and 1 need no code restructuring. Phase 2 is the only one that touches how the agents are
built, and it should be done **after** the clinical eval is expanded — that eval is what will tell
you the refactor did not change any answer.

---

## 6. Non-goals

- **Not** a rewrite of the engine/agent split, the registry, or the port convention.
- **Not** new capability. The gap here is around the code, not in it.
- **Not** multi-tenant auth or RBAC. One box, one operator; a single API key is the right size.
- **Not** deleting the unpublished doc surface — provenance has value; label it instead.
- **Not** raising test counts for their own sake. Depth where output is clinical (R25), not evenness.

---

## 7. Open decisions

| # | Decision | Blocks | Why it matters |
|---|---|---|---|
| **H-D1** | **Is LAN exposure intended?** Bind to localhost, firewall the range, or accept it and require the API key. | S1, S2 | Determines whether this is a demo box or a service. The current state — exposed *and* unauthenticated — is the one option nobody would choose deliberately. |
| **H-D2** | **AlphaMissense is CC BY-NC-SA** (non-commercial) across 132 tracked files, while the platform is Apache-2.0 and the site welcomes commercial use. *(Carried from the completion PRD as D4; still open.)* | G1, site copy | A licence conflict at the centre of the value proposition. Decide on your terms rather than discover it. |
| **H-D3** | **Merge or close the pymilvus 3.0 PRs?** 2.6.8 works; 3.0 is a major bump. | G2 | Retrieval was broken in five services a week ago. A silent client-API change would be very hard to attribute. |
| **H-D4** | **Video in git, or LFS, or out?** 192 MB of 213 MB. | G5, H5 | Every clone and every CI run pays for it. The site needs the files; git may not. |
| **H-D5** | **What is the test-depth floor for clinical output?** Current spread 4 → 1,966. | C4 | Without a number, R25 stays open forever. |

---

## 8. Risks

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| Turning on the API key breaks the demos | **high** | demos fail mid-presentation | `run_demo.py` and the eval must send the key; test both before relying on it |
| The RAG de-duplication changes an answer | medium | silent clinical regression | expand the clinical eval **first**; it is the only detector |
| pymilvus 3.0 changes the client API again | medium | retrieval breaks in all agents | H-D3; run the eval after any merge |
| Canonical LICENSE replaces a modified one | low | losing an intentional edit | diff before replacing; the only delta should be the APPENDIX |
| Firewalling breaks the Netlify/Caddy path | low | site or portal unreachable | S2 covers 8500–8600 only; 80/443 untouched |

---

## 9. What this buys

H1–H5 turn a platform that runs into one that can be handed over: a stranger can clone it in under
a minute, GitHub tells them they are allowed to use it, a red build cannot reach `main`, an
unauthenticated stranger on the network cannot ask it for dosing advice, and a bug fixed in the
retrieval path is fixed everywhere rather than in one of twelve copies.

---

*Work the [Workbook](HARDENING_WORKBOOK.md) in order. It carries the exact commands and the
verification after each step.*
