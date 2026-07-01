---
name: regulatory-compliance-posture
description: >-
  The regulatory and compliance posture of the HCLS AI Factory — what it is, what it is not, and the
  conservative framing to hold in any doc, README, demo, or outreach that touches regulatory,
  clinical-use, or data-handling questions. Consult whenever a claim, disclaimer, or design choice
  brushes against medical-device status, patient data, reproducibility, licensing, or clinician
  review. A regulator or clinician reading the platform should find it honest and conservative.
---

# Regulatory & Compliance Posture — conservative by design

The HCLS AI Factory is a **research and reference platform, not a medical device.** Everything below
keeps that framing defensible: it is decision support with a clinician in the loop, reproducible in
the spirit of 21 CFR Part 11, ships only synthetic and public data, and pins permissively licensed
dependencies. The rule of thumb: **when a reader has to guess how conservative we are, we have not
been conservative enough.**

## What this platform is — and is not

- **Not a medical device.** It is **not FDA-cleared, CE-marked, or 510(k) cleared**, and must not be
  used to make patient-care decisions (`SECURITY.md`, "Important: not a medical device"). Never
  advertise or imply clearance. The honesty register in `verify_gate.py` `block`s exactly these
  strings (`FDA[- ]?(approved|cleared)`, `CE[- ]?marked`, `510(k) cleared`).
- **Decision support, not diagnosis or prescribing.** Outputs inform a qualified clinician; they do
  not diagnose, prescribe, or direct a patient. Diagnostic-certainty and treatment-directive
  language is blocked/warned by the gate. (See `clinical-claim-honesty`.)
- **Research and demonstration outputs only.** The governance layer labels mock outputs and blocks
  clinical overclaims *by design*, but that is a backstop — outputs are for research and
  demonstration, not clinical use.

## 21 CFR Part 11-style reproducibility (the audit spine)

Regulatory-grade reproducibility is built in, not bolted on:

- **Per-run manifests.** `lib/hcls_common/reproducibility.py` (`ReproducibilityManifest`) captures,
  for every pipeline run: software/library versions, model weights + embeddings, reference-data
  versions, hardware (incl. detected GPU/CUDA), input/output files, timing, and the stages/agents
  invoked — serialized to JSON with a stable `run_id`. This is the record that makes a run defensible.
- **Audit trail + run lineage.** `lib/hcls_common/mlops.py` (`MLOpsStore`) tracks experiments,
  params, metrics, and run status/lineage on the single box (no external service required) — the
  "who ran what, with which inputs, when, and how it turned out" record.
- **Attach a manifest to anything reproducible.** New capabilities emit a manifest and register a
  run; don't reinvent tracking — use `hcls_common`. Reproducibility is a first-class output, not an
  afterthought.

## Data governance — synthetic and public only, in-repo

- **The repository ships ONLY synthetic and public reference data.** No PHI, ever, in tracked files
  (also a size + neutrality rule — the working tree stays local; only code + docs publish).
- **Real patient data is the user's responsibility, off-repo.** Anyone applying this to real data
  owns **HIPAA / GDPR** compliance, **de-identification**, **IRB / consent**, and **data
  provenance** — under applicable law, outside this repository. Say this plainly; do not imply the
  platform handles those obligations for the user.
- **Handle data at the edge.** Patient/sample IDs and free-text are validated/sanitized via
  `lib/hcls_common/security.py`; never commit real keys or patient-adjacent data
  (`scripts/pre-commit-hook.sh` blocks secret shapes and oversized files).

## Licensing & dependency posture

- **Apache-2.0, permissively licensed throughout.** `lib/hcls_common/license_gate.py` blocks
  non-commercial licenses and flags copyleft — a non-commercial dependency is a **blocker** for an
  Apache-2.0 platform.
- **Exact-pinned dependencies.** Deterministic dependency versions are part of the reproducibility
  story (21 CFR Part 11 spirit): pin exact versions in the constraints/lockfile so a run can be
  reconstructed. Gate heavy/GPU stacks behind optional extras.

## Clinician-in-the-loop review gates

- **A qualified clinician is the decision-maker, always.** The platform surfaces evidence, flags,
  candidates, and citations for that clinician to review — it never closes the loop autonomously.
- **Governance gates are non-optional and inherited by construction.** Every service is fronted by
  the input-validation + output-honesty gates via `hcls_common.api_gate`
  (`create_governed_app` / `install_governance`); generated clinical text passes
  `honesty_flags` / `assert_publishable`. Governance being real (never `|| true`) is itself part of
  the defensible posture.

## Do / Don't

**Do:** describe the platform as a research/reference platform and decision support with a clinician
in the loop; emit a reproducibility manifest + MLOps run for reproducible work; state plainly that
real-data HIPAA/GDPR, de-identification, IRB/consent, and provenance are the user's responsibility;
ship only synthetic/public data; pin deps exactly and keep them permissively licensed; keep every
claim conservative and defensible.

**Don't:** call or imply the platform is FDA-cleared, CE-marked, or a medical device; present output
as diagnosis or prescribing; suggest the repo relieves users of data-governance obligations; commit
PHI or unpinned/non-commercial dependencies; skip the manifest on a run you want to be reproducible;
overstate maturity where a regulator or clinician would expect caution.

## Compliance checklist

- [ ] Framing is "research/reference platform, decision support" — no medical-device or clearance claim.
- [ ] Reproducible run emits a `ReproducibilityManifest` + `MLOpsStore` run (software/model/reference/
      hardware + lineage captured).
- [ ] Only synthetic/public data in-repo; real-data HIPAA/GDPR/IRB/de-id/provenance stated as the
      user's responsibility.
- [ ] Dependencies exact-pinned and permissively licensed (`license_gate` clean).
- [ ] Governance gates wired; generated clinical text is `is_publishable` and clinician-reviewable.
- [ ] Every claim reads conservative and defensible to a skeptical clinician or regulator.

## Related
- `clinical-claim-honesty` — the per-claim honesty checks this posture rests on.
- `11-security-and-secrets` — secrets, data hygiene, disclosure policy, and the governance gates.
- `hcls-core-vision-mission` — the mission this conservative framing protects (honesty is load-bearing).
- `build-housekeeping-standards` — dependency pinning, registry, and "real, never mocked" discipline.
