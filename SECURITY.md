# Security Policy

## Reporting a vulnerability

If you discover a security vulnerability in the HCLS AI Factory, please report it
**privately** — do not open a public issue.

- **Email:** adam.m.jones@mac.com
- Please include: the affected component (engine / agent / platform layer), a description,
  reproduction steps, and the potential impact.
- You will receive an acknowledgement within **5 business days**, and a status update
  (accepted / more-info-needed / declined) within **15 business days**.
- Please allow a reasonable window to release a fix before any public disclosure. We are
  happy to credit reporters who wish to be named.

## Supported versions

This is an actively developed research platform. Security fixes are applied to the
`main` branch; there is no long-term-support branch. Run the latest `main`.

## Scope

In scope: the code in this repository — the eight engines, eight intelligence agents,
disease-programs, and the shared platform layer (`lib/hcls_common/`: registry, MCP
tool-surface, workflow composer, MLOps, governance gates).

Out of scope: third-party dependencies and vendored tools (e.g. `vendor_proteinmpnn/`,
`vendor_rfdiffusion/`), which should be reported to their respective upstream projects;
and issues that require pre-existing local/root access to the host.

## Important: not a medical device

The HCLS AI Factory is a **research and reference platform**. It is **not** FDA-cleared,
CE-marked, or otherwise approved for clinical use, and it must not be used to make patient
care decisions. Its governance layer labels mock outputs and blocks clinical overclaims by
design, but outputs are for research and demonstration only. Handle any real patient data
in accordance with applicable law (HIPAA, GDPR, etc.) — the repository ships only synthetic
and public reference data.
