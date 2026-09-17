# Copyright 2026 Adam Jones
# SPDX-License-Identifier: Apache-2.0
# Part of the HCLS AI Factory — https://github.com/ajones1923/hcls-ai-factory
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""Headers for calling a governed HCLS API from a first-party client.

The API gate went fail-closed on 2026-09-16 — every clinical route requires `X-API-Key`, which
closed a real hole: a POST from the LAN with no credential was returning 200 and generating
clinical prose. The demo runner and the clinical eval were updated in the same commit.

**Six Streamlit UIs were not.** They call their own agent's API with
`headers={"Content-Type": "application/json"}` and nothing else, so from the moment the gate was
enabled every query they make returns 401. The UI still loads, still renders, and still reports
no error until you ask it something — the failure is one layer in, which is why it survived a
fleet-wide "32/32 healthy".

Use this rather than hand-rolling the header, so the next credential change is one edit.
"""
from __future__ import annotations

import os


def auth_headers(extra: dict | None = None, service: str = "") -> dict:
    """Return request headers carrying the API key when one is configured.

    Mirrors `hcls_common.api_auth.resolve_key`: a per-service key wins over the shared one, and
    when neither is set the headers are returned unchanged — so a deployment with the gate
    disabled keeps working.
    """
    headers = dict(extra or {})
    key = None
    if service:
        slug = service.upper().replace("-", "_").replace(" ", "_")
        key = os.environ.get(f"HCLS_API_KEY_{slug}")
    key = key or os.environ.get("HCLS_API_KEY")
    if key:
        headers["X-API-Key"] = key
    return headers
