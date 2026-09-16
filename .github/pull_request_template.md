## What changed, and why

<!-- The change and the reason. If it fixes a defect, say what the defect actually was. -->

## How it was verified

<!-- Commands run and what they returned. "Tests pass" is not verification; the output is. -->

```
# e.g.
# ruff check --select E9,F82,F811,F706,F707 core lib scripts
# ( cd lib/hcls_common && ../../.venv/bin/python -m pytest -q )
# .venv/bin/python scripts/run_all_tests.py | tail -2
# .venv/bin/python scripts/validate_registry.py --probe | tail -1
```

## Checklist

- [ ] The merge gate passes (lint · platform tests · registry · subject suites · docs strict)
- [ ] If a capability's status, port or I/O changed — `validate_registry.py --probe` is clean
- [ ] If clinical output changed — `scripts/run_clinical_eval.py` still passes
- [ ] If a claim changed — it is backed by something that ran, not by intent
- [ ] No data, weights, keys or `.env` in the diff

## Honesty check

<!-- This project's thesis is that what it says about itself is true. If this PR changes a
     status, a badge, a count or a claim in the docs, state what now backs it. -->
