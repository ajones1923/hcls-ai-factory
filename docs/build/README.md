# Build & Audit Documentation — batch-01 (2026-08-15)

Produced during `batch-01-hcls-aif-core-updates.txt`. Every figure in these documents was measured
by running the code, not by reading it. Where a first measurement was wrong, the correction is
stated in the document rather than removed.

| Document | What it answers |
|---|---|
| [`GAP_ANALYSIS.md`](GAP_ANALYSIS.md) | What exists vs what is needed, with evidence |
| [`PRD.md`](PRD.md) | What to build, in what order, and what is blocked |
| [`BUILD_GUIDE.md`](BUILD_GUIDE.md) | Cold clone → running platform |
| [`PORT_MAP.md`](PORT_MAP.md) | The canonical allocation and the convention behind it |
| [`BEST_PRACTICES_AUDIT.md`](BEST_PRACTICES_AUDIT.md) | Security, design, and factual accuracy |
| [`LAYOUT_REVIEW.md`](LAYOUT_REVIEW.md) | Structural conformance across the 17 subjects |
| [`SCORECARD.md`](SCORECARD.md) | Scored assessment + ranked recommendations |
| [`RESOURCE_REPORT.md`](RESOURCE_REPORT.md) | CPU / GPU / memory / disk / network / models |
| [`REPO_AND_SITE_REVIEW.md`](REPO_AND_SITE_REVIEW.md) | GitHub layout + live-site accuracy |
| [`DEPLOY_RUNBOOK.md`](DEPLOY_RUNBOOK.md) | Publishing to hcls-ai-factory.org |
| [`GATED_SOFTWARE_PRD.md`](GATED_SOFTWARE_PRD.md) | What needs credentials and what it blocks |
| [`GATED_SOFTWARE_BUILD_GUIDE.md`](GATED_SOFTWARE_BUILD_GUIDE.md) | Step-by-step gated install |
| `docx/GATED_SOFTWARE_PRD.docx` · `docx/GATED_SOFTWARE_BUILD_GUIDE.docx` | Word versions, for the final upgrade pass |

## The three things that matter most

1. **Authenticate the 11 open clinical endpoints**, and stop emitting `X-HCLS-Governed` on requests
   that were not gated. (`BEST_PRACTICES_AUDIT.md` §1.2)
2. **Install CUDA PyTorch.** The GPU sat at 0 % for 200 minutes on a GB10.
   (`GATED_SOFTWARE_PRD.md` G1)
3. **Choose one runtime model.** Two half-maintained deployment paths are why the port map drifted.
   (`LAYOUT_REVIEW.md` §4)

## Verify any claim here

```bash
.venv/bin/python scripts/run_all_tests.py        # 17 subjects · 8,402 passed
.venv/bin/python scripts/validate_registry.py    # convention + supervisor cross-check
```
