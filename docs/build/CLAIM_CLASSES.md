# Claim-integrity defect classes, and what a sweep for them found

Five hand-reviewed subjects produced a small repeating taxonomy rather than five unique
problems. `scripts/sweep_claim_classes.py` looks for those classes across all 17 at once.

**It is a triage scanner, not a gate.** It is tuned for recall: most flags are fine on
reading, and a subject with zero flags has *not* passed anything — it has not been looked at.
Only class C is decidable by machine, and that one is already enforced as a test. The scanner
is deliberately not wired into CI, because a scanner reporting PASS for what it merely failed
to pattern-match would be this program's own defect, one level up.

## The classes

| | Class | First found in |
|---|---|---|
| **A** | A value presented as derived that is invented | TSC trajectory prior; GPU telemetry (v2.2.0) |
| **B** | A computed honesty verdict that reaches only the log | precision-oncology `evidence_sufficiency` |
| **C** | A declared tier contradicting the entry's own description, either direction | ACMG panel; TSC / imaging / oncology |
| **D** | A simulated number in a field named as a measurement | clinical-imaging `confidence` |

## What the first sweep found, 2026-11-04

**Class B is systemic.** The defect found by hand in precision-oncology is present in six
more subjects, and the sweep found them in one run:

| Subject | Where | Status |
|---|---|---|
| precision-oncology | `agent.py` | **fixed** — found by hand |
| rare-disease-diagnostic | `src/agent.py` | **fixed** |
| clinical-trial | `src/agent.py` | **fixed** |
| neurology | `src/agent.py` | **fixed** |
| cardiology | `src/agent.py` | **fixed** |
| pharmacogenomics | `src/agent.py` | **fixed** — disclosure prepended to the answer, since it returns dosing guidance |
| clinical-imaging | `agent/src/agent.py` | **fixed** — a different defect, see below |
| precision-intelligence | `src/rag_engine.py` | **false positive**, see below |

| cart | `src/agent.py`, `app/cart_ui.py` | **fixed** — missed by the scanner, see below |

All closed on 2026-11-04. Eight real instances, one false positive, one false negative.

### The scanner also missed one, which is the more useful failure

**cart had the defect and was reported clean.** Class B asks whether the verdict vocabulary
reaches the response surface. In cart it does — in `app/cart_ui.py`:

```python
if quality == "insufficient" and plan.sub_questions:
```

which is the discarding pattern itself, copy-pasted into the UI. The word was present, so
nothing was flagged. **A scanner that reads a second copy of the defect as evidence of the
fix is worse than no scanner**, because it converts an unexamined subject into an apparently
examined one.

It was found by hand, reviewing cart in roster order after the sweep had cleared it. The
detector now strips `if <var> == "<verdict>"` comparison lines out of the surface text before
matching, and finds cart against the pre-fix code — verified by stashing the fix: B went 0 → 2.

cart also turned out to have two problems rather than one. The API path dropped the grade
like the others; the Streamlit UI *displayed* it, but **before** the sub-question expansion
ran, so a user saw a grade describing the evidence that triggered widening rather than the
evidence the answer was built from. A stale disclosure reads as a current one.

### Two of the eight were not the defect the scanner guessed

**clinical-imaging was worse.** `evaluate_evidence` did not drop its verdict — it had **no
production caller at all.** Eight tests in `test_agent.py` asserted its classification
boundaries, so the suite reported a working evidence check while nothing in the pipeline ever
ran it. A tested guard with no caller is worse than no guard: no guard is visibly absent,
this one was visibly present and green. Same shape as therapeutic-discovery's `assert_real`
helpers, written for a real incident and then left uncalled. Now wired into `run()`, with a
test that parses the AST to assert the call exists and is not another test.

**precision-intelligence was fine.** `_verify_answer` returns a dict attached as
`result["verify"]`, and a blocking honesty violation appends a research-use correction to the
answer text. The scanner flagged it because the literal word `"unverified"` never appears
under `api/` or `app/` — a vocabulary match, which is the stated limit of class B. It is
recorded in `REVIEWED_FALSE_POSITIVES` and **still counted in the totals**, because a scanner
that silently drops what someone once dismissed stops being evidence of anything.

Three of them — rare-disease-diagnostic, clinical-trial and neurology — are **byte-identical
copies** of the same block:

```python
quality = self.evaluate_evidence(response.results)
if quality == "insufficient" and plan.sub_questions:
    for sub_q in plan.sub_questions[:2]:
        sub_response = self.rag.search(sub_q, top_k=top_k)
        if sub_response:
            response.results.extend(sub_response)

return response          # `quality` is discarded here
```

The verdict decides whether to widen the search and is then dropped. After the expansion the
answer is returned with no record of how well-grounded it was, so an answer built on almost
nothing is indistinguishable from one built on good evidence. In rare-disease-diagnostic there
is not even a log line. These are diagnostic and therapy-selection agents.

This is the sweep paying for itself: one hand-found defect became a measured seven-subject
finding in a single run, and in roster order the last of them would have surfaced twelve
subjects from now.

## Calibration against the hand reviews

The sweep's hits on the five already-read subjects match what reading them concluded, which
is the only evidence that its silence elsewhere means anything:

- **clinical-imaging A=47, D=5** — the mock NIM clients. Read and recorded as a known
  weakness: the values are flagged `is_mock` but sit in measurement-named fields.
- **therapeutic-discovery A=2, D=2** — `MockDiffDockClient`. Read and passed: the fields are
  renamed `simulated_score` and the `score_type` is one no normalizer accepts.
- **tuberous-sclerosis A=1** — `_synthetic_nhd`. The defect found by hand, now fixed.
- **genomic-foundation A=0** — matches the hand review, which passed check 2.
- **precision-oncology A=0** — matches: every `time.time()` there is a real elapsed
  measurement.

### The first draft of this table was wrong, which is the point of writing it

A regex over raw source read the repository's own prose as code. `genomic-foundation` was
flagged for `hash(str(time.time()))` occurring in the **comment recording that the
fabrication had been removed**, and `clinical-imaging` for `np.random` inside quoted strings
in a demo-guide generator. I wrote that first hit up here as "a rate-limiter clock in
security.py" without opening it — a guess that happened to be false.

A scanner that reports the description of a fixed defect as a defect is a bad failure for
this program specifically, and it inflates exactly the numbers a reader would use to judge
the sweep. It now blanks comments and string literals with `tokenize` before matching:
genomic-foundation 1 → 0, therapeutic-discovery 4 → 2, clinical-imaging 49 → 47. Class B
still reads raw lines on purpose, because the returned string literal *is* its signal.

## Limits worth stating

- Class A cannot tell a seeded synthetic-cohort generator (intended, declared) from a
  fabricated clinical value. It marks the likely-intentional ones and leaves the judgement to
  a reader.
- Class B matches on vocabulary. A verdict surfaced under a different word is a false
  positive; a verdict with no vocabulary at all is a false negative.
- Class D only looks inside files that already know about mocking.
- Nothing here addresses **check 4**, whether a clinical constant is cited to the right
  primary source. That is not mechanizable and stays per-subject reading.
