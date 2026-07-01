"""
Validation-scorecard endpoint (PRD §5; master paper §16). Runs the evaluation harness
against the cohort's known ground truth and returns measured performance — the numbers
that turn "looks right" into "is right." Construct validity on synthetic data; prospective
clinical validation is institutional Phase-1 (the response carries that caveat).
"""
from __future__ import annotations

from fastapi import APIRouter, Request

from src.eval import evaluate

router = APIRouter(prefix="/eval", tags=["eval"])


@router.get("")
def eval_scorecard(request: Request):
    app = request.app
    return evaluate(app.state.engine, app.state.manifest, app.state.featured)
