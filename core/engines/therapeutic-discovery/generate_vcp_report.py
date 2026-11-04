#!/usr/bin/env python3
"""
VCP -> Drug Candidate Report Generator — WITHDRAWN 2026-11-04.

This file was a second, older copy of the fabricated flagship PDF report. Every number in it was
a literal: a variant at chr9:35,065,263 attributed to "NVIDIA Parabricks 4.6 on DGX Spark"; an
AlphaMissense score of 0.89 and a CADD of 28.5 from no annotator; four candidates (docking -8.62
to -10.95 kcal/mol, QED, MW, logP, composite score) under the prose "NVIDIA BioNeMo MolMIM was
used to generate novel molecules ... docked against the VCP structure using DiffDock"; and a
"Pipeline Performance Metrics" table whose Parabricks figure (30-110 minutes) contradicted the
other generator's (120-240 minutes) for the same stage. The module imported only
json/datetime/pathlib/reportlab: no NIM client, no docking code, no RDKit descriptor. Nothing in
it was computed, and it shared the demo variant's coordinates and alleles with two other files
that spell them differently.

It had no caller in the repo — only README and docs references — so rather than port the whole
layout to injected data, the fabricated report is withdrawn. The surviving generator,
`generate_vcp_report_enhanced.py`, renders an injected run and raises `ReportDataUnavailable`
when none is supplied; `run_cloud_nim_report.py` produces such a run end to end.
"""

WITHDRAWN_MESSAGE = (
    "generate_vcp_report.py was withdrawn on 2026-11-04: it emitted a hard-coded result set "
    "under the Parabricks / BioNeMo MolMIM / DiffDock names. Use "
    "generate_vcp_report_enhanced.VCPReportGeneratorEnhanced(results=<run>) with a real run "
    "(see run_cloud_nim_report.py), which refuses rather than inventing candidates."
)


class VCPReportGenerator:
    """Refuses. Kept only so an old import fails with an explanation instead of a NameError."""

    def __init__(self, *args, **kwargs):
        raise RuntimeError(WITHDRAWN_MESSAGE)


def main():
    raise SystemExit(WITHDRAWN_MESSAGE)


if __name__ == "__main__":
    main()
