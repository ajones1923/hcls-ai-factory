#!/usr/bin/env python3
"""HCLS AI Factory - Automated End-to-End Demo Script.

Programmatically drives all 5 intelligence agents through a complete
patient analysis using the DEMO-VCP-001 patient, generating reports
and a timestamped narration log.

Usage:
    python3 scripts/run_demo.py [--host localhost] [--skip-unavailable] [--output-dir demo_output]

Author: Adam Jones
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from datetime import datetime
from pathlib import Path

# Add shared lib
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "lib"))

try:
    from hcls_common.demo_data import (
        DEMO_PATIENT_ID, DEMO_PATIENT_AGE, DEMO_PATIENT_SEX,
        DEMO_BIOMARKERS, DEMO_GENOTYPES, DEMO_STAR_ALLELES,
        DEMO_ONCOLOGY, DEMO_CART, DEMO_IMAGING,
        get_demo_patient_summary,
    )
except ImportError:
    print("ERROR: hcls_common not found. Run from project root.")
    sys.exit(1)

# ANSI colors
GREEN = "\033[92m"
YELLOW = "\033[93m"
RED = "\033[91m"
CYAN = "\033[96m"
BOLD = "\033[1m"
RESET = "\033[0m"
DIM = "\033[2m"


def banner():
    print(f"""{GREEN}{BOLD}
╔═══════════════════════════════════════════════════════════════╗
║       HCLS AI Factory — Automated Demo                       ║
║       Patient DNA → Drug Candidates in <5 Hours              ║
║       Running on NVIDIA DGX Spark ($4,699)                    ║
╚═══════════════════════════════════════════════════════════════╝{RESET}
""")


def log_stage(stage: str, msg: str, elapsed: float = 0):
    ts = datetime.now().strftime("%H:%M:%S.%f")[:-3]
    elapsed_str = f" ({elapsed:.1f}s)" if elapsed else ""
    print(f"  {DIM}{ts}{RESET}  {CYAN}[{stage}]{RESET}  {msg}{YELLOW}{elapsed_str}{RESET}")


def log_result(key: str, value, indent: int = 12):
    pad = " " * indent
    print(f"{pad}{GREEN}✓{RESET} {key}: {BOLD}{value}{RESET}")


def log_error(stage: str, msg: str):
    ts = datetime.now().strftime("%H:%M:%S.%f")[:-3]
    print(f"  {DIM}{ts}{RESET}  {RED}[{stage}]{RESET}  {RED}✗ {msg}{RESET}")


def check_service(host: str, port: int, path: str = "/health") -> bool:
    """Check if a service is available."""
    try:
        import urllib.request
        url = f"http://{host}:{port}{path}"
        req = urllib.request.Request(url, method="GET")
        resp = urllib.request.urlopen(req, timeout=3)
        return resp.status == 200
    except Exception:
        return False


def api_post(host: str, port: int, path: str, payload: dict, timeout: int = 30) -> dict:
    """Make a POST request and return JSON response."""
    import urllib.request
    url = f"http://{host}:{port}{path}"
    data = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(url, data=data, method="POST",
                                  headers={"Content-Type": "application/json"})
    resp = urllib.request.urlopen(req, timeout=timeout)
    return json.loads(resp.read().decode("utf-8"))


def run_demo(host: str = "localhost", skip_unavailable: bool = True, output_dir: str = "demo_output"):
    banner()

    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    print(f"  {BOLD}Demo Patient:{RESET} {get_demo_patient_summary()}")
    print(f"  {BOLD}Output Dir:{RESET}   {out_path.absolute()}")
    print(f"  {BOLD}Host:{RESET}         {host}")
    print()

    results = {}
    total_start = time.perf_counter()

    # ─── Stage 1: Service Discovery ──────────────────────────
    print(f"{BOLD}{'═' * 60}{RESET}")
    print(f"{BOLD}  Stage 0: Service Discovery{RESET}")
    print(f"{BOLD}{'═' * 60}{RESET}")

    services = {
        "Biomarker Agent":  {"port": 8528, "path": "/v1/health"},
        "Oncology Agent":   {"port": 8527, "path": "/health"},
        "CAR-T Agent":      {"port": 8521, "path": "/health"},
        "Imaging Agent":    {"port": 8525, "path": "/health"},
        "Drug Discovery":   {"port": 8505, "path": "/healthz"},
        "RAG/Chat API":     {"port": 5001, "path": "/health"},
        "Milvus":           {"port": 19530, "path": None},
    }

    available = {}
    for name, info in services.items():
        if info["path"]:
            up = check_service(host, info["port"], info["path"])
        else:
            # TCP check for Milvus
            import socket
            try:
                s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                s.settimeout(2)
                s.connect((host, info["port"]))
                s.close()
                up = True
            except Exception:
                up = False

        status = f"{GREEN}● Online{RESET}" if up else f"{RED}● Offline{RESET}"
        print(f"    {name:20s} :{info['port']}  {status}")
        available[name] = up

    print()

    # ─── Stage 1: Biomarker Analysis ─────────────────────────
    print(f"{BOLD}{'═' * 60}{RESET}")
    print(f"{BOLD}  Stage 1: Precision Biomarker Analysis{RESET}")
    print(f"{BOLD}{'═' * 60}{RESET}")

    if available.get("Biomarker Agent"):
        try:
            t0 = time.perf_counter()

            # Full analysis
            log_stage("Biomarker", "Running full patient analysis...")
            payload = {
                "patient_id": DEMO_PATIENT_ID,
                "age": DEMO_PATIENT_AGE,
                "sex": DEMO_PATIENT_SEX,
                "biomarkers": DEMO_BIOMARKERS,
                "genotypes": DEMO_GENOTYPES,
                "star_alleles": DEMO_STAR_ALLELES,
            }
            result = api_post(host, 8528, "/v1/analyze", payload)
            elapsed = time.perf_counter() - t0

            log_stage("Biomarker", "Analysis complete", elapsed)
            log_result("Biological Age", f"{result.get('biological_age', '?')}y")
            log_result("Age Acceleration", f"{result.get('age_acceleration', '?')}y")
            log_result("Mortality Risk", result.get('mortality_risk', '?'))
            log_result("Disease Trajectories", len(result.get('disease_trajectories', [])))
            log_result("PGx Results", len(result.get('pgx_results', [])))
            log_result("Critical Alerts", len(result.get('critical_alerts', [])))

            for alert in result.get('critical_alerts', []):
                print(f"            {RED}⚠ {alert}{RESET}")

            results["biomarker"] = result

            # Save result
            with open(out_path / "biomarker_analysis.json", "w") as f:
                json.dump(result, f, indent=2)

        except Exception as e:
            log_error("Biomarker", str(e))
    elif not skip_unavailable:
        log_error("Biomarker", "Service unavailable")
        sys.exit(1)
    else:
        log_stage("Biomarker", "Skipped (service unavailable)")

    print()

    # ─── Stage 2: Biological Age (standalone) ────────────────
    print(f"{BOLD}{'═' * 60}{RESET}")
    print(f"{BOLD}  Stage 2: Biological Age Calculation{RESET}")
    print(f"{BOLD}{'═' * 60}{RESET}")

    if available.get("Biomarker Agent"):
        try:
            t0 = time.perf_counter()
            log_stage("BioAge", "Computing PhenoAge...")
            payload = {"age": DEMO_PATIENT_AGE, "biomarkers": DEMO_BIOMARKERS}
            result = api_post(host, 8528, "/v1/biological-age", payload)
            elapsed = time.perf_counter() - t0

            log_stage("BioAge", "PhenoAge computed", elapsed)
            log_result("Chronological", f"{result.get('chronological_age')}y")
            log_result("Biological", f"{result.get('biological_age')}y")
            log_result("Acceleration", f"{result.get('age_acceleration')}y")
            log_result("Mortality Risk", result.get('mortality_risk'))

            results["biological_age"] = result

        except Exception as e:
            log_error("BioAge", str(e))
    else:
        log_stage("BioAge", "Skipped (service unavailable)")

    print()

    # ─── Stage 3: Disease Risk ───────────────────────────────
    print(f"{BOLD}{'═' * 60}{RESET}")
    print(f"{BOLD}  Stage 3: Disease Trajectory Analysis{RESET}")
    print(f"{BOLD}{'═' * 60}{RESET}")

    if available.get("Biomarker Agent"):
        try:
            t0 = time.perf_counter()
            log_stage("DiseaseRisk", "Analyzing 6 disease trajectories...")
            payload = {
                "age": DEMO_PATIENT_AGE,
                "sex": DEMO_PATIENT_SEX,
                "biomarkers": DEMO_BIOMARKERS,
                "genotypes": DEMO_GENOTYPES,
            }
            result = api_post(host, 8528, "/v1/disease-risk", payload)
            elapsed = time.perf_counter() - t0

            log_stage("DiseaseRisk", "Trajectories analyzed", elapsed)
            for t in result.get("trajectories", []):
                risk = t.get("risk_level", t.get("stage", "?"))
                disease = t.get("disease", t.get("name", "?"))
                color = RED if risk in ("HIGH", "CRITICAL") else YELLOW if risk == "MODERATE" else GREEN
                log_result(disease, f"{color}{risk}{RESET}")

            results["disease_risk"] = result

        except Exception as e:
            log_error("DiseaseRisk", str(e))
    else:
        log_stage("DiseaseRisk", "Skipped (service unavailable)")

    print()

    # ─── Stage 4: PGx Mapping ───────────────────────────────
    print(f"{BOLD}{'═' * 60}{RESET}")
    print(f"{BOLD}  Stage 4: Pharmacogenomic Mapping{RESET}")
    print(f"{BOLD}{'═' * 60}{RESET}")

    if available.get("Biomarker Agent"):
        try:
            t0 = time.perf_counter()
            log_stage("PGx", "Mapping star alleles to CPIC recommendations...")
            payload = {"star_alleles": DEMO_STAR_ALLELES, "genotypes": DEMO_GENOTYPES}
            result = api_post(host, 8528, "/v1/pgx", payload)
            elapsed = time.perf_counter() - t0

            log_stage("PGx", "PGx mapping complete", elapsed)
            log_result("Genes Analyzed", result.get("total_genes", 0))
            log_result("Critical Findings", len(result.get("critical_findings", [])))

            for finding in result.get("critical_findings", []):
                print(f"            {RED}⚠ {finding}{RESET}")

            results["pgx"] = result

        except Exception as e:
            log_error("PGx", str(e))
    else:
        log_stage("PGx", "Skipped (service unavailable)")

    print()

    # ─── Stage 5: RAG Evidence Query ─────────────────────────
    print(f"{BOLD}{'═' * 60}{RESET}")
    print(f"{BOLD}  Stage 5: RAG Evidence Retrieval{RESET}")
    print(f"{BOLD}{'═' * 60}{RESET}")

    if available.get("Biomarker Agent"):
        try:
            t0 = time.perf_counter()
            log_stage("RAG", "Querying biomarker knowledge base...")
            payload = {
                "question": "What is the clinical significance of VCP p.R155H mutation in frontotemporal dementia and what biomarker patterns predict disease onset?",
                "patient_profile": {
                    "patient_id": DEMO_PATIENT_ID,
                    "age": DEMO_PATIENT_AGE,
                    "sex": DEMO_PATIENT_SEX,
                    "biomarkers": DEMO_BIOMARKERS,
                    "genotypes": DEMO_GENOTYPES,
                    "star_alleles": DEMO_STAR_ALLELES,
                },
            }
            result = api_post(host, 8528, "/v1/query", payload, timeout=60)
            elapsed = time.perf_counter() - t0

            log_stage("RAG", "Evidence retrieved", elapsed)
            log_result("Evidence Count", result.get("evidence_count", 0))
            log_result("Collections", result.get("collections_searched", 0))
            log_result("Search Time", f"{result.get('search_time_ms', 0):.0f}ms")

            # Show answer preview
            answer = result.get("answer", "")
            if answer:
                preview = answer[:200] + "..." if len(answer) > 200 else answer
                print(f"            {DIM}Answer: {preview}{RESET}")

            results["rag_query"] = result

        except Exception as e:
            log_error("RAG", str(e))
    else:
        log_stage("RAG", "Skipped (service unavailable)")

    print()

    # ─── Summary ─────────────────────────────────────────────
    total_elapsed = time.perf_counter() - total_start

    print(f"{GREEN}{BOLD}{'═' * 60}{RESET}")
    print(f"{GREEN}{BOLD}  Demo Complete{RESET}")
    print(f"{GREEN}{BOLD}{'═' * 60}{RESET}")
    print()
    print(f"  {BOLD}Patient:{RESET}        {DEMO_PATIENT_ID}")
    print(f"  {BOLD}Total Time:{RESET}     {total_elapsed:.1f}s")
    print(f"  {BOLD}Stages Run:{RESET}     {len(results)}")
    print(f"  {BOLD}Output:{RESET}         {out_path.absolute()}")
    print()

    # Save combined results
    summary = {
        "demo_patient": DEMO_PATIENT_ID,
        "timestamp": datetime.now().isoformat(),
        "total_time_seconds": round(total_elapsed, 2),
        "stages_completed": list(results.keys()),
        "results": results,
    }
    with open(out_path / "demo_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    print(f"  {GREEN}✓{RESET} Summary saved to {out_path / 'demo_summary.json'}")
    print()

    # Mission statement
    print(f"  {DIM}\"No parent should ever have to bury a child due to disease.\"{RESET}")
    print(f"  {DIM}Patient DNA → Drug Candidates in <5 Hours on a $4,699 DGX Spark{RESET}")
    print()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="HCLS AI Factory Automated Demo")
    parser.add_argument("--host", default="localhost", help="Service host")
    parser.add_argument("--skip-unavailable", action="store_true", default=True,
                        help="Skip unavailable services instead of failing")
    parser.add_argument("--output-dir", default="demo_output", help="Output directory")
    args = parser.parse_args()

    run_demo(host=args.host, skip_unavailable=args.skip_unavailable, output_dir=args.output_dir)
