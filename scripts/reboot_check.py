#!/usr/bin/env python3
"""Prove the platform survived a reboot unattended (PRD acceptance criterion A3).

A3 is the one criterion that cannot be verified by running something once — it is about what
happens when nobody is watching. This script is the check, automated, so that "did it come back?"
is one command with one verdict instead of eight commands and a judgement call.

**What actually brings the fleet back.** There are no systemd units for these 32 services. The
entire boot path is:

    cron  */5  ->  health-monitor.sh fix  ->  restarts anything it finds down
    docker       restart=unless-stopped     ->  Milvus, etcd, MinIO
    cron  @reboot                            ->  re-applies the Caddy TLS edge

So the platform is expected to be DOWN for up to five minutes after boot, and then recover with
no human involvement. Running this immediately after login will show failures that are not
failures; `--wait` blocks for the first supervisor tick instead.

**Why each check is here rather than inferred from service health.** All three known reboot
failures leave the status table reading 32/32:

    no ANTHROPIC_API_KEY   services degrade to retrieval-only and still report healthy
    collections not loaded search returns "collection not loaded"; the service is up
    GPU cannot allocate    small allocations succeed, large models fail

Liveness is not readiness, and readiness is not correctness. Each is probed directly.

Usage:
    .venv/bin/python scripts/reboot_check.py            # the check
    .venv/bin/python scripts/reboot_check.py --wait     # wait for the first supervisor tick
    .venv/bin/python scripts/reboot_check.py --eval     # also run the 27-case clinical eval
"""
from __future__ import annotations

import argparse
import json
import os
import pathlib
import re
import subprocess
import time
import urllib.error
import urllib.request

ROOT = pathlib.Path(__file__).resolve().parent.parent
PY = str(ROOT / ".venv" / "bin" / "python")
results: list[tuple[str, bool, str]] = []


def record(name: str, ok: bool, detail: str = "") -> bool:
    results.append((name, ok, detail))
    print(f"  {'PASS' if ok else 'FAIL'}  {name:38s}{detail}")
    return ok


def run(cmd: list[str], timeout: int = 900) -> str:
    try:
        r = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, timeout=timeout)
        return r.stdout + r.stderr
    except subprocess.TimeoutExpired:
        return "TIMEOUT"


def uptime_minutes() -> float:
    try:
        return float(pathlib.Path("/proc/uptime").read_text().split()[0]) / 60
    except Exception:
        return -1.0


def check_supervisor() -> None:
    out = run(["./health-monitor.sh", "status"])
    m = re.search(r"Total:\s*(\d+)/(\d+)\s*services healthy", out)
    if not m:
        record("supervisor reports all services", False, "could not parse status output")
        return
    up, total = int(m.group(1)), int(m.group(2))
    record("supervisor reports all services", up == total, f"{up}/{total}")
    g = re.search(r"GPU .*?\(([\d.]+) GiB allocatable\)", out)
    if g:
        free = float(g.group(1))
        record("GPU memory allocatable", free >= 16,
               f"{free} GiB" + ("" if free >= 16 else "  — page cache holds unified memory; "
                                "sudo sysctl -w vm.drop_caches=3"))
    elif "only" in out and "allocatable" in out:
        record("GPU memory allocatable", False,
               "below the supervisor's threshold — sudo sysctl -w vm.drop_caches=3")


def check_collections() -> None:
    out = run([PY, "scripts/load_collections.py"])
    failed = re.search(r"failed\s+(\d+)", out)
    vectors = re.search(r"([\d,]+)\s+vectors", out)
    ok = bool(failed) and failed.group(1) == "0"
    detail = f"failed {failed.group(1) if failed else '?'}"
    if vectors:
        n = int(vectors.group(1).replace(",", ""))
        detail += f", {n:,} vectors"
        ok = ok and n > 0          # a loaded-but-empty corpus answers 200 and returns nothing
    record("Milvus collections loaded", ok, detail)


def check_corpus() -> None:
    """Did the corpus SURVIVE the reboot — not, is it big enough.

    Those are different questions and the first version of this check confused them: it failed
    A3 because 44 collections are empty, which is a content gap that predates any reboot and is
    not a defect in the boot path. A3 asks whether the platform came back. So this fails only if
    the corpus is GONE — Milvus unreachable, or every vector missing — and reports the depth
    separately as information, because an empty corpus right after a boot is exactly what a lost
    volume looks like.
    """
    out = run([PY, "scripts/check_corpus.py"])
    import re as _re
    m = _re.search(r"(\d+) of (\d+) collections are empty", out)
    total = sum(int(x.replace(",", "")) for x in _re.findall(r"^\s+\S.*?\s(\d[\d,]*)\s*(?:<<|$)",
                                                             out, _re.M)) if m else 0
    if not m:
        record("corpus survived the restart", False, "could not read the corpus")
        return
    record("corpus survived the restart", total > 0,
           f"{total:,} vectors across {int(m.group(2)) - int(m.group(1))} populated collections")
    thin = _re.search(r"(\d+) subject\(s\) below", out)
    if thin or m.group(1) != "0":
        print(f"  note  corpus depth                        {m.group(1)}/{m.group(2)} collections "
              f"empty"
              + (f", {thin.group(1)} subjects under the vector floor" if thin else "")
              + "  — content gap, not a boot-path defect; scripts/check_corpus.py")


def check_registry() -> None:
    out = run([PY, "scripts/validate_registry.py", "--probe"])
    probe = next((ln.strip() for ln in out.splitlines() if "endpoint" in ln), "")
    record("every live endpoint answers", "Every live endpoint answered" in out, probe[:60])


def check_demos() -> None:
    out = run([PY, "scripts/run_demo.py", "--check-all"])
    m = re.search(r"(\d+)/(\d+) demos have their prerequisites met", out)
    if not m:
        record("demo prerequisites", False, "could not parse")
        return
    record("demo prerequisites", m.group(1) == m.group(2), f"{m.group(1)}/{m.group(2)}")


def _post(port: int, path: str, payload: dict, timeout: int = 180) -> dict | None:
    body = json.dumps(payload).encode()
    headers = {"Content-Type": "application/json"}
    if os.getenv("HCLS_API_KEY"):
        headers["X-API-Key"] = os.environ["HCLS_API_KEY"]
    req = urllib.request.Request(f"http://localhost:{port}{path}", data=body,
                                 headers=headers, method="POST")
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return json.loads(r.read().decode())
    except Exception:
        return None


def check_retrieval() -> None:
    """Retrieval, probed directly: 'collection not loaded' still returns HTTP 200."""
    d = _post(8522, "/search", {"question": "CD19 CAR-T toxicity"})
    if d is None:
        record("retrieval returns evidence", False, "no response from cart :8522")
        return
    n = d.get("evidence_count") or len(d.get("results") or d.get("evidence") or [])
    record("retrieval returns evidence", bool(n), f"{n} passages")


def check_synthesis() -> None:
    """The API key is the difference between an answer and a stub, and both return 200.

    Probed by asking, not by grepping a log: a log line proves the client was constructed, not
    that the key still works.
    """
    d = _post(8508, "/query", {"question": "Which gene determines clopidogrel response?"})
    if d is None:
        record("LLM synthesis is live", False, "no response from pharmacogenomics :8508")
        return
    a = (d.get("answer") or d.get("response") or "")
    stub = a.startswith("Search completed") or len(a) < 400
    record("LLM synthesis is live", not stub,
           f"{len(a)} chars" + ("  — looks like the retrieval-only stub; is ANTHROPIC_API_KEY "
                                "loaded?" if stub else ""))


def check_eval() -> None:
    out = run([PY, "scripts/run_clinical_eval.py"], timeout=3600)
    m = re.search(r"(\d+)/(\d+) clinically correct", out)
    if not m:
        record("clinical eval", False, "could not parse")
        return
    record("clinical eval", m.group(1) == m.group(2), f"{m.group(1)}/{m.group(2)} correct")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--wait", action="store_true",
                    help="wait up to 6 min for the first supervisor tick before checking")
    ap.add_argument("--eval", action="store_true", help="also run the clinical eval (slow, costs)")
    a = ap.parse_args()

    up = uptime_minutes()
    print(f"\n  Reboot check (A3) — uptime {up:.0f} min\n")
    if 0 <= up < 6 and not a.wait:
        print("  NOTE: cron runs the supervisor every 5 minutes and that is the ONLY thing that\n"
              "        starts these services. Under 6 minutes of uptime, use --wait.\n")

    if a.wait:
        deadline = time.time() + 400
        while time.time() < deadline:
            if re.search(r"Total:\s*(\d+)/\1\s*services healthy",
                         run(["./health-monitor.sh", "status"])):
                break
            time.sleep(20)

    check_supervisor()
    check_collections()
    check_corpus()
    check_registry()
    check_demos()
    check_retrieval()
    check_synthesis()
    if a.eval:
        check_eval()

    failed = [n for n, ok, _ in results if not ok]
    print()
    if failed:
        print(f"  A3 NOT MET — {len(failed)} check(s) failed: {', '.join(failed)}")
        print("  Each failure is a real defect in the boot path, not a step to do by hand.\n"
              "  Remedies: docs/build/REBOOT_CHECK.md")
        return 1
    # A3 is a claim about the BOOT PATH. Passing on a box that has been up for weeks proves the
    # checks work and the fleet is healthy -- it does not prove anything came back by itself.
    # Saying "A3 MET" here would be the same overclaim the platform's own gate exists to block.
    if up > 120:
        print(f"  All checks pass — but uptime is {up/60:.0f} h, so this is NOT a post-reboot run.")
        print("  A3 is met only by a clean pass shortly after a boot. Re-run then:\n"
              "      .venv/bin/python scripts/reboot_check.py --wait")
        return 0
    print(f"  A3 MET — every check passed {up:.0f} min after boot, with no manual step.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
