"""
TSC Intelligence Engine — FastAPI app + portal (PRD §2.4; port 8560).

Serves a unified portal landing page (links to the three clinician surfaces, live
status, API docs), plus JSON endpoints for events, surfaces, agents, and provenance.
At startup it builds the synthetic cohort (if needed) and enrolls every patient into an
in-process engine, so the API and surfaces all read consistent, deterministic state.

Run:  uvicorn api.main:app --host 0.0.0.0 --port 8560
"""
from __future__ import annotations

import sys
from contextlib import asynccontextmanager
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fastapi import FastAPI                                  # noqa: E402
from fastapi.middleware.cors import CORSMiddleware           # noqa: E402
from fastapi.responses import HTMLResponse                   # noqa: E402

from config.settings import settings                         # noqa: E402
from src.cohort.build import build_cohort                    # noqa: E402
from src.cohort.loader import featured_map, load_manifest, load_patient  # noqa: E402
from src.orchestrator.graph import Orchestrator              # noqa: E402
from src.orchestrator.state import SqliteEventStore          # noqa: E402
from src.utils.model_router import get_router                # noqa: E402


@asynccontextmanager
async def lifespan(app: FastAPI):
    if not (settings.COHORT_DIR / "manifest.json").exists():
        build_cohort()
    # P3-1 — durable, event-sourced persistence. The event log is the source of truth and
    # persists across restarts (was :memory:, so engine.db stayed empty). Enrollment is
    # idempotent: a patient already in the append-only log is REPLAYED (projections
    # rematerialize from events) rather than re-enrolled, so restarts never duplicate events.
    # Set TSC_SQLITE_PATH=:memory: to opt back into ephemeral per-process state.
    orch = Orchestrator(store=SqliteEventStore(str(settings.SQLITE_PATH)))
    # P3-3 — optionally dispatch enrollment through the compiled LangGraph StateGraph
    # (PRD production runtime). It wraps THIS orchestrator (same store/agents), so
    # app.state.engine stays the plain Orchestrator and every route is unchanged.
    enroller = orch
    runtime = "plain-dispatcher"
    if settings.USE_LANGGRAPH:
        from src.orchestrator.langgraph_graph import LangGraphRunner   # noqa: WPS433
        enroller = LangGraphRunner(orchestrator=orch)
        runtime = "langgraph-stategraph"
    manifest = load_manifest()
    enrolled = replayed = 0
    for p in manifest["patients"]:
        pid = p["patient_id"]
        if orch.store.events_for(pid):
            replayed += 1
        else:
            enroller.enroll(pid, load_patient(pid))
            enrolled += 1
    print(f"[TSC] engine ready — {enrolled} enrolled, {replayed} replayed "
          f"({runtime}) from {settings.SQLITE_PATH}", file=sys.stderr)
    app.state.engine = orch
    app.state.manifest = manifest
    app.state.featured = featured_map()
    yield


def create_app() -> FastAPI:
    app = FastAPI(title="TSC Intelligence Engine API", version="0.1.0",
                  description="Engine 7 of the HCLS AI Factory — SYNTHETIC demonstration data.",
                  lifespan=lifespan)
    app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])

    # P1-7 — optional API-key gate. Fail-closed only when TSC_API_KEY is configured;
    # unset (default) keeps the open trusted-LAN/synthetic-demo posture. Health, portal,
    # and docs stay public so probes and the landing page still work.
    if settings.API_KEY:
        from fastapi import Request                          # noqa: WPS433
        from fastapi.responses import JSONResponse           # noqa: WPS433
        _OPEN_PATHS = ("/health", "/", "/docs", "/openapi.json", "/redoc")

        @app.middleware("http")
        async def _require_api_key(request: Request, call_next):
            if request.url.path in _OPEN_PATHS or request.method == "OPTIONS":
                return await call_next(request)
            if request.headers.get("X-API-Key") != settings.API_KEY:
                return JSONResponse({"detail": "missing or invalid X-API-Key"}, status_code=401)
            return await call_next(request)

    from api.routes import (
        agents, cohort, eval, events, provenance, surfaces, system, variant_curator, viz,
    )
    for r in (events.router, surfaces.router, agents.router, provenance.router,
              cohort.router, eval.router, system.router, variant_curator.router, viz.router):
        app.include_router(r)

    @app.get("/health")
    def health():
        return {
            "status": "ok",
            "engine": "TSC Intelligence Engine (Engine 7)",
            "patients": app.state.manifest["n_patients"],
            "featured": app.state.featured,
            "reasoning": "live" if get_router().online else "offline (deterministic stubs)",
            "watermark": "SYNTHETIC",
        }

    @app.get("/", response_class=HTMLResponse)
    def portal():
        return PORTAL_HTML

    return app


PORTAL_HTML = """<!doctype html><html><head><meta charset="utf-8">
<title>TSC Intelligence Engine — Portal</title>
<style>
 body{font-family:Inter,system-ui,Arial,sans-serif;margin:0;background:#0f1830;color:#e8eef9}
 .wrap{max-width:1140px;margin:0 auto;padding:44px 24px}
 h1{color:#fff;font-size:30px;margin:0 0 4px} .sub{color:#9db4dc;margin:0 0 8px}
 h2{color:#cfe0ff;font-size:15px;text-transform:uppercase;letter-spacing:.05em;margin:34px 0 12px}
 .wm{display:inline-block;background:#3b2330;color:#ffb3c1;border:1px solid #6b3a4d;
     padding:4px 10px;border-radius:6px;font-size:12px;margin:12px 0 8px}
 .grid{display:grid;grid-template-columns:1fr 1fr 1fr 1fr;gap:14px}
 a.card{display:block;background:#16213f;border:1px solid #26345c;border-radius:10px;
     padding:16px 18px;text-decoration:none;color:#e8eef9}
 a.card:hover{border-color:#4f7fd6;background:#1a2750}
 a.card b{color:#7fb0ff;font-size:15px} a.card span{display:block;color:#9db4dc;font-size:12px;margin-top:6px}
 .stats{display:grid;grid-template-columns:repeat(5,1fr);gap:14px}
 .stat{background:#16213f;border:1px solid #26345c;border-radius:10px;padding:14px 16px}
 .stat .n{font-size:26px;font-weight:700;color:#fff} .stat .l{font-size:12px;color:#9db4dc;margin-top:2px}
 .stat.hi .n{color:#ffd166}
 table{width:100%;border-collapse:collapse;margin-top:10px;font-size:13px}
 th{text-align:left;color:#9db4dc;font-weight:600;padding:8px 10px;border-bottom:1px solid #26345c;position:sticky;top:0;background:#0f1830}
 td{padding:7px 10px;border-bottom:1px solid #1c2950}
 tr:hover td{background:#16213f}
 .tbl{max-height:440px;overflow:auto;border:1px solid #26345c;border-radius:10px}
 .badge{display:inline-block;padding:2px 9px;border-radius:20px;color:#fff;font-weight:600;font-size:11px}
 .pill{display:inline-block;background:#22325c;color:#bcd2ff;border:1px solid #34467a;border-radius:5px;padding:1px 7px;margin:1px 3px 1px 0;font-size:11px}
 .star{color:#ffd166;font-weight:700} .mut{font-family:ui-monospace,Menlo,monospace;color:#9fc1ff;font-size:12px}
 .row{display:flex;gap:14px;flex-wrap:wrap;margin-top:18px}
 .row a{color:#7fb0ff;font-size:13px} footer{color:#6f86b3;font-size:12px;margin-top:30px}
</style></head><body><div class="wrap">
<h1>TSC Intelligence Engine</h1>
<p class="sub">Engine 7 · HCLS AI Factory · five agents + deterministic orchestrator · one NVIDIA DGX Spark</p>
<div class="wm">SYNTHETIC demonstration data — decision support, clinician review required. Not FDA-cleared.</div>

<h2>Clinician surfaces</h2>
<div class="grid">
 <a class="card" id="briefing"><b>Pre-Visit Briefing →</b><span>One-screen, mobile-readable. The night before clinic.</span></a>
 <a class="card" id="dashboard"><b>In-Visit Dashboard →</b><span>Variant · phenotype · trajectory · TAND/therapeutics.</span></a>
 <a class="card" id="alerts"><b>Async Alerts →</b><span>Disciplined — source + dismissal on every alert.</span></a>
 <a class="card" href="/docs"><b>API & Docs →</b><span>Events · surfaces · agents · cohort · provenance.</span></a>
</div>

<h2>Population — the engine runs the whole panel</h2>
<div class="stats" id="stats"></div>
<div class="tbl"><table id="cohort"><thead><tr>
 <th>Patient</th><th>Gene</th><th>Zygosity</th><th>VAF</th><th>Variant</th>
 <th>ACMG class</th><th>Trajectory</th><th>TAND</th><th>Surveillance</th>
</tr></thead><tbody></tbody></table></div>

<h2>Validation scorecard — measured against known ground truth</h2>
<div class="stats" id="evalstats"></div>
<div class="small" id="evalcaveat" style="color:#9db4dc;font-size:12px;margin-top:8px"></div>

<h2>Anatomical digital twin — Omniverse (RunPod RTX)</h2>
<div class="small" style="color:#9db4dc;font-size:13px">The Spark authors volumetric OpenUSD scenes from the real engine output — render them with RTX in Omniverse (RunPod). <a id="twin" style="color:#7fb0ff">SEGA trajectory twin →</a> &nbsp;·&nbsp; <a id="twinm" style="color:#7fb0ff">mosaic "powers-of-ten" →</a> &nbsp;·&nbsp; <a href="/viz/population?inline=true" style="color:#7fb0ff">50-patient population →</a><br><span style="color:#6f86b3">Trajectory: the lesion grows along its forecast inside its 50/90% uncertainty envelope, crossing the threshold membrane. Mosaic: a cell field where exactly the recovered VAF carries the variant.</span></div>

<div class="row">
 <a href="/health">/health</a><a href="/cohort">/cohort</a><a href="/eval">/eval</a><a href="/system/usage">/system/usage</a><a href="/agents">/agents</a>
</div>
<footer>Apache 2.0 · runs on a single NVIDIA DGX Spark · every output traceable (model · tier · latency) · Epic/LIMS integration is institutional Phase-1 (not built).</footer>
</div><script>
 const h=location.hostname;
 briefing.href=`http://${h}:8561`; dashboard.href=`http://${h}:8562`; alerts.href=`http://${h}:8563`;
 const CC={"Pathogenic":"#e5484d","Likely Pathogenic":"#e8833a","Variant of Uncertain Significance":"#8b8f9a","—":"#5b6275"};
 fetch('/health').then(r=>r.json()).then(d=>{
   document.querySelector('.wm').textContent =
     `SYNTHETIC demo · ${d.patients} patients · reasoning: ${d.reasoning} · clinician review required · not FDA-cleared`;
 });
 fetch('/cohort').then(r=>r.json()).then(d=>{
   const H=d.highlights, cd=d.distributions.classification;
   const path=(cd['Pathogenic']||0)+(cd['Likely Pathogenic']||0);
   const cards=[
     ['n',d.n_patients,'patients enrolled'],
     ['hi',H.mosaic_recovered,'mosaic variants recovered'],
     ['n',path,'pathogenic / likely path.'],
     ['n',cd['Variant of Uncertain Significance']||0,'variants of uncertain sig.'],
     ['hi',H.lesions_at_threshold+H.lesions_crossing_window,'lesions at/approaching threshold'],
   ];
   document.getElementById('stats').innerHTML = cards.map(c=>
     `<div class="stat ${c[0]=='hi'?'hi':''}"><div class="n">${c[1]}</div><div class="l">${c[2]}</div></div>`).join('');
   const tb=document.querySelector('#cohort tbody');
   tb.innerHTML = d.patients.map(p=>{
     const star=p.featured?`<span class="star">★ ${p.featured}</span> `:'';
     const rec=p.recovered?' <span class="pill">mosaic recovered</span>':'';
     const cls=`<span class="badge" style="background:${CC[p.classification]||'#5b6275'}">${p.classification}</span>`;
     const traj=(p.lesion_flags||[]).map(f=>`<span class="pill">${f}</span>`).join('')||'—';
     const tand=(p.tand||[]).map(t=>`<span class="pill">${t}</span>`).join('')||'—';
     const surv=(p.overdue||[]).map(s=>`<span class="pill">${s} overdue</span>`).join('')||'—';
     const link=`http://${h}:8562`;
     return `<tr><td>${star}<a href="${link}" style="color:#cfe0ff">${p.patient_id}</a>${rec}</td>`+
       `<td>${p.gene}</td><td>${p.zygosity}</td><td>${p.vaf??'—'}</td>`+
       `<td class="mut">${p.variant||'—'}</td><td>${cls}</td>`+
       `<td>${traj}</td><td>${tand}</td><td>${surv}</td></tr>`;
   }).join('');
 });
 const pc=x=>x==null?'—':Math.round(100*x)+'%';
 fetch('/eval').then(r=>r.json()).then(e=>{
   const H=e.headline, det=e.diagnostic_yield.detection;
   const cards=[
     ['n',pc(H.variant_classification_accuracy),'ACMG classification accuracy'],
     ['hi','+'+Math.round(100*det.uplift_points)+' pts','diagnostic detection uplift'],
     ['hi',pc(H.mosaic_recovery_sensitivity),'sub-threshold mosaic recovery'],
     ['n',H.false_benign_on_null_variant,'truncating variants called benign'],
     ['n',pc(H.provenance_completeness),'outputs fully traceable'],
   ];
   document.getElementById('evalstats').innerHTML = cards.map(c=>
     `<div class="stat ${c[0]=='hi'?'hi':''}"><div class="n">${c[1]}</div><div class="l">${c[2]}</div></div>`).join('');
   document.getElementById('evalcaveat').textContent = '⚠ ' + e.caveat;
 });
 fetch('/health').then(r=>r.json()).then(d=>{
   const f=d.featured||{};
   if(f.B){ const t=document.getElementById('twin'); t.href=`/viz/lesion/${f.B}?inline=true`;
     t.textContent=`SEGA trajectory twin (${f.B}) →`; }
   if(f.A){ const m=document.getElementById('twinm'); m.href=`/viz/mosaic/${f.A}?inline=true`;
     m.textContent=`mosaic powers-of-ten (${f.A}) →`; }
 });
</script></body></html>"""

app = create_app()
