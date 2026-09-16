#!/usr/bin/env python3
"""Run one of the 17 per-subject demonstrations and record a transcript.

Keys are E1-E8 (engines), A1-A8 (agents), P1 (TSC program) — deliberately NOT D1-D7, which is the
established patient-story portfolio in docs/demos/index.md. See DEMO_CATALOG.md.

The rule this script exists to enforce (PRD DR-3):

    A demo declared LIVE whose service is unreachable FAILS. It does not quietly degrade to a
    canned result.

That is the same discipline that caught `singlecell-compute` and `structural-biology-engine`
registered `live` in the capability registry while nothing bound their ports.

Usage:
    .venv/bin/python scripts/run_demo.py --list
    .venv/bin/python scripts/run_demo.py E8
    .venv/bin/python scripts/run_demo.py --check-all      # prerequisites only, runs nothing
"""
from __future__ import annotations
import argparse, importlib.util, json, pathlib, socket, subprocess, sys, urllib.error, urllib.request
from datetime import datetime, timezone

ROOT = pathlib.Path(__file__).resolve().parent.parent
TRANSCRIPTS = ROOT / "demo" / "transcripts"


def _headers() -> dict:
    """Request headers, carrying the API key when the gate is on.

    The gate is fail-closed once HCLS_API_KEY is set (lib/hcls_common/api_auth.py), so a runner
    that sends no credential turns 17/17 into 17 x 401 the moment auth is enabled. Health and
    docs stay open by design; the clinical routes do not.
    """
    import os
    h = {"Content-Type": "application/json"}
    key = os.getenv("HCLS_API_KEY")
    if key:
        h["X-API-Key"] = key
    return h


LIVE, REPRESENTATIVE, BURST = "LIVE", "REPRESENTATIVE", "BURST"


class Demo:
    def __init__(self, key, subject, title, label, *, port=None, packages=(),
                 payload=None, gated=(), runner=None, endpoint=None, assertion=None):
        self.key, self.subject, self.title, self.label = key, subject, title, label
        self.port, self.packages, self.payload = port, packages, payload
        self.gated, self.runner = gated, runner
        # endpoint + assertion make an HTTP demo declarative: POST the payload, then prove the
        # response carries real content. HTTP 200 is NOT the pass condition -- see http_demo().
        self.endpoint, self.assertion = endpoint, assertion

    def missing_packages(self):
        return [p for p in self.packages if importlib.util.find_spec(p) is None]

    def port_open(self):
        if not self.port:
            return None
        with socket.socket() as s:
            s.settimeout(1.5)
            return s.connect_ex(("127.0.0.1", self.port)) == 0

    def check(self):
        """-> (ok, [reasons]). A REPRESENTATIVE demo may run without its service."""
        reasons = []
        miss = self.missing_packages()
        if miss:
            reasons.append(f"missing packages: {', '.join(miss)}")
        if self.payload and not (ROOT / self.payload).exists():
            reasons.append(f"missing payload: {self.payload}")
        if self.label == LIVE and self.port is not None and not self.port_open():
            reasons.append(f"declared LIVE but nothing is listening on :{self.port}")
        if self.gated:
            reasons.append(f"gated (informational): {', '.join(self.gated)}")
        blocking = [r for r in reasons if not r.startswith("gated (informational)")]
        return (not blocking), reasons


def http_demo(demo):
    """POST the demo's payload to its endpoint and prove the answer is real.

    A 200 is deliberately NOT the pass condition. The clinical-trial agent, with no trial corpus
    loaded, answers 200 with a single `NCT-PENDING` placeholder and `total_screened: 0` -- which
    looks like a match and is not one. Accepting status codes would let that through, and the
    whole point of this runner (PRD DR-3) is that a LIVE demo which cannot really answer must
    FAIL rather than quietly degrade to a canned result. So each demo supplies an assertion that
    has to find actual content, and raises with a specific, actionable reason when it cannot.
    """
    def run(log):
        body = json.loads((ROOT / demo.payload).read_text())
        url = f"http://localhost:{demo.port}{demo.endpoint}"
        log(f"POST          {url}")
        first = next((k for k in ("question", "diagnosis", "scale_name") if k in body), None)
        if first:
            log(f"input         {first}={str(body[first])[:88]}")
        req = urllib.request.Request(
            url, data=json.dumps(body).encode(),
            headers=_headers(), method="POST")
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                status, raw = r.status, r.read().decode()
        except urllib.error.HTTPError as e:
            detail = e.read().decode()[:200]
            raise RuntimeError(f"HTTP {e.code} from {demo.endpoint} — {detail}") from None
        log(f"status        {status}")
        data = json.loads(raw)
        for line in demo.assertion(data):      # raises if the answer is not real
            log(line)
    return run


def assert_imaging(d):
    ans = (d.get("answer") or "").strip()
    if len(ans) < 80:
        raise RuntimeError(f"answer too thin to be a real RAG response ({len(ans)} chars) — "
                           "is the imaging corpus loaded?")
    yield f"answer        {len(ans)} chars"
    cites = d.get("citations") or d.get("sources") or []
    yield f"citations     {len(cites)}"
    yield f"excerpt       {ans[:150].replace(chr(10), ' ')}"


def assert_ascvd(d):
    if "score" not in d:
        raise RuntimeError(f"no score in response: {list(d)[:8]}")
    yield f"calculator    {d.get('calculator')}"
    yield f"score         {d['score']}%  ({d.get('risk_category')})"
    interp = (d.get("interpretation") or "")[:150]
    if interp:
        yield f"interpretation {interp}"


def assert_scale(d):
    if d.get("total_score") is None:
        raise RuntimeError(f"no total_score in response: {list(d)[:8]}")
    yield f"scale         {d.get('scale_name')}"
    yield f"score         {d['total_score']} / {d.get('max_score')}"
    yield f"reading       {d.get('interpretation')} · {d.get('severity_category', '-')}"


def assert_differential(d):
    diff = d.get("differential") or []
    if not diff:
        raise RuntimeError("empty differential — the rare-disease knowledge base returned nothing")
    yield f"differential  {len(diff)} candidate diagnoses"
    for row in diff[:3]:
        yield (f"  candidate   {row.get('disease_id')} {row.get('disease_name')} "
               f"(confidence {row.get('confidence')}, overlap {row.get('phenotype_overlap')})")


def assert_trial_match(d):
    """Rejects the no-corpus placeholder rather than reporting it as a match."""
    matches = d.get("matches") or []
    screened = d.get("total_screened", 0)
    real = [m for m in matches if m.get("trial_id") and m["trial_id"] != "NCT-PENDING"]
    if not screened or not real:
        raise RuntimeError(
            f"no trial corpus loaded — screened {screened} trials and the only result is the "
            "'NCT-PENDING' placeholder. Seed the clinical-trial collections first.")
    yield f"screened      {screened} trials"
    for m in real[:3]:
        yield f"  match       {m.get('trial_id')} score={m.get('overall_score')}"



def assert_biological_age(d):
    """PhenoAge — a deterministic formula, so it must return a number, not a narrative."""
    bio = d.get("biological_age")
    if bio is None:
        raise RuntimeError(f"no biological_age in response: {list(d)[:8]}")
    yield f"chronological  {d.get('chronological_age')} y"
    yield f"biological     {bio} y  (acceleration {d.get('age_acceleration')})"
    yield f"mortality risk {d.get('mortality_risk')}"
    for drv in (d.get("top_aging_drivers") or [])[:3]:
        yield f"  driver       {drv}"


def assert_warfarin(d):
    """IWPC pharmacogenomic dosing — genotype in, mg/week out."""
    dose = d.get("predicted_weekly_dose_mg")
    if not dose:
        raise RuntimeError(f"no predicted_weekly_dose_mg in response: {list(d)[:8]}")
    yield f"algorithm      {d.get('algorithm')}"
    yield f"weekly dose    {dose} mg  ({d.get('predicted_daily_dose_mg')} mg/day)"
    yield f"category       {d.get('dose_category')}"


def assert_autoimmune(d):
    diff = d.get("differential") or []
    if not diff:
        raise RuntimeError("empty differential — autoantibody interpretation returned nothing")
    yield f"differential   {len(diff)} candidates"
    for row in diff[:3]:
        ev = "; ".join((row.get("evidence") or [])[:2])
        yield f"  {row.get('disease')} (score {row.get('score')}) — {ev[:90]}"


def assert_therapy_rank(d):
    tx = d.get("therapies") or []
    if not tx:
        raise RuntimeError("no therapies ranked — is the oncology knowledge base loaded?")
    yield f"ranked         {len(tx)} therapies"
    for t in tx[:3]:
        yield (f"  {t.get('rank')}. {t.get('drug_name')} — evidence {t.get('evidence_level')}, "
               f"from {t.get('source_gene')} {t.get('source_variant')}")


def assert_annotation(d):
    cts = d.get("cell_types") or []
    if not cts:
        raise RuntimeError("no cell types annotated — marker panel returned nothing")
    yield f"annotated      {len(cts)} cell types"
    for c in cts[:4]:
        yield (f"  {c.get('cell_type')} ({c.get('cell_ontology_id')}) "
               f"confidence {c.get('confidence')} via {','.join(c.get('markers') or [])}")



def assert_cart_evidence(d):
    """A1 — cross-collection CAR-T evidence retrieval over the seeded corpus."""
    ev = d.get("evidence") or []
    if not ev:
        raise RuntimeError(
            "no evidence returned — the CAR-T collections are empty or unsearchable; "
            "run core/agents/cart/scripts/seed_*.py")
    yield f"collections   {d.get('collections_searched')} searched"
    yield f"evidence      {len(ev)} passages in {d.get('search_time_ms', '?')} ms"
    for e in ev[:3]:
        yield (f"  [{e.get('collection')}] {e.get('id')} score={round(float(e.get('score', 0)), 3)}")
        yield f"      {str(e.get('text',''))[:120]}"
    yield "decision support for a qualified clinician, not diagnosis"


DEMOS = [
    Demo("E1", "genomic-foundation", "The variant that was always there", REPRESENTATIVE,
         packages=("duckdb", "statsmodels"), gated=("Parabricks (G2)",),
         runner="genomic_foundation"),
    Demo("E2", "precision-intelligence", "Ask the evidence layer a question", LIVE, port=5001,
         runner="precision_intelligence"),
    Demo("E3", "therapeutic-discovery", "From one protein to a hundred candidates", REPRESENTATIVE,
         packages=("rdkit",), gated=("MolMIM (G3)", "DiffDock (G4)"),
         runner="therapeutic_discovery"),
    Demo("E4", "clinical-imaging", "The scan that had already answered", LIVE, port=8524,
         payload="demo/requests/imaging_query.json",
         endpoint="/api/ask", assertion=assert_imaging),
    Demo("E5", "precision-oncology", "The molecular tumour board", LIVE, port=8527,
         payload="demo/requests/oncology_therapy_rank.json",
         endpoint="/api/therapies/rank", assertion=assert_therapy_rank),
    Demo("E6", "cardiology", "Risk that changes management", LIVE, port=8127,
         payload="demo/requests/cardiology_risk.json",
         endpoint="/v1/cardio/risk/ascvd", assertion=assert_ascvd),
    Demo("E7", "structural-biology", "The shape you have to fit", REPRESENTATIVE,
         packages=("Bio",), gated=("CUDA torch (G1)", "ESMFold (G5)"),
         runner="structural_biology"),
    Demo("E8", "single-cell", "Nine populations from one sample", LIVE,
         packages=("scanpy", "anndata"), runner="single_cell"),
    Demo("A1", "cart", "Why this construct, for this patient", LIVE, port=8522,
         payload="demo/requests/cart_query.json",
         endpoint="/search", assertion=assert_cart_evidence),
    Demo("A2", "precision-biomarker", "The marker that changes the decision", LIVE, port=8529,
         payload="demo/requests/biomarker_phenoage.json",
         endpoint="/v1/biological-age", assertion=assert_biological_age),
    Demo("A3", "pharmacogenomics", "Two patients, same dose, different outcome", LIVE, port=8508,
         payload="demo/requests/pharmacogenomics_warfarin.json",
         endpoint="/v1/pgx/dosing/warfarin", assertion=assert_warfarin),
    Demo("A4", "precision-autoimmune", "Before the third flare", LIVE, port=8532,
         payload="demo/requests/autoimmune_differential.json",
         endpoint="/differential", assertion=assert_autoimmune),
    Demo("A5", "neurology", "NIHSS, and what comes next", LIVE, port=8536,
         payload="demo/requests/neurology_nihss.json",
         endpoint="/v1/neuro/scale/calculate", assertion=assert_scale),
    Demo("A6", "clinical-trial", "The trial that was open all along", LIVE, port=8539,
         payload="demo/requests/trial_match.json",
         endpoint="/v1/trial/match", assertion=assert_trial_match),
    Demo("A7", "rare-disease-diagnostic", "Ending the odyssey", LIVE, port=8545,
         payload="demo/requests/rare_disease_diagnose.json",
         endpoint="/v1/diagnostic/diagnose", assertion=assert_differential),
    Demo("A8", "single-cell", "The engine computes, the agent interprets", LIVE, port=8541,
         payload="demo/requests/single_cell_annotate.json",
         endpoint="/v1/sc/annotate", assertion=assert_annotation),
    Demo("P1", "tuberous-sclerosis", "The whole factory, one child", REPRESENTATIVE, port=8560,
         runner="tsc_program"),
]
BY_KEY = {d.key: d for d in DEMOS}


def run_single_cell(log):
    """E8 — the only demo that is real end-to-end today with no gated software and no GPU."""
    import scanpy as sc  # noqa: F401
    src = ROOT / "core/engines/single-cell/src"
    sys.path.insert(0, str(src))
    from single_cell_compute import SingleCellAnalysis, PBMC_MARKERS

    # The bundled copy lives under data/, which the project deliberately never publishes
    # (.gitignore; "Data / weights / secrets stay local"). On a fresh clone it is absent, so fall
    # back to scanpy's own copy of the same public PBMC 3k dataset. That keeps this demo genuinely
    # reproducible by anyone who clones the repo -- which is the whole point of publishing it.
    h5ad = ROOT / "core/engines/single-cell/data/pbmc3k_raw.h5ad"
    log(f"marker panel  {len(PBMC_MARKERS)} cell types")
    if h5ad.is_file():
        log(f"dataset       {h5ad.name} ({h5ad.stat().st_size/1e6:.1f} MB), local copy")
        adata = sc.read_h5ad(h5ad)
    else:
        log("dataset       PBMC 3k via scanpy (local copy absent — data/ is not published)")
        adata = sc.datasets.pbmc3k()
    log(f"loaded        {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    log("running       QC -> normalize -> HVG -> PCA -> neighbors -> Leiden -> marker DE")
    result = SingleCellAnalysis().run(adata, resolution=1.0)
    log(f"clusters      {result.get('n_clusters', '?')}")
    for ct in result.get("cell_types", []):
        log(f"  cell type   {ct}")
    return result



def run_genomic_foundation(log):
    """E1 — the variant store, Ts/Tv QC and the ACMG secondary-findings panel, all local.

    No Parabricks and no GPU: this demo deliberately starts from an already-called VCF, which is
    exactly the honest boundary the catalogue draws. Alignment and variant calling are the gated
    part; everything below runs on a clean clone.
    """
    sys.path.insert(0, str(ROOT / "core/engines/genomic-foundation/src"))
    from variant_store import VariantStore
    import acmg_sf

    # Prefer the real GIAB HG002 genome when the local data checkout is present; fall back to the
    # tracked test fixture so the demo still runs for anyone who just cloned the repo.
    big = ROOT / "hcls-ai-factory-core-data/vcf/HG002.genome.vcf.gz"
    fixture = ROOT / "core/engines/genomic-foundation/tests/fixtures/good_qc.vcf"
    src, limit = (big, 200_000) if big.is_file() else (fixture, None)
    log(f"source        {src.name}"
        f"{' (GIAB HG002, publicly consented — never a patient)' if src is big else ' (test fixture)'}")

    store = VariantStore()
    n = store.load_vcf(src, sample="HG002", limit=limit)
    scope = f"first {n:,} records" if limit else f"all {n:,} records"
    log(f"loaded        {scope} into DuckDB")
    # HG002 ships as a gVCF, so most records are non-variant reference blocks -- a low "pass
    # rate" here is that, not poor calling.
    log(f"PASS          {store.n_pass():,} ({store.pass_rate()*100:.1f}% — the rest are gVCF "
        "reference blocks, not failures)")
    tstv = store.ts_tv()
    log(f"Ts/Tv         {tstv:.3f}")
    if limit:
        log("              NB: computed over a leading slice, not the whole genome — it is the "
            "QC signal working, not a genome-wide figure (that is ~2.0-2.1)")
    if not 1.5 <= tstv <= 3.0:
        raise RuntimeError(f"Ts/Tv {tstv:.3f} is outside any plausible range — QC signal is wrong")

    panel = acmg_sf.panel_summary()
    log(f"ACMG SF panel {panel['n_genes']} genes, {panel['n_conditions']} conditions")
    log(f"              {panel['version']}")
    # A reportable finding and a deliberately non-reportable one, to show the filter discriminates.
    probe = [
        {"gene": "BRCA1", "clinical_significance": "Pathogenic"},
        {"gene": "BRCA1", "clinical_significance": "Benign"},
        {"gene": "TTN", "clinical_significance": "Pathogenic"},
    ]
    found = acmg_sf.secondary_findings(probe)
    log(f"SF filter     {len(found)} of {len(probe)} reportable")
    for v in found:
        log(f"  reportable  {v['gene']} {v['clinical_significance']} -> {v['acmg_sf_condition']}")
    if len(found) != 1:
        raise RuntimeError("ACMG SF filter should report exactly the pathogenic on-panel variant")
    log("decision support for a qualified clinician, not diagnosis")


def run_structural_biology(log):
    """E7 — developability scoring and the guided single-point optimiser, on CPU.

    ESMFold and the CUDA path are gated; these biophysical proxies are not, and they are the part
    the registry records as verified (E31K lowers instability 36.0 -> 28.2). This reproduces that.
    """
    sys.path.insert(0, str(ROOT / "core/engines/structural-biology/src"))
    from developability import develop_metrics, develop_flags, DevelopabilityScorer

    # Human lysozyme C (P61626) mature chain — a real, well-characterised sequence.
    seq = ("KVFERCELARTLKRLGMDGYRGISLANWMCLAKWESGYNTRATNYNAGDRSTDYGIFQINSRYWCNDGKTPGAVNACHLSCSALLQDNIADAVACAKRVVRDPQGIRAWVAWRNRCQNRDVRQYVQGCGV")
    m = develop_metrics(seq)
    flags, verdict = develop_flags(m)
    log(f"sequence      human lysozyme C, {m['length']} aa, {m['molecular_weight']:.0f} Da")
    log(f"GRAVY         {m['gravy']}      instability {m['instability_index']}")
    log(f"pI            {m['isoelectric_point']}   aromaticity {m['aromaticity']}")
    log(f"verdict       {verdict}" + (f" — {'; '.join(flags)}" if flags else ""))

    # Call the engine's own optimiser rather than hand-picking a substitution. It scans every
    # position against a tolerated-substitution set and keeps only those that actually lower the
    # instability index -- so "improvement" is computed, never assumed. (An arbitrary E->K does
    # not improve this sequence: +1.35. That is the difference between a guided optimiser and a
    # guess, and it is why this demo asks the optimiser instead of asserting a result.)
    opt = DevelopabilityScorer().optimize(seq, n=3)
    log(f"optimiser     scanned {m['length']} positions x 10 substitutions — "
        f"{opt['n_proposals']} lower the instability index")
    if not opt["top"]:
        raise RuntimeError("optimiser proposed no improving substitution — its premise failed")
    for p in opt["top"]:
        log(f"  proposal    {p['mutation']}  instability {opt['baseline_instability']} -> "
            f"{p['instability_index']}  ({p['delta']:+.2f})")
    best = opt["top"][0]
    if best["delta"] >= 0:
        raise RuntimeError(f"top proposal {best['mutation']} does not improve instability")

    # Say out loud what the metric cannot see. The Guruprasad instability index is a
    # sequence-only proxy: it has no concept of disulfide bonding, so it will happily propose
    # substituting a structural cysteine. Lysozyme has four disulfide bridges, and the top
    # proposals here target one of them -- a change that improves the number and would
    # destabilise the actual protein. A developability screen narrows a design space; it does
    # not rank designs on its own.
    if best["mutation"].startswith("C"):
        log(f"caveat        top proposal substitutes a cysteine ({best['mutation']}); the "
            "instability index is sequence-only and cannot see disulfide bonds")
    log("preclinical — a research bench, not a therapeutic claim; developability proxies")
    log("              narrow a design space, they do not rank designs")



def run_tsc_program(log):
    """P1 — the flagship disease program: the whole factory composed for one child.

    Walks the surfaces a clinician would actually see, in order: the cohort, one featured
    patient's briefing, the provenance chain behind it, and the five disease-specific agents
    that produced it. Every record is SYNTHETIC and the engine says so in its own health
    payload -- this demo asserts that watermark rather than trusting it, because a disease
    program that quietly served real patient data would be the worst failure in the project.
    """
    def get(path):
        with urllib.request.urlopen(f"http://localhost:8560{path}", timeout=60) as r:
            return json.loads(r.read().decode())

    health = get("/health")
    if health.get("watermark") != "SYNTHETIC":
        raise RuntimeError(
            f"TSC engine did not declare SYNTHETIC data (watermark={health.get('watermark')!r}) — "
            "refusing to run a disease-program demo that may be serving real patient records")
    log(f"engine        {health.get('engine')}")
    log(f"data          {health['watermark']} — synthetic cohort, never a real patient")

    cohort = get("/cohort")
    n = cohort.get("n_patients", 0)
    if not n:
        raise RuntimeError("cohort is empty — the synthetic cohort has not been built")
    log(f"cohort        {n} patients")
    for k, v in list((cohort.get("distributions") or {}).get("classification", {}).items())[:4]:
        log(f"  {k:<34} {v}")

    agents = get("/agents")
    log(f"agents        {len(agents)} disease-specific agents, dependency-ordered")
    for a in agents:
        log(f"  {a.get('name'):<22} emits {a.get('emits')}")

    pid = (health.get("featured") or {}).get("A") or "TSC-0043"
    brief = get(f"/surfaces/briefing/{pid}")
    hdr = brief.get("header") or {}
    log(f"patient       {pid} — {hdr.get('genotype')} {hdr.get('variant')} "
        f"({hdr.get('classification', '?')})")

    prov = get(f"/provenance/{pid}")
    if not prov:
        raise RuntimeError(f"no provenance for {pid} — the audit chain is empty")
    log(f"provenance    {len(prov)} recorded events behind that briefing")
    for e in prov[:3]:
        rec = (e.get("records") or [{}])[0]
        log(f"  {e.get('event'):<22} by {rec.get('agent')} v{rec.get('agent_version')}")

    log("TSC1/TSC2 gene therapy is preclinical; every output is decision support for a")
    log("              qualified clinician, behind a review gate — never diagnosis")



def run_precision_intelligence(log):
    """E2 — the evidence layer: annotated variants to ranked druggable targets.

    This is the hand-off between Engine 1 and Engine 3: E1 produces the QC'd variant substrate,
    E2 annotates and reasons over it, and what comes out is the target E3 designs against. The
    portal serves that as a target register, so the demo reads it rather than asking an LLM to
    narrate -- the targets are stored decisions with a mechanism and a provenance, not generated
    prose.
    """
    def get(path):
        with urllib.request.urlopen(f"http://localhost:5001{path}", timeout=90) as r:
            return json.loads(r.read().decode())

    ready = get("/api/ready")
    checks = ready.get("checks", {})
    log(f"readiness     milvus={checks.get('milvus')} "
        f"collection_loaded={checks.get('collection_loaded')} llm={checks.get('ollama')}")
    if not checks.get("milvus") or not checks.get("collection_loaded"):
        raise RuntimeError(f"evidence layer not ready: {checks}")

    status = get("/api/status")
    coll = status.get("collection", {})
    data = status.get("data", {})
    n = coll.get("num_entities", 0)
    if not n:
        raise RuntimeError(
            f"'{coll.get('name')}' holds no vectors — the clinical evidence base is not loaded")
    log(f"evidence base {coll.get('name')}: {n:,} vectors (ClinVar / AlphaMissense)")
    log(f"variant input {'present' if data.get('vcf_exists') else 'ABSENT'} — "
        f"{data.get('vcf_size')} ({pathlib.Path(str(data.get('vcf_path'))).name})")
    if not data.get("vcf_exists"):
        raise RuntimeError("the VCF the evidence layer annotates is not on disk")

    tg = get("/api/targets")
    targets = tg.get("targets") or []
    if not targets:
        raise RuntimeError("no druggable targets registered — E2 produced nothing for E3")
    summ = tg.get("summary", {})
    log(f"targets       {summ.get('total')} registered · "
        f"{summ.get('by_confidence', {}).get('high', 0)} high-confidence · "
        f"{summ.get('by_status', {}).get('validated', 0)} validated")
    for t in targets[:4]:
        log(f"  {t.get('gene'):<8} {t.get('confidence'):<7} {t.get('status', ''):<11} "
            f"{str(t.get('mechanism', ''))[:70]}")

    vcp = next((t for t in targets if t.get("gene") == "VCP"), None)
    if vcp:
        log(f"flagship      VCP (p97) for frontotemporal dementia — {vcp.get('mechanism')}")
        log(f"              {str(vcp.get('notes', ''))[:110]}")
    log("decision support for a qualified clinician, not diagnosis")



def run_therapeutic_discovery(log):
    """E3 — target to ranked candidate molecules, on the path that runs locally.

    The full ten-stage pipeline ends in DiffDock pose prediction, and the NIMs (MolMIM, DiffDock)
    are gated. What is NOT gated is the part that matters for an R&D bench: fragment-based
    generation over a seed set, and real RDKit chemistry on every candidate. So this runs
    generation + QC + ranking for real and says plainly where the gated boundary is, rather than
    narrating a docking score nobody computed.
    """
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors, QED, Crippen, Lipinski
    RDLogger.DisableLog("rdApp.*")

    # The target comes from E2's register -- this is the actual Engine 2 -> Engine 3 hand-off.
    try:
        with urllib.request.urlopen("http://localhost:5001/api/targets", timeout=60) as r:
            targets = json.loads(r.read().decode()).get("targets") or []
        vcp = next((t for t in targets if t.get("gene") == "VCP"), None)
        if vcp:
            log(f"target        VCP (p97) from E2's register — {vcp.get('mechanism')}")
            log(f"              reference compound CB-5083; FTD (frontotemporal dementia)")
    except Exception:
        log("target        VCP (p97) — E2 register unreachable, using the flagship target")

    seeds = [
        "COc1cc2c(Nc3ccc(Br)cc3F)ncnc2cc1OCC1CCN(C)CC1",
        "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1",
        "CN1CCN(CC1)c1ccc(Nc2ncc(F)c(Nc3ccccc3)n2)cc1",
    ]
    body = json.dumps({"seeds": seeds, "n": 12}).encode()
    req = urllib.request.Request("http://localhost:8574/generate", data=body,
                                 headers={"Content-Type": "application/json"}, method="POST")
    with urllib.request.urlopen(req, timeout=240) as r:
        gen = json.loads(r.read().decode())
    mols = gen.get("molecules") or []
    if not mols:
        raise RuntimeError("molecule generator returned nothing — BRICS needs a seed set it can "
                           "fragment; check :8574")
    log(f"generation    {gen.get('backend')} backend, {len(seeds)} seeds -> {len(mols)} candidates")

    # Real chemistry on every candidate -- not a score copied from the generator.
    rows = []
    for m in mols:
        smi = m.get("smiles") if isinstance(m, dict) else m
        mol = Chem.MolFromSmiles(smi or "")
        if mol is None:
            continue
        mw, logp = Descriptors.MolWt(mol), Crippen.MolLogP(mol)
        hbd, hba = Lipinski.NumHDonors(mol), Lipinski.NumHAcceptors(mol)
        violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
        rows.append({"smiles": smi, "qed": round(QED.qed(mol), 3), "mw": round(mw, 1),
                     "logp": round(logp, 2), "ro5": violations})
    if not rows:
        raise RuntimeError("no generated molecule survived RDKit parsing — chemistry QC failed")

    passed = [r for r in rows if r["ro5"] == 0]
    log(f"chemistry QC  {len(rows)} parsed by RDKit · {len(passed)} pass Lipinski Ro5 "
        f"(0 violations)")
    rows.sort(key=lambda r: -r["qed"])
    log("ranking       by QED (drug-likeness)")
    for r in rows[:4]:
        log(f"  QED {r['qed']:<6} MW {r['mw']:<7} logP {r['logp']:<6} Ro5×{r['ro5']}  {r['smiles'][:46]}")

    log("GATED — not run here: MolMIM generation (NIM), DiffDock pose prediction (NIM),")
    log("              chemprop ADMET. No binding affinity is claimed; these are drug-LIKENESS")
    log("              scores on generated chemistry, not activity against p97.")
    log("preclinical — a research bench, not a therapeutic claim")


RUNNERS = {"single_cell": run_single_cell,
           "therapeutic_discovery": run_therapeutic_discovery,
           "precision_intelligence": run_precision_intelligence,
           "tsc_program": run_tsc_program,
           "genomic_foundation": run_genomic_foundation,
           "structural_biology": run_structural_biology}


def execute(demo, verbose=True):
    TRANSCRIPTS.mkdir(parents=True, exist_ok=True)
    lines = []

    def log(msg):
        lines.append(msg)
        if verbose:
            print(f"    {msg}")

    stamp = datetime.now(timezone.utc).isoformat(timespec="seconds")
    log(f"demo          {demo.key} · {demo.subject}")
    log(f"title         {demo.title}")
    log(f"label         {demo.label}")
    log(f"started       {stamp}")

    ok, reasons = demo.check()
    for r in reasons:
        log(f"prerequisite  {r}")
    if not ok:
        log("RESULT        BLOCKED — prerequisites not met")
        (TRANSCRIPTS / f"{demo.key}.txt").write_text("\n".join(lines) + "\n")
        return False

    fn = RUNNERS.get(demo.runner) if demo.runner else (http_demo(demo) if demo.endpoint else None)
    if fn:
        try:
            fn(log)
            log("RESULT        PASS — ran on real input")
        except Exception as e:  # noqa: BLE001
            log(f"RESULT        FAIL — {type(e).__name__}: {e}")
            (TRANSCRIPTS / f"{demo.key}.txt").write_text("\n".join(lines) + "\n")
            return False
    else:
        log("RESULT        NOT IMPLEMENTED — spec only, see docs/demos/DEMO_CATALOG.md")
        (TRANSCRIPTS / f"{demo.key}.txt").write_text("\n".join(lines) + "\n")
        return False

    (TRANSCRIPTS / f"{demo.key}.txt").write_text("\n".join(lines) + "\n")
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("key", nargs="?")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--check-all", action="store_true")
    a = ap.parse_args()

    if a.list or (not a.key and not a.check_all):
        print(f"{'key':5s}{'subject':26s}{'label':16s}title")
        for d in DEMOS:
            print(f"  {d.key:3s}{d.subject:26s}{d.label:16s}{d.title}")
        return 0

    if a.check_all:
        print(f"{'key':5s}{'label':16s}{'ready':7s}why")
        ready = 0
        for d in DEMOS:
            ok, reasons = d.check()
            ready += ok
            why = "; ".join(r for r in reasons if not r.startswith("gated")) or "-"
            print(f"  {d.key:3s}{d.label:16s}{'yes' if ok else 'NO':7s}{why}")
        print(f"\n  {ready}/{len(DEMOS)} demos have their prerequisites met")
        return 0

    d = BY_KEY.get(a.key.upper())
    if not d:
        print(f"unknown demo '{a.key}' — try --list")
        return 2
    print(f"\n{d.key} · {d.title}\n")
    return 0 if execute(d) else 1


if __name__ == "__main__":
    raise SystemExit(main())
