"""This program served clinical prose with no platform governance at all.

Every other service in the factory installs `hcls_common`'s output-honesty gate, which scans
clinician-facing prose for overclaims and appends the decision-support disclaimer. The TSC
program imported nothing from `hcls_common`, so the TAND briefings and therapeutics briefs it
serves went out past both.

`hcls_common` has never been a dependency here, so a failed import warns rather than refusing
to start. What must not happen is being ungoverned *silently*, which is what `/health`
reporting the state is for.
"""
import api.main as main_mod


def test_the_app_installs_the_gate():
    app = main_mod.create_app()
    assert app.state.governance == "hcls_common output-honesty gate"


def test_the_gate_is_actually_on_the_middleware_stack():
    app = main_mod.create_app()
    names = [getattr(m.kwargs.get("dispatch", None), "__name__", "") for m in app.user_middleware]
    assert "_output_honesty" in names
    assert "_governance" in names


def test_a_missing_dependency_warns_instead_of_refusing_to_start(monkeypatch):
    """Taking the flagship program down over a newly added dependency is the worse failure."""
    import builtins

    real_import = builtins.__import__

    def no_hcls(name, *args, **kwargs):
        if name.startswith("hcls_common"):
            raise ModuleNotFoundError("No module named 'hcls_common'")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", no_hcls)
    app = main_mod.create_app()
    assert app.state.governance.startswith("NOT INSTALLED")
    assert "ModuleNotFoundError" in app.state.governance


def test_health_reports_the_governance_state():
    """A reader must be able to tell whether the prose went through the gate."""
    app = main_mod.create_app()
    # /health reads cohort state that the lifespan hook populates; this test is about the
    # governance line, so stand the rest in.
    app.state.manifest = {"n_patients": 0}
    app.state.featured = {}
    health = next(r for r in app.routes if getattr(r, "path", "") == "/health")
    body = health.endpoint()
    assert body["governance"] == "hcls_common output-honesty gate"
