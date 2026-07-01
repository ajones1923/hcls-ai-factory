"""Per-tier token + cost ledger (src/utils/cost; NFR-COST) and its router wiring."""
import pytest

from src.utils.cost import CostCapExceeded, CostLedger, get_ledger


def test_ledger_counts_and_prices():
    led = CostLedger()
    led.record("sonnet", 1_000_000, 1_000_000, billed=True)   # $3 in + $15 out = $18
    led.record("haiku", 0, 0, billed=False)
    t = led.totals()
    assert t["cost_usd"] == pytest.approx(18.0)
    assert t["calls"] == 2 and t["billed_calls"] == 1


def test_offline_calls_count_tokens_but_never_bill():
    led = CostLedger()
    led.record("opus", 500, 800, billed=False)
    t = led.totals()
    assert t["input_tokens"] == 500 and t["output_tokens"] == 800
    assert t["cost_usd"] == 0.0 and t["billed_calls"] == 0


def test_cap_blocks_billed_overspend():
    led = CostLedger()
    led.set_cap(0.01)
    with pytest.raises(CostCapExceeded):
        led.record("opus", 1_000_000, 1_000_000, billed=True)   # way over $0.01


def test_router_meters_offline_calls():
    # running any agent offline should accumulate counted-but-unbilled token volume
    from app._engine import featured, get_engine
    get_ledger().reset()
    orch, manifest = get_engine()
    orch.store.projection(featured()["B"])   # ensures the engine has run
    t = get_ledger().totals()
    assert t["calls"] > 0
    assert t["cost_usd"] == 0.0          # offline: never billed
