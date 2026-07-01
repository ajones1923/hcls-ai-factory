"""
Per-tier token + cost ledger (PRD §5.1.3 NFR-COST-1..4; master paper §16).

The model router records every LLM call here. The ledger accumulates input/output tokens
and USD cost per tier, distinguishing **billed** calls (real API usage, online) from
**estimated** ones (offline stubs — token volume counted, cost zero, never charged). It
also enforces an optional per-process spend cap so a runaway loop cannot quietly burn the
key. Prices are approximate published Anthropic rates (USD per 1M tokens) and are labeled
as such; they are the accounting basis, not a quote.

Thread-safe (the FastAPI threadpool runs agents concurrently). A module singleton is the
ledger of record; `get_ledger()` returns it.
"""
from __future__ import annotations

import threading


class CostCapExceeded(RuntimeError):
    """Raised when a billed call would push process spend past the configured cap."""


# Approximate Anthropic list prices, USD per 1,000,000 tokens (input, output). Labeled
# approximate; update as pricing changes. Local/deterministic tiers are free.
PRICES = {
    "haiku": (1.00, 5.00),
    "sonnet": (3.00, 15.00),
    "opus": (15.00, 75.00),
    "local": (0.0, 0.0),
}


def estimate_tokens(text: str) -> int:
    """Cheap offline token estimate (~4 chars/token) for stubbed calls."""
    return max(0, len(text or "") // 4)


def _cost(tier: str, in_tok: int, out_tok: int) -> float:
    pin, pout = PRICES.get(tier, (0.0, 0.0))
    return (in_tok * pin + out_tok * pout) / 1_000_000.0


class CostLedger:
    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._cap_usd: float | None = None
        self._rows: dict[str, dict] = {}

    def set_cap(self, cap_usd: float | None) -> None:
        with self._lock:
            self._cap_usd = cap_usd if cap_usd and cap_usd > 0 else None

    def _row(self, tier: str) -> dict:
        return self._rows.setdefault(
            tier, {"tier": tier, "calls": 0, "input_tokens": 0, "output_tokens": 0,
                   "billed_calls": 0, "cost_usd": 0.0})

    def record(self, tier: str, input_tokens: int, output_tokens: int, billed: bool) -> float:
        """Record one call. Returns the USD cost charged (0.0 when not billed). Raises
        CostCapExceeded *before* recording if a billed call would exceed the cap."""
        cost = _cost(tier, input_tokens, output_tokens) if billed else 0.0
        with self._lock:
            if billed and self._cap_usd is not None:
                if self._billed_total_locked() + cost > self._cap_usd:
                    raise CostCapExceeded(
                        f"billed spend cap ${self._cap_usd:.2f} would be exceeded "
                        f"(current ${self._billed_total_locked():.4f} + ${cost:.4f})")
            r = self._row(tier)
            r["calls"] += 1
            r["input_tokens"] += input_tokens
            r["output_tokens"] += output_tokens
            if billed:
                r["billed_calls"] += 1
                r["cost_usd"] = round(r["cost_usd"] + cost, 6)
        return cost

    def _billed_total_locked(self) -> float:
        return sum(r["cost_usd"] for r in self._rows.values())

    def totals(self) -> dict:
        with self._lock:
            by_tier = [dict(r) for r in self._rows.values()]
            return {
                "by_tier": by_tier,
                "calls": sum(r["calls"] for r in by_tier),
                "input_tokens": sum(r["input_tokens"] for r in by_tier),
                "output_tokens": sum(r["output_tokens"] for r in by_tier),
                "billed_calls": sum(r["billed_calls"] for r in by_tier),
                "cost_usd": round(sum(r["cost_usd"] for r in by_tier), 6),
                "cap_usd": self._cap_usd,
                "pricing_note": "approximate Anthropic list prices (USD/1M tokens); offline "
                                "stub calls count tokens but are never billed",
            }

    def reset(self) -> None:
        with self._lock:
            self._rows.clear()


_ledger: CostLedger | None = None


def get_ledger() -> CostLedger:
    global _ledger
    if _ledger is None:
        _ledger = CostLedger()
    return _ledger
