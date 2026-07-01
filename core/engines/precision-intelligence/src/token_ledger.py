"""
HCLS token-usage ledger (COST-2, 2026-06-29).

A tiny, dependency-free SQLite ledger that records every Anthropic API call's token
usage + computed cost. Written by llm_client.AnthropicClient; read directly off-disk by
the Projects Monitor's collectors.get_claude_api_usage() (no HTTP endpoint, no service
restart needed). This closes the monitor's stubbed `"hcls": _empty_api_usage()`.

Pricing ($ / 1M tokens): cache reads bill ~0.1x input, cache writes ~1.25x input.
"""
from __future__ import annotations

import os
import sqlite3
import time
from pathlib import Path

_DB = Path(os.environ.get(
    "HCLS_TOKEN_LEDGER",
    str(Path(__file__).resolve().parents[1] / "data" / "token_usage.db"),
))

# $ per 1M tokens (input, output); matched by substring of the model id, first match wins.
_PRICING = {
    "opus": (15.0, 75.0),
    "sonnet": (3.0, 15.0),
    "haiku": (1.0, 5.0),
}


def _rates(model: str) -> tuple[float, float]:
    for key, price in _PRICING.items():
        if key in (model or "").lower():
            return price
    return (3.0, 15.0)  # default to Sonnet-tier if unknown


def _connect() -> sqlite3.Connection:
    _DB.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(_DB), timeout=5)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute(
        """CREATE TABLE IF NOT EXISTS token_usage (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            ts REAL NOT NULL,
            model TEXT NOT NULL,
            input_tokens INTEGER NOT NULL,
            output_tokens INTEGER NOT NULL,
            cache_read INTEGER NOT NULL DEFAULT 0,
            cache_create INTEGER NOT NULL DEFAULT 0,
            cost REAL NOT NULL
        )"""
    )
    return conn


def record(model: str, usage) -> None:
    """Record one call's usage. `usage` is the Anthropic Message.usage object or a dict.
    Never raises into the caller — a ledger failure must not break generation."""
    try:
        def _g(k):
            return getattr(usage, k, None) if not isinstance(usage, dict) else usage.get(k)
        inp = int(_g("input_tokens") or 0)
        out = int(_g("output_tokens") or 0)
        cr = int(_g("cache_read_input_tokens") or 0)
        cc = int(_g("cache_creation_input_tokens") or 0)
        in_rate, out_rate = _rates(model)
        cost = (inp * in_rate + cr * in_rate * 0.1 + cc * in_rate * 1.25
                + out * out_rate) / 1_000_000.0
        with _connect() as conn:
            conn.execute(
                "INSERT INTO token_usage (ts, model, input_tokens, output_tokens, "
                "cache_read, cache_create, cost) VALUES (?,?,?,?,?,?,?)",
                (time.time(), model, inp, out, cr, cc, round(cost, 6)),
            )
    except Exception:
        pass  # ledger is best-effort; never disrupt the LLM call


def totals(db_path: str | os.PathLike | None = None) -> dict:
    """Aggregate the ledger into the Projects-Monitor api-usage shape. Reads off-disk;
    safe to call from another process. Returns the empty shape if the DB is absent."""
    path = Path(db_path) if db_path else _DB
    empty = {
        "total_tokens": 0, "total_cost": 0.0, "request_count": 0,
        "by_module": [], "by_model": [],
        "today_tokens": 0, "today_cost": 0.0, "month_tokens": 0, "month_cost": 0.0,
        "source": "sqlite-empty",
    }
    if not path.exists():
        return empty
    try:
        conn = sqlite3.connect(f"file:{path}?mode=ro", uri=True, timeout=5)
        now = time.time()
        day_ago, month_ago = now - 86400, now - 30 * 86400

        def scalar(sql, args=()):
            return conn.execute(sql, args).fetchone()

        tt, tc, rc = scalar(
            "SELECT COALESCE(SUM(input_tokens+output_tokens+cache_read+cache_create),0),"
            " COALESCE(SUM(cost),0), COUNT(*) FROM token_usage")
        td_tok, td_cost = scalar(
            "SELECT COALESCE(SUM(input_tokens+output_tokens+cache_read+cache_create),0),"
            " COALESCE(SUM(cost),0) FROM token_usage WHERE ts>=?", (day_ago,))
        mo_tok, mo_cost = scalar(
            "SELECT COALESCE(SUM(input_tokens+output_tokens+cache_read+cache_create),0),"
            " COALESCE(SUM(cost),0) FROM token_usage WHERE ts>=?", (month_ago,))
        by_model = [
            {"model": m, "total_tokens": int(tok), "total_cost": round(cost, 4),
             "request_count": int(n)}
            for (m, tok, cost, n) in conn.execute(
                "SELECT model, SUM(input_tokens+output_tokens+cache_read+cache_create),"
                " SUM(cost), COUNT(*) FROM token_usage GROUP BY model ORDER BY 2 DESC")
        ]
        conn.close()
        return {
            "total_tokens": int(tt), "total_cost": round(tc, 4), "request_count": int(rc),
            "by_module": [{"module": "hcls-rag-chat", "total_tokens": int(tt),
                           "total_cost": round(tc, 4), "request_count": int(rc)}] if rc else [],
            "by_model": by_model,
            "today_tokens": int(td_tok), "today_cost": round(td_cost, 4),
            "month_tokens": int(mo_tok), "month_cost": round(mo_cost, 4),
            "source": "sqlite",
        }
    except Exception:
        return empty
