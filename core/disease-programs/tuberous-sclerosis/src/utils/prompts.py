"""
Versioned prompt-template loader (PRD §2.3: each agent has a prompts/ directory).
Templates live at src/agents/<agent>/prompts/<name>.<version>.txt and are the system
instructions for the agent's LLM steps. Returns "" if a template is missing.
"""
from __future__ import annotations

from pathlib import Path

_AGENTS_DIR = Path(__file__).resolve().parents[1] / "agents"
_CACHE: dict[tuple[str, str, str], str] = {}


def load_prompt(agent: str, name: str, version: str = "v0") -> str:
    key = (agent, name, version)
    if key not in _CACHE:
        path = _AGENTS_DIR / agent / "prompts" / f"{name}.{version}.txt"
        _CACHE[key] = path.read_text().strip() if path.exists() else ""
    return _CACHE[key]
