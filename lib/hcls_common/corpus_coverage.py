"""How much of what a retrieval step meant to search actually answered.

Every multi-collection agent here declares a weight map: for a given workflow, how much of the
retrieval budget each collection is worth. clinical-trial's PATIENT_MATCHING, for instance,
puts 0.25 on `trial_eligibility` and 0.12 on `trial_sites`.

Measured 2026-11-03: **both of those collections are empty**, along with eight more of the
fourteen. Summing the declared weights of the empty ones, **73% of the retrieval budget for
patient matching points at nothing**. The query still succeeds, the answer still arrives, and
`collections_searched` still reports 14 — because a search of an empty collection is not an
error at any layer. It returns zero rows, which looks exactly like a search that found nothing
relevant.

That distinction is the whole point of this module. "I searched and the corpus had nothing for
you" and "I searched four of the fourteen places I was built to look" are different statements,
and only the first one was ever being made.

This is deliberately pure: counts and weights in, a verdict out, no Milvus and no I/O. The
caller owns retrieval; `scripts/check_corpus_weights.py` runs the same arithmetic against the
live instance.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable, Mapping

# Below this share of the declared retrieval budget, an answer is drawn from a corpus so partial
# that saying so matters more than the answer's own confidence score. 0.75 is a judgment call,
# not a published threshold: it lets one minor collection be empty without shouting, and catches
# the case measured above (0.27) by a wide margin.
DEGRADED_BELOW = 0.75

_DISCLOSURE_HEADING = "**Partial corpus.**"


@dataclass(frozen=True)
class CoverageVerdict:
    """What a retrieval step covered, against what it declared it would cover."""

    declared: tuple[str, ...]
    """Every collection the weight map named."""

    contributing: tuple[str, ...]
    """Those that returned at least one row."""

    silent: tuple[str, ...]
    """Those that returned nothing, worst-weighted first — the ones worth naming."""

    weight_declared: float
    weight_contributing: float

    @property
    def coverage(self) -> float:
        """Share of the declared retrieval budget that actually answered, 0.0-1.0.

        Weight, not collection count: ten trivial collections going quiet matters less than
        `trial_eligibility` at 0.25 doing it, and a count cannot tell those apart.
        """
        if self.weight_declared <= 0:
            return 0.0
        return self.weight_contributing / self.weight_declared

    @property
    def is_degraded(self) -> bool:
        return self.coverage < DEGRADED_BELOW

    @property
    def is_complete(self) -> bool:
        return not self.silent

    def disclosure(self, *, max_named: int = 4) -> str:
        """One paragraph for the reader, or "" when every source answered.

        Returned rather than appended by this module: the caller knows whether it is writing
        prose, a JSON field or an SSE frame.
        """
        if self.is_complete:
            return ""
        named = ", ".join(f"`{c}`" for c in self.silent[:max_named])
        more = len(self.silent) - max_named
        if more > 0:
            named += f", and {more} more"
        return (
            f"{_DISCLOSURE_HEADING} {len(self.silent)} of {len(self.declared)} evidence sources "
            f"returned nothing, carrying {1 - self.coverage:.0%} of the retrieval weight this "
            f"workflow declares: {named}. The answer above rests on the remaining "
            f"{self.coverage:.0%}. An empty source and an irrelevant one are indistinguishable "
            f"at query time, so this is a statement about the corpus, not about the question."
        )

    def as_dict(self) -> dict[str, Any]:
        """For a response body, a log line or a metric."""
        return {
            "coverage": round(self.coverage, 4),
            "degraded": self.is_degraded,
            "collections_declared": len(self.declared),
            "collections_contributing": len(self.contributing),
            "silent_collections": list(self.silent),
            "weight_declared": round(self.weight_declared, 4),
            "weight_contributing": round(self.weight_contributing, 4),
        }


def _counts_from(results: Iterable[Any]) -> dict[str, int]:
    """Rows per collection, from search results of any of the three shapes in this repo.

    Objects with a `.collection` attribute (TrialSearchResult and its siblings), dicts carrying
    a `collection` key, and plain `(name, count)` pairs all appear in callers here.
    """
    counts: dict[str, int] = {}
    for item in results:
        if isinstance(item, Mapping):
            name = item.get("collection") or item.get("collection_name")
        elif isinstance(item, tuple) and len(item) == 2:
            name, n = item
            if name:
                counts[str(name)] = counts.get(str(name), 0) + int(n)
            continue
        else:
            name = getattr(item, "collection", None) or getattr(item, "collection_name", None)
        if name:
            counts[str(name)] = counts.get(str(name), 0) + 1
    return counts


def coverage(
    weights: Mapping[str, float],
    results: Iterable[Any] | Mapping[str, int],
) -> CoverageVerdict:
    """Compare a declared weight map against what retrieval returned.

    `results` is either a per-collection row count or the result list itself; the list form is
    what callers already have in hand, so neither needs to instrument its search loop.

    A collection with a weight of 0 is not declared for this workflow and is ignored, present or
    not: a weight map that mentions it with 0 is saying "not this time".
    """
    if isinstance(results, Mapping):
        counts = {str(k): int(v) for k, v in results.items()}
    else:
        counts = _counts_from(results)

    declared = tuple(name for name, w in weights.items() if w > 0)
    total = sum(weights[name] for name in declared)

    contributing = tuple(name for name in declared if counts.get(name, 0) > 0)
    got = sum(weights[name] for name in contributing)

    # Worst first, so a truncated list names the collections that cost the most.
    silent = tuple(sorted(
        (name for name in declared if counts.get(name, 0) <= 0),
        key=lambda name: (-weights[name], name),
    ))

    return CoverageVerdict(
        declared=declared,
        contributing=contributing,
        silent=silent,
        weight_declared=total,
        weight_contributing=got,
    )
