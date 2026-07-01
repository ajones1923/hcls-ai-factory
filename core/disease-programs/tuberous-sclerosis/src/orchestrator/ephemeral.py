"""
Ephemeral state (PRD §2.5.3). Redis holds nothing durable: run locks, in-flight sets,
short-TTL surface caches, and the alert sliding-window counter. A Redis backend and an
in-process backend share one interface; `make_ephemeral()` picks Redis when TSC_USE_REDIS=1,
else the in-memory backend (the default, runs with no server).
"""
from __future__ import annotations

import json
import time
from abc import ABC, abstractmethod

from config.settings import settings


class Ephemeral(ABC):
    @abstractmethod
    def acquire_lock(self, patient_id: str, agent: str, ttl: int = 30) -> bool: ...
    @abstractmethod
    def set_inflight(self, patient_id: str, agent: str, on: bool) -> None: ...
    @abstractmethod
    def inflight(self, patient_id: str) -> list[str]: ...
    @abstractmethod
    def cache_surface(self, kind: str, patient_id: str, payload: dict, ttl: int = 60) -> None: ...
    @abstractmethod
    def get_surface(self, kind: str, patient_id: str) -> dict | None: ...
    @abstractmethod
    def alert_window_incr(self, clinician_id: str, window_s: int = 604800) -> int: ...
    @abstractmethod
    def alert_window_count(self, clinician_id: str, window_s: int = 604800) -> int: ...


class InMemoryEphemeral(Ephemeral):
    def __init__(self) -> None:
        self._kv: dict[str, tuple[float, str]] = {}      # key -> (expires_at, json)
        self._inflight: dict[str, set[str]] = {}
        self._alerts: dict[str, list[float]] = {}

    def _live(self, key: str):
        v = self._kv.get(key)
        if v and v[0] >= time.time():
            return v[1]
        self._kv.pop(key, None)
        return None

    def acquire_lock(self, patient_id, agent, ttl=30):
        key = f"lock:agent:{patient_id}:{agent}"
        if self._live(key) is not None:
            return False
        self._kv[key] = (time.time() + ttl, "1")
        return True

    def set_inflight(self, patient_id, agent, on):
        s = self._inflight.setdefault(patient_id, set())
        s.add(agent) if on else s.discard(agent)

    def inflight(self, patient_id):
        return sorted(self._inflight.get(patient_id, set()))

    def cache_surface(self, kind, patient_id, payload, ttl=60):
        self._kv[f"surface:cache:{kind}:{patient_id}"] = (time.time() + ttl, json.dumps(payload, default=str))

    def get_surface(self, kind, patient_id):
        v = self._live(f"surface:cache:{kind}:{patient_id}")
        return json.loads(v) if v else None

    def alert_window_incr(self, clinician_id, window_s=604800):
        now = time.time()
        hits = [t for t in self._alerts.get(clinician_id, []) if t >= now - window_s]
        hits.append(now)
        self._alerts[clinician_id] = hits
        return len(hits)

    def alert_window_count(self, clinician_id, window_s=604800):
        now = time.time()
        return len([t for t in self._alerts.get(clinician_id, []) if t >= now - window_s])


class RedisEphemeral(Ephemeral):  # pragma: no cover - needs a Redis server
    def __init__(self, url: str | None = None) -> None:
        import redis  # noqa: WPS433

        self.r = redis.Redis.from_url(url or settings.REDIS_URL, decode_responses=True)

    def acquire_lock(self, patient_id, agent, ttl=30):
        return bool(self.r.set(f"lock:agent:{patient_id}:{agent}", "1", nx=True, ex=ttl))

    def set_inflight(self, patient_id, agent, on):
        key = f"inflight:{patient_id}"
        self.r.sadd(key, agent) if on else self.r.srem(key, agent)

    def inflight(self, patient_id):
        return sorted(self.r.smembers(f"inflight:{patient_id}"))

    def cache_surface(self, kind, patient_id, payload, ttl=60):
        self.r.set(f"surface:cache:{kind}:{patient_id}", json.dumps(payload, default=str), ex=ttl)

    def get_surface(self, kind, patient_id):
        v = self.r.get(f"surface:cache:{kind}:{patient_id}")
        return json.loads(v) if v else None

    def alert_window_incr(self, clinician_id, window_s=604800):
        key = f"ratelimit:alerts:{clinician_id}"
        now = time.time()
        self.r.zadd(key, {str(now): now})
        self.r.zremrangebyscore(key, 0, now - window_s)
        self.r.expire(key, window_s)
        return int(self.r.zcard(key))

    def alert_window_count(self, clinician_id, window_s=604800):
        key = f"ratelimit:alerts:{clinician_id}"
        now = time.time()
        self.r.zremrangebyscore(key, 0, now - window_s)
        return int(self.r.zcard(key))


def make_ephemeral() -> Ephemeral:
    import os

    if os.environ.get("TSC_USE_REDIS") == "1":
        return RedisEphemeral()
    return InMemoryEphemeral()
