# Workbook — Hardening the HCLS AI Factory

Companion to [HARDENING_PRD.md](HARDENING_PRD.md). Work it in order; each step carries the command
and the check that settles it. Where a step has already cost time on this machine, the trap is
written next to it.

```bash
cd ~/projects/hcls-ai-factory
export PY=.venv/bin/python
```

---

## Phase 0 — Close the door  (~1 hour)

### 0.1 See the problem yourself first

```bash
ip -4 addr show | grep -oP 'inet \K[\d.]+' | grep -v 127.0.0.1   # your LAN address
LAN=$(ip -4 addr show | grep -oP 'inet \K[\d.]+' | grep -v '^127\.' | head -1)

curl -s -m 6 -o /dev/null -w "health: %{http_code}\n" http://$LAN:8522/health
curl -s -m 20 -o /dev/null -w "clinical POST: %{http_code}\n" \
  -X POST -H 'Content-Type: application/json' -d '{"question":"test"}' http://$LAN:8522/search
```

Both return **200** today. Anyone on the network can ask a clinical agent a question, get generated
prose with dosing and management guidance, and spend your API key doing it.

### 0.2 Answer H-D1 before changing anything

Three options; pick one and write it in the PRD §7:

| Option | Do this | Cost |
|---|---|---|
| **Localhost only** | services bind `127.0.0.1` | portal unreachable from other machines |
| **Firewall the range** | deny 8500–8600 except from trusted hosts | demos still work on-box |
| **Keep exposure, require the key** | turn on `HCLS_API_KEY` | every client must send a header |

Exposed **and** unauthenticated is the one combination nobody would choose deliberately. Anything
else is fine if it is chosen.

### 0.3 Turn on the API key

```bash
sed -i 's/^# HCLS_API_KEY=/HCLS_API_KEY=/' .env
grep -c '^HCLS_API_KEY=' .env        # expect 1
./health-monitor.sh restart all      # services read it at startup
```

Verify fail-closed on a real endpoint:

```bash
curl -s -o /dev/null -w "no key:    %{http_code}\n" -X POST -H 'Content-Type: application/json' \
  -d '{"question":"x"}' http://localhost:8522/search
curl -s -o /dev/null -w "wrong key: %{http_code}\n" -H "X-API-Key: wrong" -X POST \
  -H 'Content-Type: application/json' -d '{"question":"x"}' http://localhost:8522/search
curl -s -o /dev/null -w "right key: %{http_code}\n" -H "X-API-Key: $(grep '^HCLS_API_KEY=' .env | cut -d= -f2)" \
  -X POST -H 'Content-Type: application/json' -d '{"question":"x"}' http://localhost:8522/search
```

Expect **401 / 401 / 200**.

> **Trap — this will break the demos and the eval.** `run_demo.py` and `run_clinical_eval.py` send
> no auth header. Turning the key on without updating them turns 17/17 into 0/17 and 6/6 into six
> ERRORs. Update both to read `HCLS_API_KEY` from the environment **in the same change**, then:
>
> ```bash
> for k in E4 A1 A6; do $PY scripts/run_demo.py $k >/dev/null && echo "$k ok"; done
> $PY scripts/run_clinical_eval.py | tail -2
> ```

### 0.4 Firewall the service range (if that was the decision)

```bash
sudo ufw status verbose                      # ufw reports active but is not blocking these ports
sudo ufw deny proto tcp from any to any port 8500:8600
sudo ufw reload
curl -s -m 6 -o /dev/null -w "from LAN: %{http_code}\n" http://$LAN:8522/health   # expect 000
curl -s -m 6 -o /dev/null -w "on box:  %{http_code}\n" http://localhost:8522/health # expect 200
```

Leave 80/443 alone — Caddy and the site ride on those.

### 0.5 Fix the license so GitHub can see it

```bash
gh repo view ajones1923/hcls-ai-factory --json licenseInfo -q '.licenseInfo.spdxId // "none"'
# -> none

curl -sL https://www.apache.org/licenses/LICENSE-2.0.txt -o /tmp/apache.txt
diff <(sed -n '1,190p' LICENSE) <(sed -n '1,190p' /tmp/apache.txt) | head -20
```

The body matches; the file is missing the **APPENDIX** and the copyright line placement that
`licensee` keys on. Replace with the canonical text, then re-apply the attribution:

```bash
cp LICENSE /tmp/LICENSE.bak
cp /tmp/apache.txt LICENSE
# then put "Copyright 2026 Adam Jones" into the APPENDIX placeholder at the end
git diff --stat LICENSE
```

Push, wait a minute, then:

```bash
gh repo view ajones1923/hcls-ai-factory --json licenseInfo -q '.licenseInfo.spdxId'   # Apache-2.0
```

### 🚦 Phase 0 gate

- LAN clinical POST refused (401 or 000)
- 17/17 demos and 6/6 eval still pass **with** the key on
- GitHub reports `Apache-2.0`

---

## Phase 1 — Supply chain  (~half a day)

### 1.1 Least-privilege CI tokens

Add to the top of `.github/workflows/ci.yml`, above `jobs:`:

```yaml
permissions:
  contents: read
```

No job in this workflow writes to the repo. Verify the next run still passes all six.

### 1.2 Pin actions by SHA

```bash
grep -ohE 'uses: [^ ]+' .github/workflows/*.yml | sort -u
gh api repos/actions/checkout/git/ref/tags/v4 -q .object.sha
gh api repos/actions/setup-python/git/ref/tags/v5 -q .object.sha
```

Replace `@v4` / `@v5` with the SHA, keeping the tag in a trailing comment.

### 1.3 Patch the platform venv

```bash
$PY -m pip install --upgrade pip accelerate
.venv/bin/pip-audit --local 2>/dev/null | tail -5      # expect clean
```

> Only **2** packages in the platform venv are affected. A bare `pip-audit` reports ~40 because it
> audits the OS Python too (`ufw`, `python-apt`, `louis`) — those are Ubuntu's, not the project's.
> Do not let that number drive a panic.

### 1.4 Triage the 14 dependabot PRs

```bash
gh pr list --state open --limit 20
```

Suggested handling:

| Group | Action |
|---|---|
| `github-actions` bumps | merge after 1.2 |
| minor/patch python group (#132) | merge, then run the gate |
| `anthropic >= 0.116` | merge — the SDK in use is already newer than the pin |
| `streamlit >= 1.58` | merge, then confirm the 4 UIs still serve |
| **`pymilvus 2.x → 3.0.0`** (#45, #49, #52, #54) | **hold — H-D3** |

> **Trap — pymilvus 3.0 is a major bump.** 2.6.8 is installed and working, and the agents'
> retrieval was repaired only days ago by fixing a client-API mismatch (`_CollectionManager`
> vs `MilvusClient`, `param=` vs `search_params=`, Hit-objects vs dicts). A second client-API
> change would look identical and be very hard to attribute. If you merge these, run
> `$PY scripts/run_clinical_eval.py` immediately after — it is the only detector.

After each merge:

```bash
ruff check --select E9,F82,F811,F706,F707 core lib scripts
( cd lib/hcls_common && ../../$PY -m pytest -q )
$PY scripts/run_all_tests.py | tail -2
$PY scripts/run_clinical_eval.py | tail -2
```

### 1.5 Make CI actually gate

```bash
gh api -X PUT repos/ajones1923/hcls-ai-factory/branches/main/protection \
  -f 'required_status_checks[strict]=true' \
  -f 'required_status_checks[contexts][]=Lint (real-bug rules)' \
  -f 'required_status_checks[contexts][]=Subject suites (8 engines, 8 agents, 1 program)' \
  -f 'required_status_checks[contexts][]=Capability registry validation' \
  -f 'required_status_checks[contexts][]=Docs site (strict build)' \
  -f 'required_pull_request_reviews[required_approving_review_count]=0' \
  -F 'enforce_admins=true' -F 'restrictions=null'
```

Verify by pushing a deliberately failing branch and confirming the merge button is blocked.

### 🚦 Phase 1 gate

`pip-audit --local` clean · fewer than 5 open PRs · a red run cannot reach `main`.

---

## Phase 2 — De-duplicate  (~2–3 days)

**Do not start this until the clinical eval has been expanded.** The eval is the only thing that
will tell you a refactor changed an answer; six cases is thin cover for touching 13,821 lines.

### 2.1 Extract the RAG path

Twelve `rag_engine.py` files, 13,821 LOC, 15–67% similar. The parts that are genuinely shared —
and that broke identically in five services — are:

- building the Milvus client (must be a `MilvusClient`, constructed **inside** the lifespan)
- the search call (`search_params=`, not the ORM's `param=`)
- flattening results (dicts from `MilvusClient`, Hit objects from the ORM)
- normalising a hit for the routes (`_as_row`)

Move those into `hcls_common` — follow `ingest_persist.py`, which documents each trap beside the
code that handles it. Leave the domain-specific weighting and prompts in each agent.

```bash
grep -c 'search_params' core/*/*/src/rag_engine.py core/*/*/*/src/rag_engine.py
$PY scripts/run_clinical_eval.py | tail -2      # after each service is migrated
```

### 2.2 Extract the LLM client

Eight per-service `_LLMClient` classes plus `hcls_common.llm_client`. The `temperature` removal had
to be applied in **nine** places; the next model change will too.

```bash
git grep -l '_no_sampling' -- '*/api/main.py' | wc -l    # 8 today
```

Target: agents import one client; a model or parameter change is one edit.

### 2.3 Stop shadowing the stdlib

```bash
git ls-files '*/src/collections.py' | wc -l     # 11
for f in $(git ls-files '*/src/collections.py'); do
  git mv "$f" "$(dirname $f)/vector_collections.py"
done
grep -rln 'from src.collections\|from .collections\|import collections' core/*/*/src core/*/*/api | head
```

Update the imports, then remove the harness workaround:

```bash
grep -n 'shadows_stdlib\|src/ withheld' scripts/run_all_tests.py
( cd core/agents/cart && ../../../$PY -m pytest -q )   # bare pytest should now work
```

### 🚦 Phase 2 gate

One RAG path · one LLM client · bare `pytest` works in each subject · eval unchanged.

---

## Phase 3 — Polish

### 3.1 Get the video out of git

```bash
git ls-files -z 'docs/assets/videos/*.mp4' | xargs -0 du -ch | tail -1   # 192M of a 213M repo
```

Options: Git LFS, or host the files and reference them. Either way, history still carries them —
a `--depth 1` clone is the metric that matters:

```bash
git clone --depth 1 https://github.com/ajones1923/hcls-ai-factory /tmp/clonetest 2>&1 | tail -2
du -sh /tmp/clonetest && rm -rf /tmp/clonetest
```

### 3.2 Fix the 16 broken doc links

```bash
$PY - <<'PY'
import pathlib, re
root = pathlib.Path("docs"); bad = []
for md in root.rglob("*.md"):
    for m in re.finditer(r"\]\(([^)#:]+\.md)[#)]", md.read_text(errors="ignore")):
        if not (md.parent / m.group(1)).resolve().exists():
            bad.append((str(md.relative_to(root)), m.group(1)))
for a, b in bad: print(f"{a} -> {b}")
print(len(bad), "links")
PY
```

Most point at `licensing.md`, `DATA_SETUP.md`, `demo-guide.md` — files that were never written.
Write them or drop the links. **24 of the 40 raw hits resolve to build-time generated pages
(`maturity-matrix.md`, `factory/engines/index.md`) and are not broken** — `mkdocs --strict` passing
with 0 warnings is the authority.

### 3.3 Label what is superseded

`docs/HCLS_AI_FACTORY_DEMO_GUIDE.md` still prints pre-convention ports (8529 neurology, 8128
clinical-trial, 8130 single-cell). It is superseded on ports by `docs/demos/DEMO_CATALOG.md` and
`docs/build/PORT_MAP.md`. Add a banner rather than deleting — provenance has value.

### 3.4 Agree the test-depth floor (H-D5)

```bash
$PY scripts/run_all_tests.py --json /tmp/t.json >/dev/null
$PY -c "import json;r=json.load(open('/tmp/t.json'));r.sort(key=lambda x:x['passed']);print([(x['name'],x['passed']) for x in r[:5]])"
```

Spread is 4 → 1,966. The single-cell *engine* has 4 tests and is a deterministic scanpy pipeline —
low count is defensible there. Pick a floor for subjects that emit **clinical** output and record it.

---

## Traps already paid for on this machine

1. **A count is not a cause.** `run_all_tests.py` reported "errors 36" with no traceback and cost a
   full CI cycle. It now prints the failing test ids and the pytest tail.
2. **`MilvusClient` connects eagerly** in `__init__` — construct it inside the lifespan try, never
   at module scope, or every import fails where Milvus is absent.
3. **`exec 9>&-` is load-bearing in three places** in `health-monitor.sh`. A started service
   inherits the flock fd and holds the lock for its own lifetime.
4. **`str(DataType.FLOAT)` is `"10"`, not `"FLOAT"`.** Use `.name`.
5. **Fetch is not persist.** Three ingest scripts logged "N validated" and wrote nothing.
6. **Absence of a name is not absence of behaviour.** `pharmacogenomics` was reported broken
   because it lacked a `_persist` function; it had always persisted by another route.
