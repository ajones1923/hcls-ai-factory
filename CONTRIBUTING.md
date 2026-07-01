# Contributing to HCLS AI Factory

Thank you for your interest in contributing to the HCLS AI Factory. This is a precision medicine platform that takes patient DNA from sequencing to drug candidates in under five hours. Your contributions can directly accelerate therapeutic discovery and improve patient outcomes.

## How to Contribute

### Report Bugs

Open a [GitHub issue](https://github.com/ajones1923/hcls-ai-factory/issues/new?template=bug_report.md) using the **Bug Report** template. Include the agent or component affected, steps to reproduce, and your environment details.

### Suggest Features

Open a [GitHub issue](https://github.com/ajones1923/hcls-ai-factory/issues/new?template=feature_request.md) using the **Feature Request** template. Describe the use case and any clinical impact.

### Submit Code Changes

1. Fork the repository and create a feature branch from `main`.
2. Make your changes following the code standards below.
3. Add or update tests to cover your changes.
4. Open a pull request using the PR template.

## Development Setup

### Prerequisites

- **Python 3.11+**
- **Docker** and **Docker Compose**
- **NVIDIA GPU** with CUDA 12.x (optional; required only for the genomics pipeline)

### Clone and Install

```bash
git clone https://github.com/ajones1923/hcls-ai-factory.git
cd hcls-ai-factory

# Install the pre-commit guard (blocks large files, secrets, and non-neutral material)
./scripts/install-hooks.sh
# optional, for the full framework hooks:  pip install pre-commit && pre-commit install

# Create a virtual environment
python -m venv .venv && source .venv/bin/activate

# Install shared library
pip install -e lib/hcls_common/

# Install agent dependencies (example: clinical trial agent)
cd core/agents/clinical-trial
pip install -r requirements.txt
```

### Run Tests

Each agent has its own test suite:

```bash
cd core/agents/{agent_name}
pytest tests/
```

## Architecture Overview

The platform consists of eight engines and eight intelligence agents, plus disease-programs:

**Engines** (in `core/engines/`):
1. **genomic-foundation** -- Parabricks, DeepVariant, BWA-MEM2 (FASTQ to VCF) + variant store
2. **precision-intelligence** -- Milvus vector DB + Claude for variant interpretation (RAG)
3. **therapeutic-discovery** -- MolMIM, DiffDock, RDKit + real ADMET for compound generation
4. **clinical-imaging** -- DICOM analysis (VISTA-3D / MAISI / VILA-M3)
5. **precision-oncology** -- MTB packets, therapy ranking, trial matching
6. **cardiology** -- clinical workflows + risk calculators
7. **structural-biology** -- ESMFold, ESM-2 search, ProteinMPNN, developability
8. **single-cell** -- scanpy compute → cell-type annotation

**Intelligence Agents** (in `core/agents/`):
- `cart` -- CAR-T therapy design
- `precision-biomarker` -- Biomarker analysis and reporting
- `pharmacogenomics` -- Drug-gene interactions
- `precision-autoimmune` -- Autoimmune disease analysis
- `neurology` -- Neurological variant analysis
- `clinical-trial` -- Trial matching and eligibility
- `rare-disease-diagnostic` -- Rare disease diagnosis
- `single-cell` -- Single-cell transcriptomics

**Disease-programs** (in `core/disease-programs/`): vertical solutions composing the engines and
agents for one condition -- starting with `tuberous-sclerosis`.

**Shared library:** `lib/hcls_common/` (capability registry, MCP tool-surface, workflow composer,
MLOps, governance gates, Milvus client, LLM wrappers, security utilities).

Each agent follows this directory structure:

```
agent_name/
  config/        # settings.py, Pydantic configuration models
  src/           # agent.py, knowledge.py, rag_engine.py, calculators
  api/           # FastAPI routes
  app/           # Streamlit UI
  tests/         # pytest test suite
  data/          # Reference data and knowledge bases
```

## How to Add a New Agent

1. **Create the directory structure** under `core/agents/your-agent-name/` (kebab-case) following the layout above.

2. **Implement core modules:**
   - `config/settings.py` -- Pydantic `Settings` class with agent-specific configuration
   - `src/agent.py` -- Main agent orchestrator
   - `src/knowledge.py` -- Domain knowledge base and embedding logic
   - `src/rag_engine.py` -- RAG retrieval and response generation

3. **Add the API and UI:**
   - `api/main.py` + `api/routes/` -- FastAPI app and endpoint routers
   - `app/{name}_ui.py` -- Streamlit interface for interactive use

4. **Wire in governance** (see "Governance" below) and **register a capability** (see "Register a Capability").

5. **Write tests** in `tests/` covering core logic and API routes; add a `.env.example`.

6. **Register the service** in `docker-compose.dgx-spark.yml` with a port, volume mounts, and health check.

## How to Add an Engine

Engines live under `core/engines/your-engine-name/` (kebab-case) and follow the same
self-contained layout as agents (`src/`, `api/`, `tests/`, `config/`, `README.md`,
`requirements.txt`, `Dockerfile`). Use the `cardiology` engine as the canonical template.

1. **Scaffold** the directory from the cardiology layout.
2. **Build a pre-governed FastAPI app** so the gates are inherited by construction:
   ```python
   from hcls_common.api_gate import create_governed_app
   app = create_governed_app("your-engine", capability_id="your-engine-id")
   ```
3. **Register a capability** for the engine (below) — CI will fail if the new engine
   directory has no registered capability.
4. **Add a `README.md`** describing the engine, its port, and how to run it.
5. **Add tests** and wire the service into `docker-compose.dgx-spark.yml`.

## Register a Capability

Every engine/agent must appear in the capability registry so it is discoverable, wireable
by the workflow composer, and governed. This is enforced in CI.

1. **Add an entry** to `lib/hcls_common/capabilities.json` with a unique `id`, a `type`
   (`engine`/`agent`/`model`/`nim`/`stage`/`service`), typed `inputs`/`outputs` ports, an
   `endpoint`, and a `status` (`live`/`planned`). A `live` capability may **never** be
   `mock`-served — the registry rejects it.
2. **Map the directory** in `scripts/validate_registry.py` (`COVERAGE`) to the capability id(s).
3. **Validate:** `python scripts/validate_registry.py` (also runs in CI) must print `OK`.

## Governance

The platform's input-validation and output-honesty gates live in `hcls_common.api_gate`.
Wire them into every service so governance is real, not optional:

```python
# In api/main.py — one line adds the /governance surface + request-id/timing:
from hcls_common.api_gate import install_governance
install_governance(app, service="cart", capability_id="cart-intelligence-agent")

# In POST handlers:
from hcls_common.api_gate import require_valid_input, honesty_flags
payload = require_valid_input("cart-intelligence-agent", payload)   # 422 on bad input
flags = honesty_flags(answer_text)                                  # deterministic overclaim scan
```

For send-ready clinical text, use `assert_publishable(text, llm=...)` to run the full
verify gate. `cart` is the reference implementation.

## Code Standards

### Python

- Target **Python 3.11+**.
- Use **type hints** on all function signatures.
- Use **Pydantic models** for configuration, request/response schemas, and data validation.
- Use **loguru** for logging (not the stdlib `logging` module).
- Keep functions focused; prefer composition over inheritance.

### Dependencies

- New code and the shared library follow the baseline in the root
  [`constraints.txt`](constraints.txt). Install with it applied:
  `pip install -r requirements.txt -c ../../constraints.txt`.
- Prefer compatible-release ranges (`>=X,<X+1`) in `requirements.txt`; keep exact pins in
  the constraints/lockfile, not scattered across components.
- Gate heavy/GPU stacks (torch, monai, scanpy, bionemo) behind optional extras where
  possible, as `hcls_common` does — don't make them hard requirements of light services.

### Testing

- All new features and bug fixes must include tests.
- Use `pytest` as the test runner.
- Mock external services (Milvus, Claude API, NVIDIA NIMs) in unit tests.

### Clinical Accuracy

- All clinical calculators and scoring functions **must cite published sources** in docstrings (e.g., PMID, DOI, or guideline name).
- Variant classifications must reference ClinVar, ACMG, or equivalent standards.
- Never present AI-generated results as clinically validated without appropriate disclaimers.

### Security

- **No secrets in code.** API keys, credentials, and tokens go in environment variables or `.env` files (which are git-ignored).
- Sanitize all error messages returned to users; never expose stack traces or internal paths.
- Validate and sanitize all user inputs to prevent injection attacks.

## Pull Request Process

1. Ensure your branch is up to date with `main`.
2. Verify all tests pass locally: `pytest tests/`.
3. Fill out the pull request template completely.
4. A maintainer will review your PR. Address any feedback promptly.
5. Once approved, a maintainer will merge your PR.

PRs that include clinical calculations or modify agent reasoning logic require review from at least two maintainers.

## Code of Conduct

This project follows the [Contributor Covenant Code of Conduct](CODE_OF_CONDUCT.md). By participating, you agree to uphold a welcoming, inclusive, and respectful community.

## License

This project is licensed under the **Apache License 2.0**. See [LICENSE](LICENSE) for details. By contributing, you agree that your contributions will be licensed under the same terms.
