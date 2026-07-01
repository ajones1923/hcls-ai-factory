# Contributing to HCLS AI Factory

Thank you for your interest in contributing to the HCLS AI Factory. This is a precision medicine platform that takes patient DNA from sequencing to drug candidates in under five hours. Your contributions can directly accelerate therapeutic discovery and improve patient outcomes.

## How to Contribute

### Report Bugs

Open a [GitHub issue](https://github.com/your-org/hcls-ai-factory/issues/new?template=bug_report.md) using the **Bug Report** template. Include the agent or component affected, steps to reproduce, and your environment details.

### Suggest Features

Open a [GitHub issue](https://github.com/your-org/hcls-ai-factory/issues/new?template=feature_request.md) using the **Feature Request** template. Describe the use case and any clinical impact.

### Submit Code Changes

1. Fork the repository and create a feature branch from `main`.
2. Make your changes following the code standards below.
3. Add or update tests to cover your changes.
4. Open a pull request using the PR template.

## Development Setup

### Prerequisites

- **Python 3.10+**
- **Docker** and **Docker Compose**
- **NVIDIA GPU** with CUDA 12.x (optional; required only for the genomics pipeline)

### Clone and Install

```bash
git clone https://github.com/your-org/hcls-ai-factory.git
cd hcls-ai-factory

# Create a virtual environment
python -m venv .venv && source .venv/bin/activate

# Install shared library
pip install -e lib/hcls_common/

# Install agent dependencies (example: clinical trial agent)
cd ai_agent_adds/clinical_trial_intelligence_agent
pip install -r requirements.txt
```

### Run Tests

Each agent has its own test suite:

```bash
cd ai_agent_adds/{agent_name}
pytest tests/
```

## Architecture Overview

The platform consists of three core engines and eleven intelligence agents:

**Engines:**
1. **Genomics Pipeline** -- Parabricks, DeepVariant, BWA-MEM2 (FASTQ to VCF)
2. **RAG/Chat Pipeline** -- Milvus vector DB + Claude AI for variant interpretation
3. **Drug Discovery Pipeline** -- BioNeMo MolMIM, DiffDock, RDKit for compound generation

**Intelligence Agents** (in `ai_agent_adds/`):
- `precision_biomarker_agent` -- Biomarker analysis and reporting
- `precision_oncology_agent` -- Oncology variant interpretation
- `precision_autoimmune_agent` -- Autoimmune disease analysis
- `cart_intelligence_agent` -- CAR-T therapy design
- `imaging_intelligence_agent` -- Medical imaging analysis
- `clinical_trial_intelligence_agent` -- Trial matching and eligibility
- `cardiology_intelligence_agent` -- Cardiovascular genomics
- `neurology_intelligence_agent` -- Neurological variant analysis
- `pharmacogenomics_intelligence_agent` -- Drug-gene interactions
- `rare_disease_diagnostic_agent` -- Rare disease diagnosis
- `single_cell_intelligence_agent` -- Single-cell transcriptomics

**Shared library:** `lib/hcls_common/` (configuration, Milvus client, LLM wrappers, security utilities).

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

1. **Create the directory structure** under `ai_agent_adds/your_agent_name/` following the layout above.

2. **Implement core modules:**
   - `config/settings.py` -- Pydantic `Settings` class with agent-specific configuration
   - `src/agent.py` -- Main agent orchestrator
   - `src/knowledge.py` -- Domain knowledge base and embedding logic
   - `src/rag_engine.py` -- RAG retrieval and response generation

3. **Add the API and UI:**
   - `api/routes.py` -- FastAPI endpoints for the agent
   - `app/streamlit_app.py` -- Streamlit interface for interactive use

4. **Write tests** in `tests/` covering core logic and API routes.

5. **Register the agent** in `docker-compose.dgx-spark.yml` with appropriate port, volume mounts, and health check.

## Code Standards

### Python

- Use **type hints** on all function signatures.
- Use **Pydantic models** for configuration, request/response schemas, and data validation.
- Use **loguru** for logging (not the stdlib `logging` module).
- Keep functions focused; prefer composition over inheritance.

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

This project follows the [Contributor Covenant Code of Conduct](https://www.contributor-covenant.org/version/2/1/code_of_conduct/). By participating, you agree to uphold a welcoming, inclusive, and respectful community.

## License

This project is licensed under the **Apache License 2.0**. See [LICENSE](LICENSE) for details. By contributing, you agree that your contributions will be licensed under the same terms.
