# Contributing to HCLS AI Factory

Thank you for your interest in contributing to the HCLS AI Factory project.

---

## Contributor Philosophy

Contributions to this project are guided by the same principles that motivated it to be open: clarity, responsibility, and respect for the problem space.

### Build for Understanding First

Contributions should prioritize:
- Readability over cleverness
- Explicit workflows over implicit behavior
- Documentation alongside code

If a system cannot be understood, it cannot be trusted — especially in healthcare.

### Favor Composability Over Complexity

This project values modular, extensible designs that allow others to adapt workflows without rewrites. Contributions should avoid introducing unnecessary coupling, hidden dependencies, or tightly bound assumptions.

The goal is to make the system easier to extend, not harder to reason about.

### Preserve Transparency and Traceability

All changes should maintain or improve:
- Reproducibility
- Observability
- Clear data lineage
- Explainable reasoning paths

Black-box behavior is actively discouraged.

### Respect the Human Context

This project exists because real people are affected by the systems we build. Contributions should be made with an awareness of:
- Clinical implications
- Ethical considerations
- The difference between research tooling and clinical deployment

Speed matters — but not at the expense of trust.

### Shared Credit, Shared Responsibility

Contributors retain ownership of their work and receive attribution, while agreeing to uphold the standards and intent of the project. Apache 2.0 ensures freedom to build, but responsibility to be transparent about changes.

This project grows through collaboration — but it endures through care.

---

## Project Structure

The repository is organized as four independent services plus an orchestration layer:

```
hcls-ai-factory/
├── genomics-pipeline/       # Stage 1: FASTQ → VCF (Parabricks)
│   └── web-portal/          #   Flask portal for pipeline control
├── rag-chat-pipeline/       # Stage 2: VCF → Target (Milvus + LLM)
│   ├── src/                 #   Core library (parsing, embedding, search)
│   └── portal/              #   Streamlit chat interface
├── drug-discovery-pipeline/ # Stage 3: Target → Molecules (BioNeMo)
│   ├── src/                 #   Core library (scoring, docking, generation)
│   └── app/                 #   Streamlit discovery interface
├── landing-page/            # Service health dashboard
├── hls-orchestrator/        # Nextflow DSL2 orchestrator
└── docs/                    # MkDocs documentation site
```

Each pipeline is independently deployable. They communicate through API endpoints and standard file formats (VCF, JSON, SDF) — not shared Python imports.

---

## Development Setup

### Prerequisites

- Python 3.10+
- Docker (for containerized services)
- Git

GPU access is required only for running the full pipeline end-to-end. All unit tests run on CPU.

### Per-Service Setup

Each service has its own `requirements.txt`. Install only what you need:

```bash
# Example: working on the RAG pipeline
cd rag-chat-pipeline
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
pip install pytest pytest-cov
```

### Data Setup

Before running pipelines (not needed for unit tests), download the required data:

```bash
cp .env.example .env
# Edit .env with your NGC_API_KEY and ANTHROPIC_API_KEY

# Download all data (~500 GB, one-time)
./setup-data.sh --all

# Or just Stage 2 for RAG/Chat development (~2 GB, fast)
./setup-data.sh --stage2
```

See [docs/DATA_SETUP.md](docs/DATA_SETUP.md) for troubleshooting and stage-by-stage options.

### Environment Variables

Services use environment variables for configuration. Each pipeline directory contains example configuration in its README or a `.env.example` file. Required variables vary by service. Never hardcode paths, hostnames, or credentials.

---

## Running Tests

CI runs tests for all four services on every push and pull request. Run them locally before submitting:

```bash
# RAG Chat Pipeline
cd rag-chat-pipeline
python -m pytest tests/ -v -m "not slow and not integration"

# Drug Discovery Pipeline
cd drug-discovery-pipeline
python -m pytest tests/ -v -m "not slow and not integration and not rdkit"

# Genomics Portal
cd genomics-pipeline/web-portal
python -m pytest tests/ -v -m "not slow and not integration"

# Landing Page
cd landing-page
python -m pytest tests/ -v -m "not slow and not integration"
```

### Test Markers

Use pytest markers to categorize tests appropriately:

| Marker | When to Use |
|--------|-------------|
| `@pytest.mark.slow` | Tests that take more than a few seconds (large datasets, heavy computation) |
| `@pytest.mark.integration` | Tests that require running services (Milvus, Docker, API endpoints) |
| `@pytest.mark.rdkit` | Tests that require RDKit or other GPU/chemistry libraries |

CI skips `slow`, `integration`, and `rdkit` tests. These run in local or staging environments with the full stack available.

### Test Requirements

- All new code must include tests. All tests must pass before merge.
- Security-sensitive code (input validation, query construction, authentication) requires explicit test coverage.
- Use mocks for external services — tests must run without Docker, GPUs, or network access.
- Write clear, focused tests. One assertion per behavior, descriptive test names.

---

## Making Changes

### Branch Naming

Use descriptive branch names with a category prefix:

- `feature/add-variant-filtering` — new functionality
- `fix/milvus-query-timeout` — bug fixes
- `docs/update-deployment-guide` — documentation changes
- `test/add-scoring-coverage` — test additions

### Code Standards

- **Linting**: [Ruff](https://docs.astral.sh/ruff/) is enforced in CI. Run `ruff check .` locally before pushing.
- **Type hints**: Expected on public function signatures.
- **Imports**: Follow the convention already in the codebase — stdlib, then third-party, then local.
- **No hardcoded paths**: Use environment variables or path resolution relative to `PROJECT_ROOT`.
- **Dependencies**: Pin all new dependencies to exact versions in `requirements.txt`. Use `package==x.y.z`, not `>=`.

### Commit Messages

Write clear commit messages that explain *why*, not just *what*:

```
Fix Milvus filter expression injection in gene search

Add input sanitization for gene names and chromosome identifiers
before they reach f-string filter expressions in search_by_gene()
and search_by_region().
```

Keep the subject line under 72 characters. Use the body for context when the change is not self-evident.

---

## Pull Request Process

1. **One pipeline per PR** when possible. Cross-pipeline changes (like orchestrator updates) are the exception.
2. **Fill out the PR description** — see template below.
3. **All CI checks must pass** — lint, tests, and docs build.
4. **Include screenshots** for UI changes (Streamlit, landing page, portal).
5. **Keep changes focused and atomic** — update documentation if behavior changes.
6. A maintainer will review your PR. Address feedback, then the maintainer will merge.

### PR Description Template

```markdown
## Summary
Brief description of the change and why it's needed.

## Pipeline(s) Affected
- [ ] Genomics Pipeline
- [ ] RAG Chat Pipeline
- [ ] Drug Discovery Pipeline
- [ ] Landing Page
- [ ] Orchestrator
- [ ] Documentation

## Test Plan
How was this tested? What commands to run?

## Notes
Any context the reviewer should know.
```

### Areas for Contribution

- **Documentation**: Improve guides, add examples, fix typos
- **Pipelines**: Enhance genomics, RAG, or drug discovery workflows
- **Knowledge base**: Add gene-drug-disease connections, PDB structures, pathway data
- **Integrations**: Add support for new tools, models, or databases
- **UI/UX**: Improve Streamlit interfaces and the landing page
- **Testing**: Expand test coverage, especially integration tests
- **Performance**: Optimize for different hardware configurations

---

## Architecture Principles

These principles keep the platform maintainable as it grows:

1. **Pipeline independence** — Each service deploys, scales, and fails independently. Do not create cross-pipeline Python imports. Shared contracts are API endpoints and file formats.

2. **Orchestration belongs in Nextflow** — If you need to coordinate between pipelines, that logic goes in the Nextflow orchestrator, not in application code.

3. **Validate at boundaries** — All user input, API parameters, and file paths must be validated before use. Internal function calls between trusted modules do not need redundant validation.

4. **Keep it deployable on a single node** — The platform is designed to run on a single GPU workstation. Contributions should not introduce hard dependencies on distributed infrastructure, cloud services, or multi-node clusters.

5. **Research context** — This platform is for research and education, not clinical use. Do not add features that imply clinical decision support without explicit discussion with maintainers.

---

## Security

Security is critical in healthcare and life sciences software, even in research contexts.

### Requirements

- **Never commit credentials**, API keys, tokens, or patient data. Use environment variables.
- **Validate all inputs** that reach database queries, shell commands, or file system operations. See `rag-chat-pipeline/src/milvus_client.py` for the sanitization pattern used in this project.
- **No real patient data** in tests, examples, or documentation. Use synthetic or publicly available reference data only.
- **Report vulnerabilities privately** — use [GitHub Security Advisories](https://docs.github.com/en/code-security/security-advisories) or contact the maintainers directly. Do not open public issues for security vulnerabilities.

### Sanitization Pattern

When constructing filter expressions or queries from user input, always validate first:

```python
import re

_GENE_PATTERN = re.compile(r'^[A-Za-z0-9_\-]+$')

def _sanitize_gene(gene: str) -> str:
    if not gene or not _GENE_PATTERN.match(gene):
        raise ValueError(f"Invalid gene name: {gene!r}")
    return gene

# Use before any query construction
gene = _sanitize_gene(user_input)
results = collection.query(expr=f'gene == "{gene}"', ...)
```

---

## For Organizations Forking This Project

If your team is forking the HCLS AI Factory to adapt it for your own infrastructure, datasets, or research focus, these recommendations will help you maintain a productive fork.

### Initial Setup

1. **Fork, don't clone** — maintain the upstream relationship so you can pull improvements.
2. **Run the full test suite** immediately after forking to establish a green baseline.
3. **Set up your own CI** using the existing `.github/workflows/ci.yml` as a starting point. Adapt the matrix if you add or remove services.
4. **Designate a fork maintainer** — one person should own the upstream relationship and coordinate sync decisions.

### Customization Points

The platform is designed to be adapted at these layers without modifying core library code:

| Layer | Location | What to Customize |
|-------|----------|-------------------|
| **Configuration** | Environment variables | Hostnames, ports, API keys, model endpoints |
| **Knowledge base** | `rag-chat-pipeline/src/knowledge.py` | Gene-drug-disease connections for your research focus |
| **Scoring weights** | `drug-discovery-pipeline/src/scoring.py` | Generation, docking, and QED weights for your targets |
| **Orchestration** | `hls-orchestrator/nextflow.config` | Resource allocation, container registries, executor settings |
| **Documentation** | `docs/` | Internal guides, architecture updates, branding |
| **Data sources** | `setup-data.sh`, `.env` | Reference genomes, annotation databases, PDB structures |

### Staying in Sync with Upstream

- **Rebase regularly** against upstream `main` to pick up security fixes, new tests, and knowledge base updates.
- **Keep customizations in configuration** — avoid modifying core library files when environment variables or extension points exist.
- **Contribute back** — if you fix a bug or add a broadly useful feature, open a PR upstream. The community benefits and you reduce future merge conflicts.
- **Document deviations** — if you change architecture decisions, note them in your fork's README so your team understands what diverged and why.

### Internal Development Standards

Apply the same testing, review, and security standards described in this guide. The conventions exist to keep the platform trustworthy — they matter more, not less, when the platform is adapted for a specific organization's workflows.

---

## Code of Conduct

Please review our [Code of Conduct](CODE_OF_CONDUCT.md) before contributing. This project operates at the intersection of healthcare, life sciences, and AI, and requires a higher standard of professionalism, respect, and responsibility.

---

## Questions?

- Open a [GitHub Discussion](https://github.com/NVIDIA/hcls-ai-factory/discussions) for questions or ideas.
- Open an [Issue](https://github.com/NVIDIA/hcls-ai-factory/issues) for bugs or feature requests.

---

## License

By contributing, you agree that your contributions will be licensed under the Apache License 2.0.
