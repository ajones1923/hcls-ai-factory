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

If a system cannot be understood, it cannot be trusted—especially in healthcare.

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

Speed matters—but not at the expense of trust.

### Open Does Not Mean Unstructured

Open contribution does not imply lack of rigor. Contributions are expected to:
- Follow project conventions
- Include tests or validation where appropriate
- Be well-scoped and reviewable

The aim is collective progress, not unchecked accumulation.

### Shared Credit, Shared Responsibility

Contributors retain ownership of their work and receive attribution, while agreeing to uphold the standards and intent of the project. Apache 2.0 ensures freedom to build, but responsibility to be transparent about changes.

This project grows through collaboration—but it endures through care.

---

## How to Contribute

### Reporting Issues

- Use GitHub Issues to report bugs or request features
- Include steps to reproduce for bugs
- Describe expected vs. actual behavior
- Include system information (GPU, OS, Docker version)

### Development Setup

Before developing or testing changes, download the required data:

```bash
cp .env.example .env
# Edit .env with your NGC_API_KEY and ANTHROPIC_API_KEY

# Download all data (~500 GB, one-time)
./setup-data.sh --all

# Or just Stage 2 for RAG/Chat development (~2 GB, fast)
./setup-data.sh --stage2
```

See [docs/DATA_SETUP.md](docs/DATA_SETUP.md) for troubleshooting and stage-by-stage options.

### Submitting Changes

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/your-feature`)
3. Make your changes
4. Test your changes
5. Commit with clear messages (`git commit -m "Add feature X"`)
6. Push to your fork (`git push origin feature/your-feature`)
7. Open a Pull Request

### Pull Request Guidelines

- Keep changes focused and atomic
- Update documentation if needed
- Add tests for new functionality
- Ensure existing tests pass
- Follow existing code style

### Areas for Contribution

- **Documentation**: Improve guides, add examples, fix typos
- **Pipelines**: Enhance genomics, RAG, or drug discovery pipelines
- **Integrations**: Add support for new tools or databases
- **UI/UX**: Improve Streamlit interfaces
- **Testing**: Add test coverage
- **Performance**: Optimize for different hardware configurations

---

## Code of Conduct

Please review our [Code of Conduct](CODE_OF_CONDUCT.md) before contributing. This project operates at the intersection of healthcare, life sciences, and AI, and requires a higher standard of professionalism, respect, and responsibility.

---

## Questions

Open an issue with the "question" label for general questions about the project.

---

## License

By contributing, you agree that your contributions will be licensed under the Apache License 2.0.
