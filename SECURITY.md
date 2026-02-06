# Security Policy

## Supported Versions

| Version | Supported          |
| ------- | ------------------ |
| latest  | :white_check_mark: |

## Reporting a Vulnerability

If you discover a security vulnerability in the HCLS AI Factory, please report it responsibly.

**Do not open a public GitHub issue for security vulnerabilities.**

Instead, please email **security@hcls-ai-factory.org** with:

- A description of the vulnerability
- Steps to reproduce the issue
- Potential impact assessment
- Any suggested fixes (optional)

You can expect an initial response within 72 hours. We will work with you to understand and address the issue before any public disclosure.

## Scope

This policy covers the HCLS AI Factory codebase including:

- Genomics pipeline (Parabricks integration)
- RAG/Chat pipeline (Milvus, LLM integration)
- Drug discovery pipeline (BioNeMo NIMs)
- Web portals and APIs
- Deployment configurations

## Security Considerations

This platform processes genomic data and interacts with AI services. When deploying:

- Never expose internal service ports (Milvus, Grafana, Prometheus) to the public internet without authentication
- Use environment variables for all API keys and credentials
- Review CORS settings before production deployment
- Ensure genomic data is handled in compliance with applicable regulations (HIPAA, GDPR)
- Keep all container images and dependencies up to date
