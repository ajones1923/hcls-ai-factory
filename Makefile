.PHONY: help lint format typecheck test test-rag test-drug test-genomics test-landing security docs clean

help: ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

lint: ## Run ruff linter
	ruff check .

format: ## Format code with ruff
	ruff format .
	ruff check --fix .

typecheck: ## Run mypy type checking
	mypy --ignore-missing-imports landing-page/ || true

test: ## Run all service tests
	@echo "=== RAG Chat Pipeline ==="
	cd rag-chat-pipeline && python -m pytest tests/ -m "not slow and not integration" --tb=short -q || true
	@echo "\n=== Drug Discovery Pipeline ==="
	cd drug-discovery-pipeline && python -m pytest tests/ -m "not slow and not integration and not rdkit" --tb=short -q || true
	@echo "\n=== Genomics Portal ==="
	cd genomics-pipeline/web-portal && python -m pytest tests/ -m "not slow and not integration" --tb=short -q || true
	@echo "\n=== Landing Page ==="
	cd landing-page && python -m pytest tests/ -m "not slow and not integration" --tb=short -q || true

test-rag: ## Run RAG chat pipeline tests
	cd rag-chat-pipeline && python -m pytest tests/ -v -m "not slow and not integration"

test-drug: ## Run drug discovery pipeline tests
	cd drug-discovery-pipeline && python -m pytest tests/ -v -m "not slow and not integration and not rdkit"

test-genomics: ## Run genomics portal tests
	cd genomics-pipeline/web-portal && python -m pytest tests/ -v -m "not slow and not integration"

test-landing: ## Run landing page tests
	cd landing-page && python -m pytest tests/ -v -m "not slow and not integration"

security: ## Run pip-audit security scan
	@for dir in rag-chat-pipeline drug-discovery-pipeline genomics-pipeline/web-portal landing-page; do \
		echo "=== $$dir ==="; \
		cd $$dir && pip-audit --format columns 2>/dev/null || true; \
		cd $(CURDIR); \
	done

docs: ## Build MkDocs documentation
	mkdocs build

docs-serve: ## Serve docs locally with live reload
	mkdocs serve

clean: ## Remove build artifacts and caches
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name .pytest_cache -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name .mypy_cache -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name .ruff_cache -exec rm -rf {} + 2>/dev/null || true
	find . -name "*.pyc" -delete 2>/dev/null || true
	find . -name ".coverage" -delete 2>/dev/null || true
	find . -name "coverage.xml" -delete 2>/dev/null || true
	rm -rf site/ 2>/dev/null || true

preflight: ## Run preflight checks
	./scripts/preflight-check.sh
