notebook-docs:
	IS_GENERATING_DOCS=true uv run jupyter nbconvert docs/notebooks/*.ipynb \
		--config jupyter_nbconvert_config.py \
		--NbConvertApp.output_files_dir notebook-assets \
		--to markdown

notebook-docs-debug:
	IS_GENERATING_DOCS=true uv run jupyter nbconvert docs/notebooks/*.ipynb \
		--config jupyter_nbconvert_config.py \
		--NbConvertApp.output_files_dir notebook-assets \
		--log-level=DEBUG \
		--to markdown

serve: ## Serve documentation website for local development
	uv run mkdocs serve

docs: ## Build documentation website directory
	uv run mkdocs build

build: ## Build wheel package for publishing to PyPI
	rm -rf dist/
	uv build

publish: ## Publish built package to PyPI
	uv publish

# These make tasks allow the default help text to work properly.
%:
	@true

.PHONY: help notebook-docs notebook-docs-debug serve docs build publish

help:
	@echo 'Usage: make <command>'
	@echo
	@echo 'where <command> is one of the following:'
	@echo
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

.DEFAULT_GOAL := help