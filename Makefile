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

serve:
	uv run mkdocs serve

build:
	uv run mkdocs build
