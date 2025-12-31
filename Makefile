notebook-docs:
	IS_GENERATING_DOCS=true uv run jupyter nbconvert docs/notebooks/*.ipynb \
		--config jupyter_nbconvert_config.py \
		--NbConvertApp.output_files_dir notebook-assets \
		--to markdown

serve:
	uv run mkdocs serve

build:
	uv run mkdocs build
