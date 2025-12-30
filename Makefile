notebook-docs:
	uv run jupyter nbconvert docs/notebooks/*.ipynb \
		--config jupyter_nbconvert_config.py \
		--to markdown

serve:
	uv run mkdocs serve
