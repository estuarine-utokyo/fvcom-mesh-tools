.PHONY: install test lint format clean wheel

install:
	pip install -e ".[dev]"

test:
	pytest -q

lint:
	ruff check .

format:
	ruff format .

clean:
	rm -rf build dist *.egg-info .pytest_cache .ruff_cache .mypy_cache htmlcov
	find . -type d -name __pycache__ -prune -exec rm -rf {} +

# Build a clean wheel into dist/. Depends on clean so a stale build/
# tree from a prior in-place install cannot leak removed modules
# (e.g. io/dem_subset.py, mesh_engine/depth.py) into the wheel.
wheel: clean
	pip wheel . --no-deps -w dist/
