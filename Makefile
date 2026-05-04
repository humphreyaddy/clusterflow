.PHONY: install test lint docs docker-build clean

install:
	pip install -e ".[dev]"

test:
	pytest

lint:
	ruff check clusterflow tests
	black --check clusterflow tests
	mypy clusterflow

docs:
	mkdocs build

docker-build:
	docker build -t clusterflow:latest .

clean:
	rm -rf build dist *.egg-info .pytest_cache .mypy_cache .ruff_cache .coverage
	find . -type d -name __pycache__ -exec rm -rf {} +
