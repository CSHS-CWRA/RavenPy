.PHONY: clean clean-build clean-pyc clean-test coverage dist docs help install lint lint/flake8 lint/black
.DEFAULT_GOAL := help

define BROWSER_PYSCRIPT
import os, webbrowser, sys

from urllib.request import pathname2url

webbrowser.open("file://" + pathname2url(os.path.abspath(sys.argv[1])))
endef
export BROWSER_PYSCRIPT

define PRINT_HELP_PYSCRIPT
import re, sys

for line in sys.stdin:
	match = re.match(r'^([a-zA-Z_-]+):.*?## (.*)$$', line)
	if match:
		target, help = match.groups()
		print("%-20s %s" % (target, help))
endef
export PRINT_HELP_PYSCRIPT

BROWSER := python -c "$$BROWSER_PYSCRIPT"

help:
	@python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)

clean: clean-build clean-pyc clean-test ## remove all build, test, coverage and Python artifacts

clean-build: ## remove build artifacts
	rm -fr build/
	rm -fr dist/
	rm -fr .eggs/
	find . -name '*.egg-info' -exec rm -fr {} +
	find . -name '*.egg' -exec rm -f {} +

clean-docs: ## remove docs artifacts
	rm -f docs/apidoc/ravenpy*.rst
	rm -f docs/apidoc/modules.rst
	$(MAKE) -C docs clean

clean-pyc: ## remove Python file artifacts
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -fr {} +

clean-test: ## remove test and coverage artifacts
	rm -fr .tox/
	rm -f .coverage
	rm -fr htmlcov/
	rm -fr .pytest_cache

lint/black: ## check style with black
	flake8 ravenpy
	black --check ravenpy tests

lint: lint/black ## check style

test: ## run tests quickly with the default Python
	pytest

test-all: ## run tests on every Python version with tox
	tox

coverage: ## check code coverage quickly with the default Python
	coverage run --source ravenpy -m pytest
	coverage report -m
	coverage html
	$(BROWSER) htmlcov/index.html

autodoc: clean-docs ## create sphinx-apidoc files:
	sphinx-apidoc -o docs/apidoc --private --module-first ravenpy

autodoc-custom-index: clean-docs ## create sphinx-apidoc files but with special index handling for indices and indicators
	env SPHINX_APIDOC_OPTIONS="members,undoc-members,show-inheritance,noindex" sphinx-apidoc -o docs/apidoc --private --module-first ravenpy

linkcheck: autodoc ## run checks over all external links found throughout the documentation
	$(MAKE) -C docs linkcheck

docs: autodoc-custom-index ## generate Sphinx HTML documentation, including API docs
	$(MAKE) -C docs html
ifndef READTHEDOCS
	$(BROWSER) docs/_build/html/index.html
endif

servedocs: docs ## compile the docs watching for changes
	watchmedo shell-command -p '*.rst' -c '$(MAKE) -C docs html' -R -D .

release: dist ## package and upload a release
	twine upload dist/*

dist: clean ## builds source and wheel package
	flit build
	ls -l dist

install: clean ## install the package to the active Python's site-packages
	python -m pip install --no-user .

develop: clean ## install the package and development dependencies in editable mode to the active Python's site-packages
	python -m pip install --no-user --editable ".[dev]"
