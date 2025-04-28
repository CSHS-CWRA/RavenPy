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
LOCALES := docs/locales

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
	rm -fr docs/locales/fr/LC_MESSAGES/*.mo
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

lint/flake8: ## check style with flake8
	python -m ruff check src/ravenpy tests
	python -m flake8 --config=.flake8 src/ravenpy tests
	python -m numpydoc lint src/ravenpy/**.py

lint/black: ## check style with black
	python -m black --check src/ravenpy tests
	python -m blackdoc --check src/ravenpy docs
	python -m isort --check src/ravenpy tests

lint: lint/black ## check style

test: ## run tests quickly with the default Python
	python -m pytest

test-all: ## run tests on every Python version with tox
	python -m tox

coverage: ## check code coverage quickly with the default Python
	python -m coverage run --source src/ravenpy -m pytest
	python -m coverage report -m
	python -m coverage html
	$(BROWSER) htmlcov/index.html

														autodoc: clean-docs ## create sphinx-apidoc files:
	sphinx-apidoc -o docs/apidoc --private --module-first src/ravenpy

autodoc-custom-index: clean-docs ## create sphinx-apidoc files but with special index handling for indices and indicators
	env SPHINX_APIDOC_OPTIONS="members,undoc-members,show-inheritance,noindex" sphinx-apidoc -o docs/apidoc --private --module-first src/ravenpy

linkcheck: autodoc ## run checks over all external links found throughout the documentation
	$(MAKE) -C docs linkcheck

initialize-translations: autodoc-custom-index ## initialize translations, ignoring autodoc-generated files
	${MAKE} -C docs gettext
	sphinx-intl update -p docs/_build/gettext -d docs/locales -l fr
	rm -fr docs/locales/fr/LC_MESSAGES/apidoc

docs: autodoc-custom-index ## generate Sphinx HTML documentation, including API docs
	$(MAKE) -C docs html BUILDDIR="_build/html/en"
ifneq ("$(wildcard $(LOCALES))","")
	${MAKE} -C docs gettext
	$(MAKE) -C docs html BUILDDIR="_build/html/fr" SPHINXOPTS="-D language='fr'"
endif
ifndef READTHEDOCS
	$(BROWSER) docs/_build/html/en/html/index.html
endif

servedocs: docs ## compile the docs watching for changes
	watchmedo shell-command -p '*.rst' -c '$(MAKE) -C docs html' -R -D .

dist: clean ## builds source and wheel package
	python -m flit build
	ls -l dist

release: dist ## package and upload a release
	python -m flit publish dist/*

install: clean ## install the package to the active Python's site-packages
	python -m pip install --no-user .

develop: clean ## install the package and development dependencies in editable mode to the active Python's site-packages
	python -m pip install --no-user --editable ".[all]"
	pre-commit install
