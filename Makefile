.PHONY: default black black-check flake8 test test-v test-vv install serve

default: black-check flake8

black:
	black -l 100 .

black-check:
	black -l 100 --check .

flake8:
	flake8 .

test:
	pytest

test-v:
	pytest -v

test-vv:
	pytest -vv

install:
	pip install -e .

serve:
	clsify web --debug
