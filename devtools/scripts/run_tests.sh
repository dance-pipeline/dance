#!/bin/bash
# Runs tests. Intended to be run from the root of this repository.
pytest -v --cov-report term-missing --cov=dance dance/tests/
