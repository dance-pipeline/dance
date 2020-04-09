# Development, testing, and deployment tools

This directory contains a collection of tools for running Continuous Integration
(CI) tests, conda installation, and other development tools not directly related
to the coding process.

## Manifest

### Continuous Integration

For now, we only test on Linux and OSX with Travis CI
[Travis-CI](https://about.travis-ci.com/). This directory contains files related
to such testing:

- `travis-ci`: directory containing files for Travis CI
  - `before_install.sh`: Pip/Miniconda pre-package installation script for
    Travis
  - `oe_license.txt.enc`: Encrypted Openeye license for running tests on Travis
    CI

### Conda Environment

This directory contains files to set up the Conda environment for testing
purposes.

- `conda-envs`: directory containing the YAML file(s) which fully describe Conda
  Environments, their dependencies, and those dependency provenance's
  - `test_env.yaml`: Simple test environment file with base dependencies.
    Channels are not specified here and therefore respect global Conda
    configuration

### Additional Scripts

This directory contains OS agnostic helper scripts which don't fall in any of
the previous categories

- `scripts`
  - `create_conda_env.py`: Helper program for spinning up new conda environments
    based on a starter file with Python Version and Env. Name command-line
    options
  - `run_tests.sh`: A script for running tests for `dance`. Intended to be run
    from the root of the repository.
