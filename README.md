<div style="display:block; margin: 0px auto; width:300px; text-align: center">

![dance](docs/_static/logo.png)

</div>

[![Travis Build Status](https://travis-ci.com/btjanaka/dance.svg?branch=master)](https://travis-ci.com/btjanaka/dance)
[![codecov](https://codecov.io/gh/btjanaka/dance/branch/master/graph/badge.svg)](https://codecov.io/gh/btjanaka/dance/branch/master)

Visit the [documentation](https://dance.readthedocs.io/) for complete
information on DANCE.

## Table of Contents

<!-- vim-markdown-toc GFM -->

* [Overview](#overview)
* [Design](#design)
* [Installation](#installation)
* [Sample Usage](#sample-usage)
* [Documentation](#documentation)
* [Manifest](#manifest)
* [Contributing](#contributing)
* [Contributors](#contributors)
* [Copyright](#copyright)
* [Acknowledgements](#acknowledgements)

<!-- vim-markdown-toc -->

## Overview

Taking a database such as
[eMolecules](https://www.emolecules.com/info/plus/download-database) and
generating a more manageable dataset is a challenging and important step in
improving molecular dynamics force fields. For instance, in order to improve an
angle parameter, one may need a dataset of molecules with that specific
parameter, so that they can generate QM data for training. DANCE is a pipeline
that allows computational chemists to 1) identify relevant molecules and 2)
select a diverse dataset from among those molecules.

## Design

DANCE is a pipeline that takes in a database of molecules and ultimately outputs
a custom dataset. DANCE is designed with the following goals:

- Extensibility: One should be able to use DANCE to create a wide variety of
  datasets.
- Simplicity: DANCE should require as few parameters as possible to create a
  dataset. Users should be able to focus on the chemical properties of their
  molecules, rather than nitty-gritty pipeline implementation details.

With that in mind, the following diagram shows the layout of DANCE.

![pipeline](docs/_static/pipeline.png)

The pipeline consists of the following steps:

1. **Filter**: DANCE selects relevant molecules from the database. The user
   controls the filtering by providing a **relevance function**. This is a
   function that, when passed a single molecule, decides whether that molecule
   is relevant to the dataset being generated. For instance, one could pass a
   function that returns _True_ only when the molecule has an _a13_ parameter.
   1. We provide tools that enable users to easily create relevance functions.
      To illustrate, one such tool allows one to make a function that only marks
      molecules with a certain parameter as relevant.
2. **Assign Fingerprint**: DANCE assigns a "fingerprint" to each molecule. A
   **fingerprint** contains certain features of the molecule, such as the Wiberg
   Bond Order of selected bonds or aromaticity of certain atoms. The user
   controls the fingerprinting by providing a **fingerprint function**, which
   when given a single molecule, returns a fingerprint for it.
   1. Default fingerprinting options are provided.
   1. Fingerprints need to be tuples of numeric values. This allows us to easily
      order the molecules later on in the pipeline.
3. **Select**: DANCE selects molecules with diverse fingerprints to create the
   final dataset.
   1. Specifically, this step sorts the molecules by their fingerprint and
      selects them at regular intervals (i.e. it selects every _i_-th molecule
      for some constant _i_).

Thus, DANCE is intended to be easy to use. The two main parameters users need to
provide are the **relevance function** (in the Filter step) and the
**fingerprint function** (in the Assign Fingerprint) step.

## Installation

Clone the repository, then run:

```
pip install --extra-index-url https://pypi.anaconda.org/openeye/simple -e .
```

## Sample Usage

TODO

## Documentation

The complete documentation for DANCE is available at
https://dance.readthedocs.io/. To build your own version of the documentation,
see [docs/README.md](docs/README.md).

## Manifest

- **dance** - source code
- **devtools** - tools for people developing DANCE
- **docs** - documentation (built with Sphinx)

## Contributing

Interested in contributing to DANCE? See
[CONTRIBUTING.md](.github/CONTRIBUTING.md)

## Contributors

- [Bryon Tjanaka (UCI)](https://btjanaka.net/)
- [Jessica Maat (UCI)](https://github.com/jmaat)
- [David L. Mobley (UCI)](https://github.com/davidlmobley)

## Copyright

Copyright (c) 2020, Bryon Tjanaka

## Acknowledgements

Project structure based on the
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms)
version 1.1.
