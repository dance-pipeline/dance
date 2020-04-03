# Contributing

We welcome contributions from external contributors. This document describes how
to merge code changes into DANCE.

## Table of Contents

<!-- vim-markdown-toc GFM -->

* [How to contribute changes](#how-to-contribute-changes)
* [Checklist for updates](#checklist-for-updates)
* [Documentation](#documentation)
* [Running Tests](#running-tests)
* [Versioneer Auto-version](#versioneer-auto-version)
* [Auto-Formatting](#auto-formatting)
* [Codecov Code Coverage](#codecov-code-coverage)
* [LGTM Code Analysis](#lgtm-code-analysis)
* [Devtools](#devtools)

<!-- vim-markdown-toc -->

## How to contribute changes

- Clone the repository if you have write access to the main repo, fork the
  repository if you are a collaborator.
- Make a new branch with `git checkout -b {your branch name}`
- Make changes and test your code
- Ensure that the test environment dependencies (`conda-envs`) line up with the
  build and deploy dependencies (`conda-recipe/meta.yaml`)
- Push the branch to the repo (either the main or your fork) with
  `git push -u origin {your branch name}`
  - Note that `origin` is the default name assigned to the remote, yours may be
    different
- Make a PR on GitHub with your changes
- We'll review the changes and get your code into the repo after lively
  discussion!

## Checklist for updates

- [ ] Make sure there is an/are issue(s) opened for your specific update
- [ ] Create the PR, referencing the issue
- [ ] Debug the PR as needed until tests pass
- [ ] Tag the final, debugged version
  - `git tag -a X.Y.Z [latest pushed commit] && git push --follow-tags`
- [ ] Get the PR merged in

## Documentation

For information on building the documentation (including setting up a server to
automatically reload the docs), see [here](../docs/README.md).

## Running Tests

To run the tests, run `devtools/scripts/run_tests.sh` from the root of this
repository.

## Versioneer Auto-version

[Versioneer](https://github.com/warner/python-versioneer) will automatically
infer what version is installed by looking at the `git` tags and how many
commits ahead this version is. The format follows
[PEP 440](https://www.python.org/dev/peps/pep-0440/) and has the regular
expression of:

```regexp
\d+.\d+.\d+(?\+\d+-[a-z0-9]+)
```

If the version of this commit is the same as a `git` tag, the installed version
is the same as the tag, e.g. `dance-0.1.2`, otherwise it will be appended with
`+X` where `X` is the number of commits ahead from the last tag, and then
`-YYYYYY` where the `Y`'s are replaced with the `git` commit hash.

## Auto-Formatting

Please auto-format your code with [YAPF](https://github.com/google/yapf).
To install YAPF, run `pip install yapf`, and to run it on a file, use
`yapf -i <FILE>.py`. YAPF will automatically detect the settings contained in
`setup.cfg` and format your code appropriately. If you are using an editor such
as vim, you may also want to look into using a YAPF plugin to automatically
format your code as you write it.

## Codecov Code Coverage

We use [Codecov](https://codecov.io) to automatically check test coverage on
PRs. Please make sure your code is reasonably well-covered by your tests.

## LGTM Code Analysis

We use [LGTM](https://lgtm.com) to perform static analysis on PRs and detect
various programming mistakes. Please address feedback, if any, from LGTM.

## Devtools

You may wish to take a look at the [devtools](../devtools/README.md).
