# Compiling DANCE's Documentation

The docs for this project are built with
[Sphinx](http://www.sphinx-doc.org/en/master/).

## Basic Usage

Install the dependencies with:

```bash
pip install -r docs/requirements.txt
```

Run the following command to automatically serve the documentation at
http://localhost:5500/:

```bash
sphinx-reload docs/ --watch {dance,docs}/
```

This command will also open a tab in your web browser for displaying the
documentation. As you make changes to your documentation, as well as to any code
in `dance`, the documentation will automatically update in your browser.

## Advanced Usage

After installing the dependencies, use the `Makefile` in this directory to
compile static HTML pages by:

```bash
make html
```

The compiled docs will be in the `_build` directory and can be viewed by opening
`index.html` (which may itself be inside a directory called `html/` depending on
what version of Sphinx is installed).

For more options for building the documentation, run

```bash
make
```

and visit https://www.sphinx-doc.org/en/master/man/sphinx-build.html for more
info.
