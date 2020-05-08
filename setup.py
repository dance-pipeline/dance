"""
DANCE
Pipeline for generating datasets.
"""
import sys

from setuptools import find_packages, setup

import versioneer

short_description = __doc__.split("\n")

# from https://github.com/pytest-dev/pytest-runner#conditional-requirement
needs_pytest = {'pytest', 'test', 'ptr'}.intersection(sys.argv)
pytest_runner = ['pytest-runner'] if needs_pytest else []

try:
    with open("README.md", "r") as handle:
        long_description = handle.read()
except:
    long_description = "\n".join(short_description[2:])

setup(
    # Self-descriptive entries which should always be present
    name='dance',
    author='Bryon Tjanaka',
    author_email='bryon@btjanaka.net',
    description=short_description[0],
    long_description=long_description,
    long_description_content_type="text/markdown",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license='MIT',

    # Which Python importable modules should be included with DANCE upon
    # installation.
    packages=find_packages(),

    # Optional include package data to ship with DANCE -- Customize MANIFEST.in
    # if necessary.
    include_package_data=True,

    # Allows `setup.py test` to work correctly with pytest
    setup_requires=[] + pytest_runner,

    # Required packages for DANCE
    install_requires=["openeye-toolkits", "sortedcontainers"],

    # Additional info
    url="https://github.com/btjanaka/dance/",
    project_urls={
        "Bug Tracker": "https://github.com/btjanaka/dance/issues/",
        "Documentation": "https://dance.readthedocs.io/",
        "Source Code": "https://github.com/btjanaka/dance/",
    },
)
