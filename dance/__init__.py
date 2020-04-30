"""
DANCE
Pipeline for generating molecule datasets.
"""

# Imports
import dance.run
from dance.dance_pipeline import DancePipeline

# Handle Versioneer
from ._version import get_versions

versions = get_versions()  # pylint: disable=invalid-name
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
