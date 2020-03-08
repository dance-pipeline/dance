"""
Unit and regression test for the dance package.
"""

# Import package, test suite, and other packages as needed
import dance
import pytest
import sys

def test_dance_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "dance" in sys.modules
