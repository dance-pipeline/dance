"""Tests for DancePipeline.
"""
from dance import DancePipeline
import pytest


def test_init_raises_exception_with_bad_database_type():
    with pytest.raises(RuntimeError):
        dp = DancePipeline("FOOBAR", "foobar.xyz")
