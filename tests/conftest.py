#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
Pytest configuration shared by the whole test-suite.

  * puts the repository root and tests/ on sys.path, so the test files can be
    collected from any working directory and `import _tools` always works
  * registers the `slow` marker and skips slow tests unless --runslow is given
  * provides the `silence` fixture (stdout suppression of the Fortran wrappers)

The test files stay runnable without pytest (python tests/test_xxx.py); this
file is only picked up by pytest itself.
"""
import os
import sys

_TESTS = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_TESTS)
for _p in (_ROOT, _TESTS):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import pytest


def pytest_addoption(parser):
    parser.addoption(
        '--runslow', action='store_true', default=False,
        help='also run tests marked with @pytest.mark.slow',
    )


def pytest_collection_modifyitems(config, items):
    if config.getoption('--runslow'):
        return
    skip_slow = pytest.mark.skip(reason='slow test: use --runslow to run')
    for item in items:
        if 'slow' in item.keywords:
            item.add_marker(skip_slow)


@pytest.fixture
def silence():
    """Swallow stdout of the test body, including Fortran prints (fd-level).
    Usage: def test_x(silence): ... (activated for the whole test) -- or use
    _tools.silence() as a context manager for finer scope."""
    import _tools
    with _tools.silence():
        yield
