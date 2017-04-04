"""Custom ``pytest`` configurations."""
import pytest


def pytest_addoption(parser):
    parser.addoption("--remote-data", action="store_true",
                     help="Use remote data.")


def pytest_configure(config):
    config.getini('markers').append(
        'remote_data: Run tests that require data from remote servers')


def pytest_runtest_setup(item):
    if ('remote_data' in item.keywords and
            not item.config.getvalue("remote_data")):
        pytest.skip("need --remote-data option to run")
