"""Shared pytest fixtures."""

import pytest
from nurecoil.nucleus import Nucleus


@pytest.fixture
def ge76() -> Nucleus:
    return Nucleus(Z=32, A=76)


@pytest.fixture
def si28() -> Nucleus:
    return Nucleus(Z=14, A=28)
