import pytest
from chiptools.regions import Regions, Region

@pytest.fixture
def regions():
    return Regions([0, 10, 13], [3, 12, 17], [1, -1, 1])

def test_iter(regions):
    assert list(regions) == [(0,3,1), (10,12,-1), (13, 17, 1)]
