import pytest
from bdgtools import BedGraph, Regions

@pytest.fixture
def bedgraph():
    return BedGraph([0, 10, 15, 25, 40], [0, 1, 2, 3, 4], size=50)

@pytest.fixture
def regions_10b():
    return Regions([2, 13, 17], [12, 23, 27], [1, -1, 1])
