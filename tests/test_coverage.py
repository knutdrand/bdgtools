import pytest

from bdgtools.regions import Regions
from bdgtools.bedgraph import BedGraph
from bdgtools.coverage import get_coverage

@pytest.fixture
def regions():
    return Regions([2, 2, 3, 5, 5, 7],
                   [5, 6, 7, 8, 9, 10])

@pytest.fixture
def bedgraph():
    return BedGraph([0, 2, 3, 5, 6, 8, 9, 10],
                    [0, 2, 3, 4, 3, 2, 1, 0])

def test_get_coverage(regions, bedgraph):
    print(get_coverage(regions))
    assert get_coverage(regions) == bedgraph
