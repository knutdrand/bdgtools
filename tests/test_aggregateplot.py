import pytest
import numpy as np
import pandas as pd

from bdgtools.aggregateplot import *
from bdgtools import BedGraph, Regions
from .fixtures import bedgraph, regions_10b

@pytest.fixture
def regions():
    return Regions([100, 200, 300, 400, 500], [101, 202, 303, 404, 505], [1, 1, 1, 1, 1])

@pytest.fixture
def graph():
    indices = [0, 100, 101, 200, 201, 300, 301, 400, 401, 500, 501]
    return BedGraph(indices,
                    [0, 1]*(len(indices)//2) + [0], size=600)
@pytest.fixture
def true_signal():
    return np.sum([[0,0,0,0,0,0,0,0,1,1],
                   [2,2,2,2,2,2,2,2,1,1],
                   [2,2,2,2,2,2,2,2,3,3]], axis=0)/3
@pytest.fixture
def true_matrix():
    return np.array([[0,0,0,0,0,0,0,0,1,1],
                     [2,2,2,2,2,2,2,2,1,1],
                     [2,2,2,2,2,2,2,2,3,3]])


def test_signalplot(bedgraph, regions_10b, true_signal):
    plotter = SignalPlot(10, 10, do_normalize=False)
    signal = plotter([("chr1", bedgraph)], {"chr1": regions_10b})
    assert np.all(signal["y"].values==true_signal)

def test_tssplot(bedgraph, regions_10b, true_signal):
    plotter = TSSPlot(10, 10, do_normalize=False)
    mids = (regions_10b.starts+regions_10b.ends)//2
    r = Regions(mids, mids+1, regions_10b.directions)
    signal = plotter([("chr1", bedgraph)], {"chr1": r})
    assert np.all(signal["y"].values==true_signal)

def test_averageplot(bedgraph, regions_10b, true_signal):
    plotter = AveragePlot(8, do_normalize=False)
    r = Regions(regions_10b.starts+3, regions_10b.ends-3, regions_10b.directions)
    signal = plotter([("chr1", bedgraph)], {"chr1": r})
    assert np.all(signal["y"].values==true_signal[1:-1])

def test_heatplot(bedgraph, regions_10b, true_matrix):
    plotter = HeatPlot(10, 10, do_normalize=False, aspect_ratio=3/10)
    signal = plotter([("chr1", bedgraph)], {"chr1": regions_10b})
    true = pd.DataFrame(np.array(true_matrix, dtype="float"))
    true.index=np.arange(3)+1
    true.columns=np.arange(-5, 5)
    assert signal.equals(true)

def test_vplot(bedgraph, regions_10b, true_signal):
    plotter = VPlot(12, 12, do_normalize=False)
    signal = plotter([("chr1", bedgraph)], {"chr1": regions_10b})
    assert np.all(signal.iloc[10].values[1:-1]==true_signal)
