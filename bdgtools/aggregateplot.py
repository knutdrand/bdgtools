import numpy as np
import pandas as pd
import logging
from .regions import Regions, expand
from .util import get_ranks
log = logging

class AggregatePlot:
    _figure_width=2000
    _region_size=None
    _aspect_ratio=None
    def __init__(self, figure_width=2000, region_size=None, do_normalize=True):
        self._figure_width = figure_width
        self._figure_shape = (figure_width,)
        self._do_normalize = do_normalize
        self._row_counts = 0
        self._coverage = 0
        if region_size is not None:
            self._region_size = region_size

    def __call__(self, bedgraphs, regions):
        self._diffs = np.zeros(self._figure_shape)
        self._pre_process(bedgraphs, regions)
        for chrom, bedgraph in bedgraphs:
            if self._do_normalize:
                self._coverage += bedgraph.sum()
            if chrom not in regions:
                continue
            log.info("Processing %s", chrom)
            self.update_chromosome(chrom, bedgraph, self.transform_regions(regions[chrom]))
        return self.finalize()

    def transform_regions(self, regions):
        return regions

    def _pre_process(self, bedgraphs, regions):
        pass

    def finalize(self):
        values = np.cumsum(self._diffs, axis=-1)/np.maximum(self._row_counts, 1)[:, None]
        if self._do_normalize:
            values/=(self._coverage/1000000)
        table = pd.DataFrame(values)
        table.index = self.get_y_axis()
        table.columns = self.get_x_axis()
        return table

    def get_x_axis(self):
        return np.arange(self._figure_width)*self._region_size//self._figure_width-self._region_size//2
        
class MatrixPlot(AggregatePlot):
    def __init__(self, *args, aspect_ratio=None, **kwargs):
        super().__init__(*args, **kwargs)
        if aspect_ratio is not None:
            self._aspect_ratio=aspect_ratio
        self._figure_shape = (int(self._aspect_ratio*self._figure_width), self._figure_width)
        self._row_counts = np.zeros(self._figure_shape[0], dtype="int")

    def get_y_axis(self):
        return np.arange(self._row_counts.size)
    
class VPlot(MatrixPlot):
    _aspect_ratio=1
    _region_size = 50000
    xlabel="Distance from center"
    ylabel="Domain size"

    def update_chromosome(self, chrom, bedgraph, regions):
        rows = ((regions.ends-regions.starts)/self._region_size*self._figure_shape[0]).astype("int")
        mask = rows<self._figure_shape[0]
        mids = (regions.ends[mask]+regions.starts[mask])//2
        bedgraph.extract_regions(Regions(mids-self._region_size//2, mids+self._region_size//2, regions.directions[mask])).scale_x(self._figure_width).update_dense_diffs(self._diffs, rows[mask])
        rows, counts = np.unique(rows[mask], return_counts=True)
        self._row_counts[rows] += counts

    def finalize(self):
        values = np.cumsum(self._diffs, axis=-1)/np.maximum(self._row_counts, 1)[:, None]
        marked_indices = np.flatnonzero(self._row_counts)
        print(np.flatnonzero(self._row_counts==0))
        for pre, post in zip(marked_indices[:-1], marked_indices[1:]):
            D = post-pre
            if D==1:
                continue
            values[pre:post] = ((D-np.arange(D))*values[pre]+np.arange(D)*values[post])/D
        
        if self._do_normalize:
            values/=(self._coverage/1000000)
        table = pd.DataFrame(values)
        table.index = self.get_y_axis()
        table.columns = self.get_x_axis()
        return table

    def get_y_axis(self):
        return np.arange(self._row_counts.size)*self._region_size//self._row_counts.size

class HeatPlot(MatrixPlot):
    _aspect_ratio=2
    _region_size=100000
    xlabel="Distance from center"
    ylabel="Rank(domainsize)"
    def _get_y_coords(self, regions):
        sizes = [r.ends-r.starts for _, r in regions.items()]
        offsets = np.cumsum([0]+[len(s) for s in sizes])
        ranks = get_ranks(np.concatenate(sizes))
        y_coords = (ranks*self._figure_shape[0])//ranks.size
        return {chrom: y_coords[offsets[i]:offsets[i+1]] for i, chrom in enumerate(regions)}

    def _pre_process(self, bedgraphs, regions):
        self._y_coords = self._get_y_coords(regions)

    def transform_regions(self, regions):
        mids = (regions.ends+regions.starts)//2
        return Regions(mids-self._region_size//2, mids+self._region_size//2, regions.directions)

    def update_chromosome(self, chrom, bedgraph, regions):
        y_coords = self._y_coords[chrom]
        signals = bedgraph.extract_regions(regions)
        signals.scale_x(self._figure_width).update_dense_diffs(self._diffs, y_coords)
        rows, counts = np.unique(y_coords, return_counts=True)
        self._row_counts[rows] += counts

    def get_y_axis(self):
        return np.cumsum(self._row_counts)

class SignalPlot(AggregatePlot):
    def update_chromosome(self, chrom, bedgraph, regions):
        print(regions)
        signals = bedgraph.extract_regions(regions)
        if self._region_size != self._figure_width:
            signals = signals.scale_x(self._figure_width)
        signals.sum(axis=1).update_dense_diffs(self._diffs)
        self._row_counts += regions.starts.size

    def finalize(self):
        values = np.cumsum(self._diffs, axis=-1)/np.maximum(self._row_counts, 1)
        if self._do_normalize:
            values/=(self._coverage/1000000)
        return pd.DataFrame({"x":self.get_x_axis(), "y": values})

        
class TSSPlot(SignalPlot):
    _region_size=2000
    xlabel="Distance from TSS"
    ylabel="~FPKM"
    def transform_regions(self, regions):
        return expand(regions, self._region_size//2, self._region_size//2)

class AveragePlot(SignalPlot):
    xlabel="Fraction of domain"
    ylabel="~FPKM"
    def transform_regions(self, regions):
        sizes = regions.ends-regions.starts
        return Regions(regions.starts-sizes//2, regions.ends+sizes//2, regions.directions)

    def get_x_axis(self):
        return np.linspace(-2, 2, self._figure_width)
