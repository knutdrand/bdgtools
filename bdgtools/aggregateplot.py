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
            self._update_chromosome(chrom, bedgraph, self._transform_regions(regions[chrom]))
        return self._finalize()

    def get_x_axis(self):
        return np.arange(self._figure_width)*self._region_size//self._figure_width-self._region_size//2

    def _transform_regions(self, regions):
        return regions

    def _pre_process(self, bedgraphs, regions):
        pass

    def _finalize(self):
        values = np.cumsum(self._diffs, axis=-1)/np.maximum(self._row_counts, 1)[:, None]
        if self._do_normalize:
            values/=(self._coverage/1000000)
        table = pd.DataFrame(values)
        table.index = self.get_y_axis()
        table.columns = self.get_x_axis()
        return table
        
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
    xlabel="Distance from center"
    ylabel="Domain size"
    _aspect_ratio=1
    _region_size = None


    def _pre_process(self, bedgraphs, regions):
        if self._region_size is not None:
            return
        self._region_size = max(np.max(r.sizes()) for r in regions.values())
        print(self._region_size)


    def get_y_axis(self):
        return np.arange(self._row_counts.size)*self._region_size//self._row_counts.size

    def _update_chromosome(self, chrom, bedgraph, regions):
        rows = ((regions.ends-regions.starts)/self._region_size*self._figure_shape[0]).astype("int")
        mask = rows<self._figure_shape[0]
        mids = (regions.ends[mask]+regions.starts[mask])//2
        new_regions = Regions(mids-self._region_size//2, mids+self._region_size//2, regions.directions[mask])
        signals = new_regions.get_signals(bedgraph).scale_x(self._figure_width).update_dense_diffs(self._diffs, rows[mask])
        # bedgraph.extract_regions(new_regions)
        rows, counts = np.unique(rows[mask], return_counts=True)
        self._row_counts[rows] += counts

    def _finalize(self):
        values = np.cumsum(self._diffs, axis=-1)/np.maximum(self._row_counts, 1)[:, None]
        marked_indices = np.flatnonzero(self._row_counts)
        for pre, post in zip(marked_indices[:-1], marked_indices[1:]):
            D = post-pre
            if D==1:
                continue
            values[pre:post] = ((D-np.arange(D))[:, None]*values[pre]+np.arange(D)[:, None]*values[post])/D
        
        if self._do_normalize:
            values/=(self._coverage/1000000)
        table = pd.DataFrame(values)
        table.index = self.get_y_axis()
        table.columns = self.get_x_axis()
        return table

class HeatPlot(MatrixPlot):
    _aspect_ratio=2
    _region_size=100000
    xlabel="Distance from center"
    ylabel="Rank(domainsize)"

    def get_y_axis(self):
        return np.cumsum(self._row_counts)

    def _get_y_coords(self, regions):
        sizes = [r.ends-r.starts for _, r in regions.items()]
        offsets = np.cumsum([0]+[len(s) for s in sizes])
        ranks = get_ranks(np.concatenate(sizes))
        y_coords = (ranks*self._figure_shape[0])//ranks.size
        return {chrom: y_coords[offsets[i]:offsets[i+1]] for i, chrom in enumerate(regions)}

    def _pre_process(self, bedgraphs, regions):
        self._y_coords = self._get_y_coords(regions)

    def _transform_regions(self, regions):
        mids = (regions.ends+regions.starts)//2
        return Regions(mids-self._region_size//2, mids+self._region_size//2, regions.directions)

    def _update_chromosome(self, chrom, bedgraph, regions):
        y_coords = self._y_coords[chrom]
        #signals = bedgraph.extract_regions(regions)
        signals = regions.get_signals(bedgraph)
        signals.scale_x(self._figure_width).update_dense_diffs(self._diffs, y_coords)
        rows, counts = np.unique(y_coords, return_counts=True)
        self._row_counts[rows] += counts

class SignalPlot(AggregatePlot):
    xlabel="Fraction of region"
    ylabel="~FPKM"
    def _update_chromosome(self, chrom, bedgraph, regions):
        signals = regions.get_signals(bedgraph).scale_x(self._figure_width)
        # signals = bedgraph.extract_regions(regions)
        signals.sum(axis=1).update_dense_diffs(self._diffs)
        self._row_counts += regions.starts.size

    def _finalize(self):
        values = np.cumsum(self._diffs, axis=-1)/np.maximum(self._row_counts, 1)
        if self._do_normalize:
            values/=(self._coverage/1000000)
        return pd.DataFrame({"x":self.get_x_axis(), "y": values})

    def get_x_axis(self):
        return np.linspace(-0.5, 0.5, self._figure_width)
        
class TSSPlot(SignalPlot):
    _region_size=2000
    xlabel="Distance from TSS"

    def get_x_axis(self):
        return np.arange(self._figure_width)*self._region_size//self._figure_width-self._region_size//2

    def _transform_regions(self, regions):
        return expand(regions, self._region_size//2, self._region_size//2)

class BorderPlot(TSSPlot):
    xlabel="Distance from Border"
    def _transform_regions(self, regions):
        reversed_regions = Regions(regions.starts, regions.ends, -regions.directions)
        return expand(Regions.concatenate([regions, reversed_regions]) , self._region_size//2, self._region_size//2)

class AveragePlot(SignalPlot):
    def get_x_axis(self):
        return np.linspace(-2, 2, self._figure_width)

    def _transform_regions(self, regions):
        sizes = regions.sizes()
        return Regions(regions.starts-sizes//2, regions.ends+sizes//2, regions.directions)

class MetaGenePlot(SignalPlot):
    def _pre_process(self, bedgraphs, genes_dict):
        utr_l_size, cds_size, utr_r_size = (0, 0, 0)
        for genes in genes_dict.values():
            utr_l_size += np.sum(genes._coding_regions.starts)
            cds_size += np.sum(genes._coding_regions.sizes())
            utr_r_size += np.sum(genes.sizes()-genes._coding_regions.ends)
        sizes = (utr_l_size, cds_size, utr_r_size)
        self._region_sizes = np.array(sizes)*self._figure_width//sum(sizes)
        self._region_sizes[-1]+=self._figure_width-np.sum(self._region_sizes)

    def _finalize(self):
        df = super()._finalize()
        df["region"] = ["utr_l"]*self._region_sizes[0] + ["cds"]*self._region_sizes[1] + ["utr_r"]*self._region_sizes[2]
        return df

    def _update_chromosome(self, chrom, bedgraph, regions):
        regions.get_signals(bedgraph).piecewise_scale(
            [np.zeros_like(regions._coding_regions.starts),
             regions._coding_regions.starts,
             regions._coding_regions.ends,
             regions.sizes()], self._region_sizes).sum(axis=1).update_dense_diffs(self._diffs)
        # 
        #     
        # 
        # utr_l, cds, utr_r = regions.get_signals(bedgraph)
        # 
        # for signals, diffs in zip((utr_l, cds, utr_r), self._split_diffs):
        #     print(signals)
        #     signals.scale_x(diffs.shape[-1]).sum(axis=1).update_dense_diffs(diffs)
        # self._row_counts += regions.starts.size
