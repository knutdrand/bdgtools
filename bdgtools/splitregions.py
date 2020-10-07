import numpy as np
from .regions import Regions
class SplitRegions:
    def __init__(self, regions, offsets):
        self._regions = regions
        self._offsets = offsets

    def get_signals(self, bedgraph):
        signals = bedgraph.extract_regions(self._regions)
        return signals.join_rows(self._offsets)
        
    def sizes(self):
        return np.diff(np.insert(np.cumsum(self._regions.sizes()), 0, 0)[self._offsets])
        
    @property
    def starts(self):
        return self._regions.starts[self._offsets[:-1]]

    @property
    def ends(self):
        return self._regions.ends[self._offsets[:-1]]

class Genes(SplitRegions):
    def __init__(self, *args, coding_regions):
        super().__init__(*args)
        assert np.all(coding_regions.ends<=self.sizes()), (coding_regions.ends[coding_regions.ends>=self.sizes()], self.sizes()[coding_regions.ends>=self.sizes()])
        self._coding_regions = coding_regions

    def __repr__(self):
        return f"Genes({self._regions}, {self._offsets}, {self._coding_regions})"

    def __eq__(self, other):
        t = self._regions==other._regions
        t &= np.all(self._offsets==other._offsets)
        t &= self._coding_regions == other._coding_regions
        return t
    
    # def get_signals(self, bedgraph):
    #     signals = super().get_signals(bedgraph)
    #     utr_l = Regions(np.zeros_like(self._coding_regions.starts),
    #                     self._coding_regions.starts)
    #     utr_r = Regions(self._coding_regions.ends,
    #                     self.sizes())
    #     return tuple(signals.extract_regions(regions)
    #                  for regions in (utr_l, self._coding_regions, utr_r))
