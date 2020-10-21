from collections import namedtuple
import numpy as np
Region = namedtuple("Region", ["start", "end", "direction"])

def get_holes(regions):
    starts = regions.ends[:-1]
    ends = regions.starts[1:]
    return Regions(starts, ends)

class Regions:
    def __init__(self, starts, ends, directions=1):
        self.starts = np.asanyarray(starts)
        self.ends = np.asanyarray(ends)
        assert np.all(self.starts>=0), self.starts
        assert np.all(self.ends>=0), self.ends
        assert np.all(self.ends>starts), (self.starts, self.ends)

        if isinstance(directions, int) and directions == 1:
            self.directions=np.ones_like(self.starts)
        else:
            self.directions = np.asanyarray(directions)

    def __iter__(self):
        if isinstance(self, int) and self.directions == 1:
            return (Region(s, e, 1) for s, e in zip(self.starts, self.ends))
        return (Region(*t) for t in zip(self.starts, self.ends, self.directions))

    def __eq__(self, other):
        t = np.all(self.starts==other.starts)
        t &= np.all(self.ends==other.ends)
        return t & np.all(self.directions==other.directions)

    def __repr__(self):
        return f"Regions({self.starts}, {self.ends}, {self.directions})"

    def sizes(self):
        return self.ends-self.starts

    def get_signals(self, bedgraph):
        return bedgraph.extract_regions(self)

    @classmethod
    def concatenate(cls, regions_list):
        starts = np.concatenate([r.starts for r in regions_list])
        ends = np.concatenate([r.ends for r in regions_list])
        directions = np.concatenate([r.directions for r in regions_list])
        args = starts.argsort(kind="mergesort")
        return cls(starts[args], ends[args], directions[args])

def expand(regions, upstream, downstream):
    centers = np.where(regions.directions==1, regions.starts, regions.ends-1)
    starts = centers-np.where(regions.directions==1, upstream, downstream)
    ends = centers+np.where(regions.directions==1, downstream, upstream)
    args = starts.argsort(kind="mergesort")
    return Regions(starts[args], ends[args], regions.directions[args])
