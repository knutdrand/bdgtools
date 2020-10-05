from collections import namedtuple
import numpy as np
Region = namedtuple("Region", ["start", "end", "direction"])

def get_holes(regions):
    starts = regions.ends[:-1]
    ends = regions.starts[1:]
    return Regions(starts, ends)

class Regions:
    def __init__(self, starts, ends, directions=1):
        # assert np.all(np.diff(starts)>0), starts[:-1][np.diff(starts)>0]
        # assert np.all(np.diff(ends)>0)
        self.starts = np.asanyarray(starts)
        self.ends = np.asanyarray(ends)
        if isinstance(directions, int) and directions == 1:
            self.directions=np.ones_like(self.starts)
        else:
            self.directions = np.asanyarray(directions)

    def __iter__(self):
        if isinstance(self, int) and self.directions == 1:
            return (Region(s, e, 1) for s, e in zip(self.starts, self.ends))
        return (Region(*t) for t in zip(self.starts, self.ends, self.directions))

    def __repr__(self):
        return f"Regions({self.starts}, {self.ends})"

def expand(regions, upstream, downstream):
    centers = np.where(regions.directions==1, regions.starts, regions.ends-1)
    starts = centers-np.where(regions.directions==1, upstream, downstream)
    ends = centers+np.where(regions.directions==1, downstream, upstream)
    args = starts.argsort()
    return Regions(starts[args], ends[args], regions.directions[args])
