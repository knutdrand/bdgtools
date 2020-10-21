from .bedgraph import BedGraph
import numpy as np


def get_coverage(regions):
    all_indices = np.concatenate((regions.starts, regions.ends))
    args = np.argsort(all_indices, kind="mergesort")
    diffs = np.where(args<regions.starts.size, 1, -1)
    values = np.cumsum(diffs)
    sorted_indices = all_indices[args]
    changes = np.append(sorted_indices[:-1] != sorted_indices[1:], True)
    values = values[changes]
    value_changes = np.insert(values[1:]!=values[:-1], 0, True)
    indices = sorted_indices[changes][value_changes]
    values = values[value_changes]
    if indices[0] != 0:
        indices = np.insert(indices, 0, 0)
        values = np.insert(values, 0, 0)
    
    return BedGraph(indices, values)
                    
