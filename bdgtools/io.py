import logging
import pandas as pd
import numpy as np
from .regions import Regions
from .bedgraph import BedGraph

log = logging

def _peek_line(f):
    pos = f.tell()
    line = f.readline()
    f.seek(pos)
    return line

def read_bedfile(file_obj):
    n_cols = len(_peek_line(file_obj).split("\t"))
    assert n_cols >=3, n_cols
    if n_cols < 6:
        table = pd.read_table(file_obj, names=["chrom", "start", "end"], usecols=[0, 1, 2])
    else:
        table = pd.read_table(file_obj, names=["chrom", "start", "end", "direction"], usecols=[0, 1, 2, 5])
    table = table.sort_values(["chrom", "start"])
    changes = np.flatnonzero(table["chrom"].values[:-1] != table["chrom"].values[1:])+1
    changes = np.concatenate(([0], changes, [table["chrom"].values.size]))
    chrom_split = (table.iloc[start:end] for start, end in zip(changes[:-1], changes[1:]))
    r =  {t["chrom"].iloc[0]: Regions(t["start"].values, t["end"].values,
                                      np.where(t["direction"].values=="+", 1, -1) if n_cols>=6 else 1)
          for t in chrom_split}
    return r

def read_bedgraph(file_obj, size_hint=1000000):
    cur_chrom=None
    reader = pd.read_table(file_obj, names=["chrom", "start", "end", "value"], usecols=[0, 1, 2, 3], chunksize=size_hint)
    chunks = []
    for chunk in reader:
        while chunk["chrom"].iloc[-1] != cur_chrom:
            idx = np.argmax(chunk["chrom"].values!=cur_chrom)
            chunks.append(chunk.iloc[:idx])
            if cur_chrom is not None:
                log.info("Yielding chrom %s", cur_chrom)
                yield cur_chrom, BedGraph(np.concatenate([c["start"].values for c in chunks]),
                                          np.concatenate([c["value"].values for c in chunks]),
                                          chunks[-1]["end"].values[-1])
                chunks = []
            chunk = chunk.iloc[idx:]
            cur_chrom = chunk["chrom"].iloc[0]
        chunks.append(chunk)
    yield cur_chrom, BedGraph(np.concatenate([c["start"].values for c in chunks]),
                              np.concatenate([c["value"].values for c in chunks]),
                              chunks[-1]["end"].values[-1])
