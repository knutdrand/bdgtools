import logging
from itertools import chain, groupby
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

def _get_bedgraph(chunks):
    chunks = list(chunks)
    cur_chrom = chunks[0]["chrom"].iloc[0]
    starts = np.concatenate([c["start"].values for c in chunks])
    ends = np.concatenate([c["end"].values for c in chunks])
    assert starts[0] == 0, f"Bedgraph does not start on 0 on {cur_chrom}"
    assert np.all(starts[1:] == ends[:-1]), f"Begraph is not continous on {cur_chrom}, {starts[1:]}, {ends[:-1]}\n{np.flatnonzero(starts[1:]!=ends[:-1])}, {starts.size}"
    log.info("Read chromosome", cur_chrom)
    return BedGraph(starts,
                    np.concatenate([c["value"].values for c in chunks]),
                    chunks[-1]["end"].values[-1])

def _split_chunk_on_chromosomes(chunk):
    changes = list(np.flatnonzero(chunk["chrom"].values[1:]!=chunk["chrom"].values[:-1])+1)
    return (chunk.iloc[start_idx:end_idx] for start_idx, end_idx in
            zip([0] + changes, changes + [len(chunk.index)]))

def read_bedgraph(file_obj, size_hint=1000000):
    reader = pd.read_table(file_obj, names=["chrom", "start", "end", "value"], usecols=[0, 1, 2, 3], chunksize=size_hint)
    chunks = chain.from_iterable(_split_chunk_on_chromosomes(chunk) for chunk in reader)
    return ((chrom, _get_bedgraph(group)) for chrom, group in 
            groupby(chunks, lambda x: x.iloc[0]["chrom"]))
