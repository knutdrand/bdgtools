import logging
from itertools import chain, groupby
from operator import itemgetter
import pandas as pd
import numpy as np
from .regions import Regions
from .splitregions import SplitRegions
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

def read_bedgraph(file_obj, size_hint=1000000):
    reader = pd.read_table(file_obj, names=["chrom", "start", "end", "value"], usecols=[0, 1, 2, 3], chunksize=size_hint)
    grouped = groupby(chain.from_iterable(chunk.groupby("chrom", sort=False) for chunk in reader), 
                      itemgetter(0))
    grouped = ((chrom, map(itemgetter(1),  group)) for chrom, group in grouped)
    return ((chrom, _get_bedgraph(group)) for chrom, group in grouped)


def _get_split_regions(df):
    directions = np.where(df["direction"].values=="+", 1, -1)
    starts = np.concatenate([s[::d] for s, d in zip(df["exon_starts"].values, directions)])
    #       (df["exon_starts"].values)[::d])
    ends = np.concatenate([s[::d] for s, d in zip(df["exon_ends"].values, directions)])
    # ends = np.concatenate(df["exon_ends"].values)
    lens = [len(starts) for starts in df["exon_starts"]]
    offsets=np.cumsum([0]+lens)
    all_directions = np.concatenate([[d]*l for d, l in zip(directions, lens)])
    regions = Regions(starts, ends, all_directions)
    return SplitRegions(regions, offsets)

def read_refseq(file_obj):
    get_ints = lambda x: [int(r) for r in x.split(",") if r]
    df = pd.read_table(file_obj, 
                       names=["chrom", "direction", "start", "end", "exon_starts", "exon_ends"], 
                       usecols=[2,3,4,5,9,10], 
                       converters={"exon_starts": get_ints, "exon_ends": get_ints})
    grouped = df.sort_values(["chrom", "start"]).groupby("chrom", sort=False)
    return {chrom: _get_split_regions(df) for chrom, df in grouped}
