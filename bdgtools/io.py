import logging
from itertools import chain, groupby
from operator import itemgetter
import pandas as pd
import numpy as np
from .regions import Regions
from .splitregions import SplitRegions, Genes
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


def _get_genes(df):
    s = np.array([starts[0] for starts in df["exon_starts"]])
    e = np.array([ends[-1] for ends in df["exon_ends"]])
    mask = df["cds_start"] > s
    mask &= df["cds_end"] > df["cds_start"]
    mask &= e > df["cds_end"]
    if not np.any(mask):
        return None
    df = df.loc[mask]
    directions = np.where(df["direction"].values=="+", 1, -1)
    starts = np.concatenate([s[::d] for s, d in zip(df["exon_starts"].values, directions)])
    ends = np.concatenate([s[::d] for s, d in zip(df["exon_ends"].values, directions)])
    lens = [len(starts) for starts in df["exon_starts"]]
    coding_offsets = np.array([get_coding_offsets(r["exon_starts"], r["exon_ends"],
                                                  r["cds_start"], r["cds_end"], r["direction"])
                               for _, r in df.iterrows()], dtype="int")
    
    offsets=np.cumsum([0]+lens)
    all_directions = np.concatenate([[d]*l for d, l in zip(directions, lens)])
    regions = Regions(starts, ends, all_directions)
    return Genes(regions, offsets, coding_regions=Regions(coding_offsets[:, 0], coding_offsets[:, 1]))

def get_coding_offsets(exon_starts, exon_ends, cds_start, cds_end, direction):
    #assert anno.txStart<anno.cdsStart<anno.cdsEnd<anno.txEnd, anno
    offsets = np.cumsum([0]+[end-start for start, end in zip(exon_starts, exon_ends)])
    idxs = np.searchsorted(exon_starts, [cds_start, cds_end], side="right")-1

    # for idx, pos in zip(idxs, [anno.cdsStart, anno.cdsEnd]):
    #     assert idx<len(anno.exonStarts), (idx, pos, anno)
    #     assert anno.exonStarts[idx]<=pos, (anno, idx, pos)
    #     assert idx==len(anno.exonStarts)-1 or pos < anno.exonStarts[idx+1], (anno, idx, pos)
    t =  tuple(offsets[i]+pos-exon_starts[i] for i, pos in zip(idxs, [cds_start, cds_end]))
    if direction=="-":
        size = sum(exon_ends)-sum(exon_starts)
        return (size-t[1], size-t[0])
    return t

def read_refseq(file_obj):
    get_ints = lambda x: [int(r) for r in x.split(",") if r]
    df = pd.read_table(file_obj, 
                       names=["chrom", "direction", "start", "end","cds_start","cds_end", "exon_starts", "exon_ends"], 
                       usecols=[2,3,4,5,6,7,9,10], 
                       converters={"exon_starts": get_ints, "exon_ends": get_ints})
    grouped = df.sort_values(["chrom", "start"]).groupby("chrom", sort=False)
    d =  {chrom: _get_genes(df) for chrom, df in grouped}
    return {chrom: genes for chrom,  genes in d.items() if genes is not None}
