import io

from bdgtools.io import read_bedgraph, read_bedfile
from bdgtools import BedGraph, Regions

def test_read_bedgraph():
    lines = ["chr1\t0\t10\t0",
             "chr1\t10\t25\t1",
             "chr1\t25\t35\t10",
             "chr2\t0\t5\t0",
             "chr2\t5\t10\t2"]
    f = io.StringIO("\n".join(lines))
    bedgraphs = list(read_bedgraph(f))
    assert bedgraphs == [("chr1", BedGraph([0, 10, 25], [0, 1, 10])),
                         ("chr2", BedGraph([0, 5], [0, 2]))]

def test_read_bedfile():
    lines = ["chr1\t0\t10",
             "chr1\t10\t25",
             "chr1\t25\t35",
             "chr2\t0\t5",
             "chr2\t5\t10"]
    f = io.StringIO("\n".join(lines))
    bedfile = read_bedfile(f)
    assert bedfile == {"chr1": Regions([0, 10, 25], [10, 25, 35]),
                       "chr2": Regions([0, 5], [5, 10])}
