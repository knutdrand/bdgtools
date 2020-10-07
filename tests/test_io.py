import io

from bdgtools.io import read_bedgraph, read_bedfile, read_refseq
from bdgtools import BedGraph, Regions
from bdgtools.splitregions import Genes

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

def test_read_directed_bedfile():
    lines = ["chr1\t0\t10\t.\t.\t+",
             "chr1\t10\t25\t.\t.\t-",
             "chr1\t25\t35\t.\t.\t+",
             "chr2\t0\t5\t.\t.\t-",
             "chr2\t5\t10\t.\t.\t+"]
    f = io.StringIO("\n".join(lines))
    bedfile = read_bedfile(f)
    assert bedfile == {"chr1": Regions([0, 10, 25], [10, 25, 35], [1, -1, 1]),
                       "chr2": Regions([0, 5], [5, 10], [-1, 1])}

def test_read_genes():
    lines = ["1366	NM_026243	chr10	+	102374436	102391468	102378157	102389363	4	102374436,102378151,102385005,102388221,	102374644,102378304,102385153,102391468,	0	Mgat4c	cmpl	cmpl	-1,0,0,1,",
             "171	NM_172553	chr10	-	103007846	103028777	103009187	103028606	4	103007846,103022176,103025134,103028380,	103009508,103022305,103025439,103028777,	0	Alx1	cmpl	cmpl	0,0,1,0,"]
    f = io.StringIO("\n".join(lines))
    genes = read_refseq(f)["chr10"]
    exon_starts = [102374436,102378151,102385005,102388221]+[103007846,103022176,103025134,103028380][::-1]
    exon_ends = [102374644,102378304,102385153,102391468] + [103009508,103022305,103025439,103028777][::-1]
    cd_starts = [102378157-102378151 + 102374644-102374436,
                 103028777-103028606]
    cd_ends = [509+1142, 321+831]
    true_genes = Genes(Regions(exon_starts, exon_ends, [1]*4+[-1]*4), [0, 4, 8],
                       coding_regions=Regions(cd_starts, cd_ends))
    assert genes == true_genes
