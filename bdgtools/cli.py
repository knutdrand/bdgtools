"""Console script for bdgtools."""
import sys
import click
import gzip
from pathlib import PurePath

from .io import read_bedgraph, read_bedfile, read_refseq, read_large_bedfile, write_bedgraph
from .aggregateplot import *
from .coverage import get_coverage
from .plotter import plot, join_plots
plot_types = {"v": VPlot, "average": AveragePlot, "heat": HeatPlot, "tss": TSSPlot, "signal": SignalPlot,
              "metagene": MetaGenePlot, "border": BorderPlot}

@click.command()
@click.argument("plot_type", type=click.Choice(plot_types.keys()))
@click.argument("bedgraph", type=click.Path())
@click.argument("bedfile", type=click.File("r"))
@click.option("-o", "--out_im", "out_im", type=click.File("wb"), help="Path to output figure")
@click.option("-od", "--out_data", "out_data", type=click.File("wb"), help="Path to pickle of figure")
@click.option("-w", "--width", "figure_width", default=2000, help="Figure width")
@click.option("-rs", "--regionsize", "region_size", type=int, help="Genomic region size")
def do_plot(plot_type, bedgraph, bedfile, out_im, out_data, figure_width, region_size):
    bedgraphs = read_bedgraph(bedgraph)
    regions = read_bedfile(bedfile)
    f = plot_types[plot_type](figure_width=figure_width, region_size=region_size)
    fig = f(bedgraphs, regions)
    plot(fig, f, save_path=out_im, show=out_im is None and out_data is None)
    if out_data is not None:
        fig.to_pickle(out_data)
    return 0

@click.group()
def main():
    return 0

@main.command()
@click.argument("plot_type", type=click.Choice(plot_types.keys()))
@click.argument("data_files", nargs=-1, type=click.File("rb"))
@click.option("-o", "--out_im", "out_im", type=click.File("wb"))
@click.option("-n", "--name", "name", default="")
def joinfigs(plot_type, data_files, out_im, name):
    figs = [pd.read_pickle(df) for df in data_files]
    names = [PurePath(df.name).stem for df in data_files]
    click.echo("Joining figures from %s" % " ".join(names))
    cls = plot_types[plot_type]
    name = cls.__name__+ ":" + name
    join_plots(figs, names, cls, save_path=out_im, show=out_im is None, name=name)


@main.command()
@click.argument("plot_type", type=click.Choice(plot_types.keys()))
@click.argument("bedgraph", type=click.Path())
@click.argument("genefile", type=click.Path())
@click.option("-o", "--out_im", "out_im", type=click.File("wb"), help="Path to output figure")
@click.option("-od", "--out_data", "out_data", type=click.File("wb"), help="Path to pickle of figure")
@click.option("-w", "--width", "figure_width", default=2000, help="Figure width")
@click.option("-rs", "--regionsize", "region_size", type=int, help="Genomic region size")
def geneplot(plot_type, bedgraph, genefile, out_im, out_data, figure_width, region_size):
    bedgraphs = read_bedgraph(bedgraph)
    regions = read_refseq(genefile)
    f = plot_types[plot_type](figure_width=figure_width, region_size=region_size)
    fig = f(bedgraphs, regions)
    plot(fig, f, save_path=out_im, show=out_im is None and out_data is None)
    if out_data is not None:
        fig.to_pickle(out_data)
    return 0

@main.command()
@click.argument("bedfile", type=click.Path())
@click.option("-o", "--outfile", "outfile", type=click.File("w"), help="Path to bedgraph file")
def bed2bdg(bedfile, outfile):
    bedfile = read_large_bedfile(gzip.open(bedfile, "rt"))
    bedgraphs = ((chrom, get_coverage(regions)) for chrom, regions in bedfile)
    write_bedgraph(bedgraphs, outfile)
    return 0


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
