"""Console script for bdgtools."""
import sys
import click
from pathlib import PurePath

from .io import read_bedgraph, read_bedfile
from .aggregateplot import *
from .plotter import plot, join_line_plots, join_matrix_plots
plot_types = {"v": VPlot, "average": AveragePlot, "heat": HeatPlot, "tss": TSSPlot}


@click.command()
@click.argument("plot_type", type=click.Choice(plot_types.keys()))
@click.argument("bedgraph", type=click.File("r"))
@click.argument("bedfile", type=click.File("r"))
@click.option("-o", "--out_im", "out_im", type=click.File("wb"), help="Path to output figure")
@click.option("-od", "--out_data", "out_data", type=click.File("wb"), help="Path to pickle of figure")
@click.option("-w", "--width", "figure_width", default=2000, help="Figure width")
@click.option("-rs", "--regionsize", "region_size", type=int, help="Genomic region size")
def main(plot_type, bedgraph, bedfile, out_im, out_data, figure_width, region_size):
    assert plot_type in plot_types
    bedgraphs = read_bedgraph(bedgraph)
    regions = read_bedfile(bedfile)
    f = plot_types[plot_type](figure_width=figure_width, region_size=region_size)
    fig = f(bedgraphs, regions)
    plot(fig, f, save_path=out_im, show=out_im is None and out_data is None)
    return 0

@click.command()
@click.argument("data_files", nargs=-1, type=click.File("rb"))
@click.option("-o", "--out_im", "out_im", type=click.File("wb"))
def join_figs(data_files, out_im):
    figs = [pd.read_pickle(df) for df in data_files]
    names = [PurePath(df.name).stem for df in data_files]
    if "x" in figs[0]:
        axes = join_line_plots(figs)
        plt.legend(names)
    else:
        axes = join_matrix_plots(figs)
        for ax, name in zip(axes, name):
            ax.set_title(name)
    if out_im is not None:
        plt.savefig(out_im)
    else:
        plt.show()


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
