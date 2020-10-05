import seaborn as sns
import matplotlib.pyplot as plt

from .aggregateplot import *
sns.set_theme()

def plot(df, cls, save_path=None, show=False):
    plt.figure(figsize=(10, 10))
    if isinstance(cls, SignalPlot):
        p = sns.lineplot(data=df, x="x", y="y")
    else:
        p = sns.heatmap(df, cmap="gray_r")
    p.set_xlabel(cls.xlabel)
    p.set_ylabel(cls.ylabel)
    p.set_title(f"{cls.__class__.__name__}")
    if save_path is not None:
        plt.savefig(out_im)
    if show:
        plt.show()


def join_line_plots(figs):
    plt.figure(figsize=(10, 10))
    return [sns.lineplot(data=fig, x="x", y="y") for fig in figs]

def join_matrix_plots(figs):
    f, axes = plt.subplots(1, len(figs))
    return [sns.heatmap(fig, ax=ax) for ax, fig in zip(axes, figs)]
    
