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


def join_plots(dfs, names, cls, save_path=None, show=False):
    if issubclass(cls, SignalPlot):
        plt.figure(figsize=(10, 10))
        for df in dfs:
            sns.lineplot(data=df, x="x", y="y")
        plt.xlabel(cls.xlabel)
        plt.ylabel(cls.ylabel)
        plt.legend(names)
    elif issubclass(cls, MatrixPlot):
        f, axes = plt.subplots(1, len(dfs))
        for ax, df, name in zip(axes, dfs, names):
            sns.heatmap(df, cmap="gray_r", ax=ax)
            ax.set_title(name)
            ax.set_xlabel(cls.xlabel)
        axes[0].set_ylabel(cls.ylabel)
    if save_path is not None:
        plt.savefig(save_path)
    if show:
        plt.show()

