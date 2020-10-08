import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

from .aggregateplot import *
sns.set_theme()

def plot(df, cls, save_path=None, show=False):
    if save_path is None and not show:
        return
    plt.figure(figsize=(10, 10))
    if isinstance(cls, SignalPlot):
        if "region" in df:
            p = sns.lineplot(data=df, x="x", y="y", style="region")
        else:
            p = sns.lineplot(data=df, x="x", y="y")
    else:
        p = sns.heatmap(df, cmap="gray_r")
    p.set_xlabel(cls.xlabel)
    p.set_ylabel(cls.ylabel)
    p.set_title(f"{cls.__class__.__name__}")
    if save_path is not None:
        plt.savefig(save_path)
    if show:
        plt.show()


def join_plots(dfs, names, cls, save_path=None, show=False, name=""):
    if issubclass(cls, SignalPlot):
        plt.figure(figsize=(10, 10))
        cur_hue = 0
        for df in dfs:
            if "region" in df:
                p = sns.lineplot(data=df, x="x", y="y", style="region")
            else:
                sns.lineplot(data=df, x="x", y="y")
        plt.xlabel(cls.xlabel)
        plt.ylabel(cls.ylabel)
        plt.legend(names)
        plt.title(name)
    elif issubclass(cls, MatrixPlot):
        cols = int(np.sqrt(len(dfs)-1)+1)
        rows = (len(dfs)-1)//cols+1
        f, axes = plt.subplots(rows, cols, sharex=True, sharey=True, squeeze=False)
        f.set_figheight(axes.shape[0]*5)
        f.set_figwidth(axes.shape[1]*5)
        f.suptitle(name)
        for i, (df, name) in enumerate(zip(dfs, names)):
            ax = axes[i//cols, i%cols]
            sns.heatmap(df, cmap="gray_r", ax=ax)
            ax.set_title(name)
            ax.set_xlabel(cls.xlabel)
        for a in axes[:, 0]:
            a.set_ylabel(cls.ylabel)
    if save_path is not None:
        plt.savefig(save_path)
    if show:
        plt.show()

