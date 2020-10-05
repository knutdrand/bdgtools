import numpy as np


def get_ranks(array):
    args_tmp = np.argsort(array)
    args = np.empty_like(args_tmp)
    args[args_tmp] = np.arange(len(args))
    return args
