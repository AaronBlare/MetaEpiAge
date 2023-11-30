import functools
import numpy as np


def conjunction(conditions):
    return functools.reduce(np.logical_and, conditions)


def disjunction(conditions):
    return functools.reduce(np.logical_or, conditions)
