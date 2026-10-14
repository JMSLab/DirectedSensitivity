"""Numerical helpers shared across the package."""
import numpy as np


def CentralDifference(function, point, index, relative_step):
    """Central-difference partial of ``function`` in coordinate ``index``, with a step
    relative to that coordinate's magnitude."""
    spacing = relative_step * max(1.0, abs(point[index]))
    forward, backward = point.copy(), point.copy()
    forward[index] += spacing
    backward[index] -= spacing
    return (function(forward) - function(backward)) / (2.0 * spacing)


def CopyValue(value):
    return value.copy() if isinstance(value, np.ndarray) else value
