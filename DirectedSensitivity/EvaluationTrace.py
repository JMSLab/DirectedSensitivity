"""Bookkeeping for the criterion and target evaluations a sensitivity problem performs."""
import numpy as np

from .Numerics import CopyValue

COUNTED_NAMES = ('criterion', 'target')


class EvaluationBudgetExceeded(RuntimeError):
    """Raised when a watched callable reaches its declared limit."""


class EvaluationTrace:
    """Evaluation counts of the criterion and target, with optional per-name limits."""

    def __init__(self, limits=None):
        self.counts = dict.fromkeys(COUNTED_NAMES, 0)
        self.limits = dict(limits or {})

    def Watch(self, name, function):
        """``function`` cached by point, with each first evaluation at a point counted under ``name``."""
        if function is None:
            return None
        cache = {}

        def watched(point):
            key = _Key(point)
            if key in cache:
                return CopyValue(cache[key])
            self.Charge(name)
            value = function(point)
            if key is not None:
                cache[key] = CopyValue(value)
            return value

        return watched

    def Count(self, name):
        return self.counts[name]

    def Snapshot(self):
        return dict(self.counts)

    def Restore(self, snapshot):
        """Adopt the counts in a snapshot, so one budget spans processes."""
        self.counts = {name: int(snapshot.get(name, 0)) for name in COUNTED_NAMES}
        return self

    def Charge(self, name):
        """Count one evaluation under ``name``, raising once its limit is reached."""
        limit = self.limits.get(name)
        if limit is not None and self.counts[name] >= limit:
            raise EvaluationBudgetExceeded(
                f"{name} reached its limit of {limit} evaluations")
        self.counts[name] += 1


def _Key(point):
    try:
        array = np.asarray(point, dtype=float)
    except (TypeError, ValueError):
        return None
    return array.shape, array.tobytes()
