"""Exploratory derivatives at the baseline, logged to evaluations.jsonl."""
import sys
import numpy as np
from basis import BaselinePoint, Evaluator

STEPS = {0: 0.005, 1: 0.05, 2: 0.05}


def Main():
    evaluator = Evaluator()
    base = BaselinePoint()
    points = [base]
    for i, h in STEPS.items():
        for sign in (1, -1):
            p = base.copy(); p[i] += sign * h; points.append(p)
    for scale in (5.0, 50.0, 500.0):
        for side in (3, 4):
            p = base.copy(); p[5] = p[6] = np.log(scale); p[side] = 0.9; points.append(p)
    for p, r in zip(points, evaluator.Batch(points)):
        print(np.round(p, 4), f"{r['criterion']:.7f} {r['target']:.3f} pen {r['penalty']:.5f}")
    print('evaluations', len(evaluator))


if __name__ == '__main__':
    Main()
