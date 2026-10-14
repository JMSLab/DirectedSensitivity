import sys, numpy as np
sys.path.insert(0, '.')
from basis import BaselinePoint, Evaluator
ev = Evaluator(); base = BaselinePoint()
pts=[]
for d in (-1e-3, -1e-4, -1e-6, 1e-6, 1e-4, 1e-3, 3e-3, -3e-3):
    p=base.copy(); p[0]+=d; pts.append(p)
for p,r in zip(pts, ev.Batch(pts)):
    print(f"{p[0]-base[0]:+.0e} {r['criterion']:.7f} ll {-r['log_likelihood']:.7f} pen {r['penalty']:.6f} c {r['target']:.3f} match {r['match_rate']:.5f}")
