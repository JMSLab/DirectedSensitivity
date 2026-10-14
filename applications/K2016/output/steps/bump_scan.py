import sys, numpy as np
sys.path.insert(0, '.')
from basis import Evaluator, Pad
ev = Evaluator()
x0 = Pad([0.288509110693, -2.98192892931e-05, 9.50846758198e-05, 1.0, 1.0, np.log(500.0), np.log(500.0)])
QB = ev(x0)['criterion']; C0 = ev(x0)['target']
pts=[]; tags=[]
for k in (2,6,20,60,200):
    amps=(0.1,0.3) if k<=6 else (0.05,0.15)
    for side in (0,1):
        for a in amps:
            p=x0.copy(); p[7+side]=a; p[9+side]=np.log(k); pts.append(p); tags.append((k,side,a))
for t,r in zip(tags, ev.Batch(pts)):
    print(t, f"dq {r['criterion']-QB:+.6f} dc {r['target']-C0:+.3f} pen {r['penalty']:.5f}")
