import sys, numpy as np
sys.path.insert(0, '.')
from basis import Evaluator, Pad
ev = Evaluator()
x0 = Pad([0.288509110693, -2.98192892931e-05, 9.50846758198e-05, 1.0, 1.0, np.log(500.0), np.log(500.0)])
QB = ev(x0)['criterion']; C0 = ev(x0)['target']
pts=[]; tags=[]
for k in (1,2,4):
    for a in (0.6,1.2,2.5):
        p=x0.copy(); p[8]=a; p[10]=np.log(k); pts.append(p); tags.append((k,a))
for t,r in zip(tags, ev.Batch(pts)):
    print(t, f"dq {r['criterion']-QB:+.6f} dc {r['target']-C0:+.3f} pen {r['penalty']:.5f} ll {-r['log_likelihood']:.6f}")
