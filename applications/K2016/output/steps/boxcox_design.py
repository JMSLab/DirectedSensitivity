import sys, numpy as np
from scipy.stats import qmc
sys.path.insert(0, '.')
from basis import Evaluator
ev = Evaluator()
center = np.array([0.288509110693, -2.98192892931e-05, 9.50846758198e-05, 1.0, 1.0, np.log(500.0), np.log(500.0)])
lo = np.array([-0.004, -0.08, -0.12, 0.90, 0.75]); hi = np.array([0.004, 0.08, 0.12, 1.04, 1.10])
u = qmc.LatinHypercube(d=5, seed=7).random(50)
pts=[]
for row in u:
    p=center.copy(); v=lo+row*(hi-lo); p[0:3]+=v[0:3]; p[3:5]=v[3:5]; pts.append(p)
rs=ev.Batch(pts)
print(len(ev))
