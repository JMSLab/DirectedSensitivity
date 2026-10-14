import sys, itertools, numpy as np
sys.path.insert(0, '.')
from basis import BaselinePoint, Evaluator
ev = Evaluator(); base = BaselinePoint(); H = np.array([0.005, 0.05, 0.05])
pts = []
for i, j in itertools.combinations(range(3), 2):
    p = base.copy(); p[i] += H[i]; p[j] += H[j]; pts.append(p)
ev.Batch(pts)
# quadratic fit on all A points
A = [r for r in ev.records.values() if r['point'][3] == 1.0 and r['point'][4] == 1.0]
X = np.array([r['point'][:3] for r in A]); y = np.array([r['criterion'] for r in A])
Z = (X - base[:3]) / H
cols = [np.ones(len(Z))] + [Z[:, i] for i in range(3)] + [Z[:, i] * Z[:, j] for i in range(3) for j in range(i, 3)]
M = np.column_stack(cols); coef = np.linalg.lstsq(M, y, rcond=None)[0]
print(len(A), 'resid', np.abs(M @ coef - y).max())
g = coef[1:4]; Hs = np.zeros((3, 3)); k = 4
for i in range(3):
    for j in range(i, 3):
        if i == j: Hs[i, i] = 2 * coef[k]
        else: Hs[i, j] = Hs[j, i] = coef[k]
        k += 1
print('eig', np.linalg.eigvalsh(Hs))
step = -np.linalg.solve(Hs, g); print('newton (scaled)', step, 'pred', 0.5 * g @ step)
cands = [base.copy() for _ in range(4)]
for c, t in zip(cands, (0.5, 0.8, 1.0, 1.2)):
    c[:3] += t * step * H
for c, r in zip(cands, ev.Batch(cands)):
    print(np.round(c[:3], 6), f"{r['criterion']:.7f} {r['target']:.3f}")
