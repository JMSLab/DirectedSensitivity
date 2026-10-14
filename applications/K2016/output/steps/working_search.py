import sys, numpy as np
sys.path.insert(0, '.')
from basis import BaselinePoint, Evaluator
ev = Evaluator(); base = BaselinePoint()
rng = np.random.default_rng(20261014)
pts=[]
for _ in range(40):
    p=base.copy(); p[0]+=rng.normal(0,3e-5); p[1:3]+=rng.normal(0,3e-4,2); pts.append(p)
rs=ev.Batch(pts)
q=np.array([r['criterion'] for r in rs]); c=np.array([r['target'] for r in rs])
print('min',q.min(),'q at base',2.7972913365287253,'n below',(q<2.7972913365287253).sum())
print(np.percentile(q,[0,10,50,90,100]))
i=q.argmin(); print(pts[i], q[i], c[i]); print('c range', c.min(), c.max())
A=[r for r in ev.records.values() if r['point'][3]==1 and r['point'][4]==1]
b=min(A,key=lambda r:r['criterion']); print('best A overall', b['point'], b['criterion'], b['target'])
