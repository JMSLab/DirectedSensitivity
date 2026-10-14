import sys, numpy as np
sys.path.insert(0,'.')
from basis import Evaluator
from search import ESTIMATE, BUDGET
ev=Evaluator()
S=np.array([0.002,0.04,0.08]); d=np.array([-735.76272189,423.97886236,267.65786572])*S
pts=[];tags=[]
for sign,(phi,lam) in ((1,(1.2,1.0)),(-1,(0.3,2.0))):
    for dq in (2e-4,3.5e-4,5e-4,7e-4):
        t=np.sqrt(dq/417.32)
        p=ESTIMATE.copy(); p[8]=phi; p[10]=np.log(lam); p[:3]+=sign*t*d; pts.append(p); tags.append((sign,dq))
for tag,p,r in zip(tags,pts,ev.Batch(pts)):
    print(tag, np.round(p[:3],5), f"q-qbar {r['criterion']-BUDGET:+.6f} c {r['target']:.3f}")
