import sys, json, numpy as np
sys.path.insert(0,'.')
from basis import Evaluator
from search import ESTIMATE, BUDGET
ev=Evaluator()
S=np.array([0.002,0.04,0.08]); d=np.array([-735.76272189,423.97886236,267.65786572])*S
spec=json.loads(sys.argv[1])
pts=[];tags=[]
for sign,phi,lam,dq in spec:
    t=np.sqrt(dq/417.32)
    p=ESTIMATE.copy(); p[8]=phi; p[10]=np.log(lam); p[:3]+=sign*t*d; pts.append(p); tags.append((sign,phi,lam,dq))
for tag,p,r in zip(tags,pts,ev.Batch(pts)):
    print(tag, f"q-qbar {r['criterion']-BUDGET:+.6f} c {r['target']:.3f}")
print('n', len(ev))
