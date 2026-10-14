import sys, numpy as np
sys.path.insert(0,'.')
from basis import Evaluator
from search import ESTIMATE, BUDGET
ev=Evaluator(); pts=[]; tags=[]
def P(phi,lam,dg=0,df=0,da=0):
    p=ESTIMATE.copy(); p[8]=phi; p[10]=np.log(lam); p[:3]+=[dg,df,da]; return p
for tag,args in [('uf+.01',(1.2,1,0,0.01,0)),('uf+.02',(1.2,1,0,0.02,0)),('g-.0005',(1.2,1,-5e-4,0,0)),('ua+.03',(1.2,1,0,0,0.03)),
                 ('l.5p5',(5,0.5)),('l.5p3',(3,0.5)),('l.7p4',(4,0.7)),('l1p3.5',(3.5,1)),('l1.4p1.5',(1.5,1.4)),('l1.4p2',(2,1.4))]:
    pts.append(P(*args)); tags.append(tag)
for t,r in zip(tags,ev.Batch(pts)): print(t, f"q-qbar {r['criterion']-BUDGET:+.6f} c {r['target']:.3f}")
print(len(ev))
