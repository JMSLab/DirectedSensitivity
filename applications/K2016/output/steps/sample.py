import sys
sys.path.insert(0,'.')
from basis import Evaluator
from search import Sample
ev=Evaluator()
Sample(ev, 1, 30, [1.2, 2.0, 2.5, 3.0], [0.7, 1.0], (0.0, 1.5e-4), 201)
Sample(ev, -1, 30, [0.3, 0.6], [1.0, 2.0], (1.5e-4, 3.5e-4), 202)
print(len(ev))
