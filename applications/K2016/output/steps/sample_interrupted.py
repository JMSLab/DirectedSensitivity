import sys
sys.path.insert(0,'.')
from basis import Evaluator
from search import Sample
ev=Evaluator()
Sample(ev, 1, 24, [3.0, 4.0, 5.0, 7.0], [0.3, 0.5], (0.0, 1.0e-4), 301)
Sample(ev, -1, 20, [0.3, 0.6, 1.0], [0.5, 1.0], (2.0e-4, 4.0e-4), 302)
print(len(ev))
