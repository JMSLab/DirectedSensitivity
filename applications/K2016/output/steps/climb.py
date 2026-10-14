import sys
sys.path.insert(0,'.')
from basis import Evaluator
from search import Climb
sign=int(sys.argv[1]); rounds=int(sys.argv[2])
Climb(Evaluator(), sign, rounds, seed=100+sign)
