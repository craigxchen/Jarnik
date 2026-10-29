import sys, itertools
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round27/interp_checks')
from shiftlib import *
d = 2
A = [k for k in itertools.product([-1,0,1], repeat=2) if abs(k[0]+k[1]) <= 1]
L = [l for l in itertools.product(range(-2,3), repeat=2) if l != (0,0)]
for e in range(0, 5):
    print('hexagon e', e, 'dim M_e(A), dim sum_l M_e(A_l):', shift_test(A, e, d, L), 'HF', HF(A, e, d))
cube = list(itertools.product([-1,0,1], repeat=3))
L3 = [l for l in itertools.product(range(-2,3), repeat=3) if any(l)]
for e in range(0, 7):
    print('cube3 e', e, shift_test(cube, e, 3, L3))
