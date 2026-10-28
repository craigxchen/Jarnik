"""Random sample: pairs (A,C) at D=5,6 (uniform over orbit reps x C), random lambda in Q^* (small height),
and lambda=+-1: literal formula vs prediction (-N0 for non-candidate lambda). Tests Prop 2.1/Thm 2.3/(1.1) empirically."""
import sys, json, random, itertools, math
from fractions import Fraction
sys.path.insert(0, '.')
from rcore import *
rng = random.Random(int(sys.argv[2]))
D = int(sys.argv[1]); NS = int(sys.argv[3])
reps = json.load(open('rreps_D%d.json' % D))
bad = []; done = 0; hist = {}
while done < NS:
    A = [tuple(e) for e in rng.choice(reps)]
    rest = [e for e in EDGES if e not in set(A)]
    C = sorted(rng.choice(rest) for _ in range(D))
    n0 = N0(vec(A), vec(C))
    p = rng.choice([-1, 1]) * rng.randint(1, 50); q = rng.randint(1, 50)
    if math.gcd(p, q) != 1 or abs(Fraction(p, q)) == 1:
        continue
    g = full_gain(A, C, p, q)
    hist[g + n0] = hist.get(g + n0, 0) + 1
    if g != -n0:
        bad.append([A, C, p, q, g, -n0])
    done += 1
print(json.dumps(dict(D=D, samples=NS, excess_hist=hist, deviations=bad)))
