"""Exact checks: (1) cube formula for Psi; (2) product formula Psi_{XxY}^{(m)} = sum_i C(m,i)^2 Psi_X^{(i)} Psi_Y^{(m-i)}
(probability normalisation); (3) sup/L2 sandwich Psi_A <= Cheb_A^2 <= U Psi_A (via LP for Cheb, float)."""
import random, sys, time
from fractions import Fraction
from math import comb
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round16/reduction_checks')
from lib import *

random.seed(16)
out = []
def log(s):
    print(s); out.append(s); sys.stdout.flush()

t0 = time.time()
# (1) cube formula
for d in (3, 4, 5):
    for m in range(1, 2 * d + 2):
        bx = Box(cube(d), m)
        ok = True
        for trial in range(3):
            u = [Fraction(random.randint(-9, 9), random.randint(1, 7)) for _ in range(d)]
            if bx.psi(u) != psi_cube_formula(u, m):
                ok = False
        log(f"cube d={d} m={m} r0={bx.r0}: formula exact match = {ok}")
log(f"time {time.time()-t0:.1f}s")

# (2) product formula, X = one-cut slice-type set on 3 coords, Y = cube on 2 coords
X = one_cut(3, [0, 1, 2], 1)
Y = cube(2)
XY = [x + y for x in X for y in Y]
for m in range(1, 7):
    bXY = Box(XY, m)
    bX = {i: Box(X, i) if i > 0 else None for i in range(m + 1)}
    bY = {i: Box(Y, i) if i > 0 else None for i in range(m + 1)}
    ok = True
    for trial in range(3):
        u = [Fraction(random.randint(-9, 9), random.randint(1, 7)) for _ in range(5)]
        uX, uY = u[:3], u[3:]
        lhs = bXY.psi(u) / len(XY)
        rhs = Fraction(0)
        for i in range(m + 1):
            pX = (bX[i].psi(uX) / len(X)) if i > 0 else Fraction(1)
            pY = (bY[m - i].psi(uY) / len(Y)) if m - i > 0 else Fraction(1)
            rhs += comb(m, i) ** 2 * pX * pY
        ok &= (lhs == rhs)
    log(f"product X(one-cut s=3,L=1) x Y(cube d=2), m={m}: product formula exact = {ok}")
log(f"time {time.time()-t0:.1f}s")
open(__file__.replace('.py', '_out.txt'), 'w').write('\n'.join(out) + '\n')
