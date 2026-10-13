# Sanity checks for elamp.degE: (i) V_A member -> 0 ; (ii) ruling conic of a random quadric -> <= 6.
import random, sympy as sp
from sympy import Rational as R, Matrix
from fractions import Fraction as Fr
from elamp import degE, md
from star import member_through, point_of, quadric_through, t
from arr import PLANES, CHOSEN
pts = [point_of([PLANES[i], PLANES[j], PLANES[k]]) for (i, j, k) in CHOSEN]
random.seed(4)
def rnd(): return R(random.randint(-30, 30), random.randint(1, 9))
def to_mob(data):
    out = []
    for (N, D) in data:
        Np = sp.Poly(sp.expand(N), t); Dp = sp.Poly(sp.expand(D), t)
        al, be = Np.coeff_monomial(t), Np.coeff_monomial(1); ga, de = Dp.coeff_monomial(t), Dp.coeff_monomial(1)
        out.append([md(Fr(int(sp.numer(c)), int(sp.denom(c)))) for c in (al, be, ga, de)])
    return out
x = [rnd() for _ in range(7)]
data, Q, _ = member_through(PLANES, pts, x)
print('V_A member: E_Lambda degree =', degE(to_mob(data))[0])
# ruling conic of a random quadric containing the line R_x but NOT through the 6 points: use 5 random points instead
rp = [Matrix([rnd(), rnd(), rnd(), 1]) for _ in range(6)]
data2, Q2, _ = member_through(PLANES, rp, x)
print('ruling conic of a random quadric: E_Lambda degree =', degE(to_mob(data2))[0])
