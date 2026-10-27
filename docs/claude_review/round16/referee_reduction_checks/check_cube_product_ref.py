# Independent exact check of Thm 5.1(a) (cube Psi formula) and Thm 5.2 (product formula), with a NON-slab factor Y.
import itertools, random, sys
from fractions import Fraction
from math import factorial
import flint
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round16/referee_reduction_checks')
from rlib import psi_sigma, poly_basis

random.seed(7)
def rnd():
    return Fraction(random.randint(-9, 9), random.randint(1, 5))

def cube(d):
    return list(itertools.product((-1, 0, 1), repeat=d))

def formula_cube(d, m, u):
    tot = Fraction(0)
    for a in itertools.product(range(3), repeat=d):
        if sum(a) != m: continue
        mult = factorial(m)
        for ai in a: mult //= factorial(ai)
        n1 = sum(1 for ai in a if ai == 1); n2 = sum(1 for ai in a if ai == 2)
        term = Fraction(mult) ** 2 * Fraction(2, 3) ** n1 * Fraction(2, 9) ** n2
        for ui, ai in zip(u, a): term *= ui ** (2 * ai)
        tot += term
    return tot

def fq(x):
    return Fraction(int(x.p), int(x.q))

bad = 0; cnt = 0
for d in (2, 3, 4):
    A = cube(d)
    for m in range(1, 2 * d + 2):
        B = poly_basis(A, m)
        for t in range(3):
            u = [rnd() for _ in range(d)]
            ps = fq(psi_sigma(A, m, u, B)) / len(A)
            fo = formula_cube(d, m, u) if m <= 2 * d else Fraction(0)
            cnt += 1
            if ps != fo:
                bad += 1; print('CUBE MISMATCH', d, m, u, ps, fo)
print('cube formula: %d cases, %d mismatches' % (cnt, bad))

# product formula with X = one-cut block (s=3, L=1) and Y = an irregular subset of {-1,0,1}^2 (not a slab box)
X = [k for k in cube(3) if abs(sum(k)) <= 1]
Y = [k for k in cube(2) if k not in [(1, 1), (-1, 0), (0, -1)]]
XY = [x + y for x in X for y in Y]
def psi_prob(S, m, u):
    if m == 0: return Fraction(1)
    return fq(psi_sigma(S, m, u)) / len(S)
bad = 0; cnt = 0
for m in range(1, 8):
    for t in range(3):
        u = [rnd() for _ in range(5)]
        lhs = psi_prob(XY, m, u)
        rhs = sum(Fraction(factorial(m) // (factorial(i) * factorial(m - i))) ** 2 * psi_prob(X, i, u[:3]) * psi_prob(Y, m - i, u[3:]) for i in range(m + 1))
        cnt += 1
        if lhs != rhs:
            bad += 1; print('PRODUCT MISMATCH', m, u, lhs, rhs)
print('product formula (X one-cut s=3 L=1, Y irregular 6-point set): %d cases, %d mismatches' % (cnt, bad))
