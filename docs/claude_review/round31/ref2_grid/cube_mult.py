"""Test: average multiplicity over {0,1}^m of a nonzero polynomial of degree V is <= V/2
(multiplicity Schwartz-Zippel, DKSS, grid S={0,1}). Random, structured and hill-climbed attempts
at a counterexample; plus the cap E(m_S+m_{S^c}) <= V/2 on the actual grid determinants."""
import random
from fractions import Fraction
from itertools import product
from math import comb

# polynomials as dict: exponent tuple -> int coefficient
def mul(P, Q):
    R = {}
    for a, c in P.items():
        for b, e in Q.items():
            k = tuple(x + y for x, y in zip(a, b))
            R[k] = R.get(k, 0) + c * e
    return {k: v for k, v in R.items() if v}

def add(P, Q, s=1):
    R = dict(P)
    for k, v in Q.items():
        R[k] = R.get(k, 0) + s * v
    return {k: v for k, v in R.items() if v}

def shift_min_degree(P, pt):
    # expand P(pt + y): each x_i^a with pt_i=1 -> sum_b C(a,b) y_i^b
    R = {}
    for a, c in P.items():
        parts = [[(b, comb(ai, b)) for b in range(ai + 1)] if pt[i] else [(ai, 1)] for i, ai in enumerate(a)]
        for choice in product(*parts):
            k = tuple(b for b, _ in choice)
            w = c
            for _, cc in choice:
                w *= cc
            R[k] = R.get(k, 0) + w
    R = {k: v for k, v in R.items() if v}
    return min(sum(k) for k in R)

def degree(P):
    return max(sum(k) for k in P)

def avg_mult(P, m):
    return Fraction(sum(shift_min_degree(P, pt) for pt in product((0, 1), repeat=m)), 2 ** m)

def lin(m, coeffs, const):
    P = {tuple(int(i == j) for j in range(m)): c for i, c in enumerate(coeffs) if c}
    if const:
        P[(0,) * m] = P.get((0,) * m, 0) + const
    return P

random.seed(7)
worst = Fraction(0)
count = 0
# 1. random dense polynomials
for t in range(300):
    m = random.randint(2, 4); V = random.randint(1, 5)
    P = {}
    for k in product(range(V + 1), repeat=m):
        if sum(k) <= V and random.random() < 0.5:
            P[k] = random.randint(-3, 3)
    P = {k: v for k, v in P.items() if v}
    if not P or degree(P) == 0: continue
    r = avg_mult(P, m) / degree(P); worst = max(worst, r); count += 1
print('random dense: tested', count, 'max avg/deg =', worst)

# 2. random sums of products of affine forms vanishing on many cube points
worst2 = Fraction(0)
for t in range(400):
    m = random.randint(3, 5)
    terms = []
    D = random.randint(2, 5)
    P = {}
    for s in range(random.randint(1, 3)):
        Q = {(0,) * m: 1}
        for f in range(D):
            co = [random.choice((-1, 0, 0, 1)) for _ in range(m)]
            if not any(co): co[0] = 1
            Q = mul(Q, lin(m, co, random.choice((-1, 0, 0, 1, 2))))
        P = add(P, Q, random.choice((-1, 1)))
    if not P or degree(P) == 0: continue
    r = avg_mult(P, m) / degree(P); worst2 = max(worst2, r)
print('sums of products of affine forms: max avg/deg =', worst2)

# 3. hill climbing: perturb coefficients of a polynomial to maximize avg mult / degree
def score(P, m):
    return avg_mult(P, m) / degree(P) if P and degree(P) > 0 else Fraction(0)
best = Fraction(0)
for run in range(25):
    m = 3; V = 4
    mons = [k for k in product(range(V + 1), repeat=m) if sum(k) <= V]
    P = {k: random.randint(-2, 2) for k in mons}
    P = {k: v for k, v in P.items() if v} or {mons[-1]: 1}
    s = score(P, m)
    for it in range(400):
        Q = dict(P)
        k = random.choice(mons)
        Q[k] = Q.get(k, 0) + random.choice((-1, 1))
        Q = {a: b for a, b in Q.items() if b}
        if not Q: continue
        sq = score(Q, m)
        if sq >= s:
            P, s = Q, sq
    best = max(best, s)
print('hill climbing m=3, deg<=4: best avg/deg =', best)

# 4. the actual grid determinants: cap E(m_S + m_{S^c}) <= V/2 and avg multiplicity
import sympy as sp
z = sp.symbols('z0:6')
Mc = sp.Matrix([[1] + [z[(i + h) % 6] for h in range(5)] for i in range(6)])
Fc = {k: int(v) for k, v in sp.Poly(sp.expand(Mc.det()), *z).terms()}
def cap(P, m):
    tot = 0
    for S in product((0, 1), repeat=m):
        mS = min(sum(k[i] for i in range(m) if S[i]) for k in P)
        mC = min(sum(k[i] for i in range(m) if not S[i]) for k in P)
        tot += mS + mC
    return Fraction(tot, 2 ** m)
print('cyclic det (V=5): E(m_S+m_Sc) =', cap(Fc, 6), ' avg mult =', avg_mult(Fc, 6))
# d=1 modified design, generic nodes, 7 distinct marks used (s1_0,s2_0,s3_0,s3_1,s4_0,s4_1,s5_0,s5_1 + u,v)
s = sp.symbols('w0:12')
F1 = sp.expand((s[0]-s[1])*(s[2]-s[3])*(s[4]-s[5])*(s[6]-s[7])*(s[6]+s[7]-s[8]-s[9]))
F1d = {k: int(v) for k, v in sp.Poly(F1, *s[:10]).terms()}
print('d=1 modified design (V=5): E(m_S+m_Sc) =', cap(F1d, 10), ' avg mult =', avg_mult(F1d, 10))
