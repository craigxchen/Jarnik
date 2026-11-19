"""Referee checks of (II.5), (II.6), (II.7), (II.9) on random Gaussian blocks, and of the
line-graph matrix facts used in II.5.  Exact integer arithmetic throughout."""
import random, math
from itertools import combinations
from fractions import Fraction
random.seed(31)

def gmul(a, b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def gconj(a): return (a[0], -a[1])
def gnorm(a): return a[0]*a[0]+a[1]*a[1]
def br(A, B):   # [A,B] = Im(conj(A) B)
    return A[0]*B[1]-A[1]*B[0]
def re_(A, B): return A[0]*B[0]+A[1]*B[1]

E = list(combinations(range(5), 2))
def edge(i, j): return (min(i, j), max(i, j))

ntest = 0; ok7 = 0
for trial in range(3000):
    # random blocks (not necessarily prime, not necessarily coprime: identities are algebraic)
    beta = {e: (random.randint(-10**6, 10**6), random.randint(1, 10**6)) for e in E}
    n = {e: gnorm(beta[e]) for e in E}
    P = []
    for i in range(5):
        z = (1, 0)
        for j in range(5):
            if j != i: z = gmul(z, beta[edge(i, j)])
        P.append(z)
    t = {}; s = {}
    for (i, j) in E:
        m = br(P[i], P[j]); assert m % n[(i, j)] == 0
        r = re_(P[i], P[j]); assert r % n[(i, j)] == 0
        t[(i, j)] = m // n[(i, j)]; s[(i, j)] = r // n[(i, j)]
    def T(a, b): return t[(a, b)] if a < b else -t[(b, a)]
    def S(a, b): return s[edge(a, b)]
    def N(a, b): return n[edge(a, b)]
    for (i, j) in E:
        D = 1
        for k in range(5):
            if k not in (i, j): D *= N(i, k)*N(j, k)
        assert s[(i, j)]**2 + t[(i, j)]**2 == D               # (II.5)
        # |u_i-u_j|^2 = 4 t^2 / D, u = P/conj P : check with exact rationals
        # u_i - u_j = (P_i conj P_j - P_j conj P_i)/(conj P_i conj P_j) ; |.|^2 = (2 Im(conj P_j P_i))^2 / (|P_i|^2 |P_j|^2)
        assert (2*br(P[j], P[i]))**2 * D == 4*t[(i, j)]**2 * gnorm(P[i])*gnorm(P[j])
    for (i, j, k) in combinations(range(5), 3):
        for (a, b, c) in [(i, j, k), (j, k, i), (k, i, j), (j, i, k)]:
            l, mm = [h for h in range(5) if h not in (a, b, c)]
            A = T(b, c)**2*N(b, c)*N(a, l)*N(a, mm) + T(a, c)**2*N(a, c)*N(b, l)*N(b, mm) - T(a, b)**2*N(a, b)*N(c, l)*N(c, mm)
            assert 2*T(a, c)*T(b, c)*S(a, b) == A                       # (II.6), signed version
            Dab = 1
            for h in range(5):
                if h not in (a, b): Dab *= N(a, h)*N(b, h)
            assert 4*T(a, c)**2*T(b, c)**2*Dab - A**2 == 4*(T(a, b)*T(a, c)*T(b, c))**2   # (II.9)
            # (II.7): for a divisor b0 of D_ab and signed residues r_e mod b0
            for b0 in [N(a, h) for h in range(5) if h not in (a, b)][:2] + [math.gcd(Dab, 10**9+7)*N(a, c)]:
                if Dab % b0: continue
                def rr(x): y = x % b0; return y - b0 if y > b0//2 else y
                Ar = T(b, c)**2*rr(N(b, c))*rr(N(a, l))*rr(N(a, mm)) + T(a, c)**2*rr(N(a, c))*rr(N(b, l))*rr(N(b, mm)) - T(a, b)**2*rr(N(a, b))*rr(N(c, l))*rr(N(c, mm))
                assert (Ar**2 + 4*(T(a, b)*T(a, c)*T(b, c))**2) % b0 == 0
                ok7 += 1
    # (8): sum_i P_i / prod_{j != i} t_ij = 0 (signed t_ij = [P_i,P_j]/n_ij, oriented i->j)
    tot = [Fraction(0), Fraction(0)]
    for i in range(5):
        den = 1
        for j in range(5):
            if j != i: den *= T(i, j)
        tot[0] += Fraction(P[i][0], den); tot[1] += Fraction(P[i][1], den)
    assert tot == [0, 0]
    ntest += 1
print('random block tuples checked:', ntest, ' (II.7) divisibility instances:', ok7)

# Line graph of K5: M (adjacency), 6 M^{-1} = 3I + 3M - 2J
import sympy as sp
M = sp.zeros(10, 10)
for a in range(10):
    for b in range(10):
        if a != b and len(set(E[a]) & set(E[b])) == 1: M[a, b] = 1
print('det M =', M.det(), ' eigen', M.eigenvals())
Inv6 = 6*M.inv()
print('6M^{-1} == 3I+3M-2J :', Inv6 == 3*sp.eye(10)+3*M-2*sp.ones(10, 10))
