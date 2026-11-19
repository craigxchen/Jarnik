"""Referee checks of the six-norm short-interval lemma (II.4).

(1) symbolic: Q is the polynomial square root of F to order Y^3, 16Q integral, F-Q^2 has degree <= 2.
(2) exhaustive search: for X in a range and all multisets {k_1..k_6} in [0,K] with min 0,
    find the smallest |t|>0 with F(X) = s^2 + t^2 (exact), and report min |t|/sqrt(X).
    K is taken as floor(C * X^theta) for several (C, theta), including windows far wider
    than the lemma allows, to see where small-|t| representations appear.
"""
import sys, math
from itertools import combinations_with_replacement
import sympy as sp

Y = sp.symbols('Y'); ks = sp.symbols('k1:7')
F = sp.expand(sp.prod([Y+k for k in ks]))
e = [sp.Integer(1)] + [sp.expand(sum(sp.prod(c) for c in __import__('itertools').combinations(ks, j))) for j in range(1, 7)]
Q = Y**3 + e[1]/2*Y**2 + (e[2]/2 - e[1]**2/8)*Y + e[3]/2 - e[1]*e[2]/4 + e[1]**3/16
R = sp.expand(F - Q**2)
print('deg_Y(F-Q^2) =', sp.degree(R, Y))
# homogeneity degrees of coefficients
for j in range(3):
    c = sp.Poly(R, Y).coeff_monomial(Y**j)
    print(' coeff Y^%d total degree in k:' % j, sp.Poly(c, *ks).total_degree())
print(' 16Q in Z[k][Y]:', all(sp.Poly(16*Q, Y, *ks).coeffs()[i].is_integer for i in range(len(sp.Poly(16*Q, Y, *ks).coeffs()))))

def min_t(N):
    """smallest t>0 with N = s^2+t^2, searching s downward from isqrt(N) while t stays below bound"""
    s0 = math.isqrt(N)
    best = None
    for s in range(s0, max(s0-3, -1), -1):
        r = N - s*s
        if r == 0: continue
        t = math.isqrt(r)
        if t*t == r:
            best = t if best is None else min(best, t)
    return best   # s0-3 already gives t ~ sqrt(6 sqrt N) >> sqrt X, so this is exact for t < sqrt(2 sqrtN)

def scan(Xlo, Xhi, Cc, theta, step=1):
    worst = (float('inf'), None)
    for X in range(Xlo, Xhi, step):
        K = max(1, int(Cc * X**theta))
        for ksub in combinations_with_replacement(range(0, K+1), 5):
            kk = (0,) + ksub
            N = 1
            for k in kk: N *= (X+k)
            t = min_t(N)
            if t is None: continue
            r = t / math.sqrt(X)
            if r < worst[0]: worst = (r, (X, kk, t))
    return worst

if __name__ == '__main__':
    for (Xlo, Xhi, Cc, th, st) in [(10**3, 10**3+2000, 1.0, 1/6, 1), (10**5, 10**5+300, 1.0, 1/6, 1),
                                   (10**6, 10**6+40, 1.0, 1/6, 1), (10**3, 10**3+300, 1.0, 1/4, 1),
                                   (10**4, 10**4+60, 1.0, 1/4, 1), (10**3, 10**3+100, 1.0, 1/3, 1),
                                   (200, 500, 1.0, 1/2, 25)]:
        w = scan(Xlo, Xhi, Cc, th, st)
        print('X in [%d,%d) K=%.1f X^%.3f : min |t|/sqrt X = %.4g at %s' % (Xlo, Xhi, Cc, th, w[0], w[1]), flush=True)
