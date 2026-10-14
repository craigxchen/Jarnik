"""(Mar) checks for sdc_q and pal_q pattern data: (AI), G_T = L_T/Z^n_0 (Smith form), relevance of
fractional characters for k = e_1 (sdc), and the margin sigma_c >= rho_1 ||c||_1 on integral c with small ||c||_1."""
import itertools, sympy, numpy as np
from fractions import Fraction
from rlib import *
from sympy.matrices.normalforms import smith_normal_form

def omit0(n, S):
    S = set(S); return S if 0 not in S else set(range(n)) - S

def pal_mult(q, b=8, astar=1):
    H = paley_I(q); M = q+1; Hn = H.copy(); Hn[1:, :] *= -1
    cols = []
    for a in range(1, M):
        for _ in range(b): cols.append(Hn[:, a].copy())
    for x in range(M):
        lab = astar if x == 0 else x
        c = Hn[:, lab].copy(); c[x] *= -1; cols.append(c)
    mult = {}
    for c in cols:
        S = canon(M, {y for y in range(M) if c[y] == 1}); mult[S] = mult.get(S, 0) + 1
    return mult, len(cols)

def analyse(name, n, mult, kvec, maxl1):
    pats = list(mult)
    P = [omit0(n, S) for S in pats]
    N = sympy.Matrix([[1 if x in Pi else 0 for x in range(1, n)] for Pi in P])
    rk = N.rank()
    snf = smith_normal_form(N, domain=sympy.ZZ)
    divs = [snf[i, i] for i in range(min(snf.shape)) if snf[i, i] != 0]
    print(name, 'distinct patterns', len(pats), 'r =', sum(mult.values()), '(AI) rank', rk, '= d?', rk == n-1,
          'elementary divisors', [abs(x) for x in divs if abs(x) != 1], flush=True)
    # fractional characters: lambda' = N^+ z; enumerate G_T via lattice {lambda' : N lambda' in Z}
    # brute: G_T is finite; find generators by solving N lambda' = e_j-combos is heavy; instead test candidate
    # lambda* = (2/n) 1 - 2 e_0 (exact.md Thm 4.1(ii)) and its multiples, and check |G_T| = prod divisors.
    order = 1
    for x in divs: order *= abs(x)
    print('  |G_T| =', order)
    if order > 1:
        lam = [Fraction(2, n) - (2 if x == 0 else 0) for x in range(n)]
        for j in range(1, order):
            lj = [j*v for v in lam]
            inLT = all((sum(lj[x] for x in S)).denominator == 1 for S in pats) and sum(lj) == 0
            lk = sum(lj[x]*kvec[x] for x in range(n))
            print('   j=%d: j*lambda* in L_T: %s, lambda.k = %s (relevant: %s)' % (j, inLT, lk, lk.denominator == 1))
    # margin on integral zero-sum c with ||c||_1 <= maxl1 (all of them)
    worst = {}
    M = np.array([[1 if x in S else 0 for x in range(n)] for S in pats], dtype=np.int64)
    w = np.array([mult[S] for S in pats], dtype=np.int64)
    for l1 in range(2, maxl1+1, 2):
        best = None
        # enumerate c with sum 0, ||c||_1 = l1: choose positive part and negative part
        for cpos in compositions(n, l1//2):
            for cneg in compositions(n, l1//2):
                if any(a and b for a, b in zip(cpos, cneg)): continue
                c = np.array(cpos) - np.array(cneg)
                if c[0] == 0 and False: pass
                sig2 = (w * (2*np.abs(M @ c) - 1)).sum()     # = 2 sigma
                ratio = Fraction(int(sig2), 2*l1)
                if best is None or ratio < best[0]: best = (ratio, c.tolist())
        worst[l1] = best
        print('  ||c||_1 = %d: min sigma_c/||c||_1 = %s at c = %s' % (l1, best[0], best[1]), flush=True)
    return worst

def compositions(n, s):
    # all nonnegative integer vectors of length n summing to s
    if n == 1:
        yield (s,); return
    for a in range(s+1):
        for rest in compositions(n-1, s-a): yield (a,) + rest

for q in (7, 11):
    n = q+1
    mult = {S: 1 for S in sdc(q)}
    k = [0]*n; k[1] = 1
    analyse('sdc%d' % q, n, mult, k, 6 if q == 7 else 4)
    mult, r = pal_mult(q)
    analyse('pal%d (b=8)' % q, n, mult, [0]*n, 6 if q == 7 else 4)
