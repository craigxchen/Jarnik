"""Referee (analytic lens): Appendix A (Sylvester multi-block route) of sharp.md.

(1) Lemma A.2 level-set bound: max_(v != 0) prod_(u in [-3,3], u != v) (3 + |u|)/|v - u| = 225, and the
    coefficient count sum_(v != 0) |v| = 12 (so L = 12 L_GS(225)).
(2) Coset normalisation: ||1_(x+H)||_A = 1 for every coset of every subgroup H of F_2^m (m = 4, exhaustive
    over a sample of subgroups), with fhat(a) = E_x f(x)(-1)^(a.x).
(3) 'Support 4 gives exactly 1/(M-3)': exact count over F_2^m, m = 3..6, of ordered 4-tuples of distinct
    points lying in a 2-dimensional affine subspace (x1+x2+x3+x4 = 0).
(4) The Assembly union bound at the stated threshold s = 4 L log M, t = 5 blocks, entries <= 3:
    naive per-character union   N_s (|S_L|/Mult_min)^t        with N_s = C(M,s) 6^s, Mult_min = C(M, ceil(s/4))
    multiset-class union        #multisets |S_L|^t Mult_min^(1-t)
    log|S_L| <= L log2-count of (sign x affine subspaces) = L log(2 * 3.47 (m+1) 2^(m^2/4+m)).
    Reports natural logs; negative = the union bound closes.
"""
import itertools, math
import numpy as np

vals = range(-3, 4)
best = 0
for v in vals:
    if v == 0:
        continue
    p = 1.0
    for u in vals:
        if u != v:
            p *= (3 + abs(u)) / abs(v - u)
    best = max(best, p)
print(f"(1) max_(v != 0) level-set bound = {best:.1f};  sum_(v != 0) |v| = {sum(abs(v) for v in vals)}")

m = 4
M = 1 << m
pts = list(range(M))
def dot(a, x):
    return bin(a & x).count("1") & 1
ok = True
# subgroups spanned by random small sets
rng = np.random.default_rng(1)
for trial in range(200):
    gens = rng.integers(0, M, size=rng.integers(0, 4))
    H = {0}
    for g in gens:
        H |= {h ^ int(g) for h in H}
    x = int(rng.integers(0, M))
    f = np.array([1.0 if (p ^ x) in H else 0.0 for p in pts])
    fh = np.array([sum(f[p] * (-1) ** dot(a, p) for p in pts) / M for a in pts])
    ok &= abs(np.abs(fh).sum() - 1.0) < 1e-12
print(f"(2) ||1_(x+H)||_A = 1 for 200 random cosets in F_2^4: {ok}")

for m in range(3, 7):
    M = 1 << m
    cnt = tot = 0
    for x1, x2, x3 in itertools.permutations(range(M), 3):
        x4 = x1 ^ x2 ^ x3
        tot += M - 3
        if x4 not in (x1, x2, x3):
            cnt += 1
    print(f"(3) m = {m}: Pr[4 random distinct points coplanar] = {cnt}/{tot} = {cnt / tot:.6f};  1/(M-3) = {1 / (M - 3):.6f}")

def lbinom(logn, k):
    """log C(n, k) for n = e^logn huge, k << n (stable: k log n - log k! + sum log(1 - i/n))."""
    n = math.exp(logn) if logn < 700 else float('inf')
    corr = 0.0 if n == float('inf') else sum(math.log1p(-i / n) for i in range(int(k)))
    return k * logn - math.lgamma(k + 1) + corr
t = 5
print("(4) Assembly union bound at s = 4 L log M (natural logs; negative = closes):")
for m in (40, 100, 400, 2000):
    lM = m * math.log(2)
    for Lgs in (10, 1000):
        s = int(4 * Lgs * lM)
        if math.log(2 * s) > lM:
            continue
        logSL = Lgs * math.log(2 * 3.47 * (m + 1)) + Lgs * (m * m / 4 + m) * math.log(2)
        logMult = lbinom(lM, math.ceil(s / 4))
        logN = lbinom(lM, s) + s * math.log(6)
        naive = logN + t * (logSL - logMult)
        multiset = 7 * (lM + 1e-9) + t * logSL + (1 - t) * logMult
        print(f"   m = {m:5d}, L = {Lgs:5d}: naive = {naive:14.4g};  multiset = {multiset:14.4g}")
