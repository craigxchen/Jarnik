# Random/exhaustive checks of the algebraic steps of III.7, III.8, III.10.
import random, itertools, math
from math import gcd
from fractions import Fraction
import sympy
random.seed(31)

primes = [5, 13, 17, 29, 37]

def gprod(a, i, j):
    return math.prod(p**min(a[i][k], a[j][k]) for k, p in enumerate(primes))

def reduce_frac(n, d):
    g = gcd(n, d); return n//g, d//g

# --- III.8: telescoping of crossed g-quotients and unit/valuation claims --------------
bad = 0
for trial in range(20000):
    e = [random.randint(1, 4) for _ in primes]
    a = [[random.randint(0, e[k]) for k in range(len(primes))] for _ in range(5)]
    i, j, k, l, h = range(5)
    def cd(k_, l_):
        return reduce_frac(gprod(a, i, k_)*gprod(a, j, l_), gprod(a, j, k_)*gprod(a, i, l_))
    C1, D1 = cd(k, l); C2, D2 = cd(l, h); C3, D3 = cd(k, h)
    if C1*C2*D3 != D1*D2*C3: bad += 1
    # Q = product of layers separating {i,j} from {k,l,h}
    Q = 1
    for kk, p in enumerate(primes):
        for hh in range(1, e[kk]+1):
            s = [a[x][kk] >= hh for x in range(5)]
            if s[i] == s[j] and s[k] == s[l] == s[h] and s[i] != s[k]: Q *= p
    for (kk_, ll_), (Cs, Ds) in [((k, l), (C1, D1)), ((l, h), (C2, D2)), ((k, h), (C3, D3))]:
        if gcd(Cs, Q) != 1 or gcd(Ds, Q) != 1: bad += 1
        # uncrossed reduced product divisible by Q
        g1 = gprod(a, i, j)*gprod(a, kk_, ll_); g2 = gprod(a, i, kk_)*gprod(a, j, ll_); g3 = gprod(a, i, ll_)*gprod(a, j, kk_)
        c = gcd(gcd(g1, g2), g3)
        if (g1//c) % Q: bad += 1
print("III.8 telescoping / unit / divisibility failures in 20000 random prime-power arrays:", bad)

# --- pattern-indicator expansions used in III.7 (4 marks) and III.8 (5 marks) -----------
for n, coef in [(4, Fraction(1, 8)), (5, Fraction(1, 16))]:
    ok = True
    for eps in itertools.product([1, -1], repeat=n):
        for f in itertools.product([1, -1], repeat=n):
            ind = 1 if (f == eps or f == tuple(-x for x in eps)) else 0
            s = 1
            for r in range(2, n+1, 2):
                for S in itertools.combinations(range(n), r):
                    s += math.prod(eps[x]*f[x] for x in S)
            if n == 4: val = coef*s
            else: val = coef*s
            if val != ind: ok = False
    print(f"indicator expansion on {n} marks (even-degree terms, coefficient {coef}):", ok)

# --- frame bound: Gram of sign rows with negative off-diagonal has top eigenvalue <= 2 --
import numpy as np
worst = 0
for trial in range(300):
    m = random.randint(5, 40); L = random.randint(m, 3*m)
    w = np.random.rand(L); w /= w.sum()
    F = np.where(np.random.rand(m, L) < 0.5, 1.0, -1.0)
    Gm = (F*w) @ F.T
    if ((Gm - np.diag(np.diag(Gm)))[~np.eye(m, dtype=bool)] >= 0).any():
        continue  # hypothesis (pair separation) fails; skip
    worst = max(worst, np.linalg.eigvalsh(Gm).max())
print("frame bound: worst top eigenvalue among random negative-Gram sign matrices:", worst)

# --- III.10: integer realization and primitive-column criterion ------------------------
def minors(V):
    m = len(V[0]); return {(i, j): V[0][i]*V[1][j]-V[1][i]*V[0][j] for i in range(m) for j in range(i+1, m)}
fails = 0; tested = 0
for trial in range(300):
    m = 5
    V = [[random.randint(-6, 6) for _ in range(m)] for _ in range(2)]
    D = minors(V); d = 0
    for v in D.values(): d = gcd(d, v)
    if d == 0 or abs(d) > 10 or abs(d) < 2: continue
    # saturation of row space: basis of (row space) cap Z^m via nullspace of nullspace
    Mv = sympy.Matrix(V)
    K = Mv.nullspace()
    K = sympy.Matrix.hstack(*[k*sympy.ilcm(*[x.q for x in k]) for k in K]).T
    Bs = K.nullspace()
    Bs = sympy.Matrix.hstack(*[b*sympy.ilcm(*[x.q for x in b]) for b in Bs]).T
    # LLL-free saturation: HNF of the integer row lattice then saturate via Smith form
    from sympy.matrices.normalforms import smith_normal_form
    # compute saturation: lattice = Q-span cap Z^m; use: B = rows of Bs; saturate by dividing out
    Bint = [[int(x) for x in Bs.row(r)] for r in range(2)]
    # brute saturation: minors gcd of Bint; divide by searching small rational combinations
    DB = minors(Bint); dB = 0
    for v in DB.values(): dB = gcd(dB, v)
    # saturated basis has minors gcd 1; find it by searching (Bint rows)/k
    import itertools as it
    def saturate(Bm):
        while True:
            Dm = minors(Bm); g = 0
            for v in Dm.values(): g = gcd(g, v)
            if g == 1: return Bm
            p = sympy.factorint(g); p = min(p)
            done = False
            for c0, c1 in [(1, 0), (0, 1)] + [(1, c) for c in range(p)]:
                row = [c0*Bm[0][t]+c1*Bm[1][t] for t in range(m)]
                if all(x % p == 0 for x in row):
                    row = [x//p for x in row]
                    Bm = [row, Bm[1]] if c0 else [Bm[0], row]
                    done = True; break
            assert done
    Bsat = saturate(Bint)
    DS = minors(Bsat)
    ratio = set(Fraction(D[k], DS[k]) for k in D if DS[k] != 0)
    zero_ok = all((D[k] == 0) == (DS[k] == 0) for k in D)
    tested += 1
    if not (zero_ok and len(ratio) == 1 and abs(list(ratio)[0]) == d): fails += 1
    # primitive criterion: V has all columns primitive iff (Bsat columns primitive and for p|d
    # reductions miss a point of P^1(F_p)); test on V itself by recovering its multiplier
    Bsat_cols = [(Bsat[0][t], Bsat[1][t]) for t in range(m)]
    def crit():
        if any(gcd(*c) != 1 for c in Bsat_cols): return False
        for p in sympy.factorint(abs(d)):
            pts = set()
            for c in Bsat_cols:
                x, y = c[0] % p, c[1] % p
                if y: pts.add((x*pow(y, -1, p) % p, 1))
                else: pts.add((1, 0))
            if len(pts) == p+1: return False
        return True
    # crit says SOME realization with det d has primitive columns; check by brute force over
    # multipliers M (small entries) with det = +-d
    found = False
    rng = range(-abs(d)-1, abs(d)+2)
    for a_, b_, c_, e_ in it.product(rng, repeat=4):
        if abs(a_*e_-b_*c_) != abs(d): continue
        cols = [(a_*c[0]+b_*c[1], c_*c[0]+e_*c[1]) for c in Bsat_cols]
        if all(gcd(*c) == 1 for c in cols): found = True; break
    if found != crit(): fails += 1
    if tested >= 60: break
print(f"III.10 saturation/realization and primitive-column criterion: {tested} random cases, failures {fails}")
