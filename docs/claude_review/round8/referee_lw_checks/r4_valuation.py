"""Referee check R4: Lemma 2.1 and the valuation lemma (Lemma 2.2) of lw.md, with a wider pool than
lw_checks/check_valuation.py: all four units, the ramified prime 1+i, the inert prime 3, and a split
prime 17 outside J = {5, 13}.  Representations

    unit * prod_k alpha_k^{b_k} = u_5^{v_1} u_13^{v_2}        (u_p = pi_p / conj pi_p)

are enumerated exhaustively for K = 1, 2 numbers (pool: exponents in [-1,1] at the 8 places
2+i, 2-i, 3+2i, 3-2i, 1+i, 3, 4+i, 4-i, times the 4 units; |b_k| <= 3) and for K = 3 numbers
(pool: exponents in [-1,1] at the 6 places 2+i, 2-i, 3+2i, 3-2i, 1+i, 3; |b_k| <= 2).  Heights are
naive heights of primitive minimal polynomials, computed exactly.  The minimum of sum_k log H(alpha_k)
is compared with w(supp v).
Usage: python3 r4_valuation.py
"""
import itertools, math
import numpy as np
from fractions import Fraction

out = []
def log(s):
    print(s); out.append(s)

def gmul(x, y): return (x[0]*y[0] - x[1]*y[1], x[0]*y[1] + x[1]*y[0])
def gconj(x): return (x[0], -x[1])
def gnorm(x): return x[0]*x[0] + x[1]*x[1]
UNITS = [(1, 0), (0, 1), (-1, 0), (0, -1)]
PLACES = [(2, 1), (2, -1), (3, 2), (3, -2), (1, 1), (3, 0), (4, 1), (4, -1)]
NORMS = [5, 5, 13, 13, 2, 9, 17, 17]

def number(unit, ex, places):
    num, den = unit, (1, 0)
    for p, e in zip(places, ex):
        for _ in range(abs(e)):
            if e > 0: num = gmul(num, p)
            else: den = gmul(den, p)
    return num, den

def logH(num, den):
    t = gmul(num, gconj(den)); D = gnorm(den); A, B = t
    g = math.gcd(math.gcd(A, B), D); A //= g; B //= g; D //= g
    if B == 0:
        return math.log(max(abs(A), D)), True, 0.0, 0.0
    c = [D*D, -2*A*D, A*A + B*B]
    g3 = math.gcd(math.gcd(c[0], c[1]), c[2]); c = [x//g3 for x in c]
    H = max(abs(x) for x in c)
    Mahler = Fraction(max(c[0], c[2]))   # |a| max(1,|r|)^2 = max(a, c) for conj. roots
    return math.log(H), False, float(H), float(Mahler)

# ---- Lemma 2.1 (3): M <= H <= 2M, and 2h = max(P+,P-) for every pool element
pool8 = []
for ex in itertools.product((-1, 0, 1), repeat=8):
    for u in UNITS:
        num, den = number(u, ex, PLACES)
        lh, rat, H, Mh = logH(num, den)
        Pp = sum(e*math.log(n) for e, n in zip(ex, NORMS) if e > 0)
        Pm = sum(-e*math.log(n) for e, n in zip(ex, NORMS) if e < 0)
        if not rat:
            assert Mh <= H <= 2*Mh and abs(math.log(Mh) - max(Pp, Pm)) < 1e-9
        else:
            assert abs(lh - max(Pp, Pm)/2) < 1e-9
        pool8.append((ex, u, lh))
log(f"Lemma 2.1: {len(pool8)} numbers (8 places, exponents in [-1,1], 4 units): (1)-(3) PASS")

# best log H per valuation vector (unit chosen optimally; the unit is free in the identity)
best8 = {}
for ex, u, lh in pool8:
    if ex not in best8 or lh < best8[ex]: best8[ex] = lh
w5, w13 = math.log(5), math.log(13)
def target(v): return (v[0], -v[0], v[1], -v[1], 0, 0, 0, 0)

def run(v, keys, best, K, bmax):
    tg = np.array(target(v)[:len(keys[0])])
    E = np.array(keys); LH = np.array([best[k] for k in keys])
    index = {k: i for i, k in enumerate(keys)}
    mn = math.inf
    bs = [b for b in range(-bmax, bmax+1) if b]
    if K == 1:
        for b in bs:
            r = tg / b
            if np.all(r == np.round(r)):
                k = tuple(int(x) for x in r)
                if k in index: mn = min(mn, LH[index[k]])
        return mn
    if K == 2:
        for b1 in bs:
            for b2 in bs:
                R = (tg[None, :] - b1*E)
                ok = np.all(R % b2 == 0, axis=1)
                for i in np.nonzero(ok)[0]:
                    k = tuple(int(x) for x in R[i] // b2)
                    if k in index: mn = min(mn, LH[i] + LH[index[k]])
        return mn
    # K == 3
    for b1 in bs:
        for b2 in bs:
            for b3 in bs:
                for i1 in range(len(keys)):
                    R = tg[None, :] - b1*E[i1][None, :] - b2*E
                    ok = np.all(R % b3 == 0, axis=1)
                    for i2 in np.nonzero(ok)[0]:
                        k = tuple(int(x) for x in R[i2] // b3)
                        if k in index:
                            mn = min(mn, LH[i1] + LH[i2] + LH[index[k]])
    return mn

keys8 = [k for k in best8 if any(k)]
targets = [(1, 0), (0, 1), (1, 1), (1, -1), (2, 1), (2, -1), (2, 0), (3, -2)]
for v in targets:
    ws = (w5 if v[0] else 0) + (w13 if v[1] else 0)
    m1 = run(v, keys8, best8, 1, 3); m2 = run(v, keys8, best8, 2, 3)
    log(f"v={v}: w(supp v)={ws:.4f}; min sum log H: K=1 {m1:.4f}, K=2 {m2:.4f} (8 places, |b|<=3); "
        f"Lemma 2.2 {'PASS' if min(m1, m2) >= ws - 1e-12 else 'FAIL'}")

# K = 3 with the 6-place pool
best6 = {}
for ex, u, lh in pool8:
    if ex[6] == 0 and ex[7] == 0:
        k = ex[:6]
        if k not in best6 or lh < best6[k]: best6[k] = lh
keys6 = [k for k in best6 if any(k)]
for v in [(1, 0), (1, 1), (1, -1), (2, -1)]:
    ws = (w5 if v[0] else 0) + (w13 if v[1] else 0)
    m3 = run(v, keys6, best6, 3, 2)
    log(f"v={v}: K=3 (6 places, |b|<=2): min sum log H = {m3:.4f} vs w(supp v) = {ws:.4f}: "
        f"{'PASS' if m3 >= ws - 1e-12 else 'FAIL'}")
open('r4_valuation_out.txt', 'w').write('\n'.join(out) + '\n')
