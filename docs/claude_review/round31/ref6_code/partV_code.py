"""(d) Anchored parity code of the genuine Part V clusters over F_2[t] (norm of F_4(t)/F_2(t)).
Own construction: d=r+1, m=2^r, V=F_2[t]_{<d}, H=ker(l) for forms l with l(1)=1, Y=t^L+theta,
L=(d-2)2^(d-1), z_H=prod_{v in H, v!=0}(Y+v).  Verifies norms, cluster condition, factors N over F_2,
computes parities a_{H,p}-a_{H0,p} mod 2 at split primes, the minimum weighted weight of the code and
compares with K+1, K=ceil(D/2); checks deg H - deg H' >= K+1 on every nonzero word."""
import sys, itertools
sys.path.insert(0, '.')
from gf2 import *

def run(r):
    d = r + 1
    m = 2 ** r
    L = (d - 2) * 2 ** (d - 1)
    if r == 1:
        L = d
    V = list(range(2 ** d))
    forms = [w for w in range(2 ** d) if w & 1]
    assert len(forms) == m
    hyper = [[v for v in V if bin(v & w).count('1') % 2 == 0] for w in forms]
    def f_v(v):  # Y + v in F_4[t]
        a = [0] * (L + 1)
        a[L] = 1
        for i in range(d):
            if (v >> i) & 1:
                a[i] ^= 1
        a[0] ^= 2
        return p4_trim(a)
    fv = {v: f_v(v) for v in V}
    pts = []
    for H in hyper:
        z = [1]
        for v in H:
            if v:
                z = p4_mul(z, fv[v])
        pts.append(z)
    D = len(pts[0]) - 1
    assert D == (m - 1) * L
    Nn = p4_mul(pts[0], p4_conj(pts[0]))
    for z in pts:
        assert len(z) - 1 == D and z[-1] == 1
        assert p4_mul(z, p4_conj(z)) == Nn
    N2 = p4_to_gf2(Nn)
    maxdiff = max(len(p4_add(a, b)) - 1 for a, b in itertools.combinations(pts, 2))
    assert maxdiff <= D // 2
    # factor N over F_2 via factors of the norms of f_v
    primes = set()
    for v in V:
        if v in (0, 1):
            continue
        nv = p4_to_gf2(p4_mul(fv[v], p4_conj(fv[v])))
        primes.update(factor(nv).keys())
    fac = factor(N2)
    assert set(fac) == primes
    split = [p for p in primes if deg(p) % 2 == 0]
    inert = [p for p in primes if deg(p) % 2 == 1]
    # choose pi for each split p
    pis = {}
    for p in split:
        P4 = p4_from_gf2(p)
        found = None
        for v in V:
            g = p4_gcd(P4, fv[v])
            if len(g) - 1 == deg(p) // 2:
                found = g
                break
        assert found is not None, "no half-degree gcd found"
        assert p4_mul(found, p4_conj(found)) == P4
        pis[p] = found
    # valuations
    ordv = {(p, v): p4_ord(fv[v], pis[p]) for p in split for v in V if v}
    a = [[sum(ordv[(p, v)] for v in H if v) for p in split] for H in hyper]
    # check exponent consistency: a_p + abar_p = e_p
    for i, H in enumerate(hyper):
        for j, p in enumerate(split):
            ab = sum(p4_ord(fv[v], p4_conj(pis[p])) for v in H if v)
            assert a[i][j] + ab == fac[p]
    # anchored parity patterns
    cols = {}
    for j, p in enumerate(split):
        pat = tuple((a[i][j] - a[0][j]) % 2 for i in range(1, m))
        if any(pat):
            cols.setdefault(pat, []).append(p)
    n = len(cols)
    pats = list(cols)
    weights = [sum(deg(p) for p in cols[pt]) for pt in pats]
    K = -(-D // 2)
    W = 2 * D
    # rank
    rows = [sum(pt[i] << k for k, pt in enumerate(pats)) for i in range(m - 1)]
    piv = {}
    for row in rows:
        while row:
            l = row.bit_length()
            if l not in piv:
                piv[l] = row; break
            row ^= piv[l]
    rank = len(piv)
    best = None
    worst_slack = None
    for c in range(1, 2 ** (m - 1)):
        word = 0
        for i in range(m - 1):
            if (c >> i) & 1:
                word ^= rows[i]
        if word == 0:
            continue
        wt = sum(weights[k] for k in range(n) if (word >> k) & 1)
        Hp = 1
        for k in range(n):
            if (word >> k) & 1:
                for p in cols[pats[k]]:
                    Hp = mul(Hp, p)
        dH = deriv(Hp)
        assert dH != 0
        ordinf = deg(Hp) - deg(dH)
        assert ordinf >= K + 1, (c, ordinf, K)
        if best is None or wt < best:
            best = wt
        sl = ordinf - (K + 1)
        worst_slack = sl if worst_slack is None else min(worst_slack, sl)
    print(f"m={m} d={d} L={L} D={D} W={W} K={K}: #split base primes={len(split)} (degrees "
          f"{sorted(deg(p) for p in split)[:6]}...), #inert={len(inert)}, max deg(z_i-z_j)={maxdiff}")
    print(f"   distinct nonzero anchored patterns n={n}, binary rank={rank}, min weighted weight={best}, "
          f"K+1={K+1}, W/4={W/4}; min over words of [deg H - deg H' - (K+1)] = {worst_slack}")
    print(f"   total pattern weight={sum(weights)} (<= W={W});  rank/n={rank/n:.3f}")
    return best, K

for r in (2, 3, 4):
    run(r)
