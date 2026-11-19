import sys, itertools
from f4 import *
def par(x): return bin(x).count('1') & 1
def run(r, L=None, verbose=True):
    m = 2 ** r; d = r + 1; Q = 2 ** d
    V = list(range(Q))
    forms = [c for c in range(Q) if c & 1]              # l(v) = <c,v>, l(1) = c_0 = 1
    assert len(forms) == m
    if L is None: L = (d - 2) * 2 ** (d - 1)
    Y = (1 << L, 1)
    def lin(v): return ((1 << L) ^ v, 1)              # Y + v
    H = {c: [v for v in V if par(c & v) == 0] for c in forms}
    z = {}
    for c in forms:
        p = ONE
        for v in H[c]:
            if v: p = mul4(p, lin(v))
        z[c] = p
    # common norm, computed separately
    N = 1
    for v in V:
        if v not in (0, 1): N = mul2(N, norm4(lin(v)) if False else 1) if False else N
    Nsep = ONE
    for v in V:
        if v not in (0, 1): Nsep = mul4(Nsep, lin(v))
    assert Nsep[1] == 0; N = Nsep[0]
    D = (m - 1) * L
    for c in forms:
        assert norm4(z[c]) == N and deg4(z[c]) == D
        # conj(z_H) = product over complementary side minus {1}
        p = ONE
        for v in V:
            if par(c & v) == 1 and v != 1: p = mul4(p, lin(v))
        assert conj4(z[c]) == p
    assert len(set(z.values())) == m
    assert deg2(N) == (Q - 2) * L == 2 * D
    maxd = -1; nid = 0
    for c1, c2 in itertools.combinations(forms, 2):
        J = [v for v in V if par(c1 & v) == 0 and par(c2 & v) == 0]
        PJ1 = 1
        for v in J: PJ1 = mul2(PJ1, 1 ^ v)
        PJY = (PJ1, 0)
        for v in J:
            if v: PJY = mul4(PJY, lin(v))
        diff = add4(z[c1], z[c2])
        assert diff == PJY; nid += 1
        dd = deg4(diff)
        S = sum(deg2(1 ^ v) for v in J)
        assert dd == (Q // 4 - 1) * L + S
        maxd = max(maxd, dd)
    # separability of ratios z_a/z_b: derivative of ratio nonzero
    sep = all(add4(mul4(der4(z[a]), z[b]), mul4(z[a], der4(z[b]))) != (0, 0)
              for a, b in itertools.combinations(forms, 2))
    g = 0
    for c in forms: g = gcd2(g, z[c][0]) if False else g
    if verbose:
        print(f"m={m} d={d} L={L} D={D} degN={deg2(N)} maxdiff={maxd} D/2={D/2} pairs={nid} "
              f"cluster={2*maxd<=D} ratios_separable={sep} N_sqfree={gcd2(N, der2(N))==1}")
    return maxd, D
for r in (2, 3):
    run(r)
# minimality / L-range: maxdiff - D/2 = Smax - L/2
for r in (2, 3):
    d = r + 1; Lm = (d - 2) * 2 ** (d - 1)
    for L in (Lm - 1, Lm, Lm + 1, Lm + 5):
        mx, D = run(r, L, verbose=False)
        print(f"  r={r} L={L}: maxdiff-D/2 = {mx - D/2}, Smax-L/2 = {(d-2)*2**(d-2) - L/2}, c=0 cluster: {mx <= D//2}")
