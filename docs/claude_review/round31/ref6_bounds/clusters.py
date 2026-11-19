# Build the Part V clusters and the VI.5 squarefree clusters with own code.
from f4lib import *
import itertools

def sq2(a):
    r = 0; k = 0
    while a:
        if a & 1:
            r |= 1 << (2 * k)
        a >>= 1; k += 1
    return r

def span(basis):
    out = [0]
    for b in basis:
        out = out + [x ^ b for x in out]
    return out

def trace2(x, p):
    e = deg(p); x = mod2(x, p); s = 0; y = x
    for _ in range(e):
        s ^= y
        y = mod2(clmul(y, y), p)
    assert s in (0, 1), s
    return s

def solve_gf2(rows, rhs, nvars):
    # rows: list of int bitmasks over nvars; returns a solution bitmask
    piv = {}
    for r, b in zip(rows, rhs):
        for c in sorted(piv, reverse=True):
            pass
        # reduce
        for c, (pr, pb) in piv.items():
            if (r >> c) & 1:
                r ^= pr; b ^= pb
        if r == 0:
            if b:
                return None
            continue
        c = r.bit_length() - 1
        # eliminate c from existing pivots
        for c2 in list(piv):
            pr, pb = piv[c2]
            if (pr >> c) & 1:
                piv[c2] = (pr ^ r, pb ^ b)
        piv[c] = (r, b)
    x = 0
    for c, (pr, pb) in piv.items():
        # free vars = 0; pr has pivot c and possibly other non-pivot bits (=0)
        if pb:
            x |= 1 << c
    return x

def hyperplane_marks(V, E0basis, Y):
    """marks z_H for hyperplanes H of V excluding 1, V = <1> + E0.
    A hyperplane excluding 1 <-> functional l on V with l(1)=1, H = ker l;
    l is determined by its values on E0basis (any). Returns list of (signs dict, z)."""
    k = len(E0basis)
    E0 = span(E0basis)
    marks = []
    for vals in itertools.product([0, 1], repeat=k):
        # l(sum c_i e_i) = sum c_i vals_i ; l(1) = 1
        def l_of(v):
            # decompose v = a + sum c_i e_i using E0 enumeration
            return LMAP[v]
        H = []
        LMAP = {}
        for idx in range(1 << k):
            e = 0; val = 0
            for i in range(k):
                if (idx >> i) & 1:
                    e ^= E0basis[i]; val ^= vals[i]
            LMAP[e] = val
            LMAP[e ^ 1] = val ^ 1
        H = [v for v in V if LMAP[v] == 0]
        z = prod4([add4(Y, (v, 0)) for v in H if v != 0])
        signs = {e: (+1 if LMAP[e] == 0 else -1) for e in E0 if e != 0}
        marks.append((vals, signs, z))
    return marks

def partV(r):
    m = 2 ** r; d = r + 1; L = m * (r - 1)
    E0basis = [1 << i for i in range(1, d)]
    V = span([1] + E0basis)
    Y = (1 << L, 1)
    return dict(m=m, L=L, V=V, E0basis=E0basis, Y=Y, T=1 << L, B=1)

def sqfree_VI5(d, L=None):
    m = 2 ** (d - 1)
    if L is None:
        L = 4 * m * (d - 2)
    Ub = [1 << i for i in range(1, d)]
    E0basis = [sq2(sq2(f) ^ f) for f in Ub]
    V = span([1] + E0basis)
    primes = set()
    for v in V:
        if v in (0, 1):
            continue
        for p in factor2(v):
            if deg(p) % 2 == 0:
                primes.add(p)
    primes = sorted(primes)
    R = 1
    for p in primes:
        R = clmul(R, p)
    dR = deg(R)
    E = (1 << L) ^ 2 if L % 2 == 0 else (1 << L)
    rows = []; rhs = []
    for p in primes:
        e = deg(p)
        row = 0
        for k in range(dR):
            if trace2(1 << k, p):
                row |= 1 << k
        rows.append(row)
        rhs.append((1 + e // 2 + trace2(E, p)) & 1)
    S = solve_gf2(rows, rhs, dR)
    assert S is not None
    for p, b in zip(primes, rhs):
        assert trace2(S, p) == b
    T = E ^ sq2(S)
    assert deg(T) == L
    Y = (T, 1)
    return dict(m=m, L=L, V=V, E0basis=E0basis, Y=Y, T=T, B=1, R=R, degR=dR, primes=primes)
