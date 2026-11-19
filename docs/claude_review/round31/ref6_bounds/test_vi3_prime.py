# Per-prime check of the VI.3 range argument and of depth amplification (VI.2).
import random, sys, itertools
from sympy import Poly, symbols, GF
from f4lib import *
from clusters import span
from test_vi3 import coords, parity

x = symbols('x')
def to_sym(a):
    return Poly([ (a >> k) & 1 for k in range(deg(a), -1, -1)], x, modulus=2)
def from_sym(P):
    c = P.all_coeffs(); r = 0
    for b in c:
        r = (r << 1) | (int(b) % 2)
    return r
def real_factors(a):
    _, fl = to_sym(a).factor_list()
    return [from_sym(P) for P, e in fl]

def ord4(a, f):
    k = 0
    while d4(a) >= 0:
        q, r = divmod4(a, f)
        if d4(r) >= 0:
            break
        a = q; k += 1
    return k

def instance(rng, d, maxdeg):
    while True:
        basis = [rng.randrange(1, 1 << maxdeg) for _ in range(d)]
        V = span(basis)
        if len(set(V)) == 1 << d:
            break
    co = coords(basis)
    B = rng.choice([v for v in V if v and deg(v) > 0])
    T = rng.randrange(1 << rng.randrange(1, 14))
    c = T
    for b in basis:
        c = gcd2(c, b)
    if c != 1:
        return None
    Q = 1 << d
    forms = [l for l in range(1, Q) if parity(l & co[B])]
    sides = [[v for v in V if parity(l & co[v]) == rng.randrange(2) or False] for l in forms]
    # build sides properly: S = l^{-1}(c)
    sides = []
    for l in forms:
        c = rng.randrange(2)
        sides.append([v for v in V if parity(l & co[v]) == c])
    Y = (T, B)
    lin = {v: add4(Y, (v, 0)) for v in V}
    F = [prod4([lin[v] for v in S]) for S in sides]
    G = F[0]
    for z in F[1:]:
        G = gcd4(G, z)
    Dp = d4(F[0]) - d4(G)
    Nred = mul4(prod4([lin[v] for v in sides[0]]), conj4(prod4([lin[v] for v in sides[0]])))
    GG = mul4(G, conj4(G))
    Nred, rem = divmod4(Nred, GG); assert d4(rem) < 0 and is_real(Nred)
    # pairs, variability
    pairs = []; seen = set()
    for v in V:
        if v in seen: continue
        seen |= {v, v ^ B}; pairs.append(v)
    sets = [set(S) for S in sides]
    variable = [v for v in pairs if any(v in S for S in sets) and any((v ^ B) in S for S in sets)]
    fixed = [v for v in pairs if v not in variable]
    # primes
    primes = set()
    for v in V:
        Nv = mul4(lin[v], conj4(lin[v]))[0]
        if deg(Nv) > 0:
            primes |= set(real_factors(Nv))
    res = []
    for p in primes:
        if deg(p) % 2 == 1:
            continue
        f = None
        for v in V:
            g = gcd4((p, 0), lin[v])
            if 0 < d4(g) < deg(p):
                f = g; break
        if f is None:
            continue  # p divides no Y+v properly -> inert-like (p | B and p | T+v); skip
        o = {v: ord4(lin[v], f) for v in V}
        vals = [sum(o[v] - o[v ^ B] for v in S) for S in sides]
        rng_ = max(vals) - min(vals)
        need = sum(abs(o[v] - o[v ^ B]) for v in variable)
        ep = 0; a = Nred[0]
        while True:
            q, r = divmod2(a, p)
            if r: break
            a = q; ep += 1
        res.append((rng_ >= need, 2 * ep >= rng_, rng_, need, ep))
    return dict(fixed=len(fixed), res=res, Dp=Dp)

if __name__ == '__main__':
    rng = random.Random(int(sys.argv[1]))
    n = int(sys.argv[2])
    tot = 0; bad = 0; badfix = 0; tight = 0; nz = 0
    for _ in range(n):
        d = rng.choice([3, 4])
        I = instance(rng, d, rng.choice([d, 5, 7]))
        if I is None: continue
        if I['fixed'] > 1: badfix += 1
        for ok1, ok2, r_, need, ep in I['res']:
            tot += 1
            if not (ok1 and ok2): bad += 1; print("FAIL", r_, need, ep)
            if need > 0:
                nz += 1
                if r_ == need: tight += 1
    print(f"primes checked {tot}, with nonzero required range {nz}, failures {bad}, tight (range == need) {tight}, instances with >1 fixed pair {badfix}")
    # depth amplification test
    fails = 0
    for _ in range(300):
        n_ = rng.randrange(2, 30)
        f = (rng.randrange(1 << n_) | (1 << n_), rng.randrange(1 << n_))
        g = (rng.randrange(1 << n_) | (1 << n_), rng.randrange(1 << n_))
        g = (g[0] ^ (g[1] >> n_ << n_), g[1] & ((1 << n_) - 1))
        f = (f[0], f[1] & ((1 << n_) - 1))
        if f == g or lc4(f) != 1 or lc4(g) != 1: continue
        dep = depth_ratio(f, g)
        for e in range(1, 25):
            a = (e & -e).bit_length() - 1
            if depth_ratio(pow4(f, e), pow4(g, e)) != (1 << a) * dep:
                fails += 1
    print("depth amplification failures:", fails)
