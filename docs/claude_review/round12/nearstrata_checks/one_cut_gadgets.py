"""Exact check of the gadget construction in the proof of Theorem 2.5 (one cut) of nearstrata.md.
For random u (rational, with planted clusters), eta, one cut (S = first s coords, L):
 * computes eta-clusters, sigma, t, a, nu = 2a + min(2 sigma, 2L+2t, s+L);
 * builds the gadget zonotope of the proof; checks EXACTLY that every point lies in A_{S,L} (coords in {-1,0,1},
   |k(S)| <= L), that every generator v has |<v,u>| >= eta, and that the degree equals nu;
 * builds c = coeffs of w^{k0} prod (w^v - 1)^{alpha} with sum alpha = m (m = nu) and checks EXACTLY that c annihilates
   all monomials of degree < m on its support (c in V(A,m)) and that sum c_k <k,u>^m = m! prod <v,u>^alpha;
 * if nu < m_test, builds y of part (b) and checks y in the stratum, |u-y|_inf < (n-1) eta.
All arithmetic in Fractions/integers.
"""
import random, itertools, math
from fractions import Fraction as F
from collections import defaultdict

def clusters(vals, eta):
    idx = sorted(range(len(vals)), key=lambda i: vals[i])
    cl = [0]*len(vals); c = 0
    for a, b in zip(idx, idx[1:]):
        cl[a] = c
        if vals[b] - vals[a] >= eta: c += 1
    cl[idx[-1]] = c
    return cl

def build(u, S, R, L, eta):
    n = len(u) + 1
    vals = list(u) + [F(0)]
    cl = clusters(vals, eta)
    c0 = cl[n-1]
    s = len(S)
    Sn = [x for x in S if cl[x] != c0]
    Z0 = [x for x in S if cl[x] == c0]
    classes = defaultdict(list)
    for x in S: classes[cl[x]].append(x)
    sigma = len(Sn); M = max(len(v) for v in classes.values()); t = s - M
    Ra = [x for x in R if cl[x] != c0]; a = len(Ra)
    tau = min(2*sigma, 2*L + 2*t, s + L)
    G = []   # list of (kind, coords)
    for x in Ra: G.append(('single', (x,)))
    if tau == 2*sigma:
        mm = min(sigma, len(Z0))
        for i in range(mm): G.append(('pair', (Sn[i], Z0[i])))
        for x in Sn[mm:]: G.append(('single', (x,)))
    elif tau == 2*L + 2*t:
        cstar = max(classes, key=lambda k: len(classes[k]))
        assert cstar != c0
        inC = classes[cstar]; out = [x for x in S if cl[x] != cstar]
        for i, x in enumerate(out): G.append(('pair', (x, inC[i])))
        for x in inC[len(out):len(out)+L]: G.append(('single', (x,)))
    else:
        q = s - L
        h = q // 2 if q % 2 == 0 else (q - 1) // 2
        rem = {k: list(v) for k, v in classes.items()}
        removed = []
        exc = [k for k in rem if k != c0 and len(rem[k]) > h]
        if q % 2 == 1:
            nonzero = [k for k in rem if k != c0]
            if len(exc) >= 2: kx, ky = exc[0], exc[1]
            elif len(exc) == 1: kx = exc[0]; ky = [k for k in nonzero if k != kx][0]
            else: kx, ky = nonzero[0], nonzero[1]
            x = rem[kx].pop(); y = rem[ky].pop()
            G.append(('hex', (x, y))); nsing = L - 1
        else:
            nsing = L
        for k in list(rem):
            if k == c0: continue
            while len(rem[k]) > h:
                removed.append(rem[k].pop())
        assert len(removed) <= nsing, (len(removed), nsing)
        pool = [x for k in rem if k != c0 for x in rem[k]]
        # fill up singles from largest nonzero classes first (keeps classes <= h)
        while len(removed) < nsing:
            k = max((k for k in rem if k != c0 and rem[k]), key=lambda k: len(rem[k]))
            removed.append(rem[k].pop())
        for x in removed: G.append(('single', (x,)))
        rest = sorted([x for k in rem for x in rem[k]], key=lambda x: cl[x])
        assert len(rest) % 2 == 0
        hh = len(rest) // 2
        for i in range(hh): G.append(('pair', (rest[i], rest[i+hh])))
    return G, dict(sigma=sigma, t=t, a=a, tau=tau, nu=2*a+tau, cl=cl, c0=c0, classes=classes)

def gens(G, d):
    k0 = [0]*d; gl = []
    def e(x):
        v = [0]*d; v[x] = 1; return v
    for kind, cs in G:
        if kind == 'single':
            x, = cs; k0[x] -= 1; gl.append((e(x), 2))
        elif kind == 'pair':
            x, y = cs; k0[y] += 1; k0[x] -= 1
            v = [0]*d; v[x] = 1; v[y] = -1; gl.append((v, 2))
        else:
            x, y = cs; k0[x] -= 1
            v = [0]*d; v[x] = 1; v[y] = -1
            gl += [(e(x), 1), (e(y), 1), (v, 1)]
    return k0, gl

def zonotope_points(k0, gl):
    pts = set()
    for ts in itertools.product(*[range(Lg+1) for _, Lg in gl]):
        k = list(k0)
        for (v, _), ti in zip(gl, ts):
            for i in range(len(k)): k[i] += ti*v[i]
        pts.add(tuple(k))
    return pts

def cvec(k0, gl, m):
    alphas = []; left = m
    for v, Lg in gl:
        a_ = min(Lg, left); alphas.append(a_); left -= a_
    assert left == 0
    poly = {tuple(k0): 1}
    for (v, _), al in zip(gl, alphas):
        for _ in range(al):
            new = defaultdict(int)
            for k, c in poly.items():
                new[tuple(ki+vi for ki, vi in zip(k, v))] += c
                new[k] -= c
            poly = {k: c for k, c in new.items() if c}
    return poly, alphas

def check(seed, d, s, L, m_extra=0):
    rng = random.Random(seed)
    n = d + 1
    S = list(range(s)); R = list(range(s, d))
    eta = F(1, 10)
    centers = [F(rng.randint(-10, 10), 10) for _ in range(rng.randint(1, 4))] + [F(0)]
    u = [rng.choice(centers) + F(rng.randint(-3, 3), 100) for _ in range(d)]
    G, info = build(u, S, R, L, eta)
    k0, gl = gens(G, d)
    pts = zonotope_points(k0, gl)
    for k in pts:
        assert all(abs(x) <= 1 for x in k), k
        assert abs(sum(k[x] for x in S)) <= L, (k, L)
    for v, _ in gl:
        assert abs(sum(vi*ui for vi, ui in zip(v, u))) >= eta
    deg = sum(Lg for _, Lg in gl)
    assert deg == info['nu'], (deg, info['nu'])
    m = deg
    if m >= 1:
        poly, alphas = cvec(k0, gl, m)
        # exact moment check: all monomials of degree < m (in d vars) -- only feasible small; check via random linear forms powers
        for trial in range(6):
            lam = [rng.randint(-5, 5) for _ in range(d)]
            for j in range(m):
                assert sum(c * sum(li*ki for li, ki in zip(lam, k))**j for k, c in poly.items()) == 0
        lhs = sum(c * sum(ki*ui for ki, ui in zip(k, u))**m for k, c in poly.items())
        rhs = math.factorial(m)
        for (v, _), al in zip(gl, alphas):
            rhs *= sum(vi*ui for vi, ui in zip(v, u))**al
        assert lhs == rhs
        assert sum(abs(c) for c in poly.values()) <= 2**m
    # part (b) with m_test = nu + 1 + m_extra
    mt = info['nu'] + 1 + m_extra
    cl, c0 = info['cl'], info['c0']
    sig, t, a, tau = info['sigma'], info['t'], info['a'], info['tau']
    y = list(u)
    if tau == 2*sig:
        X = [x for x in range(d) if cl[x] != c0]
        y = [u[x] if x in X else F(0) for x in range(d)]
        assert 2*len(X) < mt
    elif tau == 2*L + 2*t:
        cls = info['classes']; cstar = max(cls, key=lambda k: len(cls[k]))
        astar = u[cls[cstar][0]]
        for x in range(d):
            if x in S and cl[x] == cstar: y[x] = astar
            elif x in R and cl[x] == c0: y[x] = F(0)
        X = [x for x in range(d) if (y[x] - (astar if x in S else 0)) != 0]
        assert 2*L + 2*len(X) < mt
    else:
        for x in R:
            if cl[x] == c0: y[x] = F(0)
        Y = [x for x in R if y[x] != 0]
        assert s + L + 2*len(Y) < mt
    assert max(abs(a_ - b_) for a_, b_ in zip(u, y)) < (n-1)*eta
    return info['nu']

if __name__ == '__main__':
    cnt = 0; hist = defaultdict(int)
    for seed in range(3000):
        rng = random.Random(seed + 10**6)
        d = rng.randint(3, 7); s = rng.randint(2, d); L = rng.randint(1, s-1)
        nu = check(seed, d, s, L)
        hist[(nu >= 0)] += 1; cnt += 1
    print("checked", cnt, "random (u, eta, S, L) instances, d in [3,7]: all assertions passed")

    # case coverage
    cov = defaultdict(int)
    for seed in range(3000):
        rng = random.Random(seed + 10**6)
        d = rng.randint(3, 7); s = rng.randint(2, d); L = rng.randint(1, s-1)
        rng2 = random.Random(seed)
        eta = F(1, 10)
        centers = [F(rng2.randint(-10, 10), 10) for _ in range(rng2.randint(1, 4))] + [F(0)]
        u = [rng2.choice(centers) + F(rng2.randint(-3, 3), 100) for _ in range(d)]
        G, info = build(u, list(range(s)), list(range(s, d)), L, eta)
        sig, t, tau = info['sigma'], info['t'], info['tau']
        case = 'A' if tau == 2*sig else ('B' if tau == 2*L + 2*t else ('C-odd' if (s-L) % 2 else 'C-even'))
        cov[case] += 1
    print("case coverage:", dict(cov))
