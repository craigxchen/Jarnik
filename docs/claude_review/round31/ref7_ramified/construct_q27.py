# Independent reproduction of VII.1 with p = 3, q = 27, k = 2 in a ramified model eta^2 = Delta (or the constant
# model Delta = -1).  Labels v = eta^3 C^3 (x0 + x1 t^3), x in F_27^2.  We never leave F_3[t,eta]: each Frobenius
# orbit of labels is multiplied symbolically through the elementary symmetric functions of w = x0 + x1 t^3.
import sys, itertools, time, random
import numpy as np
from fp import *
import f27 as F

Z0 = np.zeros(0, dtype=np.int64)
ONE = (const(1), Z0)

def lin(f, x):
    return F.fadd(F.fmul(f[0] % 3, x[0]), F.fmul(f[1] % 3, x[1]))

def orbit(x):
    o = [x]
    while True:
        y = (F.frob(o[-1][0]), F.frob(o[-1][1]))
        if y == x: return o
        o.append(y)

# F_27[t] polys as lists of ints (low first), tiny
def pmul27(a, b):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            r[i + j] = F.fadd(r[i + j], F.fmul(x, y))
    return r
def padd27(a, b):
    n = max(len(a), len(b)); r = [0] * n
    for i, x in enumerate(a): r[i] = F.fadd(r[i], x)
    for i, x in enumerate(b): r[i] = F.fadd(r[i], x)
    return r

def elem_sym(ws):
    # prod (Z - w_i) = sum_k (-1)^k e_k Z^{r-k}; return e_k as F_3 polys (assert coefficients in F_3)
    poly = [[1]]  # coefficients in Z (low first), each an F_27[t] list
    for w in ws:
        negw = [F.fneg(c) for c in w]
        new = [[0]] * (len(poly) + 1)
        for i, c in enumerate(poly):
            new[i + 1] = padd27(new[i + 1], c)          # * Z
            new[i] = padd27(new[i], pmul27(c, negw))    # * (-w)
        poly = new
    r = len(ws)
    es = []
    for k in range(r + 1):
        c = poly[r - k]   # coefficient of Z^{r-k} = (-1)^k e_k
        assert all(v in (0, 1, 2) for v in c), "not in F_3"
        ek = trim(np.array(c)) if k % 2 == 0 else neg(trim(np.array(c)))
        es.append(ek)
    return es

def build(Delta, C, T, forms, signs, verbose=True):
    t0 = time.time()
    O = Order(Delta)
    U = mul(O.D, powp(C, 3))            # Y = eta * U = eta^3 C^3
    labels = list(itertools.product(range(27), repeat=2))
    E = {x: [1 + sg * F.chi(lin(f, x)) for f, sg in zip(forms, signs)] for x in labels}
    gv = {x: min(E[x]) for x in labels}
    g = sum(gv.values())
    seen = set(); orbs = []
    for x in labels:
        if x in seen: continue
        o = orbit(x); seen.update(o)
        # invariance of exponents along the orbit (descent)
        assert all(E[y] == E[x] for y in o)
        orbs.append(o)
    # powers needed
    Tp = [const(1)]
    for i in range(3): Tp.append(mul(Tp[-1], T))
    Up = [const(1)]
    for i in range(3): Up.append(mul(Up[-1], U))
    Dp = [const(1), O.D]
    def orbit_elem(o):
        ws = [[x[0], 0, 0, x[1]] for x in o]   # w = x0 + x1 t^3
        es = elem_sym(ws)
        r = len(o)
        A = Z0; B = Z0
        for k in range(r + 1):
            term = mul(mul(es[k], Up[k]), Tp[r - k])
            if k % 2 == 0:
                A = add(A, scal((-1) ** k, mul(term, powp(O.D, k // 2))))
            else:
                B = add(B, scal((-1) ** k, mul(term, powp(O.D, (k - 1) // 2))))
        return (A, B)
    oe = [orbit_elem(o) for o in orbs]
    if verbose: print('orbits', len(orbs), 'g', g, 'time', round(time.time() - t0, 1), flush=True)
    def prodtree(items):
        items = [x for x in items]
        if not items: return ONE
        while len(items) > 1:
            nxt = []
            for i in range(0, len(items) - 1, 2):
                nxt.append(O.mul(items[i], items[i + 1]))
            if len(items) % 2: nxt.append(items[-1])
            items = nxt
        return items[0]
    m = len(forms)
    zs = []
    for j in range(m):
        its = []
        for o, el in zip(orbs, oe):
            e = E[o[0]][j] - gv[o[0]]
            its += [el] * e
        zs.append(prodtree(its))
        if verbose: print('mark', j, 'built', round(time.time() - t0, 1), flush=True)
    R = prodtree(oe)
    Gits = []
    for o, el in zip(orbs, oe):
        Gits += [el] * gv[o[0]]
    G = prodtree(Gits)
    out = dict(g=g, n=len(labels), L=deg(T))
    out['R_in_F3[t]'] = len(R[1]) == 0
    Rp = R[0]
    out['degR'] = deg(Rp)
    # R mod C: all primes of C must not divide R (they are where R' = a0 T' can vanish)
    _, rc = divmod_(Rp, C)
    out['gcd(R,C)'] = deg(gcd(C, rc)) if len(rc) else deg(C)
    norms = [O.norm(z) for z in zs]
    out['equal_norms'] = all(np.array_equal(norms[0], x) for x in norms)
    N0 = norms[0]
    out['degN0'] = deg(N0)
    NG = O.norm(G)
    out['N0*Norm(G)==R^2'] = np.array_equal(mul(N0, NG), mul(Rp, Rp))
    out['D2'] = sorted(set(O.pole2(z) for z in zs))
    out['distinct'] = len({(z[0].tobytes(), z[1].tobytes()) for z in zs}) == m
    out['lead'] = sorted(set((int(z[0][-1]) if len(z[0]) else 0, deg(z[0]), deg(z[1])) for z in zs))
    w = []
    for i in range(m):
        for j in range(i + 1, m):
            w.append(O.pole2(O.sub(zs[i], zs[j])) // 2)
    out['pairdiff'] = w
    D = out['D2'][0] / 2
    out['D'] = D
    out['critical'] = max(w) <= D / 2
    # independent modular cross-check of the FFT products: reduce factors mod a random f, multiply mod f
    random.seed(1)
    f = trim(np.array([random.randrange(3) for _ in range(40)] + [1]))
    def red(a): return divmod_(a, f)[1] if deg(a) >= deg(f) else a
    def omul(x, y):
        A, B = x; Cc, Ee = y
        return (red(add(mul(A, Cc), mul(red(O.D), mul(B, Ee)))), red(add(mul(A, Ee), mul(B, Cc))))
    def redel(x): return (red(x[0]), red(x[1]))
    oer = [redel(el) for el in oe]
    ok = True
    for j in range(m):
        acc = ONE
        for o, el in zip(orbs, oer):
            for _ in range(E[o[0]][j] - gv[o[0]]):
                acc = omul(acc, el)
        zr = redel(zs[j])
        ok &= np.array_equal(acc[0], zr[0]) and np.array_equal(acc[1], zr[1])
    out['modf_crosscheck'] = ok
    me, dec = max_exponent(N0) if deg(N0) < 60000 else (None, None)
    out['N0_max_exponent_direct'] = me
    out['time'] = round(time.time() - t0, 1)
    return out, zs, Rp, N0
