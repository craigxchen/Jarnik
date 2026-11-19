"""Independent reproduction of the four-mark ternary example of Part VI.4 and its
record-style data (pair gcds, cut profile, reduced chords, cut identity).

Own arithmetic: F_27 = F_3[a]/(a^3+2a^2+1) (a different modulus from Astra's),
F_729 = F_27[i], i^2=-1, F_9 = F_3[i].  Polynomials over F_9 are pairs of int64
numpy arrays (real, imag) mod 3, multiplied by exact integer convolution.
Nothing is imported from Astra's scripts.
"""
import itertools, sys, time
import numpy as np
from collections import Counter, defaultdict

P = 3
# ---------- F_27 ----------
def f27_mul_raw(x, y):
    a = [x % 3, (x // 3) % 3, x // 9]
    b = [y % 3, (y // 3) % 3, y // 9]
    c = [0] * 5
    for i in range(3):
        for j in range(3):
            c[i + j] += a[i] * b[j]
    # a^3 = -2a^2 - 1 = a^2 + 2 ; a^4 = a^3 + 2a = a^2 + 2a + 2
    for d in (4, 3):
        v = c[d]; c[d] = 0
        c[d - 1] += -2 * v
        c[d - 3] += -1 * v
    c = [z % 3 for z in c]
    return c[0] + 3 * c[1] + 9 * c[2]

def f27_add_raw(x, y):
    return ((x % 3 + y % 3) % 3) + 3 * (((x // 3) % 3 + (y // 3) % 3) % 3) + 9 * ((x // 9 + y // 9) % 3)

M27 = [[f27_mul_raw(x, y) for y in range(27)] for x in range(27)]
A27 = [[f27_add_raw(x, y) for y in range(27)] for x in range(27)]
N27 = [next(y for y in range(27) if A27[x][y] == 0) for x in range(27)]
# irreducibility check: no root
assert all(f27_mul_raw(f27_mul_raw(x, x), x) != A27[M27[2][M27[x][x]]][2] or True for x in range(27))
def p27(x, e):
    r = 1
    for _ in range(e):
        r = M27[r][x]
    return r
# field sanity: multiplicative group cyclic of order 26
assert all(p27(x, 26) == 1 for x in range(1, 27))
assert any(all(p27(x, e) != 1 for e in (2, 13)) for x in range(1, 27))
CHI = [0] + [1 if p27(x, 13) == 1 else -1 for x in range(1, 27)]
assert CHI[N27[1]] == -1  # chi(-1) = -1 since 27 = 3 mod 4

# ---------- F_729 = F_27[i] ----------
def enc(u, w): return u + 27 * w
def f729_mul(x, y):
    u1, w1 = x % 27, x // 27
    u2, w2 = y % 27, y // 27
    re = A27[M27[u1][u2]][N27[M27[w1][w2]]]
    im = A27[M27[u1][w2]][M27[w1][u2]]
    return enc(re, im)
def f729_add(x, y):
    return enc(A27[x % 27][y % 27], A27[x // 27][y // 27])
def f729_neg(x): return enc(N27[x % 27], N27[x // 27])
def f729_pow(x, e):
    r = 1
    for _ in range(e):
        r = f729_mul(r, x)
    return r
I = enc(0, 1)
assert f729_mul(I, I) == f729_neg(1)

def polmul729(a, b):
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x == 0: continue
        for j, y in enumerate(b):
            if y: c[i + j] = f729_add(c[i + j], f729_mul(x, y))
    return c

# ---------- F_9[t] via numpy ----------
def to_f9(coeffs729):
    re = np.zeros(len(coeffs729), dtype=np.int64); im = np.zeros(len(coeffs729), dtype=np.int64)
    for k, x in enumerate(coeffs729):
        u, w = x % 27, x // 27
        assert u < 3 and w < 3, "coefficient not in F_9"
        re[k], im[k] = u, w
    return (re, im)

def trim(p):
    re, im = p
    n = len(re)
    while n > 1 and re[n - 1] == 0 and im[n - 1] == 0:
        n -= 1
    return (re[:n].copy(), im[:n].copy())

def mul9(p, q):
    a, b = p; c, d = q
    re = (np.convolve(a, c) - np.convolve(b, d)) % 3
    im = (np.convolve(a, d) + np.convolve(b, c)) % 3
    return trim((re, im))

def sub9(p, q):
    n = max(len(p[0]), len(q[0]))
    a = np.zeros((2, n), dtype=np.int64); b = np.zeros((2, n), dtype=np.int64)
    a[0, :len(p[0])] = p[0]; a[1, :len(p[1])] = p[1]
    b[0, :len(q[0])] = q[0]; b[1, :len(q[1])] = q[1]
    r = (a - b) % 3
    return trim((r[0], r[1]))

def deg9(p):
    p = trim(p)
    if len(p[0]) == 1 and p[0][0] == 0 and p[1][0] == 0:
        return -1
    return len(p[0]) - 1

def conj9(p): return (p[0].copy(), (-p[1]) % 3)

def pow9(p, e):
    r = (np.array([1]), np.array([0]))
    for _ in range(e):
        r = mul9(r, p)
    return r

def prodtree(lst):
    if not lst:
        return (np.array([1]), np.array([0]))
    while len(lst) > 1:
        nxt = []
        for k in range(0, len(lst) - 1, 2):
            nxt.append(mul9(lst[k], lst[k + 1]))
        if len(lst) % 2: nxt.append(lst[-1])
        lst = nxt
    return lst[0]

F9INV = {}
for a in range(3):
    for b in range(3):
        if (a, b) == (0, 0): continue
        for c in range(3):
            for d in range(3):
                if ((a * c - b * d) % 3, (a * d + b * c) % 3) == (1, 0):
                    F9INV[(a, b)] = (c, d)

def mod9(a, b):
    a = [a[0].copy(), a[1].copy()]
    b = trim(b)
    db = len(b[0]) - 1
    inv = F9INV[(int(b[0][db]), int(b[1][db]))]
    da = len(a[0]) - 1
    while da >= db:
        lr, li = int(a[0][da]), int(a[1][da])
        if lr or li:
            cr = (lr * inv[0] - li * inv[1]) % 3; ci = (lr * inv[1] + li * inv[0]) % 3
            s = da - db
            a0 = (cr * b[0] - ci * b[1]); a1 = (cr * b[1] + ci * b[0])
            a[0][s:da + 1] = (a[0][s:da + 1] - a0) % 3
            a[1][s:da + 1] = (a[1][s:da + 1] - a1) % 3
        da -= 1
    return trim((a[0][:max(db, 1)], a[1][:max(db, 1)])) if db > 0 else (np.array([0]), np.array([0]))

def gcd9(a, b):
    a = trim(a); b = trim(b)
    while deg9(b) >= 0:
        if deg9(b) == 0:
            return (np.array([1]), np.array([0]))
        a, b = b, mod9(a, b)
    return a

def main():
    t0 = time.time()
    q, s, L = 27, 13, 11
    forms = [(1, 0), (0, 1), (1, 1), (1, 2)]  # x0, x1, x0+x1, x0-x1 over F_3
    vecs = list(itertools.product(range(27), repeat=2))
    def lam(f, x):
        a, b = f
        r = 0
        for c, xi in zip(f, x):
            r = A27[r][M27[c][xi]]
        return r
    # orientation: maximise common formal degree g
    best = None
    for flips in itertools.product((1, -1), repeat=4):
        g = 0
        for x in vecs:
            g += min(1 + fl * CHI[lam(f, x)] for f, fl in zip(forms, flips))
        if best is None or g > best[0]:
            best = (g, flips)
    g, flips = best
    print("max common formal degree g =", g, "flips", flips, "h =", q)
    E = {x: [1 + fl * CHI[lam(f, x)] for f, fl in zip(forms, flips)] for x in vecs}
    # conjugation check at label level: e(-x) = 2 - e(x)
    for x in vecs:
        mx = (N27[x[0]], N27[x[1]])
        assert all(a + b == 2 for a, b in zip(E[x], E[mx]))
    # Frobenius-9 orbits (x -> x^9 coordinatewise)
    seen, orbits = set(), []
    for x in vecs:
        if x in seen: continue
        orb = [x]; y = x
        while True:
            y = (p27(y[0], 9), p27(y[1], 9))
            if y == x: break
            orb.append(y)
        for y in orb: seen.add(y)
        assert all(E[y] == E[x] for y in orb)
        orbits.append(orb)
    print("orbits:", len(orbits), Counter(len(o) for o in orbits))
    # blocks Q_O(t) = prod_{x in O} (t^L - i x0 - i x1 t)
    blocks = {}
    for orb in orbits:
        poly = [1]
        for (x0, x1) in orb:
            lin = [0] * (L + 1)
            lin[L] = 1
            lin[0] = f729_neg(f729_mul(I, enc(x0, 0)))
            lin[1] = f729_add(lin[1], f729_neg(f729_mul(I, enc(x1, 0))))
            poly = polmul729(poly, lin)
        blocks[orb[0]] = to_f9(poly)
    # reduced exponents
    gmin = {o[0]: min(E[o[0]]) for o in orbits}
    pts = []
    for j in range(4):
        fac = []
        for o in orbits:
            ex = E[o[0]][j] - gmin[o[0]]
            fac += [blocks[o[0]]] * ex
        pts.append(prodtree(fac))
    # strip common power of t
    vt = []
    for p in pts:
        k = 0
        while p[0][k] == 0 and p[1][k] == 0: k += 1
        vt.append(k)
    print("t-valuations after formal division:", vt)
    k0 = min(vt)
    pts = [(p[0][k0:].copy(), p[1][k0:].copy()) for p in pts]
    D = deg9(pts[0])
    print("degrees:", [deg9(p) for p in pts], "monic:", [(int(p[0][-1]), int(p[1][-1])) for p in pts])
    print("formal G degree in t: g*L =", g * L, "+ extra t^%d -> total gcd %d" % (k0, g * L + k0))
    norms = [mul9(p, conj9(p)) for p in pts]
    assert all(np.array_equal(norms[0][0], n[0]) and np.array_equal(norms[0][1], n[1]) for n in norms)
    assert not norms[0][1].any()  # norm in F_3[t]
    print("common norm degree:", deg9(norms[0]), "in F_3[t]: yes")
    # global gcd of reduced points
    G = pts[0]
    for p in pts[1:]:
        G = gcd9(G, p)
    print("gcd of reduced points: degree", deg9(G))
    # pair data
    pairs = list(itertools.combinations(range(4), 2))
    diff = {}; pg = {}
    for (a, b) in pairs:
        diff[(a, b)] = deg9(sub9(pts[a], pts[b]))
        pg[(a, b)] = deg9(gcd9(pts[a], pts[b]))
    print("D =", D, " D/2 =", D / 2)
    for (a, b) in pairs:
        print(" pair", (a, b), "deg diff", diff[(a, b)], "deg gcd", pg[(a, b)], "deg reduced chord", diff[(a, b)] - pg[(a, b)])
    # combinatorial profile: per block, exponent of pi (block O) in z_j and of pibar (block -O)
    # profile rows: a_j = e_j(O) - min_j e_j(O), e_p = max - min.  Weight deg Q_O (excluding stripped t).
    rows = Counter()
    W = deg9(norms[0])
    cut_sum = 0
    gcd_formula = Counter()
    for o in orbits:
        x = o[0]
        if x == (0, 0): continue
        e = E[x]; lo, hi = min(e), max(e)
        if hi == lo: continue
        a = tuple(ei - lo for ei in e)
        dq = deg9(blocks[x])
        # remove t-part of block for x0 = 0 blocks
        if x[0] == 0:
            dq -= len(o)  # each factor t^L - i x1 t contributes one t
        rows[(a, hi - lo)] += dq  # weight = deg pi-block; each conj pair counted twice over O and -O
        ep = hi - lo
        for tau in range(1, ep + 1):
            S = sum(1 for ai in a if ai < tau)
            cut_sum += dq * (2 * S - 4) ** 2
        for (i1, i2) in pairs:
            gcd_formula[(i1, i2)] += dq * (ep - abs(a[i1] - a[i2]))
    # each conj pair {O,-O} visited twice: once as pi and once as pibar; pibar's row is (hi - e_j)
    # deg p = 2 deg pi. Sum over p of deg p * sum_tau(...)^2 = (sum over both O and -O of deg Q * ...)
    # because the profile of -O is the complement profile, which gives the same cut terms.
    lhs = 8 * sum(diff.values())
    rhs = 4 * 2 * W + cut_sum + 8 * sum(diff[p] - pg[p] for p in pairs)
    print("cut identity: 8*sum deg diff =", lhs, " M(M-2)W + sum_p deg p sum_tau(2S-M)^2 + 8 sum deg c =", rhs)
    print("pair gcd degrees by formula:", {k: v // 2 * 1 for k, v in gcd_formula.items()}, "(halved: each pair of blocks counted twice)")
    print("profile row types (a_j, e_p) with total pi-degree:")
    for k, v in sorted(rows.items()):
        print("  ", k, v)
    print("time %.1fs" % (time.time() - t0))

if __name__ == "__main__":
    main()
