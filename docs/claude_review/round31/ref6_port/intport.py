"""Integer analogues of the ternary construction VI.4 (own code, exact Gaussian integers).

Points Z_j = prod_{x in T} (X - i x)^{e_j(x)}, e_j(x) + e_j(-x) = 2, T = -T, so all Z_j have the norm
X^{e(0)} * prod_{x>0} (X^2+x^2)^2 (times X^2 if 0 in T).  For each family:
  (a) template regime: deg_X(Z_i(X) - Z_j(X)) as polynomials in Z[i][X], against the formal common degree g;
  (b) integer regime: exact Gaussian integers at chosen X, exact Gaussian gcd of all points, width exponent
      w_ij = log|Z'_i - Z'_j| / log|Z'_i|  (a c=0-type cluster needs w <= 1/2 + O(1/log R)),
  (c) best X in a range by the angle formula theta_ij = 2 sum_{x>0} d_ij(x) arctan(x/X).
"""
import itertools, math, sys
import numpy as np

def legendre(a, p):
    a %= p
    if a == 0: return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1

# ---- Gaussian integers as (re, im) python ints ----
def gmul(a, b): return (a[0]*b[0] - a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def gnorm(a): return a[0]*a[0] + a[1]*a[1]
def gdivmod(a, b):
    n = gnorm(b)
    num = gmul(a, (b[0], -b[1]))
    # nearest integer quotient
    qr = (2*num[0] + n) // (2*n); qi = (2*num[1] + n) // (2*n)
    q = (qr, qi)
    qb = gmul(q, b)
    return q, (a[0]-qb[0], a[1]-qb[1])
def ggcd(a, b):
    while b != (0, 0):
        _, r = gdivmod(a, b)
        a, b = b, r
    return a
def gexact(a, b):
    q, r = gdivmod(a, b)
    assert r == (0, 0)
    return q

def point(X, T, e):
    z = (1, 0)
    for x in T:
        for _ in range(e[x]):
            z = gmul(z, (X, -x))
    return z

def polypoint(T, e):
    # polynomial in X with Gaussian-integer coefficients, list low->high of (re,im)
    p = [(1, 0)]
    for x in T:
        for _ in range(e[x]):
            q = [(0, 0)] * (len(p) + 1)
            for k, c in enumerate(p):
                q[k+1] = (q[k+1][0] + c[0], q[k+1][1] + c[1])
                m = gmul(c, (0, -x))
                q[k] = (q[k][0] + m[0], q[k][1] + m[1])
            p = q
    return p

def polydiffdeg(p, q):
    for k in range(len(p) - 1, -1, -1):
        if p[k] != q[k]: return k
    return -1

def log_abs(z):
    # log|z| for huge Gaussian ints
    n = gnorm(z)
    if n == 0: return -math.inf
    return 0.5 * math.log(n)

def analyse(name, T, pats, Xs, Xrange=None):
    print("=" * 78); print(name)
    n = sum(pats[0][x] for x in T)
    gform = sum(min(p[x] for p in pats) for x in T)
    m = len(pats)
    # distinct?
    keys = [tuple(p[x] for x in T) for p in pats]
    assert len(set(keys)) == m, "patterns not distinct"
    for p in pats:
        assert all(p[x] + p[-x] == 2 for x in T)
    print(" m=%d points, X-degree n=%d, formal common degree g=%d" % (m, n, gform))
    polys = [polypoint(T, p) for p in pats]
    worst = max(polydiffdeg(polys[i], polys[j]) for i, j in itertools.combinations(range(m), 2))
    best = min(polydiffdeg(polys[i], polys[j]) for i, j in itertools.combinations(range(m), 2))
    print(" template regime: pair difference X-degrees in [%d, %d]; after formal division: point degree %d, "
          "max diff %d, needed <= %.1f; template width exponent %.3f"
          % (best, worst, n - gform, worst - gform, (n - gform) / 2, (worst - gform) / (n - gform)))
    for X in Xs:
        Z = [point(X, T, p) for p in pats]
        for z in Z[1:]:
            assert gnorm(z) == gnorm(Z[0])
        G = Z[0]
        for z in Z[1:]: G = ggcd(G, z)
        Zr = [gexact(z, G) for z in Z]
        lR = log_abs(Zr[0])
        ws = []
        for i, j in itertools.combinations(range(m), 2):
            d = (Zr[i][0] - Zr[j][0], Zr[i][1] - Zr[j][1])
            ws.append(log_abs(d) / lR)
        lG = log_abs(G); lform = sum(min(p[x] for p in pats) * 0.5 * math.log(X*X + x*x) for x in T)
        print(" X=%d: log|R'|=%.1f, log|gcd|=%.2f (formal %.2f), width exponents min %.4f max %.4f"
              % (X, lR, lG, lform, min(ws), max(ws)))
    if Xrange:
        Xa = np.arange(1, Xrange + 1, dtype=np.float64)
        pos = [x for x in T if x > 0]
        At = {x: np.arctan(x / Xa) for x in pos}
        logR = sum(np.log(Xa*Xa + x*x) for x in pos) + (np.log(Xa) if 0 in T else 0)
        logG = sum(min(p[x] for p in pats) * 0.5 * np.log(Xa*Xa + x*x) for x in T)
        lRp = logR - logG
        W = None
        for i, j in itertools.combinations(range(m), 2):
            th = 2 * sum((pats[i][x] - pats[j][x]) * At[x] for x in pos)
            s = np.abs(2 * np.sin(th / 2))
            s = np.maximum(s, 1e-300)
            w = 1 + np.log(s) / lRp
            W = w if W is None else np.maximum(W, w)
        k = int(np.argmin(W))
        print(" best X in [1,%d] (formal gcd): max-pair width exponent %.4f at X=%d (needed ~0.5); median %.4f"
              % (Xrange, W[k], k + 1, float(np.median(W))))
        return k + 1
    return None

def main():
    B = 15
    T = list(range(-B, B + 1))
    # E1: one odd quadratic character, all multipliers lambda
    ell = 7
    pats = {}
    for lam in range(1, ell):
        e = {x: 1 + legendre(lam * x, ell) for x in T}
        pats.setdefault(tuple(e[x] for x in T), lam)
    print("E1: chi mod 7, multipliers 1..6 give", len(pats), "distinct patterns")
    e1 = {x: 1 + legendre(x, ell) for x in T}; e2 = {x: 1 - legendre(x, ell) for x in T}
    X = 1000
    Z1, Z2 = point(X, T, e1), point(X, T, e2)
    print("   Z_- == conj(Z_+):", Z2 == (Z1[0], -Z1[1]))
    # E2: several odd real characters of different conductors
    chars = [lambda x: legendre(x, 3), lambda x: legendre(x, 7), lambda x: legendre(x, 11),
             lambda x: legendre(x, 19), lambda x: legendre(x, 3) * legendre(x, 5)]
    for mm in (2, 3, 4, 5):
        pats = [{x: 1 + ch(x) for x in T} for ch in chars[:mm]]
        best = analyse("E2: odd real characters (conductors 3,7,11,19,15)[:%d], T=[-%d,%d]" % (mm, B, B),
                       T, pats, [10, 1000, 10**6], Xrange=200000)
    # E3: literal port: balanced base-ell digit vectors, F_ell-linear forms, Legendre symbol
    for ell, k, forms in ((3, 3, [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0)]),
                          (7, 2, [(1, 0), (0, 1), (1, 1), (1, 6)]),
                          (11, 2, [(1, 0), (0, 1), (1, 1), (1, 10), (1, 2)])):
        half = (ell**k - 1) // 2
        Tk = list(range(-half, half + 1))
        def digits(x):
            d = []
            for _ in range(k):
                r = ((x + (ell - 1) // 2) % ell) - (ell - 1) // 2
                d.append(r); x = (x - r) // ell
            assert x == 0
            return d
        dig = {x: digits(x) for x in Tk}
        for x in Tk:
            assert [(-a) for a in dig[x]] == dig[-x]
        bestset = None
        for flips in itertools.product((1, -1), repeat=len(forms)):
            pats = [{x: 1 + legendre(fl * sum(a*b for a, b in zip(f, dig[x])), ell) for x in Tk}
                    for f, fl in zip(forms, flips)]
            g = sum(min(p[x] for p in pats) for x in Tk)
            if bestset is None or g > bestset[0]: bestset = (g, pats, flips)
        g, pats, flips = bestset
        analyse("E3: digit-vector port, ell=%d, k=%d digits, %d forms, flips %s, |T|=%d"
                % (ell, k, len(forms), flips, len(Tk)), Tk, pats, [10, 1000, 10**6], Xrange=200000)

if __name__ == "__main__":
    main()
