"""Exact-modular Hilbert functions of finite point sets in Z^n.

h_A(d) = dim of {polynomials of degree <= d} restricted to A  (rank of the evaluation matrix).
reg(A) = min{d : h_A(d) = |A|}  = max order of vanishing (at w=1) of a nonzero Laurent
polynomial supported on A (duality: c annihilates Poly_{<m} on A).

Rank mod a prime q is <= rank over Q, with equality for all but finitely many q.  So
  * h_A(d) mod q = |A|  ==> h_A(d) over Q = |A|  (rigorous upper bound for reg),
  * a lower bound reg(A) >= m needs a rational certificate (see certify.py).
"""
import numpy as np, itertools

Q1, Q2 = 1000003, 999983

def monos(n, d):
    if n == 0:
        if d == 0:
            yield ()
        return
    if n == 1:
        yield (d,)
        return
    for a in range(d, -1, -1):
        for rest in monos(n - 1, d - a):
            yield (a,) + rest

def rref_rows(V, q):
    """RREF of V mod q (rows). Returns (R, pivots)."""
    V = V.copy()
    r = 0
    piv = []
    rows, cols = V.shape
    for c in range(cols):
        if r == rows:
            break
        nz = np.nonzero(V[r:, c])[0]
        if len(nz) == 0:
            continue
        i = r + nz[0]
        if i != r:
            V[[r, i]] = V[[i, r]]
        inv = pow(int(V[r, c]), q - 2, q)
        V[r] = (V[r] * inv) % q
        col = V[:, c].copy()
        col[r] = 0
        nzr = np.nonzero(col)[0]
        if len(nzr):
            V[nzr] = (V[nzr] - (col[nzr, None] * V[r][None, :]) % q) % q
        piv.append(c)
        r += 1
    return V[:r], piv

def hilbert(pts, q=Q1, dmax=None):
    """Return list hist with hist[d] = h_A(d) mod q, stopping when h = |A|."""
    A = np.array(pts, dtype=np.int64)
    size, n = A.shape
    Aq = A % q
    B = np.zeros((0, size), dtype=np.int64)
    piv = []
    hist = []
    # power table
    maxdeg = size  # upper bound
    pw = {}
    def power(j, a):
        if (j, a) not in pw:
            pw[(j, a)] = np.ones(size, dtype=np.int64) if a == 0 else (power(j, a - 1) * Aq[:, j]) % q
        return pw[(j, a)]
    d = 0
    while True:
        rowsV = []
        for m in monos(n, d):
            v = np.ones(size, dtype=np.int64)
            for j, a in enumerate(m):
                if a:
                    v = (v * power(j, a)) % q
            rowsV.append(v)
        V = np.array(rowsV, dtype=np.int64)
        if len(piv):
            coeff = V[:, piv]
            V = (V - (coeff @ B) % q) % q
        R, npiv_local = rref_rows(V, q)
        if len(npiv_local):
            # eliminate new pivot columns from B
            if len(piv):
                cB = B[:, npiv_local]
                B = (B - (cB @ R) % q) % q
            B = np.vstack([B, R])
            piv = piv + list(npiv_local)
        hist.append(len(piv))
        if len(piv) == size:
            return hist
        d += 1
        if dmax is not None and d > dmax:
            return hist

def reg(pts, q=Q1):
    return len(hilbert(pts, q)) - 1

# ---------- point sets ----------
def shell(M, p, n):
    """S^M(p,n) = {k in Z^M : |k^+|_1 = p, |k^-|_1 = n}; returned in first M-1 coords."""
    out = []
    s = p - n
    def rec(prefix, P, Nn):
        j = len(prefix)
        if j == M - 1:
            last = s - sum(prefix)
            pp = P + max(last, 0)
            nn = Nn + max(-last, 0)
            if pp == p and nn == n:
                out.append(tuple(prefix))
            return
        for a in range(-(n - Nn), p - P + 1):
            rec(prefix + [a], P + max(a, 0), Nn + max(-a, 0))
    rec([], 0, 0)
    return out

def ball(M, p, n):
    """Q^M(p,n) = {k : sum k = p-n, |k|_1 <= p+n} = union of shells S(p-j,n-j)."""
    out = []
    for j in range(0, min(p, n) + 1):
        out += shell(M, p - j, n - j)
    return out

if __name__ == '__main__':
    import sys
    M = int(sys.argv[1]); kind = sys.argv[2]; tmax = int(sys.argv[3])
    f = shell if kind == 'shell' else ball
    for t in range(1, tmax + 1):
        for n in range(0, t // 2 + 1):
            p = t - n
            pts = f(M, p, n)
            h1 = hilbert(pts, Q1); h2 = hilbert(pts, Q2)
            assert h1 == h2, (M, p, n)
            r = len(h1) - 1
            print(f"M={M} {kind} p={p} n={n} t={t} |A|={len(pts)} reg={r} ratio={r/t:.4f}", flush=True)
