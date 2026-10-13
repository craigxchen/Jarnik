"""Boolean separation-lemma test (evidence).

X = {t in {0,1}^M : |<t,u_q> - c_q| <= s for all q}, 0 in X.  Columns col_t = (t^J)_{|J|<m} (multilinear
monomials, 0/1 entries).  phi on X with phi perp Poly_{<m} and phi(0) != 0 exists  <=>  col_0 in span{col_t : t != 0}.
 * char0:   span over Q            (rank test mod a 61-bit prime; rank mod P <= rank over Q)
 * pint:    span over Z_(p)        (phi(0)=1, phi integral => unique p-adic min at 0 since ell(t)=|t|>=1)
 * padic:   span over Z_(p) of p^{1-|t|} col_t  (the exact unique-minimiser criterion)
Module membership over Z_(p) is decided by Smith reduction over Z/p^e (chain ring); valid if e exceeds the largest
elementary-divisor valuation plus the target scaling; we report the max valuation seen as a sanity check.
"""
import itertools, random, sys
import numpy as np

def subsets_upto(M, m):
    out = []
    for j in range(m):
        out.extend(itertools.combinations(range(M), j))
    return out

def build(M, m, U, c, s):
    pts = []
    for bits in itertools.product((0, 1), repeat=M):
        t = np.array(bits)
        if all(abs(int(t @ u) - cq) <= s for u, cq in zip(U, c)):
            pts.append(bits)
    return pts

def colmat(pts, Js):
    G = np.zeros((len(Js), len(pts)), dtype=np.int64)
    for j, J in enumerate(Js):
        for i, t in enumerate(pts):
            if all(t[x] for x in J):
                G[j, i] = 1
    return G

def rank_modP(G, P=2147483647):
    A = G.copy() % P
    r = 0
    rows, cols = A.shape
    for c in range(cols):
        nz = np.nonzero(A[r:, c])[0]
        if len(nz) == 0:
            continue
        piv = r + nz[0]
        A[[r, piv]] = A[[piv, r]]
        inv = pow(int(A[r, c]), P - 2, P)
        A[r] = (A[r] * inv) % P
        others = np.nonzero(A[:, c])[0]
        for i in others:
            if i != r:
                A[i] = (A[i] - A[i, c] * A[r]) % P
        r += 1
        if r == rows:
            break
    return r

def vp(x, p):
    if x == 0:
        return 10**9
    v = 0
    while x % p == 0:
        x //= p; v += 1
    return v

def member_Zp(G, y, p, e):
    """Is y in the Z_(p)-column span of G?  Smith elimination over Z/p^e."""
    Q = p ** e
    A = G.copy() % Q
    y = y.copy() % Q
    rows, cols = A.shape
    r = 0
    maxv = 0
    active_cols = list(range(cols))
    while r < rows:
        sub = A[r:, :]
        if not sub.any():
            break
        # pivot of minimal valuation
        best = None
        vals = sub[sub != 0]
        # find minimal valuation quickly
        for v in range(e):
            mask = (sub % (p ** (v + 1)) != 0)
            if mask.any():
                idx = np.argwhere(mask)[0]
                best = (v, r + idx[0], idx[1])
                break
        if best is None:
            break
        v, pi, pj = best
        maxv = max(maxv, v)
        A[[r, pi]] = A[[pi, r]]
        y[[r, pi]] = y[[pi, r]]
        A[:, [r, pj]] = A[:, [pj, r]]
        # after column swap pivot at (r, r)
        a = int(A[r, r])
        unit = a // (p ** v)
        uinv = pow(unit, -1, Q)
        # eliminate rows below and above in column r
        for i in range(rows):
            if i != r and A[i, r] != 0:
                f = (int(A[i, r]) // (p ** v)) * uinv % Q
                A[i] = (A[i] - f * A[r]) % Q
                y[i] = (y[i] - f * y[r]) % Q
        # eliminate columns right in row r (column ops do not affect y)
        for j in range(cols):
            if j != r and A[r, j] != 0:
                f = (int(A[r, j]) // (p ** v)) * uinv % Q
                A[:, j] = (A[:, j] - f * A[:, r]) % Q
        r += 1
    # now A is diagonal with entries p^{v_i} * unit for i < r, zero rows beyond
    ok = True
    for i in range(rows):
        if i < r:
            v = vp(int(A[i, i]), p)
            if vp(int(y[i]), p) < v:
                ok = False
        else:
            if y[i] % Q != 0:
                ok = False
    return ok, maxv, r

def test(M, m, U, c, s, p=3, e=19, padic=True):
    pts = build(M, m, U, c, s)
    zero = tuple([0] * M)
    assert zero in pts
    pts.remove(zero)
    pts = [zero] + pts
    Js = subsets_upto(M, m)
    G = colmat(pts, Js)
    E = len(Js)
    r_all = rank_modP(G)
    r_oth = rank_modP(G[:, 1:])
    char0 = (r_all == r_oth)          # col_0 in span of others (mod P)
    res = dict(M=M, m=m, s=s, r=len(U), X=len(pts), E=E, rank=r_all, char0=char0)
    if char0:
        ok, mv, rk = member_Zp(G[:, 1:], G[:, 0], p, e)
        res['pint'] = ok; res['pint_maxv'] = mv
        if padic and not ok:
            K = max(sum(t) for t in pts[1:])
            Gs = G[:, 1:].copy()
            for i, t in enumerate(pts[1:]):
                Gs[:, i] *= p ** (K - (sum(t) - 1))
            ok2, mv2, _ = member_Zp(Gs, G[:, 0] * p ** K, p, e)
            res['padic'] = ok2; res['padic_maxv'] = mv2; res['K'] = K
    return res

if __name__ == '__main__':
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
    trials = [(8, 3, 1, 3), (8, 4, 1, 3), (9, 4, 1, 4), (10, 4, 1, 4), (10, 5, 1, 4), (10, 4, 1, 6), (11, 4, 1, 5), (12, 4, 1, 5), (12, 5, 1, 5)]
    for (M, m, s, r) in trials:
        for rep in range(2):
            U = []
            for q in range(r):
                u = [rng.choice([-1, 0, 1]) for _ in range(M)]
                U.append(np.array(u))
            c = [0] * r
            res = test(M, m, U, c, s)
            print(res); sys.stdout.flush()
