"""Exact scan over X (exact Gaussian gcd at every X) for two integer ports, and an exhaustive
vertical-line sign-pattern (odd PTE) template count for small D."""
import itertools, math, sys
sys.path.insert(0, '.')
from intport import legendre, point, ggcd, gexact, gnorm, log_abs
import numpy as np

def exact_scan(name, T, pats, Xmax):
    m = len(pats)
    best = (9, None, None)
    for X in range(1, Xmax + 1):
        Z = [point(X, T, p) for p in pats]
        G = Z[0]
        for z in Z[1:]: G = ggcd(G, z)
        Zr = [gexact(z, G) for z in Z]
        lR = log_abs(Zr[0])
        if lR <= 0: continue
        w = max(log_abs((Zr[i][0]-Zr[j][0], Zr[i][1]-Zr[j][1])) / lR for i, j in itertools.combinations(range(m), 2))
        if w < best[0]: best = (w, X, lR)
    print("%s: exact scan X<=%d: min over X of max-pair width exponent = %.4f at X=%d (log R'=%.1f)"
          % (name, Xmax, best[0], best[1], best[2]))

def main():
    B = 15; T = list(range(-B, B + 1))
    chars = [lambda x: legendre(x, 3), lambda x: legendre(x, 7), lambda x: legendre(x, 11), lambda x: legendre(x, 19)]
    for mm in (2, 3, 4):
        pats = [{x: 1 + ch(x) for x in T} for ch in chars[:mm]]
        exact_scan("E2 m=%d" % mm, T, pats, 3000)
    ell, k = 7, 2
    half = (ell**k - 1)//2; Tk = list(range(-half, half+1))
    def digits(x):
        d = []
        for _ in range(k):
            r = ((x + 3) % 7) - 3; d.append(r); x = (x - r)//7
        return d
    forms = [(1, 0), (0, 1), (1, 1), (1, 6)]
    pats = [{x: 1 + legendre(sum(a*b for a, b in zip(f, digits(x))), ell) for x in Tk} for f in forms]
    exact_scan("E3 ell=7 k=2 m=4", Tk, pats, 3000)
    # vertical-line sign-pattern templates: roots +-i x, x in S subset {1..Xb}, |S| = D,
    # points eps in {+-1}^D; c=0 cluster iff sum (eps-eps') x^m = 0 for odd m <= L, L = D-1-floor(D/2)
    for D, Xb in ((6, 12), (7, 12), (8, 13), (9, 14), (10, 14), (11, 15), (12, 16)):
        L = D - 1 - D // 2
        odd = [m for m in range(1, L + 1) if m % 2]
        E = np.array(list(itertools.product((1, -1), repeat=D)), dtype=np.int64)
        best = (0, None)
        for S in itertools.combinations(range(1, Xb + 1), D):
            Mx = np.array([[x**m for x in S] for m in odd], dtype=np.int64).T  # D x |odd|
            vals = E @ Mx
            _, cnt = np.unique(vals, axis=0, return_counts=True)
            c = int(cnt.max())
            if c > best[0]: best = (c, S)
        print("vertical-line template D=%d (odd moments m<=%d), S in [1,%d]: max cluster %d, e.g. S=%s"
              % (D, L, Xb, best[0], best[1]))

if __name__ == "__main__":
    main()
