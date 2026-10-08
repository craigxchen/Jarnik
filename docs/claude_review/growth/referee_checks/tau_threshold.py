"""Pseudo-effective threshold tau(H,E) on Bl_Y X_M, numerically.
max order of vanishing along Y of a nonzero section of O(N) = max_s dstar(A_s),
dstar(A) = min{d : polynomials of degree <= d interpolate on A}  (= len(hist)-1).
Also the twisted Ru-Vojta constant G(b) = b + int_b^tau v / v(b), computed from the
finite-N Hilbert data: sum_{m>=bN} h0(NH - mE) / (N h0(NH - bN E)) + b."""
import sys
from beta_vanishing import points, hilbert_T

def data(n, N, p=1000003):
    # returns dict m -> h0(NH - mE) for m>=0
    h0 = {}
    maxd = 0
    per_s = []
    for s in range(0, N + 1):
        if (s - N) % 2: continue
        pts = points(n, N, s)
        T, hist = hilbert_T(pts, p)
        w = 1 if s == 0 else 2
        dstar = len(hist) - 1
        maxd = max(maxd, dstar)
        per_s.append((s, len(pts), dstar))
        size = len(pts)
        # h0 contribution for order >= m: size - h_A(m-1); h_A(-1)=0
        for m in range(0, dstar + 2):
            hA = 0 if m == 0 else (hist[m-1] if m-1 < len(hist) else size)
            h0[m] = h0.get(m, 0) + w * (size - hA)
    return h0, maxd, per_s

if __name__ == '__main__':
    for arg in sys.argv[1:]:
        n, N = map(int, arg.split(','))
        h0, maxd, per_s = data(n, N)
        h0b, maxdb, _ = data(n, N, 999983)
        assert h0 == h0b and maxd == maxdb
        print(f"M={n+1} N={N}: max vanishing order along Y = {maxd}, ratio tau_N = {maxd/N:.4f}")
        print("   per s (s,|A_s|,dstar):", per_s)
        best = None
        for bm in range(0, maxd + 1):
            if h0.get(bm, 0) == 0: break
            tail = sum(h0.get(m, 0) for m in range(bm + 1, maxd + 2))
            G = bm / N + tail / (N * h0[bm])
            print(f"   b={bm}/{N}: h0(NH-bNE)={h0[bm]}, beta(L_b,E)~{tail/(N*h0[bm]):.4f}, G(b)=b+beta={G:.4f}")
        sys.stdout.flush()
