# Exact parameter exploration for small supercritical models (evidence infrastructure).
# Criterion (LP of Thm 4.1, leading order): V(A,m) != 0 and sum_p omega_p h_A(beta_p) < m, omega_p = 1/r.
import sys, random, itertools
from common import *
random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
def run(n, r, Ts, ws):
    pts = list(range(n))
    cols = []
    while len(cols) < r:
        a = tuple(random.randint(0, 1) for _ in range(n))
        if 0 < sum(a) < n: cols.append(a)
    base = n-1
    for T in Ts:
        for w in ws:
            A = slab_box(pts, base, cols, T, w)
            box = (2*T+1)**(n-1)
            hs = [width(A, [2*x-1 for x in a]) for a in cols]
            avg = sum(hs)/r
            # largest m with rank Poly_{<m}|_A < |A|
            best = None
            for m in range(1, 40):
                if math.comb(m+n-2, n-1) >= len(A) + 200: break
                Mat, mons = moment_matrix(A, base, m)
                rk = rank_mod_p(Mat)
                if rk < len(A) and m > avg:
                    best = (m, rk)
            print(f"n={n} r={r} T={T} w={w}: |A|={len(A)} (box {box}) avg h={avg:.2f} max h={max(hs)} -> "
                  f"supercritical m: {best}", flush=True)
for n, r in [(4, 4), (4, 8), (5, 5), (5, 10)]:
    run(n, r, Ts=[2, 3, 4], ws=[4, 6, 8, 10])
