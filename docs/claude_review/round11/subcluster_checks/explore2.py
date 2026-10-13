# Exact: find small supercritical models (m = floor(LP)+1, LP = avg_p h_A(beta_p)) with dim V(A,m) > 0.
import sys, random, math, time
from engine import choose_m
random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
n = int(sys.argv[2]); r = int(sys.argv[3])
cols = []
while len(cols) < r:
    a = tuple(random.randint(0, 1) for _ in range(n))
    if 0 < sum(a) < n and a not in cols and tuple(1-x for x in a) not in cols: cols.append(a)
print("cols", cols)
for T in [int(t) for t in sys.argv[4].split(',')]:
    for w in [int(t) for t in sys.argv[5].split(',')]:
        t0 = time.time()
        m, U, dV, lp, mh = choose_m(list(range(n)), n-1, cols, T, w, [1/r]*r, 1)
        print(f"n={n} r={r} T={T} w={w}: |A|={U} box={(2*T+1)**(n-1)} LP={lp:.2f} maxh={mh} m={m} dimV={dV} binom={math.comb(m+n-2,n-1)} ({time.time()-t0:.0f}s)", flush=True)
