import sys
from math import comb
from beta_vanishing import points, hilbert_T
from beta_generic import Tgen
def run(n, N):
    tot = 0; totg = 0; totA = 0
    for s in range(0, N + 1):
        if (s - N) % 2: continue
        pts = points(n, N, s)
        T, _ = hilbert_T(pts, 1000003)
        w = 1 if s == 0 else 2
        tot += w*T; totg += w*Tgen(len(pts), n); totA += w*len(pts)
    print(f"M={n+1} N={N}: h0={totA} beta_N={tot/(N*totA):.4f} generic={totg/(N*totA):.4f} ratio={tot/totg:.4f}", flush=True)
for arg in sys.argv[1:]:
    n, N = map(int, arg.split(','))
    run(n, N)
