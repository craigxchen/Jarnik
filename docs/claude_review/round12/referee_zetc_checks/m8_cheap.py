# ZE-TC^cheap in the zonotope-free model M8 = Syl(8,2), m = 8.
import itertools, numpy as np, time
from rlib import Exact
pts = list(range(1, 8))
cols = [tuple(i for i, x in enumerate(pts) if bin(a & x).count("1") % 2) for a in range(1, 8)]
A = [k for k in itertools.product((-1,0,1), repeat=7) if all(abs(sum(k[i] for i in S)) <= 2 for S in cols)]
print("|A|", len(A), "columns", cols, flush=True)
# every coordinate +-e_x in A (so every nonzero integer direction has width >= 2)
print("+-e_x in A:", all(tuple((s if j == x else 0) for j in range(7)) in set(A) for x in range(7) for s in (1, -1)))
t = time.time(); E = Exact(A); print("rho", E.rho, "H", E.H, f"{time.time()-t:.0f}s", flush=True)
S = cols[0]
for u_S in [(1, 2, 3, 5), (2, -7, 3, 11), (1, 3, 9, 27)]:
    u = [0]*7
    for i, x in enumerate(S): u[x] = u_S[i]
    print("u supported on column", S, "u =", u, "nu =", E.nu(u), flush=True)
# a generic point: nu should be rho = 8
print("generic u nu =", E.nu([1, 3, 7, 15, 31, 63, 127]))
