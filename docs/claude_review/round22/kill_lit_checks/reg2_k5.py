import sys
sys.path.insert(0, '.')
from reg2 import lattice_points, rank_mod
k = 5
for m, start in [(1,0),(2,0),(3,0),(4,8),(5,13)]:
    S = lattice_points(k, m)
    out = []
    for p in (1000003, 998244353):
        d = start; first = True
        while True:
            rk, nE = rank_mod(S, m, d, p, coords=[1,2,3,4])
            if rk == len(S):
                out.append(d if not first or d == 0 else f"<={d}")
                break
            first = False; d += 1
    print(f"k=5 m={m} |S|={len(S)} max ord (mod p, two primes) = {out}", flush=True)
