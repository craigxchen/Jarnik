import itertools, sys
from math import comb
from beta_vanishing import points
def Tgen(size, n):
    T = 0; d = 0
    while True:
        h = min(size, comb(d + n, n))
        if h == size: break
        T += size - h; d += 1
    return T
def sizes(n, N):
    out = []
    for s in range(0, N + 1):
        if (s - N) % 2: continue
        out.append((s, len(points(n, N, s))))
    return out
if __name__ == "__main__":
    for n, N, actual in [(2,16,40340),(3,10,45670),(4,6,14631),(5,5,13676)]:
        tot = 0; totA = 0
        for s_, a in sizes(n, N):
            w = 1 if s_ == 0 else 2
            tot += w * Tgen(a, n); totA += w * a
        print(f"M={n+1} N={N}: generic sumT={tot}  actual={actual}  ratio={actual/tot:.4f}  beta_gen={tot/(N*totA):.4f}")
