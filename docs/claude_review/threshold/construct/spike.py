"""Spike polytope K_M = conv{ +-(M e_i - 1) } = {x in Z^M : sum x = 0, sum_i |x_i - med(x)| <= M}.
Count lattice points of nK, its cut widths, the symmetric adversary, and (if small) reg."""
import sys, itertools, time
from fractions import Fraction as Fr
from clib import reg2, proj, widths, sym_adversary, popc

def medl1(x):
    xs = sorted(x); M = len(x)
    c = xs[(M - 1) // 2]
    return sum(abs(v - c) for v in xs)

def nK(M, n):
    R = n * M
    out = []
    def rec(pref):
        if len(pref) == M - 1:
            last = -sum(pref)
            x = pref + [last]
            if medl1(x) <= R: out.append(tuple(x))
            return
        for a in range(-R, R + 1):
            p2 = pref + [a]
            # prune: partial median-l1 lower bound is hard; use |a| <= R
            rec(p2)
    # smarter: enumerate by bounding each coordinate in [-(M-1)n, (M-1)n]
    B = (M - 1) * n
    for pref in itertools.product(range(-B, B + 1), repeat=M - 1):
        last = -sum(pref)
        if abs(last) > B: continue
        x = pref + (last,)
        if medl1(x) <= R: out.append(x)
    return out

if __name__ == "__main__":
    M = int(sys.argv[1]); n = int(sys.argv[2]); doreg = len(sys.argv) > 3
    t0 = time.time()
    S = nK(M, n)
    w = widths(S, M)
    ws = {}
    for P, v in w.items():
        ws.setdefault(min(popc(P), M - popc(P)), set()).add(v)
    ws = {j: v.pop() for j, v in ws.items()}
    adv, mu = sym_adversary(ws, M)
    print(f"M={M} n={n}: |nK|={len(S)} widths={ws} adv={adv} 2adv={2*adv} ({time.time()-t0:.1f}s)", flush=True)
    if doreg:
        r = reg2(proj(S))
        print(f"   reg={r}  ratio reg/(2adv) = {float(Fr(r)/(2*adv)):.4f} ({time.time()-t0:.1f}s)")
