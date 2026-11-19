"""Fast exact family minimum for deg T = B >= d (all factors of degree B).
g = hB + res, res = sum_{f small} deg f * (min_i ord_f(F_i) - ord_f(prod_C)),  f over monic irreducibles of
F_4[t] of degree <= (d-1)//2 (larger f cannot divide the residual gcd: proof of (V.4), level j=1 is a singleton).
req = 2*max_pairs deg P_J(1) (exact when deg T >= d).  Cross-checked against family_min.py on deg T <= 12.
usage: family_fast.py d maxdegT"""
import sys, itertools
from f4poly import *


def monic_irreducibles(k):
    out = []
    for coeffs in itertools.product(range(4), repeat=k):
        f = tuple(coeffs) + (1,)
        if k == 1:
            out.append(f); continue
        ok = True
        for kk in range(1, k // 2 + 1):
            for c2 in itertools.product(range(4), repeat=kk):
                if not divmod_(f, tuple(c2) + (1,))[1]:
                    ok = False; break
            if not ok:
                break
        if ok:
            out.append(f)
    return out


def ordf(p, f):
    n = 0
    while True:
        q, r = divmod_(p, f)
        if r:
            return n
        p, n = q, n + 1


def run(d, T, best, smallf):
    Q = 1 << d
    B = f2deg(T)
    assert B >= d
    fs = forms_and_sides(d)
    lin = [linfac(T, v) for v in range(Q)]
    ordv = [[ordf(lin[v], f) for f in smallf] for v in range(Q)]
    degf = [deg(f) for f in smallf]
    PJ1 = {}
    for a, b in itertools.combinations(range(len(fs)), 2):
        PJ1[(a, b)] = sum(f2deg(1 ^ v) for v in fs[a][1] & fs[b][1])
    sides = []
    for fi, (w, k, c) in enumerate(fs):
        for S in (k, c):
            vec = tuple(sum(ordv[v][x] for v in S) for x in range(len(smallf)))
            sides.append((fi, S, vec))
    D0 = Q * B // 2

    def gval(mv, C):
        cv = [sum(ordv[v][x] for v in C) for x in range(len(smallf))]
        return len(C) * B + sum(df * (m - c) for df, m, c in zip(degf, mv, cv))

    def visit(fi, chosen, forms, mv, C, r):
        if fi == len(fs):
            if len(chosen) >= 2:
                pd = D0 - gval(mv, C)
                s = len(chosen)
                if s not in best or pd < best[s][0]:
                    best[s] = (pd, T, tuple(chosen), r, len(C))
            return
        visit(fi + 1, chosen, forms, mv, C, r)
        for side in (2 * fi, 2 * fi + 1):
            nmv = tuple(min(a, b) for a, b in zip(mv, sides[side][2])) if chosen else sides[side][2]
            nC = C & sides[side][1]
            nr = max([r] + [2 * PJ1[(f0, fi)] for f0 in forms])
            if chosen and nr > gval(nmv, nC):
                continue
            visit(fi + 1, chosen + [side], forms + [fi], nmv, nC, nr)

    visit(0, [], [], (), frozenset(range(Q)), -10**9)


if __name__ == "__main__":
    d, maxdeg = int(sys.argv[1]), int(sys.argv[2])
    Q = 1 << d
    smallf = [f for k in range(1, (d - 1) // 2 + 1) for f in monic_irreducibles(k)]
    best = {}
    n = 0
    for B in range(d, maxdeg + 1):
        for low in range(1 << (B - d)):
            T = (1 << B) | (low << d)
            run(d, T, best, smallf)
            n += 1
        print(f"after deg T <= {B}: " + ", ".join(f"s={s}:{best[s][0]}" for s in sorted(best)), flush=True)
    print(f"d={d}: {n} classes of T mod V with d <= deg T <= {maxdeg}")
    for s in sorted(best):
        pd, T, ch, r, h = best[s]
        print(f"s={s:2d}: min point degree {pd:5d} at deg T={f2deg(T)}, req={r}, h={h}, sides={ch}")
