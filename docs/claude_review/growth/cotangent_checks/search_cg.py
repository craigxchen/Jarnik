"""Exhaustive search in the Cilleruelo--Granville five-point regime.

Theorem B: five points with diameter <= K R^(2/5) have intrinsic scale
L_0 <= K^10/128 with 6 | L_0, and in the endpoint chart
N <= (K sqrt((A/L)^2+1)/2)^(10/3).  For each L = 6,12,...,Lmax and each least
cotangent A <= Amax we enumerate every 4-clique (four finite cotangents) whose
partial all-edge lcm stays below that bound, and report the exact diameter
ratio K5 = diam / N^(1/5) of every clique found.
Usage: python3 search_cg.py Kmax Lmax Amax
"""
import math, sys
from cot_lib import edge_norm, lcm, all_edge_lcm, primitive_tuple, gnorm, gsub
from search_fixedL import factor_table, divisors
from itertools import combinations


def run(L, Amax, Kmax):
    fac = factor_table(L, Amax)
    LL = L * L
    out = []
    for A in range(L, Amax + 1):
        T = (Kmax * math.sqrt((A / L) ** 2 + 1) / 2) ** (10 / 3)
        nA = edge_norm(A, L)
        if nA > T:
            continue
        Xmax = int(L * math.sqrt(2 * T)) + 1
        cand = []
        for d in divisors(fac[A], Xmax - A):
            X = A + d
            Q = (A * X + LL) // (A - X)
            part = lcm(lcm(nA, edge_norm(X, L)), edge_norm(Q, L))
            if part <= T:
                cand.append((X, part))
        cand.sort()
        def grow(cl, part, start):
            if len(cl) == 4:
                out.append((cl, part)); return
            for idx in range(start, len(cand)):
                X, pX = cand[idx]
                p = lcm(part, pX)
                if p > T:
                    continue
                ok = True
                for Y in cl[1:]:
                    num = X * Y + LL
                    if num % (X - Y):
                        ok = False; break
                    p = lcm(p, edge_norm(num // (X - Y), L))
                    if p > T:
                        ok = False; break
                if ok:
                    grow(cl + [X], p, idx + 1)
        grow([A], nA, 0)
    res = []
    for cl, N in out:
        rows = primitive_tuple(cl, L)
        assert gnorm(rows[0]) == N
        d = max(math.sqrt(gnorm(gsub(a, b))) for a, b in combinations(rows, 2))
        res.append((d / N ** 0.2, N, tuple(cl)))
    return res


if __name__ == '__main__':
    Kmax = float(sys.argv[1]); Lmax = int(sys.argv[2]); Amax = int(sys.argv[3])
    allres = []
    for L in range(6, Lmax + 1, 6):
        r = run(L, Amax, Kmax)
        allres += [(k, N, L, cl) for k, N, cl in r]
    allres.sort()
    print('Kmax=%.2f L<=%d (multiples of 6) A<=%d: %d cliques' % (Kmax, Lmax, Amax, len(allres)))
    seen = set()
    for k, N, L, cl in allres:
        if N in seen:
            continue
        seen.add(N)
        print('  K5=%.4f N=%d L=%d X=%s' % (k, N, L, cl))
        if len(seen) >= 15:
            break
    print('largest N among cliques found:', max([r[1] for r in allres], default=None))
