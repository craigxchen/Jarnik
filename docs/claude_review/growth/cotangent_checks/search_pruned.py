"""Pruned search: integer cotangent cliques at fixed L with small C_*.

For each A (least finite cotangent) we keep only candidates X=A+d, d|A^2+L^2,
and grow cliques while the partial all-edge lcm stays below the threshold
T(A) = (Cmax / (2 arctan(L/A)))^4, i.e. C_* <= Cmax.  Every clique of
m finite coordinates with C_* <= Cmax and least element A <= Amax is found
(the partial lcm only grows when vertices are added).

Usage: python3 search_pruned.py L Amax m Cmax
"""
import math
import sys
from itertools import combinations

from cot_lib import edge_norm, lcm
from search_fixedL import factor_table, divisors


def run(L, Amax, m, Cmax, verbose=True):
    fac = factor_table(L, Amax)
    found = []
    LL = L * L
    for A in range(max(L, 1), Amax + 1):
        T = (Cmax / (2 * math.atan(L / A))) ** 4
        nA = edge_norm(A, L)
        if nA > T:
            continue
        S = A * A + LL
        # every finite X has n_{0X} | N, n_{0X} >= (X^2+L^2)/(2L^2) so X <= L sqrt(2T)
        Xmax = int(L * math.sqrt(2 * T)) + 1
        ds = divisors(fac[A], Xmax - A)
        cand = []
        for d in ds:
            X = A + d
            Q = (A * X + LL) // (A - X)  # exact since d | S
            part = lcm(lcm(nA, edge_norm(X, L)), edge_norm(Q, L))
            if part <= T:
                cand.append((X, part))
        cand.sort()
        # depth-first growth
        def grow(cl, part, start):
            if len(cl) == m:
                found.append((2 * part ** 0.25 * math.atan(L / A), part, tuple(cl)))
                return
            for idx in range(start, len(cand)):
                X, pX = cand[idx]
                ok = True
                p = lcm(part, pX)
                if p > T:
                    continue
                for Y in cl[1:]:
                    num = X * Y + LL
                    if num % (X - Y):
                        ok = False
                        break
                    p = lcm(p, edge_norm(num // (X - Y), L))
                    if p > T:
                        ok = False
                        break
                if ok:
                    grow(cl + [X], p, idx + 1)
        grow([A], nA, 0)
    found.sort()
    if verbose:
        print('L=%d Amax=%d m=%d Cmax=%.2f: %d cliques' % (L, Amax, m, Cmax, len(found)))
        for Cs, N, xs in found[:10]:
            A = xs[0]
            print('  C*=%.4f N=%d expo=%.3f X=%s' % (Cs, N, math.log(N) / math.log(A / L) if A > L else 0, xs))
    return found


if __name__ == '__main__':
    L, Amax, m = (int(a) for a in sys.argv[1:4])
    Cmax = float(sys.argv[4])
    run(L, Amax, m, Cmax)
