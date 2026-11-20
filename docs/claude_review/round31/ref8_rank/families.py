"""Exact rank data for record clusters and families (own code).
For each N, find the largest windows (by M) with arc constant <= CMAX and print
M, omega(N'), affine rank, chord constant C', theorem S (initial segment), K_C = sum S,
and the enlarged affine rank for S = {}, {3}, {3,7}, ... (empirical number of inert primes needed),
plus the 2-adic block (z mod 2).
"""
import sys
from math import log, sqrt, exp
from sympy import factorint
from scan2 import points_of, scan_N, newstats
from gauss import affine_rank, initial_S, indicator_block, PRIMES3


def detail(N, CMAX, label, wantM=None):
    fac = factorint(N)
    rec = []
    stats = newstats()
    scan_N(fac, CMAX, stats, record=rec, Nval=N)
    if not rec:
        print(label, N, 'no windows'); return
    Mmax = max(r[1] for r in rec) if wantM is None else wantM
    best = min((r for r in rec if r[1] == Mmax), key=lambda r: r[5])
    return best


def full_detail(N, CMAX, label, wantM=None):
    """recompute details for the best window of size wantM (or max)"""
    import scan2
    fac = factorint(N)
    found = []
    orig = scan2.analyse_window

    def hook(ps, es, win, stats, CMAX_unused=None, record=None, Nval=None):
        res = orig(ps, es, win, stats, record=record, Nval=Nval)
        found.append((res, ps, es, win))
        return res
    scan2.analyse_window = hook
    try:
        scan_N(fac, CMAX, newstats(), Nval=N)
    finally:
        scan2.analyse_window = orig
    if not found:
        print(label, 'none'); return
    Mmax = max(f[0][0] for f in found) if wantM is None else wantM
    cand = [f for f in found if f[0][0] == Mmax]
    (M, omega, rk, Cp), ps, es, win = min(cand, key=lambda f: f[0][3])
    zs = [w[2] for w in win]
    A = [w[1] for w in win]
    k = len(ps)
    mins = [min(r[t] for r in A) for t in range(k)]
    maxs = [max(r[t] for r in A) for t in range(k)]
    act = [t for t in range(k) if maxs[t] > mins[t]]
    rows = [[r[t] - mins[t] for t in act] for r in A]
    S = initial_S(log(Cp))
    emp = []
    cur = [list(r) for r in rows]
    erk = affine_rank(cur)
    emp.append(erk)
    used = []
    for l in PRIMES3:
        if erk == M - 1:
            break
        blk, ncls = indicator_block(zs, lambda z, l=l: (z[0] % l, z[1] % l))
        cur = [cur[i] + blk[i] for i in range(M)]
        erk = affine_rank(cur)
        used.append(l)
        emp.append(erk)
    b2, n2 = indicator_block(zs, lambda z: (z[0] % 2, z[1] % 2))
    r2 = affine_rank([rows[i] + b2[i] for i in range(M)])
    print(f"{label}: N={N} fac={dict(fac)} M={M} omega(N')={omega} active primes={[ps[t] for t in act]} "
          f"exps={[maxs[t]-mins[t] for t in act]} affrank={rk} C'={Cp:.6f} "
          f"theorem S=first {len(S)} primes 3mod4 (max {S[-1] if S else None}), K_C={sum(S)}; "
          f"enlarged ranks along 3,7,11,...: {emp} (needed {used}); with z mod 2 block: {r2}")
    return M, omega, rk, Cp, len(S), sum(S), used


if __name__ == '__main__':
    pass
