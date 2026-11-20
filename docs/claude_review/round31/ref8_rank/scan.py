"""Test of Part VIII.1 on genuine lattice-point clusters (own code).

For each norm N (prime factors all = 1 mod 4) enumerate all lattice points, take every window of
consecutive points (start angle in [0, pi/2)) with arc constant C = arc/N^(1/4) <= CMAX, normalise by
the common Gaussian gcd and record:
  M, omega(N'), affine rank of exponent rows, chord constant C' = max chord / N'^(1/4),
  enlarged affine rank with the residue indicators of the theorem (S valid for C'), etc.
Assertions encode the theorem: enlarged rows affinely independent whenever S is valid.

usage: python3 scan.py mode args...
  mode 'all X CMAX'   : every N <= X with all prime factors = 1 mod 4
  mode 'list CMAX N1 N2 ...'
"""
import sys
from math import atan2, log, sqrt, pi
from itertools import product
from collections import defaultdict
from sympy import factorint, primerange
from gauss import (gmul, gconj, gnorm, gpow, gauss_prime, affine_rank, initial_S,
                   indicator_block, PRIMES3)

LOG2 = log(2.0)


def points_of(fac):
    ps = sorted(fac)
    gp = [gauss_prime(p) for p in ps]
    es = [fac[p] for p in ps]
    pw = [[gpow(g, k) for k in range(e + 1)] for g, e in zip(gp, es)]
    pwc = [[gpow(gconj(g), k) for k in range(e + 1)] for g, e in zip(gp, es)]
    pts = []
    for a in product(*[range(e + 1) for e in es]):
        z = (1, 0)
        for t, (ai, e) in enumerate(zip(a, es)):
            z = gmul(z, gmul(pw[t][ai], pwc[t][e - ai]))
        for u in [(1, 0), (0, 1), (-1, 0), (0, -1)]:
            pts.append((gmul(z, u), a))
    return ps, es, pts


def analyse_window(ps, es, win, stats, CMAX_unused=None, record=None, Nval=None):
    M = len(win)
    zs = [w[0] for w in win]
    A = [w[1] for w in win]
    k = len(ps)
    mins = [min(r[t] for r in A) for t in range(k)]
    maxs = [max(r[t] for r in A) for t in range(k)]
    act = [t for t in range(k) if maxs[t] > mins[t]]
    rows = [[r[t] - mins[t] for t in act] for r in A]
    epr = [maxs[t] - mins[t] for t in act]
    omega = len(act)
    Wp = sum(epr[i] * log(ps[act[i]]) for i in range(omega))
    rk = affine_rank(rows)
    maxch = 0.0
    for i in range(M):
        for j in range(i + 1, M):
            d = sqrt(gnorm((zs[i][0] - zs[j][0], zs[i][1] - zs[j][1])))
            maxch = max(maxch, d)
    # normalised chord constant: chords divided by |g| = sqrt(N/N'), radius scale N'^(1/4)
    logC = log(maxch) - 0.5 * log(Nval) + 0.25 * Wp
    from math import exp
    Cp = exp(logC)
    # theorem as stated: S valid for C'
    S = initial_S(logC)
    blocks = {}
    for l in S:
        blocks[l] = indicator_block(zs, lambda z, l=l: (z[0] % l, z[1] % l))
    enl = [rows[i] + sum((blocks[l][0][i] for l in S), []) for i in range(M)]
    erk = affine_rank(enl)
    assert erk == M - 1, ("THEOREM FAILS", Nval, M, Cp, S)
    blockdim = sum(min(l, blocks[l][1] - 1) for l in S)
    assert M - 1 - rk <= blockdim
    assert M <= 1 + omega + sum(S)
    # improved variant: free factor (1+i) | b  => extra (1/2)log 2
    S2 = initial_S(logC, extra=0.5 * LOG2)
    enl2 = [rows[i] + sum((indicator_block(zs, lambda z, l=l: (z[0] % l, z[1] % l))[0][i] for l in S2), [])
            for i in range(M)]
    assert affine_rank(enl2) == M - 1, ("VARIANT FAILS", Nval, M, Cp, S2)
    if Cp < sqrt(2):
        assert rk == M - 1, ("C'<sqrt2 but rank deficient", Nval)
    # negative control: one prime short of the stated condition
    neg = None
    if len(S) >= 1:
        Sm = S[:-1]
        enlm = [rows[i] + sum((indicator_block(zs, lambda z, l=l: (z[0] % l, z[1] % l))[0][i] for l in Sm), [])
                for i in range(M)]
        neg = affine_rank(enlm) < M - 1
    key = M
    st = stats[key]
    st['n'] += 1
    st['max_def_omega'] = max(st['max_def_omega'], M - 1 - omega)
    st['max_def_rank'] = max(st['max_def_rank'], M - 1 - rk)
    if M - 1 - rk > 0:
        if Cp < st['minC_rankdef']:
            st['minC_rankdef'] = Cp
            st['ex_rankdef'] = (Nval, M, rk, omega, round(Cp, 6))
    if M - 1 - omega > 0:
        if Cp < st['minC_omdef']:
            st['minC_omdef'] = Cp
            st['ex_omdef'] = (Nval, M, rk, omega, round(Cp, 6))
    if neg:
        st['neg_dep'] += 1
        if Cp < st['minC_negdep']:
            st['minC_negdep'] = Cp
            st['ex_negdep'] = (Nval, M, rk, omega, round(Cp, 6), S[:-1])
    if record is not None:
        record.append((Nval, M, omega, rk, erk, Cp, S, blockdim))
    return M, omega, rk, Cp


def newstats():
    return defaultdict(lambda: {'n': 0, 'max_def_omega': -99, 'max_def_rank': -99,
                                'minC_rankdef': 1e9, 'ex_rankdef': None,
                                'minC_omdef': 1e9, 'ex_omdef': None,
                                'neg_dep': 0, 'minC_negdep': 1e9, 'ex_negdep': None})


def scan_N(fac, CMAX, stats, record=None, Nval=None, minM=2, maxwin=None):
    ps, es, pts = points_of(fac)
    N = 1
    for p in fac:
        N *= p ** fac[p]
    q = N ** 0.25
    ang = [(atan2(z[1], z[0]) % (2 * pi), z, a) for z, a in pts]
    ang.sort()
    n = len(ang)
    R = sqrt(N)
    for i in range(n):
        if ang[i][0] >= pi / 2:
            break
        j = i + 1
        while j < i + n:
            t = ang[j % n][0] + (2 * pi if j >= n else 0)
            span = t - ang[i][0]
            C = R * span / q
            if C > CMAX:
                break
            M = j - i + 1
            if maxwin and M > maxwin:
                break
            if M >= minM:
                win = [(ang[t % n][1], ang[t % n][2]) for t in range(i, j + 1)]
                analyse_window(ps, es, win, stats, record=record, Nval=N)
            j += 1


def gen_split_N(X):
    P1 = [p for p in primerange(5, X + 1) if p % 4 == 1]
    out = []

    def dfs(idx, cur, fac):
        if cur > 1:
            out.append((cur, dict(fac)))
        for t in range(idx, len(P1)):
            p = P1[t]
            if cur * p > X:
                break
            c = cur
            e = 0
            while c * p <= X:
                c *= p
                e += 1
                fac[p] = e
                dfs(t + 1, c, fac)
            del fac[p]
    dfs(0, 1, {})
    return out


def report(stats, title):
    print('==', title)
    print(' M  windows  max(M-1-omega)  max(M-1-rank)  minC(rank-def) ex  minC(omega-def) ex  negctrl-dep  minC ex')
    for M in sorted(stats):
        s = stats[M]
        print(M, s['n'], s['max_def_omega'], s['max_def_rank'],
              round(s['minC_rankdef'], 5) if s['ex_rankdef'] else '-', s['ex_rankdef'],
              round(s['minC_omdef'], 5) if s['ex_omdef'] else '-', s['ex_omdef'],
              s['neg_dep'], s['ex_negdep'])
    sys.stdout.flush()


if __name__ == '__main__':
    mode = sys.argv[1]
    stats = newstats()
    if mode == 'all':
        X = int(float(sys.argv[2]))
        CMAX = float(sys.argv[3])
        Ns = gen_split_N(X)
        print('number of N:', len(Ns))
        for N, fac in Ns:
            scan_N(fac, CMAX, stats, Nval=N)
        report(stats, 'all split N <= %g, C <= %g' % (X, CMAX))
    elif mode == 'list':
        CMAX = float(sys.argv[2])
        for s in sys.argv[3:]:
            N = int(s)
            fac = factorint(N)
            assert all(p % 4 == 1 for p in fac), N
            scan_N(fac, CMAX, stats, Nval=N)
        report(stats, 'list, C <= %g' % CMAX)
