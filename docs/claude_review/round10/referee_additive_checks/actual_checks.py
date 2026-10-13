"""Referee checks of additive.md on actual clusters (independent code).

Checks, on normalised clusters (maximal windows of angular width <= C N^(-1/4)):
 P31   Prop 3.1: pairs with mu > W/4 + log C' + log(2)/2 share one b = a_x + a_y - e, and {b = b*} is a matching.
 T32a  M <= C'^2 h* + 2.
 T32b  non-exceptional pairs: |u_xy| >= N^(1/4)/(C' h*), mu <= 2 log(sqrt2 C' h*), |B|/t <= C' h*;
       exceptional pairs (chord parallel to i v*) form a matching, empty if eta(v*) > 0.
 T32c  h* <= (2 pi/C')^(1/2) N^(1/8) + 1.
 TH11  |z_x - z_y|^2/Norm g = 2 B^2 or 4 B^2 (exactly one), B integer  (round9 L1L2 Thm 1.1).
 P27   I = [Z[i]:Lambda] (HNF) = gcd of all 2 Area; I even; odd q not | N, q < 200: q | I iff <= 2 residues;
       projection index of K on coordinates 4..M equals A_123/I.
 COV   covol(K)^2 * I^2 = sum over triples of (2 Area)^2   (exact; corrects the curvature loss of the
       round-9 Prop 5.2 bound M^(3/2) C^2 N^(1/2)).
 R44   M >= 5: every p_j <= max |2 Area| <= C'^3 N^(1/4)/2.
 RW    every LLL-reduced relation c in K with v(c) != 0 has
       w(v(c)) >= 3W/2 - 6 log C' - 2 log ||c||_1 + 2 log(96/sqrt2)   (referee Prop R1(ii)).
Exact integer arithmetic for all divisibility/rank/area/parallelism statements; mpmath (60 digits) for
angles and logs; numpy floats only to screen candidate directions v for h*, then recomputed in mpmath.
"""
import sys, math, random, itertools
import numpy as np
import mpmath as mp
from flint import fmpz_mat
from glib_ref import *

mp.mp.dps = 60
LOG = mp.log


def cross(u, v):
    return u[0] * v[1] - u[1] * v[0]


def int_kernel(rows_cols):
    """Saturated integer kernel of the 3 x M integer matrix B (given as list of M columns (1, re, im))."""
    M = len(rows_cols)
    lam = 10 ** 40
    A = []
    for x in range(M):
        row = [0] * M
        row[x] = 1
        row += [lam * t for t in rows_cols[x]]
        A.append(row)
    L = fmpz_mat(A).lll()
    ker = []
    for i in range(L.nrows()):
        r = [int(L[i, k]) for k in range(M + 3)]
        if r[M:] == [0, 0, 0] and any(r[:M]):
            ker.append(r[:M])
    return ker


def hstar(theta_lo, D, Np, Cp):
    """Return (h*, v*, eta*) by exhaustive search over primitive v, |v| <= bound + 3."""
    bound = float(mp.sqrt(2 * mp.pi / Cp) * mp.power(Np, mp.mpf(1) / 8) + 1) + 3
    V = int(math.ceil(bound))
    xs, ys = np.meshgrid(np.arange(0, V + 1), np.arange(-V, V + 1), indexing='ij')
    xs = xs.ravel(); ys = ys.ravel()
    keep = ((xs > 0) | ((xs == 0) & (ys > 0))) & (np.gcd(xs, ys) == 1) & (xs * xs + ys * ys <= bound * bound)
    xs = xs[keep]; ys = ys[keep]
    nrm = np.sqrt(xs * xs + ys * ys)
    al = np.arctan2(ys, xs) % np.pi
    tlo = float(theta_lo)
    Df = float(D)
    t = (al - tlo) % np.pi
    eta = np.where(t <= Df, 0.0, np.minimum(t - Df, np.pi - t))
    h = nrm * (1 + float(mp.power(Np, mp.mpf(1) / 4)) * eta / float(Cp))
    order = np.argsort(h)[:20]
    best = None
    for k in order:
        v = (int(xs[k]), int(ys[k]))
        a = mp.atan2(v[1], v[0]) % mp.pi
        tt = (a - theta_lo) % mp.pi
        et = mp.mpf(0) if tt <= D else min(tt - D, mp.pi - tt)
        hv = mp.sqrt(v[0] ** 2 + v[1] ** 2) * (1 + mp.power(Np, mp.mpf(1) / 4) * et / Cp)
        if best is None or hv < best[0]:
            best = (hv, v, et)
    return best, bound


def is_square(n):
    if n < 0:
        return False
    r = math.isqrt(n)
    return r * r == n


def check_cluster(zs0, stats, tag):
    zs, g = normalise(zs0)
    M = len(zs)
    Np = gnorm(zs[0])
    assert all(gnorm(z) == Np for z in zs)
    W = LOG(Np)
    th = [angle(z) for z in zs]
    # order by angle within the window (window width < pi/2, handle wrap)
    base = th[0]
    rel = [((t - base + mp.pi) % (2 * mp.pi)) - mp.pi for t in th]
    order = sorted(range(M), key=lambda i: rel[i])
    zs = [zs[i] for i in order]
    rel = [rel[i] for i in order]
    theta_lo = (base + rel[0]) % (2 * mp.pi)
    D = rel[-1] - rel[0]
    Cp = D * mp.power(Np, mp.mpf(1) / 4)
    N4 = mp.power(Np, mp.mpf(1) / 4)
    primes, e, pis, a = profile(zs, Np)
    r = len(primes)
    # pair data
    pairs = list(itertools.combinations(range(M), 2))
    mu = {}
    bvec = {}
    for (x, y) in pairs:
        gxy = ggcd(zs[x], zs[y])
        ng = gnorm(gxy)
        mu[(x, y)] = W / 2 - LOG(ng)
        d = sum(abs(a[x][j] - a[y][j]) * mp.log(primes[j]) for j in range(r))
        if abs((W / 2 - LOG(ng)) - (d - W / 2)) > mp.mpf(10) ** -40:
            stats['fail'].append((tag, 'gcd-profile', x, y))
        bvec[(x, y)] = tuple(a[x][j] + a[y][j] - e[j] for j in range(r))
        u = gsub(zs[x], zs[y])
        q, rem = divmod(gnorm(u), ng)
        if rem:
            stats['fail'].append((tag, 'TH11-div', x, y))
        s2 = (q % 2 == 0 and is_square(q // 2))
        s4 = (q % 4 == 0 and is_square(q // 4))
        if s2 == s4:
            stats['fail'].append((tag, 'TH11', x, y, q))
    # Prop 3.1
    thr = W / 4 + LOG(Cp) + LOG(2) / 2
    big = [p for p in pairs if mu[p] > thr]
    stats['P31_runs'] += 1
    if big:
        bs = set(bvec[p] for p in big)
        if len(bs) > 1:
            stats['fail'].append((tag, 'P31-two-b', len(bs)))
        bstar = bvec[big[0]]
        cls = [p for p in pairs if bvec[p] == bstar]
        verts = [v for p in cls for v in p]
        if len(verts) != len(set(verts)):
            stats['fail'].append((tag, 'P31-not-matching'))
        stats['P31_big'] += len(big)
    stats['P31_maxslack'] = max(stats['P31_maxslack'], max(
        [float(mu[p] - thr) for p in pairs if not big or bvec[p] != bvec[big[0]]] or [-1e9]))
    # Theorem 3.2
    (hs, vs, es), bound = hstar(theta_lo, D, Np, Cp)
    if M > Cp ** 2 * hs + 2:
        stats['fail'].append((tag, 'T32a', M, float(Cp ** 2 * hs + 2)))
    stats['T32a_maxratio'] = max(stats['T32a_maxratio'], float((M - 2) / (Cp ** 2 * hs)))
    cbound = mp.sqrt(2 * mp.pi / Cp) * mp.power(Np, mp.mpf(1) / 8) + 1
    if hs > cbound:
        stats['fail'].append((tag, 'T32c'))
    stats['T32c_maxratio'] = max(stats['T32c_maxratio'], float(hs / cbound))
    iv = (-vs[1], vs[0])
    exc = []
    for (x, y) in pairs:
        u = gsub(zs[x], zs[y])
        if cross(u, iv) == 0:
            exc.append((x, y))
            continue
        cu = math.gcd(u[0], u[1])
        uu = (u[0] // cu, u[1] // cu)
        if mp.sqrt(gnorm(uu)) < N4 / (Cp * hs) * (1 - mp.mpf(10) ** -30):
            stats['fail'].append((tag, 'T32b-u', x, y))
        lhs = mu[(x, y)] - 2 * LOG(mp.sqrt(2) * Cp * hs)
        stats['T32b_maxslack'] = max(stats['T32b_maxslack'], float(lhs))
        if lhs > mp.mpf(10) ** -30:
            stats['fail'].append((tag, 'T32b-mu', x, y))
        # |B|/t = N^(1/4) / (|g| |gamma|) with |gamma|^2 = q/B^2
        gxy = ggcd(zs[x], zs[y])
        tval = mp.sqrt(gnorm(u)) / N4
        q = gnorm(u) // gnorm(gxy)
        B2 = q // 2 if (q % 2 == 0 and is_square(q // 2)) else q // 4
        if mp.sqrt(B2) / tval > Cp * hs * (1 + mp.mpf(10) ** -30):
            stats['fail'].append((tag, 'T32b-B', x, y))
    verts = [v for p in exc for v in p]
    if len(verts) != len(set(verts)):
        stats['fail'].append((tag, 'T32b-exc-not-matching'))
    if es > 0 and exc:
        stats['fail'].append((tag, 'T32b-exc-but-eta>0'))
    stats['T32_runs'] += 1
    if M < 4:
        return
    # Prop 2.7, covolume, R44, relation weights
    stats['K_runs'] += 1
    A = {}
    for (x, y, w) in itertools.combinations(range(M), 3):
        A[(x, y, w)] = cross(gsub(zs[y], zs[x]), gsub(zs[w], zs[x]))
    I_all = 0
    for t in A.values():
        I_all = math.gcd(I_all, abs(t))
    # HNF index of the difference lattice
    diffs = [gsub(zs[x], zs[0]) for x in range(1, M)]
    H = fmpz_mat([[d[0], d[1]] for d in diffs]).hnf()
    I_hnf = abs(int(H[0, 0]) * int(H[1, 1]))
    if I_hnf != I_all or I_all % 2:
        stats['fail'].append((tag, 'P27-1', I_hnf, I_all))
    for q in [p for p in range(3, 200) if isprime(p)]:
        if Np % q == 0:
            continue
        nres = len(set(((z[0] % q), (z[1] % q)) for z in zs))
        if (I_all % q == 0) != (nres <= 2):
            stats['fail'].append((tag, 'P27-2', q))
        if I_all % q == 0:
            stats['P27_oddq'] += 1
    ker = int_kernel([(1, z[0], z[1]) for z in zs])
    if len(ker) != M - 3:
        stats['fail'].append((tag, 'rankK', len(ker)))
        return
    G = fmpz_mat(ker) * fmpz_mat(ker).transpose()
    cov2 = int(G.det())
    sA2 = sum(t * t for t in A.values())
    if cov2 * I_all * I_all != sA2:
        stats['fail'].append((tag, 'COV', cov2, sA2, I_all))
    stats['COV_ok'] += 1
    # projection index on coordinates 3..M-1 (0-indexed) vs A_012/I
    P = fmpz_mat([row[3:] for row in ker])
    idx = abs(int(P.det()))
    if idx * I_all != abs(A[(0, 1, 2)]):
        stats['fail'].append((tag, 'P27-3', idx, A[(0, 1, 2)], I_all))
    if M >= 5:
        mA = max(abs(t) for t in A.values())
        for p in primes:
            if p > mA:
                stats['fail'].append((tag, 'R44-p', p, mA))
        if mA > Cp ** 3 * N4 / 2 * (1 + mp.mpf(10) ** -30):
            stats['fail'].append((tag, 'R44-area'))
        stats['R44_runs'] += 1
    # relation weights (referee Prop R1(ii))
    for c in ker:
        n = sum(abs(t) for t in c)
        v = [sum(c[x] * a[x][j] for x in range(M)) for j in range(r)]
        if not any(v):
            stats['RW_v0'] += 1
            continue
        w = sum(abs(v[j]) * mp.log(primes[j]) for j in range(r))
        lb = 3 * W / 2 - 6 * LOG(Cp) - 2 * LOG(n) + 2 * LOG(96 / mp.sqrt(2))
        stats['RW_n'] += 1
        stats['RW_minslack'] = min(stats['RW_minslack'], float(w - lb))
        stats['RW_min_w_over_W'] = min(stats['RW_min_w_over_W'], float(w / W))
        if w < lb - mp.mpf(10) ** -30:
            stats['fail'].append((tag, 'RW', float(w), float(lb)))


def run(Ns, Cs, stats, label):
    for (N, fac) in Ns:
        pts = points_on_circle(fac)
        P = sorted(((z, angle(z)) for z in pts), key=lambda t: t[1])
        for C in Cs:
            D = mp.mpf(C) * mp.power(N, mp.mpf(-1) / 4)
            for win in windows(P, D):
                if len(win) < 3:
                    continue
                zs = [P[k][0] for k in win]
                stats['M_hist'][len(zs)] = stats['M_hist'].get(len(zs), 0) + 1
                check_cluster(zs, stats, (label, N, C, win[0]))


if __name__ == '__main__':
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 20261013)
    nrand = int(sys.argv[2]) if len(sys.argv) > 2 else 150
    stats = dict(fail=[], M_hist={}, P31_runs=0, P31_big=0, P31_maxslack=-1e9, T32_runs=0,
                 T32a_maxratio=0.0, T32b_maxslack=-1e9, T32c_maxratio=0.0, K_runs=0, P27_oddq=0,
                 COV_ok=0, R44_runs=0, RW_n=0, RW_v0=0, RW_minslack=1e9, RW_min_w_over_W=1e9)
    recs = [(1176852625, factorint(1176852625)), (537606725, factorint(537606725)),
            (235370525, factorint(235370525))]
    run(recs, [3, 4, 5.3, 5.9, 6, 7, 8], stats, 'rec')
    maxN = int(float(sys.argv[3])) if len(sys.argv) > 3 else 10 ** 13
    Ns = [random_split_N(rng, maxN) for _ in range(nrand)]
    run(Ns, [1, 2, 3, 4, 6] if len(sys.argv) > 4 else [1, 2, 3, 4], stats, "rnd")
    for k in sorted(stats):
        if k != 'fail':
            print(k, stats[k])
    print('failures', len(stats['fail']))
    for f in stats['fail'][:30]:
        print('  FAIL', f)
