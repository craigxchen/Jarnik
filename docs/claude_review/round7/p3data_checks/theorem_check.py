"""Numerical check of Theorem 5.1 (p3data.md) and self-consistency of adapted profiles.

(1) Fbar_X(t) = P(D^ind > t) for the independent demand law, by DP on a log grid with step h:
    computed twice, with every term rounded down and rounded up, giving a rigorous bracket
    Fbar_lo(t) <= Fbar(t) <= Fbar_hi(t) (up to floating-point summation error ~1e-12).
(2) For fixed profiles (cut lists with real weights): Monte Carlo in the uniform residue model U_X
    (ell uniform mod m_q*, windings 0), estimate P(all pair inequalities), compare with the bound
    min_x 1/sum_y Fbar(mu_xy) (C = sqrt 2, c_C = 0), using Fbar_lo (so the bound is valid).
    Profiles: uniform (normalised Paley Hadamard, M-1 columns), round-5 Paley `construct` with equal
    weights, and the adapted semi-continuous profiles of q_run (fitted to one genuine D), each at
    several total weights W (weights scaled).
(3) Self-consistency with ACTUAL primes: for the realised prime profiles of q_run, the genuine
    demands of their own Gaussian primes (random conjugates, k = 0), exact via demands7; count
    violated pair inequalities; fraction of 400 conjugation choices with no violation.
Evidence only (Monte Carlo); the bracket in (1) is rigorous.
"""
import json, math, sys, random
import numpy as np
from gauss7 import odd_prime_powers_upto, mQ, primes_upto, profile_construct, paley_normalised, gaussian_primes
from demands7 import point_classes, collision_moduli, demand_matrix
from qcore import mask_to_s

H_STEP = 1e-3


def fbar_tables(X):
    byq = {}
    for (q, a, Q) in odd_prime_powers_upto(X):
        byq.setdefault(q, []).append(a)
    tot = sum(2 * math.log(q) * max(a) for q, a in byq.items())
    nb = int(tot / H_STEP) + 10
    out = []
    for rnd in (math.floor, math.ceil):
        dist = np.zeros(nb); dist[0] = 1.0
        for q, levels in byq.items():
            A = max(levels)
            # P(L_q = a) = 1/m_{q^a} - 1/m_{q^{a+1}}, L_q in 0..A
            pr = [1.0 - 1.0 / mQ(q, 1)] + [1.0 / mQ(q, a) - (1.0 / mQ(q, a + 1) if a < A else 0.0) for a in range(1, A + 1)]
            new = np.zeros(nb)
            for a, p in enumerate(pr):
                sh = int(rnd(2 * a * math.log(q) / H_STEP))
                new[sh:] += p * dist[:nb - sh]
            dist = new
        tail = dist[::-1].cumsum()[::-1]       # tail[i] = P(value >= i*h)
        out.append(tail)
    return out  # (floor-rounded, ceil-rounded)


def fbar_lo(tabs, t):
    """rigorous lower bound for P(D > t): use floor-rounded values (each <= true), P(floor > t)<=P(true>t)"""
    lo = tabs[0]
    i = int(math.floor(t / H_STEP)) + 1
    return float(lo[i]) if i < len(lo) else 0.0


def sample_D(Aprof, X, rng):
    """one draw of the unit-1 demand matrix in U_X for profile A (M x r 0/1)."""
    M, r = Aprof.shape
    D = np.zeros((M, M))
    byq = {}
    for (q, a, Q) in odd_prime_powers_upto(X):
        byq.setdefault(q, []).append(a)
    for q, levels in byq.items():
        A = max(levels); m = mQ(q, A)
        ell = rng.integers(0, m, size=r)
        kap = (Aprof @ ell) % m
        for a in levels:
            ma = mQ(q, a)
            D += 2 * math.log(q) * ((kap[:, None] - kap[None, :]) % ma == 0)
    np.fill_diagonal(D, 0)
    return D


def margins_of(S, w):
    sep = (S[:, None, :] != S[None, :, :])
    mu = (sep * w).sum(2) - w.sum() / 2
    np.fill_diagonal(mu, 0)
    return mu


def check_profile(name, S, w, X, tabs, nmc, rng, scales):
    M = S.shape[0]
    A = ((1 - S) // 2).astype(np.int64)
    res = []
    for sc in scales:
        ww = w * sc
        mu = margins_of(S, ww)
        W = ww.sum()
        b = min(1.0 / max(sum(fbar_lo(tabs, mu[x, y]) for y in range(M) if y != x), 1e-300) for x in range(M))
        b51 = 2.0 / ((M - 1) * max(fbar_lo(tabs, W / (M - 1)), 1e-300))
        ok = 0
        for _ in range(nmc):
            D = sample_D(A, X, rng)
            if (mu - D)[np.triu_indices(M, 1)].min() >= 0:
                ok += 1
        res.append(dict(name=name, W=W, phat=ok / nmc, bound=min(b, 1.0), bound51=min(b51, 1.0)))
        print('  %-10s W=%8.1f  P_hat=%.4f  bound(min_x)=%.4f  bound(5.1)=%.4f' % (name, W, ok / nmc, min(b, 1), min(b51, 1)))
        sys.stdout.flush()
    return res


def main():
    rng = np.random.default_rng(7)
    allres = []
    for M, Xf in ((12, 2), (24, 2)):
        X = Xf * M
        tabs = fbar_tables(X)
        print('M=%d X=%d: Fbar bracket check at t=2log X: [%.4g, %.4g]' % (M, X, fbar_lo(tabs, 2 * math.log(X)),
              float(tabs[1][int(math.ceil(2 * math.log(X) / H_STEP))])))
        # uniform: normalised Paley Hadamard columns 1..M-1
        Hm = paley_normalised(M - 1)
        S_u = Hm[:, 1:]
        w_u = np.ones(M - 1)
        # Paley construct b = 8, equal weights
        S_c = profile_construct(M - 1, 8)
        w_c = np.ones(S_c.shape[1])
        Wt = [float(v) for v in np.array([2, 4, 8, 16]) * M * math.log(M)]
        allres += check_profile('uniform', S_u, w_u, X, tabs, 2000, rng, [t / (M - 1) for t in Wt])
        allres += check_profile('paley-b8', S_c, w_c, X, tabs, 2000, rng, [t / S_c.shape[1] for t in Wt])
        # adapted profiles from q_run (no-sep prime profiles), if available
        try:
            recs = [json.loads(l) for l in open('q_results_M%d.jsonl' % M)]
        except FileNotFoundError:
            recs = []
        for rec in recs[:2]:
            if rec['X'] != X or 'prime_cols' not in rec:
                continue
            cols = rec['prime_cols']
            S_a = np.array([mask_to_s(m, M) for m, pl in cols for p in pl]).T
            w_a = np.array([math.log(p) for m, pl in cols for p in pl])
            W0 = w_a.sum()
            allres += check_profile('adapted-t%d' % rec['t'], S_a, w_a, X, tabs, 2000, rng, [t / W0 for t in Wt])
    json.dump(allres, open('theorem_check.json', 'w'))
    # (3) self-consistency with actual primes
    print('Self-consistency of realised adapted profiles with their own actual primes:')
    for M in (8, 12, 20, 24):
        try:
            recs = [json.loads(l) for l in open('q_results_M%d.jsonl' % M)]
        except FileNotFoundError:
            continue
        for rec in recs:
            for key in ('prime_cols', 'prime_sep_cols'):
                if key not in rec:
                    continue
                cols = rec[key]; X = rec['X']
                S = np.array([mask_to_s(m, M) for m, pl in cols for p in pl]).T
                ps = [p for m, pl in cols for p in pl]
                A = ((1 + S) // 2).astype(np.int64)
                w = np.log(np.array(ps, dtype=float))
                mu = margins_of(S, w)
                iu = np.triu_indices(M, 1)
                nok, nviol = 0, []
                prng = random.Random(1234 + M + rec['t'])
                for trial in range(400):
                    gA, gB = gaussian_primes(ps, prng)
                    cl = point_classes(A, gA, gB, ps, X)
                    m1, _ = collision_moduli(cl, X, M)
                    D = demand_matrix(m1)
                    v = int(((mu - D)[iu] < 0).sum())
                    nviol.append(v); nok += (v == 0)
                print('  M=%d t=%d X=%d %-14s r=%d W=%.1f: own genuine demands violate %.1f pairs on average (max D %.1f); feasible in %d/400 conjugation choices'
                      % (M, rec['t'], X, key, len(ps), w.sum(), np.mean(nviol), D.max(), nok))
                sys.stdout.flush()


if __name__ == '__main__':
    main()
