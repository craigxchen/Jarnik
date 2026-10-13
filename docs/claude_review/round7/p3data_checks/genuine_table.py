"""Aggregate genuine_trials.jsonl (part 1) and compare with the independent uniform residue model.

Independent model (null): for every odd prime power Q = q^a <= X, classes kappa_x uniform in Z/m_Q,
nested in a (the level-a class reduces to the level-(a-1) class), independent over q and x.
Unit-1 demand D_xy = 2 sum_{Q <= X : collision} log q.  Simulated with the same M, X.
Output: genuine_table_out.txt.  Evidence only (statistics of exact per-trial data).
"""
import json, math, sys
import numpy as np
from collections import defaultdict
from gauss7 import odd_prime_powers_upto, mQ
from demands7 import expected_mean_independent

rng = np.random.default_rng(2024)


def indep_D(M, X):
    D = np.zeros((M, M))
    byq = defaultdict(list)
    for (q, a, Q) in odd_prime_powers_upto(X):
        byq[q].append(a)
    for q, levels in byq.items():
        m1 = mQ(q, 1)
        cls = rng.integers(0, m1, size=M)
        mod = m1
        for a in sorted(levels):
            if a > 1:
                cls = cls + mod * rng.integers(0, q, size=M)   # lift: Z/m_{q^a} -> uniform over fibre
                mod = mQ(q, a)
            eq = cls[:, None] == cls[None, :]
            D += 2 * math.log(q) * eq
    np.fill_diagonal(D, 0)
    return D


def dstats(D):
    M = D.shape[0]
    iu = np.triu_indices(M, 1)
    v = D[iu]
    return dict(mean=v.mean(), max=v.max(), lam=np.linalg.eigvalsh(D)[-1], maxrow=D.sum(1).max(), vals=v)


def tail_from_hist(hist, thr):
    """number of pairs with D = 2 log m >= thr"""
    return sum(c for k, c in hist.items() if 2 * math.log(int(k)) >= thr - 1e-12)


def main():
    recs = [json.loads(l) for l in open('genuine_trials.jsonl')]
    out = []
    P = lambda s: (out.append(s), print(s))
    P("# Part 1 table: unit-1 genuine pair demands D = 2 log m_xy (exact per trial), averaged over trials")
    P("# L = log M, H = log^2 M / loglog M.  'unif/Perron' = (M-1) D_max / lambda_max(D) = (uniform-margin cost 2(M-1)D_max) / (Perron bound 2 lambda_max).")
    P("# cols: M | profile | mode | X | trials | mean D (exact indep. expectation) | D_max median [min,max] | D_max/L | D_max/H | lam/(M L) | maxrow/(M L) | unif/Perron median")
    groups = defaultdict(list)
    for r in recs:
        groups[(r['M'], r['variant'], r['mode'])].append(r)
    MS = sorted({r['M'] for r in recs})
    summary = {}
    for M in MS:
        L = math.log(M); H = L * L / math.log(L)
        for Xf in (2, 4):
            X = Xf * M
            Em = expected_mean_independent(M, X)
            for key in [(M, 'construct', 'window'), (M, 'lower', 'window'), (M, 'construct', 'iid')]:
                g = groups[key]
                st = [r['u1_X%d' % Xf] for r in g]
                means = np.array([s['mean'] for s in st]); mx = np.array([s['max'] for s in st])
                lam = np.array([s['lam'] for s in st]); mr = np.array([s['maxrow'] for s in st])
                ratio = (M - 1) * mx / lam
                P("%3d | %-9s | %-6s | %3d | %3d | %6.3f (%6.3f) | %6.2f [%6.2f,%6.2f] | %5.2f | %5.2f | %5.3f | %5.3f | %5.2f"
                  % (M, key[1], key[2], X, len(g), means.mean(), Em, np.median(mx), mx.min(), mx.max(),
                     np.median(mx) / L, np.median(mx) / H, lam.mean() / (M * L), mr.mean() / (M * L), np.median(ratio)))
                summary[(M, X, key[1], key[2])] = dict(mean=means.mean(), dmax=np.median(mx), lam=lam.mean(), ratio=np.median(ratio))
            # independent model
            ntr = 100
            sts = [dstats(indep_D(M, X)) for _ in range(ntr)]
            mx = np.array([s['max'] for s in sts]); lam = np.array([s['lam'] for s in sts])
            mr = np.array([s['maxrow'] for s in sts]); means = np.array([s['mean'] for s in sts])
            ratio = (M - 1) * mx / lam
            P("%3d | %-9s | %-6s | %3d | %3d | %6.3f (%6.3f) | %6.2f [%6.2f,%6.2f] | %5.2f | %5.2f | %5.3f | %5.3f | %5.2f"
              % (M, 'indep', 'model', X, ntr, means.mean(), Em, np.median(mx), mx.min(), mx.max(),
                 np.median(mx) / L, np.median(mx) / H, lam.mean() / (M * L), mr.mean() / (M * L), np.median(ratio)))
            summary[(M, X, 'indep', 'model')] = dict(mean=means.mean(), dmax=np.median(mx), lam=lam.mean(), ratio=np.median(ratio))
    # tails
    P("")
    P("# Tail: average number of pairs per trial with D >= t*L (unit 1), genuine (construct, window) vs independent model")
    P("# cols: M | X | t=2 | t=3 | t=4 | t=5 | t=6 | t=8   (genuine / indep)")
    for M in MS:
        L = math.log(M)
        for Xf in (2, 4):
            X = Xf * M
            g = groups[(M, 'construct', 'window')]
            sims = [indep_D(M, X) for _ in range(100)]
            iu = np.triu_indices(M, 1)
            row = []
            for t in (2, 3, 4, 5, 6, 8):
                gen = np.mean([tail_from_hist(r['u1_X%d' % Xf]['hist'], t * L) for r in g])
                ind = np.mean([(S[iu] >= t * L - 1e-12).sum() for S in sims])
                row.append("%8.2f/%-8.2f" % (gen, ind))
            P("%3d | %3d | " % (M, X) + " | ".join(row))
    # window collision counts vs Poisson
    P("")
    P("# Window (sqrt M, X] level-1 unit-1 collision counts per pair: fraction of pairs with c collisions, genuine (construct) vs exact Poisson-binomial law of independent residues")
    from gauss7 import primes_upto
    for M in MS:
        for Xf in (2, 4):
            X = Xf * M
            lamQ = sum(1.0 / mQ(q, 1) for q in primes_upto(X) if q > 2 and q > math.sqrt(M))
            g = groups[(M, 'construct', 'window')]
            tot = np.zeros(12)
            for r in g:
                c = np.array(r['cnt_X%d' % Xf], dtype=float)
                tot[:len(c)] += c[:12]
            frac = tot / tot.sum()
            # exact Poisson-binomial law of sum of independent Bernoulli(1/m_q), q in window
            pb = np.zeros(13); pb[0] = 1.0
            for q in primes_upto(X):
                if q > 2 and q > math.sqrt(M):
                    p = 1.0 / mQ(q, 1); pb[1:] = pb[1:] * (1 - p) + pb[:-1] * p; pb[0] *= (1 - p)
            pois = pb[:12]
            kmax = max(k for k in range(12) if tot[k] > 0)
            P("%3d X=%3d lambda_Q=%.3f | " % (M, X, lamQ) + " ".join("c=%d: %.2e/%.2e" % (k, frac[k], pois[k]) for k in range(kmax + 1)))
    open('genuine_table_out.txt', 'w').write("\n".join(out) + "\n")
    json.dump({str(k): v for k, v in summary.items()}, open('genuine_summary.json', 'w'))


if __name__ == '__main__':
    main()
