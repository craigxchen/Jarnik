"""Margin / reduced-chord statistics of the clusters produced by scan.py.

Usage: python3 stats.py in.pkl [Mmin]
For each normalised cluster: exact pair gcds g_xy, reduced chords beta_xy = i^{k} gamma_s B_xy
(B_xy a nonzero rational integer, checked exactly), margins mu_xy = W/2 - log Norm g_xy,
identity mu_xy = log(|beta_xy|^2 / t_xy^2) with t_xy = |z_x - z_y| N^{-1/4} (checked in floating
point; it is an exact algebraic identity since Norm g_xy = N e^{-d_xy}), uniformity ratio
K = max mu / mean mu, pair demands at odd q <= X = 3M (q not dividing N) and point totals.
"""
import sys, pickle
from math import log, sqrt, atan2, pi
import numpy as np
from glib import reduced_chord_data, norm, sub, factor, primes_upto

data = pickle.load(open(sys.argv[1], "rb"))
Mmin = int(sys.argv[2]) if len(sys.argv) > 2 else 4


def analyse(N, fac, cl, Cc):
    M = len(cl)
    W = log(N)
    X = 3 * M
    odd_q = [q for q in primes_upto(X) if q > 2 and N % q != 0]
    mu = np.zeros((M, M)); D = np.zeros((M, M)); logB = np.zeros((M, M)); t = np.zeros((M, M))
    Bs = {}
    maxerr = 0.0
    for x in range(M):
        for y in range(x + 1, M):
            g, beta, s, B = reduced_chord_data(cl[x], cl[y])
            assert B != 0
            m = W / 2 - log(norm(g))
            tt = sqrt(norm(sub(cl[x], cl[y]))) / N ** 0.25
            ident = log(norm(beta)) - 2 * log(tt)
            maxerr = max(maxerr, abs(ident - m))
            mm = 1
            for q in odd_q:
                while B % (mm * q) == 0 and (B // mm) % q == 0:
                    mm *= q
            # recompute m_xy properly: product of q^{v_q(B)}
            mxy = 1
            for q in odd_q:
                b = abs(B)
                while b % q == 0:
                    b //= q; mxy *= q
            # Lemma 1.3 check (exact integer form): (M-2)^2 <= 2 C^4 Norm(g) i.e. M <= sqrt2 C^2 |g| + 2
            if not ((M - 2) ** 2 <= 2 * (Cc ** 4) * norm(g) * (1 + 1e-12)):
                LEMMA_FAIL.append((N, x, y))
            mu[x, y] = mu[y, x] = m
            D[x, y] = D[y, x] = 2 * log(mxy)
            logB[x, y] = logB[y, x] = log(abs(B))
            t[x, y] = t[y, x] = tt
            Bs[(x, y)] = (B, s)
    iu = np.triu_indices(M, 1)
    mus = mu[iu]
    res = dict(N=N, M=M, C=Cc, W=W, A=W / (M * log(M)), mean_mu=mus.mean(), max_mu=mus.max(),
               min_mu=mus.min(), K=mus.max() / mus.mean() if mus.mean() > 0 else float('inf'),
               lam_D=np.linalg.eigvalsh(D)[-1], Dmax=D[iu].max(), Dx_max=D.sum(1).max(),
               mean_logB=logB[iu].mean(), max_logB=logB[iu].max(), err=maxerr,
               argmax_t=t[iu][np.argmax(mus)], min_t=t[iu].min(),
               argmax_logB=logB[iu][np.argmax(mus)], coprime=abs(mus.max() - W / 2) < 1e-9)
    return res


rows = []
LEMMA_FAIL = []
for (N, fac, cl, Cc) in data:
    if len(cl) < Mmin:
        continue
    rows.append(analyse(N, fac, cl, Cc))
print("clusters analysed:", len(rows), " max identity error:", max(r['err'] for r in rows), " lemma1.3 failures:", len(LEMMA_FAIL))
pickle.dump(rows, open(sys.argv[1].replace('.pkl', '_stats.pkl'), 'wb'))
byM = {}
for r in rows:
    byM.setdefault(r['M'], []).append(r)
print("M  count  minC   minA   K(med) K(min) K(max)  mean_mu(med) max_mu(med) frac(argmax=closest)")
for M in sorted(byM):
    L = byM[M]
    Ks = np.array([r['K'] for r in L])
    print(M, len(L), "%.3f" % min(r['C'] for r in L), "%.2f" % min(r['A'] for r in L),
          "%.2f %.2f %.2f" % (np.median(Ks), Ks.min(), Ks.max()),
          "%.2f %.2f" % (np.median([r['mean_mu'] for r in L]), np.median([r['max_mu'] for r in L])),
          "%.2f" % np.mean([abs(r['argmax_t'] - r['min_t']) < 1e-12 for r in L]),
          "coprime-frac %.3f" % np.mean([r['coprime'] for r in L]),
          "K quartiles", np.round(np.percentile(Ks[np.isfinite(Ks)], [25, 75, 95]), 2))
