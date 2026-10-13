"""Margin / reduced-chord data for the record clusters and the algebraic families
(golden 8-point family, Pell 6-point family) of growth/families.md.

All arithmetic exact (Python ints); logs printed in floating point.
For each cluster: M, C, W = log N', A = W/(M log M), mean/max/min margin, K = max/mean,
mean log|B|, max log|B|, which pair attains max mu (its gap t and log|B|), and the largest
number of distinct odd primes q <= M^2 (q not dividing N) dividing a single B_xy.
"""
from math import log, sqrt, atan2, pi
import itertools
from glib import (mul, conj, norm, sub, ggcd, exact_div, gpow, factor, points_of_norm,
                  reduced_chord_data, primes_upto)


def normalise(cl):
    g = cl[0]
    for z in cl[1:]:
        g = ggcd(g, z)
    cl = [exact_div(z, g) for z in cl]
    return cl, norm(cl[0])


def arc_const(cl, N):
    import mpmath
    mpmath.mp.dps = 120
    ang = sorted(mpmath.atan2(z[1], z[0]) for z in cl)
    if ang[-1] - ang[0] > mpmath.pi:
        ang = sorted((a % (2 * mpmath.pi)) for a in ang)
    return float((ang[-1] - ang[0]) * mpmath.mpf(N) ** mpmath.mpf(0.25))


def analyse(cl, label):
    cl, N = normalise(cl)
    M = len(cl)
    W = log(N)
    C = arc_const(cl, N)
    qs = [q for q in primes_upto(M * M) if q > 2 and N % q]
    mus, rows = [], []
    for x, y in itertools.combinations(range(M), 2):
        g, beta, s, B = reduced_chord_data(cl[x], cl[y])
        mu = W / 2 - log(norm(g))
        t = sqrt(norm(sub(cl[x], cl[y]))) / N ** 0.25
        om = sum(1 for q in qs if B % q == 0)
        mus.append(mu)
        rows.append((mu, t, log(abs(B)), om, B, s))
    mean = sum(mus) / len(mus)
    mx = max(rows)
    print(f"{label}: N'={N} M={M} C={C:.4f} W={W:.2f} A={W/(M*log(M)):.2f} "
          f"mean_mu={mean:.2f} max_mu={mx[0]:.2f} min_mu={min(mus):.2f} K={mx[0]/mean:.2f} "
          f"mean_logB={sum(r[2] for r in rows)/len(rows):.2f} max_logB={max(r[2] for r in rows):.2f} "
          f"argmax:(t={mx[1]:.3f},logB={mx[2]:.2f}) max_omega_small={max(r[3] for r in rows)}")
    return rows


def largest_clusters(N, Cmax):
    pts = points_of_norm(factor(N))
    n = len(pts)
    ang = [atan2(z[1], z[0]) % (2 * pi) for z in pts]
    best = []
    for i in range(n):
        j = i
        while True:
            jn = (j + 1) % n
            if jn == i or ((ang[jn] - ang[i]) % (2 * pi)) * N ** 0.25 > Cmax:
                break
            j = jn
        cl = [pts[(i + t) % n] for t in range((j - i) % n + 1)]
        best.append(cl)
    m = max(len(c) for c in best)
    return [c for c in best if len(c) == m][0]


if __name__ == "__main__":
    for N, Cm in [(1176852625, 8), (235370525, 5.3), (537606725, 5.9)]:
        analyse(largest_clusters(N, Cm), f"record N={N}")
    # golden family Z_w(n), n = 1 mod 3
    fib = [0, 1]
    for _ in range(200):
        fib.append(fib[-1] + fib[-2])
    F = (1, 2)
    for n in [4, 7, 10, 13, 22, 31, 40]:
        pts = []
        for w in itertools.product([0, 1], repeat=4):
            p = sum(w)
            if p % 2 == 0:
                continue
            z = (1, 0)
            if p == 1:
                z = mul(z, F)
            if p == 3:
                z = mul(z, conj(F))
            for s in range(4):
                J = (fib[n + s + 1], fib[n + s])
                z = mul(z, J if w[s] else conj(J))
            pts.append(z)
        analyse(pts, f"golden n={n}")
    # Pell family: a = (y1^2-143)/2, points (a,0)? use the six points (a-k, +-y_k), k=0,1,2 with y0=12
    y2, y1 = 1656, 1171
    for it in range(4):
        a = (y1 * y1 - 143) // 2
        N = a * a + 144
        pts = []
        for k, y in [(0, 12), (1, y1), (2, y2)]:
            if k == 0:
                pts += [(a, 12), (a, -12)]
            else:
                pts += [(a - k, y), (a - k, -y)]
        # rotate so that the cluster is around the x-axis: points are (x, y) with y small
        for z in pts:
            assert norm(z) == N, (z, N)
        analyse(pts, f"Pell a={a}")
        # next solution: multiply y2 + y1 sqrt2 by (3+2 sqrt2)^2 = 17 + 12 sqrt2
        y2, y1 = 17 * y2 + 24 * y1, 12 * y2 + 17 * y1
