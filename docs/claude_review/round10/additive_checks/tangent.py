"""Checks for additive.md Section 3 (lattice geometry of the tangent range).

For actual clusters (maximal runs of lattice points of x^2+y^2=N within angular width C N^(-1/4)):
  Prop 3.1  every pair with mu > W/4 + log(sqrt2*C) lies in ONE class of equal b = a_x+a_y-e
            (tested in the equivalent exact form: equal primitive chord direction), and that class
            is a matching;
  Thm 3.2   with h* = min over primitive v of |v|(1 + N^(1/4) eta(v)/C):
            (a) M <= C^2 h* + 2,  (b) mu_xy <= 2 log(2 C h*) for chords not parallel to i v*.
            Also the Minkowski bound h* <= sqrt(2 pi/C) N^(1/8) + 1.
Exact integer arithmetic for gcds, chords and parallelism; mpmath (60 digits) for angles.
Usage: python3 tangent.py
"""
import itertools, sys
from math import gcd, log, sqrt, isqrt
import mpmath
from glib import (mul, conj, norm, sub, ggcd, exact_div, factor, points_of_norm, primes_upto)

mpmath.mp.dps = 60


def prim(z):
    c = gcd(abs(z[0]), abs(z[1]))
    return (z[0] // c, z[1] // c)


def normalise(cl):
    g = cl[0]
    for z in cl[1:]:
        g = ggcd(g, z)
    cl = [exact_div(z, g) for z in cl]
    return cl, norm(cl[0])


def angles(cl):
    a = [mpmath.atan2(z[1], z[0]) for z in cl]
    lo = min(a)
    if max(a) - lo > mpmath.pi:
        a = [x % (2 * mpmath.pi) for x in a]
    return a


def dist_mod_pi(t, arc_lo, arc_hi):
    """angular distance from the angle t + pi Z to the interval [arc_lo, arc_hi]."""
    best = None
    for k in range(-6, 7):
        s = t + k * mpmath.pi
        if arc_lo <= s <= arc_hi:
            return mpmath.mpf(0)
        d = min(abs(s - arc_lo), abs(s - arc_hi))
        best = d if best is None or d < best else best
    return best


def hstar(N, lo, hi, C, vmax):
    """min over primitive v with |v| <= vmax of h(v) = |v|(1 + N^(1/4) eta/C)."""
    N4 = mpmath.mpf(N) ** mpmath.mpf(0.25)
    best, bv = None, None
    for a in range(0, vmax + 1):
        for b in range(-vmax, vmax + 1):
            if a == 0 and b <= 0:
                continue
            if a * a + b * b > vmax * vmax or gcd(a, b) != 1:
                continue
            if best is not None and sqrt(a * a + b * b) >= best:
                continue
            eta = dist_mod_pi(mpmath.atan2(b, a), lo, hi)
            h = mpmath.sqrt(a * a + b * b) * (1 + N4 * eta / C)
            if best is None or h < best:
                best, bv = float(h), (a, b)
    return best, bv


def analyse(cl, label, out):
    cl, N = normalise(cl)
    M = len(cl)
    if M < 3:
        return None
    W = log(N)
    ang = angles(cl)
    lo, hi = min(ang), max(ang)
    D = hi - lo
    C = float(D * mpmath.mpf(N) ** mpmath.mpf(0.25))
    vmax = int(sqrt(2 * 3.1416 / C) * N ** 0.125) + 3
    h, v = hstar(N, lo, hi, C, vmax)
    iv = (-v[1], v[0])
    viol = []
    # (a)
    if M > C * C * h + 2 + 1e-9:
        viol.append("a")
    if h > sqrt(2 * 3.141592653589793 / C) * N ** 0.125 + 1 + 1e-9:
        viol.append("mink")
    big = []          # pairs with mu > W/4 + log(sqrt2 C)
    exc = []          # pairs with chord parallel to i v*
    worst_b = -1e9
    for x, y in itertools.combinations(range(M), 2):
        g = ggcd(cl[x], cl[y])
        mu = W / 2 - log(norm(g))
        u = prim(sub(cl[x], cl[y]))
        assert norm(u) <= 2 * norm(g)
        par = (u[0] * iv[1] - u[1] * iv[0]) == 0
        if par:
            exc.append((x, y))
        else:
            worst_b = max(worst_b, mu - 2 * log(sqrt(2) * C * h))
        if mu > W / 4 + log(sqrt(2) * C) + 1e-9:
            big.append((x, y, u))
    if worst_b > 1e-9:
        viol.append("b")
    # Prop 3.1: all big pairs share one primitive chord direction, and form a matching
    dirs = set()
    for x, y, u in big:
        uu = u if (u[0] > 0 or (u[0] == 0 and u[1] > 0)) else (-u[0], -u[1])
        dirs.add(uu)
    pts = [p for x, y, u in big for p in (x, y)]
    if len(dirs) > 1 or len(pts) != len(set(pts)):
        viol.append("P3.1")
    for P in (exc,):
        pp = [p for e in P for p in e]
        if len(pp) != len(set(pp)):
            viol.append("matching")
    out.append((label, N, M, C, W, h, v, len(exc), len(big), worst_b, viol))
    return viol


def all_runs(N, Cmax):
    pts = points_of_norm(factor(N))
    n = len(pts)
    a = [mpmath.atan2(z[1], z[0]) % (2 * mpmath.pi) for z in pts]
    N4 = mpmath.mpf(N) ** mpmath.mpf(0.25)
    runs = []
    for i in range(n):
        j = i
        while True:
            jn = (j + 1) % n
            if jn == i or ((a[jn] - a[i]) % (2 * mpmath.pi)) * N4 > Cmax:
                break
            j = jn
        L = (j - i) % n + 1
        if L >= 3:
            runs.append([pts[(i + t) % n] for t in range(L)])
    return runs


if __name__ == "__main__":
    out = []
    nviol = 0
    ncl = 0
    split = [p for p in primes_upto(400) if p % 4 == 1]
    Ns = set()
    import random
    random.seed(1)
    for _ in range(220):
        k = random.randint(3, 6)
        ps = random.sample(split[:30], k)
        N = 1
        for p in ps:
            N *= p ** random.randint(1, 2)
        if N < 10 ** 16:
            Ns.add(N)
    Ns |= {1176852625, 235370525, 537606725}
    for N in sorted(Ns):
        for cl in all_runs(N, 4.0):
            ncl += 1
            v = analyse(cl, "run", out)
            if v:
                nviol += 1
                print("VIOLATION", N, v)
    print("clusters", ncl, "violations", nviol)
    Ms = {}
    for r in out:
        Ms[r[2]] = Ms.get(r[2], 0) + 1
    print("M histogram", sorted(Ms.items()))
    # summary: distribution of h*, exceptional class sizes, and max ratio mu/(2 log 2 C h*)
    exc_nonempty = sum(1 for r in out if r[7] > 0)
    print("clusters with nonempty exceptional class:", exc_nonempty)
    print("max of mu - 2log(sqrt2 C h*) over non-exceptional pairs:", max(r[9] for r in out))
    big = sorted(out, key=lambda r: -r[2])[:8]
    for r in big:
        print(f"N={r[1]} M={r[2]} C={r[3]:.3f} W={r[4]:.1f} h*={r[5]:.2f} v*={r[6]} "
              f"C^2h*+2={r[3]**2*r[5]+2:.1f} #exc={r[7]} #big={r[8]} slack_b={r[9]:.2f}")
