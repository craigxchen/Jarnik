#!/usr/bin/env python3
"""Judge's checks: (a) X.4 on genuine witnesses (tuples at inert primes, D(x), F2-rank, divisibility);
(b) a bounded-rank genuine parent (primes a^2 + (Q b)^2); (c) spot check of M_max(q_P, C) on q_P for P = 100, 200."""
import math, sys, itertools
import mpmath as mp

def isprime(n):
    if n < 2: return False
    i = 2
    while i * i <= n:
        if n % i == 0: return False
        i += 1
    return True

def gfac(p):
    a = 1
    while True:
        b = math.isqrt(p - a * a)
        if a * a + b * b == p: return (a, b)
        a += 1

def omega_mod(pi, p, l):
    # pi/conj(pi) = pi^2 / p  mod l, as pair in F_l[i]
    a, b = pi
    re, im = (a * a - b * b) % l, (2 * a * b) % l
    inv = pow(p, -1, l)
    return (re * inv % l, im * inv % l)

def legendre(a, l):
    return 1 if pow(a % l, (l - 1) // 2, l) == 1 else -1

def f2rank(vecs):
    rows = [int("".join("1" if v == -1 else "0" for v in vec), 2) for vec in vecs]
    rank = 0; basis = []
    for r in rows:
        for bvec in basis:
            r = min(r, r ^ bvec)
        if r:
            basis.append(r); rank += 1
    return rank

def part_a(P):
    S = [p for p in range(P, 2 * P + 1) if isprime(p) and p % 4 == 1]
    pis = {p: gfac(p) for p in S}
    inert = [l for l in range(3, 200) if isprime(l) and l % 4 == 3]
    first_full = None; viol = 0; pairs_checked = 0
    out = []
    for t in range(1, len(inert) + 1):
        ls = inert[:t]
        tup = {p: tuple(omega_mod(pis[p], p, l) for l in ls) for p in S}
        D = len(set(tup.values()))
        prodl = math.prod(ls)
        for j, k in itertools.combinations(S, 2):
            if tup[j] == tup[k]:
                pairs_checked += 1
                a, b = pis[j]; c, d = pis[k]
                im = b * c - a * d          # Im(pi_j conj pi_k)
                if im == 0 or im % prodl != 0 or prodl * prodl > j * k: viol += 1
        rk2 = f2rank([[legendre(p, l) for l in ls] for p in S])
        if t <= 8: out.append((ls[-1], D, rk2))
        if D == len(S) and first_full is None: first_full = ls[-1]
    print(f"P={P} s={len(S)}: (x, D(x), F2-rank of Legendre vectors) for first inert x: {out}")
    print(f"   first x with D(x)=s: {first_full}; equal-tuple pairs checked {pairs_checked}, "
          f"divisibility violations {viol}")

def part_b():
    Q = 3 * 7 * 11
    ps = []
    for bq in range(1, 40):
        for a in range(1, 4000, 2):
            n = a * a + (Q * bq) ** 2
            if isprime(n): ps.append((n, (a, Q * bq)))
            if len(ps) >= 12: break
        if len(ps) >= 12: break
    tuples = set(tuple(omega_mod(pi, p, l) for l in (3, 7, 11)) for p, pi in ps)
    print(f"bounded-rank genuine parent: {len(ps)} primes a^2+(231 b)^2, smallest {min(p for p,_ in ps)}, "
          f"distinct tuples at l=3,7,11: {len(tuples)} -> {tuples}")

def part_c(P, Cs=(1, 4, 8)):
    mp.mp.dps = 60
    S = [p for p in range(P, 2 * P + 1) if isprime(p) and p % 4 == 1]
    q = math.prod(S)
    th = [mp.atan2(gfac(p)[1], gfac(p)[0]) for p in S]
    s = len(S)
    # all angles: sum eps_p th_p + k pi/2  (eps = +-1), reduced mod 2pi; use floats for sort, mp for check
    import numpy as np
    thf = np.array([float(t) for t in th])
    angs = np.zeros(1)
    for t in thf:
        angs = np.concatenate([angs + t, angs - t])
    allang = np.concatenate([angs + k * math.pi / 2 for k in range(4)]) % (2 * math.pi)
    allang.sort()
    gaps = np.diff(np.concatenate([allang, [allang[0] + 2 * math.pi]]))
    q14 = mp.mpf(q) ** mp.mpf(0.25)
    g = float(gaps.min())
    print(f"P={P} s={s} points {len(allang)} log10 q={float(mp.log10(q)):.1f}: min gap * q^(1/4) = {g*float(q14):.4g}")
    for C in Cs:
        w = C / float(q14)
        j = 0; best = 0
        ext = np.concatenate([allang, allang + 2 * math.pi])
        n = len(allang)
        # two-pointer
        j = 0
        for i in range(n):
            while j < i + n and ext[j] - ext[i] <= w: j += 1
            best = max(best, j - i)
        print(f"   C={C}: M_max = {best}")

if __name__ == "__main__":
    part_a(100); part_a(1000)
    part_b()
    part_c(100); part_c(200)
