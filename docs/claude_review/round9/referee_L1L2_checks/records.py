"""Referee check of Thm 1.1, (1.1), Lemma 1.3 (bisector), Prop 3.1 and the Section-4 table on the record
clusters, with the referee's own code (refgauss.py).  All Gaussian arithmetic exact.
Usage: python3 records.py
"""
import sys, itertools
from math import log, sqrt, atan2, pi
from refgauss import *


def primes_upto(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


def profile(u):
    Np = gnorm(u[0])
    fac = factor(Np)
    assert all(p % 4 == 1 for p in fac), fac
    pis = {p: gauss_prime(p) for p in fac}
    rows, ks = [], []
    for z in u:
        a = {}
        rest = z
        for p, e in fac.items():
            pp = pis[p]
            ap = vpi(z, pp)
            ac = vpi(z, gconj(pp))
            assert ap + ac == e, (z, p, ap, ac, e)
            a[p] = ap
            rest = gexactdiv(rest, gmul(gpow(pp, ap), gpow(gconj(pp), e - ap)))
        assert rest in UNITS
        rows.append(a); ks.append(UNITS.index(rest))
    for p, e in fac.items():
        col = [r[p] for r in rows]
        assert min(col) == 0 and max(col) == e, "not normalised"
    return Np, fac, pis, rows, ks


GAMMA = {0: (0, 2), 1: (1, -1), 2: (2, 0), 3: (1, 1)}


def Dsel(A, s):
    X, Y = A
    return [Y, X - Y, X, X + Y][s]


def dlog_table(q):
    """norm-one torus of (Z[i]/q)^*: elements (a,b) mod q with a^2+b^2 = 1 mod q. Return generator data."""
    T = [(a, b) for a in range(q) for b in range(q) if (a * a + b * b) % q == 1]
    n = len(T)
    nq = q - (1 if q % 4 == 1 else -1)
    assert n == nq, (q, n, nq)
    def mulq(z, w):
        return ((z[0] * w[0] - z[1] * w[1]) % q, (z[0] * w[1] + z[1] * w[0]) % q)
    for gen in T:
        tab = {}
        x = (1, 0)
        for k in range(n):
            if x in tab:
                break
            tab[x] = k
            x = mulq(x, gen)
        if len(tab) == n:
            return tab, n, mulq
    raise Exception("not cyclic")


def invq(z, q):
    n = gnorm(z) % q
    ninv = pow(n, -1, q)
    return ((z[0] * ninv) % q, (-z[1] * ninv) % q)


def analyse(name, N, cl, C):
    M = len(cl)
    u, G = normalise(cl)
    Np, fac, pis, rows, ks = profile(u)
    W = log(Np)
    Cn = C / sqrt(sqrt(gnorm(G)))  # C' = C |G|^(-1/2)
    mus, comp, Bs, fails = {}, {}, {}, []
    bis_fail = 0
    omega_small = 0
    for x, y in itertools.combinations(range(M), 2):
        g = ggcd(u[x], u[y])
        # gcd formula: Norm g = prod p^(e - |a_x - a_y|)
        ng = 1
        for p, e in fac.items():
            ng *= p ** (e - abs(rows[x][p] - rows[y][p]))
        assert gnorm(g) == ng
        d = log(Np) - log(ng)
        mu = d - W / 2
        beta = gexactdiv(gsub(u[x], u[y]), g)
        # Thm 1.1 with the formula: A_v
        A = (1, 0)
        for p in fac:
            v = rows[x][p] - rows[y][p]
            if v > 0:
                A = gmul(A, gpow(pis[p], v))
            elif v < 0:
                A = gmul(A, gpow(gconj(pis[p]), -v))
        s = (ks[y] - ks[x]) % 4
        B = Dsel(A, s)
        assert B != 0
        gam = GAMMA[s]
        # beta = eps * gamma_s * B for some unit eps (g representative may differ by a unit)
        cand = gmul(gam, (B, 0))
        ok = any(gmul(e_, cand) == beta for e_ in UNITS)
        if not ok:
            fails.append((x, y))
        # exact (1.1) in integers: |u_x - u_y|^2 = Norm g * nu_s * B^2
        assert gnorm(gsub(u[x], u[y])) == ng * gnorm(gam) * B * B
        # bisector: u_x + u_y parallel to g (s even) or g(1+-i) (s odd)
        S = gadd(u[x], u[y])
        dirs = [g, gmul(g, (0, 1))] if s % 2 == 0 else [gmul(g, (1, 1)), gmul(g, (1, -1))]
        par = any((S[0] * dd[1] - S[1] * dd[0]) == 0 for dd in dirs)
        if not par:
            bis_fail += 1
        # Lemma 1.3: M <= sqrt2 C'^2 |g| + 2
        if not (M <= sqrt(2) * Cn ** 2 * sqrt(ng) + 2 + 1e-9):
            fails.append(("L1.3", x, y))
        # floating (1.1)
        t = sqrt(gnorm(gsub(u[x], u[y]))) / Np ** 0.25
        err = abs(mu - (log(gnorm(gam) * B * B) - 2 * log(t)))
        assert err < 1e-9, err
        # odd primes q <= M^2 not dividing Np dividing B
        om = sum(1 for q in primes_upto(M * M) if q > 2 and Np % q and B % q == 0)
        omega_small = max(omega_small, om)
        mus[(x, y)] = mu; comp[(x, y)] = (ng == 1); Bs[(x, y)] = (B, s, t)
    # Thm 1.1 congruence: for odd q not dividing N', a <= 3:
    cong_checks = 0
    for (x, y), (B, s, t) in Bs.items():
        for q in primes_upto(200):
            if q == 2 or Np % q == 0:
                continue
            for a in (1, 2, 3):
                Q = q ** a
                lhs = ((u[x][0] - u[y][0]) % Q == 0) and ((u[x][1] - u[y][1]) % Q == 0)
                lhs2 = ((cl[x][0] - cl[y][0]) % Q == 0) and ((cl[x][1] - cl[y][1]) % Q == 0)
                rhs = (B % Q == 0)
                assert lhs == rhs == lhs2, (x, y, q, a)
                cong_checks += 1
    # Prop 3.1 code check for q <= 60
    code_checks = 0
    for q in primes_upto(60):
        if q == 2 or Np % q == 0:
            continue
        tab, nq, mulq = dlog_table(q)
        ell = {}
        for p, pp in pis.items():
            rho = mulq((pp[0] % q, pp[1] % q), invq(gconj(pp), q))
            ell[p] = tab[rho]
        ui = tab[(0, 1)]
        w = [(sum(rows[x][p] * ell[p] for p in fac) + ks[x] * ui) % nq for x in range(M)]
        for (x, y), (B, s, t) in Bs.items():
            assert (w[x] == w[y]) == (B % q == 0) == (((u[x][0] - u[y][0]) % q == 0) and ((u[x][1] - u[y][1]) % q == 0))
            code_checks += 1
    vals = list(mus.values())
    mbar = sum(vals) / len(vals)
    mx = max(vals)
    ncomp = sum(comp.values())
    noncomp = [m for k_, m in mus.items() if not comp[k_]]
    print(f"{name}: N={N} M={M} C={C:.6f} |G|^2={gnorm(G)} N'={Np} W={W:.3f} A=W/(MlogM)={W/(M*log(M)):.3f}")
    print(f"   mubar={mbar:.3f} max mu={mx:.3f} (W/2={W/2:.3f}) min mu={min(vals):.3f} K={mx/mbar if mbar>0 else float('inf'):.3f} complementary={ncomp}")
    if noncomp:
        m2 = sum(noncomp) / len(noncomp)
        print(f"   K over non-complementary pairs only (max/mean both without them): {max(noncomp)/m2 if m2>0 else float('inf'):.3f};"
              f" author's K2 (max non-comp / mean of all) = {max(noncomp)/mbar:.3f}")
    print(f"   Thm1.1 formula failures={len(fails)} bisector failures={bis_fail} congruence checks={cong_checks} code checks={code_checks} max#odd q<=M^2 | B: {omega_small}")
    print(f"   sqrt2 C'^2+2 = {sqrt(2)*Cn**2+2:.2f}; min Norm g = {min(gnorm(ggcd(u[x],u[y])) for x,y in itertools.combinations(range(M),2))}, Lemma1.3 needs >= {(M-2)**2/(2*Cn**4):.4f}")
    return mus


if __name__ == "__main__":
    targets = [("9pt record", 1176852625, 9, 8.0), ("10pt record", 1176852625, 10, 8.3),
               ("7pt record", 235370525, 7, 5.3), ("8pt golden n=4", 537606725, 8, 5.9)]
    for name, N, M, Cm in targets:
        ws = windows(N, M, Cm)
        ws.sort(key=lambda t: t[1])
        seen = set()
        for cl, C in ws[:16]:
            key = tuple(sorted((abs(a), abs(b)) if abs(a) >= abs(b) else (abs(b), abs(a)) for a, b in cl))
            if key in seen:
                continue
            seen.add(key)
            analyse(name, N, cl, C)
