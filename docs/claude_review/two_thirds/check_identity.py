#!/usr/bin/env python3
"""Exact checks of every algebraic/arithmetic step of proof.md on actual lattice-point
configurations (numerical evidence for the algebra; the proofs do not depend on it).

For several norms N0 (including ones with inert factors, powers of 2 and repeated split
primes, so that the gcd normalisation is exercised), enumerate all lattice points of norm
N0, sort by angle, and for every window of M consecutive points (an arc) and also for
random subsets, check exactly:

  (N)  after dividing by the Gaussian gcd g of the configuration: N = N0/|g|^2 is odd and
       has no prime factor = 3 mod 4; for p | N, the levels a_ip = v_pi(z_i) satisfy
       v_pi(z_i) + v_pibar(z_i) = e_p, min_i a_ip = 0, max_i a_ip = e_p;
  (G)  Norm(gcd(z_i,z_j)) = prod_p p^(e_p - |a_ip - a_jp|);
  (I)  the exact identity (2.3), in its integer form (2.4),
          prod_{i<j} Norm(z_i-z_j)^4
            = N^(M(M-2)) * prod_p p^(sum_tau (2 S_ptau - M)^2) * prod_{i<j} Norm(c_ij)^4 ;
  (R)  ramified prime: (1+i) | c_ij for every pair;
  (Q)  inert l: l^a | c_ij  <=>  z_i = z_j mod l^a  (a = 1, 2), and the number of residue
       classes mod l^a of Gaussian integers of norm = N mod l^a is (l+1) l^(a-1);
  (S1) split p not | N: p^a | c_ij <=> z_i = z_j mod p^a; pi | (z_i - z_j) <=> p | (z_i - z_j);
       class count (p-1) p^(a-1);
  (S2) split p | N: same level => (p | c_ij <=> r_i = r_j mod p), r_i = z_i/(pi^a pibar^(e-a));
       different levels => pi does not divide c_ij and pibar does not divide c_ij;
  (L)  the per-prime collision counts are >= the pigeonhole minima E(.,.), and
       sum_{i<j} log|c_ij| >= (sum of the per-prime lower bounds (3.6)),
       and the exact finite inequality (5.1) holds with the actual C = L/sqrt(R).
"""
import random
import sys
from math import atan2, comb, log, pi as PI, sqrt
from gauss import (sub, mul, conj, norm, divides, exact_div, gcd_g, factor_int,
                   split_prime_factor, valuation, gpow, lattice_points, is_prime_int)


def E(n, q):
    if n <= 0:
        return 0
    b, r = divmod(n, q)
    return q * comb(b, 2) + r * b


def residue(z, m):
    return (z[0] % m, z[1] % m)


def class_count(N, m):
    """number of residue classes z mod m with Norm(z) = N mod m."""
    return sum(1 for x in range(m) for y in range(m) if (x * x + y * y - N) % m == 0)


def check_config(pts, primes_small, stats):
    M = len(pts)
    # (N) normalisation
    g = pts[0]
    for z in pts[1:]:
        g = gcd_g(g, z)
    zs = [exact_div(z, g) for z in pts]
    N = norm(zs[0])
    assert all(norm(z) == N for z in zs)
    assert N % 2 == 1
    fac = factor_int(N)
    assert all(p % 4 == 1 for p in fac), fac
    pis = {p: split_prime_factor(p) for p in fac}
    lev = {}
    for p, e in fac.items():
        pi_ = pis[p]
        a = [valuation(pi_, z) for z in zs]
        b = [valuation(conj(pi_), z) for z in zs]
        assert all(a[i] + b[i] == e for i in range(M))
        assert min(a) == 0 and max(a) == e
        lev[p] = a
    # pair data
    c = {}
    for i in range(M):
        for j in range(i + 1, M):
            gij = gcd_g(zs[i], zs[j])
            ng = 1
            for p, e in fac.items():
                ng *= p ** (e - abs(lev[p][i] - lev[p][j]))
            assert norm(gij) == ng  # (G)
            cij = exact_div(sub(zs[i], zs[j]), gij)
            assert cij != (0, 0)
            c[(i, j)] = cij
            assert divides((1, 1), cij)  # (R)
    # (I) exact identity, integer form
    lhs = 1
    for (i, j) in c:
        lhs *= norm(sub(zs[i], zs[j])) ** 4
    rhs = N ** (M * (M - 2))
    for p, e in fac.items():
        expo = 0
        for tau in range(1, e + 1):
            S = sum(1 for i in range(M) if lev[p][i] < tau)
            expo += (2 * S - M) ** 2
        rhs *= p ** expo
    for cij in c.values():
        rhs *= norm(cij) ** 4
    assert lhs == rhs
    stats['identity'] += 1
    # per-prime lower bounds
    lower_total = comb(M, 2) * log(2) / 2  # ramified prime, depth one
    sum_log_c = sum(log(norm(cij)) / 2 for cij in c.values())
    cut_slack_half = 0.0
    for q in primes_small:
        if q == 2:
            continue
        if q % 4 == 3:
            for a in (1, 2):
                m = q ** a
                cnt = 0
                classes = set()
                for (i, j), cij in c.items():
                    d1 = divides((m, 0), cij)
                    d2 = residue(zs[i], m) == residue(zs[j], m)
                    assert d1 == d2  # (Q)
                    cnt += d1
                for z in zs:
                    classes.add(residue(z, m))
                if m <= 50:
                    assert class_count(N, m) == (q + 1) * q ** (a - 1)
                assert len(classes) <= (q + 1) * q ** (a - 1)
                assert cnt >= E(M, (q + 1) * q ** (a - 1))  # (L)
                lower_total += E(M, (q + 1) * q ** (a - 1)) * log(q)
                stats['inert_pairs'] += comb(M, 2)
        else:
            pi_ = split_prime_factor(q)
            if q not in fac:
                for a in (1, 2):
                    m = q ** a
                    cnt = 0
                    for (i, j), cij in c.items():
                        d1 = divides((m, 0), cij)
                        d2 = residue(zs[i], m) == residue(zs[j], m)
                        assert d1 == d2  # (S1)
                        if a == 1:
                            dpi = divides(pi_, sub(zs[i], zs[j]))
                            assert dpi == d1
                        cnt += d1
                    if m <= 50:
                        assert class_count(N, m) == (q - 1) * q ** (a - 1)
                    assert len({residue(z, m) for z in zs}) <= (q - 1) * q ** (a - 1)
                    assert cnt >= E(M, (q - 1) * q ** (a - 1))
                    lower_total += E(M, (q - 1) * q ** (a - 1)) * log(q)
                    stats['split_free_pairs'] += comb(M, 2)
            else:
                e = fac[q]
                a = lev[q]
                Np = N // q ** e
                r = [exact_div(zs[i], mul(gpow(pi_, a[i]), gpow(conj(pi_), e - a[i])))
                     for i in range(M)]
                assert all(norm(x) == Np for x in r)
                cnt = 0
                for (i, j), cij in c.items():
                    if a[i] == a[j]:
                        d1 = divides((q, 0), cij)
                        d2 = residue(r[i], q) == residue(r[j], q)
                        assert d1 == d2  # (S2)
                        cnt += d1
                    else:
                        assert not divides(pi_, cij) and not divides(conj(pi_), cij)
                        stats['diff_level_pairs'] += 1
                levels = [sum(1 for x in a if x == t) for t in range(e + 1)]
                for t in range(e + 1):
                    cls = {residue(r[i], q) for i in range(M) if a[i] == t}
                    assert len(cls) <= q - 1
                lb = sum(E(n, q - 1) for n in levels)
                assert cnt >= lb
                lower_total += lb * log(q)
                for tau in range(1, e + 1):
                    S = sum(1 for x in a if x < tau)
                    cut_slack_half += 0.5 * log(q) * (S - M / 2) ** 2
                stats['split_div_pairs'] += comb(M, 2)
    assert sum_log_c >= lower_total - 1e-9
    stats['lower_ok'] += 1
    # exact finite inequality (5.1) with the actual arc: L = R * angular span
    R = sqrt(N)
    angs = sorted(atan2(z[1], z[0]) for z in zs)
    # smallest arc containing all points = 2 pi - largest circular gap
    gaps = [angs[k + 1] - angs[k] for k in range(M - 1)] + [angs[0] + 2 * PI - angs[-1]]
    L = R * (2 * PI - max(gaps))
    Ceff = L / sqrt(R)
    # exact identity in log form, then chord <= arc
    full_slack_half = 0.0
    for p, e in fac.items():
        for tau in range(1, e + 1):
            S = sum(1 for x in lev[p] if x < tau)
            full_slack_half += 0.5 * log(p) * (S - M / 2) ** 2
    assert cut_slack_half <= full_slack_half + 1e-12
    lhs_log = sum_log_c + full_slack_half
    rhs_exact = sum(log(norm(sub(zs[i], zs[j]))) / 2 for (i, j) in c) - M * (M - 2) / 4 * log(R)
    assert abs(lhs_log - rhs_exact) < 1e-6 * max(1, abs(rhs_exact))
    assert lower_total + cut_slack_half <= M / 4 * log(R) + comb(M, 2) * log(Ceff) + 1e-9
    stats['finite_ineq_ok'] += 1


def main():
    random.seed(12345)
    norms = [
        5 ** 3 * 13 ** 2 * 17,                 # repeated split primes
        5 ** 2 * 13 * 17 * 29,
        5 * 13 * 17 * 29 * 37,
        5 ** 4 * 13 ** 2,
        2 * 9 * 5 ** 2 * 13 * 17,              # ramified + inert factor (normalisation)
        4 * 49 * 5 * 13 ** 2 * 29,
        5 ** 2 * 13 ** 2 * 17 ** 2,
        3 ** 2 * 5 ** 3 * 13 * 37 * 41,
    ]
    primes_small = [q for q in range(2, 60) if is_prime_int(q)]
    stats = dict(identity=0, inert_pairs=0, split_free_pairs=0, split_div_pairs=0,
                 diff_level_pairs=0, lower_ok=0, finite_ineq_ok=0, configs=0)
    for N0 in norms:
        pts = lattice_points(N0)
        pts.sort(key=lambda z: atan2(z[1], z[0]))
        n = len(pts)
        configs = []
        for M in range(2, min(9, n) + 1):
            for s in range(n):
                configs.append([pts[(s + k) % n] for k in range(M)])  # arcs
        for _ in range(60):
            M = random.randint(2, min(12, n))
            configs.append(random.sample(pts, M))
        for conf in configs:
            check_config(conf, primes_small, stats)
            stats['configs'] += 1
        print(f"N0={N0:>12d}  points={n:4d}  configs checked so far={stats['configs']}")
    print(stats)
    print("ALL IDENTITY/COLLISION CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
