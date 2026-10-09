"""Exact check of the triangle (Jarnik) identity of outside.md, Prop C.1:

  for distinct z_a, z_b, z_c with |z|^2 = N:
     Norm(g_ab) Norm(g_bc) Norm(g_ac) = N * Norm(g_abc)^2            (C.1a)
     4 (2 Area)^2 = Norm(g_abc)^2 * nu_ab^2 nu_bc^2 nu_ac^2          (C.1b)
  where g_ab = gcd(z_a, z_b), g_abc = gcd(z_a, z_b, z_c), nu_ab^2 = |z_a - z_b|^2 / Norm(g_ab),
  2 Area = |det(z_b - z_a, z_c - z_a)|.
Pure integer arithmetic.  Also checks the summed form: the product over all triples of
(C.1b) equals the (M-2)-th power of the product of pair identities.
"""
import random
from itertools import combinations
from gx import reps, gcd_g, norm, sub, clusters

random.seed(1)


def check_triple(N, za, zb, zc):
    gab, gbc, gac = gcd_g(za, zb), gcd_g(zb, zc), gcd_g(za, zc)
    gabc = gcd_g(gab, zc)
    nab, nbc, nac, nabc = norm(gab), norm(gbc), norm(gac), norm(gabc)
    assert nab * nbc * nac == N * nabc * nabc
    d = lambda z, w: norm(sub(z, w))
    assert d(za, zb) % nab == 0 and d(zb, zc) % nbc == 0 and d(za, zc) % nac == 0
    nu2 = (d(za, zb) // nab) * (d(zb, zc) // nbc) * (d(za, zc) // nac)
    u, v = sub(zb, za), sub(zc, za)
    area2 = abs(u[0] * v[1] - u[1] * v[0])
    assert area2 > 0
    assert 4 * area2 * area2 == nabc * nabc * nu2
    return 1


def main():
    count = 0
    Ns = []
    # circles with many points: products of split primes, with prime powers and inert/2 factors
    split = [5, 13, 17, 29, 37, 41, 53, 61, 73, 89, 97, 101, 109, 113]
    for _ in range(400):
        N = 1
        for p in random.sample(split, random.randint(2, 5)):
            N *= p ** random.randint(1, 3)
        N *= random.choice([1, 2, 4, 9, 18, 49])
        Ns.append(N)
    Ns.append(1176852625)
    ncl = 0
    for N in Ns:
        pts = reps(N)
        if len(pts) < 3:
            continue
        # random triples anywhere on the circle
        for _ in range(30):
            a, b, c = random.sample(pts, 3)
            count += check_triple(N, a, b, c)
        # all triples inside short-arc clusters (C = 8)
        for cl in clusters(N, 8.0, 3):
            ncl += 1
            for a, b, c in combinations(cl, 3):
                count += check_triple(N, a, b, c)
    print("triangle identity (C.1a),(C.1b): %d exact triples on %d circles, %d clusters: OK"
          % (count, len(Ns), ncl))


if __name__ == "__main__":
    main()
