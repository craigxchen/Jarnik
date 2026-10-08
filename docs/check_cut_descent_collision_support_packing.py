"""Exact algebra checks for cut-descent collision supports."""

from fractions import Fraction
from itertools import product

from check_ordered_affine_circuit_conductor_index import (
    gconj,
    gdiv_exact,
    gmul,
    gnorm,
    gpower,
)


def beta(alphas, history):
    out = (1, 0)
    for a, sign in zip(alphas, history):
        out = gmul(out, a if sign == 1 else gconj(a))
    return out


def determinant(z0, z1, z2):
    return ((z1[0] - z0[0]) * (z2[1] - z0[1])
            - (z1[1] - z0[1]) * (z2[0] - z0[0]))


def main():
    # alpha_0 and alpha_1 share a rational prime. These full-cube histories
    # test the factor identities with shared primes; they need not be
    # realizable as one tuple's nested threshold-cut histories.
    pi = (2, 1)
    alphas = (pi, gpower(pi, 2), (3, -2))
    ks = tuple(gnorm(a) for a in alphas)
    q = ks[0] * ks[1] * ks[2]
    histories = list(product((-1, 1), repeat=3))
    betas = {h: beta(alphas, h) for h in histories}
    assert all(gnorm(b) == q for b in betas.values())
    assert len({(b[0] % 2, b[1] % 2) for b in betas.values()}) == 1

    factorizations = 0
    for h in histories:
        for k in histories:
            if h >= k:
                continue
            common = (1, 0)
            a = (1, 0)
            support_norm = 1
            for ell, alpha in enumerate(alphas):
                if h[ell] == k[ell]:
                    common = gmul(
                        common, alpha if h[ell] == 1 else gconj(alpha)
                    )
                else:
                    a = gmul(a, alpha if h[ell] == 1 else gconj(alpha))
                    support_norm *= ks[ell]
            assert betas[h] == gmul(common, a)
            assert betas[k] == gmul(common, gconj(a))
            assert gnorm(common) * support_norm == q
            if betas[h] != betas[k]:
                diff = (betas[h][0] - betas[k][0], betas[h][1] - betas[k][1])
                expected = gmul(
                    common, (a[0] - gconj(a)[0], a[1] - gconj(a)[1])
                )
                assert diff == expected and diff[0] % 2 == diff[1] % 2 == 0
            factorizations += 1

    # A collision triple: multiplication by a common w scales its affine
    # determinant by Norm(w), and beta congruence makes D_beta divisible by 4.
    triple_h = ((-1, -1, -1), (-1, -1, 1), (-1, 1, -1))
    triple_b = [betas[h] for h in triple_h]
    assert len(set(triple_b)) == 3
    db = determinant(*triple_b)
    assert db and db % 4 == 0
    w = (1, 2)
    triple_z = [gmul(b, w) for b in triple_b]
    assert determinant(*triple_z) == gnorm(w) * db

    # Exact fractional packings: two disjoint singleton supports have mass
    # two; the three pair supports on three cuts have mass 3/2.
    disjoint = [({0}, Fraction(1)), ({1}, Fraction(1))]
    assert sum(weight for _, weight in disjoint) == 2
    assert all(sum(weight for support, weight in disjoint if ell in support) <= 1
               for ell in range(3))
    triangle = [({0, 1}, Fraction(1, 2)),
                ({0, 2}, Fraction(1, 2)),
                ({1, 2}, Fraction(1, 2))]
    assert sum(weight for _, weight in triangle) == Fraction(3, 2)
    assert all(sum(weight for support, weight in triangle if ell in support) == 1
               for ell in range(3))

    # Actual one-cut fixture with two collisions and a rectangle identity.
    points = [(-107, -54), (-98, -69), (-91, -78),
              (-78, -91), (-69, -98), (-54, -107)]
    side = {0, 1, 3}
    alpha = (4, -1)
    descended = [
        gdiv_exact(z, alpha if j in side else gconj(alpha))
        for j, z in enumerate(points)
    ]
    assert all(gnorm(z) == 845 for z in descended)
    assert descended[0] == descended[4] == (-22, -19)
    assert descended[1] == descended[5] == (-19, -22)
    assert gmul(points[0], points[5]) == gmul(points[1], points[4])

    print(
        f"PASS: {factorizations} shared-prime history factorizations, "
        "collision-triple parity/scaling, two packings, and the six-row fixture"
    )


if __name__ == "__main__":
    main()
