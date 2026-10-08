"""Finite exact stress test for the universal F2 bilinear isotropic bound."""

from itertools import product
from random import Random


def parity(n):
    return bin(n).count("1") & 1


def images(columns):
    size = 1 << len(columns)
    out = [0] * size
    for v in range(1, size):
        low = v & -v
        out[v] = out[v ^ low] ^ columns[low.bit_length() - 1]
    return out


def isotropic_dimension_at_least(columns, target):
    if target == 0:
        return True
    im = images(columns)
    size = len(im)
    iso = [u for u in range(1, size) if not parity(u & im[u])]
    if target == 1:
        return bool(iso)

    adjacent = {
        u: {
            v for v in iso
            if not parity(u & im[v]) and not parity(v & im[u])
        }
        for u in iso
    }
    for i, u in enumerate(iso):
        for v in iso[i + 1:]:
            if v not in adjacent[u]:
                continue
            if target == 2:
                return True
            if (adjacent[u] & adjacent[v]) - {u, v, u ^ v}:
                return True
    return False


def main():
    # All 2^16 forms on F2^4 have a nonzero B-isotropic vector.
    for columns in product(range(16), repeat=4):
        assert isotropic_dimension_at_least(columns, 1)
    print("all 65536 four-dimensional forms: bound passed")

    # This invertible L has no B-isotropic plane; exact half fails.
    bad = (1, 2, 4, 9)
    im = images(bad)
    assert len(set(im)) == 16
    assert [u for u in range(1, 16)
            if not parity(u & im[u])] == [3, 5, 6, 10, 11, 12, 13]
    assert not isotropic_dimension_at_least(bad, 2)
    print("invertible four-dimensional half-bound counterexample: passed")

    rng = Random(20260915)
    for n, target, samples in ((6, 2, 100), (7, 2, 100), (8, 3, 100)):
        for _ in range(samples):
            columns = tuple(rng.randrange(1 << n) for _ in range(n))
            assert isotropic_dimension_at_least(columns, target)
        print("dimension %d, target %d, %d deterministic forms: passed" %
              (n, target, samples))


if __name__ == "__main__":
    main()
