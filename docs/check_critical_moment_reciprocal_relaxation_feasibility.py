"""Exact certificates for the bare reciprocal relaxation, not admissible root order."""

from fractions import Fraction as F
from itertools import combinations
from random import Random


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def run():
    cube = [tuple(1 if mask & (1 << i) else -1 for i in range(8))
            for mask in range(256)]
    one = (1,) * 8
    balanced = (1,) * 4 + (-1,) * 4
    base = [v for v in cube if v not in (one, balanced)]
    assert len(base) == 254
    for i in range(8):
        assert sum(v[i] for v in base) == -1 - balanced[i]
        for j in range(8):
            count = sum(v[i] * v[j] for v in base)
            assert count == 256 * (i == j) - 1 - balanced[i] * balanced[j]
            assert sum(v[i] * v[j] for v in cube) == 256 * (i == j)
            assert count + 6 * sum(v[i] * v[j] for v in cube) == (
                1792 * (i == j) - 1 - balanced[i] * balanced[j])

    a_star = [1 + F(dot(v, [1 + t for t in balanced]), 248) for v in base]
    assert all(F(9, 10) < a < F(11, 10) for a in a_star)
    assert all(sum(v[i] * a for v, a in zip(base, a_star)) == 0
               for i in range(8))
    assert F(32, 248) < 2  # No coordinate difference belongs to the row space.

    # One explicit rational generic kernel perturbation, independently of
    # the finite-hyperplane existence proof in the note.
    rng = Random(20260922)
    trial = [F(rng.randrange(-1000, 1001)) for _ in base]
    image = [sum(v[i] * a for v, a in zip(base, trial)) for i in range(8)]
    inverse_image = [image[i] / 256 + F(
        sum(image) + balanced[i] * dot(balanced, image), 256 * 248)
        for i in range(8)]
    kernel = [a - dot(v, inverse_image) for v, a in zip(base, trial)]
    assert all(sum(v[i] * a for v, a in zip(base, kernel)) == 0
               for i in range(8))
    scale = F(1) / (200 * max(map(abs, kernel)))
    magnitudes = [a + scale * t for a, t in zip(a_star, kernel)]
    assert isinstance(scale, F) and all(isinstance(a, F) for a in magnitudes)
    assert len(set(magnitudes)) == 254
    assert all(F(9, 10) < a < F(11, 10) for a in magnitudes)
    assert all(sum(v[i] * a for v, a in zip(base, magnitudes)) == 0
               for i in range(8))

    # Character orthogonality certifies all four Fourier equations in (2).
    characters = [()] + [(i,) for i in range(8)] + list(combinations(range(8), 2))
    for left in characters:
        for right in characters:
            difference = set(left) ^ set(right)
            total = sum(product(v[i] for i in difference) for v in cube)
            assert total == (256 if left == right else 0)

    length = F(10000)
    gamma = 2 * length ** 2
    beta = [F(i) / length for i in range(4)] + [
        2 * length + F(j) / length for j in range(4)]
    assert len(set(beta)) == 8
    target_cross = [gamma - (beta[i] - beta[k]) ** 2 / 2
                    for i in range(4) for k in range(4, 8)]
    assert max(map(abs, target_cross)) < 7
    assert F(25400, 81) < 314
    assert F(2540, 9) < 283
    lower_w = (2 * length ** 2 - 5450) / 256
    upper_t = (8 * length + 12 / length + 2264) / 256
    assert lower_w > upper_t ** 2
    assert 4 * length ** 2 - 628 - 9 / length ** 2 > 0
    assert all((beta[i] - beta[j]) ** 2 <= 9 / length ** 2
               for group in [range(4), range(4, 8)] for i, j in combinations(group, 2))
    lower_q = 2 * length ** 2 - 314
    assert F(100, 81) < 2
    assert F(2540000, 6561) < 388
    assert lower_q > max(2 * 5136, 16 * 2)
    assert lower_q ** 2 > 512 * 388
    assert F(3, 512) < F(1, 16)
    assert F(16, 512) + F(9, 128) == F(13, 128) < F(1, 9)

    for rho in [F(3, 2), F(5, 4), F(7, 5)]:
        for target_t in [F(2), F(-1), F(7, 3)]:
            p = (1 + rho + rho ** 2) / (rho * (1 + rho) * abs(target_t))
            q = rho * p
            sign = 1 if target_t > 0 else -1
            assert p < q < p + q
            assert sign * (p + q - (p + q)) == 0
            assert sign * (1 / p + 1 / q - 1 / (p + q)) == target_t
            assert 1 / p ** 2 + 1 / q ** 2 + 1 / (p + q) ** 2 == target_t ** 2
    # A rational example of the general two-root split: T=1,W=5.
    assert 2 + (-1) == 1 and 2 ** 2 + (-1) ** 2 == 5
    assert 1790 == 4 * 447 + 2
    print("PASS: positive distinct rational first-moment base and count Gram")
    print("PASS: Fourier completion, uniform strict bounds, triads, and within-group signs")
    print("PASS: finite separation from paired sign order and fourth reciprocal identity")


def product(values):
    out = 1
    for value in values:
        out *= value
    return out


if __name__ == "__main__":
    run()
