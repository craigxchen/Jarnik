"""Exact algebra and finite checks for the general quartic slack note."""

from fractions import Fraction
from itertools import combinations, product
from math import prod, sqrt
from random import Random


def paley(q):
    squares = {x * x % q for x in range(1, q)}

    def chi(x):
        x %= q
        return 0 if x == 0 else (1 if x in squares else -1)

    h = [[1] * (q + 1)] + [
        [1] + [-1 if x == j else -chi(x - j) for j in range(q)]
        for x in range(q)
    ]
    return h, chi


def check_weighted_identity():
    for z in product((-1, 1), repeat=4):
        assert 2 * abs(sum(z)) == (
            3 + sum(z[i] * z[k] for i, k in combinations(range(4), 2))
            - prod(z)
        )
    rng = Random(57483)
    for _ in range(100):
        s = [[rng.choice((-1, 1)) for _ in range(19)] for _ in range(4)]
        w = [Fraction(rng.randrange(1, 100), 7) for _ in range(19)]
        kappa = Fraction(rng.randrange(-10, 10), 3)
        total = sum(w)
        quartic = sum(w[j] * prod(s[x][j] for x in range(4))
                      for j in range(19))
        slack = {(x, y): kappa - sum(w[j] * s[x][j] * s[y][j]
                                    for j in range(19))
                 for x, y in combinations(range(4), 2)}
        height = sum(w[j] * abs(s[0][j] + s[1][j] - s[2][j] - s[3][j]) / 2
                     for j in range(19))
        rhs = ((total - quartic) / 4 - kappa / 2
               + (slack[0, 2] + slack[0, 3] + slack[1, 2] + slack[1, 3]
                  - slack[0, 1] - slack[2, 3]) / 4)
        assert height - total / 2 == rhs
        assert total - quartic == 2 * sum(
            w[j] for j in range(19) if prod(s[x][j] for x in range(4)) == -1)


def check_paley(q, b=5):
    h, chi = paley(q)
    m, r = q + 1, b * q
    assert all(sum(h[x][j] * h[y][j] for j in range(m))
               == (m if x == y else 0)
               for x in range(m) for y in range(m))
    masks = [sum((h[x][j + 1] == -1) << j for j in range(q))
             for x in range(m)]
    physical = [sum((h[x][j % q + 1] == -1) << j for j in range(r))
                ^ (1 << x) for x in range(m)]
    # Physical column x is flipped at row x, hence columns are distinct.
    assert all(r - 2 * bin(physical[x] ^ physical[y]).count("1") <= -b + 4
               for x, y in combinations(range(m), 2))
    count, max_abs, min_defect = 0, 0, r * 2
    for rows in combinations(range(m), 4):
        raw_mask = 0
        flipped_mask = 0
        for x in rows:
            raw_mask ^= masks[x]
            flipped_mask ^= physical[x]
        raw = q - 2 * bin(raw_mask).count("1")
        flipped = r - 2 * bin(flipped_mask).count("1")
        assert abs(flipped - b * raw) <= 8
        assert abs(raw) <= 3 * sqrt(q) + 4
        finite = [x - 1 for x in rows if x]
        # Only the diagonal entries differ from the pure character product.
        diagonal_error = sum(prod(h[x][j + 1] for x in rows) for j in finite)
        polynomial_sum = raw - diagonal_error
        assert abs(diagonal_error) <= len(finite)
        assert abs(polynomial_sum) <= (len(finite) - 1) * sqrt(q)
        if count < 20:
            direct = sum(prod(-chi(x - j) for x in finite) for j in range(q))
            assert polynomial_sum == direct
        max_abs = max(max_abs, abs(raw))
        min_defect = min(min_defect, r - flipped)
        count += 1
    print(f"q={q}: {count:,} quadruples, max |raw T|={max_abs}, "
          f"min flipped defect={min_defect}/{r}")


def main():
    check_weighted_identity()
    for q in (11, 19, 31, 43):
        check_paley(q)
    print("PASS: Boolean identity; 100 exact rational weighted identities; "
          "all Paley quadruples at four orders; distinct physical one-flips.")


if __name__ == "__main__":
    main()
