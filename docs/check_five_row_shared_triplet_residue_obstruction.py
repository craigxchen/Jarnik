"""Exact shared-triplet/core/reflection identities on literal Gaussian rows.

The fixture has no controlled endpoint arc or bounded primitive residues.
It tests the algebraic identities used under those separate hypotheses.
"""

from itertools import combinations
from math import isqrt


def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    for k in range(3, isqrt(n) + 1, 2):
        if n % k == 0:
            return False
    return True


def split_primes(count):
    result = []
    candidate = 100_001
    while len(result) < count:
        if candidate % 4 == 1 and is_prime(candidate):
            result.append(candidate)
        candidate += 4
    return result


def gaussian_ab(p):
    for a in range(1, isqrt(p) + 1):
        b2 = p - a * a
        if b2 <= 0:
            break
        b = isqrt(b2)
        if b > 0 and b * b == b2:
            return a, b
    raise AssertionError(p)


def mul(z, w):
    return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]


def add(z, w):
    return z[0] + w[0], z[1] + w[1]


def neg(z):
    return -z[0], -z[1]


def conj(z):
    return z[0], -z[1]


def scale(k, z):
    return k * z[0], k * z[1]


def norm(z):
    return z[0] * z[0] + z[1] * z[1]


def imag_mul_conj(z, w):
    return z[1] * w[0] - z[0] * w[1]


def product(values):
    result = (1, 0)
    for value in values:
        result = mul(result, value)
    return result


def int_product(values):
    result = 1
    for value in values:
        result *= value
    return result


def membership(inside, a, b, c):
    return tuple(int(i in inside) for i in (a, b, c))


def exponent_a(bits):
    a, b, c = bits
    edge = a * c
    return edge + b, edge


def exponent_b(bits):
    a, b, c = bits
    edge = b * c
    return edge + a, edge


def exponent_c(bits):
    if bits in [(1, 1, 0), (1, 0, 1), (0, 1, 1)]:
        return 1, 0
    if bits == (1, 1, 1):
        return 2, 1
    return 0, 0


def exponent_res_a(bits):
    if bits == (0, 1, 0):
        return 1, 0
    if bits == (1, 0, 1):
        return 0, 1
    return 0, 0


def exponent_res_b(bits):
    if bits == (1, 0, 0):
        return 1, 0
    if bits == (0, 1, 1):
        return 0, 1
    return 0, 0


def main():
    labels = tuple(range(5))
    cuts = [frozenset(i for i in labels if mask & (1 << i))
            for mask in range(32)]
    primes = split_primes(len(cuts))
    blocks = {cut: gaussian_ab(p) for cut, p in zip(cuts, primes)}
    n = {cut: p for cut, p in zip(cuts, primes)}
    assert all(norm(blocks[cut]) == n[cut] for cut in cuts)
    a, b, c = 0, 1, 2
    inside = frozenset({3, 4})
    assert inside not in [cut for cut in cuts if cut & {a, b, c}]

    # Literal per-cut exponent equality, including full and empty cuts.
    for cut in cuts:
        bits = membership(cut, a, b, c)
        assert exponent_a(bits) == tuple(
            x + y for x, y in zip(exponent_c(bits), exponent_res_a(bits)))
        assert exponent_b(bits) == tuple(
            x + y for x, y in zip(exponent_c(bits), exponent_res_b(bits)))
    assert sum(exponent_res_a(membership(cut, a, b, c)) != (0, 0)
               for cut in cuts) == 8
    assert sum(exponent_res_b(membership(cut, a, b, c)) != (0, 0)
               for cut in cuts) == 8

    # Every inside pair has two private oriented residual signatures;
    # no monomial is shared by two different boundary restrictions.
    signatures = []
    for pair in combinations(labels, 2):
        outside = tuple(i for i in labels if i not in pair)
        aa, bb, cc = outside
        for mode in [exponent_res_a, exponent_res_b]:
            signature = tuple((tuple(sorted(cut)), mode(membership(
                cut, aa, bb, cc))) for cut in cuts
                if mode(membership(cut, aa, bb, cc)) != (0, 0))
            assert len(signature) == 8
            signatures.append(signature)
    assert len(signatures) == len(set(signatures)) == 20

    rows = {i: product(blocks[cut] for cut in cuts if i in cut)
            for i in labels}
    z = product(blocks[cut] for cut in cuts)
    reflected = {i: product(blocks[cut] for cut in cuts if i not in cut)
                 for i in labels}
    assert all(mul(rows[i], reflected[i]) == z for i in labels)

    d = {(i, j): imag_mul_conj(rows[i], rows[j])
         for i, j in combinations(labels, 2)}
    gcd_norm = {(i, j): int_product(n[cut] for cut in cuts
                                    if i in cut and j in cut)
                for i, j in combinations(labels, 2)}
    t = {}
    for edge, value in d.items():
        assert value != 0 and value % gcd_norm[edge] == 0
        t[edge] = value // gcd_norm[edge]

    def d_signed(i, j):
        return d[(i, j)] if i < j else -d[(j, i)]

    assert add(add(scale(d_signed(a, b), rows[c]),
                   scale(d_signed(b, c), rows[a])),
               scale(-d_signed(a, c), rows[b])) == (0, 0)

    g_ac = gcd_norm[(a, c)]
    g_bc = gcd_norm[(b, c)]
    c_core = product(
        product([blocks[cut]] if membership(cut, a, b, c)
                in [(1, 1, 0), (1, 0, 1), (0, 1, 1)]
                else [blocks[cut], blocks[cut], conj(blocks[cut])])
        for cut in cuts if membership(cut, a, b, c)
        in [(1, 1, 0), (1, 0, 1), (0, 1, 1), (1, 1, 1)])
    a_res = product(blocks[cut] if membership(cut, a, b, c) == (0, 1, 0)
                    else conj(blocks[cut])
                    for cut in cuts if membership(cut, a, b, c)
                    in [(0, 1, 0), (1, 0, 1)])
    b_res = product(blocks[cut] if membership(cut, a, b, c) == (1, 0, 0)
                    else conj(blocks[cut])
                    for cut in cuts if membership(cut, a, b, c)
                    in [(1, 0, 0), (0, 1, 1)])
    assert scale(g_ac, rows[b]) == mul(c_core, a_res)
    assert scale(g_bc, rows[a]) == mul(c_core, b_res)
    assert norm(c_core) == gcd_norm[(a, b)] * g_ac * g_bc
    assert inside not in [cut for cut in cuts if exponent_c(
        membership(cut, a, b, c)) != (0, 0)]
    assert exponent_c(membership(frozenset({a, b, c}), a, b, c)) == (2, 1)

    alpha, beta = 2, 3
    v = add(scale(alpha * d_signed(a, c), rows[b]),
            scale(beta * d_signed(b, c), rows[a]))
    res = add(scale(alpha * t[(a, c)], a_res),
              scale(beta * t[(b, c)], b_res))
    assert v == mul(c_core, res)
    u_vec = scale(t[(a, c)], a_res)
    v_vec = scale(t[(b, c)], b_res)
    delta = imag_mul_conj(u_vec, v_vec)
    assert delta == -t[(a, b)] * t[(a, c)] * t[(b, c)] != 0
    real_dot = u_vec[0] * v_vec[0] + u_vec[1] * v_vec[1]
    assert (2 * real_dot) ** 2 - 4 * norm(u_vec) * norm(v_vec) == -4 * delta ** 2
    assert res != (0, 0)
    for shear in [-3, -1, 0, 1, 4]:
        sheared = add(v_vec, scale(shear, u_vec))
        assert norm(sheared) * norm(u_vec) == (
            real_dot + shear * norm(u_vec)) ** 2 + delta ** 2
    assert v == add(scale(alpha * d_signed(a, b), rows[c]),
                    scale(alpha + beta, scale(d_signed(b, c), rows[a])))

    d_ref = {(i, j): imag_mul_conj(reflected[i], reflected[j])
             for i, j in combinations(labels, 2)}
    gcd_norm_ref = {(i, j): int_product(n[cut] for cut in cuts
                                        if i not in cut and j not in cut)
                    for i, j in combinations(labels, 2)}
    for edge, value in d_ref.items():
        assert value % gcd_norm_ref[edge] == 0
        assert value // gcd_norm_ref[edge] == -t[edge]
    v_ref = add(scale(alpha * d_ref[(a, c)], reflected[b]),
                scale(beta * d_ref[(b, c)], reflected[a]))
    denominator = norm(rows[a]) * norm(rows[b]) * norm(rows[c])
    assert scale(denominator, v_ref) == neg(scale(norm(z), mul(z, conj(v))))

    print("PASS: 32 exact cut exponent identities; twenty private eight-block "
          "residual signatures;")
    print("five literal Gaussian rows, ten integer primitive pair residues, "
          "small-discriminant Gram and complement-reflection identities.")
    print("Fixture residue heights and endpoint angles are uncontrolled.")


if __name__ == "__main__":
    main()
