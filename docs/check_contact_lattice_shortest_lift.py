"""Exact checks for contact_lattice_shortest_lift.md."""
from fractions import Fraction
from itertools import combinations, product
from math import gcd


def add(a, b):
    return a[0] + b[0], a[1] + b[1]


def scale(n, a):
    return n * a[0], n * a[1]


def mul(a, b):
    return a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]


def conj(a):
    return a[0], -a[1]


def norm(a):
    return a[0] ** 2 + a[1] ** 2


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1]


def det(a, b):
    return a[0] * b[1] - a[1] * b[0]


def exact_div(a, b):
    num = mul(a, conj(b))
    den = norm(b)
    if num[0] % den or num[1] % den:
        return None
    return num[0] // den, num[1] // den


def divides(b, a):
    return exact_div(a, b) is not None


def power(a, e):
    out = (1, 0)
    for _ in range(e):
        out = mul(out, a)
    return out


def valuation(a, pi):
    e = 0
    while a != (0, 0):
        quotient = exact_div(a, pi)
        if quotient is None:
            return e
        a = quotient
        e += 1
    raise ValueError("valuation of zero")


def xgcd(a, b):
    old_r, r = a, b
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t
    if old_r == -1:
        old_r, old_s, old_t = 1, -old_s, -old_t
    assert old_r == 1
    return old_s, old_t


def reduced_partner(z):
    assert gcd(z[0], z[1]) == 1
    alpha, beta = xgcd(z[0], z[1])
    w = (-beta, alpha)
    assert det(z, w) == 1
    A = norm(z)
    # One of these shears minimizes the scalar product modulo A.
    candidates = [add(w, scale(t, z)) for t in range(-A - 1, A + 2)]
    w = min(candidates, key=lambda u: abs(dot(z, u)))
    assert 2 * abs(dot(z, w)) <= A
    return w


def check_clean_period():
    z = (2, 1)
    w0 = reduced_partner(z)
    J = power((3, 2), 2)  # Norm 13^2 and coprime to z.
    b, r = 2, -3
    period = norm(J)
    solutions = [
        t for t in range(period)
        if divides(J, add(scale(r + b * t, z), scale(b, w0)))
    ]
    assert len(solutions) == 1
    return z, w0, solutions[0], period


def check_overlap_period():
    pi = (2, 1)
    p = norm(pi)
    a, e = 1, 5
    z = power(pi, a)
    G = power(pi, e)
    Q4 = G
    w0 = reduced_partner(z)
    actual_t = 3
    w = add(w0, scale(actual_t, z))
    s4 = det(z, Q4)
    r4 = det(Q4, w)
    assert Q4 == add(scale(r4, z), scale(s4, w))

    vp_s = valuation((s4, 0), pi)
    c_exp = min(e, vp_s + valuation(z, pi))
    assert vp_s == a
    assert c_exp == min(e, 2 * a)
    C4 = power(pi, c_exp)
    reduced_G = exact_div(G, C4)
    assert reduced_G is not None
    constant = add(scale(r4, z), scale(s4, w0))
    assert divides(C4, constant)

    period = norm(reduced_G)
    solutions = [
        t for t in range(period)
        if divides(G, add(constant, scale(s4 * t, z)))
    ]
    assert solutions == [actual_t % period]
    assert period == p ** (e - 2 * a)
    return period


def check_shortest_formula(z, w0, tau, M):
    A = norm(z)
    e = dot(z, w0)
    theta = Fraction(e, A)
    candidates = [(tau + k * M, abs(Fraction(tau + k * M) + theta))
                  for k in range(-3, 4)]
    t, rho = min(candidates, key=lambda item: item[1])
    w = add(w0, scale(t, z))
    exact = Fraction(A) * rho * rho + Fraction(1, A)
    assert norm(w) == exact
    assert norm(w) == min(
        norm(add(w0, scale(tau + k * M, z))) for k in range(-3, 4)
    )


def check_median_changes():
    triples = list(combinations(range(5), 3))

    def median(bits, triple):
        return int(sum(bits[i] for i in triple) >= 2)

    for first, second in combinations(triples, 2):
        disagreements = sum(
            median(bits, first) != median(bits, second)
            for bits in product((0, 1), repeat=5)
        )
        assert disagreements == (8 if len(set(first) & set(second)) == 2 else 12)


def main():
    z, w0, tau, period = check_clean_period()
    overlap_period = check_overlap_period()
    check_shortest_formula(z, w0, tau, period * overlap_period)
    check_median_changes()
    print("Clean and correction-overlap shear periods pass exactly.")
    print("The shifted least-CRT representative gives the exact shortest partner.")
    print("All changes of five-row median triples have 8 or 12 Boolean disagreements.")


if __name__ == "__main__":
    main()
