"""Exact Hadamard character and Gaussian product checks.

Matveev and the endpoint obtuse Bessel bound are cited theorem inputs,
not numerically tested here. The literal rows have uncontrolled angles.
"""

from itertools import combinations
from math import atan2, gcd, log, pi, sqrt


K = 2**40


def hadamard(order):
    rows = [[1]]
    while len(rows) < order:
        rows = [row + row for row in rows] + [row + [-x for x in row]
                                               for row in rows]
    assert len(rows) == order
    assert all(sum(rows[i][j] * rows[i][h] for i in range(order))
               == (order if j == h else 0)
               for j in range(order) for h in range(order))
    return rows


def mul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def conj(a):
    return (a[0], -a[1])


def power(a, n):
    value = (1, 0)
    for _ in range(n):
        value = mul(value, a)
    return value


def product(values):
    result = (1, 0)
    for value in values:
        result = mul(result, value)
    return result


def norm(a):
    return a[0]**2 + a[1]**2


def exact_character_fixture():
    M = 8
    H = hadamard(M)
    core = [row[1:] for row in H]
    # A nonlinear sign column, distinct from every old Walsh column.
    extra = [(-1)**(((i >> 0) & 1) * ((i >> 1) & 1))
             for i in range(M)]
    assert all(extra != [core[i][j] for i in range(M)]
               and extra != [-core[i][j] for i in range(M)]
               for j in range(M - 1))
    w = [sum(core[i][j] * extra[i] for i in range(M))
         for j in range(M - 1)]
    collision = next((j, h) for j, h in combinations(range(M - 1), 2)
                     if w[j] == w[h])
    v = [0] * (M - 1)
    v[collision[0]], v[collision[1]] = 1, -1
    ell = sum(map(abs, v))
    lam = [sum(core[i][j] * v[j] for j in range(M - 1))
           for i in range(M)]
    assert sum(lam) == 0
    assert [sum(lam[i] * core[i][j] for i in range(M))
            for j in range(M - 1)] == [M * value for value in v]
    assert sum(lam[i] * extra[i] for i in range(M)) == 0
    L = sum(map(abs, lam))
    assert L <= M * sqrt(ell)

    blocks = [(2, 1), (3, 2), (4, 1), (5, 2), (6, 1),
              (5, 4), (7, 2), (6, 5)]
    assert all(norm(block) in (5, 13, 17, 29, 37, 41, 53, 61)
               for block in blocks)
    d = (3, 0)  # Complete common inert content cancels in the character.
    units = [power((0, 1), i % 4) for i in range(M)]
    rows = []
    for i in range(M):
        signs = core[i] + [extra[i]]
        factors = [blocks[j] if signs[j] == 1 else conj(blocks[j])
                   for j in range(M)]
        rows.append(product([d, units[i]] + factors))
    assert len({norm(row) for row in rows}) == 1

    Q_num = product(power(rows[i], lam[i]) for i in range(M)
                    if lam[i] > 0)
    Q_den = product(power(rows[i], -lam[i]) for i in range(M)
                    if lam[i] < 0)
    unit = product(power(units[i], lam[i]) for i in range(M)
                   if lam[i] > 0)
    unit_den = product(power(units[i], -lam[i]) for i in range(M)
                       if lam[i] < 0)
    # unit / unit_den is another Gaussian unit, computed without division.
    unit = mul(unit, conj(unit_den))
    raw = product(blocks[j] if v[j] > 0 else conj(blocks[j])
                  for j in range(M - 1) if v[j])
    t = M // 2
    assert mul(Q_num, power(conj(raw), t)) == mul(
        Q_den, mul(unit, power(raw, t)))
    N = norm(raw)
    middle = -2 * (raw[0]**2 - raw[1]**2)
    assert gcd(N, middle, N) == 1
    assert raw[0] and raw[1] and abs(raw[0]) != abs(raw[1])
    return M, ell, L


def effective_threshold_algebra():
    C = 2
    k = 1
    a = 1
    q = 4 * k * a
    ell = 2 * q
    M = 10**16
    assert M >= 16 * q
    assert M >= 32 * K * q * log(2 * 2.718281828459045 * M * sqrt(2 * q))
    B = M * (1 + C * sqrt(ell) / 3.141592653589793)
    gamma = 1 / 4 - 2 * K * ell * log(2.718281828459045 * B) / M
    assert gamma >= 1 / 8
    assert (M - 1) * log(5) > 8 * log(2 * M * sqrt(2 * q))
    # The large-row exclusion is a theorem consequence, not a claim
    # that a literal endpoint example of this order exists.


def cyclotomic_norm_fixture():
    # beta=pi/6, zeta=e^(i pi/3), degree two over Q(i), and
    # zeta+zeta^(-1)=1. Its relative norm is an exact Gaussian integer.
    gaussian_prime = (2, 1)
    conjugate = conj(gaussian_prime)
    relative_norm = (mul(gaussian_prime, gaussian_prime)[0]
                     + mul(conjugate, conjugate)[0]
                     - norm(gaussian_prime), 0)
    assert relative_norm == (1, 0)
    distance = abs(atan2(1, 2) - pi / 6)
    assert distance >= 1 / (4 * norm(gaussian_prime))


if __name__ == "__main__":
    M, ell, L = exact_character_fixture()
    effective_threshold_algebra()
    cyclotomic_norm_fixture()
    print(f"PASS: order-{M} exact Hadamard certificate, ell={ell}, L={L}, with one nonlinear extra column.")
    print("PASS: literal Gaussian signed row product with independent units and common inert content.")
    print("PASS: primitive composite alpha polynomial and effective threshold algebra at M=10^16.")
    print("PASS: exact quadratic torsion relative norm at beta=pi/6 and its effective phase lower bound.")
    print("Matveev, Bessel, and endpoint conditions are proof inputs, not fixture claims.")
