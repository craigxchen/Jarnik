"""Exact Gaussian degree-seven lifts and their generic old-core valuations.

The table is a generic divisor calculation, not an actual conductor bound.
Variables0..5 are P_i and variables6..11 their formal conjugates.
The bracket numerator is bar(P_i)P_j-P_i bar(P_j); its denominator2i
and the common denominator(2i)^3 of the lifts are units at odd core primes.
"""

from collections import defaultdict
from pathlib import Path
import json

ZERO = (0,) * 12


def add(*polynomials):
    result = defaultdict(int)
    for polynomial in polynomials:
        for monomial, coefficient in polynomial.items():
            result[monomial] += coefficient
    return {m: c for m, c in result.items() if c}


def scale(polynomial, scalar):
    return {m: c * scalar for m, c in polynomial.items() if c * scalar}


def mul(left, right):
    result = defaultdict(int)
    for m, a in left.items():
        for n, b in right.items():
            result[tuple(x + y for x, y in zip(m, n))] += a * b
    return {m: c for m, c in result.items() if c}


def product(*polynomials):
    result = {ZERO: 1}
    for polynomial in polynomials:
        result = mul(result, polynomial)
    return result


def variable(index):
    monomial = list(ZERO)
    monomial[index] = 1
    return {tuple(monomial): 1}


u = [variable(i) for i in range(6)]
v = [variable(i + 6) for i in range(6)]


def bracket(i, j):
    return add(mul(v[i - 1], u[j - 1]), scale(mul(u[i - 1], v[j - 1]), -1))


d = bracket
A = product(d(1, 2), d(3, 4), d(5, 6))
B = product(d(1, 2), d(3, 6), d(4, 5))
C = product(d(1, 4), d(2, 3), d(5, 6))
F = add(C, scale(B, -1))
H = add(C, B)
V = add(scale(A, 2), B, C)
Q = [add(mul(H, u[3]), scale(product(d(3, 4), d(5, 6), d(4, 1), u[1]), -2)),
     add(mul(F, u[4]), scale(product(d(3, 5), d(4, 6), d(5, 1), u[1]), -2)),
     add(mul(H, u[5]), scale(product(d(3, 6), d(4, 5), d(6, 1), u[1]), 2))]
for i, polynomial, first, second in zip(range(4, 7), Q, [H, F, H], [V, V, F]):
    assert mul(d(1, 2), polynomial) == add(product(d(i, 2), first, u[0]),
                                           scale(product(d(i, 1), second, u[1]), -1))
    assert len(polynomial) == 16

rows = []
for mask in range(64):
    inside = [i for i in range(6) if mask >> i & 1]
    a = [min(sum(m[i] for i in inside) for m in polynomial) for polynomial in Q]
    b = [min(sum(m[i + 6] for i in inside) for m in polynomial) for polynomial in Q]
    phase = [x - y for x, y in zip(a, b)]
    fixed = [int(i in inside) for i in range(3)]
    width = max([0] + fixed + phase) - min([0] + fixed + phase)
    rows.append({"S": [i + 1 for i in inside], "pi": a, "barpi": b,
                 "phase_difference": phase, "radius_width": width})

assert sum(row["radius_width"] for row in rows) == 56
assert all(row["radius_width"] == int(bool(set(row["S"]) & {1, 2, 3})) for row in rows)
for i in range(3):
    assert sum(row["pi"][i] for row in rows) == 64
    assert sum(row["barpi"][i] for row in rows) == 32
    assert sum(min(row["pi"][i], row["barpi"][i]) for row in rows) == 32

Path(__file__).with_name("joint_circle_lift_cut_valuations.json").write_text(
    json.dumps(rows, indent=2) + "\n")
print("All three degree-seven inverse-normalization identities verified.")
print("Each Gaussian numerator has16 monomials; all64 old-cut valuations checked.")
print("Generic old-core log-radius coefficient:28; each row's rational content:32.")
