"""Exact polynomial and primitive-content checks for Mukai's quadratic map."""

from itertools import combinations, product
from math import gcd

from check_near_balanced_invariant_congruences import add, mul, scale, variable


def total(polynomials):
    result = {}
    for polynomial in polynomials:
        result = add(result, polynomial)
    return result


def times(polynomials):
    result = {(0,) * 6: 1}
    for polynomial in polynomials:
        result = mul(result, polynomial)
    return result


def quadratic(xs):
    return add(total(mul(x, x) for x in xs),
               scale(total(mul(x, y) for x, y in combinations(xs, 2)), -2))


def steiner(xs):
    q = add(mul(xs[4], xs[4]), scale(quadratic(xs[:4]), -1))
    return add(mul(q, q), scale(times(xs[:4]), -64))


SIGNS = ((1, 1, 1, 1), (1, 1, -1, -1),
         (1, -1, 1, -1), (1, -1, -1, 1))
x = [variable(i) for i in range(5)]
linear = [total(scale(xi, sign) for xi, sign in zip(x, signs))
          for signs in SIGNS]
y = [mul(form, form) for form in linear]
y.append(scale(add(quadratic(x[:4]), scale(mul(x[4], x[4]), -1)), 2))
source = steiner(x)
assert steiner(y) == scale(mul(source, add(source, scale(times(linear), 4))), 16)

checked = 0
max_content = 0
for vector in product(range(-3, 4), repeat=5):
    content = 0
    for entry in vector:
        content = gcd(content, entry)
    if content != 1:
        continue
    height = max(map(abs, vector))
    forms = [sum(sign * entry for sign, entry in zip(signs, vector))
             for signs in SIGNS]
    q = sum(entry * entry for entry in vector[:4])
    q -= 2 * sum(a * b for a, b in combinations(vector[:4], 2))
    output = [form * form for form in forms] + [2 * (q - vector[4] ** 2)]
    raw_height = max(map(abs, output))
    content = 0
    for entry in output:
        content = gcd(content, entry)
    assert 1 <= content <= 16 and content & (content - 1) == 0
    assert height * height <= 64 * raw_height
    assert raw_height <= 34 * height * height
    checked += 1
    max_content = max(max_content, content)

print("Exact identity S(map(x))=16*S(x)*(S(x)+4*product(L_i)) verified.")
print("Primitive integral windows checked:", checked)
print("Largest raw common content in the window:", max_content)
