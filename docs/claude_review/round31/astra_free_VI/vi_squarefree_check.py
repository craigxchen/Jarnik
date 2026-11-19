"""Check the squarefree fixed-field construction. Run python3 vi_squarefree_check.py.

Binary integers encode polynomials over F2 in ascending powers of t. The F4
operations are those used in the Part V certificate.
"""

from v_fixed2_check import carryless, degree, multiply_factors, product


def residue(a, p):
    while a.bit_length() >= p.bit_length():
        a ^= p << (a.bit_length() - p.bit_length())
    return a


def gcd(a, b):
    while b:
        a, b = b, residue(a, b)
    return a


def square(a):
    return carryless(a, a)


def trace(a, p):
    a = residue(a, p)
    total = 0
    for _ in range(p.bit_length() - 1):
        total ^= a
        a = residue(square(a), p)
    assert total in (0, 1)
    return total


def primes_through(max_degree):
    primes = []
    for e in range(1, max_degree + 1):
        for candidate in range(1 << e, 1 << (e + 1)):
            if all(residue(candidate, p) for p in primes
                   if 2 * (p.bit_length() - 1) <= e):
                primes.append(candidate)
    return primes


def derivative(a):
    # In characteristic two each odd-position bit shifts down one position.
    digits = (a.bit_length() + 1) // 2
    even_positions = ((1 << (2 * digits)) - 1) // 3
    return (a >> 1) & even_positions


def check(d):
    q = 1 << d
    m = q // 2
    # The low bit of the index selects the constant, and the remaining bits
    # encode f in span(t,...,t^(d-1)).
    values = []
    for index in range(q):
        f = index & ~1
        values.append(square(square(f) ^ f) ^ (index & 1))
    assert len(set(values)) == q
    assert all(derivative(v) == 0 for v in values)
    # v is the square of f^2+f or f^2+f+1. This bounds the search degree.
    primes = [p for p in primes_through(2 * (d - 1))
              if (p.bit_length() - 1) % 2 == 0
              and any(residue(v, p) == 0 for v in values[1:])]
    sides = [[v for index, v in enumerate(values)
              if not ((index & form).bit_count() % 2)]
             for form in range(1, q, 2)]
    sides = [set(side) for side in sides]
    periods = [sum((1 ^ v).bit_length() - 1 for v in left & right)
               for j, left in enumerate(sides) for right in sides[:j]]
    # Impose Tr(T)=1+Tr(theta)=1+(deg(p)/2 mod 2) by T=t+S^2.
    s = 0
    modulus = 1
    for p in primes:
        e = p.bit_length() - 1
        goal = 1 ^ ((e // 2) % 2)
        if trace(2 ^ s, p) != goal:  # Tr(S)=Tr(S^2)
            j = next(j for j in range(e) if trace(modulus << j, p))
            s ^= modulus << j
        modulus = carryless(modulus, p)
        assert trace(2 ^ s, p) == goal
    # Make the specialization degree exceed every deg(v) and attain the
    # required coefficient margin, without changing the residues.
    target_deg_s = 2 * m * (d - 2)
    assert max(periods) == target_deg_s
    assert modulus.bit_length() - 1 <= target_deg_s
    s ^= modulus << (target_deg_s - (modulus.bit_length() - 1))
    t_value = 2 ^ square(s)
    specialization_degree = t_value.bit_length() - 1
    assert derivative(t_value) == 1
    parent = multiply_factors([v for v in values if v not in (0, 1)], t_value)
    assert parent[1] == 0 and gcd(parent[0], derivative(parent[0])) == 1
    marks = [multiply_factors(sorted(side - {0}), t_value) for side in sides]
    assert len(set(marks)) == m
    D = (m - 1) * specialization_degree
    assert degree(parent) == 2 * D
    for z in marks:
        a, b = z
        assert degree(z) == D and product(z, (a ^ b, b)) == parent
    gaps = [degree((x[0] ^ y[0], x[1] ^ y[1]))
            for j, x in enumerate(marks) for y in marks[:j]]
    assert all(g == (m // 2 - 1) * specialization_degree + coefficient
               for g, coefficient in zip(gaps, periods))
    assert 2 * max(gaps) <= D
    print(d, m, len(primes), modulus.bit_length() - 1,
          specialization_degree, 2 * D, max(gaps), len(gaps), sep="\t")


if __name__ == "__main__":
    print("d\tm\tCRT primes\tdeg(modulus)\tdeg(T)\tdeg(N)\tmax gap\tpairs")
    for dimension in range(3, 7):
        check(dimension)
