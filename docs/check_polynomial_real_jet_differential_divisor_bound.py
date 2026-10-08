"""Exact checks for the restricted real-jet differential divisor argument."""
from fractions import Fraction
from random import Random


def trim(a):
    a = list(a)
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return tuple(a)


def add(a, b):
    return trim([(a[i] if i < len(a) else 0) +
                 (b[i] if i < len(b) else 0)
                 for i in range(max(len(a), len(b)))])


def scale(a, n):
    return trim([n * x for x in a])


def mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return trim(out)


def sq(a):
    return mul(a, a)


def derivative(a):
    return trim([i * a[i] for i in range(1, len(a))] or [0])


def degree(a):
    return len(trim(a)) - 1


def remainder(a, b):
    a, b = trim(tuple(map(Fraction, a))), trim(tuple(map(Fraction, b)))
    assert b != (0,)
    while a != (0,) and len(a) >= len(b):
        power, coefficient = len(a) - len(b), a[-1] / b[-1]
        a = add(a, (0,)*power + scale(b, -coefficient))
    return a


def differential_polynomial(H, P):
    return add(scale(mul(sq(derivative(H)), P), 4),
               scale(sq(derivative(P)), -1))


def check_identity(H, y, c):
    S = add(add(sq(y), scale(H, 2*c)), (c*c,))
    P = add(sq(H), S)
    hp, sp, yp = derivative(H), derivative(S), derivative(y)
    K = add(add(scale(mul(sq(hp), S), 4),
                scale(mul(mul(H, hp), sp), -4)), scale(sq(sp), -1))
    assert K == add(scale(mul(sq(hp), P), 4), scale(sq(derivative(P)), -1))
    rhs = add(add(scale(mul(sq(y), sq(hp)), 4),
                  scale(mul(mul(add(H, (c,)), hp), mul(y, yp)), -8)),
              scale(mul(sq(y), sq(yp)), -4))
    assert K == rhs
    assert degree(K) <= 3*degree(H)-2


def main():
    rng = Random(90630)
    for d in range(1, 15):
        for _ in range(30):
            H = tuple(rng.randint(-7, 7) for _ in range(2*d)) + (1,)
            y = tuple(rng.randint(-7, 7) for _ in range(d+1))
            check_identity(H, y, rng.randint(-9, 9))
        # (T^d +/- i)(T^d +/- 2i): four distinct equal-norm rows.
        t = (0,)*d + (1,)
        rows = []
        for a in (-1, 1):
            for b in (-2, 2):
                rows.append((add(sq(t), (-a*b,)), scale(t, a+b)))
        norms = [add(sq(x), sq(y)) for x, y in rows]
        assert len(set(norms)) == 1 and len(set(rows)) == 4
        for x, y in rows:
            assert degree(add(x, scale(rows[0][0], -1))) <= d
            assert degree(y) <= d
        representatives = dict(rows)
        assert len(representatives) == 2
        divisor = (1,)
        for y in representatives.values():
            divisor = mul(divisor, y)
        K = differential_polynomial(rows[0][0], norms[0])
        assert K != (0,) and remainder(K, divisor) == (0,)
        # The zero-K alternative admits the familiar three-row example.
        H = add(sq(t), (-1,))
        P = sq(add(H, (2,)))
        zero_rows = [(H, scale(t, 2)), (H, scale(t, -2)), (add(H, (2,)), (0,))]
        assert len(set(zero_rows)) == 3
        assert all(add(sq(x), sq(y)) == P for x, y in zero_rows)
        assert differential_polynomial(H, P) == (0,)
    for m in range(1, 100):
        assert 2*m-1 >= m
        for n in range(m, 100):
            assert min(2*n+4*m-2, 2*n+2*m-2, 4*n-2) >= m+n
    print('Passed: 420 differential identities, 14 four-row product divisors, '
          '14 zero-K three-row fixtures, multiplicity bounds.')


if __name__ == '__main__':
    main()
