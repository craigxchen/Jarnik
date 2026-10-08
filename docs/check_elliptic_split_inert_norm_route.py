"""Exact Gaussian-rational tests for the CM denominator construction."""
from fractions import Fraction as F
from math import gcd, isqrt


def C(x, y=0):
    return F(x), F(y)


def plus(a, b):
    return a[0]+b[0], a[1]+b[1]


def neg(a):
    return -a[0], -a[1]


def minus(a, b):
    return plus(a, neg(b))


def times(a, b):
    return a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]


def conj(a):
    return a[0], -a[1]


def norm(a):
    return a[0]*a[0]+a[1]*a[1]


def divide(a, b):
    n = norm(b)
    u = times(a, conj(b))
    return F(u[0])/n, F(u[1])/n


def add(p, q):
    if p is None:
        return q
    if q is None:
        return p
    x, y = p
    u, v = q
    if x == u and y == neg(v):
        return None
    slope = (divide(minus(v, y), minus(u, x)) if x != u
             else divide(minus(times(C(3), times(x, x)), C(2)), times(C(2), y)))
    xx = minus(minus(times(slope, slope), x), u)
    return xx, minus(times(slope, minus(x, xx)), y)


P = (C(2), C(2))
IP = (C(-2), C(0, 2))


def multiple(n, p=P):
    if n < 0:
        return multiple(-n, (p[0], neg(p[1])))
    answer = None
    while n:
        if n & 1:
            answer = add(answer, p)
        p = add(p, p)
        n //= 2
    return answer


def ggcd(a, b):
    while b != (0, 0):
        u = times(a, conj(b)); d = norm(b)
        quotient = tuple((2*x+d)//(2*d) for x in u)
        a, b = b, minus(a, times(quotient, b))
    return a


def exact(a, b):
    q = divide(a, b)
    assert all(x.denominator == 1 for x in q)
    return tuple(int(x) for x in q)


def denominator_square(x):
    d = x[0].denominator*x[1].denominator//gcd(x[0].denominator, x[1].denominator)
    numerator = tuple(int(z*d) for z in x)
    return exact((d, 0), ggcd(numerator, (d, 0)))


def denominator_generator(point):
    square = denominator_square(point[0])
    n = isqrt(int(norm(square)))
    assert n*n == norm(square)
    # The denominator ideal is a square; one associate is a literal square.
    for unit in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        a, b = times(unit, square)
        if (n+a) % 2:
            continue
        u = isqrt(int((n+a)//2))
        if u*u != (n+a)//2:
            continue
        if u:
            if b % (2*u):
                continue
            candidate = (u, int(b//(2*u)))
        else:
            v = isqrt(int((n-a)//2))
            candidate = (0, v)
        if times(candidate, candidate) == (a, b):
            assert norm(candidate) == n
            return candidate
    raise AssertionError("denominator ideal not square")


def main():
    counts = {p: 1+sum((y*y-x*x*x+2*x) % p == 0 for x in range(p) for y in range(p))
              for p in (3, 5)}
    assert counts == {3: 4, 5: 10}
    assert multiple(2) == (C(F(9, 4)), C(F(-21, 8)))
    shifts = (0, 1, 3, 7, -1, -3, -7, -4, -8, -10)
    cache = {}
    for n in range(-10, 31):
        r = multiple(n)
        q = add(r, IP)
        assert times(q[1], q[1]) == minus(times(times(q[0], q[0]), q[0]), times(C(2), q[0]))
        assert add(q, (conj(q[0]), conj(q[1]))) == multiple(2*n)
        if r is not None:
            x, y = r[0][0], r[1][0]
            expected = C(-2*(x-2)*(x+1)/(x+2)**2, -4*y/(x+2)**2)
            assert q[0] == expected
            assert norm(q[0]) == 4*(x-1)**2/(x+2)**2
        b = denominator_generator(q)
        # A conjugate gcd of norm at most4 proves absence of every odd inert prime.
        common_norm = int(norm(ggcd(b, conj(b))))
        assert common_norm in (1, 2, 4)
        cache[n] = (q, b)
    for N in range(1, 13):
        blocks = [cache[N-k][1] for k in shifts]
        for j in range(len(shifts)):
            for l in range(j):
                # Fixed-difference formal-depth bounds, squared to avoid unit choices.
                fixed1 = denominator_generator(multiple(shifts[l]-shifts[j]))
                fixed2 = denominator_generator(add(multiple(shifts[l]-shifts[j]), multiple(2, IP)))
                exact(fixed1, ggcd(blocks[j], blocks[l]))
                exact(fixed2, ggcd(blocks[j], conj(blocks[l])))
    q2 = cache[2][0]
    assert q2 == (C(F(-26, 289), F(168, 289)), C(F(4452, 4913), F(-3646, 4913)))
    p4 = multiple(4)
    assert p4 == (C(F(12769, 84**2)), C(F(900271, 84**3)))
    assert all(z.denominator % p for coord in q2 for z in coord for p in (3, 7))
    assert all(p4[0][0].denominator % p == 0 for p in (3, 7))
    print("Good reduction counts and exact CM addition/trace identities pass.")
    print("41 Gaussian denominator ideals are conjugate-coprime away from2, with ramified exponent<=2.")
    print("Ten shifted ideals: all fixed same/opposite-orientation gcd bounds pass in12 windows.")
    print("Exact rational trace introduces inert denominator primes3 and7; no rational-family counterexample claimed.")


if __name__ == "__main__":
    main()
