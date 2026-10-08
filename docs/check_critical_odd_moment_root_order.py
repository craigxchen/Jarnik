"""Exact rational fixtures for the critical odd-moment sign-order theorem."""
from fractions import Fraction as F


def multiply(a, b):
    out = [F(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out


def evaluate(coefficients, x):
    value = F(0)
    for c in reversed(coefficients):
        value = value*x+c
    return value


def isolate(coefficients, left, right, steps=48):
    a, b = evaluate(coefficients, left), evaluate(coefficients, right)
    assert a*b < 0
    for _ in range(steps):
        mid = (left+right)/2
        m = evaluate(coefficients, mid)
        if m == 0:
            return mid, mid
        if a*m < 0:
            right, b = mid, m
        else:
            left, a = mid, m
    return left, right


cases = roots_checked = 0
for s in range(1, 9):
    odd = [F(0), F(1)]
    for k in range(1, s+1):
        odd = multiply(odd, [-k*k, 0, 1])
    n = 2*s+1
    assert len(odd) == n+1 and odd[-1] == 1
    epsilon = F(1, 100)
    intervals = [(F(k)-epsilon, F(k)+epsilon) for k in range(-s, s+1)]
    endpoint_min = min(abs(evaluate(odd, x)) for interval in intervals for x in interval)
    assert endpoint_min > 0
    for sign in (-1, 1):
        polynomial = odd.copy()
        polynomial[0] = sign*endpoint_min/4
        power_sums = [F(n)]
        for k in range(1, n+1):
            value = -k*polynomial[n-k]
            value -= sum(polynomial[n-j]*power_sums[k-j] for j in range(1, k))
            power_sums.append(value)
        assert all(power_sums[k] == 0 for k in range(1, 2*s, 2))
        absolute_intervals = []
        for left, right in intervals:
            lo, hi = isolate(polynomial, left, right)
            assert lo*hi > 0
            if lo > 0:
                absolute_intervals.append((lo, hi, 1))
            else:
                absolute_intervals.append((-hi, -lo, -1))
        absolute_intervals.sort()
        assert all(x[1] < y[0] for x, y in zip(absolute_intervals, absolute_intervals[1:]))
        signs = [x[2] for x in absolute_intervals]
        expected = [signs[0]*(-1)**(k//2) for k in range(n)]
        assert signs == expected
        assert sum(x != y for x, y in zip(signs, signs[1:])) == s
        assert abs(sum(signs)) == 1
        cases += 1
        roots_checked += n
assert cases == 16 and roots_checked == 160
print('PASS: 16 exact rational polynomials, 160 isolated real roots; critical paired sign order and all required odd moments.')


def check_successor(polynomial, intervals, pattern, steps=48):
    n = len(polynomial)-1
    assert n % 2 == 0 and polynomial[-1] == 1
    assert polynomial[0] and polynomial[1]
    assert all(polynomial[k] == 0 for k in range(3, n, 2))
    values = []
    for left, right in intervals:
        lo, hi = isolate(polynomial, left, right, steps)
        assert lo*hi > 0
        values.append((lo, hi, 1) if lo > 0 else (-hi, -lo, -1))
    values.sort()
    assert len(values) == n
    assert all(x[1] < y[0] for x, y in zip(values, values[1:]))
    signs = [x[2] for x in values]
    expected = [(-1)**((k+1)//2) for k in range(n)]
    if pattern == 'B':
        expected = [1, 1]+expected[:-2]
    assert signs == [signs[0]*x for x in expected]


for m in range(2, 10):
    even = [F(1)]
    for k in range(1, m+1):
        even = multiply(even, [-k*k, 0, 1])
    intervals = [(F(k)-F(1, 100), F(k)+F(1, 100))
                 for k in list(range(-m, 0))+list(range(1, m+1))]
    endpoints = [x for interval in intervals for x in interval]
    r = min(abs(evaluate(even, x))/abs(x) for x in endpoints)/4
    for sign in (-1, 1):
        polynomial = even.copy()
        polynomial[1] = sign*r
        check_successor(polynomial, intervals, 'A')

for roots in ((F(1), F(2), F(3), F(-6)),
              (F(1), F(2), F(3), F(-15, 2), F(-17, 2), F(10))):
    n = len(roots)
    assert all(sum(x**k for x in roots) == 0 for k in range(1, n-2, 2))
    ordered = sorted(roots, key=abs)
    expected = [1, 1]+[(-1)**((k+1)//2) for k in range(n-2)]
    assert [1 if x > 0 else -1 for x in ordered] == expected

# Start with t*(F(t)+r), then perturb its simple zero root while retaining
# all twelve real roots. The even part gains exactly one negative x-root.
g = [F(1)]
for k in range(1, 6):
    g = multiply(g, [-k*k, 0, 1])
odd = [F(0)]+g
intervals = [(F(k)-F(1, 100), F(k)+F(1, 100)) for k in range(-5, 6)]
r = min(abs(evaluate(odd, x)) for interval in intervals for x in interval)/4
odd[0] = r
old_intervals = [isolate(odd, a, b) for a, b in intervals]
assert all(a < b and a*b > 0 for a, b in old_intervals)
zero_radius = min(abs(x) for interval in old_intervals for x in interval)/4
all_intervals = sorted(old_intervals+[(-zero_radius, zero_radius)])
assert all(a[1] < b[0] for a, b in zip(all_intervals, all_intervals[1:]))
base = [F(0)]+odd
endpoints = [x for interval in all_intervals for x in interval]
c = min(abs(evaluate(base, x))/(4*max(F(1), abs(evaluate(g, x)))) for x in endpoints)
assert c > 0
successor = base.copy()
for i, value in enumerate(g):
    successor[i] += c*value
check_successor(successor, all_intervals, 'B', steps=96)
print('PASS: successor patterns A/B, including exact rational six-root data and an isolated degree-twelve pattern-B polynomial.')
