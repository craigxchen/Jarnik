"""Exact residual spectral-factor and orientation checks; no numerical roots."""
from fractions import Fraction as F
from itertools import combinations, product


def mul(a, b):
    out = [F(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out


def ev(p, x):
    result = F(0)
    for coefficient in reversed(p):
        result = result*x+coefficient
    return result


def sign(x):
    assert x
    return 1 if x > 0 else -1


def reflect(p):
    return [(-1)**k*c for k, c in enumerate(p)]


def first_root_sign(p, intervals):
    absolute_intervals = []
    for lo, hi in intervals:
        a, b = ev(p, lo), ev(p, hi)
        assert a*b < 0
        for _ in range(64):
            mid = (lo+hi)/2
            value = ev(p, mid)
            if value == 0:
                lo = hi = mid
                break
            if a*value < 0:
                hi, b = mid, value
            else:
                lo, a = mid, value
        assert lo*hi > 0
        absolute_intervals.append((lo, hi, 1) if lo > 0 else (-hi, -lo, -1))
    absolute_intervals.sort()
    assert len(absolute_intervals) == len(p)-1
    assert all(a[1] < b[0] for a, b in zip(absolute_intervals, absolute_intervals[1:]))
    return absolute_intervals[0][2]


def check_pair(residual, s, kind, epsilon):
    d = 2*s+1
    h = d if kind == 'critical' else d+1
    assert len(residual) == h+1 and residual[-1] == 1
    common = [F(1)]
    for j in range(2*d-h):
        common = mul(common, [-100-j, 1])
    other = [(-1)**h*x for x in reflect(residual)]
    pi, pj = mul(common, residual), mul(common, other)
    assert len(pi) == 2*d+1 and len(pj) == 2*d+1
    assert mul(pi, reflect(pi)) == mul(pj, reflect(pj))
    assert pi[d+1:] == pj[d+1:]
    r = residual[0] if kind == 'critical' else residual[1]
    assert r and pi[d]-pj[d] == 2*r
    factor = common if kind == 'critical' else [F(0)]+common
    expected = [2*r*x for x in factor]+[F(0)]*(d)
    assert [x-y for x, y in zip(pi, pj)] == expected
    parity = s if kind == 'A' else s+1
    assert sign(pi[d]-pj[d]) == (-1)**parity*epsilon
    # Common positive roots are 100,...; their remaining parity orders
    # the actual pair values at every intervening positive test point.
    common_roots = list(range(100, 100+2*d-h))
    for t in [F(1, 2), F(50)]+[F(2*x+1, 2) for x in common_roots]:
        remaining = sum(x > t for x in common_roots)
        assert sign(ev(pi, t)-ev(pj, t)) == sign(2*r)*(-1)**remaining


cases = 0
for s in range(1, 7):
    for kind in ('critical', 'A'):
        if kind == 'critical':
            base = [F(0), F(1)]
            centers = list(range(-s, s+1))
            for k in range(1, s+1):
                base = mul(base, [-k*k, 0, 1])
            changed = 0
        else:
            base = [F(1)]
            centers = list(range(-s-1, 0))+list(range(1, s+2))
            for k in range(1, s+2):
                base = mul(base, [-k*k, 0, 1])
            changed = 1
        intervals = [(F(k)-F(1, 100), F(k)+F(1, 100)) for k in centers]
        size = min(abs(ev(base, x))/abs(x)**changed
                   for interval in intervals for x in interval)/4
        for orientation in (-1, 1):
            residual = base.copy()
            residual[changed] = orientation*size
            epsilon = first_root_sign(residual, intervals)
            check_pair(residual, s, kind, epsilon)
            cases += 1

for roots in ((F(1), F(2), F(3), F(-6)),
              (F(1), F(2), F(3), F(-15, 2), F(-17, 2), F(10))):
    s = len(roots)//2-1
    assert all(sum(x**k for x in roots) == 0 for k in range(1, 2*s, 2))
    for orientation in (-1, 1):
        oriented_roots = [orientation*x for x in roots]
        residual = [F(1)]
        for x in oriented_roots:
            residual = mul(residual, [-x, 1])
        epsilon = sign(min(oriented_roots, key=abs))
        check_pair(residual, s, 'B', epsilon)
        cases += 1
assert cases == 28
print('PASS: 28 exact real-rooted residual pairs; common leading coefficients, products, differences, and kappa signs.')

tournaments = 0
for parity in (0, 1):
    for q in range(5):
        # Coefficient order is B_low, A, B_high, all ascending internally.
        groups = [tuple(range(4, 4+q)), tuple(range(4)), tuple(range(4+q, 8))]
        order = sum((list(group) for group in groups), [])
        kappa = {row: rank for rank, row in enumerate(order)}
        group_id = {row: k for k, group in enumerate(groups) for row in group}
        epsilon = {}
        for i, j in combinations(range(8), 2):
            same_group = group_id[i] == group_id[j]
            coefficient_sign = sign(kappa[i]-kappa[j])
            e = (-1)**(parity if same_group else parity+1)*coefficient_sign
            epsilon[i, j] = e
            epsilon[j, i] = -e
        outgoing = {i: sum(epsilon[i, j] == 1 for j in range(8) if j != i)
                    for i in range(8)}
        assert sorted(outgoing.values()) == list(range(8))
        initial_order = sorted(range(8), key=lambda i: -outgoing[i])
        assert all(epsilon[i, j] == 1 for a, i in enumerate(initial_order)
                   for j in initial_order[a+1:])
        for b in range(4, 8):
            assert len({epsilon[a, b] for a in range(4)}) == 1
        allowed = []
        for v in product((-1, 1), repeat=8):
            if len(set(v)) == 1:
                continue
            if all((v[i]-v[j])//2 == epsilon[i, j]
                   for i, j in combinations(range(8), 2) if v[i] != v[j]):
                allowed.append(v)
        expected = [tuple(1 if i in initial_order[:k] else -1 for i in range(8))
                    for k in range(1, 8)]
        assert len(allowed) == 7 and set(allowed) == set(expected)
        tournaments += 1
assert tournaments == 10
print('PASS: all 10 parity/group-size cases; transitive first tournament and exactly 7 nonconstant directed cuts.')

# Formula for remaining common-positive roots is combinatorial and applies
# to arbitrary signed rows; exhaust all pairs of length-six sign strings.
strings = list(product((-1, 1), repeat=6))
for left in strings:
    for right in strings:
        p = sum(x > 0 for x in left)+sum(x > 0 for x in right)
        distance = sum(x != y for x, y in zip(left, right))
        for prefix in range(7):
            u = sum(x > 0 for x in left[:prefix])+sum(x > 0 for x in right[:prefix])
            k = sum(x != y for x, y in zip(left[:prefix], right[:prefix]))
            twice_remaining = p-distance-u+k
            actual = sum(x == y == 1 for x, y in zip(left[prefix:], right[prefix:]))
            assert twice_remaining == 2*actual
print('PASS: all 28,672 length-six signed-pair prefixes; exact remaining-common-root parity formula.')
