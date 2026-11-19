"""Exact primitive near-square census over all squarefree parents in a prime pool.

The upper window d <= sqrt(N) (log N)^B is omitted, so every D_B candidate
is covered for every B. No floating-point arithmetic is used.
"""
from collections import Counter, defaultdict
from itertools import combinations
from math import gcd, isqrt, prod


def squarefree_divisors(primes, start, bound, value=1, mask=0):
    yield value, mask
    for index in range(start, len(primes)):
        next_value = value * primes[index]
        if next_value > bound:
            break
        yield from squarefree_divisors(primes, index + 1, bound,
                                       next_value, mask | (1 << index))


def phase_compatible(a, b, c):
    # A signed product of the three primitive Gaussian factors must be real.
    # Conjugation replaces a negative sign on the corresponding angle.
    x, y, z = map(isqrt, (a, b, c))
    u, v, w = isqrt(a-x*x), isqrt(b-y*y), isqrt(c-z*z)
    for v_sign in (-v, v):
        ab_re, ab_im = x*y - u*v_sign, u*y + x*v_sign
        for w_sign in (-w, w):
            if ab_re*w_sign + ab_im*z == 0:
                return True
    return False


def check(primes):
    all_norms = prod(primes)
    parents = defaultdict(list)
    for d, mask in squarefree_divisors(primes, 0, all_norms):
        s = isqrt(d)
        offset = d - s*s
        t = isqrt(offset)
        if not t or t*t != offset or gcd(s, t) != 1:
            continue
        bound = d // (256 * offset**2)
        if not bound:
            continue
        unused = [p for j, p in enumerate(primes) if not mask & (1 << j)]
        for multiplier, _ in squarefree_divisors(unused, 0, bound):
            parents[d * multiplier].append(d)
    counts = Counter(map(len, parents.values()))
    best = max(counts)
    triangles = []
    for parent, divisors in parents.items():
        for a, b, c in combinations(divisors, 3):
            product = a*b*c
            if isqrt(product)**2 == product:
                triangles.append((parent, sorted((a, b, c))))
    print('Prime pool:', primes)
    print('Parents:', 2**len(primes), 'largest parent:', all_norms)
    print('Counts for nonempty parents:', sorted(counts.items()))
    parent = min(n for n, ds in parents.items() if len(ds) == best)
    print('Smallest maximizing parent:', parent, 'divisors:', sorted(parents[parent]))
    print('Triangle count:', len(triangles), 'first witness:', min(triangles, default=None))
    distinct = {tuple(d) for _, d in triangles}
    compatible = {d for d in distinct if phase_compatible(*d)}
    print('Different triples:', len(distinct), 'phase-compatible:', len(compatible))
    print('Parent/triple incidences passing phases:',
          sum(tuple(d) in compatible for _, d in triangles))
    print('First phase-compatible triple:', min(compatible, default=None))
    cliques = []
    for parent, divisors in parents.items():
        if len(divisors) < 6:
            continue
        steps = {sum(1 << j for j, p in enumerate(primes) if d % p == 0)
                 for d in divisors}
        for a, b, c in combinations(steps, 3):
            labels = {a, b, c, a ^ b, a ^ c, b ^ c}
            if len(labels) == 6 and labels <= steps:
                cliques.append((parent, (a, b, c)))
    print('Four-vertex cliques with six distinct labels, anchored at 0:', len(cliques))


if __name__ == '__main__':
    pool = [n for n in range(5, 301, 4)
            if all(n % q for q in range(2, isqrt(n) + 1))]
    check(pool[:22])
