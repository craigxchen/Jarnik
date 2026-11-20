"""Exact maximal two-class counts B_max(q,C) for squarefree q (own code).

B_max(q,C) = max size of a set of integers in an interval of length L=floor(C*sqrt q)
occupying at most two residue classes modulo every prime p | q.
Fact used: a set occupies <=2 classes mod p iff every triple of it has a pair congruent mod p;
for squarefree q this is: q | (b-a)(c-b)(c-a) for every triple a<b<c.
Translate so that the minimum is 0.  Triples (0,a,b), 0<a<b<=L, with q | a*b*(b-a):
enumerate, for every assignment of primes to the three factors a, b, b-a, the points of the
index-q lattice {a = 0 mod q1, b = 0 mod q2, b = a mod q3} in the triangle (exact integers).
Then search cliques in the graph on candidates a with edge (a,b) iff (0,a,b) admissible,
requiring all inner triples admissible.  Also verifies the signed triangle identity (IX.0a)
and the heavy-block claim on every admissible triple found.
"""
import sys, random, itertools
from math import isqrt, prod, log
from sympy import primerange


def crt_pair(r1, m1, r2, m2):
    # x = r1 mod m1, x = r2 mod m2, coprime
    t = ((r2 - r1) * pow(m1, -1, m2)) % m2
    return (r1 + m1 * t) % (m1 * m2)


def lattice_points(q1, q2, q3, L):
    """all (a,b) with 0<a<b<=L, q1|a, q2|b, q3|(b-a)."""
    m = q2 * q3
    c = crt_pair(0, q2, q1 % q3, q3) if q3 > 1 else 0
    if q2 == 1 and q3 == 1:
        c = 0
    e1 = (q1, c)
    e2 = (0, m)
    # Gauss reduction
    def dot(u, v):
        return u[0] * v[0] + u[1] * v[1]
    if dot(e1, e1) > dot(e2, e2):
        e1, e2 = e2, e1
    while True:
        n1 = dot(e1, e1)
        mu = round(dot(e1, e2) / n1) if n1 else 0
        # exact rounding
        num = dot(e1, e2)
        mu = (2 * num + n1) // (2 * n1)
        e2 = (e2[0] - mu * e1[0], e2[1] - mu * e1[1])
        if dot(e2, e2) < dot(e1, e1):
            e1, e2 = e2, e1
        else:
            break
    det = abs(e1[0] * e2[1] - e1[1] * e2[0])
    assert det == q1 * q2 * q3
    R2 = 2 * L * L  # points satisfy a^2+b^2 <= 2L^2
    n1 = dot(e1, e1)
    # x = i e1 + j e2 ; component of x orthogonal to e1 has length |j| det/|e1|
    # so j^2 det^2 / n1 <= R2
    jmax = isqrt(R2 * n1 // (det * det)) + 1
    out = []
    for j in range(-jmax, jmax + 1):
        lo, hi = -10**40, 10**40
        # constraints: a >= 1, b - a >= 1, L - b >= 0, each alpha*i + beta >= gamma
        for al, be, ga in ((e1[0], j * e2[0], 1),
                           (e1[1] - e1[0], j * (e2[1] - e2[0]), 1),
                           (-e1[1], -j * e2[1], -L)):
            # al*i >= ga - be
            rhs = ga - be
            if al > 0:
                lo = max(lo, -((-rhs) // al))
            elif al < 0:
                hi = min(hi, rhs // al if False else (-rhs) // (-al))
            elif rhs > 0:
                lo, hi = 1, 0
        for i in range(lo, hi + 1):
            a = i * e1[0] + j * e2[0]
            b = i * e1[1] + j * e2[1]
            assert 0 < a < b <= L
            out.append((a, b))
    return out


def admissible_triples(primes, L):
    s = len(primes)
    T = set()
    for assign in itertools.product(range(3), repeat=s):
        q1 = prod(p for p, x in zip(primes, assign) if x == 0)
        q2 = prod(p for p, x in zip(primes, assign) if x == 1)
        q3 = prod(p for p, x in zip(primes, assign) if x == 2)
        for ab in lattice_points(q1, q2, q3, L):
            T.add(ab)
    return T


def check_triangle(primes, ys, C, q):
    """(IX.0a) and heavy block on the triple ys=(y1,y2,y3)."""
    y1, y2, y3 = ys
    A = {'12': 1, '23': 1, '13': 1}
    G = 1
    for p in primes:
        r = [y % p for y in ys]
        if r[0] == r[1] == r[2]:
            G *= p
        elif r[0] == r[1]:
            A['12'] *= p
        elif r[1] == r[2]:
            A['23'] *= p
        elif r[0] == r[2]:
            A['13'] *= p
        else:
            raise AssertionError('three classes')
    assert G * A['12'] * A['23'] * A['13'] == q
    t12 = (y2 - y1) // (G * A['12']); assert t12 * G * A['12'] == y2 - y1
    t23 = (y3 - y2) // (G * A['23']); assert t23 * G * A['23'] == y3 - y2
    t13 = (y3 - y1) // (G * A['13']); assert t13 * G * A['13'] == y3 - y1
    assert t12 * A['12'] + t23 * A['23'] == t13 * A['13']
    W = log(q)
    # heavy block among A13, A23 (primes separating 1 and 2)
    heavy = max(log(A['13']), log(A['23']))
    assert heavy >= W / 4 - log(C) / 2 - 1e-9
    return abs(t12 * t23 * t13)


def max_clique(primes, L, q):
    T = admissible_triples(primes, L)
    adj = {}
    for a, b in T:
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set()).add(a)
    best = [0, 1] if L >= 1 else [0]
    if T:
        best = [0, *min(T)]

    def ok(x, y, z):
        return ((y - x) * (z - y) * (z - x)) % q == 0

    sys.setrecursionlimit(10000)

    def extend(cur, cands):
        nonlocal best
        if len(cur) + len(cands) <= len(best):
            return
        if len(cur) > len(best):
            best = list(cur)
        for idx, c in enumerate(cands):
            if all(ok(x, y, c) for x, y in itertools.combinations(cur[1:], 2)) and all(
                    ok(0, x, c) for x in cur[1:]):
                nc = [d for d in cands[idx + 1:] if d in adj[c]]
                extend(cur + [c], nc)
    for a in sorted(adj):
        extend([0, a], sorted(d for d in adj[a] if d > a))
    return best, T


if __name__ == '__main__':
    seed = int(sys.argv[1]); s = int(sys.argv[2]); P = int(sys.argv[3]); C = float(sys.argv[4]); trials = int(sys.argv[5])
    rng = random.Random(seed)
    pool = list(primerange(P, 2 * P))
    K = len(pool)
    dist = {}
    worst = None
    for tr in range(trials):
        primes = sorted(rng.sample(pool, s))
        q = prod(primes)
        L = int(C * isqrt(q)) if C >= 1 else isqrt(int(C * C * q))
        L = isqrt(int(C * C * 10**6) * q // 10**6)  # floor(C sqrt q) approx exact for rational C
        best, T = max_clique(primes, L, q)
        for ab in list(T)[:2000]:
            check_triangle(primes, (0, ab[0], ab[1]), C, q)
        B = len(best)
        dist[B] = dist.get(B, 0) + 1
        if worst is None or B > worst[0]:
            worst = (B, primes, best)
    W = s * log(1.5 * P)
    print(f's={s} P={P} K={K} C={C} trials={trials} dist={dict(sorted(dist.items()))} '
          f'worst={worst} larger-sieve scale logq/loglogq={W/log(W):.2f}')
