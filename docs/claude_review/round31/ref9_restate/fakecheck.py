# Checks for the rational-direction fake of IX.6 (exact where possible).
import math, itertools, random
from sympy import Matrix, legendre_symbol, nextprime, primerange, isprime
import mpmath as mp

def profile(q, b=8, astar=0):
    chi = lambda a: 0 if a % q == 0 else legendre_symbol(a % q, q)
    rows = ['inf'] + list(range(q))
    def H(x, a):  # a in F_q label, x row
        if x == 'inf': return 1
        return -((1 if x == a else 0) + chi(a - x))
    cols = []
    for a in range(q):
        for _ in range(b): cols.append([H(x, a) for x in rows])
    pair = {x: (x if x != 'inf' else astar) for x in rows}
    for xi, x in enumerate(rows):
        c = [H(y, pair[x]) for y in rows]; c[xi] = -c[xi]; cols.append(c)
    S = Matrix(cols).T  # M x r
    return S

for q in (7, 11, 19, 23):
    S = profile(q); M, r = S.shape
    assert all(abs(e) == 1 for e in S)
    rank = S.rank()
    # twin differences
    cols = [tuple(S[:, j]) for j in range(r)]
    colset = set(cols)
    ok = True
    for i in range(M):
        found = any(tuple(c[k] if k != i else -c[k] for k in range(M)) in colset for c in cols)
        ok &= found
    print(f"q={q}: M={M}, r={r}, 9M-8={9*M-8}, rank={rank}, twins for every row: {ok}")
    # integrality: v = S^T c' integral, c' real => c = 2c' integral; S^T c/2 integral iff sum c even
    for _ in range(200):
        c = Matrix([random.randint(-3, 3) for _ in range(M)])
        v = S.T * c
        par = set(int(x) % 2 for x in v)
        assert par == {int(sum(c)) % 2}
    print("   (S^T c)_j = sum(c) mod 2 for all j: checked on 200 random c")

# union bound prod(1+2/(sqrt p -1)) - 1 for s = 9m-8 primes >= m^A
print("\nHaar union bound: prod_{s primes >= m^A}(1+2/(sqrt p-1)) - 1 (taking the s smallest primes >= m^A)")
for A in (4, 5, 8):
    for m in ((8, 12, 20, 44, 100, 1000) if A < 8 else (8, 12, 20, 44, 100)):
        s = 9*m - 8; P = m**A
        prod = mp.mpf(1); p = P; cnt = 0
        while cnt < s:
            p = nextprime(p)
            if p % 4 == 1:
                prod *= 1 + 2/(mp.sqrt(p) - 1); cnt += 1
        print(f"  A={A} m={m} s={s}: {mp.nstr(prod-1, 6)}   (2s/sqrt(P) = {mp.nstr(2*s/mp.sqrt(P), 6)})")

# projection interval length: z = R e^{i theta}, theta - beta in [m^{5/4} D, (m^{5/4}+1) D], R D^2 = 1, |b| = sqrt5
mp.mp.dps = 60
print("\nprojection interval length / m^{5/4}")
for m in (10, 100, 10**4, 10**6):
    for W in (200, 2000):
        D = mp.e ** (-mp.mpf(W)/4); R = mp.e ** (mp.mpf(W)/2)
        x0 = mp.mpf(m) ** 1.25 * D
        length = mp.sqrt(5) * R * (mp.cos(x0) - mp.cos(x0 + D))
        print(f"  m={m} W={W}: length={mp.nstr(length, 10)}, ratio={mp.nstr(length / mp.mpf(m)**1.25, 10)}")

# dist(n beta, (pi/4) Z) for n up to 10^6 (n = sum c / 2 integer): min of dist * n
mp.mp.dps = 40
beta = mp.atan(mp.mpf(1)/2)
best = []
mn = mp.inf
for n in range(1, 200001):
    x = n*beta / (mp.pi/4)
    d = abs(x - mp.nint(x)) * mp.pi/4
    if d < mn:
        mn = d; best.append((n, d, d*n))
print("\nrecord minima of dist(n*beta, (pi/4)Z), n <= 2e5: (n, dist, n*dist)")
for t in best[-8:]:
    print("  ", t[0], mp.nstr(t[1], 6), mp.nstr(t[2], 6))
