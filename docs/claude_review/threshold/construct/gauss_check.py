"""Exact check of the tropical forced-divisor theorem (Theorem C1 of construct.md) on actual Gaussian
integers of a common norm N, including an inert prime (3) and the ramified prime (1+i), and including
auxiliary forms with Gaussian-integer coefficients carrying prime content.

For F = sum_k c_k (u1 v1)^{(D-|k|_1)/2} u^{k+} v^{k-}, points z_j with |z_j|^2 = N, and a split prime
p = pi*conj(pi) with v_p(N) = e, a_j = v_pi(z_j), b_j = a_j - e/2:
    v_pi(F(z))      >= min_k ( v_pi(c_k)      + e D/2 + <k,b> )
    v_conjpi(F(z))  >= min_k ( v_conjpi(c_k)  + e D/2 - <k,b> )
For an inert prime q (q^{2f} || N) and for 1+i, every monomial has the same valuation, so the bound is
min_k v(c_k) + (valuation of N^{D/2}).  All arithmetic is exact (Python integers)."""
import random, itertools, math

def gmul(a, b): return (a[0]*b[0] - a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def gadd(a, b): return (a[0] + b[0], a[1] + b[1])
def gconj(a): return (a[0], -a[1])
def gpow(a, n):
    r = (1, 0)
    for _ in range(n): r = gmul(r, a)
    return r
def gdivexact(a, b):
    """a / b if exact in Z[i], else None"""
    n = b[0]*b[0] + b[1]*b[1]
    t = gmul(a, gconj(b))
    if t[0] % n or t[1] % n: return None
    return (t[0] // n, t[1] // n)
def gval(a, pi):
    if a == (0, 0): return math.inf
    v = 0
    while True:
        q = gdivexact(a, pi)
        if q is None: return v
        a = q; v += 1

def two_squares(p):
    for x in range(1, int(math.isqrt(p)) + 1):
        y2 = p - x*x; y = math.isqrt(y2)
        if y*y == y2 and x > y: return (x, y)
    raise ValueError

random.seed(77)
N = 2 * 9 * 25 * 13 * 17 * 29
split = {5: 2, 13: 1, 17: 1, 29: 1}
pis = {p: two_squares(p) for p in split}
pts = [(x, y) for x in range(-math.isqrt(N), math.isqrt(N) + 1) for y in [math.isqrt(N - x*x)]
       if y*y == N - x*x for y in ([y, -y] if y else [0])]
pts = sorted(set(pts))
assert all(x*x + y*y == N for x, y in pts)
print("N =", N, " lattice points:", len(pts))

def F_value(c, D, z):
    """c: dict k -> Gaussian integer coefficient; z: list of Gaussian integers."""
    M = len(z)
    Nn = gmul(z[0], gconj(z[0]))
    tot = (0, 0)
    for k, ck in c.items():
        t1 = sum(abs(x) for x in k)
        assert (D - t1) % 2 == 0 and t1 <= D
        term = gmul(ck, gpow(Nn, (D - t1) // 2))
        for j in range(M):
            if k[j] > 0: term = gmul(term, gpow(z[j], k[j]))
            elif k[j] < 0: term = gmul(term, gpow(gconj(z[j]), -k[j]))
        tot = gadd(tot, term)
    return tot

def check(c, D, z):
    val = F_value(c, D, z)
    if val == (0, 0): return "zero"
    for p, e in split.items():
        pi = pis[p]; pib = gconj(pi)
        a = [gval(zj, pi) for zj in z]
        # 2*(e D/2 + <k,b>) = e D + sum_j k_j (2 a_j - e)
        lo_pi = min(2 * gval(ck, pi) + e * D + sum(kj * (2 * aj - e) for kj, aj in zip(k, a)) for k, ck in c.items())
        lo_pib = min(2 * gval(ck, pib) + e * D - sum(kj * (2 * aj - e) for kj, aj in zip(k, a)) for k, ck in c.items())
        assert 2 * gval(val, pi) >= lo_pi, ("pi", p)
        EQ[0] += (2 * gval(val, pi) == lo_pi); EQ[1] += 1
        assert 2 * gval(val, pib) >= lo_pib, ("pibar", p)
    # inert 3 (3^2 || N, so v_3(z_j) = 1) and ramified 1+i (v_{1+i}(z_j) = 1)
    for g, vz in (((3, 0), 1), ((1, 1), 1)):
        lo = min(gval(ck, g) for ck in c.values()) + vz * D
        assert gval(val, g) >= lo, ("inert/ramified", g)
    return "ok"

def rand_gauss_content():
    t = (random.choice([1, -1, 2, 3]), random.choice([0, 1, -1]))
    for p in split:
        if random.random() < 0.3: t = gmul(t, pis[p])
        if random.random() < 0.3: t = gmul(t, gconj(pis[p]))
    if t == (0, 0): t = (1, 0)
    return t

stats = {"ok": 0, "zero": 0}
EQ = [0, 0]
for trial in range(400):
    M = random.choice([2, 3, 4])
    z = random.sample(pts, M)
    kind = random.choice(["alt", "rand", "content"])
    if kind == "alt":
        r = (M - 1) // 2
        v = tuple(range(-r, M - r))
        c = {}
        for perm in itertools.permutations(range(M)):
            sgn = 1
            for i in range(M):
                for j in range(i + 1, M):
                    if perm[i] > perm[j]: sgn = -sgn
            c[tuple(v[perm[i]] for i in range(M))] = (sgn, 0)
        D = sum(abs(x) for x in v)
    else:
        s = random.randint(-1, 1); D = random.randint(2, 4)
        if (D - s) % 2: D += 1
        c = {}
        for _ in range(random.randint(1, 5)):
            k = [random.randint(-2, 2) for _ in range(M - 1)]
            k.append(s - sum(k))
            if sum(abs(x) for x in k) <= D:
                c[tuple(k)] = rand_gauss_content() if kind == "content" else (random.randint(-3, 3) or 1, 0)
        if not c: continue
    stats[check(c, D, z)] += 1
print("tropical bound verified exactly:", stats, " pi-bound attained with equality in", EQ[0], "of", EQ[1], "prime checks")
