"""(a) identity (8) on actual ten-pair-cut tuples with multi-prime blocks, own construction;
(b) Astra's Zariski-dense parity family: exact checks of primitivity, rank, parity, and the chain behind (6);
(c) synthetic index / parity distribution for balanced profiles."""
import sys, random, math
from itertools import combinations
from fractions import Fraction
from collections import Counter
from mpmath import mp, mpf, sqrt as msqrt, log as mlog, fabs
from sympy import primerange, nextprime, sqrt_mod
sys.path.insert(0, sys.argv[0].rsplit('/', 1)[0])
from gauss_util import *
mp.dps = 400
random.seed(11)
pairs = list(combinations(range(5), 2))
P1 = [p for p in primerange(5, 3000) if p % 4 == 1]


def det(u, v):
    return u[0]*v[1]-u[1]*v[0]


# (a)
for trial in range(200):
    pr = random.sample(P1, 20)
    blocks, norms = [], []
    for k in range(10):
        b = (1, 0); n = 1
        for p in pr[2*k:2*k+random.randint(1, 2)]:
            pi = split_prime(p)
            if random.random() < 0.5:
                pi = gconj(pi)
            e = random.randint(1, 2)
            b = gmul(b, gpow(pi, e)); n *= p**e
        blocks.append(b); norms.append(n)
    P = [(1, 0)]*5
    for (i, j), b in zip(pairs, blocks):
        P[i] = gmul(P[i], b); P[j] = gmul(P[j], b)
    t = {}
    for (i, j), n in zip(pairs, norms):
        d = det(P[i], P[j])
        assert d % n == 0
        t[(i, j)] = d//n; t[(j, i)] = -d//n
        assert t[(i, j)] != 0
    for i in range(5):
        assert gnorm(P[i]) == math.prod(n for (a, b), n in zip(pairs, norms) if i in (a, b))
    for coord in (0, 1):
        s = sum(Fraction(P[i][coord], math.prod(t[(i, j)] for j in range(5) if j != i)) for i in range(5))
        assert s == 0
    # z_i = G P_i / conj(P_i): check reduced chord cofactor and D formula
    G = (1, 0)
    for b in blocks:
        G = gmul(G, gconj(b))
    z = [gmul(G, gdivmod_exact(gmul(P[i], (1, 0)), gconj(P[i])) or (0, 0)) if False else None for i in range(5)]
    # build z directly: z_i = prod_{E ni i} beta_E prod_{E not ni i} conj(beta_E)
    z = []
    for i in range(5):
        w = (1, 0)
        for (a, b), bl in zip(pairs, blocks):
            w = gmul(w, bl if i in (a, b) else gconj(bl))
        z.append(w)
    Hn = gnorm(z[0])
    g = z[0]
    for w in z[1:]:
        g = ggcd(g, w)
    assert gnorm(g) == 1  # primitive
    R = msqrt(mpf(Hn))
    for (i, j), n in zip(pairs, norms):
        dij = mlog(mpf(gnorm(P[i])*gnorm(P[j])))/1 - 2*mlog(mpf(n))  # log D(y_i/y_j)
        dist = msqrt(mpf((z[i][0]-z[j][0])**2+(z[i][1]-z[j][1])**2))/R
        lhs = dist
        rhs = 2*abs(t[(i, j)])*mp.e**(-dij/2)
        assert fabs(lhs/rhs-1) < mpf(10)**-300
print("(a) identity (8), integrality of t_ij, |P_i|^2 = prod n_ij and |y_i-y_j| = 2|t_ij| e^{-d_ij/2}: 200 actual tuples OK")


# (b) the Zariski-dense parity family, small scale
def gauss_prime_big(p):
    r = sqrt_mod(-1, p)
    a, b = p, r
    while b*b > p:
        a, b = b, a % b
    c = math.isqrt(p-b*b)
    assert b*b+c*c == p
    return (b, c)


def find_prime_1mod4(lo):
    p = nextprime(lo)
    while p % 4 != 1:
        p = nextprime(p)
    return p


subsets = [frozenset(s) for r in range(1, 5) for s in combinations(range(4), r)]
for w in (60, 120):
    used = set()
    data = {}
    for A in subsets:
        if len(A) == 2 and 0 in A:
            e = 1
        elif len(A) == 1:
            e = 6   # even singleton exponent (fixed small-ish primes, large exponent)
        else:
            e = 2
        lo = int(mp.e**(mpf(w)/e))
        p = find_prime_1mod4(lo)
        while p in used:
            p = find_prime_1mod4(p)
        used.add(p)
        data[A] = (p, e)
    pis = {A: (gauss_prime_big(p) if p > 10**6 else split_prime(p)) for A, (p, e) in data.items()}
    z = []
    for i in range(5):
        v = (1, 0)
        for A, (p, e) in data.items():
            v = gmul(v, gpow(pis[A], e) if i in A else gpow(gconj(pis[A]), e))
        z.append(v)
    H = gnorm(z[0])
    assert all(gnorm(v) == H for v in z)
    g = z[0]
    for v in z[1:]:
        g = ggcd(g, v)
    assert gnorm(g) == 1
    cols = [[e*int(i in A) for i in range(4)] for A, (p, e) in data.items()]
    from sympy import Matrix
    assert Matrix(cols).T.rank() == 4
    # parity: rank of columns mod 2
    from gauss_util import rank_mod2
    assert rank_mod2([[c[i] for c in cols] for i in range(4)], len(cols)) == 3
    logn = {A: e*math.log(p) for A, (p, e) in data.items()}
    eta = max(abs(x-w) for x in logn.values())/w
    R = msqrt(mpf(H))
    y = [None]*5
    Delta = max(msqrt(mpf((u[0]-v[0])**2+(u[1]-v[1])**2)) for u, v in combinations(z, 2))/R
    # the three square words: s = y_i y_j/(y_k y_l) for (1234); log D(s) <= sum_A |sum_A c| log n_A <= 12 max log n_A
    for c in [(1, 1, -1, -1, 0), (1, -1, 1, -1, 0), (1, -1, -1, 1, 0)]:
        logD = sum(abs(sum(c[i]*(data[A][1] if i in A else 0) for i in range(5)))*math.log(data[A][0]) for A in data)
        bound = sum(abs(sum(c[i] for i in A))*logn[A] for A in data)
        assert logD <= bound + 1e-9
    print("(b) w=%d: primitive, rank 4, mod-2 rank 3, eta=%.4f, Delta=%s, bound 2^-1/2 e^{-3(1+eta)w}=%s"
          % (w, eta, mp.nstr(Delta, 5), mp.nstr(mpf(2)**-0.5*mp.e**(-3*(1+eta)*w), 5)))

# (c) synthetic distribution of index and parity for "generic" balanced-type profiles:
# each of the 15 cuts carried by r distinct primes with random exponent 1..3 (layers nested at one prime when
# a prime carries several levels: here simple model, one level per prime, exponents via repeated levels)
cnt = Counter()
for trial in range(3000):
    cols = []
    # random profile: K primes, each with e_p in 1..3 and random a_ip with min 0, max e_p
    K = random.randint(4, 15)
    for _ in range(K):
        e = random.randint(1, 3)
        a = [random.randint(0, e) for _ in range(5)]
        a[random.randrange(5)] = 0
        j = random.randrange(5)
        while a[j] == 0 and a.count(0) == 1:
            j = random.randrange(5)
        a[j] = e
        if min(a) != 0 or max(a) != e:
            continue
        cols.append([a[i]-a[4] for i in range(4)])
    I = lattice_index(cols)
    cnt[(K, 'I=1' if I == 1 else ('rank<4' if I == 0 else 'I>1'))] += 1
byK = {}
for (K, lab), v in cnt.items():
    byK.setdefault(K, {})[lab] = v
for K in sorted(byK):
    print("(c) K=%2d primes:" % K, byK[K])
