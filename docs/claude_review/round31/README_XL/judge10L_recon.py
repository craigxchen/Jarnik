"""[J] Checks of X.4 (own code, exact integers).
For q squarefree, delta, C: an admissible magnitude d (1 <= d <= L = floor(C sqrt q)) has a pattern
x in F_2^s (x_p = 0 if p | d, 1 if d^2 = delta_p^2 mod p; otherwise d is inadmissible).
A sign-free table on signatures V exists iff every difference u+v (u != v in V) is a realised pattern,
so B_free = clique number of Cay(F_2^s, X).  B_int = max number of integers in a window of length L
occupying residues {0, eps_p delta_p} mod p (eps free; translation = moving the window).
Checks: (1) statement 1, B_free <= 2 max(4, B_int) (when some p > C^2) and <= max(4,B_int) at C = 1;
(2) uniqueness of magnitudes on patterns with h > 1 when all p > C^2 (referee recon);
(3) on maximum cliques with random magnitudes on multi-valued patterns: line-embeddability,
    and embeddability after deleting one endpoint of every complementary pair;
(4) in every genuine configuration found, #complementary pairs <= 1 + floor(C^2) and
    = number of distinct values of (y_i - y_j)^2 (all primes balanced or not);
(5) the five-mark example and the 8-clique at q = 59*67*71*83, C = 6."""
import random
from itertools import combinations, product
from math import isqrt, prod
import networkx as nx

random.seed(20261120)


def patterns(primes, delta, L):
    X = {}
    for d in range(1, L + 1):
        x = 0
        for k, p in enumerate(primes):
            r = d % p
            if r == 0:
                continue
            if (r * r - delta[k] * delta[k]) % p == 0:
                x |= 1 << k
            else:
                x = -1
                break
        if x > 0:
            X.setdefault(x, []).append(d)
    return X


def crt(res, primes):
    q = prod(primes)
    y = 0
    for r, p in zip(res, primes):
        M = q // p
        y += r * M * pow(M, -1, p)
    return y % q


def b_int(primes, delta, L):
    q = prod(primes)
    s = len(primes)
    best, bestconf = 0, None
    for eps in product((1, -1), repeat=s - 1):
        eps = (1,) + eps
        ys = sorted(crt([b * e * dl for b, e, dl in zip(bits, eps, delta)], primes)
                    for bits in product((0, 1), repeat=s))
        ext = ys + [y + q for y in ys]
        j = 0
        for i in range(len(ys)):
            j = max(j, i)
            while j + 1 < len(ext) and ext[j + 1] - ext[i] <= L:
                j += 1
            if j - i + 1 > best:
                best, bestconf = j - i + 1, ext[i:j + 1]
    return best, bestconf


def clique_number(s, X):
    G = nx.Graph()
    G.add_nodes_from(range(1 << s))
    for u in range(1 << s):
        for x in X:
            v = u ^ x
            if v > u:
                G.add_edge(u, v)
    cl = max(nx.find_cliques(G), key=len)
    return len(cl), cl, G


def embeds(V, d):
    if len(V) <= 2:
        return True
    u, v = max(combinations(V, 2), key=lambda e: d[frozenset(e)])
    x = {w: (0 if w == u else d[frozenset((u, w))]) for w in V}
    return all(abs(x[i] - x[j]) == d[frozenset((i, j))] for i, j in combinations(V, 2))


def table_tests(s, X, clique, trials, full):
    """random magnitudes on a clique; return (#embeddable, #embeddable after complementary deletion)"""
    e1 = e2 = 0
    for _ in range(trials):
        d = {frozenset((i, j)): random.choice(X[i ^ j]) for i, j in combinations(clique, 2)}
        e1 += embeds(clique, d)
        V = list(clique)
        for i in list(V):
            if i in V and (full ^ i) in V:
                V.remove(full ^ i)
        e2 += embeds(V, d)
    return e1, e2


def genuine_checks(conf, primes, C):
    q = prod(primes)
    sig = lambda y: tuple(0 if y % p == conf[0] % p else 1 for p in primes)
    S = [sig(y) for y in conf]
    full = tuple(1 for _ in primes)
    comp = [(a, b) for a, b in combinations(range(len(conf)), 2)
            if tuple(x ^ y for x, y in zip(S[a], S[b])) == full]
    sq = {(conf[a] - conf[b]) ** 2 for a, b in comp}
    sums = {conf[a] + conf[b] for a, b in comp}
    assert len(comp) <= 1 + int(C * C)
    assert len(sums) <= 1 and len(sq) == len(comp)
    return len(comp)


def small_primes(lo, hi):
    return [p for p in range(lo, hi) if all(p % r for r in range(2, isqrt(p) + 1))]


out = []
# (5) the five-mark example and its modulus
primes = [59, 67, 71, 83]
q = prod(primes)
delta = [45, 61, 2, 73]
marks = [0, 11786, 11859, 23645, 28542]
for p, dl in zip(primes, delta):
    assert {y % p for y in marks} == {0, dl}
d = {frozenset((i, j)): abs(marks[i] - marks[j]) for i, j in combinations(range(5), 2)}
d[frozenset((0, 3))] = 73
H = lambda a, b, c: (a + b + c) * (-a + b + c) * (a - b + c) * (a + b - c)
h145 = H(73, 28542, 28542 - 23645)
out.append("five-mark: H(1,4,5)/q^2 = %s, embeds=%s, embeds without mark1=%s, without mark4=%s" % (
    h145 // q ** 2 if h145 % q ** 2 == 0 else "nonint", embeds(range(5), d),
    embeds([1, 2, 3, 4], d), embeds([0, 1, 2, 4], d)))
for C2 in (36, 49, 58):
    L = isqrt(C2 * q)
    X = patterns(primes, delta, L)
    multi = {x: v for x, v in X.items() if len(v) > 1}
    bf, cl, _ = clique_number(4, X)
    bi, conf = b_int(primes, delta, L)
    e1, e2 = table_tests(4, X, cl, 400, 15)
    out.append("q=%d C^2=%d: patterns=%d, multi-magnitude patterns=%s, B_free=%d clique=%s, B_int=%d; "
               "random tables embeddable %d/400, after complementary deletion %d/400" % (
                   q, C2, len(X), {format(x, '04b'): v for x, v in multi.items()}, bf,
                   [format(c, '04b') for c in cl], bi, e1, e2))

# (1)-(4) random moduli, primes > C^2, random delta
viol1 = viol2 = 0
cases = 0
worst = []
emb_fail = 0
emb_tot = 0
for C in (1, 2, 3):
    for trial in range(60):
        s = random.choice([4, 5, 6])
        lo = max(C * C + 1, 5)
        pool = small_primes(lo, 4 * lo + 40)
        primes = sorted(random.sample(pool, s))
        q = prod(primes)
        L = isqrt(C * C * q)
        bestd = None
        for _ in range(30):
            delta = [random.randrange(1, p) for p in primes]
            X = patterns(primes, delta, L)
            if bestd is None or len(X) > bestd[0]:
                bestd = (len(X), delta, X)
        _, delta, X = bestd
        cases += 1
        for x, v in X.items():
            h = prod(p for k, p in enumerate(primes) if not (x >> k) & 1)
            if h > 1 and len(v) > 1:
                viol2 += 1
        bf, cl, _ = clique_number(s, X)
        bi, conf = b_int(primes, delta, L)
        bound = max(4, bi) if C == 1 else 2 * max(4, bi)
        if bf > bound:
            viol1 += 1
        worst.append((bf - bi, C, primes, bf, bi))
        if conf and len(conf) >= 2:
            genuine_checks(conf, primes, C)
        if bf >= 5:
            e1, e2 = table_tests(s, X, cl, 50, (1 << s) - 1)
            emb_tot += 50
            emb_fail += 50 - e2
worst.sort(reverse=True)
out.append("random moduli: %d cases; statement-1 violations %d; uniqueness violations (h>1, two magnitudes) %d"
           % (cases, viol1, viol2))
out.append("largest B_free - B_int: %s" % worst[:5])
out.append("cliques with B_free>=5: tables not embeddable after complementary deletion: %d/%d" % (emb_fail, emb_tot))
print("\n".join(out))
