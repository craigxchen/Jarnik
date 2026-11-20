"""X.5(3): the four-point table at q=17*23*29 and the five-mark table at C=6."""
from itertools import combinations, product, permutations
from math import gcd, prod, isqrt

def H(a, b, c):
    return (a+b+c)*(-a+b+c)*(a-b+c)*(a+b-c)

def sign_free_check(n, D, primes, Cbound=None):
    """D: dict (i,j)->positive length. Check: at each p the zero edges are exactly the internal edges of a
    bipartition (cut), the nonzero ones have one common square mod p; gcd(d,q)=h. Returns the cut bits."""
    q = prod(primes)
    cuts = {}
    for p in primes:
        found = None
        for bits in product((0, 1), repeat=n-1):
            b = (0,) + bits
            if all((D[i, j] % p == 0) == (b[i] == b[j]) for i, j in combinations(range(n), 2)):
                found = b; break
        assert found is not None, ("no cut", p)
        sq = {D[i, j]**2 % p for i, j in combinations(range(n), 2) if found[i] != found[j]}
        assert len(sq) == 1, (p, sq)
        cuts[p] = found
    for i, j in combinations(range(n), 2):
        h = prod(p for p in primes if cuts[p][i] == cuts[p][j])
        assert gcd(D[i, j], q) == h
        if Cbound: assert D[i, j]**2 <= Cbound**2 * q
    return cuts

def line_embeddable(n, D):
    # brute force: fix vertex 0 at 0, each other vertex at +-D[0,j]
    for s in product((1, -1), repeat=n-1):
        x = [0] + [s[j-1]*D[0, j] for j in range(1, n)]
        if all(abs(x[i]-x[j]) == D[i, j] for i, j in combinations(range(n), 2)):
            return x
    return None

# (a) four-point table
q = 17*23*29
D = {(0, 1): 17, (0, 2): 46, (0, 3): 29, (1, 2): 29, (1, 3): 46, (2, 3): 17}
cuts = sign_free_check(4, D, [17, 23, 29])
print("q =", q, "sqrt q =", isqrt(q), "max length", max(D.values()), "46^2 < q:", 46**2 < q)
print("cuts", cuts)
sig = [tuple(cuts[p][i] for p in (17, 23, 29)) for i in range(4)]
print("signatures", sig)
xors = {}
for i, j in combinations(range(4), 2):
    xors.setdefault(tuple(a ^ b for a, b in zip(sig[i], sig[j])), []).append((i, j))
print("pair xors", xors, "Sidon:", all(len(v) == 1 for v in xors.values()))
print("Heron values of triples", [H(D[i, j], D[j, k], D[i, k]) for i, j, k in combinations(range(4), 3)])
print("line embeddable:", line_embeddable(4, D))
for p in (17, 23, 29):
    print("  p", p, "separated square", {D[i, j]**2 % p for i, j in combinations(range(4), 2) if cuts[p][i] != cuts[p][j]})

# (b) five-mark example
a, b = 83, 59
u, v = a*(a+b), b*(a+2*b)
pos = [0, u, v, u+v, (a+b)*(a+2*b)]
pr = [59, 67, 71, 83]
q5 = prod(pr)
print("\nu, v =", u, v, "u+v =", u+v, "|u-v| =", abs(u-v), "positions", pos, "q =", q5, "isqrt(q) =", isqrt(q5))
print("u factors divisible by", [p for p in pr if u % p == 0], "v by", [p for p in pr if v % p == 0])
for p in pr:
    print("  p", p, "classes", sorted({y % p for y in pos}), "residues", [y % p for y in pos])
orig = {(i, j): abs(pos[j]-pos[i]) for i, j in combinations(range(5), 2)}
c0 = sign_free_check(5, orig, pr, Cbound=6)
print("original: genuine, L =", max(pos)-min(pos), "L/sqrt q =", (max(pos)-min(pos))/q5**0.5)
sig5 = [tuple(c0[p][i] for p in pr) for i in range(5)]
print("signatures", sig5, "distinct:", len(set(sig5)) == 5)
alt = dict(orig); alt[0, 3] = abs(u-v)
c1 = sign_free_check(5, alt, pr, Cbound=6)
print("altered table is sign-free, same cuts:", c1 == c0, "max length", max(alt.values()), "6 sqrt q =", 6*q5**0.5)
for i, j, k in combinations(range(5), 3):
    h = H(alt[i, j], alt[j, k], alt[i, k])
    g = prod(p for p in pr if c1[p][i] == c1[p][j] == c1[p][k])
    assert h % (q5*g)**2 == 0
    print("  triple", (i+1, j+1, k+1), "lengths", (alt[i, j], alt[j, k], alt[i, k]), "g", g, "H/q^2 =", h//q5**2, "nonzero" if h else "")
print("altered line embeddable:", line_embeddable(5, alt))
# minimal deletion: which 4-subsets are line embeddable
for S in combinations(range(5), 4):
    sub = {(x, y): alt[S[x], S[y]] for x, y in combinations(range(4), 2)}
    print("  4-subset", [s+1 for s in S], "embeddable" if line_embeddable(4, sub) else "not")
