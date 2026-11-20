"""X.5(2b): the m=43 squarefree relaxation of IX.4, rebuilt from its specification."""
import sys
from itertools import combinations
from math import prod, gcd, isqrt
from collections import Counter

def primes_in(a, b):
    s = bytearray([1])*(b+1); s[0] = s[1] = 0
    for i in range(2, isqrt(b)+1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [p for p in range(a, b+1) if s[p]]

def H(a, b, c):
    return (a+b+c)*(-a+b+c)*(a-b+c)*(a+b-c)

def choose_semiprimes(m, P, ratio_den):
    pr = [p for p in primes_in(P, 2*P) if p % 4 == 1]
    T = (3*P//2)**2
    used, Q = set(), []
    prset = set(pr)
    for p in pr:
        if p in used: continue
        r = -(-T // p)
        while r <= 2*P and (p*r - T)*ratio_den <= T:
            if r in prset and r != p and r not in used:
                used |= {p, r}; Q.append((p, r)); break
            r += 1
        if len(Q) == m: break
    assert len(Q) == m
    return Q

def run(m, P, ratio_den):
    fac = choose_semiprimes(m, P, ratio_den)
    Qv = [p*r for p, r in fac]
    q = prod(Qv)
    QR = {(x*x) % m for x in range(1, m)}
    S = [[(1 if (v-i) % m in QR else -1) if v != i else -1 for v in range(m)] for i in range(m)]
    # Gram checks
    for i in range(m):
        for j in range(m):
            ip = sum(S[i][v]*S[j][v] for v in range(m))
            assert ip == (m if i == j else -1)
    k = P // (4*m)
    w = [k*i + (1 if i == m-1 else 0) for i in range(m)]
    d, h = {}, {}
    sep = {p: [] for pr in fac for p in pr}
    for i, j in combinations(range(m), 2):
        hij = prod(Qv[v] for v in range(m) if S[i][v] == S[j][v])
        dij = hij*(w[j]-w[i])
        assert dij*dij <= q and gcd(dij, q) == hij
        d[i, j] = dij; h[i, j] = hij
        for v in range(m):
            if S[i][v] != S[j][v]:
                for p in fac[v]:
                    sep[p].append(dij*dij % p)
    c = Counter()
    maxH = 0
    for i, j, l in combinations(range(m), 3):
        const = [v for v in range(m) if S[i][v] == S[j][v] == S[l][v]]
        assert len(const) == (m-3)//4
        g = prod(Qv[v] for v in const)
        assert g == gcd(gcd(h[i, j], h[j, l]), h[i, l])
        x = H(d[i, j], d[j, l], d[i, l])
        maxH = max(maxH, abs(x))
        c["triples"] += 1
        c["H nonzero"] += x != 0
        c["fail q^2 g^2 | H"] += x % (q*g)**2 != 0
        c["signed identity fails"] += (d[i, j] + d[j, l] != d[i, l])
    npairs = {len(v) for v in sep.values()}
    ndist = [len(set(v)) for v in sep.values()]
    import math
    print(f"m={m} P={P} primes={2*m} log q={math.log(q):.1f} max|d|/sqrt q={max(abs(x) for x in d.values())/math.isqrt(q):.4f}")
    print("  ", dict(c), "max|H|/q^2 =", maxH/q**2)
    print("   separated pairs per prime", npairs, "distinct d^2 mod p: min", min(ndist), "max", max(ndist))

if __name__ == "__main__":
    for P in map(int, sys.argv[1:]):
        run(43, P, 8*43)
