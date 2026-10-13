# Exact check of README Prop. H: greedy parity-respecting labels with no norm-<=K monomial colliding
# (any unit, level 1) at a modulus q > Q_K = 8K(2r+1)^K + 1. Labels live in Z/m, m = q - chi(q),
# chi = chi_{-4}; units = (m/4)Z/m; parity nu_j(q) = [(p_j/q) = -1].
import itertools, random, sys
from sympy import isprime, legendre_symbol, primerange

def monomials(j, K):
    # all v in Z^(j+1) with v_j != 0, ||v||_1 <= K (support in 0..j)
    out = []
    def rec(i, rem, cur):
        if i == j:
            for vj in range(-rem, rem + 1):
                if vj != 0:
                    out.append(cur + [vj])
            return
        for vi in range(-rem, rem + 1):
            rec(i + 1, rem - abs(vi), cur + [vi])
    rec(0, K, [])
    return out

def greedy(q, ps, K, rng):
    m = q - (1 if q % 4 == 1 else -1)
    assert m % 4 == 0
    U = {0, m // 4, m // 2, 3 * m // 4}
    nu = [1 if legendre_symbol(p % q, q) == -1 else 0 for p in ps]
    ell = []
    for j in range(len(ps)):
        mons = monomials(j, K)
        cand = [x for x in range(m) if x % 2 == nu[j]]
        rng.shuffle(cand)
        ok = None
        for x in cand:
            lab = ell + [x]
            if all((sum(v[i] * lab[i] for i in range(j + 1)) % m) not in U for v in mons):
                ok = x; break
        if ok is None:
            return None
        ell.append(ok)
    return ell, m, nu

def verify(ell, m, nu, K):
    U = {0, m // 4, m // 2, 3 * m // 4}
    r = len(ell)
    assert all(ell[j] % 2 == nu[j] for j in range(r))
    for v in itertools.product(range(-K, K + 1), repeat=r):
        n1 = sum(abs(t) for t in v)
        if 1 <= n1 <= K and (sum(a * b for a, b in zip(v, ell)) % m) in U:
            return False
    return True

rng = random.Random(1)
ps = [p for p in primerange(5, 400) if p % 4 == 1]
for r, K in [(4, 2), (4, 3), (6, 2)]:
    P = ps[:r]
    QK = 8 * K * (2 * r + 1) ** K + 1
    # (a) at the first 3 primes above Q_K greedy must succeed (proof) -- verify exhaustively
    qs = [q for q in primerange(QK + 1, QK + 2000) if q % 2 and q not in P][:3]
    res = []
    for q in qs:
        g = greedy(q, P, K, rng)
        res.append((q, g is not None and verify(*g, K)))
    # (b) empirical: smallest prime q > 5 (not dividing N) from which greedy succeeds at every later q up to Q_K
    first_fail = None
    for q in primerange(7, QK + 1):
        if q in P: continue
        if greedy(q, P, K, rng) is None:
            first_fail = q
    print(f"r={r} K={K} primes={P} Q_K={QK}: proof-range checks {res}; largest q<=Q_K where greedy failed: {first_fail}")
