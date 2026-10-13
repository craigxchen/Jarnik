"""Random-point-label control for check_pairs_fast.py (same parities, labels uniform in parity class).
Fast version of check_pairs.py via discrete logarithms in T_Q (exact): labels ell_j = log_g t_j mod m_Q,
point labels kappa = A ell, pair (x,y) collides at q^a with unit 1 iff kappa_x = kappa_y mod m_(q^a), with some
unit iff kappa_x - kappa_y = 0 mod m_(q^a)/4. Cross-checked against check_pairs.py on q0 = 43, 59."""
import random, sys, math
import numpy as np
from gauss import *
from check_T3 import profile

def factor(n):
    f, d = [], 2
    while d*d <= n:
        if n % d == 0:
            f.append(d)
            while n % d == 0: n //= d
        d += 1
    if n > 1: f.append(n)
    return f

def T_generator_and_log(q, a, rng):
    Q = q**a; m = tq_order(q, a); pf = factor(m)
    while True:
        x = (rng.randrange(Q), rng.randrange(Q))
        if norm(x, Q) % q == 0: continue
        g = mul(x, inv(conj(x, Q), Q), Q)
        if all(pw(g, m//p, Q) != (1, 0) for p in pf): break
    table, e = {}, (1, 0)
    for k in range(m):
        table[e] = k; e = mul(e, g, Q)
    return g, table, m

def run(q0, b=6, seed=5):
    rng = random.Random(seed); rng2 = random.Random(seed + 1000); rng3 = random.Random(seed + 2000)
    A, idx, twins = profile(q0, b); M, r = A.shape
    X = 3*M
    ps = []
    while len(ps) < r:
        p = rng.randrange(M**4, 2*M**4)
        if p % 4 == 1 and is_prime(p): ps.append(p)
    hi = 4*X
    while True:
        zs = gaussian_primes_norm_range(X, hi)
        if len(zs) >= r: zs = zs[:r]; break
        hi *= 2
    D1 = np.zeros((M, M)); Dg = np.zeros((M, M))
    A64 = A.astype(np.int64)
    for q in primes_upto(X)[1:]:
        Aq = 1
        while q**(Aq+1) <= X: Aq += 1
        Q = q**Aq; f, _ = sylow2_generator(q, Aq, rng)
        g, table, m = T_generator_and_log(q, Aq, rng2)
        ell = np.zeros(r, dtype=np.int64)
        for j, (x, y) in enumerate(zs):
            eps = 0 if legendre(x*x+y*y, q) == legendre(ps[j], q) else 1
            z = (x % Q, y % Q)
            t = mul(mul(z, inv(conj(z, Q), Q), Q), pw(f, eps, Q), Q)
            ell[j] = table[t]
            assert ell[j] % 2 == (0 if legendre(ps[j], q) == 1 else 1)   # forced parity
        kap = (A64 @ ell) % m
        kap = (kap % 2) + 2*np.array([rng3.randrange(m//2) for _ in range(M)], dtype=np.int64)
        diff = (kap[:, None] - kap[None, :]) % m
        lq = math.log(q)
        for lev in range(1, Aq+1):
            ml = tq_order(q, lev)
            D1 += lq * (diff % ml == 0)
            Dg += lq * (diff % (ml//4) == 0)
    iu = np.triu_indices(M, 1)
    lM = math.log(M); d1, dg = D1[iu], Dg[iu]
    Dsym = np.triu(D1, 1); Dsym = Dsym + Dsym.T
    lam = np.linalg.eigvalsh(Dsym).max()
    return M, r, d1.max()/lM, d1.mean()/lM, dg.max()/lM, dg.mean()/lM, lam/(M*lM)

if __name__ == "__main__":
    qs = [int(a) for a in sys.argv[1:]] or [43, 59]
    print("  M     r  | unit-1: max  mean | some unit: max mean | lambda_max(D1)/(M logM) | logM/loglogM")
    for q0 in qs:
        M, r, a, b_, c, d, lam = run(q0)
        print(f"{M:5d} {r:6d} |   {a:6.2f} {b_:6.3f}  |   {c:6.2f} {d:6.3f}     |   {lam:6.3f}               | "
              f"{math.log(M)/math.log(math.log(M)):5.2f}", flush=True)
