"""Lemma 3.1/3.2 (root-of-unity parity repair): exhaustive-ish exact checks over odd prime powers Q<=Qmax."""
import random, sys
from gauss import *

def sqrt_mod_prime_power(n, q, a):
    Q = q**a; n %= Q
    # Tonelli-Shanks mod q
    r = None
    for x in range(q):
        if (x*x - n) % q == 0: r = x; break
    if r is None or r == 0: return None
    # Hensel lift
    for e in range(1, a):
        mod = q**(e+1)
        r = (r - (r*r - n) * pow(2*r, -1, mod)) % mod
    return r % Q

def run(Qmax=400, trials=60, seed=1):
    rng = random.Random(seed)
    nQ = ncoll = ndeg = 0
    for q in primes_upto(Qmax)[1:]:
        a = 1
        while q**a <= Qmax:
            Q = q**a; m = tq_order(q, a); k = v2(m)
            f, kk = sylow2_generator(q, a, rng)
            assert kk == k and k >= 2
            # f is a non-square (odd label):
            assert pw(f, m//2, Q) == ((-1) % Q, 0)
            # classes: k=2 iff q=+-3 mod 8, k=3 iff q = 7,9 mod 16
            assert (k == 2) == (q % 8 in (3, 5)) and (k == 3) == (q % 16 in (7, 9))
            # if k == 2, f is +-i (labels of i are odd)
            if k == 2: assert f in ((0, 1), (0, Q-1))
            for _ in range(trials):
                cols = []
                for j in range(5):
                    while True:
                        z = (rng.randrange(-30, 31), rng.randrange(-30, 31))
                        if z[1] != 0 and (z[0]*z[0]+z[1]*z[1]) % q != 0: break
                    while True:
                        p = rng.randrange(5, 10**6)
                        if p % 4 == 1 and p != q and is_prime(p): break
                    eps = 0 if legendre(z[0]**2+z[1]**2, q) == legendre(p, q) else 1
                    zz = (z[0] % Q, z[1] % Q)
                    t = mul(mul(zz, inv(conj(zz, Q), Q), Q), pw(f, eps, Q), Q)
                    # norm-one and parity = (p/q)
                    assert mul(t, conj(t, Q), Q) == (1, 0)
                    assert pw(t, m//2, Q) == ((1 if legendre(p, q) == 1 else Q-1), 0)
                    # realise r with r/conj r = t and Norm r = p mod Q
                    r0 = None
                    for c0 in range(1, 60):          # Hilbert 90: r0 = c + t*conj(c)
                        c = (c0, rng.randrange(Q))
                        cand = mul(t, conj(c, Q), Q); cand = ((cand[0]+c[0]) % Q, (cand[1]+c[1]) % Q)
                        if norm(cand, Q) % q: r0 = cand; break
                    assert mul(r0, inv(conj(r0, Q), Q), Q) == t
                    s = sqrt_mod_prime_power(p * pow(norm(r0, Q), -1, Q), q, a)
                    assert s is not None
                    r = (r0[0]*s % Q, r0[1]*s % Q)
                    assert norm(r, Q) == p % Q and mul(r, inv(conj(r, Q), Q), Q) == t
                    cols.append((z, t))
                # random monomials (and a forced degenerate one)
                for rep in range(40):
                    v = [rng.randrange(-2, 3) for _ in cols]
                    if not any(v): continue
                    T = (1, 0)
                    al, be = 1+0j, 1+0j
                    A = (1, 0); B = (1, 0)
                    for (z, t), e in zip(cols, v):
                        T = mul(T, pw(t, e, Q), Q)
                        if e > 0:
                            for _ in range(e): A = (A[0]*z[0]-A[1]*z[1], A[0]*z[1]+A[1]*z[0])
                        elif e < 0:
                            for _ in range(-e): B = (B[0]*z[0]-B[1]*z[1], B[0]*z[1]+B[1]*z[0])
                    g = (A[0]*B[0]+A[1]*B[1], A[1]*B[0]-A[0]*B[1])   # alpha * conj(beta)
                    G = (1, 0)
                    for _ in range(1 << k): G = (G[0]*g[0]-G[1]*g[1], G[0]*g[1]+G[1]*g[0])
                    if pw(T, 4, Q) == (1, 0):       # collision with some unit at level a
                        ncoll += 1
                        assert G[1] % Q == 0, (q, a, v)
            nQ += 1
            a += 1
    # degenerate monomial check at k=2: Z_v = 1 => collides at every q=+-3 mod 8 (i-repair)
    for q in primes_upto(Qmax)[1:]:
        if q % 8 not in (3, 5): continue
        for e in range(-3, 4):
            T = pw((0, 1), e, q)
            assert pw(T, 4, q) == (1, 0); ndeg += 1
    print("prime powers checked:", nQ, " collisions tested (=> Q | Im(gamma^(2^k))):", ncoll,
          " degenerate i-repair checks:", ndeg)
    print("ALL ROOT-REPAIR CHECKS PASSED")

if __name__ == "__main__":
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 400)
