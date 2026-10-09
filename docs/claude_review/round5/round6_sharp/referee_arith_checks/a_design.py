"""Referee check (arith lens): independent class-level implementation of the design D_A of
sharp.md Def. 3.2 on the canonical Paley profile (q0 = 107, M = 108, b = 33, r = 3531 actual primes
= 1 mod 4 above 10^6), at every odd prime power <= X, for X = 3M (A = 1) and X = 3M^2 (A = 2).
(The class <-> residue dictionary itself is checked in a_cor22.py.)

 D1  kappa^(Q) = sigma_q mod 2 at the top power; kappa^(q^a) = kappa^(q^(a+1)) mod m_(q^a);
 D2  unit-1 pair collisions at q^a happen EXACTLY when q < 4M, a <= a*(q), sigma_q(x) = sigma_q(y),
     p'(q) | d, q^(a-1) | d  (Lemma 3.3(1)); none at large q or above a*;
 D3  per pair: log m_(xy,1) <= Lambda'(d) + log d, and max over pairs / log M;
 D4  units != 1 (iota = m/4, an element of order 4; parity iota = [q = +-3 mod 8]): pair collisions
     at large q, level 1: observed count vs the bound sum 12/m_q * #pairs (Lemma 3.3(4)).
"""
import math, random
import numpy as np

rng = random.Random(5)
q0, b = 107, 33
M = q0 + 1


def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


# Paley profile (canonical assignment, alpha(inf) = a* = 0)
chi = -np.ones(q0, dtype=np.int64); chi[0] = 0
for t in range(1, q0):
    chi[(t * t) % q0] = 1
H = np.zeros((M, q0), dtype=np.int64)
for x in range(q0):
    for a in range(q0):
        H[x, a] = -chi[(a - x) % q0] - (1 if a == x else 0)
H[q0, :] = 1
S = np.repeat(H, b, axis=1).copy()
for x in range(q0):
    S[x, x * b] *= -1
S[q0, 0 * b + 2] *= -1
A = (1 + S) // 2
r = A.shape[1]
ps = []
n = 10 ** 6 + 1
while len(ps) < r:
    if n % 4 == 1 and is_prime(n):
        ps.append(n)
    n += 4
ps = np.array(ps, dtype=np.int64)

# assignment p'
LIM = 4 * 3 * M * M + 10
sv = np.ones(LIM + 1, dtype=bool); sv[:2] = False
for i in range(2, int(LIM ** 0.5) + 1):
    if sv[i]:
        sv[i * i::i] = False
pr = np.nonzero(sv)[0]
pp = {3: 2, 5: 2, 7: 3}
for k in range(3, 40):
    if 2 ** k > LIM // 2:
        break
    B = pr[(pr >= 2 ** k) & (pr < 2 ** (k + 1)) & (pr > 7)]
    T = pr[(pr >= 2 ** (k - 2)) & (pr < 2 ** (k - 1))]
    for i, q in enumerate(B.tolist()):
        pp[q] = int(T[(i * len(T)) // len(B)])


def mq(q, a):
    return (q - (1 if q % 4 == 1 else -1)) * q ** (a - 1)


def legendre_vec(q):
    return np.array([0 if pow(int(p) % q, (q - 1) // 2, q) == 1 else 1 for p in ps], dtype=np.int64)


def run(X):
    ok = {'D1': True, 'D2': True, 'D3': True}
    rows = np.arange(M)
    dmat = np.abs(rows[:, None] - rows[None, :])
    iu = np.triu_indices(M, 1)
    dem = np.zeros(len(iu[0]))
    off_obs, off_bound = 0, 0.0
    npairs = len(iu[0])
    for q in pr[(pr > 2) & (pr <= X)].tolist():
        Aq = 1
        while q ** (Aq + 1) <= X:
            Aq += 1
        nu = legendre_vec(q)
        sig = (A @ nu) % 2
        m1 = mq(q, 1)
        small = q < 4 * M
        if small:
            p1 = pp[q]
            k1 = (sig + 2 * (rows % p1)) % m1
            ast = 1
            while p1 * q ** (ast - 1) < M:
                ast += 1
        else:
            y = np.array(rng.sample(range(m1 // 2), M))
            k1 = (sig + 2 * y) % m1
            ast = 1
        kap = {1: k1}
        for a in range(2, Aq + 1):
            ma, mp = mq(q, a), mq(q, a - 1)
            if small and a <= ast:
                qa1 = q ** (a - 1)
                # CRT: z = k1 mod m1, z = x mod q^(a-1)
                u = pow(m1, -1, qa1)
                kap[a] = (k1 + m1 * (((rows % qa1) - k1) * u % qa1)) % ma
            else:
                kap[a] = (kap[a - 1] + mp * np.array([rng.randrange(q) for _ in range(M)])) % ma
        # D1
        ok['D1'] &= bool(np.all((kap[Aq] - sig) % 2 == 0))
        for a in range(1, Aq):
            ok['D1'] &= bool(np.all(kap[a + 1] % mq(q, a) == kap[a]))
        # D2 / D3
        for a in range(1, Aq + 1):
            ma = mq(q, a)
            coll = (kap[a][:, None] - kap[a][None, :]) % ma == 0
            if small and a <= ast:
                pred = (sig[:, None] == sig[None, :]) & (dmat % p1 == 0) & (dmat % q ** (a - 1) == 0)
            else:
                pred = np.zeros((M, M), dtype=bool)
            np.fill_diagonal(pred, True)
            ok['D2'] &= bool(np.array_equal(coll, pred))
            dem += math.log(q) * coll[iu]
        # D4: units != 1 at large q, level 1
        if not small:
            iota = m1 // 4
            diff = (kap[1][:, None] - kap[1][None, :]) % m1
            for s in (1, 2, 3):
                off_obs += int(np.sum(diff[iu] == (iota * s) % m1))
            off_bound += 12 / m1 * npairs
    # D3: bound Lambda'(d) + log d
    Lam = np.zeros(M)
    for q, p1 in pp.items():
        if p1 < M:
            Lam[p1::p1] += math.log(q)
    dd = dmat[iu]
    bound = Lam[dd] + np.log(dd)
    ok['D3'] &= bool(np.all(dem <= bound + 1e-9))
    print(f"X = {X}: D1..D3 {ok};  max pair demand = {dem.max():.2f} = {dem.max()/math.log(M):.3f} log M;"
          f"  unit!=1 large-q pair collisions: observed {off_obs}, bound {off_bound:.1f}")
    return all(ok.values())


if __name__ == "__main__":
    allok = run(3 * M) & run(3 * M * M)
    print("ALL DESIGN (CLASS-LEVEL) CHECKS PASSED" if allok else "DESIGN FAILURE")
