"""Referee check of sharp.md Def. 3.2 / Lemma 3.3(1),(4) at the class level, on a larger profile than
the authors' M = 44 test, with genuine Legendre colourings.

Own Paley-I builder (q0 = 3 mod 4, M = q0 + 1), b copies, canonical one-flip assignment; twin
structure a_(.,f(x)) - a_(.,s(x)) = eps_x e_x is asserted.  r = b(M-1) actual primes p_j = 1 mod 4
above X.  For every odd prime q <= X = 3 M^A: sigma_q = A nu(q) mod 2, the design D_A (small q < 4M:
sigma + 2 (x mod p'), CRT levels up to a*, random lifts above; large q: sigma + 2 y, y random
injection; random lifts), iota = m/4 or 3m/4 by q mod 8 (any choice with the right parity class;
the collision pattern with unit 1 does not depend on it).  Then for every pair and every level:
  * unit-1 collisions == predicted (small q, a <= a*, sigma equal, p'|d, q^(a-1)|d);
  * the exact pair demand log m_(xy,1) and its max over pairs vs mu_M log M and the worst case;
  * unit != 1 pair collisions (m^off) recorded.
Only class arithmetic is used here; the residue-level realisation is ref_cor22.py Part B."""
import math, random, sys
import numpy as np

q0 = int(sys.argv[1]) if len(sys.argv) > 1 else 59
b = 33
A = float(sys.argv[2]) if len(sys.argv) > 2 else 2.0
rng = random.Random(7)
M = q0 + 1
X = int(3 * M ** A)


def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    i = 3
    while i * i <= n:
        if n % i == 0:
            return False
        i += 2
    return True


def leg(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


# Paley-I: rows 0..q0-1 finite, q0 = infinity; labels a in F_q0
H = np.zeros((M, q0), dtype=np.int64)
for x in range(q0):
    for a in range(q0):
        H[x, a] = -leg(a - x, q0) - (1 if a == x else 0)
H[q0, :] = 1
full = np.hstack([np.ones((M, 1), dtype=np.int64), H])
assert (full.T @ full == M * np.eye(M, dtype=np.int64)).all()
r = b * q0
col_label = np.repeat(np.arange(q0), b)
S = H[:, col_label].copy()
flip, sib = {}, {}
astar = 0
for x in range(q0):
    flip[x], sib[x] = x * b + 0, x * b + 1
flip[q0], sib[q0] = astar * b + 2, astar * b + 3
for x in range(M):
    S[x, flip[x]] *= -1
Amat = (1 + S) // 2
eps = {}
for x in range(M):
    dcol = Amat[:, flip[x]] - Amat[:, sib[x]]
    assert np.count_nonzero(dcol) == 1 and dcol[x] in (1, -1)
    eps[x] = int(dcol[x])
    assert eps[x] == -H[x, x if x < q0 else astar]
assert len(set(flip.values()) | set(sib.values())) == 2 * M
assert np.linalg.matrix_rank(Amat.astype(float)) == M

# primes p_j > X, = 1 mod 4
plist = []
n = X + 1
while len(plist) < r:
    if n % 4 == 1 and is_prime(n):
        plist.append(n)
    n += 1

# assignment p'(q)
def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0].tolist()


Pall = primes_upto(max(X, 64) * 2)
pp = {3: 2, 5: 2, 7: 3}
k = 3
while 2 ** k <= 4 * M:
    Bk = [q for q in Pall if 2 ** k <= q < 2 ** (k + 1) and q > 7]
    Tk = [p for p in Pall if 2 ** (k - 2) <= p < 2 ** (k - 1)]
    for i, q in enumerate(Bk):
        pp[q] = Tk[(i * len(Tk)) // len(Bk)]
    k += 1


def mT(q, a):
    return (q - (1 if q % 4 == 1 else -1)) * q ** (a - 1)


oddq = [q for q in Pall if 2 < q <= X]
logm1 = np.zeros((M, M))
logmoff = np.zeros((M, M))
ok = True
nlev = 0
xs = np.arange(M)
for q in oddq:
    Aq = 1
    while q ** (Aq + 1) <= X:
        Aq += 1
    nu = np.array([1 if leg(p, q) == -1 else 0 for p in plist], dtype=np.int64)
    sig = (Amat @ nu) % 2
    m1 = mT(q, 1)
    small = q < 4 * M
    if small:
        p_ = pp[q]
        k1 = (sig + 2 * (xs % p_)) % m1
        ast = 1
        while p_ * q ** (ast - 1) < M:
            ast += 1
    else:
        y = np.array(rng.sample(range(m1 // 2), M))
        k1 = (sig + 2 * y) % m1
        ast = 1
    kap = {1: k1}
    for a in range(2, Aq + 1):
        ma, mp = mT(q, a), mT(q, a - 1)
        if small and a <= ast:
            qa1 = q ** (a - 1)
            inv = pow(m1, -1, qa1)
            t = ((xs - k1) * inv) % qa1
            kap[a] = (k1 + m1 * t) % ma
        else:
            kap[a] = (kap[a - 1] + mp * np.array([rng.randrange(q) for _ in range(M)])) % ma
        ok &= bool(np.all(kap[a] % mp == kap[a - 1]))
    ok &= bool(np.all((kap[Aq] - sig) % 2 == 0))
    lq = math.log(q)
    for a in range(1, Aq + 1):
        ma = mT(q, a)
        iota = ma // 4 if q % 8 in (1, 7) else (ma // 4 if (ma // 4) % 2 else 3 * ma // 4)
        iota = ma // 4 if ((ma // 4) % 2 == (1 if q % 8 in (3, 5) else 0)) else 3 * ma // 4
        D = (kap[a][:, None] - kap[a][None, :]) % ma
        c1 = D == 0
        coff = (D == iota % ma) | (D == (2 * iota) % ma) | (D == (3 * iota) % ma)
        dd = np.abs(xs[:, None] - xs[None, :])
        if small and a <= ast:
            pred = (sig[:, None] == sig[None, :]) & (dd % pp[q] == 0) & (dd % q ** (a - 1) == 0)
        else:
            pred = np.zeros((M, M), dtype=bool)
        np.fill_diagonal(pred, True)
        ok &= bool(np.array_equal(c1, pred))
        logm1 += np.where(c1, lq, 0.0)
        logmoff += np.where(coff, lq, 0.0)
        nlev += 1
np.fill_diagonal(logm1, 0)
np.fill_diagonal(logmoff, 0)
mu = 10 + 26 / math.log(math.log(M))
print(f"q0 = {q0}, M = {M}, b = {b}, r = {r}, X = 3 M^{A:g} = {X}, odd primes <= X: {len(oddq)}, levels: {nlev}")
print(f"unit-1 pair collisions exactly as Lemma 3.3(1) predicts at every level: {ok}")
print(f"max pair demand log m_(xy,1) = {logm1.max():.2f} = {logm1.max()/math.log(M):.3f} log M  (mu_M log M = {mu*math.log(M):.1f});"
      f" mean = {logm1[np.triu_indices(M,1)].mean()/math.log(M):.3f} log M")
print(f"max unit != 1 pair modulus log m^off = {logmoff.max():.2f}; (W/4 - 1 with tau ~ 2 log M: {r*2*math.log(M)/4 - 1:.0f})")
print("ALL DESIGN-PAIR CHECKS PASSED" if ok else "FAILURE")
