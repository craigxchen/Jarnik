"""Evidence for Sections A and C: demand statistics of three residue models.

Models (all with unit-1 pair collisions at odd primes q <= X = 3M and their powers <= X):
  RAND : independent uniform classes mod m_Q (random lifts at prime powers);
  DA   : the round-6 design D_A (sharp.md Def. 3.2, Lemma 3.1 prime assignment p'(q)), with an
         independent uniformly random parity colouring sigma_q per prime (worst case is similar);
  GEN  : genuine coherent residues: e = 1 profile with random +-1 rows (r = 4M columns, rows
         distinct), genuine Gaussian primes pi_j with p_j = 1 mod 4 random in [10^6, 2*10^6],
         points zhat_x = prod pi_j^(a_xj) conj(pi_j)^(1-a_xj), k = 0; collision mod q iff
         zhat_x = zhat_y mod q (level 1 only; computed in Z[i]/q exactly).
Reported: D_max/log M (max pair demand), lam_max(D)/(M log M), max_x D_x/(M log M),
avg_x D_x/(M log M), and for the window Q = primes in (M^(1/2), M]: max pair collision count
c_max and the ratio of factorial moments  F_k = sum_{x<y} binom(c_xy, k)  to the independent
expectation  binom(M,2) e_k(1/m_q : q in Q)  (Hypothesis H(eta=1/2, delta) asks F_k >= delta^k E_k).
Evidence only.
"""
import math, sys
import numpy as np
from math import comb

rng = np.random.default_rng(31415)

def primes_upto(n):
    s = bytearray([1]) * (n + 1); s[0:2] = b'\x00\x00'
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]

def mq(q):
    return q - 1 if q % 4 == 1 else q + 1

def pprime_map(X):
    """sharp.md Lemma 3.1."""
    P = primes_upto(max(64, X))
    mp = {3: 2, 5: 2, 7: 3}
    k = 3
    while 2 ** k <= X:
        qs = [q for q in P if 2 ** k <= q < 2 ** (k + 1) and q > 7]
        ts = [t for t in P if 2 ** (k - 2) <= t < 2 ** (k - 1)]
        n, T = len(qs), len(ts)
        for i, q in enumerate(qs):
            mp[q] = ts[(i * T) // n]
        k += 1
    return mp

def esym(vals, kmax):
    e = [1.0] + [0.0] * kmax
    for v in vals:
        for k in range(kmax, 0, -1):
            e[k] += e[k - 1] * v
    return e

def classes_rand(M, q, X):
    """list over levels a (q^a <= X) of class arrays."""
    out = []
    cls = rng.integers(0, mq(q), size=M)
    Q = q
    while Q <= X:
        out.append(cls.copy())
        cls = cls * q + rng.integers(0, q, size=M)
        Q *= q
    return out

def classes_DA(M, q, X, pp):
    x = np.arange(M)
    sigma = rng.integers(0, 2, size=M)
    p = pp[q]
    out = []
    lvl1 = sigma + 2 * (x % p)
    Q, a = q, 1
    astar = 1
    while p * q ** (astar - 1) < M:
        astar += 1
    cls = lvl1.copy()
    while Q <= X:
        if a == 1:
            cls = lvl1.copy()
        elif a <= astar:
            cls = lvl1 * (q ** (a - 1)) + (x % (q ** (a - 1)))
        else:
            cls = cls * q + rng.integers(0, q, size=M)
        out.append(cls.copy())
        a += 1; Q *= q
    return out

# ---- Gaussian integers mod q (pairs of ints) ----
def gmul(a, b, q):
    return ((a[0] * b[0] - a[1] * b[1]) % q, (a[0] * b[1] + a[1] * b[0]) % q)

def gauss_prime(p):
    # p = 1 mod 4: find x^2 + y^2 = p
    for x in range(1, int(math.isqrt(p)) + 1):
        y2 = p - x * x
        y = math.isqrt(y2)
        if y * y == y2:
            return (x, y)
    raise ValueError

def classes_GEN(A, pis, q):
    M, r = A.shape
    out = []
    for x in range(M):
        z = (1, 0)
        for j in range(r):
            a, b = pis[j]
            f = (a % q, b % q) if A[x, j] else (a % q, (-b) % q)
            z = gmul(z, f, q)
        out.append(z[0] * q + z[1])
    return [np.array(out)]

def stats(M, model, X, A=None, pis=None, pp=None):
    D = np.zeros((M, M))
    cnt = np.zeros((M, M), dtype=int)
    lo = math.sqrt(M)
    Qw = []
    for q in primes_upto(X):
        if q == 2:
            continue
        if model == 'RAND':
            levels = classes_rand(M, q, X)
        elif model == 'DA':
            levels = classes_DA(M, q, X, pp)
        else:
            levels = classes_GEN(A, pis, q)
        for a, cls in enumerate(levels):
            same = (cls[:, None] == cls[None, :])
            D += 2 * math.log(q) * same
            if a == 0 and lo < q <= M:
                cnt += same
        if lo < q <= M:
            Qw.append(q)
    np.fill_diagonal(D, 0); np.fill_diagonal(cnt, 0)
    L = math.log(M)
    lam = max(np.linalg.eigvalsh(D))
    rows = D.sum(axis=1)
    iu = np.triu_indices(M, 1)
    c = cnt[iu]
    cmax = int(c.max())
    kmax = max(cmax, 1) + 1
    e = esym([1.0 / mq(q) for q in Qw], kmax)
    fm = []
    for k in range(1, kmax + 1):
        Fk = int(sum(comb(int(v), k) for v in c))
        Ek = comb(M, 2) * e[k]
        fm.append((k, Fk, Ek))
    return dict(Dmax=D.max() / L, lam=lam / (M * L), rowmax=rows.max() / (M * L), rowavg=rows.mean() / (M * L),
                cmax=cmax, fm=fm)

def fmt(name, M, s):
    fms = " ".join(f"k={k}:{F}/{E:.2f}" for (k, F, E) in s['fm'])
    return (f"{name:4s} M={M:5d} D_max/logM={s['Dmax']:6.2f} lam/(MlogM)={s['lam']:6.3f} "
            f"maxD_x/(MlogM)={s['rowmax']:6.3f} avgD_x/(MlogM)={s['rowavg']:6.3f} | window c_max={s['cmax']} | F_k/E_k: {fms}")

def main():
    for M in (64, 128, 256, 512, 1024):
        X = 3 * M
        pp = pprime_map(X)
        for model in ('RAND', 'DA'):
            s = stats(M, model, X, pp=pp)
            print(fmt(model, M, s)); sys.stdout.flush()
    # genuine coherent residues (smaller M; exact Gaussian arithmetic)
    P = [p for p in primes_upto(2_000_000) if p >= 1_000_000 and p % 4 == 1]
    for M in (32, 64, 128):
        r = 4 * M
        while True:
            A = rng.integers(0, 2, size=(M, r))
            if len({tuple(row) for row in A}) == M and np.all(A.min(axis=0) == 0) and np.all(A.max(axis=0) == 1):
                break
        ps = rng.choice(P, size=r, replace=False)
        pis = [gauss_prime(int(p)) for p in ps]
        s = stats(M, 'GEN', 3 * M, A=A, pis=pis)
        print(fmt('GEN', M, s)); sys.stdout.flush()
    print("done (evidence only)")

if __name__ == '__main__':
    main()
