"""Evidence for round6/sharp.md §6 (class P3-all: factored residues, all monomials strengthened).

Part A (exact instance, Paley q0 = 43, M = 44, b = 33, accessible moduli q^a with m < 2M):
  structured parity design for the points (Def. 3.1) + uniformly random fibre ell = nu + 2(ell'_0 + k),
  k uniform in K' = {k : A k in (Z/(m/2)) 1} (parametrised through the twins).  For every column j
  the single-prime collision modulus m_(e_j) = prod q^(a_q) (any line: <e_j, ell> in the 4-torsion
  {0, m/4, m/2, 3m/4}); also all 2-column monomials e_i +- e_j.  Prints the distribution, and checks
  that the point classes A ell are unchanged (pairs and characters unaffected by the fibre).
Part B (abstract model at larger M): r = 33 M columns, each column hit independently at each odd
  prime q < 2M with probability 2/(q - chi(q)) (two admissible lines of the right parity out of
  m/2 values); max over columns of log m_(e_j) / log M.  Shows the Poisson tail that forces
  log P >> log M for this route (Prop. 6.4: m_(e_j) <= sqrt(2 p_j) is necessary).
Part C (genuine-coherent route): actual Gaussian primes pi_j over random p_j = 1 mod 4 near 10^9,
  point residues z_x = prod pi_j^(a_xj) conj(pi_j)^(1-a_xj); pair demand sum of log q over
  accessible q^a with z_x = z_y mod q^a; compared with the structured design."""
import math, random
import numpy as np
from gres import primes_upto, m_of, legendre, is_prime, mul, conj
from profile import build
from check_design import assignment

rng = random.Random(7)
np_rng = np.random.default_rng(7)
q0, b = 43, 33
H, S, A, flip, sib, alpha = build(q0, b)
M, r = A.shape
pp = assignment(10)
plist = []
n = 10 ** 7 + 1
while len(plist) < r:
    if n % 4 == 1 and is_prime(n):
        plist.append(n)
    n += 4
fset = set(flip.values()) | set(sib.values())
others = [j for j in range(r) if j not in fset]

acc = []                                   # accessible prime powers (q, a): m_(q^a) < 2M
for q in primes_upto(4 * M):
    if q == 2:
        continue
    a = 1
    while m_of(q, a) < 2 * M:
        acc.append((q, a))
        a += 1
top = {}
for (q, a) in acc:
    top[q] = max(top.get(q, 0), a)

logm_single = np.zeros(r)
logm_pairs_plus = np.zeros((r, r))
ok = True
for q, Atop in top.items():
    mt = m_of(q, Atop)
    half = mt // 2
    nu = np.array([1 if legendre(p, q) == -1 else 0 for p in plist], dtype=np.int64)
    sigma = (A @ nu) % 2
    # structured point design at the top accessible level (CRT with x mod q^(a-1))
    m1 = m_of(q, 1)
    k1 = [(int(sigma[x]) + 2 * (x % pp[q])) % m1 for x in range(M)]
    qa1 = q ** (Atop - 1)
    kap = [(k1[x] + m1 * ((x - k1[x]) * pow(m1, -1, qa1) % qa1)) % mt if Atop > 1 else k1[x] for x in range(M)]
    yv = [((kap[x] - int((A[x] @ nu))) % mt) // 2 for x in range(M)]
    # particular solution via twins
    l0 = np.zeros(r, dtype=np.int64)
    for x in range(M):
        h = int(H[x, alpha[x]])
        l0[flip[x]] = (l0[flip[x]] - h * yv[x]) % half
        l0[sib[x]] = (l0[sib[x]] + h * yv[x]) % half
    # uniform k in K' through the twin parametrisation
    k = np.zeros(r, dtype=np.int64)
    k[others] = np_rng.integers(0, half, size=len(others))
    u = np_rng.integers(0, half, size=M)
    t = int(np_rng.integers(0, half))
    base = (A[:, others] @ k[others]) % half
    for x in range(M):
        k[sib[x]] = 0
    sib_contrib = np.array([sum(int(A[x, sib[y]]) * int(u[y]) for y in range(M)) for x in range(M)]) % half
    for x in range(M):
        h = int(H[x, alpha[x]])
        k[flip[x]] = (h * (int(base[x]) + int(sib_contrib[x]) - t)) % half
        k[sib[x]] = (int(u[x]) - int(k[flip[x]])) % half
    ok &= all(int((A[x] @ k) % half) == t for x in range(M))
    ell = (nu + 2 * ((l0 + k) % half)) % mt
    kap_check = (A @ ell) % mt
    ok &= all(int(kap_check[x]) == (kap[x] + 2 * t) % mt for x in range(M))
    for a in range(1, Atop + 1):
        ma = m_of(q, a)
        tors = {0, ma // 4, ma // 2, 3 * ma // 4}
        la = ell % ma
        hit = np.array([int(v) in tors for v in la])
        logm_single[hit] += math.log(q)
        # 2-column monomials e_i + e_j and e_i - e_j (any line)
        s = (la[:, None] + la[None, :]) % ma
        d = (la[:, None] - la[None, :]) % ma
        hp = np.isin(s, list(tors)) | np.isin(d, list(tors))
        logm_pairs_plus[hp] += math.log(q)

lm = logm_single / math.log(M)
twin_cols = sorted(fset)
print(f"Part A: M = {M}, r = {r}, accessible prime powers: {len(acc)}")
print(f"  single columns: mean log m/log M = {lm.mean():.2f}, max = {lm.max():.2f} (twin columns max {lm[twin_cols].max():.2f}),"
      f" fraction with log m/log M > 3: {(lm > 3).mean():.4f}")
iu = np.triu_indices(r, 1)
lp = logm_pairs_plus[iu] / math.log(M)
print(f"  2-column monomials e_i +- e_j ({len(lp)}): mean {lp.mean():.2f}, max {lp.max():.2f}, 99.99% quantile {np.quantile(lp, 0.9999):.2f}")
print(f"  point classes unchanged by the fibre: {ok}")

# Part B
print("Part B (abstract independent model, r = 33 M columns):")
for Mb in (10 ** 2, 10 ** 3, 10 ** 4, 3 * 10 ** 4):
    rb = 33 * Mb
    acc_log = np.zeros(rb)
    for q in primes_upto(2 * Mb):
        if q == 2:
            continue
        mq = q - 1 if q % 4 == 1 else q + 1
        p_hit = min(1.0, 2.0 / (mq / 2))
        acc_log += (np_rng.random(rb) < p_hit) * math.log(q)
    print(f"  M = {Mb:6d}: max_j log m_(e_j)/log M = {acc_log.max() / math.log(Mb):.2f};  mean = {acc_log.mean() / math.log(Mb):.2f};"
          f"  log M/loglog M = {math.log(Mb) / math.log(math.log(Mb)):.2f}")

# Part C: genuine coherent residues, pair demands
def gaussian_prime(p):
    for x in range(1, int(p ** 0.5) + 1):
        y2 = p - x * x
        y = math.isqrt(y2)
        if y * y == y2:
            return (x, y)


gp = []
while len(gp) < r:
    p = rng.randrange(10 ** 9, 2 * 10 ** 9)
    if p % 4 == 1 and is_prime(p):
        gp.append(gaussian_prime(p))
gen_dem = np.zeros((M, M))
str_dem = np.zeros((M, M))
for (q, a) in acc:
    Q = q ** a
    pis = [(x % Q, y % Q) for (x, y) in gp]
    z = []
    for x in range(M):
        u = (1, 0)
        for j in range(r):
            u = mul(u, pis[j] if A[x, j] else conj(pis[j], Q), Q)
        z.append(u)
    for x in range(M):
        for y in range(x + 1, M):
            if z[x] == z[y]:
                gen_dem[x, y] += math.log(q)
            d = y - x
            if d % pp[q] == 0 and d % q ** (a - 1) == 0:
                str_dem[x, y] += math.log(q)      # structured design, worst case (all same parity)
print(f"Part C: genuine-coherent pair demand: max {gen_dem.max() / math.log(M):.2f} log M, mean {gen_dem[np.triu_indices(M,1)].mean() / math.log(M):.2f} log M;"
      f"  structured design (parity ignored, worst case): max {str_dem.max() / math.log(M):.2f} log M")
print("P3-ALL EVIDENCE DONE" if ok else "FAILURE")
