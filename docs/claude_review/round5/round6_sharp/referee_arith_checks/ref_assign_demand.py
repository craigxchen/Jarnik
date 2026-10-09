"""Referee check of sharp.md Lemma 3.1 (dyadic prime assignment), Lemma 3.3(1)-(2) (pair demand,
structured modulus).  Independent numpy implementation.

Part 1  exact: p'(q) for all odd primes q < 2^(KMAX+1); 2p' <= m_q, q < 8p', p' <= (q-1)/2 (q>=5);
        multiplicity per target; exact n/t per block.
Part 2  R(k) (Rosser-Schoenfeld) in 50-digit Decimal for 20 <= k <= 5000, and the k > 2000 tail.
Part 3  Lambda'(d) <= 9 log d + 18.72 omega(d) for all d <= DMAX (exact assignment), and the
        monotonicity of the Robin form used in Lemma 3.3(1).
Part 4  pair demand D(d) = sum_{q<4M, p'(q)|d} (1 + min(v_q(d), a*(q)-1)) log q, max over d < M, for
        M = 2^6..2^17 (worst case over parity colourings); variant restricted to q^a <= 3M (A = 1).
Part 5  log m_struct = sum_{q<4M} a*(q) log q against pi(4M) log(32 M^2) and 11 M."""
import math, sys
from decimal import Decimal, getcontext
import numpy as np

KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 25
LIM = 2 ** (KMAX + 1)
sieve = np.ones(LIM, dtype=bool)
sieve[:2] = False
for i in range(2, int(LIM ** 0.5) + 1):
    if sieve[i]:
        sieve[i * i::i] = False
P = np.nonzero(sieve)[0]
ok = True

# ---------------- Part 1
pp = {3: 2, 5: 2, 7: 3}
blk_ratio = []
for k in range(3, KMAX + 1):
    Bk = P[(P >= 2 ** k) & (P < 2 ** (k + 1)) & (P > 7)]
    Tk = P[(P >= 2 ** (k - 2)) & (P < 2 ** (k - 1))]
    n, t = len(Bk), len(Tk)
    idx = (np.arange(n, dtype=np.int64) * t) // n
    tgt = Tk[idx]
    pp.update(zip(Bk.tolist(), tgt.tolist()))
    vals, cnts = np.unique(tgt, return_counts=True)
    blk_ratio.append((k, n, t, n / t, int(cnts.max())))
qs = np.array(sorted(pp))
pv = np.array([pp[q] for q in qs])
mq = np.where(qs % 4 == 1, qs - 1, qs + 1)
ok &= bool(np.all(2 * pv <= mq)) and bool(np.all(qs < 8 * pv)) and bool(np.all(pv < qs))
ok &= bool(np.all(pv[qs >= 5] <= (qs[qs >= 5] - 1) // 2))
vals, cnts = np.unique(pv, return_counts=True)
mult = dict(zip(vals.tolist(), cnts.tolist()))
print(f"Part 1: odd primes q < 2^{KMAX+1}: {len(qs)}; max multiplicity {cnts.max()} "
      f"(targets with mult = max: {[int(v) for v, c in zip(vals, cnts) if c == cnts.max()][:8]}...)")
print("   multiplicity of 2, 3:", mult.get(2), mult.get(3), "; preimages of 2:", [int(q) for q in qs if pp[q] == 2])
print("   per block (k, n, t, n/t, max mult):")
for row in blk_ratio[-6:]:
    print("     ", row[0], row[1], row[2], f"{row[3]:.4f}", row[4])

# ---------------- Part 2
getcontext().prec = 50
ln2 = Decimal(2).ln()


def R(k):
    eps = ln2 / (Decimal(2) ** k)
    num = Decimal("2.51012") / (k + 1) - Decimal(1) / k + eps
    den = Decimal("0.5") / (k - 1) - Decimal("0.313765") / (k - 2) - eps
    return num / den


Rmax, kmax_ = max((R(k), k) for k in range(20, 5001))
print(f"Part 2: max_(20<=k<=5000) R(k) = {float(Rmax):.5f} at k = {kmax_};  R(2000) = {float(R(2000)):.5f};"
      f"  limit 1.51012/0.186235 = {1.51012/0.186235:.5f}")
# tail: for k > 2000 numerator <= 1.51013/k ; denominator >= (0.186235 - 1/k - tiny)/k
ok &= float(Rmax) < 9 and 1.51013 / (0.186235 - 1 / 2000) < 9
# sanity: R(k) really bounds the exact n/t for 20 <= k <= KMAX
for (k, n, t, rat, mm) in blk_ratio:
    if k >= 20:
        ok &= rat <= float(R(k))
print("   R(k) >= exact n/t for 20 <= k <=", KMAX, ":", all(r[3] <= float(R(r[0])) for r in blk_ratio if r[0] >= 20))

# ---------------- Part 3
DMAX = 10 ** 6
Lam = np.zeros(DMAX + 1)
for q, p in pp.items():
    if p <= DMAX:
        Lam[p::p] += math.log(q)
omega = np.zeros(DMAX + 1, dtype=np.int64)
for p in P[P <= DMAX]:
    omega[p::p] += 1
d = np.arange(2, DMAX + 1)
lhs = Lam[2:]
rhs = 9 * np.log(d) + 18.72 * omega[2:]
ok &= bool(np.all(lhs <= rhs + 1e-9))
print(f"Part 3: Lambda'(d) <= 9 log d + 18.72 omega(d) for all 2 <= d <= {DMAX}: {bool(np.all(lhs <= rhs + 1e-9))};"
      f"  max Lambda'(d)/log d = {np.max(lhs / np.log(d)):.3f} at d = {int(d[np.argmax(lhs/np.log(d))])}")
# Robin-form function used in Lemma 3.3(1): f(d) = (10 + 26/loglog d) log d is NOT increasing
f = lambda x: (10 + 26 / math.log(math.log(x))) * math.log(x)
print(f"   Robin-form f(d) = (10 + 26/loglog d) log d:  f(3) = {f(3):.1f}, f(5) = {f(5):.1f}, f(16) = {f(16):.1f},"
      f"  f(1000) = mu_1000 log 1000 = {f(1000):.1f}  (so 'f(d) <= f(M) for d < M' fails)")
# but the direct statement max_{d<M} (Lambda'(d) + log d) <= mu_M log M holds:
direct = Lam[1:] + np.log(np.arange(1, DMAX + 1))
runmax = np.maximum.accumulate(direct)
viol = [M for M in range(4, DMAX + 1) if runmax[M - 2] > f(M)]
print(f"   direct: max_(d<M)(Lambda'(d) + log d) <= (10 + 26/loglog M) log M for all 4 <= M <= {DMAX}: {len(viol) == 0}"
      f"{'' if not viol else ' first violations ' + str(viol[:5])}")

# ---------------- Part 4
def vq(x, q):
    v = 0
    while x % q == 0:
        x //= q
        v += 1
    return v


print("Part 4: max pair demand (unit 1), worst case over parity colourings")
print(f"   {'M':>7} {'max D':>8} {'at d':>7} {'D/logM':>7} {'(A=1: q^a<=3M) D/logM':>22} {'mu_M':>7}")
for M in [64, 256, 1024, 4096, 16384, 65536, 131072]:
    qsm = [q for q in pp if q < 4 * M]
    D = np.zeros(M)
    D1 = np.zeros(M)
    for q in qsm:
        p = pp[q]
        a = 1
        while p * q ** (a - 1) < M:
            a += 1
        astar = a
        lq = math.log(q)
        for lev in range(1, astar):           # level a* never collides (p' q^(a*-1) >= M > d)
            step = p * q ** (lev - 1)
            if step >= M:
                break
            D[step::step] += lq                # d multiple of p' q^(lev-1), d >= 1
            if q ** lev <= 3 * M:
                D1[step::step] += lq
    j = int(np.argmax(D[1:])) + 1
    mu = 10 + 26 / math.log(math.log(M))
    print(f"   {M:>7} {D[j]:>8.2f} {j:>7} {D[j]/math.log(M):>7.3f} {D1[1:].max()/math.log(M):>22.3f} {mu:>7.2f}")
    ok &= D[j] <= mu * math.log(M)

# ---------------- Part 5
print("Part 5: structured modulus")
for M in [44, 100, 1000, 10 ** 4, 10 ** 5]:
    qsm = [q for q in pp if q < 4 * M]
    tot = 0.0
    for q in qsm:
        a = 1
        while pp[q] * q ** (a - 1) < M:
            a += 1
        tot += a * math.log(q)
    bnd = len(qsm) * math.log(32 * M * M)
    print(f"   M = {M:>6}: log m_struct = {tot:10.1f};  pi(4M) log(32M^2) = {bnd:10.1f};  11M = {11*M}")
    ok &= tot <= bnd
print("ALL ASSIGNMENT/DEMAND CHECKS PASSED" if ok else "FAILURE")
