"""Referee check (arith lens) of sharp.md Lemma 3.1, Lemma 3.3(1)-(2) and the pair-demand claims.
Independent numpy code (does not import check_design.py).

Part A (Lemma 3.1, exact up to 2^(KMAX+1)): the dyadic assignment q -> p'(q); per block k the
   exact maximum multiplicity; 2 p'(q) <= m_q, q < 8 p'(q), p'(q) <= (q-1)/2 (q >= 5), p' != q;
   total preimages per target prime (all blocks + special cases).
Part B (Lemma 3.1(2) for k >= 20): R(k) recomputed, and the R-S inputs pi(x) > x/ln x (x >= 17),
   pi(x) < 1.25506 x/ln x tested on the sieve range.
Part C (Lemma 3.1(3)): Lambda'(d) <= 9 log d + 18.72 omega(d) for every 2 <= d <= DMAX; Robin's
   omega(d) <= 1.3841 log d/loglog d (d >= 3) on the same range; and the chain used in Lemma 3.3(1):
   max_{d < M} (Lambda'(d) + log d) <= mu_M log M,  mu_M = 10 + 26/loglog M.
Part D (pair demand, exact): for M in a list and X = 3M (A = 1), X = 3M^2 (A = 2) and X = infinity,
   D(d) = sum_{q < 4M, q <= X, p'(q) | d} log q * min(a*(q), A_q, 1 + v_q(d)),
   a*(q) = min{a >= 1 : p'(q) q^(a-1) >= M}, worst case over parity colourings (sigma ignored);
   report max_d D(d)/log M.
Part E (Lemma 3.3(2)): log m_struct = sum_{q < 4M} a*(q) log q  vs  pi(4M) log(32 M^2)  vs 11 M,
   and max q^(a*(q)) / M^2 (claimed < 32).
"""
import math, sys
import numpy as np

KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 26
DMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 2_000_000

LIM = 2 ** (KMAX + 1)
sieve = np.ones(LIM + 1, dtype=bool); sieve[:2] = False
for i in range(2, int(LIM ** 0.5) + 1):
    if sieve[i]:
        sieve[i * i::i] = False
primes = np.nonzero(sieve)[0].astype(np.int64)
print(f"sieve to 2^{KMAX+1}: {len(primes)} primes")

# ---------------- Part A
pp_q, pp_p = [3, 5, 7], [2, 2, 3]
blockmax = {}
okA = True
for k in range(3, KMAX + 1):
    B = primes[(primes >= 2 ** k) & (primes < 2 ** (k + 1)) & (primes > 7)]
    T = primes[(primes >= 2 ** (k - 2)) & (primes < 2 ** (k - 1))]
    n, t = len(B), len(T)
    idx = (np.arange(n, dtype=np.int64) * t) // n
    tgt = T[idx]
    cnt = np.bincount(idx, minlength=t)
    blockmax[k] = (int(cnt.max()), n / t, int(math.ceil(n / t)))
    okA &= int(cnt.max()) <= math.ceil(n / t)
    pp_q.append(B); pp_p.append(tgt)
Q = np.concatenate([np.array(pp_q[:3])] + pp_q[3:])
P = np.concatenate([np.array(pp_p[:3])] + pp_p[3:])
mq = np.where(Q % 4 == 1, Q - 1, Q + 1)
okA &= bool(np.all(2 * P <= mq)) and bool(np.all(Q < 8 * P)) and bool(np.all(P != Q))
okA &= bool(np.all(P[Q >= 5] <= (Q[Q >= 5] - 1) // 2))
# every odd prime covered exactly once
okA &= len(Q) == len(primes) - 1 and len(set(Q.tolist())) == len(Q)
tot = np.bincount(P)
mult_all = int(tot.max())
print("Part A: exact max multiplicity per block k (max, n/t, ceil(n/t)):")
print("   ", {k: blockmax[k][0] for k in blockmax})
print(f"   overall max preimages of one target prime (incl. special cases 3,5->2, 7->3): {mult_all}"
      f"  [prime 2: {tot[2]}, prime 3: {tot[3]}]")
print(f"   2p' <= m_q, q < 8p', p' <= (q-1)/2 (q>=5), every odd prime < 2^{KMAX+1} assigned once: {okA}")
print(f"   max n/t over blocks: {max(v[1] for v in blockmax.values()):.4f}")

# ---------------- Part B
okB = True
Rmax, kR = 0, None
for k in range(20, 5001):
    eps = math.log(2) * 2.0 ** (-k)
    num = 2 * 1.25506 / (k + 1) - 1 / k + eps
    den = 0.5 / (k - 1) - 0.25 * 1.25506 / (k - 2) - eps
    if num / den > Rmax:
        Rmax, kR = num / den, k
lim = (2 * 1.25506 - 1) / (0.5 - 0.25 * 1.25506)
print(f"Part B: max_(20<=k<=5000) R(k) = {Rmax:.4f} at k = {kR};  limit k->inf = {lim:.4f};  all < 9: {max(Rmax, lim) < 9}")
xs = np.unique(np.round(np.logspace(math.log10(17), math.log10(LIM), 4000)).astype(np.int64))
pic = np.cumsum(sieve)
lo = bool(np.all(pic[xs] > xs / np.log(xs)))
hi = bool(np.all(pic[xs[xs > 1]] < 1.25506 * xs[xs > 1] / np.log(xs[xs > 1])))
# also at all powers of two
p2 = np.array([2 ** j for j in range(5, KMAX + 2)])
lo &= bool(np.all(pic[p2] > p2 / np.log(p2)))
hi &= bool(np.all(pic[p2] < 1.25506 * p2 / np.log(p2)))
print(f"   R-S inputs on the sieve range: pi(x) > x/ln x: {lo};  pi(x) < 1.25506 x/ln x: {hi}")
okB &= lo and hi and max(Rmax, lim) < 9

# ---------------- Part C
F = np.zeros(DMAX + 1)
sel = P <= DMAX
np.add.at(F, P[sel], np.log(Q[sel].astype(float)))
# targets <= DMAX need all q < 8 DMAX assigned
assert 8 * DMAX < LIM
Lam = np.zeros(DMAX + 1)
omega = np.zeros(DMAX + 1, dtype=np.int64)
for p in primes[primes <= DMAX]:
    Lam[p::p] += F[p]
    omega[p::p] += 1
d = np.arange(DMAX + 1, dtype=float)
okC = True
dd = d[2:]
bound = 9 * np.log(dd) + 18.72 * omega[2:]
viol = np.nonzero(Lam[2:] > bound + 1e-9)[0]
okC &= len(viol) == 0
ratio = Lam[2:] / np.log(dd)
print(f"Part C: Lambda'(d) <= 9 log d + 18.72 omega(d) for 2 <= d <= {DMAX}: {len(viol) == 0}"
      f"  (max Lambda'(d)/log d = {ratio.max():.3f} at d = {int(np.argmax(ratio)) + 2})")
d3 = d[3:]
rob = omega[3:] <= 1.3841 * np.log(d3) / np.log(np.log(d3)) + 1e-12
print(f"   Robin omega(d) <= 1.3841 log d/loglog d on 3..{DMAX}: {bool(np.all(rob))}")
okC &= bool(np.all(rob))
# the 'per-d' form (9 + 26/loglog d) log d is NOT monotone in d; compare with max over d < M
print("   M, max_{d<M}(Lambda'(d)+log d)/log M, mu_M, max_{3<=d<M} (9+26/loglog d) log d / log M:")
for M in [2 ** j for j in range(4, 21)]:
    if M > DMAX:
        break
    val = (Lam[1:M] + np.log(d[1:M])).max() / math.log(M)
    mu = 10 + 26 / math.log(math.log(M))
    perd = ((9 + 26 / np.log(np.log(d[3:M]))) * np.log(d[3:M])).max() / math.log(M)
    okC &= val <= mu
    print(f"   M = {M:8d}: {val:6.3f}  mu_M = {mu:6.3f}   per-d Robin form: {perd:9.3f}")

mo = np.maximum.accumulate(omega)
Ms = np.arange(30, DMAX + 1)
robmax = bool(np.all(mo[Ms - 1] <= 1.3841 * np.log(Ms) / np.log(np.log(Ms))))
print(f"   max_(d<M) omega(d) <= 1.3841 log M/loglog M for all 30 <= M <= {DMAX}: {robmax}")
bad = [M for M in range(3, 30) if mo[M - 1] > 1.3841 * math.log(M) / math.log(math.log(M))]
print(f"   ... fails for M in {bad}")
okC &= robmax

# ---------------- Part D
pp = dict(zip(Q.tolist(), P.tolist()))


def demand(M, X):
    Dv = np.zeros(M)            # index d = 0..M-1
    for q in primes[(primes > 2) & (primes < 4 * M)].tolist():
        if q > X:
            break
        p1 = pp[q]
        if p1 >= M:
            continue
        lq = math.log(q)
        Aq = 1
        while Aq < 64 and q ** (Aq + 1) <= X:
            Aq += 1
        ast = 1
        while p1 * q ** (ast - 1) < M:
            ast += 1
        top = min(ast, Aq)
        # level 1: p' | d
        Dv[p1::p1] += lq
        a = 2
        while a <= top:
            step = p1 * q ** (a - 1)
            if step >= M:
                break
            Dv[step::step] += lq
            a += 1
    Dv[0] = 0
    j = int(np.argmax(Dv))
    return Dv[j], j


print("Part D: max pair demand / log M (worst case over parity colourings):")
print("        M     A=1 (X=3M)        A=2 (X=3M^2)      X=inf (check_design convention)")
for M in [44, 64, 256, 1024, 4096, 16384, 65536, 262144]:
    if 4 * M >= LIM:
        break
    r1 = demand(M, 3 * M); r2 = demand(M, 3 * M * M); r3 = demand(M, float('inf'))
    L = math.log(M)
    print(f"  {M:8d}   {r1[0]/L:6.3f} (d={r1[1]:7d})   {r2[0]/L:6.3f} (d={r2[1]:7d})   {r3[0]/L:6.3f} (d={r3[1]:7d})")

# ---------------- Part E
print("Part E: log m_struct (sum over q<4M of a*(q) log q), pi(4M) log(32M^2), 11M, max q^a*/M^2:")
for M in [44, 256, 4096, 65536, 1048576]:
    if 4 * M >= LIM:
        break
    s, mx = 0.0, 0.0
    qs = primes[(primes > 2) & (primes < 4 * M)].tolist()
    for q in qs:
        p1 = pp[q]
        ast = 1
        while p1 * q ** (ast - 1) < M:
            ast += 1
        s += ast * math.log(q)
        mx = max(mx, q ** ast / M ** 2)
    pi4 = len(qs) + 1
    print(f"   M = {M:8d}: log m_struct = {s:12.1f} = {s/M:6.3f} M;  pi(4M) log(32M^2) = {pi4*math.log(32*M*M)/M:6.3f} M;  max q^a*/M^2 = {mx:6.3f}")

ok = okA and okB and okC
print("ALL ASSIGNMENT CHECKS PASSED" if ok else "ASSIGNMENT FAILURE")
