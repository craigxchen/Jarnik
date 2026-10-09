"""Merge referee, miscellaneous independent checks for round6/sharp.md.

G   Lemma 6.2 scope: explicit counterexample with general exponents (M = 2).
C   §7.3: exponent of #{c in Z^M : ||c||_1 <= 4M}; Mertens factor sum_(e^y<q<=e^(2y)) 1/q.
MU  Lemma 3.3(1) demand line: per-d chain vs max-omega chain vs mu_M log M; exact max demand.
D   Exact worst-case pair demand of D_A (parity-blind) for every M in [20, 2000], A = 1, 2.
K   Coherent example zeta_x = n_x + i (n_x^2+1 prime > X): unit-1 vs other-unit pair moduli.
"""
import math, itertools
import numpy as np
from fractions import Fraction

def sieve(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i*i::i] = False
    return s

LIM = 1 << 16
SV = sieve(LIM)
PR = np.nonzero(SV)[0].tolist()

# ---------------------------------------------------------------- G
def check_G():
    # rows x = 1, 2; columns: f(1)=(2,0), s(1)=(1,0), f(2)=(0,2), s(2)=(0,1), j=(G,0)
    for G in (10, 100):
        a1 = np.array([2, 1, 0, 0, G]); a2 = np.array([0, 0, 2, 1, 0])
        cols = np.stack([a1, a2])
        ok_norm = all(cols[:, j].min() == 0 and cols[:, j].max() >= 1 for j in range(5))
        twins = (list(cols[:, 0] - cols[:, 1]) == [1, 0]) and (list(cols[:, 2] - cols[:, 3]) == [0, 1])
        v = a1 - a2 - G * np.eye(5, dtype=int)[4]
        # L = Z (a1 - a2) since e_x - e_y spans the zero-sum lattice for M = 2
        base = a1 - a2
        inL = any((v == t * base).all() for t in range(-3, 4)) or (v == 0).all()
        # v in L only if v proportional to base with integral t; base has entry G at j, v has 0 -> t = 0 -> v = 0
        print(f"G: G={G}: normalised={ok_norm}, twins distinct={twins}, v={v.tolist()}, ||v||_1={np.abs(v).sum()}, "
              f"v in L? {inL}; v in L + G Z^r by construction; Lemma 6.2 would need G <= 2||v||_1 = {2*np.abs(v).sum()}: "
              f"{'VIOLATED' if G > 2*np.abs(v).sum() else 'ok'}")

# ---------------------------------------------------------------- C
def check_C():
    def logcount(M, n):
        # log sum_k 2^k binom(M,k) binom(n,k)
        terms = [k * math.log(2) + math.lgamma(M+1) - math.lgamma(k+1) - math.lgamma(M-k+1)
                 + math.lgamma(n+1) - math.lgamma(k+1) - math.lgamma(n-k+1) for k in range(0, M + 1)]
        mx = max(terms)
        return mx + math.log(sum(math.exp(t - mx) for t in terms))
    for M in (100, 1000, 10000):
        print(f"C: M={M}: (1/M) log #{{||c||_1 <= 4M}} = {logcount(M, 4*M)/M:.4f};  crude 2^n binom(M+n-1,n) exponent "
              f"{(4*M*math.log(2) + math.lgamma(5*M) - math.lgamma(4*M+1) - math.lgamma(M))/M:.4f}")
    H = lambda p: -p*math.log(p) - (1-p)*math.log(1-p)
    best = max((x*math.log(2) + H(x) + 4*H(x/4), x) for x in np.linspace(1e-4, 1-1e-4, 200001))
    print(f"C: limit exponent max_x [x log2 + H(x) + 4H(x/4)] = {best[0]:.4f} at x = {best[1]:.4f}")
    big = sieve(1 << 24)
    for y in (6, 8):
        lo, hi = math.exp(y), math.exp(2*y)
        if hi > (1 << 24): break
        qs = np.nonzero(big[int(lo)+1:int(hi)+1])[0] + int(lo) + 1
        print(f"C: y={y}: sum_(e^y<q<=e^2y) 1/q = {np.sum(1.0/qs):.4f} (log 2 = {math.log(2):.4f}); 'log2/y' = {math.log(2)/y:.4f}")

# ---------------------------------------------------------------- assignment p'(q) (Lemma 3.1)
def assignment(qmax):
    pp = {3: 2, 5: 2, 7: 3}
    k = 3
    while (1 << k) <= qmax:
        qs = [q for q in PR if (1 << k) <= q < (1 << (k+1)) and q > 7]
        ts = [t for t in PR if (1 << (k-2)) <= t < (1 << (k-1))]
        n, t = len(qs), len(ts)
        for i, q in enumerate(qs):
            pp[q] = ts[(i * t) // n]
        k += 1
    return pp

PP = assignment(LIM - 1)

def omega_table(n):
    om = np.zeros(n + 1, dtype=np.int64)
    for p in PR:
        if p > n: break
        om[p::p] += 1
    return om

def lambda_prime_table(n):
    lam = np.zeros(n + 1)
    for q, p in PP.items():
        lam[p::p] += math.log(q)
    return lam

def check_MU():
    n = 60000
    om = omega_table(n); lam = lambda_prime_table(n)
    d = np.arange(n + 1); logd = np.log(np.maximum(d, 1))
    # Lemma 3.1(3)
    bad = np.nonzero(lam[2:] > 9*logd[2:] + 18.72*om[2:] + 1e-9)[0]
    print(f"MU: Lambda'(d) <= 9 log d + 18.72 omega(d) for 2<=d<={n}: {'ok' if len(bad)==0 else 'FAIL at '+str(bad[:5]+2)}")
    perd = np.full(n + 1, np.nan)
    m = d >= 3
    perd[m] = (10 + 26/np.log(logd[m])) * logd[m]
    print("MU:      M   max_(d<M)(L'(d)+log d)/logM   per-d chain max/logM   max-omega chain/logM   mu_M")
    for M in (16, 44, 100, 1000, 10000, 60000):
        L = math.log(M); mu = 10 + 26/math.log(L)
        true = np.max(lam[1:M] + logd[1:M]) / L
        pdc = np.nanmax(perd[3:M]) / L
        omc = (10*L + 18.72*om[1:M].max()) / L
        print(f"MU: {M:6d}   {true:8.3f}   {pdc:10.3f}   {omc:8.3f}   {mu:7.3f}")
    # max omega(d), d < M, against Robin form 1.3841 log M/loglog M
    worst = min((1.3841*math.log(M)/math.log(math.log(M)) - om[1:M].max(), M) for M in range(16, n + 1, 7))
    print(f"MU: min over M in [16,{n}] of 1.3841 logM/loglogM - max_(d<M) omega(d) = {worst[0]:.3f} (at M={worst[1]})")

# ---------------------------------------------------------------- D
def demand_max(M, A):
    X = 3 * M ** A
    L = math.log(M)
    qs = [q for q in PR if 2 < q < 4 * M and q <= X]
    lam = np.zeros(M)
    d = np.arange(M)
    for q in qs:
        p = PP[q]
        Aq = int(math.floor(math.log(X) / math.log(q) + 1e-12))
        astar = 1
        while p * q ** (astar - 1) < M: astar += 1
        cap = min(astar, Aq)
        sel = d[p::p]
        # level a counts if q^(a-1) | d
        lev = np.ones(len(sel), dtype=np.int64)
        for a in range(2, cap + 1):
            lev += (sel % q ** (a - 1) == 0)
        lam[p::p] += lev * math.log(q)
    j = int(np.argmax(lam[1:]) + 1)
    return lam[j] / L, j

def check_D():
    for A in (1, 2):
        vals = [(demand_max(M, A), M) for M in range(20, 2001)]
        top = max(vals)
        print(f"D: A={A}: max over 20<=M<=2000 of max_d demand/log M = {top[0][0]:.3f} at M={top[1]} (d={top[0][1]})")
        for M in (23, 44, 64, 111, 256, 1024):
            v, dd = demand_max(M, A)
            print(f"D:   A={A} M={M}: {v:.3f} log M (d={dd})")
        big = [(demand_max(M, A)[0], M) for M in (4096, 16384)]
        print(f"D:   A={A}: M=4096,16384: {[f'{v:.3f}' for v, _ in big]}")

# ---------------------------------------------------------------- K
def oddsmooth(n, X):
    n = abs(n)
    while n and n % 2 == 0: n //= 2
    out = 1; q = 3; rem = n
    while q * q <= rem and q <= X:
        while rem % q == 0: rem //= q; out *= q
        q += 2
    if 1 < rem <= X: out *= rem
    return out

def check_K(X=10 ** 6, Mpts=20):
    ns = []
    n = int(math.isqrt(X)) + 1
    while len(ns) < Mpts:
        if n % 2 == 0:
            v = n * n + 1
            if all(v % p for p in PR if p * p <= v):
                ns.append(n)
        n += 1
    m1 = m_other = 0.0
    for x, y in itertools.combinations(ns, 2):
        P = x * y + 1; Q = y - x       # Omega = zeta_x conj(zeta_y) = P + iQ
        m1 = max(m1, math.log(oddsmooth(Q, X)))
        m_other = max(m_other, max(math.log(oddsmooth(t, X)) for t in (P, P - Q, P + Q)))
    Y = max(nn * nn + 1 for nn in ns)
    print(f"K: X={X}, M={Mpts}, zeta_x = n_x + i with n_x^2+1 prime in ({X}, {Y}]: max log m(pair, unit 1) = {m1:.2f}; "
          f"max log m(pair, unit != 1) = {m_other:.2f}; Lemma 4.1 pair bound 4 log Y = {4*math.log(Y):.2f}")

if __name__ == "__main__":
    check_G(); check_C(); check_MU(); check_D(); check_K(); check_K(X=10 ** 8)
