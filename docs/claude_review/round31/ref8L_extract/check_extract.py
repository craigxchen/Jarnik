#!/usr/bin/env python3
"""Referee checks for Astra's VIII.1 (extraction lemma).  Own code; exact integers/Fractions
where the claim is an identity, exact integer pattern counts for the sampled tuples."""
import itertools, random, math, sys
from fractions import Fraction as Fr
import numpy as np

rng = random.Random(20261119)
SUBS7 = list(range(1, 128))           # nonempty A' subset [7] <-> even A subset {0..7}

# ---------------------------------------------------------------- Part A: exact identities
def partA(trials=2):
    out = []
    for t in range(trials):
        M, n = 4, 6                                  # tiny so everything is exact
        mass = [Fr(rng.randint(1, 9)) for _ in range(n)]
        Wt = sum(mass); p = [m / Wt for m in mass]
        B = [[rng.choice((-1, 1)) for _ in range(n)] for _ in range(M)]
        # Gram and D
        G = [[sum(p[u] * B[x][u] * B[y][u] for u in range(n)) for y in range(M)] for x in range(M)]
        D = [[sum(mass[u] for u in range(n) if B[x][u] != B[y][u]) for y in range(M)] for x in range(M)]
        assert all(G[x][y] == 1 - 2 * D[x][y] / Wt for x in range(M) for y in range(M))
        EG = sum(G[x][y] for x in range(M) for y in range(M)) / M**2
        assert EG == sum(p[u] * (Fr(sum(B[x][u] for x in range(M)), M))**2 for u in range(n))
        assert EG >= 0
        # E_x m_A^2 = E_{u,v} g^{|A|}, all 8-tuples in X^8 (repeats allowed), |A| in {2,4,6,8}
        g = [[Fr(sum(B[x][u] * B[x][v] for x in range(M)), M) for v in range(n)] for u in range(n)]
        EGG = sum(G[x][y]**2 for x in range(M) for y in range(M)) / M**2
        for k in (2, 4, 6, 8):
            lhs = Fr(0)
            for xs in itertools.product(range(M), repeat=k):
                mA = sum(p[u] * math.prod(B[x][u] for x in xs) for u in range(n))
                lhs += mA * mA
            lhs /= M**k
            rhs = sum(p[u] * p[v] * g[u][v]**k for u in range(n) for v in range(n))
            assert lhs == rhs, (k, lhs, rhs)
            assert lhs <= EGG
        # pattern-mass expansion on all 8-tuples of a smaller sample
        for _ in range(15):
            xs = [rng.randrange(M) for _ in range(8)]
            pat = [Fr(0)] * 128
            for u in range(n):
                T = sum(1 << (i - 1) for i in range(1, 8) if B[xs[i]][u] != B[xs[0]][u])
                pat[T] += p[u]
            for T in range(128):
                s = Fr(1)
                for Ap in SUBS7:
                    A = [0] + [i for i in range(1, 8) if Ap >> (i - 1) & 1] if bin(Ap).count('1') % 2 else \
                        [i for i in range(1, 8) if Ap >> (i - 1) & 1]
                    mA = sum(p[u] * math.prod(B[xs[j]][u] for j in A) for u in range(n))
                    sign = (-1) ** bin(Ap & T).count('1')
                    s += sign * mA
                assert pat[T] == s / 128
        out.append("trial %d ok: E G=%s, E G^2=%s" % (t, EG, EGG))
    return out

# ---------------------------------------------------------------- Part B: sign systems
def walsh(v):
    v = v.astype(np.float64).copy(); h = 1
    while h < len(v):
        for i in range(0, len(v), 2 * h):
            a = v[i:i + h].copy(); b = v[i + h:i + 2 * h].copy()
            v[i:i + h] = a + b; v[i + h:i + 2 * h] = a - b
        h *= 2
    return v

def analyse(name, B, mass, a_list=(0.05, 0.1, 0.2), samples=3000):
    """B: M x n int8 matrix of +-1, mass: length-n positive integers."""
    M, n = B.shape
    Wt = int(mass.sum())
    Bf = B.astype(np.int64)
    Gnum = (Bf * mass) @ Bf.T                       # exact integers, G = Gnum / Wt
    off = Gnum[~np.eye(M, dtype=bool)]
    maxG = off.max() / Wt
    lam = max(maxG, 0) / 4
    EG2 = float((Gnum.astype(np.float64) ** 2).sum()) / Wt**2 / M**2
    EG = Gnum.sum() / Wt / M**2
    bound = 127 * math.sqrt(2 / M + 8 * lam)
    sums, devs, Dviol = [], [], 0
    for _ in range(samples):
        xs = [rng.randrange(M) for _ in range(8)]
        R = B[xs]
        diff = (R[1:] != R[0]).astype(np.int64)        # 7 x n
        code = (diff * (1 << np.arange(7))[:, None]).sum(0)
        pat = np.bincount(code, weights=None, minlength=128) if False else \
            np.bincount(code, weights=mass, minlength=128).astype(np.int64)   # exact integer masses
        assert pat.sum() == Wt
        dev = np.abs(pat * 128 - Wt).max() / Wt          # max_T |mass_T/w - 1|
        mhat = walsh(pat) / Wt                           # mhat[A'] = m_A
        s = np.abs(mhat[1:]).sum()
        assert dev <= s + 1e-12                          # deviation <= sum |m_A|
        sums.append(s); devs.append(dev)
        # 64-pattern statement: D_ij <= (1+dev) W/2 for every slot pair
        for i in range(8):
            for j in range(i + 1, 8):
                Dij = int(mass[R[i] != R[j]].sum())
                if Dij * 2 > (1 + dev) * Wt + 1e-9: Dviol += 1
    sums = np.array(sums); devs = np.array(devs)
    res = dict(name=name, M=M, n=n, maxG=round(maxG, 5), lam=round(lam, 5), EG=round(EG, 6),
               EG2=round(EG2, 6), bound_ES=round(bound, 4), mean_sum=round(sums.mean(), 4),
               mean_dev=round(devs.mean(), 4), Dviol=Dviol)
    for a in a_list:
        res['P(sum<=%g)' % a] = round(float((sums <= a).mean()), 4)
        res['Markov lower(%g)' % a] = round(max(0.0, 1 - min(bound, 127 * math.sqrt(EG2)) / a), 4)
    return res

def hadamard(k):
    n = 2 ** k
    idx = np.arange(n)
    pc = np.vectorize(lambda z: bin(z).count('1') & 1)
    H = 1 - 2 * pc(idx[:, None] & idx[None, :])
    return H.astype(np.int8)

def paley(p):
    sq = np.zeros(p, dtype=np.int8) - 1
    for t in range(p): sq[t * t % p] = 1          # chi(0) := +1
    u = np.arange(p)
    return sq[(u[None, :] + u[:, None]) % p]

def partB():
    out = []
    H = hadamard(11)                                 # M = n = 2048
    out.append(analyse('Sylvester Hadamard 2^11', H, np.ones(H.shape[1], dtype=np.int64)))
    P = paley(2003)
    out.append(analyse('Paley chi(u+x), p=2003', P, np.ones(2003, dtype=np.int64)))
    nrng = np.random.default_rng(7)
    Rnd = (1 - 2 * nrng.integers(0, 2, size=(2000, 8192))).astype(np.int8)
    out.append(analyse('uniform random M=2000,n=8192', Rnd, np.ones(8192, dtype=np.int64)))
    wmass = nrng.integers(1, 50, size=8192).astype(np.int64)
    out.append(analyse('uniform random, random masses', Rnd, wmass))
    # negative controls: hypothesis G <= 4 lambda fails with lambda large
    mu = 0.3
    Bias = np.where(nrng.random((2000, 8192)) < (1 + mu) / 2, 1, -1).astype(np.int8)
    out.append(analyse('biased mu=0.3 (G~0.09)', Bias, np.ones(8192, dtype=np.int64)))
    K = 6
    centres = (1 - 2 * nrng.integers(0, 2, size=(K, 8192))).astype(np.int8)
    lab = nrng.integers(0, K, size=2000)
    flip = np.where(nrng.random((2000, 8192)) < 0.02, -1, 1).astype(np.int8)
    Clu = (centres[lab] * flip).astype(np.int8)
    out.append(analyse('6 clusters, 2% noise', Clu, np.ones(8192, dtype=np.int64)))
    return out

# ---------------------------------------------------------------- Part C: Schwartz-Zippel / Z escape
def partC():
    out = []
    # Z = harmonic locus cr(x0,x1;x2,x3) = -1 union cr(x4,x5;x6,x7)=2 ; PGL2-invariant, proper.
    def F(x):
        f1 = (x[0]-x[2])*(x[1]-x[3]) + (x[0]-x[3])*(x[1]-x[2])
        f2 = (x[4]-x[6])*(x[5]-x[7]) - 2*(x[4]-x[7])*(x[5]-x[6])
        return f1 * f2                               # total degree 4
    for X in (list(range(1, 41)), sorted(rng.sample(range(-10**6, 10**6), 40)),
              [Fr(1, k) for k in range(1, 41)], [k * k for k in range(40)]):
        M = len(X)
        # exact count for f1 on X^4 (degree 2); full F on X^8 equals union, prob <= P1+P2
        c1 = sum(1 for x in itertools.product(X, repeat=4)
                 if (x[0]-x[2])*(x[1]-x[3]) + (x[0]-x[3])*(x[1]-x[2]) == 0)
        c2 = sum(1 for x in itertools.product(X, repeat=4)
                 if (x[0]-x[2])*(x[1]-x[3]) - 2*(x[0]-x[3])*(x[1]-x[2]) == 0)
        P = 1 - (1 - c1 / M**4) * (1 - c2 / M**4)    # exact, factors in disjoint variables
        out.append("X=%s..., M=%d: P(F=0)=%.5f <= d/M=%.5f : %s"
                   % (str(X[:3]), M, P, 4 / M, P <= 4 / M))
    return out

if __name__ == '__main__':
    for s in partA(): print('A', s)
    for s in partC(): print('C', s)
    for r in partB(): print('B', r)
