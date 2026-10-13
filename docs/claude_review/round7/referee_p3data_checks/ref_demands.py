"""Referee: independent recomputation of genuine unit-1 / any-unit collision moduli by direct
Gaussian arithmetic mod q^A (no discrete logs), from the stored primes, for the construct/lower profiles
rebuilt from construct.md section 2 (own code). Compares with the author's stored mod_*.npy."""
import json, sys, os, itertools
import numpy as np
from sympy import isprime, primerange
DATA = '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round7/p3data_checks/data'

def legendre(a, q):
    a %= q
    if a == 0: return 0
    return 1 if pow(a, (q-1)//2, q) == 1 else -1

def paley(q):
    n = q + 1
    H = np.zeros((n, n), dtype=np.int64)
    # rows/cols: 0 = inf, 1+i = F_q element i. construct.md: H(inf,.)=1, H(i,inf)=-1, H(i,j)=delta+chi(j-i); negate rows i in F_q
    H[0, :] = 1
    for i in range(q):
        H[1+i, 0] = 1  # -(-1)
        for j in range(q):
            H[1+i, 1+j] = -((1 if i == j else 0) + legendre(j - i, q))
    assert (H.T @ H == n*np.eye(n, dtype=np.int64)).all()
    return H

def prof_construct(q, b=8):
    H = paley(q); M = q + 1
    cols = [H[:, a] for a in range(1, M) for _ in range(b)]
    assign = [1] + list(range(1, M))
    for x in range(M):
        c = H[:, assign[x]].copy(); c[x] = -c[x]; cols.append(c)
    return np.array(cols).T

def gmul(a, b, Q):
    return ((a[0]*b[0] - a[1]*b[1]) % Q, (a[0]*b[1] + a[1]*b[0]) % Q)

def moduli(S, A, B, X):
    M, r = S.shape
    a = (1 + S)//2
    m1 = np.ones((M, M), dtype=object); many = np.ones((M, M), dtype=object)
    units = [(1, 0), (0, 1), (-1, 0), (0, -1)]
    for q in primerange(3, X+1):
        Amax = 1
        while q**(Amax+1) <= X: Amax += 1
        Q = q**Amax
        z = []
        for x in range(M):
            v = (1, 0)
            for j in range(r):
                g = (A[j], B[j]) if a[x, j] == 1 else (A[j], -B[j])
                v = gmul(v, g, Q)
            z.append(v)
        for x in range(M):
            for y in range(M):
                if x == y: continue
                # unit 1 level
                lev = 0
                for L in range(1, Amax+1):
                    QL = q**L
                    if (z[x][0]-z[y][0]) % QL == 0 and (z[x][1]-z[y][1]) % QL == 0: lev = L
                    else: break
                m1[x, y] *= q**lev
                # any unit: z_x = u z_y mod q^L for some unit u
                levany = 0
                for L in range(1, Amax+1):
                    QL = q**L; ok = False
                    for u in units:
                        w = gmul(u, z[y], QL)
                        if (z[x][0]-w[0]) % QL == 0 and (z[x][1]-w[1]) % QL == 0: ok = True
                    if ok: levany = L
                    else: break
                many[x, y] *= q**levany
    return m1, many

if __name__ == '__main__':
    bad = 0; n = 0
    for M in [8, 12, 20, 24, 32]:
        q0 = M - 1
        S = prof_construct(q0)
        assert S.shape[1] == 9*M - 8
        for t in range(12):
            f = f'{DATA}/primes_construct_M{M}_t{t}.json'
            if not os.path.exists(f): continue
            d = json.load(open(f)); ps, A, B = d['ps'], d['A'], d['B']
            assert len(ps) == S.shape[1] and len(set(ps)) == len(ps)
            for p, x, y in zip(ps, A, B):
                assert x*x + y*y == p and p % 4 == 1 and isprime(p)
            for k, X in (('X2', 2*M), ('X4', 4*M)):
                m1, many = moduli(S, A, B, X)
                s1 = np.load(f'{DATA}/mod_construct_M{M}_t{t}_u1_{k}.npy', allow_pickle=True)
                sa = np.load(f'{DATA}/mod_construct_M{M}_t{t}_any_{k}.npy', allow_pickle=True)
                ok1 = all(int(m1[i, j]) == int(s1[i, j]) for i in range(M) for j in range(M) if i != j)
                oka = all(int(many[i, j]) == int(sa[i, j]) for i in range(M) for j in range(M) if i != j)
                Dst = np.load(f'{DATA}/D_construct_M{M}_t{t}_u1_{k}.npy')
                okD = max(abs(Dst[i, j] - 2*np.log(float(m1[i, j]))) for i in range(M) for j in range(M) if i != j) < 1e-9
                n += 1
                if not (ok1 and oka and okD):
                    bad += 1; print('MISMATCH', M, t, k, ok1, oka, okD)
        print('M', M, 'done', flush=True)
    print('instances checked', n, 'mismatches', bad)
