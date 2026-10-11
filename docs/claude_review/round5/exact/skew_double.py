"""Skew double core (Theorem 4.1 of exact.md).

K skew-Hadamard of order M (K + K^T = 2I, K K^T = M I), normalised so that K(x,0) = 1 for all x
(D K D with D = diag(K(x,0) K(0,0))).  Cuts:  S_a = {x : K(x,a) = -1},  S'_a = S_a xor {0, a}
(a = 1..M-1).  This script
  (1) builds the Paley skew-Hadamard matrix of order q+1 (q = 3 mod 4 prime), or a Kronecker
      product skew matrix of order 2^t, and checks skewness and the Hadamard property exactly;
  (2) checks that the 2(M-1) cuts are distinct and balanced, that every pair of rows is separated
      by exactly M cuts (Plotkin-tight), and computes G_T (exact diagonal form);
  (3) checks the key inequality  Phi(c) = ||K^T c||_1 + ||K c~||_1 >= 2(M-1) + ||c||_1  (c~ = c with
      c_0 negated) exhaustively on all ternary c with support <= smax and on hill-climbed vectors
      with entries in [-2, 2]  (integer arithmetic);
  (4) optionally runs the exact verifier (inner.analyse, rho1 = 1, unit weights) on every class.
usage: python3 skew_double.py q|s<t> smax nclimb exact(0/1)
"""
import sys, itertools, random, time
import numpy as np
from fractions import Fraction
from core import pattern_of, rows_of, Group


def normalise(K):
    M = K.shape[0]
    d = K[:, 0] * K[0, 0]
    K = (d[:, None] * K) * d[None, :]
    assert np.all(K[:, 0] == 1)
    assert np.all(K + K.T == 2 * np.eye(M, dtype=np.int64)), 'not skew'
    assert np.all(K @ K.T == M * np.eye(M, dtype=np.int64)), 'not Hadamard'
    return K


def paley_skew(q):
    chi = lambda a: 0 if a % q == 0 else (1 if pow(a % q, (q - 1) // 2, q) == 1 else -1)
    M = q + 1
    S = np.zeros((M, M), dtype=np.int64)
    for i in range(M):
        for j in range(M):
            if i == j:
                continue
            if i == 0:
                S[i, j] = -1
            elif j == 0:
                S[i, j] = 1
            else:
                S[i, j] = chi((j - 1) - (i - 1))
    return normalise(np.eye(M, dtype=np.int64) + S)


def kron_skew(t):
    """Skew-Hadamard of order 2^t:  K_1 = [[1,1],[-1,1]],  K_{t+1} = [[K_t, K_t],[-K_t^T... ]] via the
    standard doubling  K' = [[K, K], [K - 2I, -K + 2I]]  (K skew => K' skew Hadamard)."""
    K = np.array([[1, 1], [-1, 1]], dtype=np.int64)
    for _ in range(t - 1):
        n = K.shape[0]
        I = np.eye(n, dtype=np.int64)
        K = np.block([[K, K], [K - 2 * I, -K + 2 * I]])
    return normalise(K)


def cuts(K):
    M = K.shape[0]
    S = [frozenset(x for x in range(M) if K[x, a] == -1) for a in range(1, M)]
    S2 = [s ^ {0, a} for a, s in zip(range(1, M), S)]
    return S, S2


def Phi(K, c):
    ct = c.copy(); ct[0] = -ct[0]
    return int(np.abs(K.T @ c).sum() + np.abs(K @ ct).sum())


def main(spec, smax, nclimb, exact):
    K = kron_skew(int(spec[1:])) if spec.startswith('s') else paley_skew(int(spec))
    M = K.shape[0]
    S, S2 = cuts(K)
    T = [pattern_of(M, sorted(s)) for s in S + S2]
    ok = len(set(T)) == 2 * (M - 1) and all(len(s) == M // 2 for s in S + S2)
    sep = set()
    for x, y in itertools.combinations(range(M), 2):
        sep.add(sum(1 for s in S + S2 if (x in s) != (y in s)))
    G = Group(M, T)
    print(f'[{spec}] M={M}: 2(M-1) distinct balanced cuts: {ok}; pair separation values {sorted(sep)} '
          f'(Plotkin-tight iff [{M}]); G_T invariants {G.invs}', flush=True)
    # sigma identity check on random integer c:  2 sigma_c = Phi(c) - 2(M-1)
    P = np.array([[1 if x in s else 0 for x in range(M)] for s in S + S2], dtype=np.int64)
    rng = random.Random(1)
    for _ in range(200):
        c = np.array([rng.randint(-3, 3) for _ in range(M)], dtype=np.int64)
        c[0] -= c.sum()
        two_sigma = 2 * int(np.abs(P @ c).sum()) - 2 * (M - 1)
        assert two_sigma == Phi(K, c) - 2 * (M - 1), 'identity 2 sigma = Phi - 2(M-1) fails'
    print('  identity 2 sigma_c = Phi(c) - 2(M-1): checked on 200 random integer c', flush=True)
    worst = None
    cnt = 0
    for h in range(2, smax + 1, 2):
        for supp in itertools.combinations(range(M), h):
            for pos in itertools.combinations(supp, h // 2):
                if supp[0] not in pos:
                    continue
                c = np.zeros(M, dtype=np.int64)
                for x in supp:
                    c[x] = 1 if x in pos else -1
                val = Phi(K, c) - (2 * (M - 1) + h)
                cnt += 1
                if worst is None or val < worst[0]:
                    worst = (val, h, supp)
    print(f'  exhaustive ternary, support <= {smax}: {cnt} vectors (one of each +-pair), '
          f'min Phi - (2M-2+n) = {worst[0]} (at n = {worst[1]})', flush=True)
    best = None
    for t in range(nclimb):
        c = np.zeros(M, dtype=np.int64)
        k = rng.randrange(2, M + 1)
        for x in rng.sample(range(M), k):
            c[x] = rng.choice((-1, 1))
        c[rng.randrange(M)] -= c.sum()
        if not c.any():
            continue
        f = lambda v: Phi(K, v) - (2 * (M - 1) + int(np.abs(v).sum()))
        cur = f(c)
        for step in range(600):
            x, y = rng.sample(range(M), 2)
            d = rng.choice((-1, 1))
            c2 = c.copy(); c2[x] += d; c2[y] -= d
            if not c2.any() or np.abs(c2).max() > 2:
                continue
            v = f(c2)
            if v <= cur:
                c, cur = c2, v
        if best is None or cur < best[0]:
            best = (cur, c.copy())
    if best:
        print(f'  hill climbing ({nclimb} starts, entries in [-2,2]): min Phi - (2M-2+n) = {best[0]} '
              f'(at n = {int(np.abs(best[1]).sum())})', flush=True)
    if exact:
        from inner import analyse, find_hom
        t0 = time.time()
        A = analyse(M, T, [Fraction(1)] * len(T), 0, 1)
        eta = find_hom(G, list(A['bad'].keys()))
        print(f'  exact verifier (rho1 = 1, unit weights, all classes of G_T): integral violations '
              f'{len(A["int_viol"])}, bad rational classes {len(A["bad"])}, avoiding homomorphism {eta}, '
              f'min slack over all classes {A["min_slack"][0]} ({time.time()-t0:.0f}s)', flush=True)


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
