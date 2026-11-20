#!/usr/bin/env python3
"""Referee check (block fakedefs): definitional properties of Astra's VIII.3 profile.

Own code; exact integer arithmetic (numpy int64 on small integer matrices, sympy for ranks mod q / Smith form).
For q = 3 mod 4 prime: K_ij = chi(j - i) on F_q, class A = rows with K_{i0} = +1, S subset of A (density rho),
a = 1 + K on I = F_q minus {0}, e = 2.  Checks:
  (a) K^T = -K, K K^T = qI - J; every pair of rows of A has l1 distance q+1 on I, coefficient sum
      sum_{j in I}(|K_ij - K_i'j| - 1) = 2, pair support (q+3)/2 on I;
  (b) rank_Q K_S = m; rank of K_S^T (on I) mod q equals m  => every rational zero-sum lambda with
      K_S^T lambda integral is integral (Astra's forced-character lemma); Smith form of the row differences
      on I (small q): all invariants 1 (saturated character lattice);
  (c) free dimension m_free = r - rank L, r = q - 1, rank L = m - 1;
  (d) normalisation: columns of I on which the values over S miss -1 or +1;
  (e) s = 4 characters lambda = e_a + e_b - e_c - e_d: min ||K_S^T lambda||_1 and min support on I,
      against Astra's 1.2 q and (q-1)/2.
"""
import itertools, random, sys
import numpy as np
from sympy import Matrix, isprime
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ


def chi(x, q):
    x %= q
    if x == 0:
        return 0
    return 1 if pow(x, (q - 1) // 2, q) == 1 else -1


def rank_mod_p(M, p):
    A = [list(map(lambda t: t % p, row)) for row in M]
    rows, cols = len(A), len(A[0])
    r = 0
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if A[i][c] % p:
                piv = i
                break
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        inv = pow(A[r][c], p - 2, p)
        A[r] = [(x * inv) % p for x in A[r]]
        for i in range(rows):
            if i != r and A[i][c]:
                f = A[i][c]
                A[i] = [(x - f * y) % p for x, y in zip(A[i], A[r])]
        r += 1
        if r == rows:
            break
    return r


def run(q, rho, seed, do_snf, do_s4):
    assert isprime(q) and q % 4 == 3
    K = np.array([[chi(j - i, q) for j in range(q)] for i in range(q)], dtype=np.int64)
    J = np.ones((q, q), dtype=np.int64)
    ok_a = bool((K.T == -K).all() and (K @ K.T == q * np.eye(q, dtype=np.int64) - J).all())
    A = [i for i in range(q) if K[i, 0] == 1]
    rnd = random.Random(seed)
    S = [i for i in A if rnd.random() < rho]
    m = len(S)
    I = list(range(1, q))
    KS = K[np.ix_(S, I)]
    # (a) pairs
    l1set, csum, supp = set(), set(), set()
    for x, y in itertools.combinations(range(m), 2):
        d = np.abs(KS[x] - KS[y])
        l1set.add(int(d.sum()))
        csum.add(int((d - 1).sum()))
        supp.add(int((d != 0).sum()))
    # (b) ranks
    rQ = rank_mod_p(KS.tolist(), 1000000007)  # rank mod a large prime <= rank over Q <= m
    rq = rank_mod_p(KS.T.tolist(), q)
    snf_ok = None
    if do_snf:
        Bd = Matrix((KS[1:] - KS[0]).tolist())
        D = smith_normal_form(Bd, domain=ZZ)
        diag = [D[i, i] for i in range(min(D.shape)) if D[i, i] != 0]
        snf_ok = (len(diag), sorted(set(abs(int(t)) for t in diag)))
    # (d) normalisation
    vals = KS  # values of K over S, columns I
    miss = int(sum(1 for j in range(q - 1) if not ((vals[:, j] == 1).any() and (vals[:, j] == -1).any())))
    # (e) s = 4 characters
    s4 = None
    if do_s4:
        h_min, sup_min = 10**9, 10**9
        idx = range(m)
        for a, b in itertools.combinations(idx, 2):
            pos = KS[a] + KS[b]
            for c, d in itertools.combinations([t for t in idx if t not in (a, b)], 2):
                v = pos - KS[c] - KS[d]
                h = int(np.abs(v).sum())
                sp = int((v != 0).sum())
                if h < h_min:
                    h_min = h
                if sp < sup_min:
                    sup_min = sp
        s4 = (h_min, round(h_min / q, 3), sup_min, (q - 1) // 2)
    r = q - 1
    m_free = r - (m - 1)
    print(f"q={q} |A|={len(A)} m=|S|={m} (a) KKt ok={ok_a} pair l1 {sorted(l1set)} coefsum {sorted(csum)} "
          f"pair supp {sorted(supp)} (=(q+3)/2={(q+3)//2}) | (b) rankQ={rQ} rank mod q={rq} SNF={snf_ok} | "
          f"(c) r={r} m_free={m_free} m_free/m={m_free/m:.3f} | (d) columns missing an endpoint={miss} | (e) s=4: {s4}")
    sys.stdout.flush()


if __name__ == "__main__":
    for q, snf, s4 in [(43, True, True), (67, False, True), (103, False, True), (127, False, True),
                       (163, False, False), (383, False, False)]:
        run(q, 0.99, 12345 + q, snf, s4)
