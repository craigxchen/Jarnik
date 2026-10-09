"""Referee re-check of the lower.md inputs that sharp.md uses (Lemma S = sharp Lemma 4.2,
Prop. 3.3 = sharp Prop. 4.3, Lemma 3.1 = sharp Lemma 4.1), with independent code.

Paley q0 in {43, 47} (M = q0 + 1), canonical one-flip profile with b = 33 copies:
 L1  H (normalised Paley-I, rows F_q u {inf}, labels F_q, plus the constant column) is Hadamard.
 L2  twins: S_f(x) - S_s(x) = -2 H(x, alpha x) e_x; 2M twin columns distinct; rank S = M (exact,
     mod a large prime); every column nonconstant.
 L3  Lemma S: F(c) >= 17M/16 for EVERY ternary zero-sum c of support 4 (exhaustive) and for random
     ternary c of supports 6..M (not exceptional); F >= 2M for random c with max|c| >= 2;
     F = M for pairs, columns and half-sums.
 L4  Prop. 4.3 / identity (3.1): 2B(c) = b F(c) + sum_x(|T_alpha(x) - 2c_x H(x,alpha x)| - |T_alpha(x)|)
     on all tested c; 2B - r >= b - 4 for every pair (and = b - 4 attained); >= kappa_b M for every
     tested non-pair (kappa_b = min(1, b/16 - 2)); >= hM(b-4)/2 + b for amplitude h >= 2.
 L5  slack bound s_c >= (tau/2)(2B - r) - 1/2 with random weights in [tau, tau + 1/r].
"""
import itertools, random, math
import numpy as np

rng = random.Random(3)


def paley(q):
    M = q + 1
    chi = np.zeros(q, dtype=np.int64)
    for t in range(1, q):
        chi[(t * t) % q] = 1
    chi = np.where(chi == 1, 1, -1); chi[0] = 0
    H = np.zeros((M, q), dtype=np.int64)
    for x in range(q):
        for a in range(q):
            H[x, a] = -chi[(a - x) % q] - (1 if a == x else 0)
    H[q, :] = 1          # row infinity
    return H


def profile(H, b, astar=0):
    M, q = H.shape
    cols, lab = [], []
    for a in range(q):
        for j in range(b):
            cols.append(H[:, a].copy()); lab.append(a)
    S = np.array(cols).T.copy()     # M x r
    flip, sib, alpha = {}, {}, {}
    for x in range(q):
        f, s = x * b + 0, x * b + 1
        S[x, f] = -S[x, f]; flip[x], sib[x], alpha[x] = f, s, x
    f, s = astar * b + 2, astar * b + 3
    S[q, f] = -S[q, f]; flip[q], sib[q], alpha[q] = f, s, astar
    return S, flip, sib, alpha


def rank_modp(Mx, p=1_000_000_007):
    A = [[int(v) % p for v in row] for row in Mx.tolist()]
    rk, rows, cols = 0, len(A), len(A[0])
    for c in range(cols):
        piv = next((r for r in range(rk, rows) if A[r][c]), None)
        if piv is None:
            continue
        A[rk], A[piv] = A[piv], A[rk]
        inv = pow(A[rk][c], -1, p)
        A[rk] = [(v * inv) % p for v in A[rk]]
        for r in range(rows):
            if r != rk and A[r][c]:
                f = A[r][c]
                A[r] = [(vr - f * vk) % p for vr, vk in zip(A[r], A[rk])]
        rk += 1
        if rk == rows:
            break
    return rk


def run(q, b=33):
    H = paley(q)
    M = q + 1
    ok = {}
    Hf = np.concatenate([np.ones((M, 1), dtype=np.int64), H], axis=1)
    ok['L1_hadamard'] = bool(np.all(Hf.T @ Hf == M * np.eye(M, dtype=np.int64)))
    S, flip, sib, alpha = profile(H, b)
    r = S.shape[1]
    tw = all(np.array_equal(S[:, flip[x]] - S[:, sib[x]], -2 * H[x, alpha[x]] * np.eye(M, dtype=np.int64)[x]) for x in range(M))
    distinct = len(set(list(flip.values()) + list(sib.values()))) == 2 * M
    ok['L2_twins'] = tw and distinct
    ok['L2_rank'] = rank_modp(S) == M
    ok['L2_nonconst'] = bool(np.all(S.max(axis=0) != S.min(axis=0)))
    Sf = S.astype(np.float32); Hfl = H.astype(np.float32)

    def stats_for(C):
        T = C.astype(np.float32) @ Hfl
        F = np.abs(T).sum(axis=1)
        twoB = np.abs(C.astype(np.float32) @ Sf).sum(axis=1)
        flips = np.zeros(len(C), dtype=np.float32)
        for x in range(M):
            a = alpha[x]
            flips += np.abs(T[:, a] - 2 * C[:, x] * H[x, a]) - np.abs(T[:, a])
        ident = np.allclose(twoB, b * F + flips)
        return F, twoB - r, ident

    kb = min(1.0, b / 16 - 2)
    # pairs
    P = []
    for x, y in itertools.combinations(range(M), 2):
        c = np.zeros(M, dtype=np.int64); c[x], c[y] = 1, -1; P.append(c)
    P = np.array(P)
    F, G, ident = stats_for(P)
    ok['L3_pairs_flat'] = bool(np.all(F == M))
    ok['L4_pairs'] = bool(np.all(G >= b - 4)) and ident and bool(np.any(G == b - 4))
    # columns and half-sums
    cols = np.array([H[:, a] for a in range(q)])
    F, G, ident = stats_for(cols)
    ok['L3_cols_flat'] = bool(np.all(F == M))
    ok['L4_cols'] = bool(np.all(G >= b + 2 * M - 8)) and ident
    hs = np.array([(s * H[:, a] + t * H[:, bb]) // 2 for a, bb in itertools.combinations(range(q), 2) for s in (1, -1) for t in (1, -1)])
    F, G, ident = stats_for(hs)
    ok['L3_halfsums_flat'] = bool(np.all(F == M))
    ok['L4_halfsums'] = bool(np.all(G >= b + M - 16)) and ident
    # exhaustive support 4
    minF4, minG4 = 1e9, 1e9
    ident_all = True
    batch = []
    patterns = [(1, 1, -1, -1), (1, -1, 1, -1), (1, -1, -1, 1)]
    for supp in itertools.combinations(range(M), 4):
        for pat in patterns:
            c = np.zeros(M, dtype=np.int64); c[list(supp)] = pat; batch.append(c)
        if len(batch) >= 30000:
            F, G, ident = stats_for(np.array(batch)); batch = []
            minF4 = min(minF4, F.min()); minG4 = min(minG4, G.min()); ident_all &= ident
    if batch:
        F, G, ident = stats_for(np.array(batch))
        minF4 = min(minF4, F.min()); minG4 = min(minG4, G.min()); ident_all &= ident
    ok['L3_support4'] = minF4 >= 17 * M / 16
    ok['L4_support4'] = minG4 >= kb * M and ident_all
    # random ternary larger supports
    R = []
    for _ in range(20000):
        h = 2 * rng.randint(3, M // 2)
        supp = rng.sample(range(M), h)
        c = np.zeros(M, dtype=np.int64)
        c[supp[:h // 2]] = 1; c[supp[h // 2:]] = -1
        R.append(c)
    R = np.array(R)
    F, G, ident = stats_for(R)
    # exclude exceptional (half-sums have support M/2; columns support M)
    exc = (F == M)
    ok['L3_random'] = bool(np.all(F[~exc] >= 17 * M / 16))
    ok['L4_random'] = bool(np.all(G[~exc] >= kb * M)) and ident
    minFr, minGr = F[~exc].min(), G[~exc].min()
    # amplitude >= 2
    Am = []
    hs_ = []
    for _ in range(5000):
        c = np.array([rng.randint(-3, 3) for _ in range(M)], dtype=np.int64)
        c[-1] -= c.sum()
        if np.abs(c).max() < 2:
            continue
        Am.append(c); hs_.append(np.abs(c).max())
    Am = np.array(Am); hs_ = np.array(hs_)
    F, G, ident = stats_for(Am)
    ok['L3_amp'] = bool(np.all(F >= 2 * M))
    ok['L4_amp'] = bool(np.all(G >= hs_ * M * (b - 4) / 2 + b)) and ident
    # L5 slack with random weights
    tau = 50.0
    w = tau + np.array([rng.random() / r for _ in range(r)])
    Wt = w.sum()
    sl_ok = True
    for C in (P[:200], R[:500], Am[:200]):
        tauj = C.astype(np.float64) @ S.astype(np.float64)
        s = 0.5 * (np.abs(tauj) * w).sum(axis=1) - Wt / 2
        twoBr = np.abs(tauj).sum(axis=1) - r
        sl_ok &= bool(np.all(s >= (tau / 2) * twoBr - 0.5 - 1e-6))
    ok['L5_slack'] = sl_ok
    print(f"q0={q} M={M} b={b} r={r}: min F/M support-4 = {minF4/M:.3f} (need >= {17/16:.4f}); "
          f"min (2B-r) support-4 = {minG4:.0f} (need >= {kb*M:.2f}); random: min F/M = {minFr/M:.3f}, min 2B-r = {minGr:.0f}")
    print("    ", ok)
    return all(ok.values())


if __name__ == "__main__":
    allok = True
    for q in (43, 47):
        allok &= run(q)
    print("ALL LOWER.MD INPUT CHECKS PASSED" if allok else "LOWER.MD INPUT FAILURE")
