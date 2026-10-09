"""Independent (round6) checks of round5/lower.md Lemma S and Prop. 3.3, which round6/sharp.md
re-derives and uses.  Evidence only; the proofs are in sharp.md §4.

For Paley q0 in {43, 47, 59} (M = q0 + 1), F(c) = sum over the q0 nonconstant labels of |(H^T c)_a|.
 (1) exhaustive: every ternary zero-sum c of support 4 has F >= (17/16) M   (no exceptions exist
     at support 4: pairs have support 2, columns and half-sums support >= M/2);
 (2) random ternary c of support 6..M, random amplitude-2 vectors, all signed columns and
     half-sums: non-exceptional ones have F >= (17/16) M, amplitude >= 2 has F >= 2M;
 (3) canonical one-flip profile with b = 33 (profile.py): 2B - r >= b - 4 (pairs, exhaustive),
     >= (b/16 - 2) M + b (all other tested characters; Prop. 3.3 of lower.md, from Lemma S),
     where 2B = sum_j |c^T S_j|.
Prints the minima found and ALL LEMMA S CHECKS PASSED."""
import itertools, random
import numpy as np
from profile import paley_H, build

rng = random.Random(11)
ok = True
for q0 in (43, 47, 59):
    H = paley_H(q0)                      # M x q0, nonconstant labels
    M = q0 + 1
    thr = 17 * M / 16
    # (1) exhaustive support 4
    mins4 = 10 ** 9
    subsets = np.array(list(itertools.combinations(range(M), 4)), dtype=np.int64)
    for signs in ((1, 1, -1, -1), (1, -1, 1, -1), (1, -1, -1, 1)):
        s = np.array(signs, dtype=np.int64)
        for start in range(0, len(subsets), 20000):
            Sb = subsets[start:start + 20000]
            T = (H[Sb] * s[None, :, None]).sum(axis=1)        # (batch, q0)
            Fv = np.abs(T).sum(axis=1)
            mins4 = min(mins4, int(Fv.min()))
    ok &= mins4 >= thr
    # (2) random larger supports, amplitude 2, columns, half-sums
    minrand, minamp = 10 ** 9, 10 ** 9
    for _ in range(20000):
        h = rng.choice(range(6, M + 1, 2))
        rows = rng.sample(range(M), h)
        c = np.zeros(M, dtype=np.int64)
        for k, x in enumerate(rows):
            c[x] = 1 if k < h // 2 else -1
        Fv = int(np.abs(c @ H).sum())
        T = c @ H
        is_col = (np.abs(T) == M).sum() == 1 and (np.abs(T) == 0).sum() == q0 - 1
        is_half = (np.abs(T) == M // 2).sum() == 2 and (np.abs(T) == 0).sum() == q0 - 2
        if not (is_col or is_half):
            minrand = min(minrand, Fv)
    for _ in range(5000):
        c = np.array([rng.choice([-2, -1, 0, 0, 1, 2]) for _ in range(M)], dtype=np.int64)
        c[rng.randrange(M)] -= c.sum()
        if np.abs(c).max() >= 2 and c.any():
            minamp = min(minamp, int(np.abs(c @ H).sum()))
    ok &= minrand >= thr and minamp >= 2 * M
    colF = [int(np.abs(H[:, a] @ H).sum()) for a in range(q0)]
    ok &= all(f == M for f in colF)
    # (3) profile with b = 33
    b = 33
    Hh, S, A, flip, sib, alpha = build(q0, b)
    r = S.shape[1]
    minpair = 10 ** 9
    for x in range(M):
        for y in range(x + 1, M):
            c = np.zeros(M, dtype=np.int64)
            c[x], c[y] = 1, -1
            minpair = min(minpair, int(np.abs(c @ S).sum()) - r)
    ok &= minpair >= b - 4
    minother = 10 ** 9
    tests = []
    for a in range(q0):
        tests.append(Hh[:, a].copy())
    for _ in range(300):
        a1, a2 = rng.sample(range(q0), 2)
        tests.append((Hh[:, a1] + rng.choice([1, -1]) * Hh[:, a2]) // 2)
    for x in subsets[rng.sample(range(len(subsets)), 3000)]:
        c = np.zeros(M, dtype=np.int64)
        c[x] = (1, 1, -1, -1)
        tests.append(c)
    for _ in range(3000):
        h = rng.choice(range(6, M + 1, 2))
        rows = rng.sample(range(M), h)
        c = np.zeros(M, dtype=np.int64)
        for k, x in enumerate(rows):
            c[x] = 1 if k < h // 2 else -1
        tests.append(c)
    for c in tests:
        if np.abs(c).sum() == 0:
            continue
        minother = min(minother, int(np.abs(c @ S).sum()) - r)
    ok &= minother >= (b / 16 - 2) * M + b
    print(f"q0 = {q0:3d} M = {M}: min F support 4 = {mins4} (17M/16 = {thr:.1f}); min F random non-exc = {minrand}; "
          f"min F amplitude>=2 = {minamp} (2M = {2*M}); b=33 profile: min 2B-r pairs = {minpair} (b-4 = {b-4}), "
          f"others = {minother} (proved bound (b/16-2)M+b = {(b/16-2)*M+b:.1f})")
print("ALL LEMMA S CHECKS PASSED" if ok else "FAILURE")
