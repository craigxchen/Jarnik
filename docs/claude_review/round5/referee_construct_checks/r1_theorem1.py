"""Referee check R1: Theorem 1 of construct.md (character margins of skew-paired fully flipped profiles).

Independent implementation (does not import construct's lib.py).  All arithmetic is integer.

Profile P(H, pi, b, a*): H normalised Hadamard (column 0 all ones), labels 1..M-1; b unflipped copies of
every label; for every row x one flipped copy of label a(x) (entry at row x negated), a(x)=pi(x) for
x != x_inf, a(x_inf)=a*.  G(c) = sum_j (|(S^T c)_j| - 1).
Claims: (i) G >= b-4 for pairs; (ii) G >= 2n+b-8 for +-1 characters with n>=4; (iii) G >= (b-3)n/2+b
if ||c||_inf >= 2.  Also: pair minimum with a skew pairing is exactly b-2 (Remark 3.3).

Matrices tested:
  * Paley-I, q = 7, 11, 19, 23 (identity pairing after normalising), built from scratch;
  * a NON-Paley skew-Hadamard matrix of order 16, obtained by the doubling
    K = [[A, A], [A - 2I, -A + 2I]] from the Paley skew matrix A of order 8 (K + K^T = 2I, K K^T = 16 I);
    Sylvester-16 is the only order-16 class that cannot be skew (construct Remark 2.3), so this is a
    genuinely different core.
Exhaustive over all +-1 zero-sum characters up to the stated n; exhaustive over small-amplitude
||c||_inf >= 2 characters; adversarial hill-climbing at M = 44, 68 over integer c."""
import itertools, sys, time
import numpy as np

def legendre(a, q):
    a %= q
    return 0 if a == 0 else (1 if pow(a, (q - 1) // 2, q) == 1 else -1)

def paley_skew(q):
    """K = I + J, J skew, order q+1 (rows/cols: 0 = infinity, 1+i = i in F_q)."""
    n = q + 1
    K = np.zeros((n, n), dtype=np.int64)
    for i in range(n):
        for j in range(n):
            if i == j:
                K[i, j] = 1
            elif i == 0:
                K[i, j] = 1
            elif j == 0:
                K[i, j] = -1
            else:
                K[i, j] = legendre((j - 1) - (i - 1), q)
    assert (K + K.T == 2 * np.eye(n, dtype=np.int64)).all()
    assert (K @ K.T == n * np.eye(n, dtype=np.int64)).all()
    return K

def double_skew(A):
    n = A.shape[0]
    I = np.eye(n, dtype=np.int64)
    K = np.block([[A, A], [A - 2 * I, -A + 2 * I]])
    m = 2 * n
    assert (K + K.T == 2 * np.eye(m, dtype=np.int64)).all()
    assert (K @ K.T == m * np.eye(m, dtype=np.int64)).all()
    assert set(np.unique(K)) == {-1, 1}
    return K

def normalise_with_pairing(K):
    """Normalise column 0 to all ones by row signs; the skew pairing is x -> x (x != 0), x_inf = 0."""
    H = K * K[:, [0]]
    M = H.shape[0]
    assert (H[:, 0] == 1).all() and (H.T @ H == M * np.eye(M, dtype=np.int64)).all()
    # check skew-pairing minors for all distinct x,y != 0
    for x in range(1, M):
        for y in range(x + 1, M):
            assert H[x, x] * H[y, y] != H[x, y] * H[y, x]
    return H

def profile(H, b, astar):
    M = H.shape[0]
    cols = []
    for a in range(1, M):
        for _ in range(b):
            cols.append(H[:, a].copy())
    for x in range(M):
        a = astar if x == 0 else x
        col = H[:, a].copy(); col[x] = -col[x]
        cols.append(col)
    S = np.array(cols, dtype=np.int64).T
    assert np.linalg.matrix_rank(S.astype(float)) == M
    return S

def bound(n, kinf, b):
    if kinf >= 2:
        return (b - 3) * n / 2 + b
    if n == 2:
        return b - 4
    return 2 * n + b - 8

def pm1_chars(M, n):
    """all zero-sum +-1 characters with ||c||_1 = n, as an array (batch)."""
    h = n // 2
    out = []
    for supp in itertools.combinations(range(M), n):
        for plus in itertools.combinations(range(n), h):
            if 0 not in plus:   # c and -c give the same G; keep the one with +1 at the first support row
                continue
            c = np.zeros(M, dtype=np.int64)
            c[list(supp)] = -1
            c[[supp[i] for i in plus]] = 1
            out.append(c)
    return np.array(out, dtype=np.int64)

def Gbatch(S, C):
    Y = C @ S
    return np.abs(Y).sum(axis=1) - S.shape[1]

results = []
def run_exhaustive(name, H, b, astar, nmax):
    M = H.shape[0]
    S = profile(H, b, astar)
    worst = {}
    for n in range(2, nmax + 1, 2):
        # chunk by support to bound memory
        total = 0; mins = None
        for first in range(M):
            rows = []
            for rest in itertools.combinations(range(first + 1, M), n - 1):
                supp = (first,) + rest
                for plus in itertools.combinations(range(1, n), n // 2 - 1):
                    c = np.full(M, 0, dtype=np.int64)
                    c[list(supp)] = -1
                    c[first] = 1
                    c[[supp[i] for i in plus]] = 1
                    rows.append(c)
            if not rows:
                continue
            Cm = np.array(rows, dtype=np.int64)
            G = Gbatch(S, Cm)
            slack = G - bound(n, 1, b)
            k = int(np.argmin(slack))
            total += len(rows)
            if mins is None or slack[k] < mins:
                mins = int(slack[k])
        worst[n] = (total, mins)
        assert mins >= 0, (name, b, n, mins)
    return worst

def run_amp2(name, H, b, astar, nmax, K=2):
    """all zero-sum c with entries in [-K,K], ||c||_inf >= 2, ||c||_1 <= nmax."""
    M = H.shape[0]
    S = profile(H, b, astar)
    vals = [v for v in range(-K, K + 1) if v]
    cnt = 0; mins = None
    for k in range(2, nmax + 1):
        for supp in itertools.combinations(range(M), k):
            batch = []
            for vv in itertools.product(vals, repeat=k):
                if sum(vv) != 0 or max(abs(v) for v in vv) < 2 or sum(abs(v) for v in vv) > nmax:
                    continue
                c = np.zeros(M, dtype=np.int64); c[list(supp)] = vv
                batch.append(c)
            if batch:
                Cm = np.array(batch); G = Gbatch(S, Cm)
                n = np.abs(Cm).sum(axis=1)
                sl = G - ((b - 3) * n / 2 + b)
                cnt += len(batch)
                m = float(sl.min())
                mins = m if mins is None or m < mins else mins
    assert mins >= 0
    return cnt, mins

def hillclimb(H, b, astar, iters, rng):
    M = H.shape[0]
    S = profile(H, b, astar)
    def score(c):
        n = int(np.abs(c).sum()); kinf = int(np.abs(c).max())
        return int(np.abs(c @ S).sum()) - S.shape[1] - bound(n, kinf, b)
    best = None
    for it in range(iters):
        # random start: +-1 of random even support, or column/half-sum, or amplitude-2
        mode = it % 4
        c = np.zeros(M, dtype=np.int64)
        if mode == 0:
            h = 2 * int(rng.integers(2, M // 2 + 1))
            supp = rng.choice(M, size=h, replace=False)
            c[supp[: h // 2]] = 1; c[supp[h // 2:]] = -1
        elif mode == 1:
            a = int(rng.integers(1, M)); c = H[:, a].copy()
        elif mode == 2:
            a, bb = rng.choice(np.arange(1, M), size=2, replace=False)
            c = (H[:, a] + H[:, bb]) // 2
        else:
            supp = rng.choice(M, size=4, replace=False)
            c[supp[0]] = 2; c[supp[1]] = -1; c[supp[2]] = -1
        if not c.any() or c.sum() != 0:
            continue
        s = score(c)
        improved = True
        while improved:
            improved = False
            for _ in range(400):
                x, y = rng.choice(M, size=2, replace=False)
                d = int(rng.choice([1, -1]))
                c2 = c.copy(); c2[x] += d; c2[y] -= d
                if not c2.any() or np.abs(c2).max() > 3:
                    continue
                s2 = score(c2)
                if s2 < s:
                    c, s = c2, s2; improved = True
        if best is None or s < best[0]:
            best = (s, int(np.abs(c).sum()), int(np.abs(c).max()))
    assert best[0] >= 0
    return best

t0 = time.time()
print("R1: Theorem 1 (construct.md), independent implementation")
mats = []
for q in (7, 11, 19, 23):
    mats.append(("Paley-I q=%d (M=%d)" % (q, q + 1), normalise_with_pairing(paley_skew(q))))
K16 = double_skew(paley_skew(7))
mats.append(("doubled skew-Hadamard (M=16, non-Paley, non-Sylvester)", normalise_with_pairing(K16)))
plan = {8: 8, 12: 12, 16: 8, 20: 6, 24: 6}
for name, H in mats:
    M = H.shape[0]
    for b in (3, 5, 8):
        w = run_exhaustive(name, H, b, 1, plan[M])
        s = "  ".join("n=%d:#%d min slack %d" % (n, t, m) for n, (t, m) in w.items())
        print("%-55s b=%d exhaustive +-1:  %s" % (name, b, s))
    pair_min = w[2][1] + (8 - 4)   # b=8 last run: G_min on pairs
    print("%-55s     pair minimum G (b=8) = %d  (= b-2 claimed)" % (name, pair_min))
    assert pair_min == 6
for name, H in mats[:3] + mats[4:]:
    M = H.shape[0]
    nm = 8 if M <= 12 else 6
    cnt, mn = run_amp2(name, H, 5, 1, nm)
    print("%-55s b=5 amplitude>=2, entries in [-2,2], ||c||_1<=%d: %d chars, min slack %.1f" % (name, nm, cnt, mn))
rng = np.random.default_rng(12345)
for q in (43, 67):
    H = normalise_with_pairing(paley_skew(q))
    for b in (3, 8):
        best = hillclimb(H, b, 1, 120, rng)
        print("Paley-I q=%d (M=%d) b=%d: hill-climb (120 starts) min slack %d at n=%d, ||c||_inf=%d"
              % (q, q + 1, b, best[0], best[1], best[2]))
print("ALL R1 CHECKS PASSED  (%.0fs)" % (time.time() - t0))
