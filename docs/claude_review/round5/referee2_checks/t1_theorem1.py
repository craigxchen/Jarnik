"""Referee-2 check T1: Theorem 1 and Lemma 2.2 of construct.md, independent code (no import of construct's lib).

Theorem 1: H normalised Hadamard (column 0 all ones) with a skew pairing (x_inf, pi); profile = b unflipped
copies of every label 1..M-1 plus, for every row x, one copy of label a(x) negated in row x
(a(x) = pi(x), a(x_inf) = a*).  G(c) = sum_j (|(S^T c)_j| - 1).  Claims, for nonzero zero-sum integer c,
n = ||c||_1:
   (i)  G >= b-4 (pairs);  (ii) G >= 2n+b-8 (+-1 valued, n >= 4);  (iii) G >= (b-3)n/2+b (||c||_inf >= 2).
Lemma 2.2: G = (b+1)(F-M) + b + |T_a*| + sum_x Delta_x.

What is tested (all integer arithmetic):
 A. Paley-I built directly from the construct's definition H(inf,.)=1, H(i,inf)=-1, H(i,j)=delta_ij+chi(j-i),
    rows i negated.  ALL skew pairings (x_inf, pi) of Paley-8 are enumerated (8 * 7! candidates); for every one of
    them, every a*, b in {3,8}: exhaustive over all +-1 zero-sum c, and over amplitude-2 c with ||c||_1 <= 8.
 B. Paley-12: 40 random skew pairings found by randomised backtracking (most are not the identity); exhaustive
    over all +-1 zero-sum c (73,788) for b = 3 and 8; amplitude-2 c with ||c||_1 <= 6.
 C. Necessity of the skew hypothesis: Sylvester-16 (which has no skew pairing) with several bijections; we
    exhibit +-1 characters with G < 2n+b-8 (so (ii) genuinely fails without the hypothesis).
 D. Lemma 2.2 identity on 3,000 random integer c for Paley-20/24 with random skew-or-not assignments.
"""
import itertools, random, sys, time
import numpy as np

t0 = time.time()
rng = random.Random(20261009)

def chi(a, q):
    a %= q
    return 0 if a == 0 else (1 if pow(a, (q - 1) // 2, q) == 1 else -1)

def paley_I(q):
    n = q + 1
    H = np.zeros((n, n), dtype=np.int64)       # index 0 = infinity, index 1+i = element i of F_q
    for J in range(n):
        H[0, J] = 1
    for I in range(1, n):
        H[I, 0] = -1
        for J in range(1, n):
            i, j = I - 1, J - 1
            H[I, J] = (1 if i == j else 0) + chi(j - i, q)
    assert (H @ H.T == n * np.eye(n, dtype=np.int64)).all()
    Hn = H * H[:, [0]]                            # negate rows i in F_q
    assert (Hn[:, 0] == 1).all()
    return Hn

def sylvester(t):
    H = np.array([[1]], dtype=np.int64)
    for _ in range(t):
        H = np.block([[H, H], [H, -H]])
    return H

def is_skew_pairing(H, xinf, pi):
    rows = [x for x in range(H.shape[0]) if x != xinf]
    for x, y in itertools.combinations(rows, 2):
        if H[x, pi[x]] * H[y, pi[y]] == H[x, pi[y]] * H[y, pi[x]]:
            return False
    return True

def profile(H, b, assign):
    M = H.shape[0]
    cols = []
    for a in range(1, M):
        for _ in range(b):
            cols.append(H[:, a].copy())
    for x in range(M):
        col = H[:, assign[x]].copy(); col[x] = -col[x]; cols.append(col)
    return np.array(cols, dtype=np.int64).T

def pm1_chars(M, nmax):
    out = []
    for n in range(2, nmax + 1, 2):
        half = n // 2
        sgn = []
        for pos in itertools.combinations(range(n), half):
            s = -np.ones(n, dtype=np.int8); s[list(pos)] = 1; sgn.append(s)
        sgn = np.array(sgn)
        for supp in itertools.combinations(range(M), n):
            c = np.zeros((len(sgn), M), dtype=np.int8); c[:, supp] = sgn
            out.append(c)
    return np.concatenate(out)

def amp2_chars(M, nmax):
    out = []
    vals = [-2, -1, 1, 2]
    for k in range(2, nmax + 1):
        for supp in itertools.combinations(range(M), k):
            for vv in itertools.product(vals, repeat=k):
                if sum(vv) == 0 and sum(map(abs, vv)) <= nmax and max(map(abs, vv)) == 2:
                    c = np.zeros(M, dtype=np.int8); c[list(supp)] = vv; out.append(c)
    return np.array(out)

def bound(n, amp, b):
    if amp >= 2:
        return (b - 3) * n / 2 + b
    return b - 4 if n == 2 else 2 * n + b - 8

def check(S, C, b):
    """return (min slack, worst c) over the character matrix C"""
    Y = C.astype(np.int64) @ S
    G = np.abs(Y).sum(axis=1) - S.shape[1]
    n = np.abs(C).sum(axis=1).astype(np.int64)
    amp = np.abs(C).max(axis=1)
    bd = np.where(amp >= 2, (b - 3) * n / 2 + b, np.where(n == 2, b - 4, 2 * n + b - 8))
    sl = G - bd
    i = int(np.argmin(sl))
    return float(sl[i]), C[i], int(G[i])

ALL_OK = True
# ---------------- A: all skew pairings of Paley-8 ----------------
H8 = paley_I(7); M = 8
pairings = []
for xinf in range(M):
    rows = [x for x in range(M) if x != xinf]
    for perm in itertools.permutations(range(1, M)):
        pi = dict(zip(rows, perm))
        if is_skew_pairing(H8, xinf, pi):
            pairings.append((xinf, pi))
print("A. Paley-8: %d skew pairings among %d candidates (identity pairing included: %s)"
      % (len(pairings), M * 5040, any(xi == 0 and all(p[x] == x for x in p) for xi, p in pairings)))
C8 = pm1_chars(8, 8); A8 = amp2_chars(8, 8)
worstA = None; pairmin = None
for xinf, pi in pairings:
    for astar in range(1, M):
        assign = [astar if x == xinf else pi[x] for x in range(M)]
        for b in (3, 8):
            S = profile(H8, b, assign)
            for CC in (C8, A8):
                sl, cw, G = check(S, CC, b)
                if worstA is None or sl < worstA[0]:
                    worstA = (sl, list(cw), b)
                if sl < 0:
                    ALL_OK = False
            # pair minimum
            P = np.array([np.eye(M, dtype=np.int8)[x] - np.eye(M, dtype=np.int8)[y] for x in range(M) for y in range(M) if x != y])
            Gp = (np.abs(P.astype(np.int64) @ S).sum(axis=1) - S.shape[1]).min()
            pm = int(Gp) - (b - 2)
            pairmin = pm if pairmin is None else min(pairmin, pm)
print("   all pairings x all a* x b in {3,8}: %d +-1 chars + %d amplitude-2 chars each; min slack over bound = %s (at c=%s, b=%d)"
      % (len(C8), len(A8), worstA[0], worstA[1], worstA[2]))
print("   pair minimum G - (b-2) over all cases = %d  (Remark 3.3 / previous referee: exactly b-2)" % pairmin)

# ---------------- B: random skew pairings of Paley-12 ----------------
H12 = paley_I(11); M = 12
def random_skew_pairing(H, tries=2000):
    M = H.shape[0]
    for _ in range(tries):
        xinf = rng.randrange(M)
        rows = [x for x in range(M) if x != xinf]; rng.shuffle(rows)
        labels = list(range(1, M))
        pi = {}
        def bt(i):
            if i == len(rows):
                return True
            x = rows[i]
            cand = [a for a in labels if a not in pi.values()]
            rng.shuffle(cand)
            for a in cand:
                if all(H[x, a] * H[y, pi[y]] != H[x, pi[y]] * H[y, a] for y in rows[:i]):
                    pi[x] = a
                    if bt(i + 1):
                        return True
                    del pi[x]
            return False
        if bt(0):
            return xinf, dict(pi)
    return None
C12 = pm1_chars(12, 12); A12 = amp2_chars(12, 6)
found = []; nonid = 0
while len(found) < 40:
    res = random_skew_pairing(H12)
    if res is None:
        break
    xinf, pi = res
    assert is_skew_pairing(H12, xinf, pi)
    found.append(res)
    if not (xinf == 0 and all(pi[x] == x for x in pi)):
        nonid += 1
worstB = None
for xinf, pi in found:
    astar = rng.randrange(1, M)
    assign = [astar if x == xinf else pi[x] for x in range(M)]
    for b in (3, 8):
        S = profile(H12, b, assign)
        for CC in (C12, A12):
            sl, cw, G = check(S, CC, b)
            if worstB is None or sl < worstB[0]:
                worstB = (sl, list(cw), b)
            if sl < 0:
                ALL_OK = False
print("B. Paley-12: %d random skew pairings (%d non-identity); exhaustive %d +-1 chars and %d amplitude-2 chars each, b in {3,8}: min slack %s (c=%s, b=%d)"
      % (len(found), nonid, len(C12), len(A12), worstB[0], worstB[1], worstB[2]))

# ---------------- C: Sylvester-16, no skew pairing: (ii) fails ----------------
H16 = sylvester(4); M = 16
C16 = pm1_chars(16, 8)
viol = []
for trial in range(6):
    if trial == 0:
        assign = [1] + list(range(1, 16))          # x -> x (x != 0), x_inf = 0 -> label 1
    else:
        perm = list(range(1, 16)); rng.shuffle(perm)
        assign = [rng.randrange(1, 16)] + perm
    for b in (8,):
        S = profile(H16, b, assign)
        sl, cw, G = check(S, C16, b)
        viol.append((sl, int(np.abs(cw).sum()), G))
print("C. Sylvester-16 (no skew pairing), b=8, all +-1 chars with n<=8: min slack per assignment =", [v[0] for v in viol],
      " (negative = Theorem 1 (ii) bound violated; e.g. n=%d, G=%d)" % (min(viol)[1], min(viol)[2]))

# ---------------- D: Lemma 2.2 identity ----------------
cnt = 0
for q in (19, 23):
    H = paley_I(q); M = q + 1
    for _ in range(1500):
        assign = [rng.randrange(1, M) for _ in range(M)]
        b = rng.choice([3, 5, 8])
        S = profile(H, b, assign)
        c = np.array([rng.randint(-3, 3) for _ in range(M)], dtype=np.int64)
        c[-1] -= c.sum()
        if not c.any():
            continue
        T = H.T @ c; assert T[0] == 0
        F = int(np.abs(T[1:]).sum())
        Dsum = sum(abs(int(T[assign[x]]) - 2 * int(c[x]) * int(H[x, assign[x]])) - abs(int(T[assign[x]])) for x in range(M))
        # |T_a*| term: a appears twice in assign iff it's the x_inf label; general form: sum_x |T_a(x)| - F
        extra = sum(abs(int(T[assign[x]])) for x in range(M)) - F
        G = int(np.abs(c @ S).sum()) - S.shape[1]
        rhs = (b + 1) * (F - M) + b + extra + Dsum + (F - sum(abs(int(T[a])) for a in set(assign)) if False else 0)
        # For a general assignment hitting labels with multiplicities m_a: G = b(F-M+1) + sum_x(|T_a(x)|-1) + sum Delta
        rhs2 = b * (F - M + 1) + sum(abs(int(T[assign[x]])) - 1 for x in range(M)) + Dsum
        assert G == rhs2, (G, rhs2)
        if len(set(assign)) == M - 1:   # a bijection onto labels plus one repeat: Lemma 2.2's exact form
            astar_vals = [a for a in set(assign) if assign.count(a) == 2]
            if len(astar_vals) == 1:
                assert G == (b + 1) * (F - M) + b + abs(int(T[astar_vals[0]])) + Dsum
        cnt += 1
print("D. Lemma 2.2 (general-assignment form G = b(F-M+1) + sum_x(|T_a(x)|-1) + sum Delta) on %d random integer c: OK" % cnt)

print("ALL T1 CHECKS PASSED" if ALL_OK else "T1: A FAILURE OCCURRED", " (%.0fs)" % (time.time() - t0))
