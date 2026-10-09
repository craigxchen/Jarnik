"""Exact checks for the fully flipped Paley profile of Theorem 1 (equal-weight form) and for the
residue design of Lemma 4.

Profile: Paley H of order M = q+1, b columns per nonconstant label a (all oriented as H_a);
row x is flipped at one column f(x) of label alpha(x), alpha at most 2-to-1, each flipped column
has an unflipped sibling in its label.  S is the M x r sign matrix (r = b(M-1)).
For an integral zero-sum c put  G(c) = sum_j (|c^T S_j| - 1)  (equal-weight slack: s_c = (tau/2) G(c)).

Checked exactly:
 (a) identity  G(c) = b * sum_a (|T_a| - 1) + sum_x (|T_alpha(x) - 2 c_x H(x,alpha(x))| - |T_alpha(x)|);
 (b) pairs: G >= b - 4;  columns: G >= b + 2M - 8;  half-sums: G >= b + M - 16;
 (c) every other tested c: G >= b * max(F - M + 1, n - M + 1) - 2n   (n = ||c||_1);
 (d) saturation: S_f(x) - S_sib(x) = -2 H(x,alpha(x)) e_x for every row (twins), rank S = M;
 (e) residue design: rho_x = beta * w_x / conj(w_x) mod q^a with Gaussian primes w_x above distinct
     split primes; pair collision at q^a  <=>  q^a | Im(w_x conj w_y); for random characters c the
     generous collision modulus m_c (product of q^a with prod rho^c == unit mod q^a) divides
     Im(Omega^4)/4, Omega = prod_{c>0} w^c prod_{c<0} conj(w)^|c|.
"""
import random
from itertools import combinations
import numpy as np
from paley import paley, transform, F, exceptional, is_prime

random.seed(5)
np.random.seed(5)
FAIL = []


def check(cond, msg):
    if not cond:
        FAIL.append(msg)
        print("FAIL:", msg)


def build_profile(q, b):
    H = paley(q)
    M = q + 1
    labels = list(range(1, M))
    # alpha: rows -> labels, at most 2-to-1, random
    slots = labels + labels
    random.shuffle(slots)
    alpha = slots[:M]
    cols = []          # (label, flipped_row or None)
    flipcol = {}
    sib = {}
    for a in labels:
        rows_a = [x for x in range(M) if alpha[x] == a]
        start = len(cols)
        for t in range(b):
            cols.append((a, rows_a[t] if t < len(rows_a) else None))
        for t, x in enumerate(rows_a):
            flipcol[x] = start + t
            sib[x] = start + len(rows_a) + t   # an unflipped column of the same label
    r = len(cols)
    S = np.zeros((M, r), dtype=np.int64)
    for j, (a, x) in enumerate(cols):
        S[:, j] = H[:, a]
        if x is not None:
            S[x, j] = -S[x, j]
    return H, S, alpha, flipcol, sib


def G(S, c):
    return int((np.abs(c @ S) - 1).sum())


def run_profile(q, b):
    H, S, alpha, flipcol, sib = build_profile(q, b)
    M = q + 1
    # (d) twins and rank
    for x in range(M):
        d = S[:, flipcol[x]] - S[:, sib[x]]
        e = np.zeros(M, dtype=np.int64); e[x] = -2 * H[x, alpha[x]]
        check((d == e).all(), "twin q=%d x=%d" % (q, x))
    check(np.linalg.matrix_rank(S.astype(float)) == M, "rank")

    def ident(c):
        T = transform(H, c)
        rhs = b * int((np.abs(T) - 1).sum())
        for x in range(M):
            Ta = int(T[alpha[x] - 1])
            rhs += abs(Ta - 2 * int(c[x]) * int(H[x, alpha[x]])) - abs(Ta)
        return rhs

    mins = {"pair": 10 ** 9, "column": 10 ** 9, "halfsum": 10 ** 9, "other": 10 ** 9}
    # pairs
    for x, y in combinations(range(M), 2):
        c = np.zeros(M, dtype=np.int64); c[x] = 1; c[y] = -1
        g = G(S, c)
        check(g == ident(c), "identity pair")
        check(g >= b - 4, "pair bound")
        mins["pair"] = min(mins["pair"], g)
    # columns and half-sums
    for a in range(1, M):
        c = H[:, a].copy()
        g = G(S, c)
        check(g == ident(c) and g >= b + 2 * M - 8, "column bound")
        mins["column"] = min(mins["column"], g)
    for a, bb in combinations(range(1, M), 2):
        for s in (1, -1):
            c = (H[:, a] + s * H[:, bb]) // 2
            g = G(S, c)
            check(g == ident(c) and g >= b + M - 16, "halfsum bound")
            mins["halfsum"] = min(mins["halfsum"], g)
    # others: random structured + exhaustive support 4 (+-1)
    from check_lemmaS import structured
    tested = 0
    worst_ratio = 10 ** 9
    vecs = list(structured(H, 3000))
    for supp in random.sample(list(combinations(range(M), 4)), 3000):
        c = np.zeros(M, dtype=np.int64); c[list(supp)] = [1, 1, -1, -1]
        vecs.append(c)
    for c in vecs:
        if exceptional(H, c):
            continue
        n = int(np.abs(c).sum())
        g = G(S, c)
        Fc = F(H, c)
        check(g == ident(c), "identity other")
        bound = b * max(Fc - M + 1, n - M + 1) - 2 * n
        check(g >= bound, "other bound")
        mins["other"] = min(mins["other"], g)
        worst_ratio = min(worst_ratio, g / (b * M))
        tested += 1
    print("q=%3d M=%3d b=%3d: min G  pairs %d (>= b-4=%d), columns %d (>= %d), half-sums %d (>= %d), "
          "others %d over %d vectors (min G/(bM) = %.3f)"
          % (q, M, b, mins["pair"], b - 4, mins["column"], b + 2 * M - 8, mins["halfsum"],
             b + M - 16, mins["other"], tested, worst_ratio))


# ------------------------------------------------------------------ residue design (e)
def gmul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def gconj(a):
    return (a[0], -a[1])


def gmod(a, m):
    return (a[0] % m, a[1] % m)


def ginv_mod(a, m):
    # inverse of a mod m (a coprime to m): a^{-1} = conj(a) / Norm(a)
    nrm = (a[0] * a[0] + a[1] * a[1]) % m
    t = pow(nrm, -1, m)
    c = gconj(a)
    return gmod((c[0] * t, c[1] * t), m)


def gpow_mod(a, e, m):
    if e < 0:
        a = ginv_mod(a, m)
        e = -e
    r = (1, 0)
    while e:
        if e & 1:
            r = gmod(gmul(r, a), m)
        a = gmod(gmul(a, a), m)
        e >>= 1
    return r


def gaussian_prime_above(p):
    for u in range(1, int(p ** 0.5) + 2):
        v2 = p - u * u
        if v2 <= 0:
            break
        v = int(round(v2 ** 0.5))
        if v * v == v2:
            return (u, v)
    raise ValueError


def residue_design_check(M=24, X=60):
    # w_x Gaussian primes above the first M split primes > X
    ells = []
    p = X + 1
    while len(ells) < M:
        if is_prime(p) and p % 4 == 1:
            ells.append(p)
        p += 1
    w = [gaussian_prime_above(l) for l in ells]
    qs = [q for q in range(3, X + 1) if is_prime(q)]
    units = [(1, 0), (0, 1), (-1, 0), (0, -1)]
    Nfake = 1
    for l in ells:
        pass
    tested_pairs = 0
    tested_chars = 0
    for q in qs:
        A = 1
        while q ** (A + 1) <= X:
            A += 1
        mod = q ** A
        # beta: any element of norm == N mod q^A; norm is surjective, take N = 5^? -> choose beta=(2,1)
        beta = (2, 1)
        rho = [gmod(gmul(beta, gmul(wx, ginv_mod(gconj(wx), mod))), mod) for wx in w]
        for x, y in combinations(range(M), 2):
            im = gmul(w[x], gconj(w[y]))[1]
            for a in range(1, A + 1):
                m = q ** a
                coll = gmod(rho[x], m) == gmod(rho[y], m)
                check(coll == (im % m == 0), "pair collision criterion q^a=%d" % m)
                tested_pairs += 1
    # generous collision modulus for random characters
    for _ in range(400):
        s = random.randint(2, 6)
        supp = random.sample(range(M), s)
        vals = [random.choice([1, -1, 2]) for _ in range(s - 1)]
        vals.append(-sum(vals))
        if vals[-1] == 0:
            continue
        Om = (1, 0)
        for x, v in zip(supp, vals):
            f = w[x] if v > 0 else gconj(w[x])
            for _ in range(abs(v)):
                Om = gmul(Om, f)
        O2 = gmul(Om, Om)
        O4 = gmul(O2, O2)
        target = abs(O4[1]) // 4
        check(O4[1] % 4 == 0 and target != 0, "Im Omega^4 divisible by 4, nonzero")
        mgen = 1
        for q in qs:
            A = 1
            while q ** (A + 1) <= X:
                A += 1
            best = 0
            for a in range(1, A + 1):
                mq = q ** a
                prod = (1, 0)
                for x, v in zip(supp, vals):
                    prod = gmod(gmul(prod, gpow_mod(gmul(w[x], ginv_mod(gconj(w[x]), mq)), v, mq)), mq)
                if any(gmod(u, mq) == prod for u in units):
                    best = a
            mgen *= q ** best
        check(target % mgen == 0, "m_gen divides Im(Omega^4)/4")
        tested_chars += 1
    print("residue design: %d pair/level collision tests, %d character tests (X=%d, M=%d)"
          % (tested_pairs, tested_chars, X, M))


if __name__ == "__main__":
    for q, b in ((43, 8), (47, 8), (59, 6), (67, 12), (83, 8)):
        run_profile(q, b)
    residue_design_check()
    print("ALL CHECKS PASSED" if not FAIL else "FAILURES: %d" % len(FAIL))
