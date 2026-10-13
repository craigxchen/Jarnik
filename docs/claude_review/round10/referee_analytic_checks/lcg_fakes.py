"""Referee check: test the claimed fake-violated input LCG (analytic.md Def. 5.1, Thm D) on explicit
realisations of (i) the round-5 construct.md fake (fully flipped Paley-I profile, b copies) and
(ii) the skew double core of exact.md / lw.md Thm 6.1 (one prime per cut, e = 1, windings k = e_1).

For each sample:
  * primes p_j = 1 mod 4 chosen as in the source theorem (window, size), deterministic Miller-Rabin;
  * delta uniform in [0, Delta]^M (random rationals), realisation phi = particular + random kernel element,
    computed with mpmath at high precision; the realisation (1.1) is re-verified;
  * EVIDENCE that the sample is a fake: Lemma 0.1 gaps for every v = e_j, every e_j +- e_l, every
    character c (zero-sum, entries in [-2,2], ||c||_1 <= 4), and 20000 random monomials;
  * LCG test: for each j, min over alpha in Z[i], 0 < Norm alpha <= p_j, of |Im(alpha e^{-i phi_j})| * sqrt(p_j).
    LCG_{p_j}(phi_j; p_j) fails iff this ratio is < 1 (and > 0).  For the genuine angle arg(pi_j) the same
    ratio is |Im(alpha conj pi_j)| >= 1 for every non-parallel alpha (exact integers);
  * fibre integrality (Prop. A remark): with e = 1, the fibre sum N^{2k} S_k(N)/4 = prod_j 2 Re(pi_j^{4k}) is an
    integer for actual clusters; we print its distance to Z at the fake angles (k = 1).
usage: python3 lcg_fakes.py
"""
import random, math, sys
import numpy as np
from mpmath import mp, mpf, matrix, lu_solve, cos, sin, pi as mpi, asin, exp, log, sqrt, nint, fabs, atan2

def is_prime(n):
    if n < 2: return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0: return n == p
    d = n - 1; s = 0
    while d % 2 == 0: d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True

def legendre(a, q):
    a %= q
    if a == 0: return 0
    return 1 if pow(a, (q - 1) // 2, q) == 1 else -1

def paley1(q):
    Q = np.array([[legendre(j - i, q) for j in range(q)] for i in range(q)], dtype=np.int64)
    n = q + 1; S = np.zeros((n, n), dtype=np.int64)
    S[0, 1:] = 1; S[1:, 0] = -1; S[1:, 1:] = Q
    H = np.eye(n, dtype=np.int64) + S
    H = H * H[:, [0]]
    assert (H.T @ H == n * np.eye(n, dtype=np.int64)).all()
    return H

def construct_profile(q, b):
    H = paley1(q); M = q + 1
    assign = [1] + list(range(1, q + 1))
    cols = []
    for a in range(1, M):
        for _ in range(b): cols.append(H[:, a].copy())
    for x in range(M):
        col = H[:, assign[x]].copy(); col[x] = -col[x]; cols.append(col)
    return np.array(cols, dtype=np.int64).T          # M x r, entries +-1 = 2a - e

def skew_double_profile(q):
    chi = lambda a: 0 if a % q == 0 else legendre(a, q)
    M = q + 1; S = np.zeros((M, M), dtype=np.int64)
    for i in range(M):
        for j in range(M):
            if i == j: continue
            if i == 0: S[i, j] = -1
            elif j == 0: S[i, j] = 1
            else: S[i, j] = chi((j - 1) - (i - 1))
    K = np.eye(M, dtype=np.int64) + S
    d = K[:, 0] * K[0, 0]; K = (d[:, None] * K) * d[None, :]
    assert (K + K.T == 2 * np.eye(M, dtype=np.int64)).all() and (K @ K.T == M * np.eye(M, dtype=np.int64)).all()
    cuts = [frozenset(x for x in range(M) if K[x, a] == -1) for a in range(1, M)]
    cuts2 = [c ^ {0, a} for a, c in zip(range(1, M), cuts)]
    cols = [[1 if x in c else -1 for x in range(M)] for c in cuts + cuts2]
    return np.array(cols, dtype=np.int64).T

def primes_window(P, r, logwidth):
    """r primes = 1 mod 4 in [P', P' e^logwidth] with P' >= P (first window found)."""
    out = []; n = P + (1 - P) % 4
    while True:
        if is_prime(n):
            out.append(n)
            while math.log(out[-1] / out[0]) > logwidth: out.pop(0)
            if len(out) == r: return out
        n += 4

def two_squares(p):
    # Cornacchia / brute: a^2 + b^2 = p, a > b > 0
    for b in range(1, math.isqrt(p // 2) + 1):
        a2 = p - b * b; a = math.isqrt(a2)
        if a * a == a2: return max(a, b), min(a, b)

def realise(Smat, ps, C0, kvec, rng, dps):
    mp.dps = dps
    M, r = Smat.shape
    W = sum(log(mpf(p)) for p in ps)
    Delta = C0 * exp(-W / 4)
    delta = [Delta * mpf(rng.randrange(1, 10**60)) / 10**60 for _ in range(M)]
    # unknowns (phi_1..phi_r, theta0):  sum_j S_xj phi_j - theta0 = delta_x + (pi/2) k_x
    A = np.hstack([Smat.astype(float), -np.ones((M, 1))])
    # basis of kernel (float, then exact rational-ish combination is fine: any real kernel vector is allowed)
    _, sv, Vt = np.linalg.svd(A)
    rank = int((sv > 1e-9).sum())
    kern = Vt[rank:]
    # particular solution: pick M independent columns
    cols = []
    for j in range(r + 1):
        if np.linalg.matrix_rank(A[:, cols + [j]]) > len(cols): cols.append(j)
        if len(cols) == rank: break
    assert rank == M, rank
    B = matrix([[int(A[x, j]) for j in cols] for x in range(M)])
    rhs = matrix([delta[x] + mpi / 2 * kvec[x] for x in range(M)])
    sol = lu_solve(B, rhs)
    vec = [mpf(0)] * (r + 1)
    for i, j in enumerate(cols): vec[j] = sol[i]
    # random kernel element: Haar-generic angles.  Kernel vectors are computed in float; to keep the
    # realisation exact we project them onto the exact kernel by re-solving for the pivot coordinates.
    free = [j for j in range(r + 1) if j not in cols]
    for j in free:
        vec[j] = mpf(rng.randrange(-10**60, 10**60)) / 10**58
    # re-solve pivots with free coordinates fixed
    rhs2 = matrix([delta[x] + mpi / 2 * kvec[x] - sum(int(A[x, j]) * vec[j] for j in free) for x in range(M)])
    sol2 = lu_solve(B, rhs2)
    for i, j in enumerate(cols): vec[j] = sol2[i]
    phi = [vec[j] % (2 * mpi) for j in range(r)]
    theta0 = vec[r]
    # verify (1.1) modulo 2 pi:  sum S_xj phi_j - theta0 - delta_x - (pi/2)k_x in 2 pi Z
    worst = max(fabs(((sum(int(Smat[x, j]) * phi[j] for j in range(r)) - theta0 - delta[x] - mpi / 2 * kvec[x]) / (2 * mpi))
                     - nint((sum(int(Smat[x, j]) * phi[j] for j in range(r)) - theta0 - delta[x] - mpi / 2 * kvec[x]) / (2 * mpi)))
                for x in range(M))
    return phi, delta, W, Delta, worst

def gap_ok(v, phi, logp):
    ang = sum(int(v[j]) * phi[j] for j in range(len(v)) if v[j])
    w = sum(abs(int(v[j])) * logp[j] for j in range(len(v)) if v[j])
    t = ang % (mpi / 2)
    d_ax = min(t, mpi / 2 - t)
    t2 = (ang - mpi / 4) % (mpi / 2)
    d_dg = min(t2, mpi / 2 - t2)
    e = exp(-w / 2)
    return min(d_ax / asin(e), d_dg / asin(e / sqrt(2)))

def lcg_ratio(phi, n):
    """min over alpha, 0 < Norm alpha <= n, of |Im(alpha e^{-i phi})| * sqrt(n), and the minimiser."""
    c, s = float(cos(phi)), float(sin(phi))
    R = math.isqrt(n)
    best = (float('inf'), None)
    if abs(c) >= abs(s):
        xs = np.arange(1, R + 1, dtype=np.float64)
        ys = np.rint(xs * s / c)
        val = np.abs(ys * c - xs * s)
        cand = np.argsort(val)[:20]
        pts = [(int(xs[i]), int(ys[i])) for i in cand]
    else:
        ys = np.arange(1, R + 1, dtype=np.float64)
        xs = np.rint(ys * c / s)
        val = np.abs(ys * c - xs * s)
        cand = np.argsort(val)[:20]
        pts = [(int(xs[i]), int(ys[i])) for i in cand]
    for (x, y) in pts:
        if x * x + y * y > n or (x == 0 and y == 0): continue
        im = fabs(y * cos(phi) - x * sin(phi))
        if im == 0: continue
        rr = im * sqrt(n)
        if rr < best[0]: best = (rr, (x, y))
    return best

def run(label, Smat, ps, C0, kvec, seed, dps, nrand=20000):
    rng = random.Random(seed)
    M, r = Smat.shape
    phi, delta, W, Delta, worst = realise(Smat, ps, C0, kvec, rng, dps)
    dps_gap = int(float(W) / (4 * 2.302585)) + 60
    mp.dps = dps_gap
    logp = [log(mpf(p)) for p in ps]
    print("=== %s: M=%d r=%d primes %d..%d, W=%.2f, W/(M log M)=%.2f, realisation residual %.1e" %
          (label, M, r, ps[0], ps[-1], float(W), float(W) / (M * math.log(M)), float(worst)))
    # evidence: gaps
    g1 = min(gap_ok([1 if i == j else 0 for i in range(r)], phi, logp) for j in range(r))
    g2 = mpf('inf')
    for j in range(r):
        for l in range(j + 1, r):
            for sgn in (1, -1):
                v = [0] * r; v[j] = 1; v[l] = sgn
                g2 = min(g2, gap_ok(v, phi, logp))
    # characters: zero-sum c with entries in [-2,2], ||c||_1 <= 4 ; v(c) = S^T c / 2
    gc = mpf('inf'); nchar = 0
    import itertools
    for supp in itertools.combinations(range(M), 2):
        for a in (1, 2):
            c = np.zeros(M, dtype=np.int64); c[supp[0]] = a; c[supp[1]] = -a
            for sg in (1, -1):
                v = (Smat.T @ (sg * c)) // 2
                if np.any(v): gc = min(gc, gap_ok(list(v), phi, logp)); nchar += 1
    for supp in itertools.combinations(range(M), 4):
        for signs in itertools.product((1, -1), repeat=4):
            if sum(signs) != 0: continue
            c = np.zeros(M, dtype=np.int64)
            for x, sg in zip(supp, signs): c[x] = sg
            v = (Smat.T @ c) // 2
            if np.any(v): gc = min(gc, gap_ok(list(v), phi, logp)); nchar += 1
    for supp in itertools.combinations(range(M), 3):
        for signs in ((2, -1, -1), (-2, 1, 1), (-1, 2, -1), (1, -2, 1), (-1, -1, 2), (1, 1, -2)):
            c = np.zeros(M, dtype=np.int64)
            for x, sg in zip(supp, signs): c[x] = sg
            v = (Smat.T @ c) // 2
            if np.any(v): gc = min(gc, gap_ok(list(v), phi, logp)); nchar += 1
    rng2 = random.Random(seed + 1); gr = mpf('inf')
    for _ in range(nrand):
        k = rng2.randint(1, 6); v = [0] * r
        for j in rng2.sample(range(r), k): v[j] = rng2.choice((-2, -1, 1, 2))
        gr = min(gr, gap_ok(v, phi, logp))
    print("  EVIDENCE fake gaps (min dist/required, >1 means satisfied): singles %.3g, pairs %.3g, "
          "characters(%d) %.3g, random monomials %.3g" % (float(g1), float(g2), nchar, float(gc), float(gr)))
    # LCG at fake angles
    rat = []
    for j in range(r):
        rr, al = lcg_ratio(phi[j], ps[j])
        rat.append((float(rr), al, ps[j]))
    mx = max(rat, key=lambda t: t[0]); mn = min(rat, key=lambda t: t[0])
    nviol = sum(1 for t in rat if t[0] < 1)
    print("  LCG at fake angles: violated at %d/%d primes; ratio min %.4f (alpha=%s), max %.4f (alpha=%s, Norm=%d <= p=%d)" %
          (nviol, r, mn[0], mn[1], mx[0], mx[1], mx[1][0]**2 + mx[1][1]**2, mx[2]))
    # LCG at genuine angles (exact): ratio = |Im(alpha conj pi)|, minimised over the same candidate search
    gmin = None
    for j, p in enumerate(ps):
        a, b = two_squares(p)
        R0 = math.isqrt(p)
        # exact: for every alpha = x+iy with Norm <= p (x dominant coordinate scan), m = Im(pi conj alpha) = b x - a y;
        # the minimal nonzero |m| over the 3 nearest y per x is the LCG ratio of the genuine angle.
        xs = np.arange(-R0, R0 + 1, dtype=np.int64)
        ms = []
        for dy in (-1, 0, 1):
            ys = (b * xs) // a + dy
            ok = xs * xs + ys * ys <= p
            m = np.abs(b * xs[ok] - a * ys[ok])
            ms.append(m[m > 0].min())
        m = int(min(ms))
        assert m >= 1
        gmin = m if gmin is None else min(gmin, m)
    print("  LCG at genuine angles arg(pi_j): best candidate has |Im(alpha conj pi)| = %d >= 1 (exact) for all j" % gmin)
    # fibre integrality, k = 1:  prod_j 2 p_j^2 cos(4 phi_j) must be an integer for actual clusters
    mp.dps = dps
    prodf = mpf(1)
    for j, p in enumerate(ps): prodf *= 2 * mpf(p) ** 2 * cos(4 * phi[j])
    dist = fabs(prodf - nint(prodf))
    single = [fabs(2 * mpf(p) ** 2 * cos(4 * phi[j]) - nint(2 * mpf(p) ** 2 * cos(4 * phi[j]))) for j, p in enumerate(ps)]
    print("  fibre integrality N^2 S_1(N)/4 = prod_j 2 p_j^2 cos(4 phi_j): distance to Z = %.3f (|value| ~ 10^%d); "
          "per-prime distances min %.3f" % (float(dist), int(float(log(fabs(prodf), 10))), float(min(single))))
    sys.stdout.flush()

if __name__ == '__main__':
    # (i) construct.md fake: Paley-I q=11, b=8, primes near 10^7 (as instance.py 11 8 1 2 7), window <= 1/(4r)
    q, b = 11, 8
    S1 = construct_profile(q, b); M1, r1 = S1.shape
    ps1 = primes_window(10**7, r1, 1 / (4 * r1))
    assert min(ps1) >= 21 * r1 * r1
    digits = int(sum(math.log10(2.0 * p * p) for p in ps1)) + 60
    run("construct.md Paley-I flipped q=11 b=8", S1, ps1, mpf(1), [0] * M1, 11, digits, 5000)
    # (ii) skew double core, lw.md Thm 6.1 parameters: P >= (40M/C0)^4, window log-width 1/(8q), C0 = 1, k = e_1
    for q in (3, 7, 11):
        S2 = skew_double_profile(q); M2, r2 = S2.shape
        P = (40 * M2) ** 4
        ps2 = primes_window(P, r2, 1 / (8 * q))
        kv = [0] * M2; kv[1] = 1
        digits = int(sum(math.log10(2.0 * p * p) for p in ps2)) + 60
        run("skew double core q=%d (lw.md 6.1, P=(40M)^4)" % q, S2, ps2, mpf(1), kv, 100 + q, digits)
