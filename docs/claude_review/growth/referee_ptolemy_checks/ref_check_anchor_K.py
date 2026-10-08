"""Referee checks:
 (L5.2) anchoring lemma, independent implementation (cluster discs with anchor at infinity);
 (K)    canonical class of Mbar_{0,n}: boundary formula vs Kapranov formula, and K.beta closed form;
 (D)    the numerical constants of Steps 3-5 of Theorem D, in exact rational arithmetic.
"""
import math, itertools, random
from fractions import Fraction as F
from ref_common import factor, vp, det2

random.seed(7)
# ---------------- (L5.2)
cnt = 0; pchecks = 0
for trial in range(300):
    k = random.randint(4, 8); m = k - 1
    pts = set()
    while len(pts) < m:
        x, y = random.randint(-10**4, 10**4), random.randint(1, 10**4)
        g = math.gcd(x, y); pts.add((x // g, y // g))
    w = [(1, 0)] + list(pts)
    # optionally plant p-adic clustering to create deep trees
    D = 1
    for j in range(1, k):
        D = D * abs(w[j][1]) // math.gcd(D, abs(w[j][1]))
    P = []
    for (x, y) in w:
        u, v = D * x, y
        g = math.gcd(u, v); P.append((u // g, v // g))
    assert P[0] == (1, 0)
    for j in range(1, k):
        assert abs(det2(P[0], P[j])) == 1
    # cross-ratios unchanged
    for a_, b_, c_, d_ in itertools.combinations(range(k), 4):
        cr = lambda Q: F(det2(Q[a_], Q[b_]) * det2(Q[c_], Q[d_]), det2(Q[a_], Q[c_]) * det2(Q[b_], Q[d_]))
        assert cr(w) == cr(P)
    primes = set()
    for i, j in itertools.combinations(range(k), 2):
        primes |= set(factor(det2(w[i], w[j])))
    bs = [frozenset(T) for s in range(2, m) for T in itertools.combinations(range(1, k), s)]
    for p in primes:
        # cluster tree of x_i = X_i/Y_i (anchor 0 = infinity) from ORIGINAL w
        def vd(i, j):
            return vp(det2(w[i], w[j]), p) - vp(w[i][1], p) - vp(w[j][1], p)
        ell = {}
        for T in bs:
            rT = min(vd(i, j) for i, j in itertools.combinations(sorted(T), 2))
            sT = max(vd(i, j) for i in T for j in range(1, k) if j not in T)
            ell[T] = max(0, rT - sT)
        deltas = set()
        for i, j in itertools.combinations(range(1, k), 2):
            deltas.add(vp(det2(P[i], P[j]), p) - sum(ell[T] for T in bs if i in T and j in T))
        assert len(deltas) == 1 and min(deltas) >= 0, (p, deltas)
        pchecks += 1
    cnt += 1
print(f"(L5.2) anchoring lemma: {cnt} configurations, {pchecks} prime checks: OK")

# ---------------- (K)
def c(s, n):
    return F(s * (n - s), n - 1) - 2
for n in range(5, 13):
    # boundary splits: unordered {S,S^c}, 2<=|S|<=n-2
    splits = []
    for s in range(2, n // 2 + 1):
        for S in itertools.combinations(range(n), s):
            if s == n - s and 0 not in S:
                continue
            splits.append(frozenset(S))
    assert len(splits) == 2**(n - 1) - n - 1
    Kbeta = sum(c(len(S), n) for S in splits)
    assert Kbeta == 2**(n - 3) * (n - 8) + n + 2
    assert sum(c(len(S), n) + 1 for S in splits) == 2**(n - 3) * (n - 4) + 1
    # Kapranov: marked point n-1 is psi-point; E_I (I subset [n-1], 1<=|I|<=n-4) = D_{I u {n}},
    # hyperplane through p_I (|I| = n-3) = H - sum_{J subsetneq I} E_J = D_{I u {n}}.
    last = n - 1
    coefH = F(0); coefE = {}
    for S in splits:
        Sn = S if last in S else frozenset(range(n)) - S
        I = Sn - {last}
        if len(I) == n - 3:
            coefH += c(len(S), n)
            for r in range(1, n - 3):
                for J in itertools.combinations(sorted(I), r):
                    coefE[frozenset(J)] = coefE.get(frozenset(J), 0) - c(len(S), n)
        else:
            coefE[I] = coefE.get(I, 0) + c(len(S), n)
    assert coefH == -(n - 2)
    for I, v in coefE.items():
        assert v == n - 3 - len(I), (n, I, v)
    print(f"(K) n={n}: boundary formula == Kapranov K = -(n-2)H + sum (n-3-|I|)E_I ;  K.beta = {Kbeta}")

# ---------------- (D) constants
n = 8; m = 7
pos = sum(c(len(S), n) for S in splits if False)  # placeholder
cs = {2: c(2, 8), 3: c(3, 8), 4: c(4, 8)}
num = {2: 28, 3: 56, 4: 35}
P_ = cs[3] * num[3] + cs[4] * num[4]; N_ = -cs[2] * num[2]
assert (P_, N_) == (18, 8)
up = 2**(m - 2)  # 32
# h_K >= P_((1-eta)w - 2H) - N_((1+eta+up*eta)w + 3H)
w_coef = P_ - N_; eta_coef = P_ + N_ * (1 + up); H_coef = 2 * P_ + 3 * N_
assert (w_coef, eta_coef, H_coef) == (10, 282, 60)
Asum = sum((cs[s] + 1) * num[s] for s in cs)
assert Asum == 129
eps = F(1, 43)
rhs_w = eps * Asum; rhs_eta = eps * Asum * (1 + up); rhs_H = eps * Asum * 3
assert (rhs_w, rhs_eta, rhs_H) == (3, 99, 9)
assert (w_coef - rhs_w, eta_coef + rhs_eta, H_coef + rhs_H) == (7, 381, 69)
# contradiction threshold: (6w - C')/69 <= 3584 w/(M-1), with w >= C' -> 5w/69 <= 3584 w/(M-1)
assert F(3584 * 69, 5) == F(247296, 5)
assert 28 * 128 == 3584
eta_need = (2 * math.sqrt(2) * 255 * 381) ** 2 + 7
print(f"(D) Step 3-5 constants verified exactly; M threshold for eta<=1/381: {eta_need:.3e}")
