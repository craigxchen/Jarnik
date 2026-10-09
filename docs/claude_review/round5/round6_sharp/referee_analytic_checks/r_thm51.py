"""Referee checks for sharp.md Theorem 5.1 (P3-char) and Lemmas 3.1, 3.3, 5.1a.

 (a) own implementation of the prime assignment p'(q) (Lemma 3.1) and of a*(q); exact
     log m_struct = sum_{q<4M} a*(q) log q against 11 M and against pi(4M) log(32 M^2);
     exact max pair demand max_d sum_{q: p'|d} (1 + min(v_q(d), a*-1)) log q against
     mu_M log M = (10 + 26/loglog M) log M.
 (b) Rosser-Schoenfeld ratio R(k) (Lemma 3.1(2)) re-evaluated.
 (c) the small-q moment E[q^(theta (a-a*)^+)] <= 1 + 2 q^(theta-1) claimed in Thm 5.1 step 4:
     exact worst case (Pr[>= j] = q^(-j)) at q = 3, 5 and theta = 1/2.
 (d) where the union bound of Thm 5.1 (as written: constants 11 M, 66 Sigma_theta, mu_M) closes:
     smallest log M at which each family's total is < 1/4, for b = 33, A = 1, 2, C = 1.
 (e) Monte Carlo of Lemma 3.3(3), large q, level 1: Pr[sum c_x kappa_x = iota s mod m] <= 4 c_min/m.
 (f) Monte Carlo of Lemma 5.1a: Pr[dist(c.delta/2,(pi/4)Z) < arcsin eps] <= 4 eps (n + pi/Delta).
"""
import math, random
import numpy as np

rng = random.Random(5)
nprng = np.random.default_rng(5)
out = []
ok = True


def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


P = primes_upto(2 ** 20)
pp = {3: 2, 5: 2, 7: 3}
for k in range(3, 19):
    B = [q for q in P if 2 ** k <= q < 2 ** (k + 1) and q > 7]
    T = [p for p in P if 2 ** (k - 2) <= p < 2 ** (k - 1)]
    n, t = len(B), len(T)
    for i, q in enumerate(B):
        pp[q] = T[(i * t) // n]


def m1(q):
    return q - 1 if q % 4 == 1 else q + 1


assert all(2 * pp[q] <= m1(q) and q < 8 * pp[q] for q in pp)

# (a) m_struct and pair demand
rows = []
for M in (44, 128, 512, 2048, 8192, 32768, 120000):
    qs = [q for q in P if 2 < q < 4 * M]
    astar = {}
    for q in qs:
        a = 1
        while pp[q] * q ** (a - 1) < M:
            a += 1
        astar[q] = a
    lms = sum(astar[q] * math.log(q) for q in qs)
    piB = len(qs) * math.log(32 * M * M)
    ok &= lms <= 11 * M
    if M <= 32768:
        worst = 0.0
        for d in range(1, M):
            D = 0.0
            for q in qs:
                if d % pp[q] == 0:
                    v, dd = 0, d
                    while dd % q == 0:
                        dd //= q
                        v += 1
                    D += (1 + min(v, astar[q] - 1)) * math.log(q)
            worst = max(worst, D)
        mu = 10 + 26 / math.log(math.log(M))
        ok &= worst <= mu * math.log(M)
        rows.append(f"   M={M:6d}: log m_struct = {lms:9.1f} = {lms/M:.2f} M (claim <= 11 M; pi(4M)log(32M^2) = {piB/M:.2f} M); "
                    f"max pair demand = {worst/math.log(M):.2f} log M (mu_M = {mu:.2f})")
    else:
        rows.append(f"   M={M:6d}: log m_struct = {lms:9.1f} = {lms/M:.2f} M")
out.append("(a) structured part and pair demand:")
out += rows
# robin-route monotonicity remark: (10 + 26/loglog d) log d at small d
vals = [(d, (10 + 26 / math.log(math.log(d))) * math.log(d)) for d in (3, 4, 5, 8, 16)]
out.append("   note: the d-wise Robin form (10+26/loglog d) log d at d = 3,4,5,8,16: " +
           ", ".join(f"{v:.0f}" for d, v in vals) + "  (not monotone; the true demand there is O(1))")

# (b) R(k)
Rmax = 0
for k in range(20, 2001):
    eps = math.log(2) * 2.0 ** (-k)
    num = 2.51012 / (k + 1) - 1 / k + eps
    den = 0.5 / (k - 1) - 0.313765 / (k - 2) - eps
    Rmax = max(Rmax, num / den)
out.append(f"(b) max_(20<=k<=2000) R(k) = {Rmax:.4f}; ceil = {math.ceil(Rmax)}; asymptote 1.51012/0.186235 = {1.51012/0.186235:.4f}")
ok &= Rmax < 9

# (c) small-q moment
for q in (3, 5, 7):
    th = 0.5
    exact = 1 + sum((q ** (th * j) - q ** (th * (j - 1))) * q ** (-j) for j in range(1, 200))
    claim = 1 + 2 * q ** (th - 1)
    out.append(f"(c) q={q}, theta=1/2: worst-case E[q^(theta (a-a*)^+)] = {exact:.4f} vs claimed bound {claim:.4f}  ({'OK' if exact <= claim + 1e-12 else 'VIOLATED'})")

# (d) union bound thresholds (log space), b = 33, C = 1
def thresholds(A, b=33, C=1.0):
    """log-space: M = e^L.  Families normalised by M where they are exponential in M."""
    th = 1.0 / (A + 1)
    kb = min(1.0, b / 16 - 2)
    res = {}
    L = 4.0
    while L < 10 ** 6 and len(res) < 3:
        logr = math.log(b) + L
        tau = max(math.log(3) + A * L, math.log(21) + 2 * logr) + 2 * math.log(L)
        sig_over_M = math.exp(th * math.log(3) + (th * A - 1) * L) / th
        invM = math.exp(-L)
        # ternary: total = exp(M * tern_rate)
        tern_rate = math.log(3) - th * (tau * kb / 4 - 11) + (th * (math.log(4) + L + math.log1p(math.pi / C * invM)) + th / 4) * invM + 66 * sig_over_M
        # amplitude h = 2 (dominant): total = exp(M * amp_rate)
        h = 2
        amp_rate = h * (-(3.5 * th * tau) + math.log(2 * h + 1) / h + 11 * th / h) + 66 * h * sig_over_M + th * (math.log(4 * h) + L) * invM
        # pairs: log total
        mu = 10 + 26 / math.log(L)
        pairs = 2 * L - math.log(2) + math.log(4 * (2 + math.pi / C)) + mu * L - (b - 4) * tau / 4 + 0.25
        if "ternary" not in res and tern_rate < 0:
            res["ternary"] = L
        if "amplitude2" not in res and amp_rate < 0:
            res["amplitude2"] = L
        if "pairs" not in res and pairs < math.log(0.25):
            res["pairs"] = L
        L += 1.0 if L < 2000 else 50.0
    return res


for A in (1, 2):
    t = thresholds(A)
    out.append(f"(d) A={A}, b=33, C=1: union-bound families < 1/4 from log M >= " +
               ", ".join(f"{k}: {v:.0f}" for k, v in t.items()))

# (e) Lemma 3.3(3) Monte Carlo
worst_ratio = 0
for trial in range(300):
    q = rng.choice([p for p in P if 400 < p < 3000])
    m = m1(q)
    Mrows = rng.randint(10, min(60, m // 4))
    sigma = [rng.randint(0, 1) for _ in range(Mrows)]
    cmax = rng.choice([1, 2, 3, 6])
    c = [rng.randint(-cmax, cmax) for _ in range(Mrows)]
    c[0] -= sum(c)
    if not any(c):
        continue
    cmin = min(abs(v) for v in c if v)
    iota = m // 4 if (q % 8 in (1, 7)) else m // 4  # any fixed target; bound is per target
    targets = [(iota * s) % m for s in range(4)]
    hits = [0] * 4
    N = 3000
    for _ in range(N):
        y = rng.sample(range(m // 2), Mrows)
        kap = sum(cx * (sx + 2 * yx) for cx, sx, yx in zip(c, sigma, y)) % m
        for s in range(4):
            if kap == targets[s]:
                hits[s] += 1
    ratio = max(hits) / N / (4 * cmin / m)
    worst_ratio = max(worst_ratio, ratio)
out.append(f"(e) Lemma 3.3(3) level-1 MC: max over trials of empirical Pr / (4 c_min/m) = {worst_ratio:.3f} (small counts, noisy; superseded by r_lemma33_mc.py, which gives 0.55)")

# (f) Lemma 5.1a Monte Carlo
worst = 0
for trial in range(200):
    Mrows = rng.randint(2, 12)
    c = [rng.randint(-3, 3) for _ in range(Mrows)]
    c[0] -= sum(c)
    if not any(c):
        continue
    n = sum(abs(v) for v in c)
    Delta = rng.choice([0.01, 0.1, 0.5, 1.0])
    eps = rng.choice([1e-3, 1e-2, 5e-2])
    d = nprng.random((40000, Mrows)) * Delta
    val = d @ np.array(c, dtype=float) / 2
    dist = np.abs(val - np.round(val / (math.pi / 4)) * (math.pi / 4))
    pr = float((dist < math.asin(min(1, eps))).mean())
    bound = 4 * eps * (n + math.pi / Delta)
    worst = max(worst, pr / bound)
out.append(f"(f) Lemma 5.1a MC: max empirical Pr / bound = {worst:.3f} (<= 1 expected)")
out.append("REFEREE THM 5.1 CHECKS DONE" + ("" if ok else " (FAILURE in (a)/(b))"))
print("\n".join(out))
