"""Referee checks for IX.3 (factor2): nonpair bounds, real correction, spread windings.
Exact integer / Fraction arithmetic.  python3 paley_checks.py q seed
"""
import sys, random
from fractions import Fraction as Fr
from itertools import combinations

q = int(sys.argv[1]) if len(sys.argv) > 1 else 43
seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1
random.seed(seed)
sq = {(i*i) % q for i in range(1, q)}
chi = lambda t: 0 if t % q == 0 else (1 if t % q in sq else -1)
K = [[chi(j - i) for j in range(q)] for i in range(q)]
Acls = [i for i in range(1, q) if K[i][0] == 1]
I = list(range(1, q))

def run(S, label):
    m = len(S)
    # (1) nonpair inequalities on random and structured characters
    worst = {}
    viol = 0
    def test(c):
        nonlocal viol
        n = sum(abs(t) for t in c.values())
        if n == 2: return   # pairs are treated separately in IX.3
        v = [sum(c[i]*K[i][j] for i in c) for j in range(q)]
        assert v[0] == 0
        h = sum(abs(t) for t in v[1:])
        c2 = sum(t*t for t in c.values())
        ok = (h >= 1.2*q - 1e-9) and ((n+1)*h >= q*c2 + n*q/2) and (h*m >= q*n)
        if not ok: viol += 1
        key = 'n<=q/10' if n <= q/10 else ('q/10<n<=m' if n <= m else 'n>m')
        val = (h - n)/q
        if key not in worst or val < worst[key][0]: worst[key] = (val, n, h)
    # exhaustive n = 4 with entries +-1
    cnt = 0
    for a, b in combinations(S, 2):
        for cc, dd in combinations(S, 2):
            if len({a, b, cc, dd}) < 4 or (a, b) > (cc, dd): continue
            test({a: 1, b: 1, cc: -1, dd: -1}); cnt += 1
    for a in S:      # 2e_a - e_b - e_c
        for b, cc in combinations([s for s in S if s != a], 2):
            test({a: 2, b: -1, cc: -1}); cnt += 1
    # random characters of many lengths, and balanced +-1 vectors on subsets
    for _ in range(3000):
        k = random.randint(2, m)
        sub = random.sample(S, k)
        c = {i: random.choice([-3, -2, -1, 1, 2, 3]) for i in sub}
        s = sum(c.values())
        # fix the sum to zero
        i0 = sub[0]; c[i0] -= s
        if c[i0] == 0: del c[i0]
        if not c: continue
        test(c); cnt += 1
    for _ in range(500):
        k = 2*random.randint(1, m//2)
        sub = random.sample(S, k)
        c = {i: (1 if t < k//2 else -1) for t, i in enumerate(sub)}
        test(c); cnt += 1
    print(label, "q", q, "m", m, "characters tested", cnt, "violations of the three VIII.3 bounds:", viol)
    print("   min (h - n)/q by range:", {k: (round(v[0], 3), v[1], v[2]) for k, v in worst.items()})

    # (3) real correction: t_j = (2q)^{-1} sum_{i in S} K_ij eps_i ; check 2 K_S t = eps - (2/q) sum eps
    eps = {i: Fr(random.randint(0, 10**6), 10**6) for i in S}
    t = {j: sum(K[i][j]*eps[i] for i in S)/(2*q) for j in I}
    tot = sum(eps.values())
    for x in S:
        assert 2*sum(K[x][j]*t[j] for j in I) == eps[x] - Fr(2, q)*tot
    assert all(abs(t[j]) <= max(eps.values()) / 2 for j in I)
    # character evaluation: v.t = c.eps/2 for zero-sum c
    for _ in range(50):
        sub = random.sample(S, 6); c = {sub[0]: 2, sub[1]: -1, sub[2]: 1, sub[3]: -1, sub[4]: 1, sub[5]: -2}
        v = {j: sum(c[i]*K[i][j] for i in c) for j in I}
        assert sum(v[j]*t[j] for j in I) == sum(c[i]*eps[i] for i in c)/2
    print("   real correction: 2K_S t = eps - (2/q)sum(eps) and v.t = c.eps/2 verified exactly; max|t| <= max eps/2")

    # (3') spread windings: b_j = 1/12 + (j/q + z_j)/3 with z built as in VIII.4 (T random small)
    tau1 = [sum(K[i][j]*j for j in range(q)) for i in range(q)]
    assert all(t_ % q == 0 for t_ in tau1)
    tau1 = [t_ // q for t_ in tau1]
    Sset = set(S); Sc = [i for i in range(q) if i not in Sset]
    T = random.sample(I, 3)
    # solve over F_3: choose u_S = tau1_S, u_{S^c} with (K^T u)_T = 0 mod 3
    rows = [[K[i][j] % 3 for i in Sc] for j in T]
    rhs = [(-sum(K[i][j]*tau1[i] for i in S)) % 3 for j in T]
    # Gaussian elimination mod 3
    nvar = len(Sc); M = [r + [b] for r, b in zip(rows, rhs)]; piv = []; r0 = 0
    for col in range(nvar):
        pr = next((r for r in range(r0, len(M)) if M[r][col] % 3), None)
        if pr is None: continue
        M[r0], M[pr] = M[pr], M[r0]
        inv = 1 if M[r0][col] % 3 == 1 else 2
        M[r0] = [(x*inv) % 3 for x in M[r0]]
        for r in range(len(M)):
            if r != r0 and M[r][col]:
                f = M[r][col]; M[r] = [(x - f*y) % 3 for x, y in zip(M[r], M[r0])]
        piv.append(col); r0 += 1
        if r0 == len(M): break
    assert r0 == len(T), "rank deficiency on T"
    uSc = [0]*nvar
    for r, col in enumerate(piv): uSc[col] = M[r][-1]
    u = [0]*q
    for i in S: u[i] = tau1[i]
    for k, i in enumerate(Sc): u[i] = uSc[k]
    KTu = [sum(K[i][j]*u[i] for i in range(q)) for j in range(q)]
    assert all(KTu[j] % 3 == 0 for j in T)
    qinv3 = pow(q, -1, 3)
    z = [(-qinv3*KTu[j]) % 3 for j in range(q)]
    assert all(z[j] == 0 for j in T)
    b = {j: Fr(1, 12) + (Fr(j, q) + z[j])/3 for j in I}
    Kb = {x: sum(K[x][j]*b[j] for j in I) for x in S}
    x0 = S[0]
    diffs = [Kb[x] - Kb[x0] for x in S]
    assert all(d.denominator == 1 for d in diffs), "Kb differences not integral"
    print("   spread: (Kb)_x - (Kb)_x0 integral for all marks => windings 2Kb have EVEN differences; T =", T,
          "; sample diffs", [int(d) for d in diffs[:6]])

run(Acls, "S = A (rho = 1)")
S = [i for i in Acls if random.random() < 0.9]
run(S, "S random, density 0.9")
