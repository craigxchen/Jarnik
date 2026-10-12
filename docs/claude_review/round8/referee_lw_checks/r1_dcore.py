"""Referee check R1 (independent of lw_checks/dcore.py): the Paley skew double core.

Builds K from scratch (border index 'inf' stored as index q, F_q elements as 0..q-1), and checks exactly:
  (a) K skew-Hadamard, K(x,inf) = 1;
  (b) the identities (5.1): <c,s1_a> = (c_a - c_inf) + Sigma_a,  <c,s2_a> = -(c_a - c_inf) + Sigma_a;
  (c) Sigma_a = R(a) mod q with R(X) = sum_y c_y (X-y)^d, d = (q-1)/2;
  (d) Theorem 5.1's bound s(c) >= M/2 on random c AND on adversarial 'Case 2' vectors c built so that
      p_k = sum_y c_y y^k = 0 mod q for k = 0..d (c_y = g(y) + q*noise, deg g <= (q-3)/2), where the
      proof claims s >= q+3 when c is primitive and c != const mod q;
  (e) structure of the column character: values of v(K(.,a)) and its support size M/2+1;
      label pairs and row pairs.
Exact integer arithmetic only.  Usage: python3 r1_dcore.py
"""
import random, sys
random.seed(20261012)
out = []
def log(s):
    print(s); out.append(s)

def chi(a, q):
    a %= q
    if a == 0: return 0
    return 1 if pow(a, (q-1)//2, q) == 1 else -1

def build(q):
    M = q + 1; INF = q
    K = [[0]*M for _ in range(M)]
    for x in range(q):
        for y in range(q):
            K[x][y] = (1 if x == y else 0) + chi(y - x, q)
        K[x][INF] = 1      # K(x, inf) = 1  (first column ones)
        K[INF][x] = -1     # K(inf, y) = -1
    K[INF][INF] = 1
    for x in range(M):
        for y in range(M):
            assert K[x][y] + K[y][x] == (2 if x == y else 0)
            assert sum(K[x][k]*K[y][k] for k in range(M)) == (M if x == y else 0)
    s1 = [[K[x][a] for x in range(M)] for a in range(q)]                 # columns a of K (K(x,inf)=1)
    s2 = [[K[a][x]*K[INF][x] for x in range(M)] for a in range(q)]       # K^T(x,a) K^T(x,inf)
    for a in range(q):  # s2 = s1 with rows inf and a negated
        for x in range(M):
            assert s2[a][x] == (-s1[a][x] if x in (INF, a) else s1[a][x])
    return K, s1, s2

def supp(c, types):
    return sum(1 for s in types if sum(ci*si for ci, si in zip(c, s)) != 0)

def check_q(q, ntrials):
    M = q + 1; INF = q; d = (q-1)//2
    K, s1, s2 = build(q)
    types = s1 + s2
    # (b), (c)
    for _ in range(ntrials):
        c = [random.randint(-5, 5) for _ in range(M)]
        c[INF] -= sum(c)
        for a in range(q):
            Sig = sum(chi(a - y, q)*c[y] for y in range(q))
            u = sum(c[x]*s1[a][x] for x in range(M)); v = sum(c[x]*s2[a][x] for x in range(M))
            assert u == (c[a] - c[INF]) + Sig and v == -(c[a] - c[INF]) + Sig
            R = sum(c[y]*pow(a - y, d, q) for y in range(q)) % q
            assert (Sig - R) % q == 0
    # (d) random c of various sparsities
    mins = M*10
    for _ in range(ntrials):
        k = random.randint(2, M)
        rows = random.sample(range(M), k)
        c = [0]*M
        for r in rows: c[r] = random.choice([-3, -2, -1, 1, 2, 3])
        c[rows[0]] -= sum(c)
        if all(x == 0 for x in c): continue
        s = supp(c, types); mins = min(mins, s)
        assert s >= (q+1)//2, (q, c, s)
    # (d') adversarial Case-2 vectors
    case2_min = M*10; n2 = 0
    for _ in range(ntrials):
        deg = random.randint(0, (q-3)//2)
        g = [random.randint(0, q-1) for _ in range(deg+1)]
        c = [0]*M
        for y in range(q):
            c[y] = sum(gk*pow(y, k, q) for k, gk in enumerate(g)) % q + q*random.randint(-1, 1)
        c[INF] = -sum(c[:q])
        # check p_k = 0 mod q for k <= d
        assert all(sum(c[y]*pow(y, k, q) for y in range(q)) % q == 0 for k in range(d+1))
        from math import gcd
        gg = 0
        for x in c: gg = gcd(gg, abs(x))
        if gg == 0: continue
        c = [x//gg for x in c]
        nonconst = any((c[y] - c[INF]) % q for y in range(q))
        s = supp(c, types)
        assert s >= (q+1)//2
        if nonconst and gg == 1:
            n2 += 1; case2_min = min(case2_min, s)
            assert s >= q + 3, (q, c, s)
    # (e) column characters, label pairs, row pairs
    col_s = set(); col_vals = set()
    for a in range(q):
        c = [K[x][a] for x in range(M)]
        img = [sum(ci*si for ci, si in zip(c, s)) for s in types]
        col_s.add(sum(1 for t in img if t)); col_vals.add(tuple(sorted(set(abs(t)//2 for t in img if t))))
    lp = set()
    for a in range(q):
        for b in range(a+1, q):
            c = [(K[x][a] - K[x][b])//2 for x in range(M)]
            lp.add(supp(c, types))
    rp = set()
    for x in range(M):
        for y in range(x+1, M):
            c = [0]*M; c[x] = 1; c[y] = -1
            rp.add(supp(c, types))
    log(f"q={q} M={M}: K skew-Hadamard OK; (5.1) and Sigma=R(a) mod q OK on {ntrials} random c; "
        f"min s over {ntrials} random sparse c = {mins} (bound M/2={M//2}); "
        f"Case-2 adversarial primitive nonconst c: {n2} tested, min s = {case2_min} (claimed >= q+3={q+3}); "
        f"columns: s={sorted(col_s)}, |v| values={sorted(col_vals)}; label pairs s={sorted(lp)}; row pairs s={sorted(rp)}")

if __name__ == '__main__':
    for q, nt in ((7, 400), (11, 400), (19, 300), (23, 300), (31, 200), (43, 150), (47, 100), (59, 60)):
        check_q(q, nt)
    open('r1_dcore_out.txt', 'w').write('\n'.join(out) + '\n')
