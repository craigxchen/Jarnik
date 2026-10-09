"""Sweep of the pair-inequality LP (see ref_pair_lp.py) over many assignments and capacities.

For each instance: rational simplex optimum rho = max W_F/W subject to all pair inequalities G_xy <= 0.
Also returns the exact dual certificate and K = sum of pair multipliers, so that with the true right-hand
side gamma = 4 log C - D the certificate gives   W_F <= rho W + K gamma   (exact LP duality).
The certificate is re-verified independently: sum_pairs mu_xy * coef_xy  >=  (l-coeff 1, U-coeff 0) - rho*(1,...,1).
"""
from fractions import Fraction as Fr
import random, sys
sys.path.insert(0, '.')
from ref_pair_lp import paley, sylv

def simplex(A, b, c):
    m, n = len(A), len(c)
    T = [[Fr(v) for v in A[i]] + [Fr(int(i == k)) for k in range(m)] + [Fr(b[i])] for i in range(m)]
    obj = [Fr(-v) for v in c] + [Fr(0)]*m + [Fr(0)]
    basis = [n + i for i in range(m)]
    while True:
        e = next((j for j in range(n + m) if obj[j] < 0), None)
        if e is None: break
        best, r = None, None
        for i in range(m):
            if T[i][e] > 0:
                ratio = T[i][-1] / T[i][e]
                if best is None or ratio < best or (ratio == best and basis[i] < basis[r]):
                    best, r = ratio, i
        piv = T[r][e]; T[r] = [v / piv for v in T[r]]
        for i in range(m):
            if i != r and T[i][e] != 0:
                f = T[i][e]; T[i] = [vi - f*vr for vi, vr in zip(T[i], T[r])]
        if obj[e] != 0:
            f = obj[e]; obj = [vo - f*vr for vo, vr in zip(obj, T[r])]
        basis[r] = e
    return obj[-1], obj[n:n+m]

def build(H, assign):
    M = len(H); nv = (M - 1) + M; rows = []
    for x in range(M):
        for y in range(x + 1, M):
            coef = [0]*nv
            for a in range(1, M):
                h = H[x][a]*H[y][a]; coef[a - 1] += h
                for x2 in range(M):
                    if assign[x2] == a: coef[M - 1 + x2] += h
            coef[M - 1 + x] += -2*H[x][assign[x]]*H[y][assign[x]]
            coef[M - 1 + y] += -2*H[y][assign[y]]*H[x][assign[y]]
            rows.append(coef)
    return rows

def solve(H, assign):
    M = len(H); nv = 2*M - 1
    rows = build(H, assign)
    val, dual = simplex(rows + [[1]*nv], [0]*len(rows) + [1], [0]*(M - 1) + [1]*M)
    mu, t = dual[:-1], dual[-1]
    assert t == val
    # independent verification of the certificate
    comb = [sum(mu[i]*rows[i][k] for i in range(len(rows))) for k in range(nv)]
    target = [Fr(0)]*(M - 1) + [Fr(1)]*M
    assert all(mu_i >= 0 for mu_i in mu)
    assert all(comb[k] + t >= target[k] for k in range(nv)), "certificate fails"
    return val, sum(mu)

def rand_assign(M, B, rng, surj=True):
    if surj and B >= 2:
        labels = list(range(1, M)) + [rng.randrange(1, M)]; rng.shuffle(labels); return labels[:M]
    cnt = {}; lab = []
    for x in range(M):
        ch = [a for a in range(1, M) if cnt.get(a, 0) < B]; a = rng.choice(ch); lab.append(a); cnt[a] = cnt.get(a, 0) + 1
    return lab

if __name__ == '__main__':
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 11)
    plan = [('Walsh-8', sylv(3), [(2, True, 40), (3, False, 40), (4, False, 40)]),
            ('Paley-12', paley(11), [(2, True, 40), (3, False, 30), (5, False, 30)]),
            ('Walsh-16', sylv(4), [(2, True, 8), (3, False, 6)]),
            ('Paley-20', paley(19), [(2, True, 6), (4, False, 4)]),
            ('Paley-24', paley(23), [(2, True, 3)])]
    overall = Fr(0)
    for name, H, specs in plan:
        M = len(H)
        # canonical capacity-two assignment
        v, K = solve(H, [1] + list(range(1, M)))
        print(f"{name} canonical: rho={v} ({float(v):.4f}), K={float(K):.3f}", flush=True)
        overall = max(overall, v)
        for B, surj, n in specs:
            vals = []
            for _ in range(n):
                v, K = solve(H, rand_assign(M, B, rng, surj)); vals.append((v, K))
            mx = max(v for v, _ in vals); overall = max(overall, mx)
            print(f"{name} B={B} {'surjective' if surj else 'random'} x{n}: rho in [{float(min(v for v,_ in vals)):.4f}, {float(mx):.4f}],"
                  f" max K={float(max(K for _,K in vals)):.3f}", flush=True)
    print("max rho over all instances:", overall, float(overall))
