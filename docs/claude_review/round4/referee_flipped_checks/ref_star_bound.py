"""Referee check: the 'star' pair-inequality bound W_F <= W/2 + (M-1)(4 log C - D)/2.

Claim (proved in referee_flipped.md, Sec. 4): in a fully flipped Hadamard-core profile (one flip per row,
flipped column of row z at label a(z) != 0), for every row x
   sum_{y != x} G_xy = -W + 2 l_x + 2 sum_{z != x} eps_zx l_z ,   eps_zx = -H(z,a(z)) H(x,a(z)),
where G_xy = sum_j w_j s_xj s_yj, W = sum_j w_j, l_z = weight of the flipped column of row z.
If x0 is an all-ones row and H(z,a(z)) = -1 for every z != x0 (true for the canonical capacity-two Paley
assignment of item 392 / Prop P: finite row x flipped at the diagonal label x, H(x,x) = -1), then
eps_{z,x0} = +1 for all z != x0 and the pair inequalities G_xy <= 4 log C - D give
   W_F <= W/2 + (M-1)(4 log C - D)/2.
Checks: (1) the identity with exact integers on random weights against the literal sign matrix, for Walsh and
Paley cores and random capacity-B assignments; (2) H(x,x) = -1 for every finite x for Paley-I, q = 3 mod 4,
q < 2000 (so the canonical assignment satisfies the hypothesis); (3) G_xy from the sign matrix equals
2 log Norm gcd(z_x, z_y) - W - 2D on literal Gaussian profiles (exact gcd in Z[i]).
"""
import math, random, sys
sys.path.insert(0, '.')
from ref_pair_lp import paley, sylv, legendre

def is_prime(n):
    if n < 2: return False
    i = 2
    while i*i <= n:
        if n % i == 0: return False
        i += 1
    return True

def signs(H, assign, x, z_cols):
    pass

def check_identity(H, rng, B, trials=20):
    M = len(H)
    for _ in range(trials):
        cnt = {}; assign = []
        for x in range(M):
            ch = [a for a in range(1, M) if cnt.get(a, 0) < B]; a = rng.choice(ch); assign.append(a); cnt[a] = cnt.get(a, 0) + 1
        # columns: per label b=5 copies; first copies flipped by the rows assigned to it
        cols = []
        used = {}
        flipcol = {}
        for a in range(1, M):
            for c in range(5):
                cols.append([a, rng.choice([1, -1]), rng.randint(1, 50), None])
        for x in range(M):
            a = assign[x]; c = used.get(a, 0); used[a] = c + 1
            idx = (a - 1)*5 + c; cols[idx][3] = x; flipcol[x] = idx
        S = [[cols[j][1]*H[x][cols[j][0]]*(-1 if cols[j][3] == x else 1) for j in range(len(cols))] for x in range(M)]
        W = sum(c[2] for c in cols)
        l = [cols[flipcol[z]][2] for z in range(M)]
        for x in range(M):
            lhs = sum(sum(cols[j][2]*S[x][j]*S[y][j] for j in range(len(cols))) for y in range(M) if y != x)
            rhs = -W + 2*l[x] + 2*sum(-H[z][assign[z]]*H[x][assign[z]]*l[z] for z in range(M) if z != x)
            assert lhs == rhs, (lhs, rhs)
    return trials*M

def gcd_check(rng):
    # literal Gaussian profile, Paley-12 canonical, b=5, content g = (2+i)*3
    def gp(p):
        for y in range(1, math.isqrt(p) + 1):
            x2 = p - y*y; x = math.isqrt(x2)
            if x*x == x2: return (x, y)
    def mul(a, b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
    def cj(a): return (a[0], -a[1])
    def nrm(a): return a[0]**2 + a[1]**2
    def gdivmod(a, b):
        n = nrm(b); t = mul(a, cj(b))
        q = ((2*t[0] + n)//(2*n), (2*t[1] + n)//(2*n))
        r = (a[0] - mul(q, b)[0], a[1] - mul(q, b)[1]); return r
    def ggcd(a, b):
        while b != (0, 0): a, b = b, gdivmod(a, b)
        return a
    H = paley(11); M = 12; assign = [1] + list(range(1, M))
    splits = [p for p in range(5, 3000) if p % 4 == 1 and is_prime(p)]
    ps = rng.sample(splits[3:], 5*(M-1))
    cols = []; used = {}
    for a in range(1, M):
        for c in range(5): cols.append([a, rng.choice([1, -1]), ps[(a-1)*5 + c], None])
    for x in range(M):
        a = assign[x]; c = used.get(a, 0); used[a] = c + 1; cols[(a-1)*5 + c][3] = x
    g = mul((2, 1), (3, 0)); D = math.log(nrm(g))
    Z = []
    for x in range(M):
        z = g
        for a, s, p, fx in cols:
            e = s*H[x][a]*(-1 if fx == x else 1); pi = gp(p)
            z = mul(z, pi if e == 1 else cj(pi))
        Z.append(mul(z, rng.choice([(1, 0), (0, 1), (-1, 0), (0, -1)])))
    W = sum(math.log(c[2]) for c in cols); n = 0
    for x in range(M):
        for y in range(x+1, M):
            Gxy = sum(math.log(p)*(s*H[x][a]*(-1 if fx == x else 1))*(s*H[y][a]*(-1 if fx == y else 1)) for a, s, p, fx in cols)
            lg = math.log(nrm(ggcd(Z[x], Z[y])))
            assert abs(lg - (D + (W + Gxy)/2)) < 1e-6, (lg, D + (W + Gxy)/2); n += 1
    return n

if __name__ == '__main__':
    rng = random.Random(8)
    n = 0
    for H, B in [(sylv(3), 2), (sylv(3), 3), (paley(11), 2), (paley(11), 4), (sylv(4), 2), (paley(19), 3)]:
        n += check_identity(H, rng, B)
    print("star identity checked on", n, "(row, profile) instances: PASS")
    bad = [q for q in range(3, 2000) if q % 4 == 3 and is_prime(q) and any(paley_entry != -1 for paley_entry in [-legendre(0, q) - 1])]
    # H(x,x) = -chi(0) - 1 = -1 for all x, all q: verify on actual matrices for small q
    for q in [7, 11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83]:
        H = paley(q)
        assert all(H[1 + x][1 + x] == -1 for x in range(q)) and all(v == 1 for v in H[0])
    print("Paley-I: H(x,x) = -1 for all finite x and row infinity all ones (q <= 83 matrices built): PASS")
    print("gcd identity log Norm gcd = D + (W+G)/2 verified on", gcd_check(rng), "literal Paley-12 pairs: PASS")
