# E_Lambda . C (degree of t -> Q_{x(t)}) modulo a prime by evaluation/interpolation (fast).
# C = Mobius curve (pure, or with one triple event T at t = 0). Arrangement and event points from arr.py.
import sys, random, itertools
from fractions import Fraction as Fr
import flint
from star import point_of, quadric_through, LAB
from arr import PLANES, CHOSEN
P = 2305843009213693951  # 2^61-1


def md(q):
    q = Fr(q)
    return q.numerator % P * pow(q.denominator % P, P - 2, P) % P


def det_mod(M):
    M = [row[:] for row in M]
    n = len(M)
    d = 1
    for c in range(n):
        piv = next((r for r in range(c, n) if M[r][c] % P), None)
        if piv is None:
            return 0
        if piv != c:
            M[c], M[piv] = M[piv], M[c]
            d = -d
        d = d * M[c][c] % P
        inv = pow(M[c][c], P - 2, P)
        for r in range(c + 1, n):
            f = M[r][c] * inv % P
            if f:
                for k in range(c, n):
                    M[r][k] = (M[r][k] - f * M[c][k]) % P
    return d % P


pl = [[md(c) for c in H] for H in PLANES]
pts = [point_of([PLANES[i], PLANES[j], PLANES[k]]) for (i, j, k) in CHOSEN]
qb = [[[md(Q[i, j]) for j in range(4)] for i in range(4)] for Q in quadric_through(pts)]
assert len(qb) == 4


def quad(Q, u, v):
    return sum(u[i] * Q[i][j] * v[j] for i in range(4) for j in range(4)) % P


def kernel_at(mob, t0):
    rows = []
    for a in range(7):
        al, be, ga, de = mob[a]
        N = (al * t0 + be) % P
        D = (ga * t0 + de) % P
        rows.append([N * pl[a][i] % P for i in range(4)] + [D * pl[a][i] % P for i in range(4)])
    return [(-1) ** j * det_mod([[r[k] for k in range(8) if k != j] for r in rows]) % P for j in range(8)]


def newton_interp(xs, ys):
    n = len(xs)
    coef = list(ys)
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            coef[i] = (coef[i] - coef[i - 1]) * pow((xs[i] - xs[i - j]) % P, P - 2, P) % P
    poly = flint.nmod_poly([coef[-1]], P)
    X = flint.nmod_poly([0, 1], P)
    for i in range(n - 2, -1, -1):
        poly = poly * (X - xs[i]) + coef[i]
    return poly


def degE(mob, N=40):
    xs = list(range(3, 3 + N))            # kernel entries have degree <= 7
    kers = [kernel_at(mob, x) for x in xs]
    K = [newton_interp(xs, [k[j] for k in kers]) for j in range(8)]
    g = K[0]
    for k in K[1:]:
        g = g.gcd(k)
    K = [k // g for k in K]
    dk = max(k.degree() for k in K)
    xs2 = list(range(5, 5 + 6 * dk + 10))   # minors have degree <= 6 dk
    vals = []
    for x in xs2:
        L1 = [int(K[i](x)) for i in range(4)]
        L2 = [int(K[4 + i](x)) for i in range(4)]
        C = [[quad(Q, u, v) for Q in qb] for (u, v) in [(L1, L1), (L1, L2), (L2, L2)]]
        vals.append([(-1) ** j * det_mod([[r[k] for k in range(4) if k != j] for r in C]) % P for j in range(4)])
    m = [newton_interp(xs2, [v[j] for v in vals]) for j in range(4)]
    g = m[0]
    for k in m[1:]:
        g = g.gcd(k)
    m = [k // g for k in m]
    return max(k.degree() for k in m), dk


def rnd():
    while True:
        n = random.randint(-60, 60)
        if n:
            return Fr(n, random.randint(1, 25))


def good_mob(T):
    while True:
        c = [[rnd() for _ in range(4)] for a in range(7)]
        if T:
            a0 = T[0]
            z0 = c[a0][1] / c[a0][3]
            for b in T[1:]:
                c[b][1] = z0 * c[b][3]
        if any(ci[0] * ci[3] - ci[1] * ci[2] == 0 for ci in c):
            continue
        return [[md(v) for v in ci] for ci in c]


if __name__ == '__main__':
    random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
    d0 = [degE(good_mob(()))[0] for _ in range(2)]
    print('mu', d0, flush=True)
    tot = 0
    for T in itertools.combinations(range(7), 3):
        ds = [degE(good_mob(T))[0] for _ in range(2)]
        print([LAB[a] for a in T], ds, flush=True)
        tot += ds[0]
    print('E.mu =', d0[0], ' sum_T E.mu_T =', tot, ' E.beta = sum - 32 E.mu =', tot - 32 * d0[0])
    print('lambda_A(beta) = 7 - 2 E.beta =', 7 - 2 * (tot - 32 * d0[0]))
