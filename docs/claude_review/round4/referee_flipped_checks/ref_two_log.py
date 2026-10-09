"""Referee check (independent code) of Lemma 1.1 / Props 2.1-2.2 of round4/flipped.md.

For random Hadamard-core profiles with ARBITRARY modification sets F_j (several rows per column,
modified label-0 columns), a nontrivial common content g, random orientations and row units:
  (i)  exact Gaussian identity  prod_x zeta_x^{H(x,a)} = +- (G_a/conj G_a)^M (conj P_a/P_a)^2
       (tested as: A_a * conj(G_a^M conj(P_a)^2) is real or purely imaginary), and the pair version
       with (G_a conj G_b)^M conj(Q_ab)^4;
  (ii) the angle bookkeeping: with real lifts phi_j, theta' and the decomposition
       (S phi)_x = theta' + delta_x + (pi/2) k_x, the quantity
       (M Log u_a - 2 Log v_a - 2 i eta_a)/(i pi/2) is an integer, eta_a = sum_x H(x,a) delta_x;
       and the same for pairs with eta_ab = sum_x (H(x,a)-H(x,b)) delta_x;
  (iii) nonvanishing: conj-primitive part of Gamma_a is off the axes and diagonals whenever some
       column of label a has |F_j| < M/2 (pairs: < M/4);
  (iv) minimal polynomial N X^2 - 2 Re(gamma^2) X + N of gamma/conj(gamma) is primitive, H <= 2N.
Own implementation; does not import the author's fcommon.py.
"""
import math, random, cmath

def is_prime(n):
    if n < 2: return False
    i = 2
    while i * i <= n:
        if n % i == 0: return False
        i += 1
    return True

def gp(p):
    for y in range(1, math.isqrt(p) + 1):
        x2 = p - y * y; x = math.isqrt(x2)
        if x * x == x2: return (x, y)
def mul(a, b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def cj(a): return (a[0], -a[1])
def nrm(a): return a[0]*a[0]+a[1]*a[1]
def pw(a, e):
    r = (1, 0)
    if e < 0: a, e = cj(a), -e
    for _ in range(e): r = mul(r, a)
    return r
def prod(l):
    r = (1, 0)
    for t in l: r = mul(r, t)
    return r

def legendre(a, q):
    a %= q
    if a == 0: return 0
    return 1 if pow(a, (q-1)//2, q) == 1 else -1

def paley(q):
    M = q + 1
    H = [[1]*M] + [[1] + [-legendre(a - x, q) - (a == x) for a in range(q)] for x in range(q)]
    return H

def sylv(t):
    M = 2**t
    return [[(-1)**bin(x & a).count('1') for a in range(M)] for x in range(M)]

def hadcheck(H):
    M = len(H)
    for x in range(M):
        assert H[x][0] == 1
        for y in range(M):
            assert sum(H[x][k]*H[y][k] for k in range(M)) == (M if x == y else 0)

def ang(z):
    k = max(abs(z[0]).bit_length(), abs(z[1]).bit_length()) - 900
    if k > 0:
        return math.atan2(z[1] >> k if z[1] >= 0 else -((-z[1]) >> k), z[0] >> k if z[0] >= 0 else -((-z[0]) >> k))
    return math.atan2(z[1], z[0])

def run(H, rng, ntr):
    M = len(H); hadcheck(H)
    splits = [p for p in range(5, 4000) if p % 4 == 1 and is_prime(p)]
    nid = nang = nnz = 0
    for tr in range(ntr):
        ps = rng.sample(splits, 2*(M-1) + 4)
        cols = []
        for i, p in enumerate(ps):
            a = (i % (M-1)) + 1 if i < 2*(M-1) else 0
            # modification sets: mostly small, sometimes large
            k = rng.choice([0, 0, 1, 1, 2, 3]) if rng.random() < 0.9 else rng.randint(0, M//2)
            F = set(rng.sample(range(M), k))
            cols.append(dict(a=a, s=rng.choice([1, -1]), pi=gp(p), p=p, F=F))
        g = prod([gp(13), gp(17), (3, 0)])        # nontrivial content (split and inert parts)
        units = [rng.randrange(4) for _ in range(M)]
        U = [(1,0),(0,1),(-1,0),(0,-1)]
        def sgn(x, c):
            v = c['s'] * H[x][c['a']]
            return -v if x in c['F'] else v
        Z = []
        for x in range(M):
            z = mul(g, U[units[x]])
            for c in cols: z = mul(z, pw(c['pi'], sgn(x, c)))
            Z.append(z)
        N0 = nrm(Z[0]); assert all(nrm(z) == N0 for z in Z)
        # blocks
        def c_of(c, a): return sum(H[x][a]*H[x][c['a']] for x in c['F'])
        G = {a: prod([pw(c['pi'], c['s']) for c in cols if c['a'] == a]) for a in range(1, M)}
        P = {a: prod([pw(c['pi'], c['s']*c_of(c, a)) for c in cols]) for a in range(1, M)}
        # angles
        phi = [ang(c['pi']) for c in cols]
        Sphi = [sum(sgn(x, c)*phi[j] for j, c in enumerate(cols)) for x in range(M)]
        theta = Sphi[0]
        kx = [round((Sphi[x]-theta)/(math.pi/2)) for x in range(M)]
        dx = [Sphi[x]-theta-(math.pi/2)*kx[x] for x in range(M)]
        for a in range(1, M):
            A = prod([Z[x] if H[x][a] == 1 else cj(Z[x]) for x in range(M)])
            Gam = mul(pw(G[a], M), pw(cj(P[a]), 2))
            t = mul(A, cj(Gam)); assert t[0] == 0 or t[1] == 0; nid += 1
            # angle bookkeeping with principal logs
            Lu = 2j*ang(G[a]); Lu = complex(0, math.remainder(Lu.imag, 2*math.pi))
            Lv = 2j*ang(P[a]) if P[a] != (1, 0) else 0j; Lv = complex(0, math.remainder(Lv.imag, 2*math.pi))
            eta = sum(H[x][a]*dx[x] for x in range(M))
            r = (M*Lu - 2*Lv - 2j*eta).imag/(math.pi/2)
            assert abs(r - round(r)) < 1e-6, r; nang += 1
            # nonvanishing
            if any(c['a'] == a and len(c['F']) < M/2 for c in cols):
                ex = {}
                core = [(c['pi'], c['s']*((M if c['a'] == a else 0) - 2*c_of(c, a))) for c in cols]
                gam = prod([pw(pi, e) for pi, e in core if e != 0])
                assert any(e != 0 for _, e in core)
                assert gam[0] != 0 and gam[1] != 0 and abs(gam[0]) != abs(gam[1]); nnz += 1
                # Gamma_a = positive integer * gam (up to unit): check Gamma_a / gam is rational integer * unit
                Gm = Gam; n2 = nrm(gam)
                q = mul(Gm, cj(gam)); assert q[0] % n2 == 0 and q[1] % n2 == 0
                q = (q[0]//n2, q[1]//n2); assert q[0] == 0 or q[1] == 0
        for _ in range(8):
            a, b = rng.sample(range(1, M), 2)
            d = [(c_of(c, a)-c_of(c, b)) for c in cols]; assert all(v % 2 == 0 for v in d)
            Q = prod([pw(c['pi'], c['s']*d[j]//2) for j, c in enumerate(cols)])
            Uab = mul(G[a], cj(G[b]))
            Aa = prod([Z[x] if H[x][a] == 1 else cj(Z[x]) for x in range(M)])
            Ab = prod([Z[x] if H[x][b] == 1 else cj(Z[x]) for x in range(M)])
            t = mul(mul(Aa, cj(Ab)), cj(mul(pw(Uab, M), pw(cj(Q), 4)))); assert t[0] == 0 or t[1] == 0; nid += 1
            Lu = complex(0, math.remainder(2*ang(Uab), 2*math.pi))
            Lw = complex(0, math.remainder(2*ang(Q), 2*math.pi)) if Q != (1, 0) else 0j
            eta = sum((H[x][a]-H[x][b])*dx[x] for x in range(M))
            r = (M*Lu - 4*Lw - 2j*eta).imag/(math.pi/2)
            assert abs(r - round(r)) < 1e-6; nang += 1
            if any(c['a'] in (a, b) and len(c['F']) < M/4 for c in cols):
                core = [(c['pi'], c['s']*((M if c['a'] == a else -M if c['a'] == b else 0) - 4*(d[j]//2))) for j, c in enumerate(cols)]
                gam = prod([pw(pi, e) for pi, e in core if e != 0])
                assert gam[0] != 0 and gam[1] != 0 and abs(gam[0]) != abs(gam[1]); nnz += 1
    return nid, nang, nnz

def minpoly_check(rng, n=3000):
    splits = [p for p in range(5, 400) if p % 4 == 1 and is_prime(p)]
    for _ in range(n):
        ps = rng.sample(splits, rng.randint(1, 4))
        gam = prod([pw(gp(p), rng.choice([1, -1])*rng.randint(1, 3)) for p in ps])
        N = nrm(gam); re2 = mul(gam, gam)[0]
        assert math.gcd(math.gcd(N, 2*re2), N) == 1 and abs(2*re2) <= 2*N
        # u = gam/conj gam is a root of N X^2 - 2 re2 X + N
        u = complex(*gam)/complex(*cj(gam))
        assert abs(N*u*u - 2*re2*u + N) < 1e-6*N
    return n

if __name__ == '__main__':
    rng = random.Random(77)
    tot = [0, 0, 0]
    for name, H, ntr in [('Walsh-8', sylv(3), 30), ('Paley-12', paley(11), 30), ('Walsh-16', sylv(4), 12), ('Paley-20', paley(19), 8)]:
        r = run(H, rng, ntr)
        print(name, 'identities', r[0], 'angle-bookkeeping', r[1], 'nonvanishing', r[2])
        tot = [u+v for u, v in zip(tot, r)]
    print('minpoly checks', minpoly_check(rng))
    print('TOTAL', tot, 'PASS')
