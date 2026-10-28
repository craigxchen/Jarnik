"""Theorem L1 check: no class-beta_7 curve over E_t (any genus).

Claim checked (over F_p, which implies it over Qbar by properness, see lift.md Sec. 4):
for EVERY line bundle L of degree 8 on E_t, no nonzero sigma in H^0(E_t x P^1, L (x) O(2)) vanishes at
the 25 points (B_A, v_A).  Parametrisation: L_q = 3H - r1 - r2 - r3 - q, q in E (bijection E -> Pic^8);
H^0(L_q) = cubic forms (normal monomials mod u^2, w^2) vanishing at r1, r2, r3, q.
sigma(., v) = sum_j l_j(v) sigma_j (Lagrange at x1, x2, x3); sigma_j in H^0(L(-sum B_j)) = <tau_j>;
nonzero sigma exists  <=>  rank [tau_j(B_{4,k})]_{k<=4, j<=3} <= 2.

Usage: python3 beta7lift.py p x1 x2 x3 x4
"""
import sys
from fractions import Fraction as Fr
from math import lcm
import flint
from etdata import etdata

MONS = [(i, e, d, 3 - e - d - i) for (e, d) in [(0, 0), (1, 0), (0, 1), (1, 1)] for i in range(3 - e - d + 1)]
assert len(MONS) == 12

def mp(x, p):
    x = Fr(x)
    return (x.numerator % p) * pow(x.denominator % p, -1, p) % p

def emon_h(P, one):
    """monomials at homogeneous point P = (a,u,w,h) (field elements or ints)."""
    a, u, w, h = P
    return [a ** i * u ** e * w ** d * h ** k for (i, e, d, k) in MONS]

def grad_dot(P, T, p):
    """d/dtau of each monomial at affine P (h=1) along tangent T; ints mod p."""
    a, u, w = P
    out = []
    for (i, e, d, k) in MONS:
        da = i * pow(a, i - 1, p) * pow(u, e, p) * pow(w, d, p) if i > 0 else 0
        du = e * pow(a, i, p) * pow(w, d, p) if e > 0 else 0
        dw = d * pow(a, i, p) * pow(u, e, p) if d > 0 else 0
        out.append((da * T[0] + du * T[1] + dw * T[2]) % p)
    return out

class Ring:
    """F_p[a][u,w]/(u^2 - f, w^2 - g); elements are 4-tuples (A,B,C,D) = A + B u + C w + D uw."""
    def __init__(self, f, g, p):
        self.f, self.g, self.p = f, g, p
        self.fg = f * g
    def mul(self, X, Y):
        A1, B1, C1, D1 = X; A2, B2, C2, D2 = Y
        f, g = self.f, self.g
        return (A1 * A2 + f * B1 * B2 + g * C1 * C2 + self.fg * D1 * D2,
                A1 * B2 + B1 * A2 + g * (C1 * D2 + D1 * C2),
                A1 * C2 + C1 * A2 + f * (B1 * D2 + D1 * B2),
                A1 * D2 + D1 * A2 + B1 * C2 + C1 * B2)
    def add(self, X, Y):
        return tuple(x + y for x, y in zip(X, Y))
    def neg(self, X):
        return tuple(-x for x in X)
    def norm(self, X):
        A, B, C, D = X
        Y = self.mul(X, (A, -B, C, -D))
        Z = self.mul(Y, (Y[0], Y[1], -Y[2], -Y[3]))
        assert Z[1] == 0 and Z[2] == 0 and Z[3] == 0
        return Z[0]

def setup(xs, p):
    E = etdata(xs)
    kb, s, kc, tp = (mp(E[k], p) for k in ('kb', 's', 'kc', 'tp'))
    red = lambda P: tuple(mp(c, p) for c in P)
    B = {j: [red(P) for lab, A, jj, P in E['pts'] if jj == j] for j in (1, 2, 3, 4)}
    R = [red(P) for P in E['aux']]
    X = [mp(x, p) for x in E['xs']]
    allp = [P for j in B for P in B[j]] + R
    assert len(set(allp)) == 28, 'points collide mod p'
    assert len(set(X)) == 4 and all(x != 0 for x in X)
    assert kb and kc and s and tp and (s - tp) % p, 'bad reduction'
    for P in allp:
        a, u, w = P
        assert (u * u - kb * (a * a - s)) % p == 0 and (w * w - kc * (a * a - tp)) % p == 0
    return E, dict(p=p, kb=kb, s=s, kc=kc, tp=tp, B=B, R=R, X=X)

def lagr4(D):
    p = D['p']; x1, x2, x3, x4 = D['X']
    l1 = (x4 - x2) * (x4 - x3) * pow((x1 - x2) * (x1 - x3), -1, p) % p
    l2 = (x4 - x1) * (x4 - x3) * pow((x2 - x1) * (x2 - x3), -1, p) % p
    l3 = (x4 - x1) * (x4 - x2) * pow((x3 - x1) * (x3 - x2), -1, p) % p
    return [l1, l2, l3]

def direct_nullity(D, base):
    """base: dict point(affine tuple or ('inf',hom tuple)) -> multiplicity; returns nullity of the
    36-unknown system for L = 3H - base.  Affine points only get derivative conditions (mult <= 2)."""
    p = D['p']; kb, kc = D['kb'], D['kc']
    lam = lagr4(D)
    rows = []
    def conds(P, order):
        """list of 12-vectors: vanish to order >= order at P."""
        if order <= 0:
            return []
        if P[0] == 'inf':
            assert order == 1
            return [[x % p for x in emon_h(P[1], 1)]]
        a, u, w = P
        v0 = [x % p for x in emon_h((a, u, w, 1), 1)]
        if order == 1:
            return [v0]
        assert order == 2
        T = (4 * u * w % p, 4 * kb * a * w % p, 4 * kc * a * u % p)
        assert any(T)
        return [v0, grad_dot(P, T, p)]
    def place(vec, j, coef=1):
        r = [0] * 36
        for n in range(12):
            r[12 * j + n] = vec[n] * coef % p
        return r
    for j in range(3):
        for P, m in base.items():
            for v in conds(P, m):
                rows.append(place(v, j))
        for P in D['B'][j + 1]:
            for v in conds(P, base.get(P, 0) + 1):
                rows.append(place(v, j))
    for P in D['B'][4]:
        for v in conds(P, base.get(P, 0) + 1):
            r = [0] * 36
            for j in range(3):
                for n in range(12):
                    r[12 * j + n] = (r[12 * j + n] + lam[j] * v[n]) % p
            rows.append(r)
    M = flint.nmod_mat(rows, p)
    return 36 - M.rank()

def symbolic(D):
    p = D['p']
    P = lambda c: flint.nmod_poly(c, p)
    f = P([(-D['kb'] * D['s']) % p, 0, D['kb']]); g = P([(-D['kc'] * D['tp']) % p, 0, D['kc']])
    ring = Ring(f, g, p)
    Rr = [[x % p for x in emon_h(r + (1,), 1)] for r in D['R']]
    B4 = [[x % p for x in emon_h(b + (1,), 1)] for b in D['B'][4]]
    unit = lambda n: [1 if i == n else 0 for i in range(12)]
    G = [[None] * 3 for _ in range(4)]
    for j in range(3):
        fixed = Rr + [[x % p for x in emon_h(b + (1,), 1)] for b in D['B'][j + 1]]
        K = [[int(flint.nmod_mat(fixed + [unit(n), unit(m)], p).det()) for n in range(12)] for m in range(12)]
        for k in range(4):
            gam = [sum(B4[k][m] * K[m][n] for m in range(12)) % p for n in range(12)]
            comp = {(0, 0): [0] * 4, (1, 0): [0] * 4, (0, 1): [0] * 4, (1, 1): [0] * 4}
            for n, (i, e, d, kk) in enumerate(MONS):
                comp[(e, d)][i] = (comp[(e, d)][i] + gam[n]) % p
            G[k][j] = tuple(P(comp[c]) for c in [(0, 0), (1, 0), (0, 1), (1, 1)])
    def det3(rows):
        (a, b, c), (d, e, f_), (g_, h, i) = rows
        t1 = ring.mul(a, ring.add(ring.mul(e, i), ring.neg(ring.mul(f_, h))))
        t2 = ring.mul(b, ring.add(ring.mul(d, i), ring.neg(ring.mul(f_, g_))))
        t3 = ring.mul(c, ring.add(ring.mul(d, h), ring.neg(ring.mul(e, g_))))
        return ring.add(ring.add(t1, ring.neg(t2)), t3)
    minors = [det3([G[r] for r in range(4) if r != k]) for k in range(4)]
    norms = [ring.norm(m) for m in minors]
    return G, minors, norms, ring

def eval_G_at(D, G, pt, F):
    """evaluate the 4x3 matrix at a point (a,u,w) with coordinates in field F (fq_default)."""
    a, u, w = pt
    M = [[None] * 3 for _ in range(4)]
    for k in range(4):
        for j in range(3):
            A, B, C, Dd = G[k][j]
            ev = lambda poly: sum((F(int(c)) * a ** i for i, c in enumerate(poly.coeffs())), F(0))
            M[k][j] = ev(A) + ev(B) * u + ev(C) * w + ev(Dd) * u * w
    return M

def rank_le2(M, F):
    import itertools
    for rows in itertools.combinations(range(4), 3):
        (a, b, c), (d, e, f), (g, h, i) = [M[r] for r in rows]
        det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
        if not det.is_zero():
            return False
    return True

def run(xs, p, verbose=True):
    E, D = setup(xs, p)
    G, minors, norms, ring = symbolic(D)
    for N in norms:
        assert not N.is_zero()
    g = norms[0]
    for N in norms[1:]:
        g = g.gcd(N)
    fac = g.factor()[1] if g.degree() > 0 else []
    special = set(D['R']) | set(P for j in D['B'] for P in D['B'][j])
    special_a = {P[0] for P in special}
    genuine = []
    examined = 0
    for phi, mult in fac:
        r = phi.degree()
        F = flint.fq_default_ctx(p, 2 * r)
        PF = flint.fq_default_poly_ctx(F)
        phiF = PF([F(int(c)) for c in phi.coeffs()])
        for a0, _ in phiF.roots():
            fa = F(int(ring.f.coeffs()[0])) + F(int(ring.f.coeffs()[2])) * a0 ** 2
            ga = F(int(ring.g.coeffs()[0])) + F(int(ring.g.coeffs()[2])) * a0 ** 2
            us = [x for x, _ in PF([-fa, F(0), F(1)]).roots()] if not fa.is_zero() else [F(0)]
            ws = [x for x, _ in PF([-ga, F(0), F(1)]).roots()] if not ga.is_zero() else [F(0)]
            for u0 in us:
                for w0 in ws:
                    examined += 1
                    M = eval_G_at(D, G, (a0, u0, w0), F)
                    if rank_le2(M, F):
                        # spurious?  only if the point is one of the 28 rational special points
                        is_sp = False
                        if r == 1:
                            for P in special:
                                if F(P[0]) == a0 and F(P[1]) == u0 and F(P[2]) == w0:
                                    is_sp = True
                        if not is_sp:
                            genuine.append((str(phi), str(a0), str(u0), str(w0)))
    # spurious points and points at infinity: direct computation
    sp_null = {}
    for q in special:
        base = {r: 1 for r in D['R']}
        base[q] = base.get(q, 0) + 1
        sp_null[q] = direct_nullity(D, base)
    inf_null = []
    kb, kc = D['kb'], D['kc']
    sb = [x for x in range(p) if x * x % p == kb]; sc = [x for x in range(p) if x * x % p == kc]
    if sb and sc:
        for ub in sb:
            for wc in sc:
                base = {r: 1 for r in D['R']}
                base[('inf', (1, ub, wc, 0))] = 1
                inf_null.append(direct_nullity(D, base))
    else:
        inf_null = None
    # cross-check symbolic vs direct at random rational non-special affine points
    import random
    rnd = random.Random(1)
    checks = 0
    while checks < 6:
        a0 = rnd.randrange(p)
        fa = (D['kb'] * (a0 * a0 - D['s'])) % p; ga = (D['kc'] * (a0 * a0 - D['tp'])) % p
        su = [x for x in range(p) if x * x % p == fa] if p < 5000 else []
        sw = [x for x in range(p) if x * x % p == ga] if p < 5000 else []
        if not su or not sw:
            if p >= 5000:
                break
            continue
        q = (a0, su[0], sw[0])
        if q in special:
            continue
        base = {r: 1 for r in D['R']}; base[q] = 1
        Fq = flint.fq_default_ctx(p, 1)
        M = eval_G_at(D, G, tuple(Fq(c) for c in q), Fq)
        assert (direct_nullity(D, base) > 0) == rank_le2(M, Fq)
        checks += 1
    res = dict(p=p, gcd_deg=g.degree(), factors=[(str(f_), m_) for f_, m_ in fac], examined=examined,
               genuine=genuine, spurious_nullities=sorted(sp_null.values()), inf_nullities=inf_null,
               norm_degs=[N.degree() for N in norms], crosschecks=checks)
    if verbose:
        for k, v in res.items():
            print(k, ':', v)
    return res

if __name__ == '__main__':
    p = int(sys.argv[1]) if len(sys.argv) > 1 else 10007
    xs = [Fr(v) for v in sys.argv[2:6]] if len(sys.argv) >= 6 else [Fr(2), Fr(3), Fr(-5), Fr(1)]
    run(xs, p)
