"""Theorem L1 (classification of class-beta_7 lifts of E_t), exact over Q.

Base points r1, r2, r3 = the non-contact points over y_4 with S = {a}, {b}, {c};
L_q = 3H - r1 - r2 - r3 - q; q* := residual point of the plane <r1,r2,r3> on E, so L_{q*} = 2H = T^*O(1).
G(q) = [tau_j(B_{4,k})] (4x3), entries cubic in q; m_i = 3x3 minors (degree-9 forms on E, 36 zeros).
Forced zeros (identically in t):  r_i: order >= 3 (all columns vanish);  B_{j,k}, j<=3: order >= 1 (column j);
B_{4,k}: order >= 1 for the three minors containing row k;  q*: order >= 2 (rank <= 1 there: the pencil).
rho_i = N(m_i) / prod over forced points (a - a(P))^ord.   gcd(rho_1..rho_4) = 1 and no common zero at
infinity  ==> the only q with rank <= 2 are special points or q*;  then direct nullities at the 28 special
points (0) and at q* (2) finish the proof at this t; the same data give a Zariski-open set of t (lift.md 4.3).

Usage: python3 beta7lift2.py x1 x2 x3 x4
"""
import sys, itertools
from fractions import Fraction as Fr
import flint
from etdata import etdata
from beta7lift import MONS, emon_h

Q = flint.fmpq
CHOICE = 'aux'

def fq(x):
    x = Fr(x); return Q(x.numerator, x.denominator)

class Ring:
    def __init__(self, f, g):
        self.f, self.g, self.fg = f, g, f * g
    def mul(self, X, Y):
        A1, B1, C1, D1 = X; A2, B2, C2, D2 = Y
        f, g = self.f, self.g
        return (A1 * A2 + f * B1 * B2 + g * C1 * C2 + self.fg * D1 * D2,
                A1 * B2 + B1 * A2 + g * (C1 * D2 + D1 * C2),
                A1 * C2 + C1 * A2 + f * (B1 * D2 + D1 * B2),
                A1 * D2 + D1 * A2 + B1 * C2 + C1 * B2)
    def sub(self, X, Y):
        return tuple(x - y for x, y in zip(X, Y))
    def add(self, X, Y):
        return tuple(x + y for x, y in zip(X, Y))
    def norm(self, X):
        A, B, C, D = X
        Y = self.mul(X, (A, -B, C, -D))
        Z = self.mul(Y, (Y[0], Y[1], -Y[2], -Y[3]))
        assert Z[1] == 0 and Z[2] == 0 and Z[3] == 0
        return Z[0]

def ev_mons(P):
    a, u, w = P
    return [a ** i * u ** e * w ** d for (i, e, d, k) in MONS]

def grad_dot(P, T):
    a, u, w = P
    out = []
    for (i, e, d, k) in MONS:
        da = i * a ** (i - 1) * u ** e * w ** d if i > 0 else 0
        du = e * a ** i * w ** d if e > 0 else 0
        dw = d * a ** i * u ** e if d > 0 else 0
        out.append(da * T[0] + du * T[1] + dw * T[2])
    return out

def rank_Q(rows):
    if not rows:
        return 0
    M = flint.fmpq_mat(len(rows), len(rows[0]), [fq(x) for r in rows for x in r])
    return M.rank()

def residual_point(E, r):
    """4th intersection of the plane through r1,r2,r3 with E (exact)."""
    import sympy as sp
    a, u, w = sp.symbols('a u w')
    M = sp.Matrix([[sp.Rational(str(c)) for c in P] + [1] for P in r])
    ns = M.nullspace()[0]          # plane alpha a + beta u + gamma w + delta = 0
    al, be, ga, de = ns
    kb, s, kc, tp = (sp.Rational(str(E[k])) for k in ('kb', 's', 'kc', 'tp'))
    # eliminate: solve plane for one variable present
    eqs = [u ** 2 - kb * (a ** 2 - s), w ** 2 - kc * (a ** 2 - tp), al * a + be * u + ga * w + de]
    sols = sp.solve(eqs, [a, u, w], dict=True)
    known = {tuple(sp.Rational(str(c)) for c in P) for P in r}
    out = [tuple(sl[v] for v in (a, u, w)) for sl in sols]
    rest = [P for P in out if P not in known]
    assert len(rest) == 1 and all(c.is_rational for c in rest[0]), (out,)
    return tuple(Fr(int(c.p), int(c.q)) for c in rest[0])

def run(xs, verbose=True):
    E = etdata(xs)
    kb, s, kc, tp = E['kb'], E['s'], E['kc'], E['tp']
    X = E['xs']
    B = {j: [P for lab, A, jj, P in E['pts'] if jj == j] for j in (1, 2, 3, 4)}
    # non-contact points over y_4: all 8 sign combinations minus the 4 contacts
    x4 = X[3]
    U4 = abs(B[4][0][1]); W4 = abs(B[4][0][2])
    over4 = [(sa * x4, su * U4, sw * W4) for sa in (1, -1) for su in (1, -1) for sw in (1, -1)]
    labels = {}
    for P in over4:
        aa, bb, cc = E['to_abc'](*P)
        S = frozenset(n for n, val in (('a', aa), ('b', bb), ('c', cc)) if val == x4)
        labels[P] = S
    if CHOICE == 'y4':
        r = [next(P for P in over4 if labels[P] == frozenset({m})) for m in 'abc']
    else:
        r = list(E['aux'])
    special = r + [P for j in B for P in B[j]]
    assert len(set(special)) == 28
    qs = residual_point(E, r)
    qs_special = qs in special
    if qs_special:
        assert qs in B[4], 'q* special but not in B_4'
    # Lagrange coefficients
    x1, x2, x3, x4 = X
    lam = [(x4 - x2) * (x4 - x3) / ((x1 - x2) * (x1 - x3)), (x4 - x1) * (x4 - x3) / ((x2 - x1) * (x2 - x3)),
           (x4 - x1) * (x4 - x2) / ((x3 - x1) * (x3 - x2))]

    def conds(P, order):
        if order <= 0:
            return []
        v0 = ev_mons(P)
        if order == 1:
            return [v0]
        a, u, w = P
        T = (4 * u * w, 4 * kb * a * w, 4 * kc * a * u)
        return [v0, grad_dot(P, T)]

    def direct_nullity(base):
        rows = []
        def place(vec, j):
            rr = [0] * 36
            for n in range(12):
                rr[12 * j + n] = vec[n]
            return rr
        for j in range(3):
            for P, m in base.items():
                for v in conds(P, m):
                    rows.append(place(v, j))
            for P in B[j + 1]:
                for v in conds(P, base.get(P, 0) + 1):
                    rows.append(place(v, j))
        for P in B[4]:
            for v in conds(P, base.get(P, 0) + 1):
                rr = [0] * 36
                for j in range(3):
                    for n in range(12):
                        rr[12 * j + n] += lam[j] * v[n]
                rows.append(rr)
        return 36 - rank_Q(rows)

    # symbolic G(q)
    P = lambda c: flint.fmpq_poly([fq(x) for x in c])
    f = P([-kb * s, 0, kb]); g = P([-kc * tp, 0, kc])
    ring = Ring(f, g)
    Rr = [ev_mons(P_) for P_ in r]
    B4 = [ev_mons(b) for b in B[4]]
    unit = lambda n: [1 if i == n else 0 for i in range(12)]
    G = [[None] * 3 for _ in range(4)]
    for j in range(3):
        fixed = Rr + [ev_mons(b) for b in B[j + 1]]
        K = [[None] * 12 for _ in range(12)]
        for m in range(12):
            for n in range(12):
                rows = fixed + [unit(n), unit(m)]
                K[m][n] = Fr(str(flint.fmpq_mat(12, 12, [fq(x) for rr in rows for x in rr]).det()))
        for k in range(4):
            gam = [sum(B4[k][m] * K[m][n] for m in range(12)) for n in range(12)]
            comp = {c: [Fr(0)] * 4 for c in [(0, 0), (1, 0), (0, 1), (1, 1)]}
            for n, (i, e, d, kk) in enumerate(MONS):
                comp[(e, d)][i] += gam[n]
            G[k][j] = tuple(P(comp[c]) for c in [(0, 0), (1, 0), (0, 1), (1, 1)])
    def det3(rows):
        (a, b, c), (d, e, f_), (g_, h, i) = rows
        t1 = ring.mul(a, ring.sub(ring.mul(e, i), ring.mul(f_, h)))
        t2 = ring.mul(b, ring.sub(ring.mul(d, i), ring.mul(f_, g_)))
        t3 = ring.mul(c, ring.sub(ring.mul(d, h), ring.mul(e, g_)))
        return ring.add(ring.sub(t1, t2), t3)
    minors = {k: det3([G[rr] for rr in range(4) if rr != k]) for k in range(4)}
    rhos = {}
    info = {}
    for k, m in minors.items():
        N = ring.norm(m)
        forced = {}
        def addf(Pt, o):
            forced[Pt[0]] = forced.get(Pt[0], 0) + o
        for P_ in r: addf(P_, 3)
        for j in (1, 2, 3):
            for P_ in B[j]: addf(P_, 1)
        for kk, P_ in enumerate(B[4]):
            if kk != k and P_ != qs: addf(P_, 1)
        addf(qs, 2)
        rem = N
        for av, o in forced.items():
            lin = P([-av, 1]) ** o
            qq, rr = divmod(rem, lin)
            assert rr == 0, ('forced order fails', k, av, o)
            rem = qq
        rhos[k] = rem
        info[k] = dict(normdeg=N.degree(), forced=sum(forced.values()), rhodeg=rem.degree(),
                       bin_deg=36 - sum(forced.values()))
    gg = rhos[0]
    for k in (1, 2, 3):
        gg = gg.gcd(rhos[k])
    inf_common = all(info[k]['rhodeg'] < info[k]['bin_deg'] for k in range(4))
    # direct nullities
    sp_null = []
    for qq_ in special:
        base = {P_: 1 for P_ in r}
        base[qq_] = base.get(qq_, 0) + 1
        sp_null.append(direct_nullity(base))
    base = {P_: 1 for P_ in r}; base[qs] = base.get(qs, 0) + 1
    qs_null = direct_nullity(base)
    res = dict(xs=[str(x) for x in X], choice=CHOICE, qstar=[str(c) for c in qs], qstar_special=qs_special, minor_info=info,
               gcd_rho=str(gg), gcd_deg=gg.degree(), common_zero_at_infinity=inf_common,
               special_nullities=sorted(sp_null), qstar_nullity=qs_null)
    if verbose:
        for k_, v in res.items():
            print(k_, ':', v)
    return res

if __name__ == '__main__':
    xs = [Fr(v) for v in sys.argv[1:5]] if len(sys.argv) >= 5 else [Fr(2), Fr(3), Fr(-5), Fr(1)]
    if len(sys.argv) >= 6: CHOICE = sys.argv[5]
    run(xs)
