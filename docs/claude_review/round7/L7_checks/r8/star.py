# Exact construction: star arrangement A of 7 planes in P^3 over Q, member of the conic family through x.
import sympy as sp, itertools, random
from sympy import Rational as R, Matrix, symbols, Poly, gcd, QQ
t = symbols('t')
LAB = ['O','a1','a2','b1','b2','c1','c2']
def kernel(M):
    ns = Matrix(M).nullspace()
    return ns
def star_arrangement(rnd):
    # H_O : X3 = 0 ; lines l_i = H_O cap {m_i = 0}; planes through l_i: m_i + alpha * X3
    lO = Matrix([[0,0,0,1]])
    ms = [Matrix([[rnd(), rnd(), rnd(), 0]]) for i in range(3)]
    planes = [lO]
    for i in range(3):
        for j in range(2):
            alpha = rnd()
            planes.append(ms[i] + alpha*lO)
    return planes  # order: O, a1, a2, b1, b2, c1, c2
def point_of(planes_sub):
    ns = kernel(Matrix.vstack(*planes_sub))
    assert len(ns) == 1, 'not a point'
    return ns[0]
def line_through_x(planes, x):
    # unknowns L1 (4), L2 (4): l_a(x_a L1 + L2) = 0
    rows = []
    for a, l in enumerate(planes):
        rows.append(list(x[a]*l) + list(l))
    ns = kernel(Matrix(rows))
    return ns
def quadric_through(points):
    # symmetric 4x4: unknowns q00,q01,q02,q03,q11,q12,q13,q22,q23,q33
    idx = [(i,j) for i in range(4) for j in range(i,4)]
    rows = []
    for p in points:
        rows.append([p[i]*p[j]*(1 if i==j else 2) for (i,j) in idx])
    ns = kernel(Matrix(rows))
    out = []
    for v in ns:
        Q = sp.zeros(4,4)
        for k,(i,j) in enumerate(idx):
            Q[i,j] = v[k]; Q[j,i] = v[k]
        out.append(Q)
    return out
def other_line(Q, p, Rpts):
    # T_p Q cap Q = R_x cup R'; find R' (two points spanning it) rationally
    tp = (p.T*Q)  # 1x4
    basis = kernel(tp)  # 3 vectors spanning plane T_pQ (as subspace)
    B = Matrix.hstack(*basis)  # 4x3
    q3 = B.T*Q*B  # ternary form, rank <= 2
    u = symbols('u0:3')
    form = sp.expand((Matrix([u]).T if False else Matrix(u)).T*q3*Matrix(u))[0]
    fl = sp.factor_list(form)
    lins = []
    for f, e in fl[1]:
        if sp.Poly(f, *u).total_degree() == 1:
            lins += [f]*e
    assert len(lins) == 2, ('factor', fl)
    # each linear form defines a line in plane; points of line = kernel of coefficient row
    lines = []
    for f in lins:
        c = Matrix([[sp.Poly(f, *u).coeff_monomial(v) for v in u]])
        ker = kernel(c)  # 2 vectors in u-space
        pts = [B*k for k in ker]
        lines.append(pts)
    # identify which is R_x: contains the Rpts
    def on_line(pts, q):
        return Matrix.hstack(*(pts+[q])).rank() == 2
    rest = [L for L in lines if not all(on_line(L, q) for q in Rpts)]
    assert len(rest) == 1
    return rest[0]

def member_through(planes, events_pts, x):
    """Return Mobius data [(N_a, D_a)] (linear polys in t) of the member through x, plus Q."""
    ns = line_through_x(planes, x)
    assert len(ns) == 1, ('line system rank deficient', len(ns))
    v = ns[0]; L1 = v[:4, :]; L2 = v[4:, :]
    assert Matrix.hstack(L1, L2).rank() == 2
    Qs = quadric_through([L1, L2, L1+L2] + events_pts)
    assert len(Qs) == 1, ('quadric not unique', len(Qs))
    Q = Qs[0]
    p, pp = L1, L2
    Rp = other_line(Q, p, [L1, L2]); Rpp = other_line(Q, pp, [L1, L2])
    # w on R'_p not p, w'' on R'_pp not pp
    w = [q for q in Rp if Matrix.hstack(q, p).rank() == 2][0]
    w2 = [q for q in Rpp if Matrix.hstack(q, pp).rank() == 2][0]
    q_t = p + t*w
    den = (q_t.T*Q*w2)[0]; num = (w.T*Q*pp)[0]
    q2_t = den*pp - t*num*w2
    data = []
    for l in planes:
        Na = sp.expand(-(l*q2_t)[0]); Da = sp.expand((l*q_t)[0])
        data.append((Na, Da))
    return data, Q, (L1, L2)
def analyze(data, verbose=True):
    """Exact contact analysis of t -> (N_a/D_a). Returns dict cluster -> event polynomial."""
    n = len(data)
    P = {}
    for a, b in itertools.combinations(range(n), 2):
        Na, Da = data[a]; Nb, Db = data[b]
        P[(a,b)] = Poly(sp.expand(Na*Db - Nb*Da), t, domain=QQ)
    # contacts at t = infinity: degree drop
    infinity = [ab for ab in P if P[ab].degree() < 2]
    ev = {}
    # triple gcds
    for T in itertools.combinations(range(n), 3):
        g = P[(T[0],T[1])]
        for ab in itertools.combinations(T, 2): g = g.gcd(P[ab])
        if g.degree() > 0: ev[T] = g.monic()
    for T in itertools.combinations(range(n), 4):
        g = P[(T[0],T[1])]
        for ab in itertools.combinations(T, 2): g = g.gcd(P[ab])
        if g.degree() > 0: ev[T] = g.monic()
    # pure pair part
    for ab in P:
        g = P[ab]
        for T, e in ev.items():
            if set(ab) <= set(T) and len(T) == 3:
                q, r = g.div(e); assert r.is_zero; g = q
        if g.degree() > 0: ev[ab] = g.monic()
    return ev, infinity, P
