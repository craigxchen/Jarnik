# E_Lambda . C = degree of t -> Q_{x(t)} (quadric through R_{x(t)} and the 6 chosen points), exact over Q.
import sys, random, itertools
import sympy as sp
from sympy import Matrix, Poly, QQ, symbols, Rational as R
from star import star_arrangement, point_of, quadric_through, LAB
t = symbols('t')
def basis_quadrics(pts):
    return quadric_through(pts)  # list of 4 symmetric matrices
def deg_EL(planes, qbasis, mob):
    """mob: list of 7 pairs (N_a, D_a) linear polys in t. Returns degree of t -> Q(t) in P^3."""
    rows = []
    for a, l in enumerate(planes):
        Na, Da = mob[a]
        rows.append([sp.expand(Na*l[i]) for i in range(4)] + [sp.expand(Da*l[i]) for i in range(4)])
    M = Matrix(rows)
    ker = []
    for j in range(8):
        sub = M[:, [k for k in range(8) if k != j]]
        ker.append(sp.expand((-1)**j * sub.det(method='berkowitz')))
    g = 0
    for k in ker: g = sp.gcd(g, k)
    ker = [sp.cancel(k/g) for k in ker]
    L1 = Matrix(ker[:4]); L2 = Matrix(ker[4:])
    conds = []
    for (u, v) in [(L1, L1), (L1, L2), (L2, L2)]:
        conds.append([sp.expand((u.T*Q*v)[0]) for Q in qbasis])
    C = Matrix(conds)  # 3 x 4
    minors = []
    for j in range(4):
        sub = C[:, [k for k in range(4) if k != j]]
        minors.append(sp.expand((-1)**j * sub.det(method='berkowitz')))
    g = 0
    for m in minors: g = sp.gcd(g, m)
    minors = [sp.cancel(m/g) for m in minors]
    return max(Poly(m, t).degree() if m != 0 else -1 for m in minors), [Poly(L1[i], t).degree() for i in range(4)]
def random_mobius(rnd, event=None):
    mob = []
    for a in range(7):
        mob.append([rnd(), rnd(), rnd(), rnd()])  # (al, be, ga, de): (al t + be)/(ga t + de)
    if event is not None:
        a0 = event[0]; z0 = mob[a0][1]/mob[a0][3]  # value at t=0
        for b in event[1:]:
            mob[b][1] = z0*mob[b][3]
    return [(c[0]*t + c[1], c[2]*t + c[3]) for c in mob]
if __name__ == '__main__':
    random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 3)
    rnd = lambda: R(random.randint(-9, 9), random.randint(1, 4))
    planes = star_arrangement(rnd)
    chosen = [(1,3,5),(1,3,6),(1,4,5),(1,4,6),(2,3,5),(2,4,5)]
    pts = [point_of([planes[i], planes[j], planes[k]]) for (i,j,k) in chosen]
    qb = basis_quadrics(pts); assert len(qb) == 4
    d, dl = deg_EL(planes, qb, random_mobius(rnd))
    print('pure Mobius mu: E_Lambda deg', d, 'L1 degs', dl, flush=True)
    tot = 0
    for T in itertools.combinations(range(7), 3):
        d3, _ = deg_EL(planes, qb, random_mobius(rnd, T))
        tot += d3
        print([LAB[a] for a in T], d3, flush=True)
    print('sum over triples', tot, ' E.beta = sum - 32*E.mu =', tot - 32*d)
