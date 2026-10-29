# Judge check J2 (own code): killer's Proposition K, steps 4-6.
# Five lines L_a: w = rho_a*y + kappa_a in the (y,w)-plane, generic (distinct slopes, no 3 concurrent, none through O).
# Claim: (i) quartics through the 10 nodes form a 5-dim space spanned by the products of 4 of the lines;
#        (ii) no nonzero quartic through the 10 nodes has order >= 3 at O;
#        (iii) if 4 of the lines are concurrent, such quartics exist (the proof's only escape).
# Also an end-to-end check of steps 1-4: random q, Psi, rho, kappa -> the curve phi_a = Psi + kappa_a q/(s-rho_a);
# r_ab = common roots of N_a, N_b exist only if G vanishes at the nodes.
import random, itertools, sys
import sympy as sp
from fractions import Fraction as F
random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 3)
y, w = sp.symbols('y w')
monos = [y**i * w**j for i in range(5) for j in range(5) if i + j <= 4]   # all monomials of degree <= 4 (affine quartics)
low = [m for m in monos if sp.Poly(m, y, w).total_degree() <= 2]          # order >= 3 at O kills these coefficients

def nodes(rho, kap):
    P = []
    for a, b in itertools.combinations(range(5), 2):
        yy = -(kap[a] - kap[b]) / (rho[a] - rho[b])
        P.append((yy, rho[a] * yy + kap[a]))
    return P

def space_dims(rho, kap):
    P = nodes(rho, kap)
    rows = [[m.subs({y: p[0], w: p[1]}) for m in monos] for p in P]
    A = sp.Matrix(rows)
    d_nodes = len(monos) - A.rank()
    B = sp.Matrix(rows + [[1 if m == l else 0 for m in monos] for l in low])
    d_trip = len(monos) - B.rank()
    # products of 4 lines lie in the node space and are independent
    lines = [w - rho[a] * y - kap[a] for a in range(5)]
    prods = [sp.expand(sp.prod([lines[b] for b in range(5) if b != c])) for c in range(5)]
    C = sp.Matrix([[sp.Poly(pr, y, w).coeff_monomial(m) for m in monos] for pr in prods])
    ok_in = all(sp.simplify(pr.subs({y: p[0], w: p[1]})) == 0 for pr in prods for p in P)
    return d_nodes, d_trip, C.rank(), ok_in

def rnd():
    return sp.Rational(random.randint(-30, 30), random.randint(1, 9))

res = []
for inst in range(8):
    while True:
        rho = [rnd() for _ in range(5)]; kap = [rnd() for _ in range(5)]
        if len(set(rho)) < 5 or 0 in kap: continue
        P = nodes(rho, kap)
        if len(set(P)) < 10: continue   # no three concurrent
        break
    res.append(('generic', space_dims(rho, kap)))
    print(res[-1], flush=True)
# four concurrent lines through a point (y0,w0) != O, fifth generic
for inst in range(3):
    y0, w0 = rnd(), rnd()
    rho = [rnd() for _ in range(5)]
    kap = [w0 - rho[a] * y0 for a in range(4)] + [rnd()]
    if len(set(rho)) < 5 or 0 in kap: continue
    P = nodes(rho, kap)
    rows = [[m.subs({y: p[0], w: p[1]}) for m in monos] for p in P]
    B = sp.Matrix(rows + [[1 if m == l else 0 for m in monos] for l in low])
    print(('4-concurrent', 'triple-point space dim', len(monos) - B.rank()), flush=True)

# end-to-end: build phi_a from random q (monic quartic), Psi (cubic), rho, kappa; check resultants of N_a,N_b
s = sp.symbols('s')
for inst in range(3):
    q = sp.expand(sp.prod([s - rnd() for _ in range(4)]))
    Psi = sum(rnd() * s**i for i in range(4))
    rho = [rnd() for _ in range(5)]; kap = [rnd() for _ in range(5)]
    N = [sp.expand(Psi * (s - rho[a]) + kap[a] * q) for a in range(5)]
    nres = sum(1 for a, b in itertools.combinations(range(5), 2) if sp.resultant(N[a], N[b], s) == 0)
    print(('end-to-end random data: #pairs (a,b) with common root of N_a,N_b (need 10 for the class)', nres), flush=True)
