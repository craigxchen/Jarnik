# Exact verification of Theorems A/B (G(2,4)-model, family V_A) for the arrangement in arr.py.
import sys, random, itertools
import sympy as sp
from sympy import Matrix, Rational as R, Poly, QQ
from star import point_of, line_through_x, quadric_through, member_through, analyze, kernel, LAB, t
from arr import PLANES, CHOSEN
ok = True
# (1) arrangement checks
for S in itertools.combinations(range(7), 6):
    if Matrix.vstack(*[PLANES[a] for a in S]).rank() < 4: print('point on 6 planes', S); ok = False
coax = [T for T in itertools.combinations(range(7), 3) if Matrix.vstack(*[PLANES[a] for a in T]).rank() == 2]
print('coaxial triples:', [[LAB[a] for a in T] for T in coax])
five = [S for S in itertools.combinations(range(7), 5) if Matrix.vstack(*[PLANES[a] for a in S]).rank() == 3]
print('points on 5 planes:', [[LAB[a] for a in S] for S in five])
pts = [point_of([PLANES[i], PLANES[j], PLANES[k]]) for (i, j, k) in CHOSEN]
print('event points:', [list(p.T) for p in pts])
Lam = quadric_through(pts); print('dim of quadrics through the 6 points:', len(Lam)); ok &= len(Lam) == 4
# (2) members through random rational points
random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 7)
def rnd(): return R(random.randint(-40, 40), random.randint(1, 15))
def cr(a, b, c, d): return (a - c)*(b - d)/((a - d)*(b - c))
for trial in range(int(sys.argv[2]) if len(sys.argv) > 2 else 3):
    x = []
    while len(x) < 7:
        c = rnd()
        if c not in x: x.append(c)
    data, Q, (L1, L2) = member_through(PLANES, pts, x)
    okQ = Q.rank() == 4
    # passes through x at t=0: compare cross-ratios
    v0 = [sp.Rational(sp.sympify(N).subs(t, 0)) / sp.Rational(sp.sympify(D).subs(t, 0)) for (N, D) in data]
    okx = all(cr(v0[0], v0[1], v0[2], v0[a]) == cr(x[0], x[1], x[2], x[a]) for a in range(3, 7))
    ev, inf, P = analyze(data)
    prod = Poly(1, t, domain=QQ)
    for e in ev.values(): prod = prod*e
    sqf = prod.sqf_part().degree() == prod.degree()
    four = [T for T in ev if len(T) == 4]
    pair = sum(e.degree() for T, e in ev.items() if len(T) == 2)
    trip = sum(e.degree() for T, e in ev.items() if len(T) == 3)
    cls = {tuple(LAB[a] for a in T): e.degree() for T, e in ev.items()}
    good = okQ and okx and sqf and not inf and not four and pair == 6 and trip == 12
    ok &= good
    print('x =', x, '| Q smooth', okQ, '| through x', okx, '| contacts distinct', sqf, '| none at inf', not inf,
          '| 4-collisions', four, '| pair', pair, 'triple', trip)
    print('   class:', cls)
print('ALL CHECKS PASSED' if ok else 'SOME CHECK FAILED')
