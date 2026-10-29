# Killer: net of sections of -3K with order >= 4 at P, and members with square tangent cone (own code).
import sympy as sp
x, y, z, e, s = sp.symbols('x y z e s')
f = x**3*y + 2*x**2*y*z - 2*x*y**3 - 2*x*y*z**2 + 2*x*z**3 + 2*y**4 + 2*y**3*z + y**2*z**2 - y*z**3 + z**4
g = sp.expand(f.subs({z: 1}).subs({x: e*x, y: e*y}))
wser = sp.expand(sp.series(sp.sqrt(g), e, 0, 7).removeO())   # w(ex,ey) up to e^6
# basis of H^0(-3K): cubic monomials in x,y,z, and s*{x,y,z}
cub = [x**i*y**j*z**(3-i-j) for i in range(4) for j in range(4-i)]
basis = cub + [s*x, s*y, s*z]
cs = sp.symbols('c0:%d' % len(basis))
G = sum(c*m for c, m in zip(cs, basis))
# local expansion at P: z=1, x->e x, y->e y, s->w
loc = sp.expand(G.subs({s: sp.Symbol('W')}).subs({z: 1}).subs({x: e*x, y: e*y}).subs({sp.Symbol('W'): wser}))
eqs = []
for d in range(4):
    cd = sp.expand(loc.coeff(e, d))
    eqs += sp.Poly(cd, x, y).coeffs() if cd != 0 else []
sol = sp.linsolve(eqs, cs)
sol = list(sol)[0]
free = sorted(set().union(*[sp.sympify(v).free_symbols for v in sol]) & set(cs), key=str)
print("dimension of sections of -3K with ord_P >= 4:", len(free))
Gen = sp.expand(sum(v*m for v, m in zip(sol, basis)))
gens = [sp.expand(Gen.coeff(fr)) if True else None for fr in free]
gens = [sp.expand(sp.diff(Gen, fr)) for fr in free]
# tangent cones (degree-4 parts)
TC = []
for gg in gens:
    lg = sp.expand(gg.subs({s: sp.Symbol('W')}).subs({z: 1}).subs({x: e*x, y: e*y}).subs({sp.Symbol('W'): wser}))
    for d in range(4):
        assert sp.expand(lg.coeff(e, d)) == 0
    TC.append(sp.expand(lg.coeff(e, 4)))
print("tangent cones:"); [print("  ", t) for t in TC]
# span of tangent cones: 5-vectors
vec = lambda q4: [sp.Poly(q4, x, y).coeff_monomial(x**(4-i)*y**i) for i in range(5)]
Mt = sp.Matrix([vec(t) for t in TC])
print("rank of tangent-cone span:", Mt.rank())
# annihilators of the span
ann = Mt.nullspace()
print("number of linear conditions for q^2 in span:", len(ann))
al, be, ga = sp.symbols('al be ga')
q2 = sp.expand((al*x**2 + be*x*y + ga*y**2)**2)
v2 = sp.Matrix(vec(q2))
conics = [sp.expand((a.T*v2)[0]) for a in ann]
conics = [sp.Poly(sp.together(c)*sp.denom(sp.together(c)), al, be, ga).as_expr() for c in conics]
print("conics:", conics)
# intersect the two conics in P^2: chart ga=1, and ga=0
r = sp.resultant(conics[0].subs(ga, 1), conics[1].subs(ga, 1), be)
print("eliminant in al (chart ga=1):", sp.factor_list(r))
r2 = sp.resultant(conics[0].subs(al, 1), conics[1].subs(al, 1), be)
print("eliminant in ga (chart al=1):", sp.factor_list(r2))
inf = [sp.factor(c.subs(ga, 0)) for c in conics]
print("chart ga=0 (conics restricted):", inf, " common root?", sp.gcd(conics[0].subs({ga: 0, be: 1}), conics[1].subs({ga: 0, be: 1})))
# the members xD, yD: D = s - T
T = -x**2/2 - x*y/2 + x*z + sp.Rational(3, 8)*y**2 - y*z/2 + z**2
for lin in (x, y):
    h = sp.expand(lin*(s - T))
    ok = sp.Matrix([sp.Poly(h, x, y, z, s).coeff_monomial(m) for m in basis]).T
    # check membership in span of gens
    A = sp.Matrix([[sp.Poly(gg, x, y, z, s).coeff_monomial(m) for m in basis] for gg in gens])
    print("lin*D in net:", A.T.rank() == A.T.row_join(ok.T).rank())
