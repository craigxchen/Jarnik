"""Independent referee check of split.md Theorem 2 (i)-(iii),(v).
Symbolic in a,b,c; exact rank of the Jacobian at rational points via a rational parametrisation of the conic."""
import sympy as sp, itertools, random
a, b, c, x, z, t = sp.symbols('a b c x z t')
d = a + b + c
sq = [a**2, b**2, c**2, d**2]
S = sum(sq)
e2 = sum(p*q for p, q in itertools.combinations(sq, 2))
e3 = sum(p*q*r for p, q, r in itertools.combinations(sq, 3))
m = a*b*c*d
B = sp.expand((e3 - S*m)/4)
conic = m*x**2 + B  # z^2 = conic
F0 = (x-a)*(x-b)*(x-c)*(x+d)
F1 = x**4 - S/2*x**2 + m + 2*x*z
F1b = x**4 - S/2*x**2 + m - 2*x*z
F2 = (x**2-b**2)*(x**2-c**2); F3 = (x**2-a**2)*(x**2-c**2); F4 = (x**2-a**2)*(x**2-b**2)
n = sp.prod([x**2 - r for r in sq])
# (i) norm identity
N1 = sp.expand(sp.expand(F1*F1b).subs(z**2, conic))
print("(i) F1*F1bar - n == 0 on conic:", sp.simplify(N1 - sp.expand(n)) == 0)
print("    x^2-coefficient of conic minus m:", sp.factor(sp.expand((S**2/4 + 2*m - e2)/4 - m)))
for r in [a, b, c, d]:
    val = sp.factor(sp.expand(m*r**2 + B))
    print("    m r^2 + B at r =", r, ":", val)
# (iii) differences F_j - F0: coefficients in x
for j, Fj in enumerate([F1, F2, F3, F4], 1):
    P = sp.Poly(sp.expand(Fj - F0), x)
    print(f"(iii) F{j}-F0 coeffs (x^4,x^3,x^2):", [sp.factor(P.coeff_monomial(x**k)) for k in (4, 3, 2)])
# (ii) valuation vectors at the 16 places, exactly at several parameter points
def profile(av, bv, cv):
    sub = {a: av, b: bv, c: cv}
    Fs = [sp.expand(F.subs(sub)) for F in (F0, F1, F2, F3, F4)]
    con = sp.expand(conic.subs(sub)); dv = av + bv + cv
    places = []
    for r in [av, -av, bv, -bv, cv, -cv, dv, -dv]:
        zz2 = con.subs(x, r); zr = sp.sqrt(zz2)
        assert zr.is_rational and zr != 0, (r, zz2)
        for s in (zr, -zr):
            # order of vanishing in local parameter (x - r) since z != 0 there: z is a unit-analytic fn of x
            vec = []
            for F in Fs:
                # substitute z = s*sqrt(conic/zz2) series in (x-r)
                h = sp.symbols('h')
                zser = sp.series(s*sp.sqrt(con.subs(x, r+h)/zz2), h, 0, 6).removeO()
                ser = sp.expand(F.subs(z, zser).subs(x, r+h))
                k = 0
                while sp.simplify(ser.coeff(h, k)) == 0 and k < 6: k += 1
                vec.append(k)
            places.append(((r, s), tuple(vec)))
    return places
cuts_all = set()
for A in itertools.chain.from_iterable(itertools.combinations(range(5), k) for k in (1, 2)):
    cuts_all.add(frozenset(A))
def cut_of(vec):
    mn = min(vec); v = [w - mn for w in vec]
    assert max(v) <= 1, vec
    A = frozenset(i for i in range(5) if v[i] == 1)
    if not A: return None
    Ac = frozenset(range(5)) - A
    return A if len(A) <= 2 else Ac
for (av, bv, cv) in [(1, 2, 4), (1, 3, 5), (2, 3, 7), (1, 5, 6), (sp.Rational(1), sp.Rational(3, 2), sp.Rational(5, 3)), (3, 5, 11)]:
    pl = profile(sp.Rational(av), sp.Rational(bv), sp.Rational(cv))
    cs = [cut_of(v) for _, v in pl]
    const = sum(1 for q in cs if q is None)
    hit = [q for q in cs if q is not None]
    ok = (const == 1 and len(hit) == 15 and set(hit) == cuts_all)
    print(f"(ii) (a,b,c)=({av},{bv},{cv}): constant places={const}, cuts hit={len(set(hit))}/15, each once={ok}")
# (v) exact Jacobian rank at rational points: parametrise the conic through (a, z_a)
def jac_rank(av, bv, cv, tv):
    sub = {a: av, b: bv, c: cv}
    con = sp.expand(conic.subs(sub))
    za = sp.sqrt(con.subs(x, av))
    # line through (av, za) with slope t: z = za + t (x - av); second intersection
    X = sp.symbols('X')
    eq = sp.expand((za + t*(X - av))**2 - con.subs(x, X))
    q, r_ = sp.div(sp.Poly(eq, X), sp.Poly(X - av, X))
    assert r_.is_zero
    xs = sp.solve(q.as_expr(), X)[0]
    return xs
# Use generic symbolic parametrisation in (a,b,c,t): x(t), z(t) rational functions
X = sp.symbols('X')
za_sym = sp.sqrt(sp.factor(conic.subs(x, a)))
print("z_a^2 factor:", sp.factor(conic.subs(x, a)))
