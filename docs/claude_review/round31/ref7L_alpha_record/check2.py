"""(a) exact Jacobian of seed map at fixed t; (b) resultant bound for gcd; (c) fixed T, n -> infinity; (d) linear family."""
import sympy as sp
from math import log
import sys
sys.path.insert(0, '.')
from check import family, to_expr, gdeg, gsub, ggcd_int

t = sp.symbols('t')
S = sp.symbols('s1:5')

def cayley(m, z0):
    # q = m/z0 has |q|=1; Cayley coordinate u = Im q/(1+Re q) = Im(m conj z0)/(N + Re(m conj z0))
    mr, mi = m; zr, zi = z0
    re = mr*zr + mi*zi; im = mi*zr - mr*zi
    N = zr**2 + zi**2
    return im/(N + re)

def jac(n, T0, pt):
    marks, _ = family(n, list(S))
    # marks are Polys in t with symbolic seeds -> domain issue; rebuild symbolically
    return None

def family_sym(n, seeds, tval):
    A = sp.Matrix([[tval**2+1, tval], [tval, 1]])
    An = A**n
    G = []
    for s in seeds:
        v = An*sp.Matrix([1, s]); G.append((sp.expand(v[0]), sp.expand(v[1])))
    mul = lambda a, b: (sp.expand(a[0]*b[0]-a[1]*b[1]), sp.expand(a[0]*b[1]+a[1]*b[0]))
    cj = lambda a: (a[0], -a[1])
    beta = (tval, 2)
    def prodc(idx):
        r = (1, 0)
        for l in idx: r = mul(r, cj(G[l]))
        return r
    z0 = mul(beta, prodc(range(4)))
    marks = [z0]
    for (j, k) in [(0, 1), (0, 2), (1, 2), (0, 3)]:
        marks.append(mul(mul(cj(beta), mul(G[j], G[k])), prodc([l for l in range(4) if l not in (j, k)])))
    return marks

print('(a) exact Jacobian determinant of seeds -> Cayley coords of 4 ratios, at fixed t')
for n in (1, 2):
    for T0 in (1, 3):
        marks = family_sym(n, list(S), sp.Integer(T0))
        us = [cayley(m, marks[0]) for m in marks[1:]]
        Jm = sp.Matrix([[sp.diff(u, s) for s in S] for u in us])
        pt = {S[0]: sp.Rational(1, 2), S[1]: -2, S[2]: sp.Rational(7, 3), S[3]: 4}
        d = sp.nsimplify(Jm.subs(pt).det())
        print(f'  n={n} t={T0}: det = {sp.N(d, 6)} (exact rational, nonzero: {d != 0})')

print('(b) resultant bound for the Gaussian gcd (seeds 1,2,3,5; integral)')
for n in (1, 2):
    marks, _ = family(n, [1, 2, 3, 5])
    e0 = sp.expand(to_expr(marks[0])); e1 = sp.expand(to_expr(marks[1]))
    lc0 = sp.Poly(e0, t).LC(); lc1 = sp.Poly(e1, t).LC()
    Rz = sp.expand(sp.resultant(e0, e1, t))
    NR = sp.expand(Rz*sp.conjugate(Rz))
    print(f'  n={n}: lc = {lc0},{lc1}; N(Res(z0,z12)) = {sp.factorint(NR)}')
    for T in (10, 999, 10**6+7, 10**15+37):
        vals = [(int(m[0].eval(T)), int(m[1].eval(T))) for m in marks]
        g = vals[0]
        for v in vals[1:]: g = ggcd_int(g, v)
        gn = g[0]**2+g[1]**2
        assert NR % gn == 0
    print('   gcd norms divide N(Res) at all tested T')

print('(c) fixed T, n -> infinity: distance of ratio class from e')
for T0 in (1, 2, 3):
    for n in (5, 10, 20, 40):
        marks = family_sym(n, [1, 2, 3, 5], sp.Integer(T0))
        z0 = complex(*[complex(x) for x in marks[0]]) if False else complex(float(marks[0][0]), float(marks[0][1]))
        dist = max(abs(complex(float(m[0]), float(m[1]))/z0 - 1) for m in marks[1:])
        print(f'  T={T0} n={n}: max|q-1| = {dist:.6f}')

print('(d) linear star family P_j = t + i c_j: exponent')
c = [1, 2, 3, 5]
tt = sp.symbols('tt')
for T in (10**3, 10**6, 10**12):
    Pj = [(T, cj) for cj in c]
    mul = lambda a, b: (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
    cj_ = lambda a: (a[0], -a[1])
    def prodc(idx):
        r = (1, 0)
        for l in idx: r = mul(r, cj_(Pj[l]))
        return r
    z = [prodc(range(4))] + [mul(Pj[j], prodc([l for l in range(4) if l != j])) for j in range(4)]
    g = z[0]
    for v in z[1:]: g = ggcd_int(g, v)
    R2 = z[0][0]**2+z[0][1]**2; H = R2//(g[0]**2+g[1]**2)
    diam2 = max((z[a][0]-z[b][0])**2+(z[a][1]-z[b][1])**2 for a in range(5) for b in range(a+1, 5))
    print(f'  T={T:.0e}: gamma = {log(H)/(-0.5*(log(diam2)-log(R2))):.4f}')
