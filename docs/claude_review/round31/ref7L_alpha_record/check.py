"""Independent exact check of the polynomial Pell five-mark family and its integer specialisations.
Gaussian polynomials are pairs (re, im) of sympy Polys over QQ in t (t real)."""
import sys
from fractions import Fraction as F
from math import log
import sympy as sp

t = sp.symbols('t')
P = lambda e: sp.Poly(e, t, domain='QQ')

def gmul(a, b):
    return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def gconj(a):
    return (a[0], -a[1])
def gsub(a, b):
    return (a[0]-b[0], a[1]-b[1])
def gdeg(a):
    return max(a[0].degree() if not a[0].is_zero else -1, a[1].degree() if not a[1].is_zero else -1)

def family(n, seeds):
    A = sp.Matrix([[t**2+1, t], [t, 1]])
    An = A**n
    G = []
    for s in seeds:
        v = An*sp.Matrix([1, s])
        G.append((P(sp.expand(v[0])), P(sp.expand(v[1]))))
    beta = (P(t), P(2))
    one = (P(1), P(0))
    def prodc(idx):
        r = one
        for l in idx:
            r = gmul(r, gconj(G[l]))
        return r
    z0 = gmul(beta, prodc(range(4)))
    marks = [z0]
    for (j, k) in [(0, 1), (0, 2), (1, 2), (0, 3)]:
        r = gmul(gmul(gconj(beta), gmul(G[j], G[k])), prodc([l for l in range(4) if l not in (j, k)]))
        marks.append(r)
    return marks, G

def to_expr(a):
    return a[0].as_expr() + sp.I*a[1].as_expr()

def ggcd_int(a, b):
    # Gaussian integer gcd, a,b as python complex-int tuples
    def mod(x, y):
        # x - q y with q nearest
        nr = x[0]*y[0]+x[1]*y[1]; ni = x[1]*y[0]-x[0]*y[1]; d = y[0]**2+y[1]**2
        q = ((2*nr+d)//(2*d), (2*ni+d)//(2*d))
        return (x[0]-(q[0]*y[0]-q[1]*y[1]), x[1]-(q[0]*y[1]+q[1]*y[0]))
    while b != (0, 0):
        a, b = b, mod(a, b)
    return a

def run(n, seeds, Ts):
    marks, G = family(n, seeds)
    D = gdeg(marks[0])
    norms = [(m[0]**2 + m[1]**2) for m in marks]
    assert all(N == norms[0] for N in norms), 'norms differ'
    degs = [gdeg(m) for m in marks]
    pd = [gdeg(gsub(marks[a], marks[b])) for a in range(5) for b in range(a+1, 5)]
    # gcd over Q(i)[t]
    ex = [sp.expand(to_expr(m)) for m in marks]
    g = ex[0]
    for e in ex[1:]:
        g = sp.gcd(g, e, extension=sp.I)
    print(f'n={n} seeds={seeds} D={D} degs={set(degs)} pairdegs={sorted(set(pd))} gcd={g} 2D/h={F(2*D, D-max(pd))}')
    assert D == 8*n+1 and set(pd) == {4*n+1} and sp.degree(g, t) == 0
    # Bezout constant: s0*z0 + s1*z1 = 1 over Q(i)[t] (if z0,z1 coprime)
    s_, s1_, h_ = sp.gcdex(ex[0], ex[1], t, extension=sp.I) if False else (None, None, None)
    # denominators of seeds
    den = 1
    for s in seeds:
        den = sp.ilcm(den, sp.fraction(sp.Rational(s))[1])
    for T in Ts:
        vals = []
        for m in marks:
            re = m[0].eval(T)*den**4; im = m[1].eval(T)*den**4
            assert re.q == 1 and im.q == 1
            vals.append((int(re), int(im)))
        g0 = vals[0]
        for v in vals[1:]:
            g0 = ggcd_int(g0, v)
        gn = g0[0]**2+g0[1]**2
        R2 = vals[0][0]**2+vals[0][1]**2
        assert all(v[0]**2+v[1]**2 == R2 for v in vals)
        assert R2 % gn == 0
        H = R2//gn  # primitive norm R''^2
        diam2 = max((vals[a][0]-vals[b][0])**2+(vals[a][1]-vals[b][1])**2 for a in range(5) for b in range(a+1, 5))
        mind2 = min((vals[a][0]-vals[b][0])**2+(vals[a][1]-vals[b][1])**2 for a in range(5) for b in range(a+1, 5))
        assert mind2 > 0
        # Delta = diam/R ; exponent gamma = log H / log(1/Delta)
        lD = 0.5*(log(diam2)-log(R2))
        gam = log(H)/(-lD)
        arcexp = 0.5*log(diam2)/(0.5*log(H))  # diameter as power of R''
        print(f'   T={T:.0e} N(gcd)={gn} log10 H={log(H)/log(10):.1f} gamma={gam:.4f} diam=R\'\'^{arcexp:.4f}')
    return marks

if __name__ == "__main__" and False:
    for n in (1, 2, 3):
        for seeds in ([sp.Rational(1, 2), -2, sp.Rational(7, 3), 4], [1, 2, 3, 5], [sp.Rational(-3, 5), sp.Rational(11, 7), 6, sp.Rational(-1, 9)]):
            run(n, seeds, [10, 10**3, 10**6, 10**12, 10**30])
