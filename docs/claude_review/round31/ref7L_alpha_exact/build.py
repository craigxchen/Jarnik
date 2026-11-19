"""Exact check of the five-mark Pell family (own code).
A(t)=[[t^2+1,t],[t,1]], beta=t+2i, seeds s_j=p_j/c_j, integer vectors v_j=(c_j,p_j).
G_j = X_j + i Y_j, (X_j,Y_j)^T = A^n v_j  (in Z[t]).
marks: z0 = beta*prod conj(G_l); z_jk = conj(beta) G_j G_k prod_{l!=j,k} conj(G_l), jk in 12,13,23,14.
conj = coefficientwise complex conjugation (t real).
"""
import sys, math, json
from fractions import Fraction
from sympy import symbols, I, Poly, QQ_I, ZZ_I, Matrix, expand, conjugate, ilcm, factorint, Rational

t = symbols('t', real=True)

def Gpolys(n, seeds):
    A = Matrix([[t**2 + 1, t], [t, 1]])
    An = A**n
    out = []
    for s in seeds:
        s = Fraction(s)
        v = Matrix([s.denominator, s.numerator])
        w = (An * v).applyfunc(expand)
        out.append(expand(w[0] + I * w[1]))
    return out

def cj(e):
    return expand(conjugate(e))

EDGES = [(0, 1), (0, 2), (1, 2), (0, 3)]
NAMES = ['z0', 'z12', 'z13', 'z23', 'z14']

def marks(n, seeds):
    G = Gpolys(n, seeds)
    Gc = [cj(g) for g in G]
    beta = t + 2 * I
    z = []
    e = beta
    for l in range(4):
        e *= Gc[l]
    z.append(expand(e))
    for (j, k) in EDGES:
        e = (t - 2 * I) * G[j] * G[k]
        for l in range(4):
            if l not in (j, k):
                e *= Gc[l]
        z.append(expand(e))
    return z

def P(e):
    return Poly(e, t, domain=QQ_I)

def coeffs_int(e):
    """list of (re,im) python ints, lowest degree first; requires Z[i] coefficients"""
    p = Poly(e, t)
    d = p.degree()
    out = []
    for k in range(d + 1):
        c = p.coeff_monomial(t**k)
        re, im = c.as_real_imag()
        assert re.is_Integer and im.is_Integer
        out.append((int(re), int(im)))
    return out

def poly_checks(n, seeds):
    z = marks(n, seeds)
    Pz = [P(e) for e in z]
    D = Pz[0].degree()
    degs = [p.degree() for p in Pz]
    norms = [expand(e * cj(e)) for e in z]
    one_norm = all(expand(norms[i] - norms[0]) == 0 for i in range(5))
    g = Pz[0]
    for p in Pz[1:]:
        g = g.gcd(p)
    diff = {}
    for a in range(5):
        for b in range(a + 1, 5):
            diff[NAMES[a] + '-' + NAMES[b]] = (Pz[a] - Pz[b]).degree()
    # pairwise gcd degrees too
    pg = {}
    for a in range(5):
        for b in range(a + 1, 5):
            pg[NAMES[a] + ',' + NAMES[b]] = Pz[a].gcd(Pz[b]).degree()
    # Bezout: sum a_i z_i = 1 over Q(i)[t]
    h = Pz[0]
    co = [Poly(1, t, domain=QQ_I)] + [Poly(0, t, domain=QQ_I)] * 4
    for i in range(1, 5):
        s_, t_, h2 = h.gcdex(Pz[i])
        co = [s_ * c for c in co]
        co[i] = co[i] + t_
        h = h2
    assert h.degree() == 0
    hc = h.LC()
    co = [c * Poly(1 / hc, t, domain=QQ_I) for c in co]
    tot = Poly(0, t, domain=QQ_I)
    for c, p in zip(co, Pz):
        tot = tot + c * p
    assert tot == Poly(1, t, domain=QQ_I)
    # common denominator: lcm of denominators of real and imaginary parts of all coefficients
    dens = []
    for c in co:
        e = c.as_expr()
        for k in Poly(e, t).all_coeffs():
            re, im = k.as_real_imag()
            dens += [Rational(re).q, Rational(im).q]
    N = int(ilcm(*dens))
    return dict(D=D, degs=degs, one_norm=one_norm, gcd_deg=g.degree(), diff_deg=diff, pair_gcd_deg=pg, bezout_N=N), z

if __name__ == '__main__':
    n = int(sys.argv[1]); seeds = [Fraction(x) for x in sys.argv[2].split(',')]
    res, z = poly_checks(n, seeds)
    res['bezout_N_factor'] = {str(k): v for k, v in factorint(res['bezout_N']).items()}
    res['coeffs'] = [coeffs_int(e) for e in z]
    res['n'] = n; res['seeds'] = [str(s) for s in seeds]
    out = sys.argv[3]
    json.dump(res, open(out, 'w'))
    r2 = dict(res); del r2['coeffs']
    print(json.dumps(r2, indent=1))
