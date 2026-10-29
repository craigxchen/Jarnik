"""Fibre sections on E_l: F in L(k S_inf) = span{m^i c^{k-i} (i<=k), s m^j c^{k-2-j} (j<=k-2)} with
ord_P F >= a, ord_P' F >= b, a+b = 2k-1 (a >= k+1). In t = m/c: F/c^k = sum alpha_i t^i + w(t) sum beta_j t^j, and
at P' the branch is -w(t). The coefficients (alpha, beta) are maximal minors; check
  deg_{(a',b')} beta_j = -mu + j,  mu = (2k^2 + 1 - a^2 - b^2)/2,
and that the minor vector is not identically zero. Also checks Gamma_k classes are (-1)-sections on Bl_{P,P'}S."""
import sympy as sp
from sympy import symbols, expand, Poly, Matrix, sqrt, series
from geom import default_f
x, y, z, t, A, B = symbols('x y z t A B')
f = default_f()

def wcoeffs(N):
    g = expand(f.subs({x: t*A, y: t*B, z: 1}))
    ser = expand(series(sp.sqrt(g), t, 0, N + 1).removeO())
    return [expand(ser.coeff(t, n)) for n in range(N + 1)]

W = wcoeffs(16)
for n in range(1, 9):
    P_ = Poly(W[n], A, B)
    assert all(sum(mon) == n for mon in P_.monoms()), n
print('w_n homogeneous of degree n in (a\',b\') for n<=8: ok')

def fibre_section(k, a, b):
    # unknowns alpha_0..alpha_k, beta_0..beta_{k-2}; Taylor coefficient of t^i at P: alpha_i + sum_j beta_j w_{i-j}
    # at P': alpha_i - sum_j beta_j w_{i-j}
    nal, nbe = k + 1, k - 1
    rows = []
    for i in range(a):
        rows.append([1 if (ii == i) else 0 for ii in range(nal)] + [W[i - j] if i - j >= 0 else 0 for j in range(nbe)])
    for i in range(b):
        rows.append([1 if (ii == i) else 0 for ii in range(nal)] + [-W[i - j] if i - j >= 0 else 0 for j in range(nbe)])
    M = Matrix(rows)
    ns = M.nullspace()
    assert len(ns) == 1, (k, a, b, len(ns))
    v = ns[0]
    # clear denominators -> polynomial vector in A,B
    den = sp.lcm([sp.fraction(sp.together(e))[1] for e in v])
    v = [sp.factor(sp.together(e*den)) for e in v]
    return v

for (k, a, b) in [(3, 5, 0), (3, 4, 1), (4, 7, 0), (4, 6, 1), (5, 9, 0)]:
    mu = sp.Rational(2*k*k + 1 - a*a - b*b, 2)
    v = fibre_section(k, a, b)
    # remove common polynomial factor
    g = v[0]
    for e in v[1:]:
        g = sp.gcd(g, e)
    v = [sp.cancel(e/g) for e in v]
    degs_beta = [Poly(sp.numer(sp.together(v[k + 1 + j])), A, B).total_degree() if v[k + 1 + j] != 0 else None
                 for j in range(k - 1)]
    degs_alpha = [Poly(sp.numer(sp.together(v[i])), A, B).total_degree() if v[i] != 0 else None for i in range(k + 1)]
    print('k=%d a=%d b=%d: mu=%s, predicted deg beta_j = %s; found beta degs %s, alpha degs %s, common factor %s'
          % (k, a, b, mu, [-mu + j for j in range(k - 1)], degs_beta, degs_alpha, sp.factor(g)))

# Gamma_k classes on Bl_{P,P'} S: (d,a,b) = (k^2-k, k^2-1, k^2-2k); L^2=2, E^2=-1, K = -L+E_P+E_P'
kk = symbols('k')
d_, a_, b_ = kk**2 - kk, kk**2 - 1, kk**2 - 2*kk
print('Gamma_k^2 =', sp.expand(2*d_**2 - a_**2 - b_**2), '; K.Gamma_k =', sp.expand(-2*d_ + a_ + b_),
      '; F.Gamma_k =', sp.expand(2*d_ - a_ - b_))
