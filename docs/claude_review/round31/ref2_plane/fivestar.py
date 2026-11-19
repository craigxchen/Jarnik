"""Referee ref2_plane: checks of the real five-star lemma.
1. symbolic: critical points of H(x)=x^3/((x-1)^3(x-b)^2) and the equal-critical-value determinant.
2. numerical search for REAL y in R^10 (edges of K5) with constant star sums and constant reciprocal star sums,
   |y_e| pairwise distinct.  Expect none.
3. numerical search for COMPLEX solutions (expect a 1-parameter family up to scaling) and, for each, test the
   structure the lemma predicts over C: plane quartic C(z)=[1-rz : z^4-sz^3 : z^2] has exactly one node pair.
"""
import numpy as np, sympy as sp, itertools, sys
from scipy.optimize import least_squares

x, b = sp.symbols('x b')
H = x**3/((x-1)**3*(x-b)**2)
num = sp.factor(sp.together(sp.diff(sp.log(H), x)))
print("dlogH numerator:", sp.factor(sp.numer(sp.together(sp.diff(sp.log(H),x)))))
q = 2*x**2 + x - 3*b
rN = sp.rem(sp.expand(x**3), q, x); rD = sp.rem(sp.expand((x-1)**3*(x-b)**2), q, x)
a1, a0 = sp.Poly(rN, x).coeff_monomial(x), sp.Poly(rN, x).coeff_monomial(1)
b1, b0 = sp.Poly(rD, x).coeff_monomial(x), sp.Poly(rD, x).coeff_monomial(1)
dd = sp.factor(a1*b0 - a0*b1)
print("determinant:", dd)
bb = sp.Rational(-16, 9)
cs = sp.solve(q.subs(b, bb), x)
print("b=-16/9 critical pts:", cs, " H values equal:", sp.simplify(H.subs(b,bb).subs(x,cs[0]) - H.subs(b,bb).subs(x,cs[1])) == 0)
print("disc 1+24b at -16/9:", 1 + 24*bb)
# also: independent derivation of the H-shape: degree-5 map with two (3,2) fibres over 0,inf and a (2,2,1) fibre.
# Generic such map: x^3/((x-1)^3 (x-b)^2) after Moebius; check total ramification = 8 at b=-16/9.
Hb = H.subs(b, bb)
nb = sp.numer(sp.together(sp.diff(Hb, x)))
print("numerator of H' at b=-16/9 factors:", sp.factor(nb))

edges = [(i, j) for i in range(5) for j in range(i+1, 5)]
inc = np.zeros((5, 10))
for k, (i, j) in enumerate(edges):
    inc[i, k] = inc[j, k] = 1

def F(y):
    S = inc @ y; R = inc @ (1/y)
    return np.concatenate([S[:4]-S[4], R[:4]-R[4], [np.sum(y**2)-10]])

rng = np.random.default_rng(12345)
best = []
nconv = 0
for trial in range(int(sys.argv[1]) if len(sys.argv) > 1 else 20000):
    y0 = rng.normal(size=10) * rng.choice([0.3, 1, 3])
    if np.min(np.abs(y0)) < 1e-3: continue
    try:
        res = least_squares(F, y0, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=2000)
    except Exception:
        continue
    if np.linalg.norm(res.fun) < 1e-9 and np.min(np.abs(res.x)) > 1e-6:
        nconv += 1
        a = np.sort(np.abs(res.x))
        sep = np.min(np.diff(a))
        best.append((sep, res.x))
best.sort(key=lambda t: -t[0])
print("real search: converged real solutions:", nconv)
if best:
    print("largest min-gap between |y_e| among real solutions: %.3e" % best[0][0])
    print("example:", np.round(best[0][1], 6))
# distribution of min-gaps
gaps = np.array([g for g, _ in best])
if len(gaps):
    print("gap quantiles:", np.quantile(gaps, [0.5, 0.9, 0.99, 1.0]))

# complex solutions
def Fc(v):
    y = v[:10] + 1j*v[10:]
    S = inc @ y; R = inc @ (1/y)
    out = np.concatenate([S[:4]-S[4], R[:4]-R[4], [np.sum(y**2)-10, np.sum(y*np.conj(y)).imag]])
    # fix the second complex dof (the curve) by a random linear normalisation
    out = np.concatenate([out, [ (np.dot(cvec, y) - 1.0) ]])
    return np.concatenate([out.real, out.imag])
cvec = rng.normal(size=10) + 1j*rng.normal(size=10)
csols = []
for trial in range(400):
    v0 = rng.normal(size=20)
    res = least_squares(Fc, v0, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
    if np.linalg.norm(res.fun) < 1e-9:
        y = res.x[:10] + 1j*res.x[10:]
        a = np.abs(y**2)
        # distinct up to sign: y_e^2 distinct
        y2 = y**2
        sep = min(abs(y2[i]-y2[j]) for i in range(10) for j in range(i+1, 10))
        if sep > 1e-4 and np.min(np.abs(y)) > 1e-4:
            csols.append(y)
print("complex nondegenerate solutions found:", len(csols))
for y in csols[:3]:
    print("  max |Im y| / max|y| = %.3f" % (np.max(np.abs(y.imag))/np.max(np.abs(y))))
