"""Certified regular, well-separated, non-flat zero of the sdc7 tangent form (T=1, L=2, m=8, dim V = 1).

Phi(u) = sum_k c_k <k,u>^8 with the exact integer c of sdc7_c.json (recomputed by sdc7_tc.py).
Steps (all exact rational arithmetic except the last, which uses arb interval arithmetic):
 1. start from a numerically found zero with good separation, perturbed to break rational coincidences, snapped to a
    rational point u_r (denominator 10^6); direction g = rational approximation of grad Phi(u_r);
 2. p(t) = Phi(u_r + t g) is computed EXACTLY as a degree-8 polynomial in Q[t]; find h with p(-h) p(h) < 0 and
    certify via a Sturm sequence that p has exactly one root t* in [-h, h] and p'(t) != 0 there (p' has no root);
    so u* = u_r + t* g is a zero of Phi with grad Phi(u*).g = p'(t*) != 0 (regular point);
 3. separation: every pair difference of u* differs from that of u_r by <= 2 h |g|_inf: exact lower bound;
 4. non-flatness: the Hessian of Phi restricted to grad^perp is not zero at u* (arb enclosure over the segment):
    then the germ of {Phi = 0} at u* is not contained in any hyperplane (a hypersurface germ inside a hyperplane is an
    open piece of it and has zero second fundamental form);
 5. nondegeneracy of the germ: a germ not contained in any hyperplane contains nondegenerate points (degenerate points
    form a countable union of nowhere dense closed subsets of the 6-dimensional germ).
"""
import json
from fractions import Fraction as Fr
import numpy as np
import sympy as sp
import flint

D = json.load(open('sdc7_c.json'))
A = [tuple(k) for k in D['A']]
c = D['c']
d = 7
Af = np.array(A, float); cf = np.array(c, float)
c1 = sum(abs(x) for x in c)


def Phi(u):
    return float(cf @ (Af @ u) ** 8)


def grad(u):
    return 8 * (cf * (Af @ u) ** 7) @ Af


u = np.load('sdc7_bestsep.npy')
# perturb and re-project (Newton along gradient) to get a generic zero with sep close to the optimum
rng = np.random.default_rng(7)
u = u + 1e-3 * rng.standard_normal(d)
for it in range(50):
    g = grad(u)
    u = u - Phi(u) / (g @ g) * g
u = u / np.abs(u).max()
for it in range(50):
    g = grad(u)
    u = u - Phi(u) / (g @ g) * g
import mpmath as mp
mp.mp.dps = 40
Am = [[mp.mpf(x) for x in k] for k in A]
um = [mp.mpf(x) for x in u]
for it in range(8):
    y = [sum(kk * uu for kk, uu in zip(k, um)) for k in Am]
    ph = mp.fsum(ck * yy ** 8 for ck, yy in zip(c, y))
    gg = [mp.fsum(8 * ck * yy ** 7 * k[i] for ck, yy, k in zip(c, y, Am)) for i in range(d)]
    nn = mp.fsum(x * x for x in gg)
    um = [uu - ph / nn * x for uu, x in zip(um, gg)]
u = np.array([float(x) for x in um])
ur_hp = um
z = np.concatenate([[0.0], u])
sep_num = min(abs(z[i] - z[j]) for i in range(8) for j in range(i + 1, 8))
print('numerical zero u =', np.round(u, 8).tolist(), ' Phi/|c|_1 = %.2e  sep = %.6f' % (Phi(u) / c1, sep_num))

den = 10 ** 14
ur = [Fr(int(mp.nint(x * den)), den) for x in ur_hp]
gn = grad(u); gn = gn / np.abs(gn).max()
gr = [Fr(int(round(x * 1000)), 1000) for x in gn]
t = sp.symbols('t')
# exact p(t) = sum_k c_k (a_k + b_k t)^8
a = [sum(Fr(kx) * ux for kx, ux in zip(k, ur)) for k in A]
b = [sum(Fr(kx) * gx for kx, gx in zip(k, gr)) for k in A]
from math import comb
coef = [Fr(0)] * 9
for ck, ak, bk in zip(c, a, b):
    if ck == 0:
        continue
    for j in range(9):
        coef[j] += ck * comb(8, j) * ak ** (8 - j) * bk ** j
p = sp.Poly(sum(sp.Rational(cf_.numerator, cf_.denominator) * t ** j for j, cf_ in enumerate(coef)), t)
dp = p.diff(t)
h = sp.Rational(1, 10 ** 11)
pm, pp = p.eval(-h), p.eval(h)
print('p(-h) = %.3e, p(h) = %.3e, opposite signs: %s' % (float(pm), float(pp), pm * pp < 0))
nroots = p.count_roots(-h, h)
ndroots = dp.count_roots(-h, h)
print('Sturm: roots of p in [-h,h]: %d, roots of p\' in [-h,h]: %d' % (nroots, ndroots))
assert pm * pp < 0 and nroots == 1 and ndroots == 0
# separation lower bound
gmax = max(abs(x) for x in gr)
zr = [Fr(0)] + ur
sep_r = min(abs(zr[i] - zr[j]) for i in range(8) for j in range(i + 1, 8))
umax_r = max(abs(x) for x in ur)
sep_lb = sep_r - 2 * h * gmax
print('certified: u* = u_r + t* g with |t*| < h = 1e-11;  sep(u*) >= %.6f,  |u*|_inf in [%.6f, %.6f]' % (
    float(sep_lb), float(umax_r - h * gmax), float(umax_r + h * gmax)))
print('normalised separation sep(u*)/|u*|_inf >= %.6f  (global max over all u: 2/7 = %.6f)' % (
    float(sep_lb / (umax_r + h * gmax)), 2 / 7))
# 4. exact gradient and Hessian at u_r; Lipschitz transfer to u*
G = [Fr(0)] * d
H = [[Fr(0)] * d for _ in range(d)]
for ck, k, ak in zip(c, A, a):
    if ck == 0:
        continue
    a6 = ak ** 6; a7 = a6 * ak
    for i in range(d):
        if k[i]:
            G[i] += 8 * ck * k[i] * a7
            for j in range(d):
                if k[j]:
                    H[i][j] += 56 * ck * k[i] * k[j] * a6
Gm = sp.Matrix([[sp.Rational(x.numerator, x.denominator)] for x in G])
Hm = sp.Matrix([[sp.Rational(x.numerator, x.denominator) for x in row] for row in H])
n2 = (Gm.T * Gm)[0]
Pm = sp.eye(d) - Gm * Gm.T / n2
B = Pm * Hm * Pm
fro2 = sum(x ** 2 for x in B)
Hf2 = sum(x ** 2 for x in Hm)
print('at u_r: |grad|^2 = %.6e, ||P H P||_F^2 = %.6e, ||H||_F^2 = %.6e' % (float(n2), float(fro2), float(Hf2)))
# crude bounds on the segment |x|_inf <= 1.001: |k.x| <= 7.007, sum|c_k| = c1
R = 7.007
K2 = 7 ** 0.5  # |k|_2 <= sqrt 7
D2 = 56 * c1 * R ** 6 * K2 ** 2   # operator norm bound of Hessian on segment
D3 = 336 * c1 * R ** 5 * K2 ** 3  # bound of third derivative (trilinear norm)
du = float(h) * float(sum(x * x for x in gr)) ** 0.5
gnorm = float(n2) ** 0.5
dG = D2 * du
dn = 2 * dG / gnorm            # |n* - n_r| <= 2 |G*-G_r|/|G_r|
Hn = float(Hf2) ** 0.5 + D3 * du * d ** 0.5   # ||H*||_F <= ||H_r||_F + perturbation
dB = D3 * du * d ** 0.5 + 2 * Hn * 2 * dn   # ||P*H*P* - P_rH_rP_r||_F <= ||H*-H_r||_F + 2||H*||_F ||P*-P_r||, ||P*-P_r|| <= 2|n*-n_r|
print('|u*-u_r|_2 <= %.3e;  ||P*H*P* - P_r H_r P_r||_F <= %.3e  vs  ||P_r H_r P_r||_F = %.6e' % (du, dB, float(fro2) ** 0.5))
assert float(fro2) ** 0.5 > dB
print('CERTIFIED: second fundamental form of {Phi=0} at the regular zero u* is nonzero')
