"""T5 calibration: points of D approaching P along the real branch (slopes = convergents of the real root theta of
F_3(1,t)*8 = 8 + 20t - 10t^2 + 19t^3 ... in the slope t = y/x). Check: u = 0, x = log(1/X)/log c -> 1/2,
nu_+ = log m_+/log c -> 1/4, and the Gamma_3 inequality x <= 3/4 - nu_+ - (3/8) nu_- is (asymptotically) tight."""
import math
from fractions import Fraction
from math import gcd
import mpmath as mp
import sympy as sp
from usearch import f, T8

F3 = lambda a, b: 8*a**3 + 20*a**2*b - 10*a*b**2 + 19*b**3          # 8 F_3
F4p = lambda a, b: -16*a**4 + 32*a**3*b + 8*a**2*b**2 - 104*a*b**3 + 119*b**4   # 64 F_4'
mp.mp.dps = 60
t = sp.symbols('t')
roots = sp.Poly(8 + 20*t - 10*t**2 + 19*t**3, t).nroots(n=50)
theta = [r for r in roots if abs(sp.im(r)) < 1e-30][0]
theta = mp.mpf(str(sp.re(theta)))
print('real root of F_3(1,t):', theta)
# continued fraction convergents of theta: slope b'/a' = p/q
cf = []
xx = theta
for _ in range(30):
    a0 = int(mp.floor(xx)); cf.append(a0); xx = 1/(xx - a0)
pm, qm, p0, q0 = 1, 0, cf[0], 1
conv = [(p0, q0)]
for a0 in cf[1:]:
    pm, qm, p0, q0 = p0, q0, a0*p0 + pm, a0*q0 + qm
    conv.append((p0, q0))
for (pp, qq) in conv[3:16]:
    a1, b1 = qq, pp                     # direction (a',b') = (1, t) scaled
    A = a1*F3(a1, b1) * 8               # use 64F_3 * a' etc: point (F_3 a' : F_3 b' : -F_4') with common scaling
    # pi(D) on the line (x,y)=m(a',b'): z F_3 + m F_4' = 0  ->  (m : z) = (F_3 : -F_4') ; use integer forms 8F_3, 64F_4'
    m_, z_ = 8*F3(a1, b1), -F4p(a1, b1)
    a, b, c = m_*a1, m_*b1, z_
    g = gcd(gcd(a, b), c); a, b, c = a//g, b//g, c//g
    if c < 0: a, b, c = -a, -b, -c
    val = f(a, b, c); s = math.isqrt(val); assert s*s == val
    if (T8(a, b, c) - 8*s) != 0: s = -s
    u8 = 8*s - T8(a, b, c)
    m = gcd(a, b)
    mp_ = 1
    for p, e in sp.factorint(m).items():
        if p != 2 and (s - c*c) % p == 0: mp_ *= p**e
    X = max(abs(a), abs(b))/c
    lc = math.log(c)
    x = math.log(1/X)/lc; nup = math.log(mp_)/lc if mp_ > 1 else 0.0
    print('c=%.3e  8u=%d  x=%.4f  nu_+=%.4f  Gamma_3 bound 3/4-nu_+=%.4f' % (c, u8, x, nup, 0.75 - nup))
