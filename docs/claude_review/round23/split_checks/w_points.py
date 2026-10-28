"""Rational A-points on the Pell-balanced curves C_{abc} (split LM_5), beyond the integer Pell case.
Conic z^2 = m x^2 + B (m = abc(a+b+c)); rational point P0 = (a, z_a). Line through P0 with slope s = p/q (convergents of
sqrt m) meets the conic again at x(s) -> infinity. Evaluate F_0..F_4, clear denominators, make primitive, compute
H = (max/min) lcm and dist = (max-min)/min exactly; report gamma = log H / (-log dist) (-> 15/4 predicted)."""
from fractions import Fraction as Fr
from math import gcd, lcm, log, isqrt
import sys
def run(a, b, c, nconv=40):
    d = a + b + c
    sq = [a*a, b*b, c*c, d*d]
    S = sum(sq)
    e2 = sum(sq[i]*sq[j] for i in range(4) for j in range(i+1, 4))
    e3 = sum(sq[i]*sq[j]*sq[k] for i in range(4) for j in range(i+1, 4) for k in range(j+1, 4))
    m = a*b*c*d
    A = Fr(S*S, 4) + 2*m - e2; A /= 4
    B = Fr(e3 - S*m, 4)
    assert A == m
    x0 = Fr(a); z0 = -(x0**4 - Fr(S, 2)*x0**2 + m)/(2*x0)
    assert z0*z0 == A*x0*x0 + B
    # continued fraction of sqrt(m)
    r = isqrt(m); assert r*r != m
    # convergents via standard algorithm
    mm, dd, aa = 0, 1, r
    p0, p1, q0, q1 = 1, r, 0, 1
    out = []
    for it in range(nconv):
        mm = dd*aa - mm; dd = (m - mm*mm)//dd; aa = (r + mm)//dd
        p0, p1 = p1, aa*p1 + p0; q0, q1 = q1, aa*q1 + q0
        s = Fr(p1, q1)
        # second intersection: z = z0 + s(x - x0); (z0 + s(x-x0))^2 = A x^2 + B -> quadratic in x with root x0
        # coefficients: (s^2 - A) x^2 + 2 s (z0 - s x0) x + (z0 - s x0)^2 - B = 0; product of roots = const/(s^2-A)
        k = z0 - s*x0
        if s*s == A: continue
        x1 = ((k*k - B)/(s*s - A))/x0
        z1 = z0 + s*(x1 - x0)
        assert z1*z1 == A*x1*x1 + B
        X, Z = x1, z1
        F = [(X-a)*(X-b)*(X-c)*(X+d), X**4 - Fr(S, 2)*X**2 + m + 2*X*Z, (X*X-b*b)*(X*X-c*c), (X*X-a*a)*(X*X-c*c), (X*X-a*a)*(X*X-b*b)]
        if any(f == 0 for f in F): continue
        sign = 1 if F[0] > 0 else -1
        F = [sign*f for f in F]
        if any(f <= 0 for f in F): continue
        den = 1
        for f in F: den = lcm(den, f.denominator)
        u = [int(f*den) for f in F]
        g = 0
        for v in u: g = gcd(g, v)
        u = [v//g for v in u]
        L = 1
        for v in u: L = lcm(L, v)
        logH = log(L) + log(max(u)) - log(min(u))
        logdist = log(max(u) - min(u)) - log(min(u))
        if logdist < log(1e-3):
            out.append((float(X), logH, logH/(-logdist)))
    return out
for (a, b, c) in [(1, 2, 4), (1, 3, 5), (2, 3, 7), (1, 5, 6)]:
    res = run(a, b, c)
    print((a, b, c), "m =", a*b*c*(a+b+c), " last 4 (x, log H, gamma):", [("%.3g" % t[0], "%.0f" % t[1], "%.4f" % t[2]) for t in res[-4:]])
print("log(dist^{15/4} H) along convergents (bounded <=> alpha_C = 15/4 attained with bounded constant):")
for (a, b, c) in [(1, 2, 4), (1, 5, 6), (2, 3, 7)]:
    res = run(a, b, c, nconv=60)
    vals = [t[1] - 3.75*(t[1]/t[2]) for t in res]
    print((a, b, c), "n_pts", len(res), "min %.2f max %.2f" % (min(vals), max(vals)), "last 6:", ["%.2f" % v for v in vals[-6:]])
