"""Evidence: rational points of S near P with small u = s - T(a,b,c) (near-tight Liouville).
For each (a,b) with |a|,|b| <= Amax and U = 8u in [-Umax, Umax] minus 0, solve 64 f(a,b,c) = (8T(a,b,c) + U)^2,
a quadratic in c (exactly), keep integer roots c > max(|a|,|b|), primitive (a,b,c), s = (8T+U)/8 integer.
Report rho = X^{3/2} H (Liouville ratio, >= const by Liouville), x = log(1/X)/log c, m = gcd(a,b), m_+, m_-.
Also D-points (u = 0) from slopes near the real root of F3, to check x -> 1/2, nu_+ -> 1/4 (section 7, T5)."""
import math, sys
from math import gcd, isqrt
import sympy as sp

def f(a, b, c):
    x, y, z = a, b, c
    return (x**3*y + 2*x**2*y*z - 2*x*y**3 - 2*x*y*z**2 + 2*x*z**3 + 2*y**4 + 2*y**3*z + y**2*z**2 - y*z**3 + z**4)

def T8(a, b, c):
    x, y, z = a, b, c
    return 8*z*z + 8*x*z - 4*y*z - 4*x*x - 4*x*y + 3*y*y

def quad_coeffs(a, b, U):
    vals = [64*f(a, b, c) - (T8(a, b, c) + U)**2 for c in (0, 1, 2, 3)]
    # degree <= 2 in c ? check third difference
    d1 = [vals[i+1] - vals[i] for i in range(3)]
    d2 = [d1[i+1] - d1[i] for i in range(2)]
    assert d2[0] == d2[1], (a, b, U)
    A2 = d2[0] // 2
    A1 = d1[0] - A2
    A0 = vals[0]
    return A2, A1, A0

def int_roots(A2, A1, A0):
    if A2 == 0:
        if A1 != 0 and A0 % A1 == 0:
            return [-A0 // A1]
        return []
    disc = A1*A1 - 4*A2*A0
    if disc < 0:
        return []
    r = isqrt(disc)
    if r*r != disc:
        return []
    out = []
    for sgn in (1, -1):
        num = -A1 + sgn*r
        if num % (2*A2) == 0:
            out.append(num // (2*A2))
    return out

def mpm(a, b, c, s):
    m = gcd(a, b)
    mp = mm = 1
    for p, e in sp.factorint(m).items():
        if p == 2:
            continue
        if (s - c*c) % p == 0:
            mp *= p**e
        elif (s + c*c) % p == 0:
            mm *= p**e
    return m, mp, mm

if __name__ == '__main__':
    Amax = int(sys.argv[1]); Umax = int(sys.argv[2])
    found = []
    for a in range(-Amax, Amax + 1):
        for b in range(-Amax, Amax + 1):
            if a == 0 and b == 0:
                continue
            for U in range(-Umax, Umax + 1):
                if U == 0:
                    continue
                A2, A1, A0 = quad_coeffs(a, b, U)
                for c in int_roots(A2, A1, A0):
                    if c <= max(abs(a), abs(b)) or gcd(gcd(a, b), c) != 1:
                        continue
                    num = T8(a, b, c) + U
                    if num % 8:
                        continue
                    s = num // 8
                    assert s*s == f(a, b, c)
                    X = max(abs(a), abs(b)) / c
                    found.append((X**1.5*c, math.log(1/X)/math.log(c), a, b, c, s, U))
    found.sort()
    print('points with 0 < |8u| <= %d, |a|,|b| <= %d, c > max(|a|,|b|): %d' % (Umax, Amax, len(found)))
    for rho, xx, a, b, c, s, U in found[:25]:
        m, mp, mm = mpm(a, b, c, s)
        print('rho=X^1.5*H=%.4f  x=%.4f  (a,b,c)=(%d,%d,%d)  8u=%d  m=%d m+=%d m-=%d' % (rho, xx, a, b, c, U, m, mp, mm))
