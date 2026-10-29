# Killer: all points of S near P (c > max(|a|,|b|)) with max(|a|,|b|) <= A and 0 < |8u| <= Umax, via the quadratic in c.
# u = s - T(a,b,c); u(2T+u) = c*F3(a,b) + F4'(a,b). Multiply by 64: with U = 8u,
#   U(2 T8 + U) = 64 c F3 + 64 F4'  where T8 = 8T = 8c^2 + 8c L + 8Q.
import math, sys
from fractions import Fraction as Fr
A = int(sys.argv[1]) if len(sys.argv) > 1 else 400
Umax = int(sys.argv[2]) if len(sys.argv) > 2 else 24
def fz(a, b, c): return a**3*b + 2*a**2*b*c - 2*a*b**3 - 2*a*b*c**2 + 2*a*c**3 + 2*b**4 + 2*b**3*c + b**2*c**2 - b*c**3 + c**4
def T8(a, b, c): return 8*c*c + c*(8*a - 4*b) + (-4*a*a - 4*a*b + 3*b*b)
found = []
for a in range(-A, A + 1):
    for b in range(-A, A + 1):
        if a == 0 and b == 0: continue
        G3 = 8*a**3 + 20*a*a*b - 10*a*b*b + 19*b**3          # 8 F3
        G4 = -16*a**4 + 32*a**3*b + 8*a*a*b*b - 104*a*b**3 + 119*b**4   # 64 F4'
        L8 = 8*a - 4*b; Q8 = -4*a*a - 4*a*b + 3*b*b
        for U in range(-Umax, Umax + 1):
            if U == 0: continue
            # U(2 (8c^2 + L8 c + Q8) + U) = 8 c G3 + G4  (64 u(2T+u) = 64 c F3 + 64 F4')
            # 16 U c^2 + (2 U L8 - 8 G3) c + (2 U Q8 + U^2 - G4) = 0
            qa = 16*U; qb = 2*U*L8 - 8*G3; qc = 2*U*Q8 + U*U - G4
            D = qb*qb - 4*qa*qc
            if D < 0: continue
            r = math.isqrt(D)
            if r*r != D: continue
            for num in (-qb + r, -qb - r):
                if num % (2*qa): continue
                c = num // (2*qa)
                if c <= max(abs(a), abs(b)): continue
                if math.gcd(math.gcd(a, b), c) != 1: continue
                s8 = T8(a, b, c) + U
                if s8 % 8: continue
                s = s8 // 8
                assert s*s == fz(a, b, c)
                X = max(abs(a), abs(b))/c
                found.append((a, b, c, s, U, X**1.5*c, -math.log(X)/math.log(c)))
for row in sorted(set(found)):
    print(row)
print("total", len(set(found)))
