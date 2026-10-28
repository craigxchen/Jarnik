# Judge J2 (independent code): Theorem 2 of split.md at NEW parameter points, and the m-czech divisors.
import sympy as sp
from math import isqrt, gcd
from fractions import Fraction as Fr
from itertools import combinations
x, z = sp.symbols('x z')
def family(a,b,c):
    a,b,c = map(sp.Rational,(a,b,c)); d = a+b+c
    sq = [a**2,b**2,c**2,d**2]
    S = sum(sq); e3 = sum(sp.prod(t) for t in combinations(sq,3))
    m = a*b*c*d; B = (e3 - S*m)/4
    F = [(x-a)*(x-b)*(x-c)*(x+d), x**4-(S/2)*x**2+m+2*x*z,
         (x**2-b**2)*(x**2-c**2), (x**2-a**2)*(x**2-c**2), (x**2-a**2)*(x**2-b**2)]
    return a,b,c,d,m,B,F
def check(a,b,c):
    a,b,c,d,m,B,F = family(a,b,c)
    n = sp.prod([x**2-r**2 for r in (a,b,c,d)])
    norm = sp.expand((x**4-(sp.Rational(1,2)*sum([a**2,b**2,c**2,d**2]))*x**2+m)**2 - 4*x**2*(m*x**2+B) - n)
    places = []
    for x0 in (a,-a,b,-b,c,-c,d,-d):
        z2 = m*x0**2+B; z0 = sp.sqrt(z2)
        assert z0.is_rational and z0 != 0, (x0, z2)
        for s in (1,-1):
            v = []
            for j,f in enumerate(F):
                if j == 1:
                    val = f.subs({x:x0, z:s*z0}); v.append(0 if val != 0 else 1)  # simple zero (n has simple roots)
                else:
                    v.append(sp.Poly(f,x).as_expr().subs(x,x0) == 0 and 1 or 0)
            places.append(((x0,s), tuple(v)))
    cuts = []; const = 0
    for p,v in places:
        A = frozenset(i for i in range(5) if v[i])
        if len(A) in (0,5): const += 1; continue
        cuts.append(min(A, frozenset(range(5))-A, key=lambda s: (len(s), sorted(s))))
    BC = sum(max(v)-min(v) for p,v in places)
    print((a,b,c), "m=%s B=%s" % (m,B), "norm identity:", norm == 0, " const places:", const,
          " distinct cuts:", len(set(cuts)), "of", len(cuts), " B.C =", BC,
          " m square?", sp.sqrt(m).is_rational)
for t in [(1,2,4),(2,5,9),(1,4,6),(3,4,10),(2,7,8)]:
    check(*t)
# m-czech: z^2 = 56 x^2 + 65 ; F_1 = x^4 - 35 x^2 + 56 + 2 x z ; window (sqrt n, sqrt n + C n^{1/4})
def icbrt4(N):  # floor(N^(1/4))
    r = isqrt(isqrt(N))
    while (r+1)**4 <= N: r += 1
    while r**4 > N: r -= 1
    return r
X, Z = 1, 11
print("m-czech Pell points (exact): x, #divisors in window, max offset / n^(1/4) (approx), offsets/x^2")
for k in range(6):
    n = (X**2-1)*(X**2-4)*(X**2-16)*(X**2-49)
    Fs = [(X-1)*(X-2)*(X-4)*(X+7), X**4-35*X**2+56+2*X*Z, (X**2-4)*(X**2-16), (X**2-1)*(X**2-16), (X**2-1)*(X**2-4)]
    if X > 7:
        r = isqrt(n); assert r*r < n
        ok = all(n % f == 0 and f > r for f in Fs)
        lo = min(Fs); g = 0
        for f in Fs: g = gcd(g, f)
        C = Fr(max(Fs) - isqrt(n), icbrt4(n))
        print(X, ok, float(C), [float(Fr(f - Fs[0], X*X)) for f in Fs], "gcd", g)
    X, Z = 2*Z+15*X, 15*Z+112*X
