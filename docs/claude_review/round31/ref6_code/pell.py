"""(e) The integer analogue of the PARITY (mod 2) code condition is false, already for 2 points.
z1 = gamma^2 * pi, z0 = conj(z1), pi = 2+i (p = 5), gamma = x+iy from X^2-5y^2 = +-1, x = X-2y.
Then z1 - z0 = 2i*Im(gamma^2 pi) = 2i*(x^2+4xy-y^2) = +-2i.  The anchored parity word of c=(1) is supported
on p=5 only (exponent difference 2*ord + 1 at 5, even at all other primes), so its 'odd-support weight' is log 5,
while W = log N -> infinity.  The Z-support (radical) weight is ~W/2 >= W/4, as abc predicts.
Also a mod-3 analogue: z1 = gamma^3 pi with |Im(gamma^3 pi)| small, found by search."""
from math import log, gcd, isqrt
from sympy import factorint  # only for reporting the factorisation of N(gamma); checked by multiplication

def gmul(a, b):
    return (a[0]*b[0] - a[1]*b[1], a[0]*b[1] + a[1]*b[0])

def conj(a):
    return (a[0], -a[1])

def gauss_val(z, pi):
    # valuation of Gaussian integer z at Gaussian prime pi
    n = pi[0]**2 + pi[1]**2
    k = 0
    while True:
        q = gmul(z, conj(pi))
        if q[0] % n or q[1] % n:
            return k
        z = (q[0]//n, q[1]//n)
        k += 1

def gauss_prime_over(p):
    # find pi = a+bi with a^2+b^2 = p
    for a in range(1, isqrt(p)+1):
        b2 = p - a*a
        b = isqrt(b2)
        if b*b == b2:
            return (a, b)

pi = (2, 1)
sols = []
X, y = 2, 1  # 2^2-5 = -1
fund = (9, 4)  # 9^2-5*16 = 1
for k in range(12):
    sols.append((X, y))
    X, y = X*fund[0] + 5*y*fund[1], X*fund[1] + y*fund[0]
print("l = 2 (parity) family")
print(f"{'X':>14} {'y':>14} {'W=log N':>10} {'W/4':>8} {'odd-supp wt':>11} {'Z-supp wt':>10} {'|z1-z0|':>7} {'C=|z1-z0|/R^(1/2)':>18}")
for X, y in sols:
    for sgn in (1, -1):
        x = X - 2*y
        g = (x, sgn*y) if sgn == 1 else (x, y)
        if sgn == -1:
            continue
        g = (x, y)
        if g[0]**2 + g[1]**2 <= 1:
            continue
        z1 = gmul(gmul(g, g), pi)
        z0 = conj(z1)
        N = z1[0]**2 + z1[1]**2
        assert N == z0[0]**2 + z0[1]**2 and N % 2 == 1
        assert gcd(z1[0], z1[1]) == 1  # primitive => coprime to conj (N odd)
        diff = (z1[0]-z0[0], z1[1]-z0[1])
        dist = (diff[0]**2 + diff[1]**2) ** 0.5
        fac = factorint(N)
        odd_w, z_w = 0.0, 0.0
        for p, e in fac.items():
            assert p % 4 == 1
            P = gauss_prime_over(p)
            a1, a0 = gauss_val(z1, P), gauss_val(z0, P)
            assert a1 + gauss_val(z1, conj(P)) == e
            if (a1 - a0) % 2:
                odd_w += log(p)
            if a1 != a0:
                z_w += log(p)
        W = log(N)
        R = N ** 0.5
        print(f"{X:>14} {y:>14} {W:10.2f} {W/4:8.2f} {odd_w:11.3f} {z_w:10.2f} {dist:7.1f} {dist/R**0.5:18.3e}")

print()
print("l = 3 analogue: z1 = gamma^3 pi, need |z1 - z0| <= N^(1/4) (C = 1); odd-mod-3 support = {5}")
found = 0
best = []
for y0 in range(1, 4000):
    for x0 in range(-3*y0, 3*y0+1):
        if gcd(x0, y0) != 1 or (x0 + y0) % 2 == 0:
            continue
        g = (x0, y0)
        z1 = gmul(gmul(gmul(g, g), g), pi)
        im = z1[1]
        n = x0*x0 + y0*y0
        N = n**3 * 5
        if 4*im*im <= N ** 0.5 and n % 5:
            best.append((n, x0, y0, im))
    if len(best) >= 6:
        break
for n, x0, y0, im in best[:6]:
    N = n**3*5
    z1 = gmul(gmul(gmul((x0, y0), (x0, y0)), (x0, y0)), pi)
    assert gcd(z1[0], z1[1]) == 1
    print(f"gamma={x0}+{y0}i  N(gamma)={n}  |z1-z0|=2|Im|={2*abs(im)}  N^(1/4)={N**0.25:.1f}  W={log(N):.2f}  W/4={log(N)/4:.2f}  "
          f"mod-3 odd support weight = log 5 = {log(5):.3f}")
