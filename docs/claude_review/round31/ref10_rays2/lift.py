# Referee check of the transfer of X.1 (n = 2) used for (X.3) on a 5 mod 8 parent.
# For r = pi*gamma/conj(gamma) = alpha/b (lowest terms), alpha^2 = x + i y:
#  (a) x-y, x+y are odd, nonzero, prime to every parent prime and to J = p b^2;
#  (b) their prime-power factors are <= sqrt2 * p * Norm(gamma)  (cutoff (*) at n=2);
#  (c) the two forced sine lower bounds at units i, -i square-sum to b^4 (> 1 iff b > 1).
from sympy import primerange, factorint
from math import gcd, isqrt
import random
random.seed(7)

def gfactor(p):
    for x in range(1, isqrt(p)+1):
        y2 = p - x*x; y = isqrt(y2)
        if y*y == y2: return (x, y)
def gm(a, b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])

P = 1000
S = [l for l in primerange(P, 2*P+1) if l % 8 == 5]
print("parent: primes 5 mod 8 in [P,2P], P=1000:", len(S), "labels")
checked = nonint = 0; worst = 0.0
for trial in range(20000):
    p = random.choice(S); pi = gfactor(p)
    if random.random() < 0.5: pi = (pi[0], -pi[1])
    while True:
        g = (random.randint(-60, 60), random.randint(-60, 60))
        if g != (0, 0): break
    Ng = g[0]**2 + g[1]**2
    num = gm(pi, gm(g, g))               # r = pi gamma^2 / Norm gamma
    d = gcd(gcd(num[0], num[1]), Ng)
    alpha = (num[0]//d, num[1]//d); b = Ng//d
    # make lowest terms with positive b; d removes common rational factors
    assert gcd(gcd(alpha[0], alpha[1]), b) == 1
    J = alpha[0]**2 + alpha[1]**2
    assert J == p*b*b
    if gcd(alpha[0], alpha[1]) != 1:
        # common factor of coordinates is a unit/2-power issue only if b even; X.1 says impossible
        raise SystemExit("gcd(x,y) != 1")
    x, y = gm(alpha, alpha)
    for w in (x-y, x+y):
        assert w % 2 == 1 and w != 0
        assert gcd(w, J) == 1
        for l in S: assert w % l != 0, ("parent prime divides", l)
        mx = int(max(q**e for q, e in factorint(abs(w)).items())) if abs(w) > 1 else 1
        worst = max(worst, mx/(2**0.5*p*b*b)); assert abs(w) <= 2**0.5*p*b*b
        assert mx <= 2**0.5*p*Ng*Ng
    lhs = ((x-y)**2 + (x+y)**2)            # = 2 * (forced sine^2 sum) * p^2
    assert lhs == 2*J*J and J*J == p*p*b**4
    checked += 1; nonint += (b > 1)
print(f"{checked} random lifts (|gamma coords| <= 60), {nonint} nonintegral: all pass (a),(b),(c);"
      f" max ratio (largest prime power of x+-y)/(sqrt2 p b^2), b <= Norm gamma = {worst:.4f}")
