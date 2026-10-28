"""m-czech (forum, July 2026, via formal-conjectures PR #7218) five-divisor Pell family for Erdos #887, checked exactly.
n=(x^2-1)(x^2-4)(x^2-16)(x^2-49), z^2=56x^2+65. Reports divisors of n in (sqrt n, sqrt n + C n^{1/4}), the actual C,
and the split-torus data of the 5-tuple: primitive u, H=(max/min)lcm(u), dist=(max-min)/min, gamma=log H/(-log dist)."""
import sympy as sp
from math import gcd, lcm, isqrt, log
x, z = sp.symbols('x z')
n_poly = sp.expand((x**2-1)*(x**2-4)*(x**2-16)*(x**2-49))
f5 = x**4 - 35*x**2 + 56 + 2*x*z
f5b = x**4 - 35*x**2 + 56 - 2*x*z
prod = sp.expand((f5*f5b).subs(z**2, 56*x**2+65))
prod = sp.expand(sp.expand(f5*f5b).subs(z**2, 56*x**2+65))
print("f5*f5bar mod (z^2-56x^2-65) =", sp.factor(prod), "; equals n?", sp.expand(prod - n_poly) == 0)
# Pell solutions of z^2 - 56 x^2 = 65
sols = set()
for X in range(0, 200000):
    t = 56*X*X+65; r = isqrt(t)
    if r*r == t: sols.add((X, r))
base = sorted(sols)[:12]
print("small solutions (x,z):", base)
# fundamental unit of Z[sqrt56]: 15+2sqrt56 (15^2-56*4=1)
def nxt(X, Z): return (15*X + 2*Z, 15*Z + 56*2*X)
allsol = set()
for (X, Z) in base:
    a, b = X, Z
    for _ in range(8):
        allsol.add((a, b)); a, b = nxt(a, b)
for (X, Z) in sorted(allsol)[:40]:
    if X < 10: continue
    n = (X*X-1)*(X*X-4)*(X*X-16)*(X*X-49)
    assert Z*Z == 56*X*X+65
    s = isqrt(n)
    q4 = n**0.25
    cand = set()
    import itertools
    lin = [X+1, X-1, X+2, X-2, X+4, X-4, X+7, X-7]
    for r in range(0, 9):
        for c in itertools.combinations(lin, r):
            d = 1
            for v in c: d *= v
            cand.add(d)
    cand.add(X**4 - 35*X*X + 56 + 2*X*Z)
    ds = sorted(d for d in cand if n % d == 0 and s < d and d <= s + 40*q4 + 1)
    if len(ds) < 5: print(X, "fewer than 5:", len(ds)); continue
    Cact = (max(ds) - s) / q4
    d5 = ds[:5]
    g = 0
    for d in d5: g = gcd(g, d)
    u = [d//g for d in d5]
    L = 1
    for v in u: L = lcm(L, v)
    H = L * max(u) / min(u)
    dist = (max(u) - min(u)) / min(u)
    print(f"x={X}: #div in window(C=40)={len(ds)}, C_actual={Cact:.3f}, gcd={g}, log n={log(n):.1f}, log L/log n={log(L)/log(n):.4f}, gamma={log(H)/(-log(dist)):.4f}")
