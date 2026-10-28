#!/usr/bin/env python3
"""Independent checks of toric.md Prop 1.1, Prop 2.3, Thm 2.1 numbers, alpha_bal(k).
Own code (nothing imported from toric_checks/)."""
import itertools
from fractions import Fraction
from math import comb, factorial

def pts(k, m):
    """lattice points of m*R_k: x in Z^k, sum x = 0, ||x||_1 <= 2m"""
    out = []
    for x in itertools.product(range(-2*m, 2*m+1), repeat=k-1):
        last = -sum(x)
        y = x + (last,)
        if sum(abs(t) for t in y) <= 2*m:
            out.append(y)
    return out

def subsets(k):
    for r in range(1, k):
        for A in itertools.combinations(range(k), r):
            yield A

def lagrange_poly(xs, ys):
    """exact interpolation; return coefficient list (ascending)"""
    n = len(xs)
    coeffs = [Fraction(0)]*n
    for i in range(n):
        # basis polynomial
        b = [Fraction(1)]
        den = Fraction(1)
        for j in range(n):
            if j == i: continue
            b = [Fraction(0)] + b  # multiply by x
            for t in range(len(b)-1):
                b[t] -= xs[j]*b[t+1]
            den *= (xs[i]-xs[j])
        for t in range(n):
            coeffs[t] += ys[i]*b[t]/den
    return coeffs

for k in (3, 4, 5):
    n = k-1
    counts = [len(pts(k, m)) for m in range(0, n+2)]
    co = lagrange_poly(list(range(0, n+2)), counts)
    deg_lead = co[n]
    vol_norm = deg_lead*factorial(n)
    print(f"k={k}: counts m=0..{n+1}: {counts}; Ehrhart coeffs {[str(c) for c in co]}; B^{n} = {vol_norm}; binom(2k-2,k-1) = {comb(2*k-2,k-1)}")
    # check predicted counts for higher m
    for m in range(n+2, n+4):
        pred = sum(co[t]*m**t for t in range(len(co)))
        actual = len(pts(k, m))
        assert pred == actual, (k, m, pred, actual)
    # vertices: for permutation s, solve <m,1_{s1..st}> = -1 t<k, sum m=0 -> m = e_{s_k}-e_{s_1}
    P1 = pts(k, 1)
    for s in itertools.permutations(range(k)):
        sols = [x for x in P1 if all(sum(x[j] for j in s[:t]) == -1 for t in range(1, k))]
        exp = tuple((1 if j == s[-1] else 0) - (1 if j == s[0] else 0) for j in range(k))
        assert sols == [exp], (s, sols)
    # P_B = {<x,1_A> <= 1} equals ||x||_1 <= 2 on sum 0 (check on box)
    for x in itertools.product(range(-3, 4), repeat=k-1):
        y = x + (-sum(x),)
        c1 = all(sum(y[j] for j in A) <= 1 for A in subsets(k))
        c2 = sum(abs(t) for t in y) <= 2
        assert c1 == c2
    # width along integer v: 2(max v - min v)
    for v in itertools.product(range(-2, 3), repeat=k):
        vals = [sum(a*b for a, b in zip(x, v)) for x in P1]
        assert max(vals)-min(vals) == 2*(max(v)-min(v))
    # Vandermonde polytope Perm(-(k-1)/2..) -- only integral when k odd; check ratio ord/m
    print(f"   vertices of P_B are the roots e_a-e_b (checked for all {factorial(k)} perms); widths 2(max v-min v) ok")

# Prop 2.3 / mean width numbers
for k in range(3, 9):
    cuts = [A for A in subsets(k) if 0 not in A]  # one side of each cut not containing 0... use A not containing index 0 -> wait use complement rule
    # each cut {A, A^c}: choose rep containing element 0
    cuts = [A for A in subsets(k) if 0 in A]
    assert len(cuts) == 2**(k-1)-1
    s = sum(len(A)*(k-len(A)) for A in cuts)
    ordV = comb(k, 2)
    # balanced curve: B.C = (2^k-2) kappa ; D_V.C = s kappa >= ordV * mult -> mult <= s/ordV kappa
    alpha_bal = Fraction(2**k-2)*ordV/s
    print(f"k={k}: sum_cuts |A||A^c| = {s} (= k(k-1)2^(k-3) = {k*(k-1)*2**(k-3)}); alpha_bal = {alpha_bal} = 4(1-2^(1-k)) = {4*(1-Fraction(1, 2**(k-1)))}; covering k-1 = {k-1}; delta_bal = 2/alpha-1/2 = {Fraction(2)/alpha_bal-Fraction(1,2)} vs 1/(2^k-2) = {Fraction(1,2**k-2)}")

# Vandermonde at k=5: polytope of prod(x_a-x_b)/(x1..x5)^2 inside 3R_5, order 10
k = 5
mons = {}
for signs in itertools.product((0, 1), repeat=10):
    pass
# expand prod_{a<b}(x_a - x_b) symbolically via dict
poly = {tuple([0]*5): 1}
for a, b in itertools.combinations(range(5), 2):
    new = {}
    for e, c in poly.items():
        ea = list(e); ea[a] += 1
        eb = list(e); eb[b] += 1
        new[tuple(ea)] = new.get(tuple(ea), 0) + c
        new[tuple(eb)] = new.get(tuple(eb), 0) - c
    poly = {e: c for e, c in new.items() if c}
supp = [tuple(t-2 for t in e) for e in poly]
print("V support size", len(supp), "max ||s||_1 =", max(sum(abs(t) for t in s) for s in supp), "(<= 2*3 => in 3R_5)")
# widths of Q_V along 1_A
for A in [(0,), (0, 1)]:
    vals = [sum(s[j] for j in A) for s in supp]
    print("  width_QV(1_A) |A|=%d:" % len(A), max(vals)-min(vals))
# 3-point vandermonde in R_5
p3 = {(0,0,0): 1}
for a, b in itertools.combinations(range(3), 2):
    new = {}
    for e, c in p3.items():
        ea = list(e); ea[a] += 1
        eb = list(e); eb[b] += 1
        new[tuple(ea)] = new.get(tuple(ea), 0) + c
        new[tuple(eb)] = new.get(tuple(eb), 0) - c
    p3 = {e: c for e, c in new.items() if c}
print("3-pt Vandermonde support shifted by -1:", sorted(tuple(t-1 for t in e) for e in p3), "-> l1 norms", sorted(set(sum(abs(t-1) for t in e) for e in p3)))
# lower bound floor(10m/3): V^a * V3^b with 3a+b=m, ord 10a+3b
for m in range(1, 10):
    best = max(10*a+3*(m-3*a) for a in range(0, m//3+1))
    print(f"  m={m}: Vandermonde-product order {best}, floor(10m/3) = {10*m//3}")
