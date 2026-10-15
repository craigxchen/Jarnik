# Referee check of Thm B1's product bound (independent): P(e_j survives) <= prod_{q=+-3 (8), Y<q<=X, q not | N} (1-4/m_q)
# compared with (log Y/log X)^2, Y = sqrt(2 p_j); effect of excluding q = p_l (moduli dividing N are not in (1.2)).
from sympy import primerange
from math import log, sqrt
def prod_bound(p, X, excl):
    Y = sqrt(2*p); P = 1.0
    for q in primerange(int(Y)+1, X+1):
        if q % 8 in (3, 5) and q not in excl and q != p:
            m = q - (1 if q % 4 == 1 else -1)
            P *= 1 - 4/m
    return P, (log(Y)/log(X))**2
split = [p for p in primerange(3, 10**6) if p % 4 == 1]
M, A = 30, 3; X = 3*M**A
for p in [101, 1009, 5009, 26993]:
    P0, ref = prod_bound(p, X, set())
    excl = set(split[:3*M])          # a prime set containing the 90 smallest split primes
    P1, _ = prod_bound(p, X, excl)
    print(f"M={M} A={A} X={X} p_j={p}: product (no exclusion) {P0:.4f}, (log Y/log X)^2 {ref:.4f}, product excluding q in first 90 split primes {P1:.4f} (ratio {P1/P0:.2f})")
