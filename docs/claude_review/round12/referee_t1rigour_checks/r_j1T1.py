# Referee: J1's p-adic point at T = 1 via the single direction 1_{V_c}.
# Identity: 1[a.x=c1]1[b.x=c2] = (1/4)(1 + e1 chi_a + e2 chi_b + e1 e2 chi_{a+b}), chi = beta for Sylvester columns.
# On the T=1 widened box A_{Q,lambda} (|k(S_q)| <= lambda m'/2 for q in Q, |k(S_q)| <= j otherwise):
#   |k(V_c)| <= (lambda m' + 4j)/4, so the family (b, a+b, 1_{V_c}) costs <= 4j + 2 floor((lambda m' + 4j)/4).
import itertools, math
from fractions import Fraction
# exact identity check for k = 3..6
fails = 0
for k in range(3, 7):
    pts = list(itertools.product((0,1), repeat=k))
    dot = lambda a, x: sum(ai*xi for ai, xi in zip(a, x)) % 2
    vecs = [a for a in pts if any(a)]
    for a in vecs:
        for b in vecs:
            if a == b: continue
            ab = tuple((x+y) % 2 for x, y in zip(a, b))
            for c1 in (0,1):
                for c2 in (0,1):
                    e1 = (-1)**c1; e2 = (-1)**c2
                    for x in pts:
                        lhs = 1 if (dot(a,x) == c1 and dot(b,x) == c2) else 0
                        chi = lambda v: (-1)**dot(v, x)
                        rhs = Fraction(1 + e1*chi(a) + e2*chi(b) + e1*e2*chi(ab), 4)
                        if lhs != rhs: fails += 1
print("indicator identity failures (k=3..6):", fails)
mmax = {12: 2631, 13: 5260, 14: 10518, 15: 21034}   # exact values (independently: n=4096 -> 2631 in r_params)
for k in (12, 13, 14, 15):
    n = 2**k; r = n - 1; w = math.sqrt(8*(n-1)*math.log(4*r)); j = int(w//2); lo = math.ceil(4*w)
    for lam in (Fraction(1,4), Fraction(1,2), Fraction(1), Fraction(2)):
        ok = [m for m in range(lo, mmax[k]+1) if 4*j + 2*math.floor((lam*m + 4*j)/4) < m]
        print(f"k={k} n={n} j={j} window=[{lo},{mmax[k]}] lambda={lam}: family (b,a+b,1_Vc) cheap for m' in",
              (ok[0], ok[-1]) if ok else None, "; author's coordinate family cheap iff m' >", n//2 + 4*j)
