"""Referee ref2_plane: the curve {[y]: [y_1^2:...:y_5^2] in Lambda} of II.7.
Lambda is cut out by the three rows of the 3x5 matrix with columns (1, A_i)/D_i.  The curve is the complete
intersection of three diagonal quadrics in P^4.  Check (i) all 3x3 minors nonzero for a sample configuration
(=> smooth CI of three quadrics => genus 5 by adjunction K=O(1), deg 8); (ii) smoothness via Groebner basis
of quadrics + 3x3 minors of the Jacobian; (iii) point counts over F_p versus the Weil bound for g=5 and
compare with g=4 / g=6 bounds as a weak sanity check; (iv) Riemann-Hurwitz recomputation.
"""
import sympy as sp, itertools, random
random.seed(3)
# sample: five rational points, no three collinear, and arbitrary nonzero D_i (in II.7 D_i=prod_j c_ij)
A = [(0,0),(1,0),(0,1),(2,3),(5,-1)]
D = [1,2,-3,5,7]
cols = [sp.Matrix([1, a[0], a[1]])/d for a, d in zip(A, D)]
M = sp.Matrix.hstack(*cols)
minors = [M[:, list(c)].det() for c in itertools.combinations(range(5), 3)]
print("all 3x3 minors nonzero:", all(m != 0 for m in minors))
y = sp.symbols('y1:6')
Q = [sum(M[k, i]*y[i]**2 for i in range(5)) for k in range(3)]
Q = [sp.expand(q*sp.ilcm(*[sp.fraction(sp.nsimplify(M[k,i]))[1] for i in range(5)])) for k, q in enumerate(Q)]
J = sp.Matrix([[sp.diff(q, v) for v in y] for q in Q])
jm = [J[:, list(c)].det() for c in itertools.combinations(range(5), 3)]
G = sp.groebner(Q + jm, *y, order='grevlex')
# check that every variable has a power in the ideal (=> projective zero set empty => smooth)
ok = all(any(sp.Poly(g, *y).is_monomial and sp.Poly(g, *y).free_symbols_in_domain == set() and
             sp.Poly(g, *y).degree(v) == sp.Poly(g, *y).total_degree() for g in G.exprs) or
         G.reduce(v**12)[1] == 0 for v in y)
print("Jacobian criterion: y_i^12 in ideal for all i ->", ok)

# point counts over F_p
def count(p):
    Qp = [[int(sp.Poly(q, *y).coeff_monomial(v**2)) % p for v in y] for q in Q]
    # enumerate projective points
    n = 0
    for pt in itertools.product(range(p), repeat=5):
        # normalise: first nonzero coordinate =1
        nz = [c for c in pt if c]
        if not nz or nz[0] != 1: continue
        if all(sum(Qp[k][i]*pt[i]*pt[i] for i in range(5)) % p == 0 for k in range(3)):
            n += 1
    return n
for p in [11, 13, 17, 19]:
    bad = any(m.p % p == 0 or m.q % p == 0 for m in minors) or any(d % p == 0 for d in D)
    if bad: continue
    c = count(p)
    print("p=%d  #C(F_p)=%d  p+1=%d  |dev|=%d  10*sqrt(p)=%.1f" % (p, c, p+1, abs(c-p-1), 10*p**0.5))
print("Riemann-Hurwitz: degree 16 cover of P^1, 5 branch points each with 8 ramification points of index 2:",
      "2g-2 =", 16*(-2) + 5*8, "=> g =", (16*(-2)+5*8+2)//2)
