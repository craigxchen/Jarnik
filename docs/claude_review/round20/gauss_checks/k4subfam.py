"""Exact rank of the rational 2-parameter subfamily (X,m) x curve parameter s -> T_4 (similarity classes) and -> P(X_4).
Subfamily: r123=i, r12=u+i(X-1), r13=u+i(x_b-1), r23=u+i(x_c-1), x_b=X(1-X)/(m^2+X), x_c=1-X-x_b, u=m x_b.
Rank 3 onto T_4 => the subfamily's lattice-point tuples have Zariski-dense similarity classes (rational params are dense in Q^2 x Q)."""
from fractions import Fraction as F
import k4dominance as K
def run(Xv, mv, sv):
    X, m, s = K.var(0, Xv), K.var(1, mv), K.var(2, sv)
    xb = X*(1-X)/(m*m+X); xc = 1 - X - xb; u = m*xb
    A, B, Cc = (u, X-1), (u, xb-1), (u, xc-1); rho = (K.D(0), K.D(1))
    roots = {(1,2,3): K.C(0,1), (1,2): K.C(*A), (1,3): K.C(*B), (2,3): K.C(*Cc),
             (1,): K.sing([rho,A,B]), (2,): K.sing([rho,A,Cc]), (3,): K.sing([rho,B,Cc])}
    assert all(K.phi(p,q).v == 0 for p,q in [(A,B),(A,Cc),(B,Cc)])
    H = {T: K.C(s, K.D(0)) - r for T, r in roots.items()}
    z = []
    for x in range(4):
        acc = K.C(1, 0)
        for T in H: acc = acc * (H[T] if x in T else H[T].conj())
        z.append(acc)
    t1 = [(z[x]*z[0].conj()) for x in (1,2,3)]; t1 = [w.im/w.re for w in t1]
    t2 = [w.im/w.re for w in z]
    return K.rank([t.g[:3] for t in t1]), K.rank([t.g[:3] for t in t2])
for p in [(F(-3),F(-3),F(7)), (F(-2),F(1,3),F(9)), (F(5,2),F(2),F(4)), (F(-7,3),F(5,4),F(-11))]:
    print("(X,m,s)=", [str(x) for x in p], " rank onto T_4 =", run(*p)[0], "(3 = dominant)", " rank onto P(X_4) =", run(*p)[1], "(max 3 < 4)")
