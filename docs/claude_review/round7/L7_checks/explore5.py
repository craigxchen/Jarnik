from fractions import Fraction as Fr
import qpoly as Q, param2 as P2, ratrec2 as RR
base = dict(b=Fr(3), tau=Fr(-2), mu=Fr(5, 2), sp=Fr(7, 3))
for var in ['sp', 'b', 'tau', 'mu']:
    def f(x):
        a = dict(base); a[var] = x
        return P2.full(a['b'], a['tau'], a['mu'], a['sp'])[1]
    p, q = RR.reconstruct(f, 24)
    print(var, 'deg num', Q.deg(p), 'deg den', Q.deg(q), 'rational roots num', sorted(Q.rational_roots(p)), flush=True)
