"""exact numerators of FINAL in each variable; strip rational roots; report remaining degree."""
from fractions import Fraction as Fr
import qpoly as Q
from ratrec import nullspace
import param3 as P
vals = dict(tau=Fr(-2), a=Fr(3), b=Fr(5, 2), c=Fr(7, 3), mu=Fr(11, 5))
degs = {('tabc','tau'):(7,6),('tabc','a'):(8,8),('tabc','b'):(9,8),('tabc','c'):(7,8),
        ('mabc','mu'):(7,7),('mabc','a'):(6,6),('mabc','b'):(6,5),('mabc','c'):(7,8)}
def exact(f, dn, dd):
    pts = []; k = 0
    while len(pts) < dn + dd + 4:
        x = Fr(k + 1, 1) + Fr(1, 3); k += 1
        try: pts.append((x, f(x)))
        except ZeroDivisionError: pass
    rows = [[x ** i for i in range(dn + 1)] + [-y * x ** j for j in range(dd + 1)] for x, y in pts]
    ns = nullspace(rows); assert len(ns) == 1, len(ns)
    v = ns[0]; return Q.norm(v[:dn + 1]), Q.norm(v[dn + 1:])
if __name__ == "__main__":
     for (name, v), (dn, dd) in degs.items():
        fn = P.final_tabc if name == 'tabc' else P.final_mabc
        vars_ = ['tau', 'a', 'b', 'c'] if name == 'tabc' else ['mu', 'a', 'b', 'c']
        f = lambda x: fn(*[x if w == v else vals[w] for w in vars_])
        num, den = exact(f, dn, dd)
        rr = Q.rational_roots(num)
        rest = num
        mult = {}
        for r in rr:
            while True:
                qq, rem = Q.divmod_(rest, Q.lin(r))
                if rem: break
                rest = qq; mult[r] = mult.get(r, 0) + 1
        print(name, v, 'num deg', Q.deg(num), 'rational roots', {str(k): m for k, m in mult.items()}, 'remaining deg', Q.deg(rest), flush=True)
