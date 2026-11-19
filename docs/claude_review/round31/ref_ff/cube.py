import sympy as sp, itertools
exec(open('canon_helper.py').read())
t=sp.symbols('t',real=True); R=sp.Rational
q=[R(-9,4),R(-2,3),R(1,4),R(8,3)]; e=[1,-1,-1,1]
r=[e[j]/q[j]+sp.I*q[j] for j in range(4)]
roots=r+[-x for x in r]           # indices 0..7
allv=roots+[sp.conjugate(x) for x in roots]
assert len(set(sp.nsimplify(x) for x in allv))==16
subsets=[(),tuple(range(8)),(0,1,4,5),(2,3,6,7),(0,2,4,6),(1,3,5,7),(0,1,2,3),(4,5,6,7)]
polys=[]
for S in subsets:
    P=sp.expand(sp.prod([t-(sp.conjugate(roots[k]) if k in S else roots[k]) for k in range(8)]))
    polys.append(sp.Poly(P,t))
N=sp.expand(polys[0].as_expr()*sp.conjugate(polys[0].as_expr()))
for P in polys: assert sp.expand(P.as_expr()*sp.conjugate(P.as_expr())-N)==0
print('degrees of differences',sorted(set(sp.Poly(sp.expand(A.as_expr()-B.as_expr()),t).degree() for A,B in itertools.combinations(polys,2))))
# sign code w.r.t. upper-half-plane reps
signs=[]
for S in subsets:
    row=[]
    for k in range(8):
        used=sp.conjugate(roots[k]) if k in S else roots[k]
        row.append(1 if sp.im(used)>0 else -1)
    signs.append(row)
print(canon_signs(signs))
