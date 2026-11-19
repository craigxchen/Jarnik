import sys, itertools, ast
from f4lib import *
from clusters import span
from test_vi3 import coords, parity
d=int(sys.argv[1]); B=int(sys.argv[2],0); T=int(sys.argv[3],0); sides=ast.literal_eval(sys.argv[4])
basis=[1<<i for i in range(d)]; V=span(basis); co=coords(basis)
forms=[l for l in range(1,1<<d) if parity(l&co[B])]
Y=(T,B); pts=[]
for l,c in zip(forms,sides):
    pts.append(prod4([add4(Y,(v,0)) for v in V if parity(l&co[v])==c]))
G=pts[0]
for z in pts[1:]: G=gcd4(G,z)
w=[]
for z in pts:
    q,r=divmod4(z,G); assert d4(r)<0; w.append(monic4(q))
D=d4(w[0]); N=mul4(w[0],conj4(w[0]))
ok=is_real(N) and all(mul4(x,conj4(x))==N for x in w) and all(d4(x)==D for x in w)
md=max(d4(add4(a,b)) for a,b in itertools.combinations(w,2))
print(f"d={d} m={len(w)} distinct={len(set(w))==len(w)} norms_equal={ok} D={D} parent={2*D} maxdiff={md} D/2={D/2} critical={md<=D//2} gcd={d4(G)}")
