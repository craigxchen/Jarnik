# exact check of the degenerate limit of the rank-5 near-hits at D=8: roots snapped to 8th roots of unity
import sympy as sp, re, ast, pickle, numpy as np, itertools
t=sp.symbols('t'); w=sp.sqrt(2)/2
pts=[w+sp.I*w, -w+sp.I*w, -w-sp.I*w, w-sp.I*w]
num=[complex(p) for p in pts]
for log,pk in (("ext8_b.txt","t8_b.pkl"),("ext8_a.txt","t8_a.pkl")):
    T=pickle.load(open(pk,'rb'))
    for line in open(log):
        if not line.startswith("HIT"): continue
        m=re.match(r"HIT type (\d+) F \[([0-9 ]+)\] degen (\S+) (.*)",line.strip())
        idx=int(m.group(1)); F=list(map(int,m.group(2).split())); q=ast.literal_eval(m.group(4))
        E=np.array(T[idx]).T; sets=[list((E[i]!=E[0]).astype(int)) for i in range(len(E))]+[F]
        qs=[pts[int(np.argmin([abs(x-p) for p in num]))] for x in q]
        polys=[]
        for S in sets:
            z=sp.expand(sp.prod([(t-sp.conjugate(qs[k])) if S[k] else (t-qs[k]) for k in range(8)]))
            polys.append(sp.Poly(sp.nsimplify(z),t))
        distinct=len(set(tuple(sp.nsimplify(c) for c in p.all_coeffs()) for p in polys))
        degs=[sp.degree(sp.simplify(polys[i].as_expr()-polys[j].as_expr()),t) for i,j in itertools.combinations(range(len(polys)),2)]
        print(log,"type",idx,"snapped roots",[str(x) for x in qs],"distinct points",distinct,"max pair deg",max(d for d in degs if d!=-sp.oo) if degs else None, "pair degs",degs)
