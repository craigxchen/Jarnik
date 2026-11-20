import numpy as np, itertools, sys
from sympy import Matrix
exec(open('paley_half.py').read().split('for m in')[0])
for m,tmax in [(11,6),(19,9),(23,11)]:
    s=S(m); miss=[i for i in range(m) if s[i,0]!=1]
    T=np.vstack([np.ones(m,dtype=int),s[miss]]); r=np.linalg.matrix_rank(T)
    for t in range(1,tmax+1):
        vecs=[]; maxh=0
        combos=list(itertools.combinations(range(m),t))
        for c0 in range(0,len(combos),20000):
            cs=np.array(combos[c0:c0+20000])
            rk=np.linalg.matrix_rank(T[:,cs].transpose(1,0,2).astype(float))
            for idx in np.nonzero(rk<t)[0]:
                c=cs[idx]; ns=Matrix(T[:,c]).nullspace()
                if len(ns)==1:
                    v=ns[0]; v=v*np.lcm.reduce([int(x.q) for x in v]); x=np.zeros(m,dtype=object); x[list(c)]=list(v)
                    vecs.append(x); maxh=max(maxh,max(abs(int(z)) for z in v))
        if vecs:
            d=Matrix(np.array(vecs,dtype=object).tolist()).rank()
            print(f"m={m} t={t}: circuits {len(vecs)}, span dim {d} of ker dim {m-r}, max coeff {maxh}"); sys.stdout.flush()
