import numpy as np, sys
from hlib import *
rng=np.random.default_rng(1)
for q in [7,11,19,23,31,43,47,59]:
    P,Pp=patterns(q); ch=chi_vec(q); bad=0
    for trial in range(400):
        f=rng.integers(-3,4,size=q)*(rng.random(q)<rng.random())
        lam=np.zeros(q+1,dtype=np.int64); lam[:q]=f; lam[q]=-f.sum()
        if sigma_direct(lam,q,P,Pp)!=sigma_formula(f,q,ch): bad+=1
    # chi*chi identity
    cc=conv_chi(ch.copy(),q,ch); ok=(cc[0]==1-q) and all(cc[1:]==1)
    print(f"q={q} formula mismatches={bad}/400 chi*chi=1-q*delta0: {ok}",flush=True)
