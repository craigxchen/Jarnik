from sympy import prime
from math import log, ceil, exp
for k in [5,6,7,8]:
    m=k-1; V=2**m-m-2; per=k-3; r=ceil(V/per); pr=prime(r)
    # sum bound
    sm=sum(log(prime(ceil(j/per))) for j in range(1,V+1))
    out=[]
    for eta in [0.1,0.25,0.5]:
        w1=log(pr)/(1+eta); w2=sm/(V*(1+eta)); w=max(w1,w2)
        out.append("eta=%.2f: e^w>=%.1f (max-bound %.1f, sum-bound %.1f); minor>=10^%.1f"%(eta,exp(w),exp(w1),exp(w2),(2**(k-3)-1)*(1-eta)*w/log(10)))
    print(k,V,per,r,pr,"|"," | ".join(out))
