# Judge's own check of the VII.1 finite-zero separation step: R = prod_v (T - v) squarefree in F_p[t].
# Uses q = p (labels x_j in F_p), so every pair {v,-v} gives a real factor T^2 - v^2.
import sys, itertools
sys.setrecursionlimit(10000)
from witness12 import pmul, padd, psub, pgcd, deriv, pdivmod, trim
def ppow(f,e,p):
    r=[1]
    for _ in range(e): r=pmul(r,f,p)
    return r
def run(p,model,L,control=False):
    k=2
    base=[0,p-1]+[0]*(p-2)+[1]          # t^p - t
    Delta=[1] if model=='const' else [1,0,2]   # ramified: Delta = 2t^2+1 (p=3)
    C=pmul(Delta,base,p) if model!='const' else base
    degC=len(C)-1; s=L-p*degC
    Cp=ppow(C,p,p)
    if s%p: T=padd([1],pmul(Cp,[0]*s+[1],p),p)
    else:   T=padd([1],pmul(Cp,padd([0]*s+[1],[0,1],p),p),p)
    if control: T=padd([1],[0]*s+[1],p)  # control: T without the C^p factor
    # eta^2 = d (d=-1 constant, d=Delta ramified); v = eta^p C^p w^p, v^2 = d^p C^{2p} w^{2p}
    d=[(p-1)] if model=='const' else Delta
    dp=ppow(d,p,p); C2p=pmul(Cp,Cp,p)
    R=T[:]; facs=[T]
    labels=[(x0,x1) for x0 in range(p) for x1 in range(p) if (x0,x1)!=(0,0)]
    seen=set()
    for x0,x1 in labels:
        if ((-x0)%p,(-x1)%p) in seen: continue
        seen.add((x0,x1))
        w=trim([x0,x1])
        wp=ppow(w,2*p,p)
        f=psub(pmul(T,T,p),pmul(pmul(dp,C2p,p),wp,p),p)
        facs.append(f); R=pmul(R,f,p)
    g=pgcd(R,deriv(R,p),p)
    gC=pgcd(R,C,p)
    print(f"p={p} model={model} L={L} s={s} control={control}: deg R={len(R)-1}, gcd(R,R') deg={len(g)-1}, gcd(R,C) deg={len(gC)-1}, R mod C ={pdivmod(R,C,p)[1]}")
run(3,"const",11); run(3,"const",12); run(3,"const",13); run(3,"const",15)
run(3,"ram",17); run(3,"ram",18); run(3,"ram",22); run(3,"ram",24)
run(3,'const',11,control=True)
run(7,"const",51); run(7,"const",57); run(7,"const",63)
