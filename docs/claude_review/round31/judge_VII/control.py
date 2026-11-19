# Controls: (i) C missing the factor t^p - t (C = t); (ii) non-p-th-power labels v = i*C*w.
from witness12 import pmul, padd, psub, pgcd, deriv, trim
def ppow(f,e,p):
    r=[1]
    for _ in range(e): r=pmul(r,f,p)
    return r
def R_of(p,T,vsq_list):
    R=T[:]
    for vs in vsq_list: R=pmul(R,psub(pmul(T,T,p),vs,p),p)
    return R
def labels(p):
    seen=set(); out=[]
    for x0 in range(p):
        for x1 in range(p):
            if (x0,x1)==(0,0) or ((-x0)%p,(-x1)%p) in seen: continue
            seen.add((x0,x1)); out.append(trim([x0,x1]))
    return out
for p in (3,7):
    bad1=bad2=0; tot=0
    for s in range(1,25):
        # (i) C = t : v = i^p t^p w^p ; v^2 = -t^{2p} w^{2p}
        C=[0,1]; Cp=ppow(C,p,p)
        T=padd([1],pmul(Cp,[0]*s+[1],p),p) if s%p else padd([1],pmul(Cp,padd([0]*s+[1],[0,1],p),p),p)
        vs=[pmul([p-1],pmul(ppow(C,2*p,p),ppow(w,2*p,p),p),p) for w in labels(p)]
        R=R_of(p,T,vs); bad1+= len(pgcd(R,deriv(R,p),p))>1
        # (ii) non-p-th-power labels with full C: v = i C w, v^2 = -C^2 w^2
        C=[0,p-1]+[0]*(p-2)+[1]; Cp=ppow(C,p,p)
        T=padd([1],pmul(Cp,[0]*s+[1],p),p) if s%p else padd([1],pmul(Cp,padd([0]*s+[1],[0,1],p),p),p)
        vs=[pmul([p-1],pmul(pmul(C,C,p),pmul(w,w,p),p),p) for w in labels(p)]
        R=R_of(p,T,vs); bad2+= len(pgcd(R,deriv(R,p),p))>1
        tot+=1
    print(f"p={p}: s=1..24: non-squarefree R with C=t: {bad1}/{tot}; with non-p-th-power labels: {bad2}/{tot}")
