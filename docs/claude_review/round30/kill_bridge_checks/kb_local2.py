# For l>=7 in P (not dividing root denominators): gain at l iff some residue a mod l has v_l(Ng(a))>0 (complete test).
# Then refine by lifting. For l in {2,3,5}: brute force as in kb_local.
import sys, itertools
from fractions import Fraction as Fr
from kb_gauss import *
from kb_local import vl, normg, badprimes, localmax
import math
def bigprime_max(rts,l,masks=MASKS):
    best=0; arg=0
    hits=[a for a in range(l) if vl(normg(rts,Fr(a),masks),l)>0]
    for a0 in hits:  # lift: search t = a0 + l*b, recursive up to depth 4
        frontier=[a0]; 
        for depth in range(1,5):
            new=[]
            for a in frontier:
                for b in range(l if l<60 else 0):
                    tt=a+b*l**depth
                    new.append(tt)
            if not new: new=frontier
            vals=[(vl(normg(rts,Fr(x),masks),l),x) for x in new]
            m=max(vals)[0]
            if m>best: best,arg=m,max(vals)[1]
            frontier=[x for v,x in vals if v==m][:3]
            if l>=60: break
        v=vl(normg(rts,Fr(a0),masks),l)
        if v>best: best,arg=v,a0
    return best,arg,len(hits)
def C_inf(rts,masks,Gam,verbose=False,P=None):
    P=P or badprimes(rts)
    lp=0; out={}
    for l in P:
        if l==2: v,a=localmax(rts,2,5,8,masks)
        elif l==3: v,a=localmax(rts,3,3,5,masks)
        elif l==5: v,a=localmax(rts,5,1,4,masks)
        else: v,a,h=bigprime_max(rts,l,masks)
        out[l]=(v,a); lp+=v*math.log(l)/2
        if verbose: print(l,v,a,flush=True)
    return Gam/math.exp(lp/2), out
if __name__=="__main__":
    y=[Fr(-9,4),Fr(-2,3),Fr(1,4),Fr(8,3)]
    rts=roots_cube(y)
    c,out=C_inf(rts,MASKS,153335/1296,verbose=True)
    print("C_inf",c)
