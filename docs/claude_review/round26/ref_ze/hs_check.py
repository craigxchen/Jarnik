"""Independent check of ze.md Prop 4.5 inequalities (exact integers), worst case L=3."""
import math, sys
from fractions import Fraction
def m_max(d):
    t=3**d; lo,hi=0,2*d
    while lo<hi:
        mid=(lo+hi+1)//2
        if 8*math.comb(mid+d-1,d)<=t: lo=mid
        else: hi=mid-1
    return lo
L=3
for n in (2966,4096):
    d=n-1; w=math.sqrt(8*d*math.log(4*n)); hw=math.floor(w/2)
    for m in (math.ceil(4*w), m_max(d)):
        hG=2*math.comb(m-1+d,d); A=(3**d+1)//2
        worst=None; ok=True
        b=1  # binom(m-1+e,e) iteratively
        for e in range(0,d+1):
            if e>0: b=b*(m-1+e)//e
            LB=max(3**min(d-e,hw), -(-A//3**e))
            lhs1=max(L+1,LB)*b
            r1=math.log(lhs1)-math.log(hG)
            if lhs1<=hG: ok=False
            cand=[(r1,e,'Vx1')]
            if e<d:
                lhs2=LB*2*b
                if lhs2<=hG: ok=False
                cand.append((math.log(lhs2)-math.log(hG),e,'VxGm'))
            for c in cand:
                if worst is None or c[0]<worst[0]: worst=c
        print(f"n={n} m={m}: all hold={ok}; tightest ln(LHS/RHS)={worst[0]:.4f} at e={worst[1]} {worst[2]}",flush=True)
