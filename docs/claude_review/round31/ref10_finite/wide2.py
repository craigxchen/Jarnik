from cut import labels, Lcm, killed
from math import gcd
S=labels(100,False)
for (p,b,al,qs) in [(157,185,(-2307,-226),(5,37)),(173,145,(-1838,-509),(5,29))]:
    lo,hi=439,2000
    def ok(X):
        L=Lcm(X,set(S))
        for q in qs:
            while L%q==0: L//=q
        return killed(al,p,L,60)!=0
    while hi-lo>1:
        mid=(lo+hi)//2
        if ok(mid): hi=mid
        else: lo=mid
    print(p,b,al,"single-n kill cutoff (zero charge at q|b):",hi, "kill n at hi:", end=" ")
    L=Lcm(hi,set(S))
    for q in qs:
        while L%q==0: L//=q
    print(killed(al,p,L,60))
