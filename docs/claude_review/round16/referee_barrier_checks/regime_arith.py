# Regime arithmetic for sdc_q: d=q, r=2q, w=sqrt(8 d ln 4r), j=floor(w/2), m'=ceil(4w).
import math, sympy
def par(q):
    d=q; r=2*q; w=math.sqrt(8*d*math.log(4*r)); j=math.floor(w/2); m=math.ceil(4*w)
    return d,r,w,j,m
first=None; first2=None; viol=[]
for q in sympy.primerange(3,200000):
    if q%4!=3: continue
    d,r,w,j,m=par(q)
    if d<316*math.log(4*r): continue
    if first is None: first=(q,w,j,m)
    if first2 is None and 2*m<q: first2=(q,w,j,m)
    if m>8*j+4: viol.append((q,round(w,3),j,m,8*j+4))
print("first regime sdc q:",first," 2m'<q there?",2*first[3]<first[0])
print("first regime q with 2m'<q:",first2)
print("regime q with m'>8j+4 (count, first 5):",len(viol),viol[:5])
print("max m'-8j over regime q<2e5:",max(par(q)[4]-8*par(q)[3] for q in sympy.primerange(3000,200000) if q%4==3))
# general d (not only sdc): m' - 8j
mx=max(math.ceil(4*math.sqrt(8*d*math.log(4*r)))-8*math.floor(math.sqrt(8*d*math.log(4*r))/2) for d in range(1000,20000) for r in [d,2*d,10*d] if d>=316*math.log(4*r))
print("max m'-8j over d in [1000,20000), r in {d,2d,10d} regime:",mx)
