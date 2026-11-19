# integer quintuples with e4=0, e3=-e1*e2, -e1*e5 a nonzero square, distinct nonzero, b_ij distinct nonzero
import itertools, math, sys
from sympy import integer_nthroot
R=int(sys.argv[1])
vals=[v for v in range(-R,R+1) if v]
found=0
for c in itertools.combinations(vals,5):
    # e4 = sum of products of 4 = e5*sum(1/x)
    e5=1
    for v in c: e5*=v
    if sum(e5//v for v in c)!=0: continue
    e1=sum(c); e2=sum(a*b for a,b in itertools.combinations(c,2)); e3=sum(a*b*d for a,b,d in itertools.combinations(c,3))
    if e3!=-e1*e2: continue
    t=-e1*e5
    if t<=0: continue
    r,ex=integer_nthroot(t,2)
    if not ex: continue
    b=[a*bb*(e1-a-bb) for a,bb in itertools.combinations(c,2)]
    if 0 in b or len(set(b))<10: continue
    print(c,'s1..s5',e1,e2,e3,0,e5,'l=',r); found+=1
print('found',found)
