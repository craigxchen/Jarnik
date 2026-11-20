# Wider class of X.1: denominators may contain exterior primes; at q | b charge NOTHING
# (the most favourable choice for survival).  Which nonintegral lifts survive single-n tests?
import sys
from math import gcd, isqrt
from cut import labels, Lcm, gprime, mul, canon, killed
from sympy import factorint
P=int(sys.argv[1]); X=int(sys.argv[2]); five=False
S=labels(P,five); H=2*P
L=Lcm(X,set(S))
R=isqrt(H); cases={}
for p in S:
    a,bb=gprime(p)
    for pi in [(a,bb),(a,-bb)]:
        for u in range(-R,R+1):
            for v in range(-R,R+1):
                Ng=u*u+v*v
                if Ng==0 or Ng>H: continue
                num=mul(pi,mul((u,v),(u,v)))
                d=gcd(gcd(*num),Ng); al=(num[0]//d,num[1]//d); b=Ng//d
                if b>1: cases[(p,b,canon(al))]=1
print("all nonintegral lifts up to sym:",len(cases))
surv=[]
for (p,b,al) in cases:
    Lb=L
    for q in factorint(b):
        while Lb%q==0: Lb//=q
    if killed(al,p,Lb,40)==0: surv.append((p,b,al))
print("not killed by single-n tests n<=40 at X =",X,":",len(surv),surv[:10])
