import json, sys, random
from math import log
from fractions import Fraction as Fr
sys.path.insert(0,'.')
import ref_phiN as P
import flint
x=P.x
def gauss(b1,b2):
    def n(v): return v[0]*v[0]+v[1]*v[1]
    if n(b1)>n(b2): b1,b2=b2,b1
    while True:
        m=round(Fr(b1[0]*b2[0]+b1[1]*b2[1], n(b1)))
        b2=(b2[0]-m*b1[0],b2[1]-m*b1[1])
        if n(b2)>=n(b1): return b1
        b1,b2=b2,b1
def lattice(conds):  # list of (alpha,beta,mod): u*beta - alpha*v = 0 mod
    M=1; r=0
    for al,be,m in conds:
        if m==1: continue
        al%=m; be%=m
        rr=al*pow(be,-1,m)%m
        # CRT r mod M, rr mod m
        r=(r+M*((rr-r)*pow(M,-1,m)%m))%(M*m); M*=m
    return M, gauss((M,0),(r,1))
Q1,Q2=P.Qj(0),P.Qj(1); D=P.add(Q2,P.neg(Q1))
for N in map(int,sys.argv[1:]):
    R=P.add(P.Q1 if False else Q1,P.mul(N,D)); a,b,c=P.lifts(R)[0]
    d,ok,cop=P.profile([x[0],x[1],x[2],a,b,c]); w=sum(log(v) for v in d.values())/25/2
    val=lambda C: x[min(C)] if min(C)<=2 else x[3]
    conds=[(val(C).numerator,val(C).denominator,d[C]) for C in P.sides if 3 in C]
    M,v=lattice(conds)
    ua,va=a.numerator,a.denominator
    print(N,'w=%.2f'%w,'#contacts with a:',len(conds),'index/w=%.2f'%(log(M)/w),'log lam1/w=%.3f'%(log(max(abs(v[0]),abs(v[1])))/w),
          'h(a)/w=%.3f'%(log(max(abs(ua),abs(va)))/w),'lam1 = +-(u_a,v_a):',(abs(v[0]),abs(v[1]))==(abs(ua),abs(va)),
          'skew/w=%.2f'%((log(M)/2-log(max(abs(v[0]),abs(v[1]))))/w),flush=True)
