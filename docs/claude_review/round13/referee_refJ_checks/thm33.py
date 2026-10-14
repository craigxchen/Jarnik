from rlib import *
import itertools
def hilb_box(d, cuts):
    A=[k for k in itertools.product((-1,0,1),repeat=d) if all(abs(sum(e*k[i] for i,e in zip(S,eps)))<=L for S,eps,L in cuts)]
    cum=[]; e=0
    while True:
        ms=[m for m in monos_le(d,e) if max(m)<=2]
        r=fmpz_mat([[evmono(m,a) for m in ms] for a in A]).rank(); cum.append(r)
        if r==len(A): break
        e+=1
    return len(A),[cum[0]]+[cum[i]-cum[i-1] for i in range(1,len(cum))]
def polymul(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): c[i+j]+=x*y
    return c
def trin(s):
    h=[1]
    for _ in range(s): h=polymul(h,[1,1,1])
    return h
def Tsl(s,L):
    h=trin(s); N=2*L+1
    return [max(h[i]-(h[i-N] if i-N>=0 else 0),0) for i in range(len(h))]
def predicted(d,cuts):
    r=[1]; used=0
    for S,eps,L in cuts: r=polymul(r,Tsl(len(S),L)); used+=len(S)
    r=polymul(r,trin(d-used))
    while r and r[-1]==0: r.pop()
    return r
tests=[(6,[((0,1,2),(1,-1,1),1),((3,4),(1,1),0)]),
       (7,[((0,1,2,3),(1,-1,-1,1),1),((4,5,6),(-1,1,1),0)]),
       (7,[((0,1),(1,-1),0),((2,3),(1,1),0),((4,5,6),(1,1,-1),1)]),
       (6,[((0,1,2,3,4,5),(1,1,1,-1,-1,1),2)])]
for d,cuts in tests:
    n,Hh=hilb_box(d,cuts); pr=predicted(d,cuts)
    print(d,cuts,"|A|",n,"H",Hh,"pred",pr,"MATCH" if Hh==pr else "MISMATCH",flush=True)
# overlapping cuts: does (x^3, l_j^{2L+1}) still give |A|? (expect not in general)
d=5; cuts=[((0,1,2),(1,1,1),0),((2,3,4),(1,1,1),0)]
n,Hh=hilb_box(d,cuts); print("overlap",cuts,"|A|",n,"H",Hh)
