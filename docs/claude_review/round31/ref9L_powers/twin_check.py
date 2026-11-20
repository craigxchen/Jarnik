#!/usr/bin/env python3
"""Twin with prime powers: exact check of the triangle division, and random checks
of the extraction inequalities (Cauchy-Schwarz steps)."""
import itertools, random
from fractions import Fraction as F
from sympy.ntheory.modular import crt

def vp(n,p):
    k=0
    while n%p==0: n//=p; k+=1
    return k

def twin(mods, classes, lo, hi):
    q=1
    for m in mods.values(): q*=m
    ys=[]
    for choice in itertools.product(*[classes[p] for p in mods]):
        r,_=crt(list(mods.values()), list(choice))
        r=int(r)
        y=r+((lo-r)//q)*q
        while y<hi:
            if y>=lo: ys.append((y,choice))
            y+=q
    ys.sort()
    return q,ys

def check(mods,classes,lo,hi):
    q,ys=twin(mods,classes,lo,hi)
    plist=list(mods)
    lowagree=0; n=0
    for (y1,c1),(y2,c2),(y3,c3) in itertools.permutations(ys,3):
        h=lambda ca,cb: eval_h(ca,cb)
        def eval_h(ca,cb):
            r=1
            for k,p in enumerate(plist):
                if ca[k]==cb[k]: r*=mods[p]
            return r
        h12,h23,h13=eval_h(c1,c2),eval_h(c2,c3),eval_h(c1,c3)
        assert (y2-y1)%h12==0 and (y3-y2)%h23==0 and (y3-y1)%h13==0
        t12,t23,t13=(y2-y1)//h12,(y3-y2)//h23,(y3-y1)//h13
        assert 0 not in (t12,t23,t13)
        G=1
        for k,p in enumerate(plist):
            if c1[k]==c2[k]==c3[k]: G*=mods[p]
        A12,A23,A13=h12//G,h23//G,h13//G
        assert t12*A12+t23*A23==t13*A13
        for k,p in enumerate(plist):
            e=vp(mods[p],p)
            vs=[vp(A,p) for A in (A12,A23,A13)]
            if c1[k]==c2[k]==c3[k]: assert vs==[0,0,0]
            else:
                assert sorted(vs)==[0,0,e], (p,vs)   # exactly one block, full exponent
                # disagreeing pairs: p^v | t where v = v_p(r-r')
                v=vp(classes[p][0]-classes[p][1],p)
                for (ca,cb,t) in ((c1,c2,t12),(c2,c3,t23),(c1,c3,t13)):
                    if ca[k]!=cb[k]:
                        assert vp(t,p)>=v
                        if v>0: lowagree+=1
        n+=1
    print(f"q={q} marks={len(ys)} ordered triples={n} disagreeing-pair events absorbing lower agreement={lowagree}")

def ineq_checks(trials=20000):
    rnd=random.Random(1)
    worst=None
    for _ in range(trials):
        e=rnd.randint(1,5); b=rnd.randint(2,12)
        a=[rnd.randint(0,e) for _ in range(b)]
        al=[F(a.count(x),b) for x in range(e+1)]
        delta=sum(x*(1-x) for x in al)
        Pexcl=sum(x*x*(1-x) for x in al)
        assert Pexcl*e>=delta*delta
        ED=sum(abs(x-y) for x in a for y in a)*F(1,b*b)
        assert ED<=e*delta
        # threshold capacity: sum over unordered pairs |a_i-a_j| <= e*b^2/4
        assert 4*sum(abs(x-y) for x,y in itertools.combinations(a,2))<=e*b*b
        # direct check P(a1=a2!=a3) with replacement
        P=F(sum(1 for x in a for y in a for z in a if x==y!=z),b**3)
        assert P==Pexcl
    print("inequality checks passed:",trials)

if __name__=="__main__":
    mods={3:9,5:25,7:49,11:11}
    classes={3:[1,4],5:[1,6],7:[2,3],11:[0,5]}   # mod 9 and mod 25 classes agree mod 3, mod 5
    check(mods,classes,0,3*9*25*49*11)
    mods={2:8,3:27,5:5}
    classes={2:[1,5],3:[2,11],5:[1,2]}
    check(mods,classes,-1080,2160)
    ineq_checks()
