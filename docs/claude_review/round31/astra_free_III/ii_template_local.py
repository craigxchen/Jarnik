"""Discovery of local contractions in integral degree-eight specializations."""
from collections import defaultdict
from itertools import combinations
from math import gcd,lcm,prod
from sympy import factorint
from ii_template_scan import m,pairs
ROOTS=[(-12,4),(-8,6),(-3,16),(0,10),(0,12),(3,16),(8,6),(12,4)]

def mul(z,w):
    return z[0]*w[0]-z[1]*w[1], z[0]*w[1]+z[1]*w[0]

def get_solutions(roots):
    vec = [m(*r) for r in roots]
    sols=[]
    for grp in pairs(vec).values():
        for ij,kl in combinations(grp,2):
            sol=([ (j,1) for j in ij ]+[ (j,-1) for j in kl ])
            sols.append(sol)
    return sols

def value(roots,sol,n):
    z=(1,0)
    for j,s in sol:
        x,y=roots[j]
        z=mul(z,(n-x,-s*y))
    return z

def data(roots,sols,n):
    fs=[(n-x)**2+y*y for x,y in roots]
    T=prod(fs)
    d=[]; g=[]
    for s in sols:
        z=value(roots,s,n)
        k=gcd(*z)
        g.append(k)
        d.append(prod(fs[j] for j,_ in s)//k**2)
    return T,d,g

if __name__=='__main__':
    import sys
    roots=ROOTS
    sols=get_solutions(roots)
    b=[abs(value(roots,s,0)[1]) for s in sols]
    L=1; pows={}
    for n in range(0,int(sys.argv[1])):
        T,d,g=data(roots,sols,n)
        x=T//lcm(*d)
        if L%x:
            for p,e in factorint(x).items():
                if p**e > pows.get(p,(0,0))[0]:pows[p]=(int(p**e),n)
            L=lcm(L,x)
    print('gains',pows,'L',L,'sqrtL',int(L**.5),'b',b)
    print('pass',[x for x in b if L>256*x**4])
