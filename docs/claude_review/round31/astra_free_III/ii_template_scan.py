"""Nonexhaustive discovery scan for rational D=8 configurations; integer moments."""
from collections import defaultdict, Counter
from itertools import combinations
from sys import argv

def m(x,y):
    return (y,x*y,3*x*x*y-y**3)

def plus(a,b): return tuple(x+y for x,y in zip(a,b))

def pairs(vec):
    pp=defaultdict(list)
    for i,j in combinations(range(len(vec)),2): pp[plus(vec[i],vec[j])].append((i,j))
    return pp

def fibres(vec):
    f=Counter({(0,0,0):1})
    for v in vec:
        f+=Counter({plus(k,v):n for k,n in f.items()})
    return max(f.values())

def scan(h,w):
    roots=[(x,y) for x in range(-w,w+1) for y in range(1,h+1)]
    vec=[m(*r) for r in roots]
    pc=pairs(vec)
    quad={frozenset(a+b) for group in pc.values() for a,b in combinations(group,2)}
    print("points, parallelograms",len(roots),len(quad), flush=True)
    by2=defaultdict(set)
    for e in quad:
        for p in combinations(sorted(e),2):
            by2[p].add(e)
    vs=set(quad)
    for step in range(4):
        nv=set()
        for s in vs:
            if len(s)==8:continue
            for q in combinations(sorted(s),2):
                for e in by2[q]:
                    z=s|e
                    if len(z)<=8 and len(z)>len(s):nv.add(z)
        vs=nv
        if not vs: break
        stats={}
        for s in vs:
            if len(s)!=8:continue
            n=sum(frozenset(c) in quad for c in combinations(s,4))
            if n<5:continue
            sol=sorted(roots[i] for i in s)
            f=fibres([m(*r) for r in sol])
            if f not in stats or n>stats[f][0]:stats[f]=(n,sol)
        print("step, candidates and records",step,len(nv),stats,flush=True)

if __name__=="__main__":scan(int(argv[1]),int(argv[2]))
