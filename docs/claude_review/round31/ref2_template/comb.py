# Independent combinatorial checks for II.2 (referee ref2_template)
from fractions import Fraction as F
from itertools import combinations, product
from collections import defaultdict, Counter
import sys

C1=[(-12,4),(-8,6),(-3,16),(0,10),(0,12),(3,16),(8,6),(12,4)]
C2=[(-23,4),(-11,10),(-8,16),(2,16),(5,10),(5,12),(13,6),(17,4)]
fr=(F(-9,4),F(-2,3),F(1,4),F(8,3))
CUBE=[(x,abs(t)) for t,e in zip(fr,(1,-1,-1,1)) for x in (e/t,-e/t)]

def gmul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def polymul(P,Q):
    R=[(0,0)]*(len(P)+len(Q)-1)
    for i,a in enumerate(P):
        for j,b in enumerate(Q):
            c=gmul(a,b); R[i+j]=(R[i+j][0]+c[0],R[i+j][1]+c[1])
    return R
def Bpoly(roots,v):
    # coefficients low->high of prod_{v=+1}(t-rho) prod_{v=-1}(t-rhobar)
    P=[(1,0)]
    for j,s in enumerate(v):
        if s==0: continue
        x,y=roots[j]
        sg=1 if s>0 else -1
        for _ in range(abs(s)): P=polymul(P,[(-x,-sg*y),(1,0)])
    return P
def kernel_vectors(roots,allow2=False):
    D=len(roots); out=[]
    rng=(-2,-1,0,1,2) if allow2 else (-1,0,1)
    for v in product(rng,repeat=D):
        if sum(abs(a) for a in v)!=4: continue
        P=Bpoly(roots,v)   # degree = sum|v| =4 (with multiplicity)
        if all(P[k][1]==0 for k in (1,2,3)) and P[0][1]!=0:
            # first nonzero entry positive -> up to reversal
            fnz=next(a for a in v if a)
            if fnz>0: out.append(v)
    return out
def kernel_vectors_allIm(roots,allow2=False):
    D=len(roots); out=[]
    rng=(-2,-1,0,1,2) if allow2 else (-1,0,1)
    for v in product(rng,repeat=D):
        if sum(abs(a) for a in v)!=4: continue
        P=Bpoly(roots,v)
        if all(P[k][1]==0 for k in (1,2,3)):
            fnz=next(a for a in v if a)
            if fnz>0: out.append(v)
    return out

def mom(r):
    x,y=r; return (y,x*y,3*x*x*y-y**3)

def analyse(name,roots,selected=None):
    print('=====',name)
    dirs=kernel_vectors_allIm(roots)
    print('ternary weight-4 directions (up to sign):',len(dirs))
    for k,v in enumerate(dirs):
        P=Bpoly(roots,v)
        print(k,v,'pos',[j for j in range(8) if v[j]>0],'neg',[j for j in range(8) if v[j]<0],'b=',P[0][1])
    all4=kernel_vectors_allIm(roots,True)
    print('all weight-4 integer kernel vectors incl. entries 2 (up to sign):',len(all4),' midpoints:',[v for v in all4 if 2 in map(abs,v)])
    # pair-sum coincidences
    h=[mom(r) for r in roots]
    ps=defaultdict(list)
    for i,j in combinations(range(8),2): ps[tuple(a+b for a,b in zip(h[i],h[j]))].append((i,j))
    print('pair-sum classes:',[g for g in ps.values() if len(g)>1])
    # fibres
    fib=defaultdict(list)
    for m in range(256):
        fib[tuple(sum(h[i][c] for i in range(8) if m>>i&1) for c in range(3))].append(m)
    hist=Counter(len(f) for f in fib.values())
    print('fibre histogram',sorted(hist.items()),'max',max(hist))
    # cliques in {0,1}^8 with differences in +-dirs (selected or all)
    def maxclique(dset,distinct=False):
        S=set()
        for v in dset: S.add(v); S.add(tuple(-a for a in v))
        verts=[tuple((m>>i)&1 for i in range(8)) for m in range(256)]
        adj={u:[w for w in verts if tuple(a-b for a,b in zip(u,w)) in S] for u in verts}
        best=[1,None]
        def lab(u,w):
            d=tuple(a-b for a,b in zip(u,w)); fnz=next(a for a in d if a)
            return d if fnz>0 else tuple(-a for a in d)
        def ext(cl,cand):
            if len(cl)>best[0]:
                if not distinct or len({lab(a,b) for a,b in combinations(cl,2)})==len(cl)*(len(cl)-1)//2:
                    best[0]=len(cl); best[1]=list(cl)
            for i,w in enumerate(cand):
                if distinct:
                    labs=[lab(a,b) for a,b in combinations(cl+[w],2)]
                    if len(set(labs))<len(labs): continue
                ext(cl+[w],[x for x in cand[i+1:] if x in adj[w]])
        for u in verts:
            ext([u],[w for w in adj[u] if w>u])
        return best
    print('max min-shell clique (all dirs):',maxclique(dirs))
    print('max min-shell clique distinct labels:',maxclique(dirs,True)[0])
    if selected is not None:
        sel=[dirs[k] for k in selected]
        print('max clique selected only:',maxclique(sel)[0])
        tri=[]
        for a,b,c in combinations(range(len(sel)),3):
            for s,t in product((1,-1),repeat=2):
                if all(x+s*y+t*z==0 for x,y,z in zip(sel[a],sel[b],sel[c])): tri.append((selected[a],selected[b],selected[c],s,t))
        print('signed triangles among selected:',tri)
        # max fibre coverage of selected labels
        supp=[sum(1<<j for j in range(8) if v[j]) for v in sel]
        cov=max(len({k for k,s in enumerate(supp) if any(a^b==s for a,b in combinations(f,2))}) for f in fib.values())
        print('max selected labels covered by one fibre:',cov)
    return dirs

if __name__=='__main__':
    analyse('cube',CUBE)
    d1=analyse('C1',C1,None)
    d2=analyse('C2',C2,None)
