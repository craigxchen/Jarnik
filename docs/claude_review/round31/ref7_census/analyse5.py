import sys, itertools, collections
sys.setrecursionlimit(10000)
from witness import *
D=int(sys.argv[2]); KL=D//2
def dec(c):
    p=[]
    while c: c,d=divmod(c,9); p.append((d%3,d//3))
    return p
def key(p): return tuple(p)
def norm(p): return pmul(p,pconj(p))
def primes_of(U):
    if deg(U)<=0: return []
    out=[]
    for f,e in factor_F3(norm(U)):
        out.append(pgcd(U,f))
    return out
parents={}
for line in open(sys.argv[1]):
    _,n,u,v,t,g=line.split(); blocks=[dec(int(x)) for x in (u,v,t,g)]
    N=[(1,0)]
    for b in blocks: N=pmul(N,norm(b))
    k=key(N)
    if k in parents: continue
    parents[k]=blocks
print(f"D={D}: distinct squarefree parents carrying a 5-cluster (from canonical triples): {len(parents)}")
profiles=collections.Counter(); sizes=collections.Counter()
for k,(U,V,T,G) in parents.items():
    prs=[]
    for b in (U,V,T,G): prs+=primes_of(b)
    z1=pmul(pmul(G,U),pmul(pconj(V),pconj(T)))
    groups=collections.defaultdict(list)
    for s in itertools.product((0,1),repeat=len(prs)):
        z=[(1,0)]
        for bit,p in zip(s,prs): z=pmul(z,p if bit else pconj(p))
        groups[tuple(z[KL+1:])].append(z)
    big=[gz for gz in groups.values() if len(gz)>=5]
    for gz in big:
        sizes[len(gz)]+=1
        m=len(gz); cuts=[]
        for p in prs:
            col=tuple(1 if not pmod(z,p) else 0 for z in gz)
            ones=[i for i in range(m) if col[i]]; zeros=[i for i in range(m) if not col[i]]
            mino=ones if len(ones)<=len(zeros) else zeros
            cuts.append((tuple(mino),deg(p)))
        # canonical form under S_m
        best=None
        for perm in itertools.permutations(range(m)):
            c=sorted((tuple(sorted(perm[i] for i in mi)),w) for mi,w in cuts)
            # merge equal minorities
            agg=collections.Counter()
            for mi,w in c: agg[mi]+=w
            c=tuple(sorted(agg.items()))
            if best is None or c<best: best=c
        profiles[best]+=1
        Nlist=[c[0] for c in k]
        even=all(Nlist[j]==0 for j in range(1,len(Nlist),2))
        profiles[('even' if even else 'noneven',best)]+=1
print("cluster sizes:",dict(sizes))
for p,c in profiles.items():
    if p[0] in ('even','noneven'): continue
    print(" profile (minority set: total upstairs degree):",p," count",c, " even-parent count", profiles[('even',p)])
