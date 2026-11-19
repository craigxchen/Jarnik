import sys, itertools, collections
from witness import *
def dec(c):
    p=[]
    while c: c,d=divmod(c,9); p.append((d%3,d//3))
    return p
def norm(p): return pmul(p,pconj(p))
D=int(sys.argv[2]); KL=D//2; n_lines=int(sys.argv[3])
for line in open(sys.argv[1]).readlines()[:n_lines]:
    _,n,*bl=line.split(); U,V,T,G=[dec(int(x)) for x in bl]
    prs=[]
    for b in (U,V,T,G):
        if deg(b)>0: prs+=[pgcd(b,f) for f,e in factor_F3(norm(b))]
    groups=collections.defaultdict(list)
    for s in itertools.product((0,1),repeat=len(prs)):
        z=[(1,0)]
        for bit,p in zip(s,prs): z=pmul(z,p if bit else pconj(p))
        groups[tuple(z[KL+1:])].append(z)
    z1=pmul(pmul(G,U),pmul(pconj(V),pconj(T)))
    cl=groups[tuple(z1[KL+1:])]
    N=norm(cl[0]); assert all(norm(z)==N for z in cl) and all(deg(z)==D for z in cl)
    assert deg(pgcd(N,deriv(N)))==0
    gaps=[deg(psub(a,b)) for a,b in itertools.combinations(cl,2)]
    print(f"D={D} claimed {n} found cluster size {len(cl)} maxgap {max(gaps)} (<= {KL}) squarefree N, primes {len(prs)}; largest class on parent {max(len(g) for g in groups.values())}")
