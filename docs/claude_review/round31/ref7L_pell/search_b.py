# exhaustive exact search over multiplicative Pell-block templates at n=1 (and n=2 for r<=3)
import sys; sys.path.insert(0,'.')
from pell_core import *
import itertools, random
from sympy import Matrix
A=[[P(1,0,1),P(0,1)],[P(0,1),P(1)]]; beta=P(G(0,2),1); betab=pconj(beta)
def powp(p,k):
    r=P(1)
    for _ in range(k): r=pmul(r,p)
    return r
def run(r,w,w0,n,seeds,maxsub=None):
    An=matpow(A,n); Gs=[block(An,s) for s in seeds]; Gb=[pconj(x) for x in Gs]
    vecs=list(itertools.product(*[range(x+1) for x in w],range(w0+1)))
    classes={}
    for v in vecs:
        ch=sum(2*v[l]-w[l] for l in range(r))+2*(2*v[r]-w0)
        classes.setdefault(ch,[]).append(v)
    cache={}
    def mark(v):
        if v not in cache:
            z=pmul(powp(beta,v[r]),powp(betab,w0-v[r]))
            for l in range(r): z=pmul(z,pmul(powp(Gs[l],v[l]),powp(Gb[l],w[l]-v[l])))
            cache[v]=z
        return cache[v]
    best={True:None,False:None}; cnt=0
    for ch,cl in classes.items():
        for sub in itertools.combinations(cl,5):
            ms=[mark(v) for v in sub]
            a=analyse(ms); cnt+=1
            if not a['distinct']: continue
            rk=Matrix([[sub[i][l]-sub[0][l] for l in range(r)] for i in range(1,5)]).rank()
            dense=(rk==4); val=2*a['h']-a['D']
            if best[dense] is None or val>best[dense][0]: best[dense]=(val,a['D'],a['h'],sub)
    return cnt,best
random.seed(1)
for (r,w,w0,n) in [(4,(1,1,1,1),1,1),(4,(1,1,1,1),2,1),(3,(1,1,1),2,1),(3,(2,1,1),1,1),(3,(2,2,1),1,1),(3,(1,1,1),1,2),(5,(1,1,1,1,1),1,1),(4,(2,1,1,1),1,1)]:
    seeds=[Fr(random.randint(-40,40),random.randint(1,9)) for _ in range(r)]
    cnt,best=run(r,w,w0,n,seeds)
    print("r",r,"w",w,"w0",w0,"n",n,"subsets",cnt,"best dense (2h-D,D,h,marks):",best[True],"\n    best non-dense:",best[False],flush=True)
