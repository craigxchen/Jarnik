import random, time, math, sys
from collections import Counter
from gtight import *
def rand_primes(r, lo, hi, rng):
    out=set()
    while len(out)<r:
        n=rng.randrange(lo,hi)|1
        if n%4==1 and is_prime(n): out.add(n)
    return sorted(out)
def run(M,b,K,seed=0,Xf=2):
    rng=random.Random(seed)
    X=Xf*M; r=b*(M-1); P=int(M**K)
    ps=rand_primes(r,P,2*P,rng)
    gp=[]
    for p in ps:
        a,bb=two_squares(p)
        if rng.random()<0.5: bb=-bb
        gp.append((a,bb))
    rels=tight_relations(gp,X)
    kinds=Counter(rel[0] for rel in rels)
    Pn,log=closure(r,rels)
    st=Counter(pinned_status(Pn,r))
    comp=components(Pn,r)
    big=max(comp.values())
    npairs=r*(r-1)//2*8
    return dict(M=M,b=b,K=K,r=r,X=X,S=kinds['S'],T=kinds['T'],Pp=kinds['P'],
                dens_pair=kinds['P']/npairs, deg=2*kinds['P']/r, free=st[0],half=st[1],exact=st[2],bigcomp=big,certs=len(log))
if __name__=='__main__':
    for K in [float(x) for x in sys.argv[1].split(',')]:
        for M in [int(x) for x in sys.argv[2].split(',')]:
            for b in [int(x) for x in sys.argv[3].split(',')]:
                t=time.time()
                d=run(M,b,K)
                print(' '.join(f'{k}={v:.3g}' if isinstance(v,float) else f'{k}={v}' for k,v in d.items()), 'time=%.0f'%(time.time()-t), flush=True)
