# Evidence: search structured multiplicatively dependent families for large min-anchor/(n h(ell)).
import sys; sys.path.insert(0,'.')
from gauss import *
from itertools import product, combinations
def prodg(gs,ex):
    r=(1,0)
    for g,e in zip(gs,ex):
        if e>=0: r=mul(r,pw(g,e))
        else: r=mul(r,pw(conj(g),-e))
    return r
def stats(u):
    n=len(u)-1; hl,G=h_line(u); s,a=minanchor(u); return n,hl,G,s
def show(name,u,rank):
    if not nonreal_pairs(u): print(name,"has R-proportional pair; skip"); return
    n,hl,G,s=stats(u)
    print(f"{name}: n={n} rank={rank} h(ell)={hl:.3f} logG={log(G):.3f} minanchor={s:.3f} ratio/n={s/hl/n:.4f} (need >1/3 to beat S_min)",flush=True)
A=10**30
for r in [2,3,4,5]:
    gs=[(A+j,1) for j in range(1,r+1)]
    u=[prodg(gs,ex) for ex in product([0,1],repeat=r)]
    show(f"subsets r={r} A_j=A+j",u,r)
for r in [2,3,4]:
    gs=[(A**(j),1) for j in range(1,r+1)]
    u=[prodg(gs,ex) for ex in product([0,1],repeat=r)]
    show(f"subsets r={r} A_j=A^j",u,r)
for r,m in [(2,2),(2,3),(3,2)]:
    gs=[(A+j,1) for j in range(1,r+1)]
    u=[prodg(gs,ex) for ex in product(range(m+1),repeat=r)]
    show(f"box r={r} exps 0..{m}",u,r)
for r in [1,2,3]:
    b=[10**12+2*j for j in range(r)]
    gs=[(bj,1) for bj in b]
    u=[]
    for t in product([-1,0,1],repeat=r):
        z=(1,0)
        for g,tt in zip(gs,t): z=mul(z,mul(pw(g,1+tt),pw(conj(g),1-tt)))
        u.append(z)
    show(f"equal-norm circle products r={r}",u,r)
for n in [2,3,4,5,6,8]:
    g=(A,1); show(f"powers 0..{n}",[pw(g,k) for k in range(n+1)],1)
for n in [2,3,4,6]:
    g=(A,1); show(f"circle GP beta^x conj(beta)^(n-x), n={n}",[mul(pw(g,x),pw(conj(g),n-x)) for x in range(n+1)],1)
