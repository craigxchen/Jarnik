# Referee: membership of the count-feasible products in J(A)_m^perp for a small model (exact rank criterion).
import itertools, numpy as np, sys
from math import factorial
from flint import fmpz_mat, nmod_mat
pre=sys.argv[1]; m=int(sys.argv[2]); PR=2147483647
A=np.load(pre+'_A.npy'); d=A.shape[1]
pool=[list(map(int,l.split())) for l in open(pre+'_pool.txt').read().split('\n')[1:] if l.strip()]
mons=[a for a in itertools.product((0,1,2),repeat=d) if sum(a)<=m]
K=[[int(np.prod([k[i]**a[i] for i in range(d)])) for k in A] for a in mons]
print('|A|',len(A),'monomials',len(mons))
rq=fmpz_mat(K).rank(); rp=nmod_mat(K,PR).rank(); print('rank_Q',rq,'rank_p',rp)
assert rq==rp
idx={a:i for i,a in enumerate(mons)}
def polyprod(zs):
    P={(0,)*d:1}
    for z in zs:
        Q={}
        for e,c in P.items():
            for x in range(d):
                if z[x]:
                    f=list(e); f[x]+=1; f=tuple(f); Q[f]=Q.get(f,0)+c*z[x]
        P=Q
    return P
def multinom(a):
    r=factorial(sum(a))
    for x in a: r//=factorial(x)
    return r
Kp=nmod_mat(K,PR)
nm=0; nn=0; nonred=0
for ln in open(pre+'_leaves.txt'):
    zs=[pool[i][1:] for i in map(int,ln.split())]
    P=polyprod(zs)
    if any(c and max(e)>2 for e,c in P.items()): nonred+=1; nn+=1; continue
    b=[0]*len(mons)
    for e,c in P.items():
        if c: b[idx[e]]=c*pow(multinom(e),-1,PR)%PR
    aug=nmod_mat([K[i]+[b[i]] for i in range(len(mons))],PR)
    if aug.rank()==rp: nm+=1
    else: nn+=1
print('leaves',nm+nn,'members mod p',nm,'non-members',nn,'(of which nonreduced coefficient)',nonred,' negatives rigorous since rank_Q==rank_p')
