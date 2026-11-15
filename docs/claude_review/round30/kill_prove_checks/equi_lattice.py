# Independent: equidistant (pairwise orthogonal) M-point codes in {+-1}^8 (c=0, L=3), and the saturated-lattice
# necessary conditions of Lemma 1.2: every nonzero a in span_R(differences) cap Z^8 has ||a||_1>=L+1 and
# |sum_k sigma_k a_k| <= ||a||_1 - L (sigma = normalisation signs), plus Gordan (exists y>0 orth to sigma*V).
import itertools, sys, numpy as np, sympy as sp
from scipy.optimize import linprog
D=8; L=3; M=int(sys.argv[1])
ones=np.ones(D,int)
W=[np.array(v) for v in itertools.product([1,-1],repeat=D) if sum(v)==0]   # orthogonal to all-ones
n=len(W)
adj=[[int(W[i]@W[j]==0) for j in range(n)] for i in range(n)]
cl=[]
def ext(c,cand):
    if len(c)==M-1: cl.append(c); return
    for j in cand:
        if j>(c[-1] if c else -1):
            ext(c+[j],[k for k in cand if k>j and adj[j][k]])
ext([],list(range(n)))
print('cliques (s1=1):',len(cl))
def canon(rows):
    best=None
    for perm in itertools.permutations(range(len(rows))):
        R=np.array([rows[p] for p in perm])
        key=tuple(sorted(tuple(R[:,k]) for k in range(D)))
        if best is None or key<best: best=key
    return best
classes={}
for c in cl:
    key=canon([W[j] for j in c])
    classes.setdefault(key,c)
print('classes up to S_8 x S_{M-1} (s1=1, column flips used):',len(classes))
# short integer vectors, l1<=7
short=[]
for j in range(1,8):
    for supp in itertools.combinations(range(D),j):
        # compositions of total <=7 into j positive parts
        for tot in range(j,8):
            for parts in itertools.product(range(1,tot+1),repeat=j):
                if sum(parts)!=tot: continue
                for sg in itertools.product([1,-1],repeat=j):
                    a=np.zeros(D,int)
                    for k,p,s in zip(supp,parts,sg): a[k]=p*s
                    short.append(a)
short=np.array(short); print('short vectors',len(short))
l1=np.abs(short).sum(1)
results=[]
for key,c in classes.items():
    pts=[ones]+[W[j] for j in c]
    diffs=sp.Matrix([list((p-ones)//2) for p in pts[1:]])  # (M-1) x D
    perp=diffs.nullspace()   # vectors orthogonal to rows
    A=np.array([[int(x) for x in (v*sp.ilcm(*[sp.fraction(e)[1] for e in v]))] for v in perp]) if perp else np.zeros((0,D),int)
    inV=np.all(short@A.T==0,axis=1) if len(A) else np.ones(len(short),bool)
    sv=short[inV]; svl1=np.abs(sv).sum(1)
    minl1=svl1.min() if len(sv) else None
    surv=[]
    for sig in itertools.product([1,-1],repeat=D):
        sig=np.array(sig)
        if len(sv) and np.any(np.abs(sv@sig) > svl1-L): continue
        if len(sv) and minl1<L+1: continue
        # Gordan: y>0 (y>=1 after scaling) with y . (sig*diff)=0 for all diffs
        Aeq=np.array(diffs.tolist(),float)*sig
        lp=linprog(np.zeros(D),A_eq=Aeq,b_eq=np.zeros(M-1),bounds=[(1,None)]*D,method='highs')
        if lp.status!=0: continue
        surv.append(tuple(sig))
    results.append((key,minl1,len(surv)))
    print('class',len(results)-1,'min l1 in saturation:',minl1,'surviving sigma:',len(surv),
          'example:',surv[0] if surv else None, flush=True)
