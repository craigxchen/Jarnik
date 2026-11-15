# are all realisations of the half-cube type at D=8 on the cube locus (even N after centering)?
import numpy as np, itertools, sys
sys.path.insert(0,'.'); from basesolve import *
D=8; L=3
U=np.array(list(itertools.product([0,1],repeat=3)))
Fs=np.array([np.ones(D,int)]+[1-U[:,i] for i in range(3)])
out=solve_base(Fs,L,starts=400,seed=7)
good=[o for o in out if o[0]<1e-10 and o[1]>1e-3]
print("nondegenerate solutions:",len(good))
asym=[]
for cost,dg,q in good:
    allr=np.concatenate([q,np.conj(q)]); allr=allr-allr.real.mean()
    # closed under negation?
    d=max(min(abs(-z-w) for w in allr) for z in allr)
    asym.append(d)
asym=np.array(asym); print("max over roots of dist(-root, rootset) after centering: min %.2e max %.2e"%(asym.min(),asym.max()))
print("count symmetric (<1e-8):",(asym<1e-8).sum())
