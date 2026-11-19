import pickle, numpy as np, sympy as sp
src=open('census_ref.py').read().split('# enumerate all cliques')[0]
exec(src)
def arank(S):
    return sp.Matrix([[x-y for x,y in zip(bits(r),bits(S[0]))] for r in S[1:]]).rank() if len(S)>1 else 0
reps=[R for R in pickle.load(open('reps6.pkl','rb')) if arank(R)<=4]
print('6-row orbits of rank<=4:',len(reps))
k=6
while reps:
    new={}
    for R in reps:
        cand=set.intersection(*[adj[r] for r in R])
        for s in cand:
            TT=R+(s,)
            if arank(TT)<=4:
                c=canon(TT); new[c]=1
    reps=sorted(new); k+=1
    print(k,len(reps),[arank(R) for R in reps])
    if k==8: pickle.dump(reps,open('reps8_rank4.pkl','wb'))
