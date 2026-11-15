# enumerate cluster types: M sign vectors in {+-1}^D (normalised), pairwise pos,neg>=h, pos+neg>=L+1,
# modulo S_M (point relabeling), S_D (column order: multisets), global negation.
import itertools, numpy as np, sys
from scipy.optimize import linprog
def enum_types(D,M,L,h,allow_const=True):
    pats=[p for p in itertools.product([1,-1],repeat=M) if allow_const or len(set(p))>1]
    pats=[p for p in pats if p[0]==1 or True]
    P=np.array(pats); npat=len(pats)
    pairs=list(itertools.combinations(range(M),2))
    pos=np.array([[1 if (p[i]==1 and p[j]==-1) else 0 for (i,j) in pairs] for p in pats])
    neg=np.array([[1 if (p[i]==-1 and p[j]==1) else 0 for (i,j) in pairs] for p in pats])
    out=set(); res=[]
    perms=list(itertools.permutations(range(M)))
    def canon(cols):
        best=None
        for pm in perms:
            for s in (1,-1):
                t=tuple(sorted(tuple(s*pats[c][q] for q in pm) for c in cols))
                if best is None or t<best: best=t
        return best
    def dfs(start, cols, ps, ns):
        left=D-len(cols)
        need=np.maximum(h-ps,0)+np.maximum(h-ns,0)
        need=np.maximum(need, (L+1)-(ps+ns))
        if need.max()>left: return
        if left==0:
            cf=canon(cols)
            if cf not in out: out.add(cf); res.append(cf)
            return
        for c in range(start,npat):
            dfs(c, cols+[c], ps+pos[c], ns+neg[c])
    dfs(0,[],np.zeros(len(pairs),int),np.zeros(len(pairs),int))
    return res
def gordan_ok(E):
    # exists y>0 with (E_i-E_0).y = 0 for all i  (j=1 moment equation, Im rho>0)
    E=np.array(E,float).T if False else np.array(E,float)
    A=(E[1:]-E[0])/2; D=E.shape[1]
    r=linprog(np.zeros(D),A_eq=A,b_eq=np.zeros(len(A)),bounds=[(1,None)]*D,method='highs')
    return r.status==0
def lattice_ok(E,h,L,cmax=2):
    E=np.array(E,int); d=(E[1:]-E[0])//2; n=len(d)
    for c in itertools.product(range(-cmax,cmax+1),repeat=n):
        if not any(c): continue
        a=np.array(c)@d
        if not a.any(): continue
        if (a[a>0].sum()<h) or (-a[a<0].sum()<h) or abs(a).sum()<L+1: return False
    return True
if __name__=="__main__":
    D,M,c=map(int,sys.argv[1:4]); L=D-1-(D+2*c)//2; h=(L+1)//2
    T=enum_types(D,M,L,h)
    print("D",D,"M",M,"L",L,"h",h,"types",len(T))
    T2=[t for t in T if gordan_ok(np.array(t).T)]
    print("after Gordan",len(T2))
    T3=[t for t in T2 if lattice_ok(np.array(t).T,h,L)]
    print("after lattice(c<=2)",len(T3))
    import pickle; pickle.dump(T3,open(f"types_D{D}_M{M}_c{c}.pkl","wb"))
