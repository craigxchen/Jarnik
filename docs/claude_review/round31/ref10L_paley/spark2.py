from itertools import combinations
from paley1 import chi, S_of
from spark import rank_mod
for m in (11,19):
    S=S_of(m); half=[i for i in range(m) if chi(0-i,m)==1]
    M=[i for i in range(m) if i not in half]; k=len(M)
    T=[[1]*m]+[S[i] for i in M]
    for u in range(1,k+2):
        defic=[U for U in combinations(range(m),u) if rank_mod([[T[r][c] for c in U] for r in range(k+1)])<u]
        if defic:
            print("m=%d k=%d spark=%d (vs generic k+2=%d); #deficient sets=%d; contain col 0: %d; example %s"%(m,k,u,k+2,len(defic),sum(0 in U for U in defic),defic[0])); break
