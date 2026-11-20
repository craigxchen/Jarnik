# IX.2 check: 7 finite marks, blocks e_T (3<=|T|<=6) pairwise coprime, cofactors c_ij=(y_i-y_j)/prod_{T>={i,j}} e_T.
# At p=2 at most one e_T is even (valuation v). Minimal number of even c_ij over all residues mod 2^(v+1)
# (the 2-adic valuation of c_ij is >=1 iff v2(y_i-y_j) > v*[i,j in T]); enumerate y mod 2^(v+2) to be safe.
from itertools import combinations, product
pairs=list(combinations(range(7),2))
def mineven(T,v,mod):
    best=99
    Ts=set(T)
    for y in product(range(mod),repeat=7):
        if T and any((y[i]-y[min(T)])%(2**v) for i in T): continue
        cnt=0
        for i,j in pairs:
            d=(y[i]-y[j])%mod
            need = v+1 if (T and i in Ts and j in Ts) else 1
            if d % (2**need)==0: cnt+=1
        best=min(best,cnt)
        if best<2: return best
    return best
res={}
res['none even']=mineven((),0,2)
for k in range(3,7):
    for v in (1,2):
        res[(k,v)]=mineven(tuple(range(k)),v,2**(v+1))
print(res)
