# independent pure-python brute force of M_2(D,0): group all monic z by (norm, coefficients above D/2)
import sys
from f2lib import *
from collections import Counter
for D in range(1,int(sys.argv[1])+1):
    K=D//2; C=Counter(); ex={}
    for a in range(1<<D):
        A=(1<<D)|a
        for B in range(1<<D):
            key=(norm4((A,B)),A>>(K+1),B>>(K+1)); C[key]+=1; ex.setdefault(key,[]).append((A,B))
    M=max(C.values()); k=[k for k in C if C[k]==M][0]
    # recheck the cluster directly
    pts=ex[k]
    assert all(norm4(p)==k[0] for p in pts) and all(2*deg4(add4(p,q))<=D for p in pts for q in pts if p!=q)
    print(D,M,len([1 for v in C.values() if v==M]),[ (bin(p[0]),bin(p[1])) for p in pts][:4])
