# Count solutions mod k of the formal k-th-power conditions in the uniqueness calculation.
from itertools import combinations, product
from math import gcd
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
marks=range(5); pairs=list(combinations(marks,2))
# variables: a_0..a_4, v_e (10)
rows=[]
for (i,j) in pairs:
    r=[0]*15; r[i]+=1; r[j]+=1; r[5+pairs.index((i,j))]+=1; rows.append(r)       # beta exponent
    r=[0]*15; r[i]-=1; r[j]-=1; r[5+pairs.index((i,j))]+=1
    for e in pairs:
        if i in e: r[5+pairs.index(e)]-=1
        if j in e: r[5+pairs.index(e)]-=1
    rows.append(r)                                                               # conj beta exponent
rows.append([1]*15)  # homogeneity (exact =0, a fortiori mod k)
M=sp.Matrix(rows)
S=smith_normal_form(M,domain=sp.ZZ)
d=[abs(S[i,i]) for i in range(min(S.shape)) if S[i,i]!=0]
print("Smith invariants:",d,"rank",len(d))
Jvec=[-2]*5+[1]*10
print("J exponents mapped:",list(M*sp.Matrix(Jvec))); assert all(x%3==0 for x in M*sp.Matrix(Jvec))
for k in range(2,25):
    cnt=1
    for x in d: cnt*=gcd(x,k)
    cnt*=k**(15-len(d))
    print(k,cnt,"expected",gcd(3,k))
    assert cnt==gcd(3,k)
