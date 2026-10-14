from rlib import *
import itertools, pickle, time
A=build_A(2); print("|A|",len(A))
# Hilbert function of R/J(A): cumulative exact ranks, monomials with exponents <=2 (x^3-x in I(A))
cum=[]
for e in range(0,15):
    ms=[m for m in monos_le(7,e) if max(m)<=2]
    r=fmpz_mat([[evmono(m,a) for m in ms] for a in A]).rank()
    cum.append(r)
    if r==len(A): break
H=[cum[0]]+[cum[i]-cum[i-1] for i in range(1,len(cum))]
print("H(R/J(A))",H, "sum",sum(H))
dimJ=[math.comb(i+6,6)-H[i] for i in range(9)]
print("dim J(A)_i i<=8",dimJ)
# width set: primitive v, |v_i|<=3, tau+1<=8
W=[]
for v in itertools.product(range(-3,4),repeat=7):
    if not any(v): continue
    g=0
    for x in v: g=math.gcd(g,abs(x))
    if g!=1: continue
    # canonical sign: first nonzero positive
    fn=[x for x in v if x][0]
    if fn<0: continue
    vals=set(sum(a*b for a,b in zip(k,v)) for k in A)
    tau=len(vals)-1; h=max(vals)-min(vals)
    if tau+1<=8: W.append((v,tau,h))
from collections import Counter
print("width dirs",len(W),Counter(h for v,t,h in W), "tau==h all:",all(t==h for v,t,h in W), "max|v|",max(max(abs(x) for x in v) for v,t,h in W))
pickle.dump((A,H,dimJ,W),open('m8_a.pkl','wb'))
