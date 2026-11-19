# Transfer of the record's pigeonhole + chain-lemma lower bound (two_thirds/proof.md) to integral fakes
# (no norm condition) over A = Z and A = F_q[t].  Upper side (exact, see report):
#    sum_{i<j} log|c_ij| + Q/2 <= m W/8 + K log C1 ,   c_ij = D_ij/g_ij.
# Lower side: prime l not dividing N: columns in P^1(A/l^a), (|l|+1)|l|^(a-1) classes;
# prime p | N: chain lemma with lambda = 1/(|p|+1)  (a valid, slightly lossy, version).
import math
from sympy import primerange, mobius, divisors
def E(n,Q):
    if Q>=n: return 0
    b,r=divmod(n,Q); return Q*b*(b-1)//2+r*b
def cchain(lam): return (math.sqrt(1+4*lam)-1)/2
def free(m,P):
    s=0;Q=P+1
    while Q<m: s+=E(m,Q); Q*=P
    return s
def Gp(m,P,can_divide=True):
    f=free(m,P)
    if not can_divide: return f
    ch=max(0.0,cchain(1/(P+1))*m*m/4-m/2)
    return min(f,ch)
def LB_Z(m,splitonly):
    tot=m*(m-1)/2*0  # no ramified gain claimed for no-norm model
    for l in primerange(2,m+1):
        cand = (not splitonly) or (l%4==1)
        tot+=math.log(l)*Gp(m,l,cand)
    return tot
def nirr(q,k): return sum(mobius(d)*q**(k//d) for d in divisors(k))//k
def LB_Fq(m,q):
    tot=0;k=1
    while q**k+1<m:
        tot+=nirr(q,k)*k*math.log(q)*Gp(m,q**k,True); k+=1
    return tot
print(" m   | Z(any N)  Z(split N)  F2[t]   F3[t]   F_{m+1}[t]   (entries: 8*LB/(m * m log m), i.e. coefficient c in W>=c m log m - ...)")
for m in [16,64,256,1024,4096]:
    row=[8*LB_Z(m,False)/m, 8*LB_Z(m,True)/m, 8*LB_Fq(m,2)/m, 8*LB_Fq(m,3)/m, 8*LB_Fq(m,m+1)/m]
    print(f"{m:5d} | "+"  ".join(f"{x/(m*math.log(m)):7.3f}" for x in row))
# trace construction check of the exact transferred inequality (degree units, q=2^d, m=q-1, C1=1):
for d in [5,6,7,10,20]:
    q=2**d;m=q-1;L=q-1;s=q//2-1
    lhs=0+L*(s-m/2)**2/2; rhs=m*L/8
    print(f"trace q={q}: sum deg c_ij + Q/2 = {lhs}  <=  mL/8 = {rhs}: {lhs<=rhs};  pigeonhole gain = {LB_Fq(m,q)}")
