# Check: on a genuine parent, the e-all collision modulus of v=e_j-e_k at unit i^s, computed from residues
# r_j = pi_j mod q^a, equals the X-truncated odd exterior part of D_s(pi_j conj pi_k) in {Y, X-Y, X, X+Y}.
import math
from sympy import primerange
def gauss(p):
    for b in range(1,math.isqrt(p)+1):
        a=math.isqrt(p-b*b)
        if a*a==p-b*b and a>=b: return (a,b)
P=100; Xc=3000
S=[p for p in primerange(P,2*P+1) if p%4==1]; pis=[gauss(p) for p in S]
def mul(z,w,m): return ((z[0]*w[0]-z[1]*w[1])%m,(z[0]*w[1]+z[1]*w[0])%m)
def inv(z,m):
    n=(z[0]*z[0]+z[1]*z[1])%m; ni=pow(n,-1,m); return ((z[0]*ni)%m,(-z[1]*ni)%m)
units=[(1,0),(0,1),(-1,0),(0,-1)]
qs=[q for q in primerange(3,Xc+1) if q not in S]
bad=0; tot=0; maxlog=0
for j in range(len(S)):
  for k in range(len(S)):
    if j==k: continue
    a,b=pis[j]; c,d=pis[k]
    A=(a*c+b*d, b*c-a*d)  # pi_j * conj(pi_k)
    Xr,Yr=A
    Ds=[Yr, Xr-Yr, Xr, Xr+Yr]
    for s in range(4):
      m_res=1; m_int=1
      for q in qs:
        A_q=int(math.log(Xc)/math.log(q)+1e-12)
        while q**(A_q+1)<=Xc: A_q+=1
        while q**A_q>Xc: A_q-=1
        # residue side
        aa=0
        for e in range(1,A_q+1):
          M=q**e
          rj=( a%M, b%M); rk=(c%M,d%M)
          tj=mul(rj,inv((a%M,(-b)%M),M),M); tk=mul(rk,inv((c%M,(-d)%M),M),M)
          t=mul(tj,inv(tk,M),M)
          u=(units[s][0]%M,units[s][1]%M)
          if t==u: aa=e
          else: break
        m_res*=q**aa
        # integer side
        v=0; D=Ds[s]
        while D%q==0 and v<A_q: D//=q; v+=1
        m_int*=q**v
      tot+=1; bad+= (m_res!=m_int)
      # sine inequality with genuine angle: |D_s| >= m
      assert abs(Ds[s])>=m_int and Ds[s]!=0
      maxlog=max(maxlog, math.log(m_int)/math.log(math.sqrt(S[j]*S[k])))
print('pairs*units',tot,'mismatches',bad,'max log m / log sqrt(pj pk)',round(maxlog,3))
