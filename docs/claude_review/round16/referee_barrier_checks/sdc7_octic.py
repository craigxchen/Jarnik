# sdc7, L=2, m=8: exact V(A,8) (integer nullspace), translation invariance of c, and value of Phi_c on Fourier circles (mpmath 50 digits).
import sys, itertools, numpy as np, flint, mpmath as mp
sys.path.insert(0,'/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round16/referee_barrier_checks')
from thm32_check import box, monos
mp.mp.dps=50
q,L,m=7,2,8
A=box(q,L); d=q
ms=monos(d,m-1)
E=[[int(np.prod([k[x] for x in c])) if c else 1 for c in ms] for k in A]
ET=flint.fmpz_mat([[E[i][j] for i in range(len(A))] for j in range(len(ms))])
X,nul=ET.nullspace()
print("|A|=",len(A)," dim V(A,8)=",nul,flush=True)
c=[int(X[i,0]) for i in range(len(A))]
from math import gcd
g=0
for v in c: g=gcd(g,abs(v))
c=[v//g for v in c]
idx={tuple(k):i for i,k in enumerate(A)}
inv=all(c[idx[tuple(np.roll(k,1))]]==c[i] for i,k in enumerate(A))
neg=all(c[idx[tuple(-np.array(k))]]==c[i] for i,k in enumerate(A))
print("c translation-invariant:",inv," c even:",neg," |c|_1=",sum(abs(v) for v in c),flush=True)
for a in [1,2,3]:
    vals=[]
    for th in [mp.mpf(0),mp.mpf('0.3'),mp.mpf('1.1'),mp.pi/7]:
        u=[mp.cos(th+2*mp.pi*a*x/q) for x in range(q)]
        s=mp.mpf(0)
        for i,k in enumerate(A):
            if c[i]: s+=c[i]*mp.fsum(int(k[x])*u[x] for x in range(q))**m
        vals.append(s)
    print(f"a={a}: Phi_c on C_a at 4 thetas:",[mp.nstr(v,20) for v in vals],flush=True)
