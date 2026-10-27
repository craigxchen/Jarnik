# Referee: EXACT check (python-flint, rational arithmetic) of Prop 4.1/4.2 on the toy sdc7, L=2, m=6 (< q=7):
# N_i = |top(khat(1)^i khat(-1)^(m-i))|^2 in L^2(mu_A) is rational, N_i = N_(m-i), and
# sum_i binom(m,i)^2 N_i / 4^m equals Psi_A(u_theta) (checked at theta=0 against an exact rational evaluation
# is not possible since u_theta is irrational; we check the identity numerically at 30 digits instead).
import itertools, flint, mpmath as mp
from math import comb
mp.mp.dps=30
q=7; L=2; m=6
def chi(a):
    a%=q
    return 0 if a==0 else (1 if pow(a,(q-1)//2,q)==1 else -1)
pats=[[a]+[x for x in range(q) if chi(a-x)==1] for a in range(q)]+[[a]+[x for x in range(q) if chi(a-x)==-1] for a in range(q)]
A=[k for k in itertools.product([-1,0,1],repeat=q) if all(abs(sum(k[x] for x in S))<=L for S in pats)]
U=len(A); print("|A|=",U)
ms=[e for e in itertools.product(range(3),repeat=q) if sum(e)<m]
def mon(k,e):
    v=1
    for x in range(q): v*=k[x]**e[x]
    return v
Bt=flint.fmpz_mat([[mon(k,e) for k in A] for e in ms])   # (#monos) x U
C=Bt.nullspace()[0]   # U x U, first 'nullity' columns span ker(Bt) = V
rk=Bt.rank(); nul=U-rk
print("rank Poly_<m|_A =",rk," dim V =",nul)
Cm=flint.fmpz_mat([[C[r,c] for c in range(nul)] for r in range(U)])
G=Cm.transpose()*Cm
Ginv=flint.fmpq_mat(G).inv()
# Z[zeta] arithmetic: polys mod Phi_7, coefficient vectors in basis 1..zeta^5
Phi=flint.fmpz_poly([1]*q)
def X(e): return flint.fmpz_poly([0]*(e%q)+[1])
kp=[]; km=[]
for k in A:
    kp.append(sum((X(x)*k[x] for x in range(q)),flint.fmpz_poly([0]))%Phi)
    km.append(sum((X(-x)*k[x] for x in range(q)),flint.fmpz_poly([0]))%Phi)
def coeffs(f):
    return [int(f[t]) for t in range(q-1)]
Nvals=[]
for i in range(m+1):
    vecs=[coeffs((kp[r]**i*km[r]**(m-i))%Phi) for r in range(U)]   # U x 6
    Mt=flint.fmpz_mat([[vecs[r][t] for t in range(q-1)] for r in range(U)])
    Y=Cm.transpose()*Mt            # nul x 6
    Gst=flint.fmpq_mat(Y).transpose()*Ginv*flint.fmpq_mat(Y)   # 6x6 rational: <m^(s), P_V m^(t)>
    # N_i = (1/U) sum_{s,t} G_st conj(zeta^s) zeta^t = (1/U) sum G_st zeta^(t-s); reduce mod Phi in Q(zeta)
    acc=[flint.fmpq(0)]*q
    for s_ in range(q-1):
        for t_ in range(q-1):
            acc[(t_-s_)%q]+=Gst[s_,t_]
    # zeta^0..zeta^6 with zeta^6 = -(1+...+zeta^5): subtract acc[6] from all
    red=[acc[e]-acc[q-1] for e in range(q-1)]
    rational = all(red[e]==0 for e in range(1,q-1))
    N=red[0]/U
    Nvals.append(N)
    print(f"i={i}: N_i rational: {rational}; N_i = {N}  (~{float(N.p)/float(N.q):.10g})",flush=True)
print("N_i = N_(m-i):",all(Nvals[i]==Nvals[m-i] for i in range(m+1)))
pred=sum(comb(m,i)**2*Nvals[i] for i in range(m+1))/4**m
print("predicted Psi_A on the circle (exact rational) =",pred," ~",float(pred.p)/float(pred.q))
# 30-digit numerical Psi_A(u_theta) via the same V basis: Psi = (1/U) y^T Ginv y, y = C^T g
GI=mp.matrix([[mp.mpf(int(Ginv[a,b].p))/mp.mpf(int(Ginv[a,b].q)) for b in range(nul)] for a in range(nul)])
for th in [0,mp.mpf('0.3'),mp.mpf('1.1')]:
    u=[mp.cos(th+2*mp.pi*x/q) for x in range(q)]
    g=[sum(k[x]*u[x] for x in range(q))**m for k in A]
    y=[mp.fsum(int(Cm[r,c])*g[r] for r in range(U)) for c in range(nul)]
    Gi=GI if "GI" in globals() else mp.matrix([[mp.mpf(int(Ginv[a,b].p))/mp.mpf(int(Ginv[a,b].q)) for b in range(nul)] for a in range(nul)])
    yv=mp.matrix(y)
    val=(yv.T*Gi*yv)[0]/U
    print(f"theta={th}: Psi_A(u_theta) = {mp.nstr(val,25)}  rel.diff to prediction = {mp.nstr(val/(mp.mpf(int(pred.p))/mp.mpf(int(pred.q)))-1,5)}",flush=True)
