# complex solutions of the sigma system (c=-1), then rebuild V_i and verify the (G1) identity
import numpy as np, mpmath as mp, sys
from scipy.optimize import root
E=[(0,1),(0,2),(0,3),(0,4),(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)]
inc=[[k for k,e in enumerate(E) if i in e] for i in range(5)]
c=-1.0
rng=np.random.default_rng(int(sys.argv[1]))
def Fc(z,K):
    s=z[:10]+1j*z[10:]; r=[]
    for i in range(5):
        r.append(np.prod(s[inc[i]])/K-1); r.append(np.sum(s[inc[i]]+c/s[inc[i]]))
    r=np.array(r); return np.concatenate([r.real,r.imag])
sol=None
for t in range(500):
    K=complex(rng.normal(),rng.normal())*2
    S=root(Fc,rng.normal(size=20),args=(K,),method='lm')
    if np.max(np.abs(Fc(S.x,K)))<1e-12:
        s=S.x[:10]+1j*S.x[10:]; rho=(c-s**2)/s
        g=min(abs(rho[a]-rho[b]) for a in range(10) for b in range(a))
        if g>1e-3: sol=(s,K); break
print('complex nondeg solution found:',sol is not None)
s,K=sol
# refine with mpmath
mp.mp.dps=50
def eqs(*v):
    ss=v; out=[]
    for i in range(5):
        p=mp.mpf(1)
        for k in inc[i]: p*=ss[k]
        out.append(p-K); 
    for i in range(5):
        out.append(sum(ss[k]+c/ss[k] for k in inc[i]))
    return out
ss=mp.findroot(eqs,[mp.mpc(x) for x in s])
ss=[ss[k] for k in range(10)]
cK=mp.mpc(K)
a2=-c**3/cK; q=c/a2; m=(1+q)/(1-q); L2=c/(m**2-1)
print('m=',mp.nstr(m,8),' L^2=',mp.nstr(L2,8))
L=mp.sqrt(L2)
# xi_i(s)= s^-3 (s^2+a^2) prod(s-sigma); P(s)=(c-s^2)/(2 m s), w=(c+s^2)/(2s)
import itertools
Ps=[mp.mpf(k)/7+mp.mpf(1)/3 for k in range(1,9)]
def s_of_P(P): # root of s^2+2mPs-c=0
    return (-2*m*P+mp.sqrt(4*m*m*P*P+4*c))/2
def xi(i,sv):
    v=sv**-3*(sv**2+a2)
    for k in inc[i]: v*=(sv-ss[k])
    return v
def u_y(i,P):
    sv=s_of_P(P); si=-c/sv
    w=(c+sv**2)/(2*sv)
    return (xi(i,sv)+xi(i,si))/2,(xi(i,sv)-xi(i,si))/(2*w)
# y must be constant
for i in range(5):
    ys=[u_y(i,P)[1] for P in Ps[:3]]
    print('vertex',i,'y spread',mp.nstr(max(abs(ys[0]-y) for y in ys),5))
def V(i,P):
    u,y=u_y(i,P); A=P*P+L2; x=(u-P*y)/A
    return (x*P+y, L*x), y
mu=m*m
maxres=0
for P in Ps:
    Vs=[V(i,P) for i in range(5)]
    for i in range(5):
        (X,Y),y=Vs[i]
        n=X*X+Y*Y-mu*y*y
        prod=1
        for j in range(5):
            if j!=i: prod*=(Vs[i][0][0]*Vs[j][0][1]-Vs[i][0][1]*Vs[j][0][0])
        Vs_ratio=n/prod
        if P==Ps[0]: r0=[None]*5
    # ratio n_i/prod_j det_ij must be independent of P
ratios={}
for P in Ps:
    Vs=[V(i,P) for i in range(5)]
    for i in range(5):
        (X,Y),y=Vs[i]; n=X*X+Y*Y-mu*y*y; prod=1
        for j in range(5):
            if j!=i: prod*=(X*Vs[j][0][1]-Y*Vs[j][0][0])
        ratios.setdefault(i,[]).append(n/prod)
for i in range(5):
    r=ratios[i]; print('i',i,'n_i/prod det: rel spread',mp.nstr(max(abs(x/r[0]-1) for x in r),5))
# determinants linear in P?
for (i,j) in [(0,1),(2,4)]:
    d=[ (lambda A,B: A[0][0]*B[0][1]-A[0][1]*B[0][0])(V(i,P),V(j,P)) for P in Ps[:4]]
    # second difference with equal spacing
    print('det',i,j,'2nd diff',mp.nstr(abs(d[2]-2*d[1]+d[0]),5),' |d|',mp.nstr(abs(d[0]),5))
