# build the curve from the rational parameter m and verify it independently (places, vectors, contact)
import numpy as np, sympy as sp, mpmath as mp, itertools, sys
mp.mp.dps=50
def build(mv, check=True):
    m=mp.mpf(mv)
    cub=m**3-11*m**2+39*m-13
    u=-(m+3)*cub/(16*(m-3)**2*(m-1)**2); v=-(m+1)*cub/(4*(m-3)**3*(m-1))
    if not (v<0 and u>0): return None,'need v<0,u>0 (u='+mp.nstr(u,5)+', v='+mp.nstr(v,5)+')'
    s=1; e2=v*s**2; e3=v**2*s**3; e4=u*v**2*s**4      # e3 sign choice +
    bs=mp.polyroots([1,0,e2,-e3,e4],maxsteps=200,extraprec=200)
    if any(abs(mp.im(b))>mp.mpf(10)**-30 for b in bs): return None,'b-quartic has non-real roots'
    bs=[mp.re(b) for b in bs]
    t=mp.sqrt(e4)
    p=[t/b+1j*b for b in bs]
    # assign roles: p_0 = rho_123 (partner rho_4), p_1..p_3 = rho_1..rho_3 (partners rho_23, rho_13, rho_12)
    return p,'ok'
def verify(p):
    # Thm 4.1 check in X, then u_4 block, then places on the conic via parametrisation, as in verify_x
    mb=lambda z: -mp.conj(z)
    rho={frozenset([1,2,3]):p[0],frozenset([1]):p[1],frozenset([2]):p[2],frozenset([3]):p[3],
         frozenset([2,3]):mb(p[1]),frozenset([1,3]):mb(p[2]),frozenset([1,2]):mb(p[3])}
    r4=mb(p[0])
    def pfr(rs):
        out=[mp.mpc(1)]
        for r in rs:
            new=[mp.mpc(0)]*(len(out)+1)
            for i,c in enumerate(out): new[i]+=c; new[i+1]-=c*r
            out=new
        return out
    thm41=max(abs(mp.im(pfr([rho[S] for S in rho if j in S])[k])) for j in (1,2,3) for k in (1,2,3))
    n=pfr(list(rho.values())+[r4]); In=[mp.im(c) for c in n]; Rn=[mp.re(c) for c in n]
    sig=-mp.fprod([abs(z)**2 for z in p])
    # ell(Y): L = Re n - Im n^2/(4 sig) + sig, even in X
    Im2=[mp.mpf(0)]*17
    for i,a in enumerate(In):
        for j,b in enumerate(In): Im2[i+j]+=a*b
    L=[Rn[k]-Im2[8+k]/(4*sig) for k in range(9)]; L[8]+=sig
    odd=max(abs(L[k]) for k in range(1,9,2))
    ell=[L[0],L[2],L[4],L[6],L[8]]  # coefficients of Y^4..Y^0
    rts=mp.polyroots(ell,maxsteps=200,extraprec=200)
    rts=sorted(rts,key=lambda z:abs(z))
    # expect roots: 0, -e (double), -D
    zero=rts[0]; others=rts[1:]
    # find double root
    pairs=sorted([(abs(a-b),a,b) for a,b in itertools.combinations(others,2)],key=lambda x:x[0])
    dbl=(pairs[0][1]+pairs[0][2])/2; single=[z for z in others if abs(z-pairs[0][1])>1e-20 and abs(z-pairs[0][2])>1e-20][0]
    e=-mp.re(dbl); D=-mp.re(single)
    print('   Thm4.1 residual',mp.nstr(thm41,3),' odd part of L',mp.nstr(odd,3),' |l(0)|',mp.nstr(abs(zero),3),' double-root gap',mp.nstr(pairs[0][0],3),' e =',mp.nstr(e,8),' D =',mp.nstr(D,8))
    return rho,r4,sig,e,D,In
if __name__=='__main__':
    for mv in sys.argv[1:]:
        p,msg=build(mv); print('m =',mv,':',msg)
        if p is None: continue
        verify(p)

def full_check(rho,r4,sig,e,D,In):
    # conic z^2 = qt(X) = -(X^2 + D)/sig ; H = alpha_r + i + beta z, alpha_r = Im n/(2 sig), beta = X^3 + e X
    Q2=-1/sig; Q1=mp.mpf(0); Q0=-D/sig
    if Q2<=0: print('   conic has no real points at infinity'); return
    r2=mp.sqrt(Q2); ar=[c/(2*sig) for c in In]; beta=[1,0,e,0]
    def wroots_over(xv): return mp.polyroots([1,2*r2*xv,-(Q0+xv*Q1)])
    def Hval(xv,zv,sgn):
        A=sum(ar[k]*xv**(len(ar)-1-k) for k in range(len(ar))); Bv=sum(beta[k]*xv**(3-k) for k in range(4))
        return A+sgn*mp.mpc(0,1)+Bv*zv
    pl={}
    def add(w,j,o):
        for kk in list(pl):
            if abs(kk-w)<mp.mpf(10)**(-15): pl[kk][j]+=o; return
        pl[w]=[0]*5; pl[w][j]+=o
    for S,r in rho.items():
        for w in wroots_over(r):
            for j in S: add(w,j-1,1)
        for w in wroots_over(mp.conj(r)):
            for j in S: add(w,j-1,-1)
    nH=0
    for r in list(rho.values())+[r4]:
        for w in wroots_over(r):
            xv=(Q0-w**2)/(2*r2*w-Q1)
            if abs(Hval(xv,r2*xv+w,1))<mp.mpf(10)**(-15): add(w,3,1); nH+=1
        for w in wroots_over(mp.conj(r)):
            xv=(Q0-w**2)/(2*r2*w-Q1)
            if abs(Hval(xv,r2*xv+w,-1))<mp.mpf(10)**(-15): add(w,3,-1)
    from collections import Counter
    cnt=Counter(); ok=True
    for w,v in pl.items():
        v=np.array(v); v=v-v.min()
        if not any(v): continue
        if max(v)!=1: ok=False
        cnt[frozenset(i+1 for i in range(5) if v[i]==1)]+=1
    ws=list(pl.keys()); sepw=min(abs(a-b) for i,a in enumerate(ws) for b in ws[i+1:])
    print('   min separation of the 30 places (w-coordinate):',mp.nstr(sepw,4))
    print('   zeros of H:',nH,' places:',sum(cnt.values()),' rays hit:',len(cnt),' all once:',all(c==1 for c in cnt.values()),' ray vectors:',ok)
    w0=Q1/(2*r2)
    def us(wv):
        xv=(Q0-wv**2)/(2*r2*wv-Q1); zv=r2*xv+wv; out=[]
        for j in (1,2,3):
            P=mp.fprod([xv-r for S,r in rho.items() if j in S]); Pb=mp.fprod([xv-mp.conj(r) for S,r in rho.items() if j in S]); out.append(P/Pb)
        out.append(Hval(xv,zv,1)/Hval(xv,zv,-1)); return out
    for name,ws,c0 in (('w0',[w0+mp.mpf(10)**-6,w0+mp.mpf(10)**-8],w0),('inf',[mp.mpf(10)**6,mp.mpf(10)**8],0)):
        a=us(ws[0]); b=us(ws[1])
        print('   contact at',name,[mp.nstr(abs(mp.log(abs(b[j]-1)/abs(a[j]-1))/mp.log(abs((ws[1]-c0)/(ws[0]-c0)))),4) for j in range(4)])
if __name__=='__main__':
    for mv in sys.argv[1:]:
        p,msg=build(mv)
        if p is None: continue
        print('full check m =',mv); full_check(*verify(p))
