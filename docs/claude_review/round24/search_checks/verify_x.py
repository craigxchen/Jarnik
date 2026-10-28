# high-precision refinement and full verification of x-structured solutions found by xstruct.py
import numpy as np, ast, itertools, sys
import mpmath as mp
mp.mp.dps=60
import xstruct as X
S3=X.S3; MEM=X.MEM
def polymul(a,b):
    out=[mp.mpc(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): out[i+j]+=x*y
    return out
def poly_from_roots(rs):
    p=[mp.mpc(1)]
    for r in rs: p=polymul(p,[mp.mpc(1),-r])
    return p
def Fmp(y,sg):
    rho=[mp.mpc(y[k],y[7+k]) for k in range(7)]; r4=mp.mpc(y[14],y[15]); sig=sg*mp.e**y[16]
    beta=[mp.mpf(1),y[17],y[18],y[19]]; q=[y[20],y[21],y[22]]
    out=[]
    for j in range(3):
        c=poly_from_roots([rho[m] for m in MEM[j]]); out+=[mp.im(c[k]) for k in (1,2,3)]
    n=poly_from_roots(rho+[r4]); out+=[mp.im(n[k]) for k in (1,2,3)]
    In=[mp.im(x) for x in n]; Rn=[mp.re(x) for x in n]
    In2=polymul(In,In)[-9:]
    lhs=[Rn[k]-mp.re(In2[k])/(4*sig) for k in range(9)]; lhs[8]+=sig
    rhs=polymul(polymul(beta,beta),q); rhs=[mp.mpf(0)]*(9-len(rhs))+rhs
    out+=[lhs[k]-mp.re(rhs[k]) for k in range(9)]
    out+=[mp.re(rho[6]), mp.im(rho[6])-1]
    return out
def refine(y,sg):
    y=[mp.mpf(v) for v in y]
    for it in range(30):
        f=Fmp(y,sg); nf=mp.sqrt(sum(v*v for v in f))
        if nf<mp.mpf(10)**(-50): break
        h=mp.mpf(10)**(-25); J=mp.matrix(23,23)
        for k in range(23):
            yp=list(y); yp[k]+=h; ym=list(y); ym[k]-=h
            fp=Fmp(yp,sg); fm=Fmp(ym,sg)
            for i in range(23): J[i,k]=(fp[i]-fm[i])/(2*h)
        d=mp.lu_solve(J,mp.matrix([-v for v in f]))
        y=[y[k]+d[k] for k in range(23)]
    return y,nf
def verify(y,sg,label):
    y,nf=refine(y,sg)
    print(f'== {label}: residual after refinement {mp.nstr(nf,5)}')
    rho=[mp.mpc(y[k],y[7+k]) for k in range(7)]; r4=mp.mpc(y[14],y[15]); sig=sg*mp.e**y[16]
    beta=[mp.mpf(1),y[17],y[18],y[19]]; q=[y[20],y[21],y[22]]
    n=poly_from_roots(rho+[r4])
    In=[mp.im(x) for x in n]
    print('   Im n (x^8..x^0):',[mp.nstr(v,6) for v in In])
    # kappa = 1, c = sig, q_true = -q/sig, alpha_r = Im(n)/(2 sig)
    c=sig; qt=[-v/sig for v in q]; ar=[v/(2*sig) for v in In]
    # check (alpha_r + i)^2 - beta^2 q_true == n/c
    a=list(ar); a[-1]+=mp.mpc(0,1)
    lhs=polymul(a,a); b2q=polymul(polymul(beta,beta),qt); b2q=[0]*(len(lhs)-len(b2q))+b2q
    lhs=[lhs[k]-b2q[k] for k in range(len(lhs))]
    ncn=[0]*(len(lhs)-len(n))+[x/c for x in n]
    print('   identity (alpha+i)^2 - beta^2 q = n/c : max err',mp.nstr(max(abs(lhs[k]-ncn[k]) for k in range(len(lhs))),5))
    Q2,Q1,Q0=qt; print('   conic z^2 = q(x), leading coeff',mp.nstr(Q2,12),' (needs >0)  m-class ~ -1/sigma =',mp.nstr(-1/sig,12))
    r2=mp.sqrt(Q2)
    # parametrisation: x = (Q0 - w^2)/(2 r2 w - Q1), z = r2 x + w ; branches at w = inf and w0 = Q1/(2 r2)
    def wroots_over(xv):  # points over x = xv: Q0 - w^2 - xv (2 r2 w - Q1) = 0 -> w^2 + 2 r2 xv w - (Q0 + xv Q1) = 0
        return mp.polyroots([1,2*r2*xv,-(Q0+xv*Q1)])
    places={}  # w -> vector
    def add(w,j,o):
        key=None
        for kk in places:
            if abs(kk-w)<mp.mpf(10)**(-12): key=kk
        if key is None: places[w]=[0]*5; key=w
        places[key][j]+=o
    for i,S in enumerate(S3):
        for w in wroots_over(rho[i]):
            for j in S: add(w,j-1,+1)      # zero of P_j
        for w in wroots_over(mp.conj(rho[i])):
            for j in S: add(w,j-1,-1)      # pole of P_j (zero of Pbar_j)
    # u_4 = H/Hbar ; H(w)*(2 r2 w - Q1)^4 polynomial in w
    def Hpoly(sgn):
        den=[2*r2,-Q1]; num_x=[-1,0,Q0]   # x = num_x/den (num_x in w: -w^2 + Q0)
        # build alpha_r(x) + sgn*i + beta(x) z, z = r2 x + w, all times den^4
        def xpow(k):  # x^k * den^4 as poly in w
            p=[mp.mpc(1)]
            for _ in range(k): p=polymul(p,num_x)
            for _ in range(4-k): p=polymul(p,den)
            return p
        tot=[mp.mpc(0)]*9
        def addp(p,coef):
            p=[0]*(9-len(p))+p
            for k in range(9): tot[k]+=coef*p[k]
        A=list(ar); A=[0]*(5-len(A))+A   # degree<=4 coefficient list high->low
        for d in range(5): addp(xpow(4-d) if False else xpow(d),A[4-d])
        addp(xpow(0),sgn*mp.mpc(0,1))
        # beta(x) z = beta(x) (r2 x + w): beta deg 3
        B=beta
        for d in range(4):
            coef=B[3-d]
            # x^d * (r2 x + w) * den^4 = r2 x^{d+1} den^4 + w x^d den^4
            addp(xpow(d+1),coef*r2)
            p=xpow(d) if d<=4 else None
            addp(polymul(p,[1,0])[-9:] if len(polymul(p,[1,0]))>9 else polymul(p,[1,0]),coef)
        return tot
    def Hval(xv,zv,sgn):
        A=sum(ar[k]*xv**(len(ar)-1-k) for k in range(len(ar))); Bv=sum(beta[k]*xv**(3-k) for k in range(4))
        return A+sgn*mp.mpc(0,1)+Bv*zv
    nH=0; nHb=0
    for r in rho+[r4]:
        for w in wroots_over(r):
            xv=(Q0-w**2)/(2*r2*w-Q1); zv=r2*xv+w
            if abs(Hval(xv,zv,+1))<mp.mpf(10)**(-20): add(w,3,+1); nH+=1
        for w in wroots_over(mp.conj(r)):
            xv=(Q0-w**2)/(2*r2*w-Q1); zv=r2*xv+w
            if abs(Hval(xv,zv,-1))<mp.mpf(10)**(-20): add(w,3,-1); nHb+=1
    print('   zeros of H found over the 8 x-values:',nH,' of Hbar:',nHb,' (H has exactly 8 zeros on the conic)')
    vecs=[]
    for w,v in places.items():
        v=np.array(v,int); v=v-v.min(); vecs.append((w,tuple(v)))
    nz=[(w,v) for w,v in vecs if any(v)]
    print('   number of boundary places:',len(nz))
    from collections import Counter
    cnt=Counter()
    ok=True
    for w,v in nz:
        A=frozenset(i+1 for i in range(5) if v[i]==1)
        if max(v)!=1: ok=False; print('   non-ray vector',v)
        cnt[A]+=1
    cuts=Counter()
    for A,k in cnt.items():
        Ac=frozenset(range(1,6))-A; cuts[min(tuple(sorted(A)),tuple(sorted(Ac)))]+=k
    print('   all ray vectors:',ok,' distinct rays hit:',len(cnt),' each ray once:',all(k==1 for k in cnt.values()))
    print('   cut weights:',sorted(set(cuts.values())),' number of cuts:',len(cuts))
    # branch points: w = inf and w0 = Q1/(2 r2): check they are not places
    w0=Q1/(2*r2); print('   min distance of places to branch w0:',mp.nstr(min(abs(w-w0) for w,_ in nz),5), ' max |w|:',mp.nstr(max(abs(w) for w,_ in nz),5))
    # contact: u_j - 1 order at branches, numerically
    def u_all(wv):
        xv=(Q0-wv**2)/(2*r2*wv-Q1); zv=r2*xv+wv
        us=[]
        for j in (1,2,3):
            P=mp.mpc(1)
            for i,S in enumerate(S3):
                if j in S: P*= (xv-rho[i])
            Pb=mp.mpc(1)
            for i,S in enumerate(S3):
                if j in S: Pb*= (xv-mp.conj(rho[i]))
            us.append(P/Pb)
        A=sum(ar[k]*xv**(len(ar)-1-k) for k in range(len(ar))); Bv=sum(beta[k]*xv**(3-k) for k in range(4))
        us.append((A+mp.mpc(0,1)+Bv*zv)/(A-mp.mpc(0,1)+Bv*zv))
        return us
    for name,ws in (('w0',[w0+mp.mpf(10)**(-k) for k in (6,8)]),('inf',[mp.mpf(10)**k for k in (6,8)])):
        a=u_all(ws[0]); b=u_all(ws[1])
        slopes=[mp.log(abs(b[j]-1)/abs(a[j]-1))/mp.log(abs((ws[1]-(w0 if name=='w0' else 0))/(ws[0]-(w0 if name=='w0' else 0)))) for j in range(4)]
        print('   contact orders at',name,':',[mp.nstr(abs(s),6) for s in slopes])
    # phases: u_j(branch) -> 1 ?
    print('   u_j at w = w0+1e-20:',[mp.nstr(v,8) for v in u_all(w0+mp.mpf(10)**(-20))])
    return y
if __name__=='__main__':
    for line in open(sys.argv[1]):
        if 'NONDEG' not in line: continue
        y=np.array(ast.literal_eval(line[line.index('['):].strip().replace('np.float64','')))
        sg=None
        for s_ in (1.0,-1.0):
            X.SIGN=s_
            if np.linalg.norm(X.F(y))<1e-9: sg=s_; break
        verify(list(y),sg,'solution '+line.split()[0])
