import mpmath as mp, itertools, sys
mp.mp.dps=int(sys.argv[2]) if len(sys.argv)>2 else 60
def build(m):
    m=mp.mpf(m)
    c=m**3-11*m**2+39*m-13
    u=-(m+3)*c/(16*(m-3)**2*(m-1)**2); v=-(m+1)*c/(4*(m-3)**3*(m-1))
    bs=mp.polyroots([1,0,v,-v**2,u*v**2],maxsteps=500,extraprec=200)
    bs=[mp.re(b) if abs(mp.im(b))<mp.mpf(10)**(-mp.mp.dps+10) else b for b in bs]
    return u,v,bs
def curve(m):
    u,v,bs=build(m)
    assert all(isinstance(b,mp.mpf) for b in bs), 'b not real'
    T2=bs[0]*bs[1]*bs[2]*bs[3]; t=mp.sqrt(T2)
    p=[t/b+1j*b for b in bs]
    lab={(1,2,3):p[0],(4,):-mp.conj(p[0]),(1,):p[1],(2,):p[2],(3,):p[3],(2,3):-mp.conj(p[1]),(1,3):-mp.conj(p[2]),(1,2):-mp.conj(p[3])}
    def pmul(a,b):
        r=[0]*(len(a)+len(b)-1)
        for i,x in enumerate(a):
            for j,y in enumerate(b): r[i+j]+=x*y
        return r
    n=[mp.mpc(1)]
    for r in range(4): n=pmul(n,[1,-p[r]]); n=pmul(n,[1,mp.conj(p[r])])
    sig=-mp.fprod([abs(x)**2 for x in p])
    ren=[mp.re(x) for x in n]; imn=[mp.im(x) for x in n]   # coefficients highest first, degree 8
    im2=pmul(imn,imn)[-9:]
    L=[ren[i]-im2[i]/(4*sig) for i in range(9)]; L[-1]+=sig
    # L even: coefficients of X^8,X^6,..; L = Y(Y+e)^2(Y+D)
    cub=[L[0],L[2],L[4],L[6]]
    rts=mp.polyroots(cub,maxsteps=500,extraprec=200)
    rts=sorted(rts,key=lambda z:mp.re(z))
    # find double root
    best=None
    for i,j in itertools.combinations(range(3),2):
        d=abs(rts[i]-rts[j])
        if best is None or d<best[0]: best=(d,i,j)
    d,i,j=best; k=3-i-j
    e=-(rts[i]+rts[j])/2; D=-rts[k]
    return dict(m=m,u=u,v=v,bs=bs,t=t,p=p,lab=lab,n=n,sig=sig,imn=imn,L=L,e=e,D=D,dbl=d,odd=max(abs(L[1]),abs(L[3]),abs(L[5]),abs(L[7]),abs(L[8])))
def peval(cf,x):
    return mp.polyval(cf,x)
def funcs(C):
    lab,sig,e,D,imn=C['lab'],C['sig'],C['e'],C['D'],C['imn']
    kap=-1/sig; sk=mp.sqrt(kap)
    def Xz(w): return (kap*D/w-w)/(2*sk),(kap*D/w+w)/2
    def P(j,X,conj=False):
        r=mp.mpc(1)
        for S,x in lab.items():
            if j in S: r*= (X-(mp.conj(x) if conj else x))
        return r
    def H(X,z,conj=False):
        return peval(imn,X)/(2*sig)+(-1j if conj else 1j)+X*(X**2+e)*z
    def U(w):
        X,z=Xz(w)
        return [P(j,X)/P(j,X,True) for j in (1,2,3)]+[H(X,z)/H(X,z,True)]
    return Xz,P,H,U,kap,sk
def places(C):
    Xz,P,H,U,kap,sk=funcs(C)
    D=C['D']; xs=[]
    for x in C['p']: xs+= [x,-mp.conj(x),mp.conj(x),-x]
    out=[]
    for x0 in xs:
        for w in mp.polyroots([1,2*sk*x0,-kap*D],extraprec=100):
            X,z=Xz(w)
            o=[]
            for j in (1,2,3):
                o.append(int(abs(P(j,X))<mp.mpf(10)**(-mp.mp.dps//2))-int(abs(P(j,X,True))<mp.mpf(10)**(-mp.mp.dps//2)))
            o.append(int(abs(H(X,z))<mp.mpf(10)**(-mp.mp.dps//2))-int(abs(H(X,z,True))<mp.mpf(10)**(-mp.mp.dps//2)))
            vec=o+[0]
            if max(vec)==min(vec): continue
            out.append((w,vec))
    return out
def ray(vec):
    lo=min(vec); hi=max(vec)
    assert hi-lo==1 and all(x in (lo,hi) for x in vec), vec
    return tuple(i+1 for i,x in enumerate(vec) if x==hi)
def check(m):
    C=curve(m)
    pl=places(C)
    rays=[ray(v) for w,v in pl]
    allA=[A for r in range(1,5) for A in itertools.combinations(range(1,6),r)]
    ok_rays = sorted(rays)==sorted(allA)
    ws=[w for w,v in pl]
    sep=min(abs(a-b) for a,b in itertools.combinations(ws,2))
    realdist=min(abs(mp.im(w)) for w in ws)
    Xz,P,H,U,kap,sk=funcs(C)
    cont=[]
    for W0 in (mp.mpf(10)**8,mp.mpf(10)**16, mp.mpf(10)**-8, mp.mpf(10)**-16):
        cont.append([abs(x-1) for x in U(W0)])
    sl_inf=[float(mp.log10(cont[1][j]/cont[0][j])/8) for j in range(4)]
    sl_0=[float(mp.log10(cont[3][j]/cont[2][j])/(-8)) for j in range(4)]
    # also pairs
    Uinf=U(mp.mpf(10)**16); U0=U(mp.mpf(10)**-16)
    pairs_inf=min(float(-mp.log10(abs(Uinf[a]/Uinf[b]-1))/16) for a,b in itertools.combinations(range(4),2))
    return C,pl,dict(m=float(m),nplaces=len(pl),rays_ok=ok_rays,sep=mp.nstr(sep,5),dist_to_R=mp.nstr(realdist,5),
                     slope_inf=sl_inf,slope_0=sl_0,pair_inf_min=pairs_inf,dbl=mp.nstr(C['dbl'],3),odd=mp.nstr(C['odd'],3),
                     D=mp.nstr(C['D'],8),e=mp.nstr(C['e'],8),sig=mp.nstr(C['sig'],8))
if __name__=='__main__':
    for m in sys.argv[1].split(','):
        C,pl,r=check(mp.mpf(m)); print(r)
