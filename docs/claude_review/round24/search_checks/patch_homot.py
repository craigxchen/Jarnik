s=open('homot.py').read()
new_track='''def track(x0,eps,rng,out,maxsteps=6000):
    EN0,_=EJ(x0,eps,KN); sl=rng.normal(size=30); s0=sl@x0
    z=np.concatenate([x0,[0.0]]); h=0.01; prev=None
    def HJ(z):
        H,Jx,Jt=Hfun(z[:30],z[30],eps,EN0,sl,s0); return H,np.hstack([Jx,Jt[:,None]])
    for step in range(maxsteps):
        H,J=HJ(z)
        _,_,Vt=np.linalg.svd(J); t=Vt[-1]
        if prev is None:
            if t[30]<0: t=-t
        elif t@prev<0: t=-t
        zp=z+h*t; ok=False
        for it in range(10):
            H,J=HJ(zp)
            A=np.vstack([J,t[None,:]]); b=-np.concatenate([H,[t@(zp-z)-h]])
            try: d=np.linalg.solve(A,b)
            except np.linalg.LinAlgError: break
            zp=zp+d
            if np.linalg.norm(d)<1e-11*(1+np.linalg.norm(zp)): ok=True; break
        if ok and np.all(np.isfinite(zp)) and np.linalg.norm(HJ(zp)[0])<1e-9:
            if zp[30]>=1.0 and z[30]<1.0:
                x=z[:30]+(zp[:30]-z[:30])*(1.0-z[30])/(zp[30]-z[30])
                for it in range(20):
                    Hh,Jx,Jt=Hfun(x,1.0,eps,EN0,sl,s0); d=np.linalg.solve(Jx,-Hh); x=x+d
                    if np.linalg.norm(d)<1e-13: break
                return 'reached',x,1.0
            prev=t; z=zp; h=min(h*1.4,0.05)
            r=z[:15]+1j*eps*np.exp(z[15:30]); d=diag(r)
            if d['sep']<1e-6 or d['imin']<1e-7 or d['rmax']>1e6 or d['rmin']<1e-6: return 'degenerate',z[:30],z[30]
            if z[30]<-0.5: return 'turned_back',z[:30],z[30]
        else:
            h/=3
            if h<1e-11: return 'stalled',z[:30],z[30]
    return 'maxsteps',z[:30],z[30]
'''
i=s.index('def track('); j=s.index("if __name__=='__main__':")
s=s[:i]+new_track+s[j:]
open('homot.py','w').write(s)
