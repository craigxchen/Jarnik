# Independent re-check of Prop 4.3 / Theorem 4.4 (bit strata), Sylvester n=2^k, r=n-1.
import mpmath as mp, math, sys
mp.mp.dps=50
def params(n,r):
    d=n-1
    w=mp.sqrt(8*d*mp.log(4*r)); j=int(mp.floor(w/2)); m=int(mp.ceil(4*w))
    return d,w,j,m
def lnfact_sum(s):
    # sum_{t=1}^{s-1} (s-t) ln t  exactly via mpmath (s up to 2^22 -> use loggamma: sum_{t=1}^{s-1} ln t! ... )
    # sum_{t=1}^{s-1}(s-t)ln t = sum_{u=1}^{s-1} ln(u!) = sum_{u=1}^{s-1} loggamma(u+1)
    # use Barnes G: ln G(s+1) = sum_{u=1}^{s-1} ln(u!)
    return mp.log(mp.barnesg(s+1)) if s<200 else mp.loggamma(1)*0 + mp.mpf(mp.log(mp.barnesg(mp.mpf(s+1))))
def run(k, theta_kind, Wfac=None, Wmin=False):
    n=2**k; r=n-1; d,w,j,m=params(n,r)
    W = 3*n*mp.log(n) - (14*n if Wmin else 0) if Wfac is None else Wfac*n*mp.log(n)
    eps0 = d*mp.e**(-W/4)/(m+1)
    abar = W/(2*d); N2=mp.mpf(n)*(n-1)/2
    if theta_kind=='half': one_minus=mp.mpf(1)/2
    elif theta_kind=='t0_author': one_minus=mp.mpf(1)/4*(4*r)**mp.mpf(-0.5)
    elif theta_kind=='t0_correct': one_minus=mp.mpf(1)/2*(4*r)**mp.mpf(-0.5)
    elif theta_kind=='tS': one_minus=2*r*mp.e**(-3*w**2/(8*n))
    out=[]
    for q in range(1,k+1):
        # minimal c>=0 with q*2*exp(-3(c+1)^2/(2n)) <= one_minus
        c=0
        while q*2*mp.e**(-3*mp.mpf(c+1)**2/(2*n))>one_minus: c+=1
        rho=2*c*q
        if rho>m-1: continue
        o=m-rho; K=2**q; s=n//K
        ln_eta=(mp.log(eps0)-m*mp.log(4))/o
        eta=mp.e**ln_eta
        okEta = eta<=mp.mpf(1)/(2*(K-1)) if K>1 else True
        S = K*(mp.mpf(s)*(s-1)/2*mp.log(s*(1+eta)/eta) - lnfact_sum(s)) + (N2-K*mp.mpf(s)*(s-1)/2)*mp.log(2*(1+eta)*(K-1) if K>1 else 1)
        Gam = 2*S/N2
        out.append((q,c,rho,o,float(ln_eta),bool(okEta),float(Gam),float(abar), bool(okEta and Gam<=abar)))
    return dict(n=n,m=m,j=j,w=float(w),rows=out)
if __name__=='__main__':
    kinds=sys.argv[1].split(',')
    ks=[int(x) for x in sys.argv[2].split(',')]
    Wfac=None if sys.argv[3]=='3' else float(sys.argv[3])
    Wmin = (len(sys.argv)>4 and sys.argv[4]=='Wmin')
    for th in kinds:
        for k in ks:
            res=run(k,th,Wfac if Wfac!=3 else None,Wmin)
            fails=[r[0] for r in res['rows'] if r[-1]]
            adm=[r[0] for r in res['rows']]
            print(th,'k=',k,'m=',res['m'],'j=',res['j'],'admissible q:',adm,'failing q (in MIN region):',fails, flush=True)
