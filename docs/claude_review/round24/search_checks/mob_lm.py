# LM in the Moebius/gamma formulation. zhat_5 = id, zhat_a = rotation about P_a by theta_a (a=1..4).
import numpy as np, itertools, sys
from scipy.optimize import least_squares
from mob import rot_about, PAIRS
def mats(x):
    Ms=[rot_about(complex(x[3*a],np.exp(x[3*a+1])),x[3*a+2]) for a in range(4)]+[np.eye(2)]
    return Ms
def coinc_pts(Ma,Mb):
    # roots of det(Ma v, Mb v), v=(s,1): quadratic a s^2 + b s + c
    # det([Ma v, Mb v]) = (Ma v)_0 (Mb v)_1 - (Ma v)_1 (Mb v)_0
    a0,a1=Ma[0],Ma[1]; b0,b1=Mb[0],Mb[1]
    # (a0[0]s+a0[1])(b1[0]s+b1[1]) - (a1[0]s+a1[1])(b0[0]s+b0[1])
    c2=a0[0]*b1[0]-a1[0]*b0[0]; c1=a0[0]*b1[1]+a0[1]*b1[0]-a1[0]*b0[1]-a1[1]*b0[0]; c0=a0[1]*b1[1]-a1[1]*b0[1]
    return np.roots([c2,c1,c0])
def build(x,signs):
    Ms=mats(x); rows=[]; rhos=[]
    for (a,b),sg in zip(PAIRS,signs):
        r=coinc_pts(Ms[a],Ms[b])
        if len(r)<2: r=np.array([r[0],np.conj(r[0])]) if len(r)==1 else np.array([1j,-1j])
        rho=r[0] if ((r[0].imag>0)==(sg>0)) else r[1]
        v=Ms[a]@np.array([rho,1.0])
        pw=np.array([rho**k for k in range(4,-1,-1)])
        rows.append(np.concatenate([pw*v[1],-pw*v[0]])); rhos.append(rho)
    return Ms,np.array(rows),np.array(rhos)
def solve_g(A):
    # normalise coefficient of s^4 in g1 = 1: A[:,0]*1 + A[:,1:] c = 0
    c,*_=np.linalg.lstsq(A[:,1:],-A[:,0],rcond=None)
    g=np.concatenate([[1.0],c]); res=A@g
    return g[:5],g[5:],res
def Delta(g1,g2):
    p=np.polymul(g1,np.conj(g2)); return p.imag   # deg 8, 9 coeffs (high->low)
LSIGN=1.0
def residual(y,signs,L,c0):
    x=y[:12]; ll,s1,s2=y[12:15]; lam=LSIGN*np.exp(ll)
    Ms,A,rhos=build(x,signs); g1,g2,res=solve_g(A)
    D=Delta(g1,g2); target=lam*np.polynomial.polynomial.polypow([s2,-s1,1.0],4)[::-1]
    D=np.concatenate([np.zeros(9-len(D)),D])
    return np.concatenate([res.real,res.imag,(D-target),0.1*(L@y-c0)])
def places(y,signs):
    x=y[:12]; Ms,A,rhos=build(x,signs); g1,g2,res=solve_g(A)
    verts=[]
    for a in range(5):
        M=Ms[a]  # zhat_a(s) = M (s,1); h_a = det(zhat_a, g) = (M0 . (s,1)) g2 - (M1 . (s,1)) g1
        h=np.polysub(np.polymul([M[0,0],M[0,1]],g2),np.polymul([M[1,0],M[1,1]],g1))
        rts=np.roots(h)
        edge=[rhos[i] for i,(p,q) in enumerate(PAIRS) if a in (p,q)]
        rem=[r for r in rts if min(abs(r-e) for e in edge)>1e-6]
        verts.append(rem)
    return rhos,verts,g1,g2
if __name__=='__main__':
    N=int(sys.argv[1]); rng=np.random.default_rng(int(sys.argv[2])); out=open(sys.argv[3],'w')
    hits=0
    for it in range(N):
        signs=rng.choice([1,-1],size=10)
        x0=np.concatenate([[rng.normal()*0.3,rng.normal()*0.3,rng.uniform(0.3,6.0)] for _ in range(4)])
        y0=np.concatenate([x0,[rng.normal(),rng.normal()*2,rng.normal()-2]])
        globals()['LSIGN']=rng.choice([1.0,-1.0])
        L=rng.normal(size=(4,15)); c0=L@y0
        try:
            sol=least_squares(residual,y0,args=(signs,L,c0),method='lm',xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=4000)
        except Exception as e: continue
        r=residual(sol.x,signs,L,c0); rn=np.linalg.norm(r[:-4])
        if not np.isfinite(rn) or rn>1e-9: continue
        ll,s1,s2=sol.x[12:15]; lam=LSIGN*np.exp(ll); disc=s1*s1-4*s2
        rhos,verts,g1,g2=places(sol.x,signs)
        vflat=[v for vv in verts for v in vv]
        allp=np.concatenate([rhos,np.array(vflat)])
        ok_v=all(len(v)==1 for v in verts)
        allv=np.concatenate([allp,np.conj(allp)])
        iu=np.triu_indices(len(allv),1); sep=np.min(np.abs(allv[:,None]-allv[None,:])[iu]) if len(allv)>1 else 0
        imin=np.min(np.abs(allp.imag))
        good= ok_v and disc>1e-8 and sep>1e-4 and imin>1e-4 and abs(lam)>1e-8
        hits+=good
        print(it,f"res={rn:.1e} disc={disc:.2e} lam={lam:.2e} sep={sep:.1e} imin={imin:.1e} nverts={[len(v) for v in verts]}",'GOOD' if good else '', (list(sol.x),list(signs)) if good else '',file=out,flush=True)
    print('hits',hits,file=out)
