# Moebius formulation: 5 real Moebius maps zhat_a(s) = (A_a s + B_a) with A_a,B_a in R^2 (columns); gamma=[g1:g2] deg 4.
import numpy as np, itertools
PAIRS=list(itertools.combinations(range(5),2))
def coinc(Ma,Mb):
    # Ma = 2x2 real matrix acting on (s,1): zhat(s) = Ma @ (s,1). coincidences: det(Ma v, Mb v)=0, v=(s,1): real quadratic
    # det([Ma v, Mb v]) = q(s)
    def q(s):
        v=np.array([s,1.0]); return np.linalg.det(np.column_stack([Ma@v,Mb@v]))
    # fit quadratic
    S=np.array([-1.0,0.0,1.0]); Y=np.array([q(x) for x in S]); c=np.polyfit(S,Y,2)
    return np.roots(c), c
def interp_matrix(Ms,signs):
    rows=[]; pts=[]
    for (a,b),sg in zip(PAIRS,signs):
        r,c=coinc(Ms[a],Ms[b])
        if abs(r[0].imag)<1e-12: return None,None
        rho=r[0] if (r[0].imag>0)==(sg>0) else r[1]
        v=Ms[a]@np.array([rho,1.0])   # complex 2-vector, direction w = v[0]/v[1]
        # gamma(rho) || v : g1(rho) v[1] - g2(rho) v[0] = 0 ; unknowns: coeffs of g1 (5), g2 (5)
        pw=np.array([rho**k for k in range(4,-1,-1)])
        rows.append(np.concatenate([pw*v[1],-pw*v[0]])); pts.append(rho)
    return np.array(rows),pts
if __name__=='__main__':
    rng=np.random.default_rng(0)
    for trial in range(5):
        Ms=[rng.normal(size=(2,2)) for _ in range(5)]
        ranks=[]
        for signs in itertools.product([1,-1],repeat=10):
            A,pts=interp_matrix(Ms,signs)
            if A is None: break
            sv=np.linalg.svd(A,compute_uv=False); ranks.append(sv[-1]/sv[0])
        if ranks: print(trial,'min rel sigma_min over sign choices',min(ranks),'max',max(ranks))
        else: print(trial,'real coincidences present')

def rot_about(P,theta):
    # real Moebius matrix: rotation by theta about P in upper half plane
    x,y=P.real,P.imag
    A=np.array([[np.sqrt(y),x/np.sqrt(y)],[0,1/np.sqrt(y)]])   # maps i -> P
    R=np.array([[np.cos(theta/2),np.sin(theta/2)],[-np.sin(theta/2),np.cos(theta/2)]])
    return A@R@np.linalg.inv(A)
def elliptic_all(Ms):
    for a,b in PAIRS:
        M=np.linalg.inv(Ms[b])@Ms[a]; tr=np.trace(M); dt=np.linalg.det(M)
        if tr*tr>=4*dt: return False
    return True
def sample_elliptic(rng,spread=0.3):
    while True:
        Ms=[rot_about(complex(rng.normal()*spread, np.exp(rng.normal()*spread)), rng.uniform(0.3,2*np.pi-0.3)) for _ in range(5)]
        if elliptic_all(Ms): return Ms
