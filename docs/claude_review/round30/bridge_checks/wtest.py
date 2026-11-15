# Numerical membership test: do points of the cube 5-subset curves lie on round 23's hypersurface W (any S5 labelling)?
# W = closure of {[F_0:..:F_4](a,b,c,x)}: F_0=(x-a)(x-b)(x-c)(x+d), d=a+b+c, F_1=x^4-(S/2)x^2+m+2xz, z^2=m x^2+B,
# S=a^2+b^2+c^2+d^2, m=abcd, B=(e_3-S m)/4 with e_3 = e_3(a^2,b^2,c^2,d^2) (round23/split.md Thm 2).
# Validation: points of W itself must give residual ~0.
import sys, itertools, random, cmath
sys.path.insert(0,'/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round30/bridge_checks')
import numpy as np
from scipy.optimize import least_squares
from cube_basic import cube
from glib import peval, G
from fractions import Fraction as Fr

def Fvec(a,b,c,x,branch):
    d=a+b+c; sq=[a*a,b*b,c*c,d*d]; S=sum(sq); m=a*b*c*d
    e3=sq[0]*sq[1]*sq[2]+sq[0]*sq[1]*sq[3]+sq[0]*sq[2]*sq[3]+sq[1]*sq[2]*sq[3]
    B=(e3-S*m)/4; z=branch*cmath.sqrt(m*x*x+B)
    F0=(x-a)*(x-b)*(x-c)*(x+d); F1=x**4-(S/2)*x*x+m+2*x*z
    F2=(x*x-b*b)*(x*x-c*c); F3=(x*x-a*a)*(x*x-c*c); F4=(x*x-a*a)*(x*x-b*b)
    return np.array([F0,F1,F2,F3,F4])

def residual(p, tries=30, seed=1):
    rng=random.Random(seed); best=1e9
    p=np.array(p,dtype=complex); p=p/p[2]
    for _ in range(tries):
        x0=[rng.uniform(-3,3) for _ in range(6)]
        for branch in (1,-1):
            def f(v):
                a=v[0]+1j*v[1]; b=v[2]+1j*v[3]; c=v[4]+1j*v[5]
                F=Fvec(a,b,c,1.0,branch)
                if abs(F[2])<1e-14: return np.ones(8)*1e3
                r=F/F[2]-p
                r=np.delete(r,2)
                return np.concatenate([r.real,r.imag])
            try:
                s=least_squares(f,x0,xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=4000)
                best=min(best,np.linalg.norm(s.fun))
            except Exception: pass
        if best<1e-10: break
    return best

if __name__=="__main__":
    # validation on W
    F=Fvec(1.3+0.2j,0.7-0.4j,2.1+0.1j,1.0,1)
    print("validation residual (point of W):",residual(F,tries=40))
    F=np.array([1,1.3+0.2j,0.8-0.1j,1.7+0.5j,0.4+0.9j])
    print("random point residual:",residual(F,tries=40))
    y=[Fr(-9,4),Fr(-2,3),Fr(1,4),Fr(8,3)]
    roots,P=cube(y); names=list(P)
    ts=[1.37,-2.21,0.6+0.8j]
    worst=1e9; cnt=0
    for S in itertools.combinations(names,5):
        for perm in itertools.permutations(S):
            if not (perm[2]<perm[3]<perm[4]): continue   # W is symmetric in F2,F3,F4 (S3 on a,b,c)
            rmin=1e9
            for t in ts:
                p=[complex(float(peval(P[n][0],G(Fr(t.real).limit_denominator(1000))).a),0) if False else None for n in perm]
                tt=t
                vals=[]
                for n in perm:
                    coeffs=[complex(float(c.a),float(c.b)) for c in P[n][0]]
                    vals.append(sum(cc*tt**k for k,cc in enumerate(coeffs)))
                r=residual(vals,tries=6)
                rmin=min(rmin,r)
                if r>1e-6: break     # one point off W suffices for this labelling
            worst=min(worst,rmin); cnt+=1
            if cnt%100==0: print('progress',cnt,'min so far',worst,flush=True)
            if rmin<1e-8: print("ON W?",S,perm,rmin,flush=True)
    print("labellings tested:",cnt," smallest residual over all labellings (first point that failed):",worst)
