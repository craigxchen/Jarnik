# Floating-point sanity check of the two atom bounds for lattice Gaussians (not a proof).
import numpy as np, math, random
def dist(A,cut=12):
    n=np.arange(-int(cut*A),int(cut*A)+1); p=np.exp(-n**2/(2*A*A)); return n,p/p.sum()
def atom(w,A):
    n,p=dist(A); acc={0:1.0}
    lo,hi=0,0; vec=np.array([1.0]); off=0
    for wi in w:
        # distribution of wi*xi on support
        m=abs(wi); L=len(n); arr=np.zeros(m*(L-1)+1); arr[::m]=p
        vec=np.convolve(vec,arr); off+= -m*n[-1]
    return vec.max(), vec
for A in (3.0,6.0):
    for w in ([1]*5,[1]*20,[2,4,6],[2,2,2,2],[1,3],[3,5,7,9],[1,-1,2,-2,1]):
        mx,_=atom(w,A); n0=len(w); nrm=math.sqrt(sum(x*x for x in w)); g=math.gcd(*w)
        print(f"A={A} w={w}: atom={mx:.4f}  1/(A sqrt(2pi n0))={1/(A*math.sqrt(2*math.pi*n0)):.4f}  g/(A|w|sqrt(2pi))={g/(A*nrm*math.sqrt(2*math.pi)):.4f}")
