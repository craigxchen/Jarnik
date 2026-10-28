"""Check of Prop 2.2 (spectral LP cap): for S in a slice of Z^n and any profile,
   rho(S) <= Omega(S) <= sum_p om_p h_S(beta_p) |G^+ bt_p|_2 <= (sqrt(n)/lam) Lam(S).
rho(S) = max{j : V(S,j) != 0} computed EXACTLY (Fraction Gaussian elimination);
Omega(S) exact rational (enumeration of eps); lam, |G^+ bt_p| floating (evidence)."""
import itertools, random, math
from fractions import Fraction
import numpy as np

def rank_exact(rows):
    M=[[Fraction(x) for x in r] for r in rows]
    rk=0; ncol=len(M[0]) if M else 0
    for c in range(ncol):
        piv=None
        for i in range(rk,len(M)):
            if M[i][c]!=0: piv=i;break
        if piv is None: continue
        M[rk],M[piv]=M[piv],M[rk]
        for i in range(len(M)):
            if i!=rk and M[i][c]!=0:
                f=M[i][c]/M[rk][c]
                M[i]=[a-f*b for a,b in zip(M[i],M[rk])]
        rk+=1
    return rk

def monos(d,deg):
    out=[]
    for t in range(deg+1):
        for comb in itertools.combinations_with_replacement(range(d),t):
            out.append(comb)
    return out

def rho(S,d):
    # largest j with rank(Poly_{<j}|_S) < |S|
    j=0
    while True:
        # check V(S,j+1) != 0 : Poly_{<j+1} = degree <= j
        ms=monos(d,j)
        rows=[[math.prod(k[i] for i in mono) for k in S] for mono in ms]
        if rank_exact(rows)<len(S): j+=1
        else: return j
        if j>40: return j

def width(S,v):
    vals=[sum(a*b for a,b in zip(k,v)) for k in S]
    return max(vals)-min(vals)

random.seed(26)
worst=0.0; nt=0; fails=0
for trial in range(60):
    n=random.choice([4,5,6,7])
    r=random.randint(n-1,2*n)
    betas=[[random.choice([-1,1]) for _ in range(n)] for _ in range(r)]
    om=[random.random() for _ in range(r)]; s=sum(om); om=[x/s for x in om]
    P=np.eye(n)-np.ones((n,n))/n
    bt=[P@np.array(b,float) for b in betas]
    G=sum(o*np.outer(b,b) for o,b in zip(om,bt))
    ev=np.linalg.eigvalsh(G); ev=sorted(ev)[1:]  # drop the 1-direction (eigenvalue 0)
    lam=ev[0]
    if lam<1e-9: continue
    Gp=np.linalg.pinv(G)
    # random support S in slice sum=0 with |k_x|<=2 for x<n
    box=[k for k in itertools.product(range(-1,2),repeat=n-1)]
    size=random.randint(6,min(40,len(box)))
    Sd=random.sample(box,size)
    S=[tuple(list(k)+[-sum(k)]) for k in Sd]
    rh=rho([k[:-1] for k in S],n-1)
    Om=Fraction(0)
    for eps in itertools.product([-1,1],repeat=n):
        Om+=width(S,eps)
    Om/=2**n
    Lam=sum(o*width(S,b) for o,b in zip(om,betas))
    mid=sum(o*width(S,b)*np.linalg.norm(Gp@bb) for o,b,bb in zip(om,betas,bt))
    cap=math.sqrt(n)/lam*Lam
    ok = rh<=Om+1e-12 and float(Om)<=mid+1e-9 and mid<=cap+1e-9
    nt+=1; fails+= (not ok)
    worst=max(worst, rh/ (cap) if cap>0 else 0)
    if trial<12 or not ok:
        print(f"n={n} r={r} |S|={len(S)} rho={rh} Omega={float(Om):.3f} mid={mid:.3f} cap={cap:.3f} lam={lam:.3f} ok={ok}")
print("tests",nt,"failures",fails,"max rho/cap",round(worst,4))
# Hadamard/Sylvester profile: lam = n/(n-1) exactly
for k in [2,3,4]:
    n=2**k
    pts=list(itertools.product([0,1],repeat=k))
    cols=[[(-1)**(sum(a*b for a,b in zip(av,x))%2) for x in pts] for av in pts if any(av)]
    om=1/(n-1)
    P=np.eye(n)-np.ones((n,n))/n
    G=sum(om*np.outer(P@np.array(c,float),P@np.array(c,float)) for c in cols)
    ev=sorted(np.linalg.eigvalsh(G))[1:]
    print(f"Sylvester n={n}: lam_min on 1-perp = {ev[0]:.6f}, n/(n-1) = {n/(n-1):.6f}")
