"""Exact checks of Prop 2.2 (Hilbert 90 binary-form model) on random primitive 5-tuples.
Gaussian integers as Python complex-free pairs (a,b)."""
import random, math
from fractions import Fraction
def mul(z,w): return (z[0]*w[0]-z[1]*w[1], z[0]*w[1]+z[1]*w[0])
def conj(z): return (z[0],-z[1])
def N(z): return z[0]**2+z[1]**2
def divmod_g(z,w):
    n=N(w); t=mul(z,conj(w))
    q=((2*t[0]+n)//(2*n),(2*t[1]+n)//(2*n))
    r=(z[0]-mul(q,w)[0], z[1]-mul(q,w)[1]); return q,r
def ggcd(z,w):
    while w!=(0,0):
        _,r=divmod_g(z,w); z,w=w,r
    return z
def exact_div(z,w):
    q,r=divmod_g(z,w); assert r==(0,0); return q
# Gaussian primes above p = 1 mod 4
def pi_of(p):
    for a in range(1,p):
        b2=p-a*a
        b=int(math.isqrt(b2)) if b2>0 else 0
        if b*b==b2 and b>0: return (a,b)
PR=[5,13,17,29,37,41,53,61]
PI={p:pi_of(p) for p in PR}
def vpi(z,pi):
    v=0
    while True:
        q,r=divmod_g(z,pi)
        if r!=(0,0) or z==(0,0): return v
        z=q; v+=1
def randP():
    z=(1,0)
    if random.random()<0.3: z=mul(z,(1,1))
    for p in PR:
        k=random.choice([0,0,1,2])
        g=PI[p] if random.random()<0.5 else conj(PI[p])
        for _ in range(k): z=mul(z,g)
    u=random.choice([(1,0),(0,1),(-1,0),(0,-1)])
    return mul(z,u)
def primitive(z): return math.gcd(z[0],z[1])==1
random.seed(7); ok=0; tests=0; t_int_same_class=0; t_sqrt_half=0
for trial in range(300):
    Ps=[]
    while len(Ps)<5:
        z=randP()
        if primitive(z) and z not in Ps and conj(z) not in Ps and (-z[0],-z[1]) not in Ps: Ps.append(z)
    L=(1,0)
    for P in Ps: L=mul(L,conj(P))
    zeta=[exact_div(mul(P,L),conj(P)) for P in Ps]   # = P_j * prod_{k!=j} conj(P_k)
    Nz=[N(z) for z in zeta]; assert len(set(Nz))==1
    g=zeta[0]
    for z in zeta[1:]: g=ggcd(g,z)
    red=[exact_div(z,g) for z in zeta]; R2=N(red[0])
    eps={p:[vpi(P,PI[p])-vpi(P,conj(PI[p])) for P in Ps] for p in PR}
    H=1
    for p in PR: H*=p**(max(eps[p])-min(eps[p]))
    # (iii) height
    assert H==R2, (H,R2)
    for j,P in enumerate(Ps):
        e=1 if (P[0]-P[1])%2==0 else 0
        assert N(P)==2**e*math.prod(p**abs(eps[p][j]) for p in PR)
    # (ii),(iv)
    for a in range(5):
        for b in range(a+1,5):
            Pa,Pb=Ps[a],Ps[b]; br=Pa[0]*Pb[1]-Pa[1]*Pb[0]
            c=math.prod(p**max(0,min(eps[p][a],eps[p][b]),min(-eps[p][a],-eps[p][b])) for p in PR)
            assert br%c==0
            # |x_a/x_b - 1|^2 = 4 br^2/(N Pa N Pb)  (exact rational check via zeta)
            za,zb=red[a],red[b]; d2=N((za[0]-zb[0],za[1]-zb[1]))
            assert Fraction(d2,R2)==Fraction(4*br*br,N(Pa)*N(Pb))
            # t_ab: |zeta_a-zeta_b| = 2|t| R'' e^{-d/2}, e^{-d}=prod p^{-|eps_a-eps_b|}
            ed=math.prod(p**abs(eps[p][a]-eps[p][b]) for p in PR)
            t2=Fraction(d2*ed,4*R2)   # = t_ab^2
            ea=1 if (Pa[0]-Pa[1])%2==0 else 0; eb=1 if (Pb[0]-Pb[1])%2==0 else 0
            assert t2==Fraction(br*br, c*c*2**(ea+eb))
            if ea==eb: assert t2.denominator==1; t_int_same_class+=1
            else: t_sqrt_half+=1
            tests+=1
    ok+=1
print(f"trials {ok}, pairs {tests}: H_B=R''^2, N(P_j) formula, [P_a,P_b] divisible by prod p^c_ab, dist identity,")
print(f"t_ab^2 = [P_a,P_b]^2/(prod p^(2c_ab) 2^(e_a+e_b)) all exact. pairs with e_a=e_b: {t_int_same_class} (t_ab integer), e_a!=e_b: {t_sqrt_half} (t_ab^2 in Z/2)")
