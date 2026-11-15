import math, cmath, itertools
from sympy import factorint
def pts(N):
    r=math.isqrt(N); out=[]
    for a in range(-r,r+1):
        b2=N-a*a; b=math.isqrt(b2)
        if b*b==b2: out+= [(a,b)] if b==0 else [(a,b),(a,-b)]
    return out
# (1) lambda_k product formula; include cases with 2 and p=3 mod 4
def gauss_prime(p):
    for a in range(1,math.isqrt(p)+1):
        b=math.isqrt(p-a*a)
        if a*a+b*b==p: return complex(a,b)
for N in [5**2*13*17, 5*13**3*29, 2**3*5**2*13, 3**2*5*17, 2*9*5*13]:
    P=pts(N); f=factorint(N)
    for k in [1,2,3,7]:
        lam=sum((complex(a,b)/abs(complex(a,b)))**(4*k) for a,b in P)
        prod=4
        for p,e in f.items():
            if p%4==1:
                ph=cmath.phase(gauss_prime(p)); x=4*k*ph
                prod*= (e+1) if abs(math.sin(x))<1e-12 else math.sin((e+1)*x)/math.sin(x)
            elif p==2: prod*=(-1)**(k*e)
            else: prod*= (1 if e%2==0 else 0)
        print(f"N={N} k={k} lambda={lam.real:+.6f}{lam.imag:+.1e}i product(with 2,3mod4 factors)={prod:+.6f}")
# (2) Prop 2.6: a^2=a'^2 mod N forced equal for a,a' in an interval of length l=C sqrt R inside (0,R]
for N,C in [(10**8+1,4),(5**4*13**2*17**2,4),(2576450045,8)]:
    R=math.sqrt(N); l=int(C*math.sqrt(R)); a0=int(R*0.6)
    J=[a for a in range(a0,a0+l) if math.gcd(a,N)==1]
    sq={}
    coll=0
    for a in J:
        s=a*a%N; coll+= sq.get(s,0); sq[s]=sq.get(s,0)+1
    print(f"N={N} R>16C^2:{R>16*C*C} |J|={len(J)} nontrivial collisions a^2=a'^2 mod N: {coll}")
