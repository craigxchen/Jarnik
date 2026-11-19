# Judge's own check of the D=12 squarefree five-mark witness (no imports from Astra/referees).
# F_9 = F_3[i]/(i^2+1), element (a,b) = a+b*i ; encoding digit = a+3b, little-endian base 9 (Astra's table format)
P=3
def dec(x):
    c=[]
    while x: x,d=divmod(x,9); b,a=divmod(d,3); c.append((a,b))
    return c
def trim(f):
    while f and f[-1]==0: f.pop()
    return f
def pmul(f,g,p=P):
    if not f or not g: return []
    r=[0]*(len(f)+len(g)-1)
    for i,a in enumerate(f):
        if a:
            for j,b in enumerate(g): r[i+j]=(r[i+j]+a*b)%p
    return trim(r)
def padd(f,g,p=P):
    n=max(len(f),len(g)); r=[((f[i] if i<len(f) else 0)+(g[i] if i<len(g) else 0))%p for i in range(n)]
    return trim(r)
def psub(f,g,p=P): return padd(f,[(-x)%p for x in g],p)
def pdivmod(f,g,p=P):
    f=f[:]; q=[0]*max(1,len(f)-len(g)+1); inv=pow(g[-1],p-2,p)
    while len(f)>=len(g) and f:
        c=f[-1]*inv%p; d=len(f)-len(g); q[d]=c
        for j,b in enumerate(g): f[d+j]=(f[d+j]-c*b)%p
        trim(f)
    return trim(q),f
def pgcd(f,g,p=P):
    while g: f,g=g,pdivmod(f,g,p)[1]
    if f: inv=pow(f[-1],p-2,p); f=[x*inv%p for x in f]
    return f
def deriv(f,p=P): return trim([(i*f[i])%p for i in range(1,len(f))])
def ppow_mod(b,e,m,p=P):
    r=[1]; b=pdivmod(b,m,p)[1]
    while e:
        if e&1: r=pdivmod(pmul(r,b,p),m,p)[1]
        b=pdivmod(pmul(b,b,p),m,p)[1]; e>>=1
    return r
def ddf(f,p=P):
    out=[]; d=1; h=[0,1]
    while len(f)>1:
        if 2*d>len(f)-1: out.append((len(f)-1,f)); break
        h=ppow_mod(h,p,f,p); g=pgcd(f,psub(h,[0,1],p),p)
        if len(g)>1: out.append((d,g)); f=pdivmod(f,g,p)[0]; h=pdivmod(h,f,p)[1]
        d+=1
    return out
line="391358742215,391358571511,391360438826,391357359206,391358930855"
Z=[dec(int(s)) for s in line.split(",")]
A=[[c[0] for c in z] for z in Z]; B=[[c[1] for c in z] for z in Z]
for a,b in zip(A,B): trim(a); trim(b)
deg=[len(z)-1 for z in Z]; print("degrees",deg,"monic",[z[-1] for z in Z])
N=[padd(pmul(a,a),pmul(b,b)) for a,b in zip(A,B)]
print("norms equal:",all(n==N[0] for n in N))
target=[0]*25
for e,c in [(24,1),(20,1),(18,2),(16,1),(10,2),(2,2),(0,2)]: target[e]=c
print("N == t^24+t^20+2t^18+t^16+2t^10+2t^2+2:",N[0]==target)
print("gcd(N,N')=",pgcd(N[0],deriv(N[0])))
print("DDF degrees:",[(d,len(g)-1) for d,g in ddf(N[0])])
# distinctness and gaps (difference in F_9[t]: componentwise)
gaps=[]
for i in range(5):
    for j in range(i):
        da=psub(A[i],A[j]); db=psub(B[i],B[j]); gaps.append(max(len(da),len(db))-1)
print("pair gap degrees",sorted(gaps),"max",max(gaps),"floor(D/2)=",12//2, "distinct",min(gaps)>=0)
# conjugates/unit multiples: z_i = u * conj(z_j) or u*z_j with u in F_9^*?
units=[(a,b) for a in range(3) for b in range(3) if (a,b)!=(0,0)]
def umul(u,z): a,b=u; return [((a*x-b*y)%3,(a*y+b*x)%3) for x,y in z]
def conj(z): return [(x,(-y)%3) for x,y in z]
bad=0
for i in range(5):
    for j in range(5):
        for u in units:
            if i!=j and umul(u,Z[j])==Z[i]: bad+=1
            if umul(u,conj(Z[j]))==Z[i]: bad+=1
print("unit/conjugate coincidences:",bad)
# infinity inert: constant extension, deg N even, leading coeff 1 (norm of monic) -> fine
