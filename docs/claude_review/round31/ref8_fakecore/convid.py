# Exact check (in Q(sqrt(-q))) of the identity f * (f o n) = (sum f)^2/q for the finite part f of an
# eigenvector y=(C+g)v of the Paley skew conference matrix, n a nonsquare; and of the SL_2 signed-permutation
# symmetry commuting with C.
import sympy as sp, random
rng=random.Random(2)
for q in [11,19,23,43]:
    g=sp.sqrt(-q); sq={i*i%q for i in range(1,q)}; chi=lambda x: 0 if x%q==0 else (1 if x%q in sq else -1)
    n=q+1; C=sp.zeros(n,n)
    for i in range(q):
        for j in range(q): C[i,j]=chi(j-i)
        C[i,q]=1; C[q,i]=-1
    assert C*C==-q*sp.eye(n)
    v=sp.Matrix([rng.randint(-3,3) for _ in range(n)]); y=(C+g*sp.eye(n))*v
    assert sp.simplify(C*y-g*y)==sp.zeros(n,1)
    f=[sp.expand(y[i]) for i in range(q)]; nn=next(x for x in range(2,q) if x not in sq)
    conv=[sp.expand(sum(f[t]*f[(nn*(x-t))%q] for t in range(q))) for x in range(q)]
    c=sp.expand(sum(f)**2/q)
    ok=all(sp.simplify(cv-c)==0 for cv in conv)
    # SL2 symmetry: gamma=[[a,b],[c,d]] acting on s_i=(1,i), s_inf=(0,1) (row vectors, gamma*s as column)
    def act(gm):
        a,b,c_,d=gm; perm={}; sign={}
        for idx in list(range(q))+['inf']:
            s=(1,idx) if idx!='inf' else (0,1)
            t=((a*s[0]+b*s[1])%q,(c_*s[0]+d*s[1])%q)
            if t[0]!=0: lam=t[0]; img=(t[1]*pow(t[0],-1,q))%q
            else: lam=t[1]; img='inf'
            perm[idx]=img; sign[idx]=chi(lam)
        Pm=sp.zeros(n,n); pos=lambda z: q if z=='inf' else z
        for idx in perm: Pm[pos(perm[idx]),pos(idx)]=sign[idx]
        return Pm
    comm=True
    for gm in [(1,1,0,1),(2,0,0,pow(2,-1,q)),(0,1,q-1,0),(3,5,1,2*pow(1,1,q))]:
        a,b,c_,d=gm
        if (a*d-b*c_)%q!=1: d=((1+b*c_)*pow(a,-1,q))%q if a%q else d; 
        if (a*d-b*c_)%q!=1: continue
        Pm=act((a,b,c_,d)); comm&= (Pm*C==C*Pm)
    print("q=%d convolution identity: %s ; SL2 signed permutations commute with C: %s"%(q,ok,comm))
