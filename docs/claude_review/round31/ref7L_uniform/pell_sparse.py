"""Pell five-mark family pulled back by t=u-1/u: bounded-term quotient presentation? height W' and contact h'."""
import sympy as sp
t,u=sp.symbols('t u'); I=sp.I
def blocks(n,seeds):
    A=sp.Matrix([[t**2+1,t],[t,1]]); An=A**n
    out=[]
    for s in seeds:
        v=An*sp.Matrix([1,s]); out.append((sp.expand(v[0]+I*v[1]),sp.expand(v[0]-I*v[1])))
    return out
def nterms(e):
    e=sp.expand(e); return len(sp.Add.make_args(e))
for n in [1,2,3,4]:
    seeds=[1,3,5,8]
    B=blocks(n,seeds)
    sub=u-1/u
    counts=[]
    for G,Gb in B:
        Gu=sp.expand(sp.cancel(G.subs(t,sub)*(u+1/u)))
        counts.append(nterms(Gu))
    beta_u=sp.expand((sub+2*I)*u)  # = (u+i)^2
    # ratio map in u: q_jk = g_j g_k/g_0
    edges=[(0,1),(0,2),(1,2),(0,3)]
    g=[sp.cancel(G.subs(t,sub)/Gb.subs(t,sub)) for G,Gb in B]
    g0=sp.cancel((sub+2*I)/(sub-2*I))
    q=[sp.cancel(g[j]*g[k]/g0) for j,k in edges]
    # intrinsic height W = sum_v (max - min) over places incl 0 and inf; compute via projective degrees
    nums=[];dens=[]
    for x in q:
        a,b=sp.fraction(sp.cancel(x)); nums.append(sp.Poly(a,u,extension=I)); dens.append(sp.Poly(b,u,extension=I))
    L=dens[0]
    for d in dens[1:]: L=sp.lcm(L,d)
    F=[L]+[sp.Poly(sp.cancel(nums[i].as_expr()*L.as_expr()/dens[i].as_expr()),u,extension=I) for i in range(4)]
    g_=F[0]
    for f in F[1:]: g_=sp.gcd(g_,f)
    F=[sp.div(f,g_)[0] for f in F]
    dplus=max(f.degree() for f in F)
    # reciprocal
    Lr=nums[0]
    for d in nums[1:]: Lr=sp.lcm(Lr,d)
    Gl=[Lr]+[sp.Poly(sp.cancel(dens[i].as_expr()*Lr.as_expr()/nums[i].as_expr()),u,extension=I) for i in range(4)]
    g2=Gl[0]
    for f in Gl[1:]: g2=sp.gcd(g2,f)
    Gl=[sp.div(f,g2)[0] for f in Gl]
    dminus=max(f.degree() for f in Gl)
    # contact at u=inf
    h=min(sp.degree(sp.denom(sp.cancel(x-1)),u)-sp.degree(sp.numer(sp.cancel(x-1)),u) for x in q)
    # contact at u=0
    h0=min(sp.Poly(sp.numer(sp.cancel(x-1)),u).monoms()[-1][0] - sp.Poly(sp.denom(sp.cancel(x-1)),u).monoms()[-1][0] for x in q)
    print((n, 'blockterms',counts,'W_u',dplus+dminus,'2W_t',2*(16*n+2),'h_inf',h,'h_0',h0,'4h_inf',4*h))
