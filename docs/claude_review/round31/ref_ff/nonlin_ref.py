# Independent nonlinear exclusion of the 21 survivors.
import pickle, itertools, sympy as sp, sys
a,b=sp.symbols('a b')
surv=pickle.load(open('survivors.pkl','rb'))
def bits(r): return [(r>>j)&1 for j in range(8)]
def Hmat(S):
    M=sp.Matrix([[x-y for x,y in zip(bits(r),bits(S[0]))] for r in S[1:]])
    return sp.Matrix.hstack(*M.nullspace())
def collinear4(H):
    rows=[tuple(H.row(j)) for j in range(8)]
    for q in itertools.combinations(range(8),4):
        R=[sp.Matrix(rows[j]) for j in q]
        if len(set(rows[j] for j in q))<4: continue
        M=sp.Matrix.hstack(R[1]-R[0],R[2]-R[0],R[3]-R[0])
        if M.rank()==1: return q
    return None
def proportional(H):
    for j,k in itertools.combinations(range(8),2):
        if sp.Matrix.hstack(H.row(j).T,H.row(k).T).rank()<=1: return (j,k)
    return None
res={}
j0=int(sys.argv[1]) if len(sys.argv)>1 else 7   # gauge row (choose a different one from Astra's)
for idx,S in enumerate(surv):
    H=Hmat(S)
    p=proportional(H)
    if p: res[idx]=('prop',p); print(idx,'prop',p,flush=True); continue
    c=collinear4(H)
    if c: res[idx]=('col4',c); print(idx,'col4',c,flush=True); continue
    # change basis so that row g becomes e3 : H -> H G with h_g G = e3
    g=j0 if any(H.row(j0)) else next(j for j in range(8) if any(H.row(j)))
    hg=H.row(g)
    # complete hg to basis: G^{-1} has third row hg
    for e1,e2 in itertools.combinations([sp.Matrix([[1,0,0]]),sp.Matrix([[0,1,0]]),sp.Matrix([[0,0,1]])],2):
        Gi=sp.Matrix.vstack(e1,e2,hg)
        if Gi.det()!=0: break
    G=Gi.inv(); H2=H*G
    assert list(H2.row(g))==[0,0,1]
    u=sp.Matrix([a,b,1])
    y=[ (H2.row(j)*u)[0] for j in range(8)]
    mon=lambda h:[h[0]**2,h[1]**2,h[2]**2,2*h[0]*h[1],2*h[0]*h[2],2*h[1]*h[2]]
    A=sp.Matrix([mon(list(H2.row(j))) for j in range(8)])
    assert A.rank()==6, (idx,A.rank())
    lk=A.T.nullspace()          # z with z^T A = 0
    comp=[sp.expand(sum(z[j]*y[j]**4 for j in range(8))) for z in lk]
    # solve Q from all equations (least squares exact = any 6 independent rows)
    sel=None
    for T in itertools.combinations(range(8),6):
        if A.extract(list(T),list(range(6))).rank()==6: sel=T;break
    q=A.extract(list(sel),list(range(6))).LUsolve(sp.Matrix([y[j]**4 for j in sel]))
    Q=sp.Matrix([[q[0],q[3],q[4]],[q[3],q[1],q[5]],[q[4],q[5],q[2]]])
    F=sp.Matrix([[1,0,-a],[0,1,-b]])
    det=sp.expand(sp.together((F*Q*F.T).det()))
    det=sp.numer(sp.together(det))
    comp=[sp.numer(sp.together(c)) for c in comp]
    Gb=sp.groebner(comp+[det],a,b,order='grevlex',domain=sp.QQ)
    if list(Gb.exprs)==[1]:
        res[idx]=('unit',); print(idx,'unit',flush=True); continue
    # otherwise: positivity: product of y_j and |rho_j-rho_k|^2 y_j^2 y_k^2
    prod=sp.Integer(1)
    for j in range(8): prod=Gb.reduce(sp.expand(prod*y[j]))[1]
    for j,k in itertools.combinations(range(8),2):
        z=(y[k]*H2.row(j)-y[j]*H2.row(k)).T
        Bjk=(y[j]-y[k])**2*y[j]**2*y[k]**2+sp.Rational(1,3)*(z.T*Q*z)[0]
        Bjk=sp.numer(sp.together(sp.expand(Bjk)))  # Q entries are polynomials (denominators constants)
        prod=Gb.reduce(sp.expand(prod*Bjk))[1]
        if prod==0: break
    res[idx]=('pos0' if prod==0 else 'OPEN',Gb.exprs); print(idx,res[idx][0],Gb.exprs if prod!=0 else '',flush=True)
pickle.dump(res,open(f'nonlin_res_{j0}.pkl','wb'))
