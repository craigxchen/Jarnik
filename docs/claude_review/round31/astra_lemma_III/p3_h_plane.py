"""Exact equations on three charts of the projective midpoint-shape plane.

Run p3_h_certificate.py for the published elimination certificates.
"""
from math import lcm
from itertools import combinations
from sympy import Matrix
import flint

edges=list(combinations(range(5),2))
tris=list(combinations(range(5),3))
R=lambda pair:edges.index(tuple(sorted(pair)))
E=Matrix([[int(i in e)-int(4 in e) for e in edges] for i in range(4)])
basis=Matrix.hstack(*[v*lcm(*[int(a.q) for a in v]) for v in E.nullspace()])
assert all(a.is_Integer for a in basis)
def determinant(A):
    if len(A)==1:return A[0][0]
    return sum((1 if j%2==0 else -1)*a*determinant([r[:j]+r[j+1:] for r in A[1:]]) for j,a in enumerate(A[0]))
def minors(m,ij):
    return [[r[j] for j in ij] for r in m]
def calc(p=101,H=(-5,-3,0,2,6),chart=2):
    assert sum(H)==0
    from sympy.polys.matrices import DomainMatrix
    from sympy import GF
    assert DomainMatrix.from_Matrix(basis).convert_to(GF(p)).rank()==6
    assert DomainMatrix.from_Matrix(E).convert_to(GF(p)).rank()==4
    ctx=flint.nmod_mpoly_ctx.get(("u","v"),p)
    u,v=ctx.gens()
    if chart==2:a,b,c=u,v,1
    if chart==1:a,b,c=u,1,v
    if chart==0:a,b,c=1,u,v
    X=[a*(H[i]+H[j])+b*(H[i]*H[i]+H[j]*H[j]-H[i]*H[j])+c*H[i]*H[j]*(H[i]+H[j]) for i,j in edges]
    ll=[[sum(int(E[i,j])*X[j]*int(basis[j,k]) for j in range(10)) for k in range(6)]for i in range(4)]
    dets={ij:determinant(minors(ll,ij)) for ij in combinations(range(6),4)}
    piv,D=next((ij,d) for ij,d in dets.items() if d)
    Y=[]
    for l in sorted(set(range(6))-set(piv)):
        a=[0]*6;a[l]=D
        for j,k in enumerate(piv):
            I=list(piv);I[j]=l
            a[k]=-determinant(minors(ll,I))
        Y.append([sum(int(basis[i,j])*a[j] for j in range(6)) for i in range(10)])
    # Twenty real triangle equations. All are linear in y0^2, y0*y1, y1^2.
    rows=[]
    mix=lambda k,i,j:Y[0][i]*Y[0][j] if k==1 else (Y[1][i]*Y[1][j] if k==3 else Y[0][i]*Y[1][j]+Y[1][i]*Y[0][j])
    for d in (2,3):
        for tri in tris:
            l,m=sorted(set(range(5))-set(tri))
            row=[0]*4
            for k in range(3):
                i,j,n=tri[k],tri[(k+1)%3],tri[(k+2)%3]
                e,f,g=R((i,j)),R((n,l)),R((n,m));h=H[j]-H[i]
                if d==2:
                    row[0]+=h*(X[e]*X[f]+X[e]*X[g]+X[f]*X[g])
                    for z in (1,2,3):row[z]+=h*(mix(z,e,f)+mix(z,e,g)-mix(z,f,g))
                else:
                    row[0]+=h*X[e]*X[f]*X[g]
                    for z in (1,2,3):row[z]+=h*(X[g]*mix(z,e,f)+X[f]*mix(z,e,g)-X[e]*mix(z,f,g))
            rows.append(row)
    eqs=[]
    choices=[(0,1,2), (2,6,9), (3,7,10), (15,17,19),
             (10,1,6), (5,0,18), (11,6,16)]
    extra_rows=iter([(3,17,10,4),(2,8,19,4),(8,18,4,15),
                     (14,13,9,10),(6,17,10,1),(15,9,1,11)])
    Delta=None
    for pick in choices:
        tr=[rows[i] for i in pick]
        det=determinant([r[1:] for r in tr])
        if not det:continue
        if Delta is None:Delta=det
        cr=[]
        for j in (1,2,3):
            cols=[1,2,3];cols[j-1]=0
            cr.append(-determinant(minors(tr,cols)))
        use=range(len(rows)) if not eqs else next(extra_rows)
        for i in use:
            r=rows[i]
            f=det*r[0]+sum(cr[j]*r[j+1] for j in range(3))
            if f:
                assert not f%(D**3)
                eqs.append(f//(D**3))
        vr=cr[1]**2-cr[0]*cr[2]
        # The discriminant of the adjugate binary symmetric tensor has
        # relative weight four under GL_2 basis change.
        assert not vr%(D**4)
        eqs.append(vr//(D**4))
    return ctx, X,dets,eqs,D,Delta

def elim(eqs,axis):
    base=eqs[0]
    X=None
    for f in eqs[1:]:
        v=base.resultant(f,axis)
        if not v:continue
        X=v if X is None else X.gcd(v)
        if X.is_constant():break
    return X
