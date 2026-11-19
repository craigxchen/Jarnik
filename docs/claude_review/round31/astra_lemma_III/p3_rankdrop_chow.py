"""Boundary-format determinant detecting possible rank loss of the 4-by-6 matrix."""
from itertools import combinations
from sympy import Matrix
from math import lcm
from flint import nmod_mat,nmod_poly
edges=list(combinations(range(5),2))
E=Matrix([[int(i in e)-int(4 in e) for e in edges] for i in range(4)])
Z=Matrix.hstack(*[v*lcm(*[int(c.q) for c in v]) for v in E.nullspace()])
W=[[[int(E[i,e]*Z[e,j]) for e in range(10)] for j in range(6)] for i in range(4)]
def mon(d):return [(a,b,d-a-b) for a in range(d+1) for b in range(d+1-a)]
three,four=mon(3),mon(4)
positions={(i,j,k):(i*15+four.index(tuple(z+int(m==k) for m,z in enumerate(mm))), j*10+l) for i in [] for j in [] for l,mm in enumerate(three) for k in range(3)}
def value(h,p):
    h=[(5*z-sum(h))%p for z in h]
    X=[[(h[i]+h[j])%p,(h[i]**2+h[j]**2-h[i]*h[j])%p, (h[i]*h[j]*(h[i]+h[j]))%p] for i,j in edges]
    M=[[0]*60 for _ in range(60)]
    for i in range(4):
        for j in range(6):
            L=[sum(c*x[k] for c,x in zip(W[i][j],X))%p for k in range(3)]
            for l,mm in enumerate(three):
                for k in range(3):
                    row=i*15+four.index(tuple(z+int(a==k) for a,z in enumerate(mm)))
                    M[row][j*10+l]=L[k]
    return int(nmod_mat(M,p).det())
def interpolate(values,p):
    base=nmod_poly([1],p);f=nmod_poly([],p)
    for a,v in enumerate(values):
        k=((v-int(f(a)))*pow(int(base(a)),-1,p))%p
        f+=base*k
        base*=nmod_poly([-a,1],p)
    return f
if __name__=='__main__':
    h=[-5,-3,0,2];p=1009
    v=interpolate([value(h+[a],p) for a in range(121)],p)
    assert int(v(582))==value(h+[582],p)
    print('specialized degree',v.degree(),'factors',[(str(f) if f.degree()<=1 else f.degree(),e) for f,e in v.factor()[1]])
