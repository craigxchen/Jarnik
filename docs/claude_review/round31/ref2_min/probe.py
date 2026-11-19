import sys, pickle, itertools, sympy as sp, time
from realize import setup, a, b, s, t, z
exec(open('solve2d.py').read().split('def satur')[0].split("GENT,out")[0])
GENT,out=pickle.load(open('c_mid12.out.classes.pkl','rb'))
def lin(e):
    p=sp.Poly(sp.expand(e),s,t); r=[0,0,0,0]
    for (i,j),c in p.terms(): r[{(2,0):0,(1,1):1,(0,2):2,(0,0):3}[(i,j)]]+=c
    return r
def build(n):
    c=out[n]; E,H,S,Eb=setup(GENT,c['rep'])
    u=sp.Matrix([a,b,1]); v=sp.Matrix([s,t,0])
    Y=[sp.expand((H[i,:]*u)[0]) for i in range(8)]; V=[sp.expand((H[i,:]*v)[0]) for i in range(8)]
    rows=[]
    for k in range(Eb.rows):
        m=Eb[k,:]; L=sp.Mul(*[Y[i] for i in range(8) if m[i]!=0])
        rows.append(lin(sum(m[i]*sp.cancel((3*V[i]**2-Y[i]**4)*L/Y[i]) for i in range(8) if m[i]!=0)))
    return H,Y,V,rows
if __name__=='__main__':
    n=int(sys.argv[1]); H,Y,V,rows=build(n)
    print('H',[tuple(H[i,:]) for i in range(8)]); print('Y',Y)
    M=sp.Matrix([r[:3] for r in rows]); cc=sp.Matrix([r[3] for r in rows])
    print('deg of entries',[[sp.Poly(M[k,q],a,b).total_degree() if M[k,q]!=0 else -1 for q in range(3)]+[sp.Poly(cc[k],a,b).total_degree()] for k in range(5)])
    for R in itertools.combinations(range(5),3):
        MR=M.extract(list(R),[0,1,2]); dR=sp.expand(MR.det())
        if dR!=0: break
    t0=time.time()
    adj=MR.adjugate(); Xn=[sp.expand(x) for x in -(adj*cc.extract(list(R),[0]))]
    conds=[sp.expand(sum(M[k,q]*Xn[q] for q in range(3))+cc[k]*dR) for k in range(5) if k not in R]+[sp.expand(Xn[0]*Xn[2]-Xn[1]**2)]
    print('R',R,'cond degrees',[sp.Poly(c,a,b).total_degree() for c in conds], time.time()-t0)
    g=sp.gcd_list(conds); print('gcd factored:',sp.factor(g), time.time()-t0)
    print('dR factored:',sp.factor(dR))
