# Two-variable certificate.  Equations are linear in X=(s^2,st,t^2): M(a,b) X + c(a,b) = 0 (5 rows).
# For every 3-row subset R with det M_R != 0 identically, on the open set det M_R != 0 we solve X by Cramer
# and get polynomial conditions in (a,b) only; we saturate by det_R * prod Y_i and by each D_ij numerator.
# Unit ideal for every R => no real realization with rank M = 3.  The locus rank M <= 2 is treated by
# the full 6-variable ideal plus all 3x3 minors of M (small).
import sys, pickle, itertools, signal, sympy as sp
from realize import setup, a, b, s, t, z
w,X1,X2,X3=sp.symbols('w X1 X2 X3')
GENT,out=pickle.load(open(sys.argv[1],'rb'))
idx=[int(x) for x in sys.argv[2].split(',')]
def lin(e):
    p=sp.Poly(sp.expand(e),s,t); r=[0,0,0,0]
    for (i,j),c in p.terms(): r[{(2,0):0,(1,1):1,(0,2):2,(0,0):3}[(i,j)]]+=c
    return r
def satur(G, polys, vars_):
    for D in polys:
        if D==0: continue
        GB=sp.groebner(G+[w*D-1],w,*vars_,order='lex')
        G=[g for g in GB.exprs if w not in g.free_symbols]
        if G==[1]: return G
    return G
for n in idx:
    c=out[n]; print('class',n,'nq',c['nq'],'nmid',c['nm'],'rep',c['rep'][:5]); sys.stdout.flush()
    E,H,S,Eb=setup(GENT,c['rep'])
    u=sp.Matrix([a,b,1]); v=sp.Matrix([s,t,0])
    Y=[sp.expand((H[i,:]*u)[0]) for i in range(8)]; V=[sp.expand((H[i,:]*v)[0]) for i in range(8)]
    rows=[]
    for k in range(Eb.rows):
        m=Eb[k,:]; L=sp.Mul(*[Y[i] for i in range(8) if m[i]!=0])
        rows.append(lin(sum(m[i]*sp.cancel((3*V[i]**2-Y[i]**4)*L/Y[i]) for i in range(8) if m[i]!=0)))
    M=sp.Matrix([r[:3] for r in rows]); cc=sp.Matrix([r[3] for r in rows])
    Dl=[lin((Y[j]*V[i]-Y[i]*V[j])**2+(Y[i]-Y[j])**2*Y[i]**2*Y[j]**2) for i,j in itertools.combinations(range(8),2)]
    PY=sp.Mul(*Y); allunit=True
    for R in itertools.combinations(range(5),3):
        MR=M.extract(list(R),[0,1,2]); dR=sp.expand(MR.det())
        if dR==0: continue
        adj=MR.adjugate(); Xn=-(adj*cc.extract(list(R),[0]))   # X = Xn/dR
        Xn=[sp.expand(x) for x in Xn]
        conds=[sp.expand(sum(M[k,q]*Xn[q] for q in range(3))+cc[k]*dR) for k in range(5) if k not in R]
        conds.append(sp.expand(Xn[0]*Xn[2]-Xn[1]**2))
        bad=sp.factor_list(sp.expand(dR*PY))[1]; badf=[f for f,_ in bad if f.free_symbols]
        while True:
            g=sp.gcd_list(conds); hit=False
            for f in badf:
                if sp.rem(g,f,a)==0 if f.has(a) else sp.rem(g,f,b)==0:
                    conds=[sp.cancel(cnd/f) for cnd in conds]; hit=True
            if not hit: break
        g=sp.gcd_list(conds)
        if sp.Poly(g,a,b).total_degree()>0: print('   residual common factor (curve component):',sp.factor(g))
        G=list(sp.groebner(conds+[z*dR*PY-1],z,a,b,order='lex').exprs)
        Ds=[sp.expand(sum(d[q]*Xn[q] for q in range(3))+d[3]*dR) for d in Dl]
        if G!=[1]: G=satur(G,Ds,(z,a,b))
        print('  R',R,'unit' if G==[1] else ('NOT unit: '+str([sp.factor(g) for g in G]))); sys.stdout.flush()
        allunit &= (G==[1])
    # rank-drop locus
    minors=[sp.expand(M.extract(list(R),[0,1,2]).det()) for R in itertools.combinations(range(5),3)]
    G=list(sp.groebner([m for m in minors if m!=0]+[z*PY-1],z,a,b,order='lex').exprs)
    if G==[1]: print('  rank<=2 locus empty (Y!=0)')
    else:
        print('  rank<=2 locus:',[sp.factor(g) for g in G if z not in g.free_symbols])
        eqs=[sum(r[q]*x for q,x in enumerate((s**2,s*t,t**2)))+r[3] for r in rows]
        G2=list(sp.groebner(eqs+[g for g in G if z not in g.free_symbols]+[z*PY-1],z,s,t,a,b,order='lex').exprs)
        Dfull=[(Y[j]*V[i]-Y[i]*V[j])**2+(Y[i]-Y[j])**2*Y[i]**2*Y[j]**2 for i,j in itertools.combinations(range(8),2)]
        G2=satur(G2,[sp.expand(d) for d in Dfull],(z,s,t,a,b))
        print('  rank<=2 full system after saturation:', 'unit' if G2==[1] else [sp.factor(g) for g in G2])
        allunit &= (G2==[1])
    print('  ALL UNIT' if allunit else '  needs manual look'); sys.stdout.flush()
