# Same certificate as saturate.py but with X=(s^2,st,t^2) as unknowns (the equations and every D_ij are
# linear in X), plus X1*X3 = X2^2.  Unit ideal => no real admissible realization.
import sys, pickle, itertools, signal, sympy as sp
from realize import setup, a, b, s, t, z
w,X1,X2,X3=sp.symbols('w X1 X2 X3')
GENT,out=pickle.load(open(sys.argv[1],'rb'))
idx=[int(x) for x in sys.argv[2].split(',')]; tl=int(sys.argv[3])
class TO(Exception): pass
def h(*_): raise TO()
signal.signal(signal.SIGALRM,h)
def lin(e):
    p=sp.Poly(sp.expand(e),s,t); r=0
    for (i,j),c in p.terms():
        r+=c*{(2,0):X1,(1,1):X2,(0,2):X3,(0,0):1}[(i,j)]
    return sp.expand(r)
VARS=(w,z,X1,X2,X3,a,b)
for n in idx:
    c=out[n]
    if c['rank']!=5: continue
    print('class',n,'nq',c['nq'],'nmid',c['nm'],'rep',c['rep'][:5]); sys.stdout.flush()
    signal.alarm(tl)
    try:
        E,H,S,Eb=setup(GENT,c['rep'])
        u=sp.Matrix([a,b,1]); v=sp.Matrix([s,t,0])
        Y=[sp.expand((H[i,:]*u)[0]) for i in range(8)]; V=[sp.expand((H[i,:]*v)[0]) for i in range(8)]
        eqs=[]
        for k in range(Eb.rows):
            m=Eb[k,:]; L=sp.Mul(*[Y[i] for i in range(8) if m[i]!=0])
            eqs.append(lin(sum(m[i]*sp.cancel((3*V[i]**2-Y[i]**4)*L/Y[i]) for i in range(8) if m[i]!=0)))
        G=list(sp.groebner(eqs+[X1*X3-X2**2, z*sp.Mul(*Y)-1],*VARS[1:],order='lex').exprs)
        done = G==[1]
        if done: print('  unit before pair saturation')
        for i,j in itertools.combinations(range(8),2):
            if done: break
            D=lin((Y[j]*V[i]-Y[i]*V[j])**2+(Y[i]-Y[j])**2*Y[i]**2*Y[j]**2)
            GB=sp.groebner(G+[w*D-1],*VARS,order='lex')
            newG=[g for g in GB.exprs if w not in g.free_symbols]
            if newG==[1] or list(GB.exprs)==[1]:
                print('  unit after saturating by pair',(i,j)); done=True; break
            if set(map(sp.expand,newG))!=set(map(sp.expand,G)):
                G=newG; print('  saturated by',(i,j),'-> (a,b)-part',[sp.factor(g) for g in G if not (g.free_symbols&{z,X1,X2,X3})]); sys.stdout.flush()
        if not done: print('  NOT unit; remaining GB:',[sp.factor(g) for g in G])
    except TO:
        print('  TIMEOUT')
    signal.alarm(0); sys.stdout.flush()
