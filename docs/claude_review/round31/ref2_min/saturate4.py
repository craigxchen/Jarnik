# As saturate2.py, but the z-saturation also removes Y_i = Y_j for every midpoint triple (j midpoint of i,k):
# equal heights on a midpoint triple put three collinear points at one height, impossible on the
# fixed-height parabola 3V^2 = Y W + Y^4 (checked in the report).
import sys, pickle, itertools, signal, sympy as sp
from realize import setup, a, b, s, t, z
import realize
w=sp.symbols('w')
GENT,out=pickle.load(open(sys.argv[1],'rb'))
idx=[int(x) for x in sys.argv[2].split(',')]; tl=int(sys.argv[3])
class TO(Exception): pass
def h(*_): raise TO()
signal.signal(signal.SIGALRM,h)
for n in idx:
    c=out[n]; print('class',n,'nq',c['nq'],'nmid',c['nm'],'rep',c['rep'][:5]); sys.stdout.flush()
    signal.alarm(tl)
    try:
        E,H,S,Eb=setup(GENT,c['rep'])
        u=sp.Matrix([a,b,1]); v=sp.Matrix([s,t,0])
        Y=[sp.expand((H[i,:]*u)[0]) for i in range(8)]; V=[sp.expand((H[i,:]*v)[0]) for i in range(8)]
        # check the midpoint relations hold for formal rows, then collect height differences
        hd=[]
        for r in c['rel']:
            if r[0]=='m':
                j,ik=r[1]; i,k=sorted(ik); assert H[i,:]+H[k,:]==2*H[j,:]; hd.append(sp.expand(Y[i]-Y[j]))
        eqs=[]
        for k in range(Eb.rows):
            m=Eb[k,:]; L=sp.Mul(*[Y[i] for i in range(8) if m[i]!=0])
            eqs.append(sp.expand(sum(m[i]*sp.cancel((3*V[i]**2-Y[i]**4)*L/Y[i]) for i in range(8) if m[i]!=0)))
        G=list(sp.groebner(eqs+[z*sp.Mul(*Y)*sp.Mul(*hd)-1],z,s,t,a,b,order='lex').exprs)
        done = G==[1]
        if done: print('  unit before pair saturation')
        for i,j in itertools.combinations(range(8),2):
            if done: break
            D=sp.expand((Y[j]*V[i]-Y[i]*V[j])**2+(Y[i]-Y[j])**2*Y[i]**2*Y[j]**2)
            GB=sp.groebner(G+[w*D-1],w,z,s,t,a,b,order='lex')
            newG=[g for g in GB.exprs if w not in g.free_symbols]
            if newG==[1] or list(GB.exprs)==[1]: print('  unit after saturating by pair',(i,j)); done=True; break
            if set(map(sp.expand,newG))!=set(map(sp.expand,G)):
                G=newG; print('  saturated by',(i,j),'-> (a,b)-part',[sp.factor(g) for g in G if not (g.free_symbols&{z,s,t})]); sys.stdout.flush()
        if not done: print('  NOT unit; remaining GB:',[sp.factor(g) for g in G])
    except TO: print('  TIMEOUT')
    signal.alarm(0); sys.stdout.flush()
