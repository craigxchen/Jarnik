# Saturate the realization ideal by D_ij=(Y_jV_i-Y_iV_j)^2+(Y_i-Y_j)^2Y_i^2Y_j^2 (= (Y_iY_j)^2|rho_i-rho_j|^2, >0 on
# real configurations with distinct roots).  Unit ideal => no admissible real realization.
import sys, pickle, itertools, sympy as sp
from realize import test, a, b, s, t, z
w=sp.symbols('w')
GENT,out=pickle.load(open(sys.argv[1],'rb'))
for c in out:
    if c['rank']!=5: continue
    r=test(GENT,c['rep']); Y,V=r['Y'],r['V']
    G=list(r['G'].exprs)
    print('class nq',c['nq'],'nmid',c['nm'],'rep',c['rep'][:5]); sys.stdout.flush()
    done=False
    # order pairs: try pairs whose Y's can coincide on elimination factors first (cheap heuristic: all pairs)
    for i,j in itertools.combinations(range(8),2):
        D=sp.expand((Y[j]*V[i]-Y[i]*V[j])**2+(Y[i]-Y[j])**2*Y[i]**2*Y[j]**2)
        if sp.groebner(G,z,s,t,a,b,order='lex').reduce(D)[1]==0 and False: pass
        GB=sp.groebner(G+[w*D-1],w,z,s,t,a,b,order='lex')
        newG=[g for g in GB.exprs if w not in g.free_symbols]
        if newG==[1] or list(GB.exprs)==[1]:
            print('  unit after saturating by pair',(i,j)); done=True; break
        if set(map(sp.expand,newG))!=set(map(sp.expand,G)):
            G=newG
            print('  saturated by',(i,j),'-> (a,b)-part',[sp.factor(g) for g in G if not (g.free_symbols&{z,s,t})]); sys.stdout.flush()
    if not done: print('  NOT unit; remaining GB:',[sp.factor(g) for g in G])
    sys.stdout.flush()
