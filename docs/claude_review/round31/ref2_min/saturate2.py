import sys, pickle, itertools, signal, sympy as sp
from realize import test, a, b, s, t, z
w=sp.symbols('w')
GENT,out=pickle.load(open(sys.argv[1],'rb'))
idx=[int(x) for x in sys.argv[2].split(',')]; tl=int(sys.argv[3])
class TO(Exception): pass
def h(*_): raise TO()
signal.signal(signal.SIGALRM,h)
for n in idx:
    c=out[n]
    if c['rank']!=5: continue
    print('class',n,'nq',c['nq'],'nmid',c['nm'],'rep',c['rep'][:5]); sys.stdout.flush()
    signal.alarm(tl)
    try:
        r=test(GENT,c['rep']); Y,V=r['Y'],r['V']; G=list(r['G'].exprs)
        elim=[g for g in G if not (g.free_symbols&{z,s,t})]
        # heuristic order: pairs with Y_i-Y_j dividing... first pairs whose Y_i-Y_j vanishes on some elimination factor
        facs=set()
        for g in elim:
            for f,_ in sp.factor_list(g)[1]: facs.add(f)
        def score(p):
            i,j=p; d=sp.expand(Y[i]-Y[j])
            return 0 if any(sp.rem(d,f,a)==0 or sp.rem(d,f,b)==0 for f in facs if f.free_symbols) else 1
        pairs=sorted(itertools.combinations(range(8),2),key=score)
        done=False
        for i,j in pairs:
            D=sp.expand((Y[j]*V[i]-Y[i]*V[j])**2+(Y[i]-Y[j])**2*Y[i]**2*Y[j]**2)
            GB=sp.groebner(G+[w*D-1],w,z,s,t,a,b,order='lex')
            newG=[g for g in GB.exprs if w not in g.free_symbols]
            if newG==[1] or list(GB.exprs)==[1]:
                print('  unit after saturating by pair',(i,j)); done=True; break
            if set(map(sp.expand,newG))!=set(map(sp.expand,G)):
                G=newG; print('  saturated by',(i,j),'-> (a,b)-part',[sp.factor(g) for g in G if not (g.free_symbols&{z,s,t})]); sys.stdout.flush()
        if not done: print('  NOT unit; remaining GB:',[sp.factor(g) for g in G])
    except TO:
        print('  TIMEOUT')
    signal.alarm(0); sys.stdout.flush()
