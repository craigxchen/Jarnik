# Referee realizability test for planar (rank-5) flats on the moment surface 3V^2 = Y W + Y^4, Y>0.
# Actual (Y,V,W) columns lie in ker E = col span(H); normalise H so three independent rows are e1,e2,e3,
# then u=(Y_s1,Y_s2,Y_s3) has positive entries; dilate to Y_s3=1; translate x so x_s3=0 (V_s3=0).
# W is eliminated: W_i=(3V_i^2-Y_i^4)/Y_i must be orthogonal to E.
import sys, pickle, itertools, sympy as sp
a,b,s,t,z = sp.symbols('a b s t z')
def gvec(T):
    v=[0]*8
    if T[0]=='q': i,j,k,l=T[1]; v[i]+=1; v[l]+=1; v[j]-=1; v[k]-=1
    else: i,j,k=T[1]; v[i]+=1; v[k]+=1; v[j]-=2
    return v
def setup(GENT, ids):
    E=sp.Matrix([gvec(GENT[i]) for i in ids]); r=E.rank(); assert r==5
    N=sp.Matrix.hstack(*E.nullspace())
    for S in itertools.combinations(range(8),3):
        NS=N.extract(list(S),[0,1,2])
        if NS.det()!=0: break
    H=(N*NS.inv()).applyfunc(sp.nsimplify)
    Eb=sp.Matrix.vstack(*[sp.Matrix([list(row)]) for row in E.rowspace()])
    return E,H,S,Eb
def test(GENT, ids, verbose=True):
    E,H,S,Eb=setup(GENT,ids)
    u=sp.Matrix([a,b,1]); v=sp.Matrix([s,t,0])
    Y=[sp.expand((H[i,:]*u)[0]) for i in range(8)]; V=[sp.expand((H[i,:]*v)[0]) for i in range(8)]
    eqs=[]
    for k in range(Eb.rows):
        m=Eb[k,:]
        L=sp.Mul(*[Y[i] for i in range(8) if m[i]!=0])
        e=sp.expand(sum(m[i]*sp.cancel((3*V[i]**2-Y[i]**4)*L/Y[i]) for i in range(8) if m[i]!=0))
        eqs.append(e)
    Yall=sp.Mul(*[Y[i] for i in range(8)])
    G=sp.groebner(eqs+[z*Yall-1], z,s,t,a,b, order='lex')
    res=dict(S=S,H=H,Y=Y,V=V,eqs=eqs,G=G)
    return res
if __name__=='__main__':
    GENT,out=pickle.load(open(sys.argv[1],'rb'))
    for c in out:
        if c['rank']!=5: continue
        r=test(GENT,c['rep'])
        print('class nq',c['nq'],'nmid',c['nm'],'rep',c['rep'][:5],'S',r['S'])
        print('  H rows',[tuple(r['H'][i,:]) for i in range(8)])
        print('  Y',r['Y'])
        if list(r['G'].exprs)==[1]: print('  UNIT IDEAL: no complex solution with all Y_i != 0')
        else:
            print('  GB (a,b)-elimination part:',[sp.factor(g) for g in r['G'].exprs if not (g.free_symbols & {z,s,t})])
            print('  GB size',len(r['G'].exprs))
        sys.stdout.flush()
