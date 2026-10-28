# convert a verified x-structured solution to the s-form places rho_T (T subset [4]) and test the general system (3.1)
import numpy as np, ast, itertools, mpmath as mp, sys
import xstruct as X
from verify_x import refine, poly_from_roots
mp.mp.dps=60
def places_of(y,sg):
    rho=[mp.mpc(y[k],y[7+k]) for k in range(7)]; r4=mp.mpc(y[14],y[15]); sig=sg*mp.e**y[16]
    beta=[mp.mpf(1),y[17],y[18],y[19]]; q=[y[20],y[21],y[22]]
    n=poly_from_roots(rho+[r4]); In=[mp.im(x) for x in n]
    qt=[-v/sig for v in q]; ar=[v/(2*sig) for v in In]
    Q2,Q1,Q0=qt; r2=mp.sqrt(Q2); w0=Q1/(2*r2)
    def wroots_over(xv): return mp.polyroots([1,2*r2*xv,-(Q0+xv*Q1)])
    def Hval(xv,zv,sgn):
        A=sum(ar[k]*xv**(len(ar)-1-k) for k in range(len(ar))); Bv=sum(beta[k]*xv**(3-k) for k in range(4))
        return A+sgn*mp.mpc(0,1)+Bv*zv
    pl={}
    def add(w,j,o):
        for kk in list(pl):
            if abs(kk-w)<mp.mpf(10)**(-20): pl[kk][j]+=o; return
        pl[w]=[0]*5; pl[w][j]+=o
    for i,S in enumerate(X.S3):
        for w in wroots_over(rho[i]):
            for j in S: add(w,j-1,1)
        for w in wroots_over(mp.conj(rho[i])):
            for j in S: add(w,j-1,-1)
    for r in rho+[r4]:
        for w in wroots_over(r):
            xv=(Q0-w**2)/(2*r2*w-Q1)
            if abs(Hval(xv,r2*xv+w,1))<mp.mpf(10)**(-20): add(w,3,1)
        for w in wroots_over(mp.conj(r)):
            xv=(Q0-w**2)/(2*r2*w-Q1)
            if abs(Hval(xv,r2*xv+w,-1))<mp.mpf(10)**(-20): add(w,3,-1)
    rhoT={}
    for w,v in pl.items():
        v=np.array(v); v=v-v.min()
        if not any(v): continue
        A=frozenset(i+1 for i in range(5) if v[i]==1)
        if 5 not in A: rhoT[A]=w-w0
    return rhoT,w0
SUBS=[frozenset(c) for k in range(1,5) for c in itertools.combinations(range(1,5),k)]
def eqs(rT):
    out=[]
    for j in range(1,5):
        c=poly_from_roots([rT[T] for T in SUBS if j in T])
        out+=[mp.im(c[k]) for k in (1,2,3,5,6,7,8)]
    return out
if __name__=='__main__':
    lines=[l for l in open(sys.argv[1]) if 'NONDEG' in l]
    for l in lines:
        y=np.array(ast.literal_eval(l[l.index('['):].strip().replace('np.float64','')))
        sg=None
        for s_ in (1.0,-1.0):
            X.SIGN=s_
            if np.linalg.norm(X.F(y))<1e-9: sg=s_
        y,nf=refine(list(y),sg)
        if nf>1e-40: print(l.split()[0],'not refined; skip'); continue
        rT,w0=places_of(y,sg)
        print('== solution',l.split()[0],': places found',len(rT))
        if len(rT)!=15: continue
        e=eqs(rT); print('   max |eq (3.1)| =',mp.nstr(max(abs(v) for v in e),5))
        # Jacobian rank of the 28 equations w.r.t. 30 real unknowns
        h=mp.mpf(10)**(-25); cols=[]
        for T in SUBS:
            for d in (1,1j):
                rp=dict(rT); rp[T]=rp[T]+h*d; rm=dict(rT); rm[T]=rm[T]-h*d
                cols.append([(a-b)/(2*h) for a,b in zip(eqs(rp),eqs(rm))])
        J=np.array([[float(c[i]) for c in cols] for i in range(28)])
        sv=np.linalg.svd(J,compute_uv=False)
        print('   singular values (smallest 4):',sv[-4:],' rank:',int(np.sum(sv>1e-10*sv[0])),' -> local dim incl. gauge',30-int(np.sum(sv>1e-10*sv[0])))
        print('   rho_T (s-form):')
        for T in SUBS: print('     ',''.join(map(str,sorted(T))),mp.nstr(rT[T],20))
