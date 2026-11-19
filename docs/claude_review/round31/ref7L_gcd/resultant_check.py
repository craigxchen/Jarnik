# Check step (e): R = Res_X(X^-a A(Z X^v), X^-a C(Z X^v)) is nonzero, one v-weight, and in ideal (A,C) (Groebner, r=2,3)
import random, itertools, sympy as sp
random.seed(3)
def run(r,v,trials=3):
    Z=sp.symbols('z1:%d'%(r+1)); X=sp.Symbol('X')
    pts=list(itertools.product(range(3),repeat=r))
    for _ in range(trials):
        S=random.sample(pts,4)
        A=sum(random.randint(-5,5) or 1 for _ in [0]) and sum((random.randint(1,7))*sp.prod([z**e for z,e in zip(Z,u)]) for u in S)
        C=sum((random.randint(-7,-1) if i%2 else random.randint(1,7))*sp.prod([z**e for z,e in zip(Z,u)]) for i,u in enumerate(S))
        wts=[sum(a*b for a,b in zip(u,v)) for u in S]; a=min(wts)
        sub={z:z*X**vj for z,vj in zip(Z,v)}
        At=sp.expand(A.subs(sub,simultaneous=True)/X**a); Ct=sp.expand(C.subs(sub,simultaneous=True)/X**a)
        assert sp.expand(At.subs(X,1)-A)==0
        R=sp.expand(sp.resultant(At,Ct,X))
        P=sp.Poly(R,*Z); wset={sum(a*b for a,b in zip(m,v)) for m in P.monoms()}
        G=sp.groebner([A,C],*Z,order='grevlex')
        rem=G.reduce(R)[1]
        d=sp.degree(At,X)
        print('r',r,'v',v,'d',d,'R!=0',R!=0,'#terms',len(P.monoms()),'v-weights',wset,'pred',d*d+2*a*d,'R in (A,C):',sp.expand(rem)==0)
run(2,(1,3)); run(2,(2,-1)); run(3,(1,2,-1),2)
