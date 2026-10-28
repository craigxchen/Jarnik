import sympy as sp, pickle, flint, time
b0,b1,b2,X = sp.symbols('b0 b1 b2 X')
Ln,Ld=pickle.load(open('c1_L.pkl','rb'))
LP=sp.Poly(Ln,X)
ctx=flint.fmpz_mpoly_ctx.get(('b0','b1','b2'),'lex')
def tof(e):
    P=sp.Poly(sp.expand(e),b0,b1,b2)
    return ctx.from_dict({m:int(c) for m,c in P.terms()})
A={}
for m,c in zip(LP.monoms(),LP.coeffs()):
    A[m[0]//2-1]=tof(c)
A3,A2,A1,A0=A[3],A[2],A[1],A[0]
print('Ld factor:', [ (str(f)[:60],e) for f,e in tof(Ld).factor()[1]])
for name,Aa in (('A3',A3),('A2',A2),('A1',A1),('A0',A0)):
    c,ff=Aa.factor(); print(name,c,[(str(f)[:70],e) for f,e in ff])
t0=time.time()
D = A2**2*A1**2-4*A3*A1**3-4*A2**3*A0-27*A3**2*A0**2+18*A3*A2*A1*A0
print('disc nterms',len(D),'deg',D.total_degree(), time.time()-t0)
c,facs=D.factor()
print('factor time',time.time()-t0, 'content',c)
res=[]
for f,e in facs:
    print('factor total degree',f.total_degree(),'mult',e,'nterms',len(f), str(f)[:150])
    res.append((str(f),e))
pickle.dump(res,open('c1b_factors.pkl','wb'))
