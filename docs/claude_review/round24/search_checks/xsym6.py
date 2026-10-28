# express the degree-24 factor in S_4-invariants e2,e3,e4 (e1 = sum b = 0), exactly
import sympy as sp, pickle, random, itertools
b1,b2=sp.symbols('b1 b2')
disc,cub,den=pickle.load(open('xsym4.pkl','rb'))
fl=sp.factor_list(disc)[1]
f=[g for g,e in fl if sp.Poly(g,b1,b2).total_degree()==24][0]
pickle.dump(f,open('f24.pkl','wb'))
b0=sp.Integer(1); b3=-1-b1-b2; bs=[b0,b1,b2,b3]
E2=sp.expand(sum(bs[i]*bs[j] for i,j in itertools.combinations(range(4),2)))
E3=sp.expand(sum(bs[i]*bs[j]*bs[k] for i,j,k in itertools.combinations(range(4),3)))
E4=sp.expand(b0*b1*b2*b3)
mons=[(i,j,k) for i in range(13) for j in range(9) for k in range(7) if 2*i+3*j+4*k==24]
print('monomials',len(mons))
cs=sp.symbols('c0:%d'%len(mons))
random.seed(3); eqs=[]
for trial in range(len(mons)+15):
    u=sp.Rational(random.randint(-30,30),random.randint(1,9)); v=sp.Rational(random.randint(-30,30),random.randint(1,9))
    e2,e3,e4=[E.subs({b1:u,b2:v}) for E in (E2,E3,E4)]
    eqs.append(sum(c*e2**i*e3**j*e4**k for c,(i,j,k) in zip(cs,mons))-f.subs({b1:u,b2:v}))
sol=sp.linsolve(eqs,cs)
sol=list(sol)
print('solutions:',len(sol))
if sol:
    s=sol[0]; e2,e3,e4=sp.symbols('e2 e3 e4')
    F=sp.expand(sum(c*e2**i*e3**j*e4**k for c,(i,j,k) in zip(s,mons)))
    print('F nonzero terms',len(F.as_ordered_terms()))
    pickle.dump(F,open('F24.pkl','wb'))
    print(sp.factor(F))
