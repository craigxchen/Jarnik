import sympy as sp, itertools, math
P,q=sp.symbols('P q')
x=[-27,-270,135,120,40]; s1,s2,s3,s5=-2,-53925,-5042250,4723920000; l=97200; b0=3*s5//2
L=l*q; B=b0*q**2
U=[P+q*(xi**3-s1*xi**2+s2*xi) for xi in x]
V=[(P*U[i]+B*x[i], L*U[i]) for i in range(5)]
J={(i,j):P+q*x[i]*x[j]*(s1-x[i]-x[j]) for i,j in itertools.combinations(range(5),2)}
n=[sp.expand(math.prod(J[tuple(sorted((i,j)))] for j in range(5) if j!=i)) for i in range(5)]
a,b,c=sp.symbols('a b c')
sol=sp.solve([a*V[i][0]**2+b*V[i][0]*V[i][1]+c*V[i][1]**2-n[i] for i in range(3)],[a,b,c],dict=True)[0]
sol={k:sp.factor(v) for k,v in sol.items()}
print('Q* coefficients:',sol)
print('exact at marks 4,5:',[sp.simplify(sol[a]*V[i][0]**2+sol[b]*V[i][0]*V[i][1]+sol[c]*V[i][1]**2-n[i]) for i in (3,4)])
D=sp.factor(sol[b]**2-4*sol[a]*sol[c])
print('disc(Q*) =',D)
Ex=sp.factor((1-sol[a])*(1-sol[c])-(sol[b]/2)**2); print('det E =',Ex)
