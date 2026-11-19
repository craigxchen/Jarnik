# Split-form ansatz: V_i=(X_i,1), n_i = X_i^2 - w^2 = L_i^+ L_i^-,
# L_i^+ = prod of the two C+ blocks at i, L_i^- = prod of the two C- blocks, L^+ - L^- = 2w.
import sympy as sp, itertools
E=[(1,2),(1,3),(1,4),(1,5),(2,3),(2,4),(2,5),(3,4),(3,5),(4,5)]
b={e:sp.Symbol('b%d%d'%e) for e in E}
w=sp.Symbol('w')
def ed(i,j): return (min(i,j),max(i,j))
Cp=[ed(1,2),ed(2,3),ed(3,4),ed(4,5),ed(5,1)]
Cm=[e for e in E if e not in Cp]
def nb(i,C): return [e for e in C if i in e]
lin=[];quad=[]
for i in range(1,6):
    p=nb(i,Cp); m=nb(i,Cm)
    lin.append(b[p[0]]+b[p[1]]-b[m[0]]-b[m[1]])
    quad.append(b[p[0]]*b[p[1]]-b[m[0]]*b[m[1]]-2*w)
print('Cm',Cm)
M=sp.Matrix([[sp.diff(l,b[e]) for e in E] for l in lin]); print('rank lin',M.rank())
# normalise translation: b12=0
sol=sp.solve(lin,[b[ed(1,3)],b[ed(2,4)],b[ed(3,5)],b[ed(1,4)],b[ed(2,5)]],dict=True)
print(sol)
s=sol[0]
Q=[sp.expand((q+2*w).subs(s)) for q in quad]
for q in Q: print(sp.factor(q))
D=[sp.factor(Q[i]-Q[0]) for i in range(1,5)]
print('differences:');[print(d) for d in D]
fr=[b[e] for e in Cp]
G=sp.groebner([Q[i]-Q[0] for i in range(1,5)],*fr,order='lex')
print('GB',G.exprs)
sols=sp.solve(G.exprs,fr,dict=True)
for so in sols:
    full={e:sp.simplify(b[e].subs(s).subs(so)) for e in E}
    vals=list(full.values())
    distinct=len(set(vals))
    print(so,' #distinct b:',distinct, ' w=',sp.factor(Q[0].subs(so)/2))
