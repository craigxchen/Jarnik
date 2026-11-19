import sys; sys.path.insert(0,'.')
from circ import *
def prod(L): return reduce(lambda a,b:a*b,L)
T=((-3,-1,-1,1,2),(-2,-1,-1,-1,3),(-3,-2,1,1,1))
marks=[P(prod([t-sp.I*s for s in m])) for m in T]
for m in marks: print(m.as_expr())
print('pairwise diff degrees',[ (a-b).degree() for a,b in combinations(marks,2)])
res,tri=analyze(marks,'imag-axis triple')
print('2h=',2*res['h'],' D=',res['D'],' (circle bound for triples needs 2h<=D+2u; here u=',tri[0][3],')')
T5=((-5,-1,-1,2,3,4),(-4,-3,1,1,2,5),(-3,-2,-1,-1,4,5),(-4,-2,-1,1,3,5),(-5,-2,1,1,3,4))
m5=[P(prod([t-sp.I*s for s in m])) for m in T5]
res,tri=analyze(m5,'imag-axis five')
print('20h=',20*res['h'],' 9D=',9*res['D'])
