import sys; sys.path.insert(0,'.')
from circ import *
R=sp.Rational
def prod(L): return reduce(lambda a,b:a*b,L)
# cube family (round 29 Theorem C)
y=[R(-9,4),R(-2,3),R(1,4),R(8,3)]; eps=[1,-1,-1,1]
r=[e/v+sp.I*v for e,v in zip(eps,y)]
roots=[r[0],-r[0],r[1],-r[1],r[2],-r[2],r[3],-r[3]]
Q=P(prod([t-v for v in roots]))
cl=[]
for mask in range(256):
    z=P(prod([t-sp.conjugate(v) if (mask>>k)&1 else t-v for k,v in enumerate(roots)]))
    if (z-Q).is_zero or (z-Q).degree()<=4: cl.append(z)
print('cube cluster size',len(cl))
worst=None
for S in combinations(range(len(cl)),5):
    res,_=analyze([cl[i] for i in S],f'cube {S}',triples=False,verbose=False)
    if worst is None or res['slack']<worst['slack']: worst=res
print('cube worst five-subset', worst)
analyze([cl[i] for i in (0,1,2,3,4)],'cube sample with triples')
# D=4 six-point cluster: roots a_k +- i c
for a,c in [([0,1,3,7],1),([R(-3),R(-1),R(1),R(3)],2),([0,1,2,5],R(1,2))]:
    marks=[]
    for S in combinations(range(4),2):
        marks.append(P(prod([t-a[k]-sp.I*c if k in S else t-a[k]+sp.I*c for k in range(4)])))
    print('D4 cluster h=',contact(marks))
    for S in combinations(range(6),5):
        analyze([marks[i] for i in S],f'D4 a={a} c={c} {S}')
