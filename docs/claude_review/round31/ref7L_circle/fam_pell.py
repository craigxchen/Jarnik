import sys; sys.path.insert(0,'.')
from circ import *
def pell_blocks(n, slopes):
    out=[]
    for s in slopes:
        x,y=P(1),P(s)
        for _ in range(n):
            x,y=P(t*t+1)*x+P(t)*y, P(t)*x+y
        out.append(x+P(sp.I)*y)
    return out
def words(blocks,edges):
    beta=P(t+2*sp.I)
    first=beta*reduce(lambda a,b:a*b,[bar(g) for g in blocks])
    return [first]+[bar(beta)*reduce(lambda a,b:a*b,[g if j in e else bar(g) for j,g in enumerate(blocks)]) for e in edges]
tri=[(0,1),(0,2),(1,2)]
for n in [1,2,3]:
    res,T=analyze(words(pell_blocks(n,[1,3,5,8]),tri+[(0,3)]),f'pell5 n={n}')
    res,T=analyze(words(pell_blocks(n,[1,3,5]),tri),f'pell4 n={n}')
for d in [1,3,6]:
    bl=[P(t**d-a+sp.I) for a in [1,3,8]]
    analyze([reduce(lambda a,b:a*b,[g if k==j else bar(g) for k,g in enumerate(bl)]) for j in range(3)],f'three d={d}')
    bl=[P(t**d-a+sp.I) for a in [1,-1,3]]
    analyze([reduce(lambda a,b:a*b,[g if k==j else bar(g) for k,g in enumerate(bl)]) for j in range(3)],f'three(+-) d={d}')
