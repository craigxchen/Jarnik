import sys; sys.path.insert(0,'.')
from circ import *
from fam_pell import pell_blocks, words, tri
def prod(L): return reduce(lambda a,b:a*b,L)
def coeffmat(marks):
    D=max(m.degree() for m in marks)
    return sp.Matrix([[ (m.as_expr().expand().coeff(t,k)) for k in range(D+1)] for m in marks]).T
print('--- linear independence of four-mark Pell families (2h>D)')
for n in [1,2,3]:
    ms=words(pell_blocks(n,[1,3,5]),tri)
    print('n',n,'D',ms[0].degree(),'h',contact(ms),'rank',coeffmat(ms).rank())
    ms5=words(pell_blocks(n,[1,3,5,8]),tri+[(0,3)])
    print('   five-mark: D',ms5[0].degree(),'h',contact(ms5),'rank',coeffmat(ms5).rank())
print('--- D=4 cluster (2h=D boundary): relation among five marks, reality, conj-stability')
a=[0,1,3,7];c=1
marks=[P(prod([t-a[k]-sp.I*c if k in S else t-a[k]+sp.I*c for k in range(4)])) for S in combinations(range(4),2)]
five=marks[:5]
ns=coeffmat(five).nullspace()
for v in ns:
    v=v/[x for x in v if x!=0][0]
    print('relation',list(v))
    cv=[sp.conjugate(x) for x in v]
    s=sum((ci*m.as_expr() for ci,m in zip(cv,five)),0)
    print('  conj relation applied:',sp.expand(s))
print('--- tight triple with r=0: G1=t^d-w, G2=t^d+w, marks G1G2, conj(G1G2), G1*conj(G2)')
w=sp.Rational(3,5)+sp.I*sp.Rational(4,5)
for d in [1,2,3]:
    G1=P(t**d-w);G2=P(t**d+w)
    ms=[G1*G2,bar(G1)*bar(G2),G1*bar(G2)]
    res,T=analyze(ms,f'tight triple d={d}')
