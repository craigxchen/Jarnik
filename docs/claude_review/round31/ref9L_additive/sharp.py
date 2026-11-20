#!/usr/bin/env python3
"""Sharpness of IX.1 (exponent 2/3) and of the parallelogram remark (exponent 1/2):
q = Q1 Q2 Q3 with Q_i ~ q^{1/3} pairwise coprime; y_x = x1 Q2Q3 + x2 Q1Q3 + x3 Q1Q2
has at most two classes mod each p^f || q, signature = affine image of F_2^3
(an affine 3-space), span = Q2Q3+Q1Q3+Q1Q2 ~ 3 q^{2/3}. Also check Walsh divisibility."""
import itertools, math
from sympy import factorint, nextprime
def sigs(ys, q):
    pps = sorted(p**f for p,f in factorint(q).items())
    cls = []
    for m in pps:
        r = sorted(set(y % m for y in ys)); assert len(r) <= 2; cls.append(r)
    return pps, [tuple(int(y % m != r[0]) for m, r in zip(pps, cls)) for y in ys]
for (Q1,Q2,Q3) in [(101,103,107),(2**7,3**4,5**3),(1009,1013,1019*7)]:
    q=Q1*Q2*Q3
    X=list(itertools.product((0,1),repeat=3))
    ys=[x[0]*Q2*Q3+x[1]*Q1*Q3+x[2]*Q1*Q2 for x in X]
    pps,S=sigs(ys,q)
    span=max(ys)-min(ys)
    F={u:sum((-1)**(sum(a*b for a,b in zip(u,x)))*y for x,y in zip(X,ys)) for u in X if any(u)}
    print('q=%d span=%d span/q^(2/3)=%.3f q^(2/3)/4=%.0f distinct sigs=%d'%(q,span,span/q**(2/3),q**(2/3)/4,len(set(S))),
          'F_u!=0 for', [u for u in F if F[u]!=0], 'q/Q_u|F_u:', all(F[u]%(q//[Q1,Q2,Q3][u.index(1)])==0 for u in F if F[u]!=0 and sum(u)==1))
# parallelogram at span ~ 2 sqrt q
for (Q1,Q2) in [(101,103),(1009,1013)]:
    q=Q1*Q2; ys=[0,Q2,Q1,Q1+Q2]; pps,S=sigs(ys,q)
    print('q=%d parallelogram span=%d =%.3f sqrt q; sigs'%(q,Q1+Q2,(Q1+Q2)/math.sqrt(q)),S)
