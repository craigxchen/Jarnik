"""Random genuine five-tuples over Q(i)[t] with contact h=2 or 3 imposed exactly by linear algebra.
Roots a_k +/- i b_k; mark j picks a_k + i s_jk b_k. p1 equal <=> sum_k s_jk b_k equal; p2 equal additionally
<=> sum_k s_jk a_k b_k equal. Also clusters with small integer data to make circles likely."""
import sys, random; sys.path.insert(0,'.')
from circ import *
R=sp.Rational
def prod(L): return reduce(lambda a,b:a*b,L)
random.seed(int(sys.argv[1]) if len(sys.argv)>1 else 0)
results=[]
trials=0
while len(results)<int(sys.argv[2] if len(sys.argv)>2 else 30):
    trials+=1
    D=random.randint(5,8); want=random.choice([2,3])
    S=[[1]*D]+[[random.choice([1,-1]) for _ in range(D)] for _ in range(4)]
    if len({tuple(s) for s in S})<5: continue
    M=sp.Matrix([[S[j][k]-S[0][k] for k in range(D)] for j in range(1,5)])
    ker=M.nullspace()
    if not ker: continue
    b=sum((random.randint(-3,3)*v for v in ker), sp.zeros(D,1))
    if any(x==0 for x in b): continue
    b=b*sp.lcm([sp.Rational(x).q for x in b])
    if want==3:
        M2=sp.Matrix([[ (S[j][k]-S[0][k])*b[k] for k in range(D)] for j in range(1,5)])
        ker2=M2.nullspace()
        if not ker2: continue
        a=sum((random.randint(-3,3)*v for v in ker2), sp.zeros(D,1))
        a=a*sp.lcm([sp.Rational(x).q for x in a])
    else:
        a=sp.Matrix([random.randint(-3,3) for _ in range(D)])
    marks=[P(prod([t-a[k]-sp.I*S[j][k]*b[k] for k in range(D)])) for j in range(5)]
    if len({m.as_expr() for m in marks})<5: continue
    g=reduce(lambda x,y:x.gcd(y),marks)
    if g.degree()>0: marks=[m.exquo(g) for m in marks]
    if marks[0].degree()<2: continue
    try:
        res,tri=analyze(marks,f'rand',verbose=False)
    except AssertionError as e:
        continue
    results.append(res)
    print(res['D'],res['h'],res['r'],'slack',res['slack'],'min triple slack',res['min_triple_slack'],flush=True)
print('min slack',min(r['slack'] for r in results),' min triple slack',min(r['min_triple_slack'] for r in results))
