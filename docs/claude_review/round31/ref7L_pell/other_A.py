import sys; sys.path.insert(0,'.')
from pell_core import *
def A_of(a): return [[padd(pmul(a,a),P(1)),a],[a,P(1)]]
def run(a,ns,seeds,edges,beta):
    A=A_of(a); r=len(seeds)
    Gs=[block(matpow(A,ns[l]),seeds[l]) for l in range(r)]; Gb=[pconj(x) for x in Gs]
    marks=[pmul(beta,prod(Gb))]+[pmul(pconj(beta),prod([Gs[l] if l in e else Gb[l] for l in range(r)])) for e in edges]
    return analyse(marks)
S=[Fr(1,2),Fr(1,3),Fr(2,7),Fr(-5,3)]; E=[(0,1),(0,2),(1,2),(0,3)]
for name,a in [("t",P(0,1)),("t^2",P(0,0,1)),("t+1",P(1,1)),("t^2+t-3",P(-3,1,1)),("t^3",P(0,0,0,1))]:
    beta=padd(a,P(G(0,2)))
    for ns in [(1,1,1,1),(2,2,2,2),(1,2,2,2),(2,2,2,1),(2,1,1,1)]:
        x=run(a,ns,S,E,beta)
        print(f"a={name:8} n={ns} D={x['D']:3} h={x['h']:3} 2h-D={2*x['h']-x['D']:3} gcd={x['gcddeg']} eqnorm={x['eqnorm']} distinct={x['distinct']}")
