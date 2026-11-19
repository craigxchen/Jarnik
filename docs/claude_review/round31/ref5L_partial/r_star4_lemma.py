# Four-block star at a NEW shape (-9,-5,7,8), specialised exactly as in the lemma (dilation D*L, T=1 mod L).
import sys; sys.path.insert(0,'.')
from glib import *
import sympy as sp, math, random
from itertools import combinations
T=sp.symbols('T')
w=[sp.Integer(1)]+[sp.Integer(a) for a in (-9,-5,7,8)]
W=w[1:]; S=[sum(sp.prod(c) for c in combinations(W,m)) for m in (1,2,3,4)]
assert S[2]==S[0]*S[1]
u=[a*a for a in w]; E=list(combinations(range(5),2))
G=sp.expand(sp.prod([T-sp.I*a for a in W])); A=-T/S[0]+sp.I
F=[sp.expand(G+(ui-1)/ui*(ui+T**2)*A) for ui in u]
coef=[sp.Poly(f,T).all_coeffs() for f in F]
dens=set()
for cs in coef:
    for c in cs: dens|={sp.re(c).q,sp.im(c).q}
D=1
for d in dens: D=sp.ilcm(D,d)
D=int(D)
cof0={e:S[3]*(1/u[e[0]]-1/u[e[1]]) for e in E}
bad=set(sp.primerange(2,24))|set(sp.factorint(D))
for c in cof0.values(): bad|=set(sp.factorint(abs(c.p)))|set(sp.factorint(c.q))
for a,b in combinations(E,2): bad|=set(sp.factorint(abs(int(u[a[0]]*u[a[1]]-u[b[0]]*u[b[1]]))))
for e in E: bad|=set(sp.factorint(abs(int(w[e[0]]*w[e[1]]))))
L=math.prod(int(p) for p in bad); D=D*L
print('bad primes',sorted(bad),' D digits',len(str(D)))
cre=[[(int(sp.re(c)*D**k),int(sp.im(c)*D**k)) for k,c in enumerate(cs)] for cs in coef]
def col(i,m): return (sum(a*m**(4-k) for k,(a,b) in enumerate(cre[i])), sum(b*m**(4-k) for k,(a,b) in enumerate(cre[i])))
icof={e:sp.Integer(D)**6*cof0[e] for e in E}; assert all(c.is_Integer for c in icof.values())
random.seed(3)
import mpmath as mp; mp.mp.dps=50
for SC in (10**3,10**6):
  for trial in range(500):
    m=1+L*(SC*D//L+random.randrange(1,10**6))
    V=[col(i,m) for i in range(5)]
    blocks={e:m*m+(D*int(w[e[0]]*w[e[1]]))**2 for e in E}
    ok=all(math.gcd(*v)==1 for v in V) and all(math.gcd(blocks[a],blocks[b])==1 for a,b in combinations(E,2)) \
       and all(math.gcd(blocks[e],abs(int(icof[f])))==1 for e in E for f in E) \
       and all(blocks[e]%(p*p) for e in E for p in sp.primerange(2,10**5)) \
       and all(math.gcd(m, D*int(w[e[0]]*w[e[1]]))==1 for e in E)
    if ok: break
  N=math.prod(blocks.values()); n=[math.prod(blocks[e] for e in E if i in e) for i in range(5)]
  for (i,j) in E:
    assert math.gcd(n[i],n[j])==blocks[i,j] and det(V[i],V[j])==int(icof[i,j])*blocks[i,j]
  iso=[(e,i) for e in E for i in e if nrm(V[i])%blocks[e]==0]
  Dl2=max(chord2(V[i],V[j]) for i,j in E)
  print('m/D ~',SC,'trial',trial,'parent digits',len(str(N)),'\n isotropic:',iso,
        '\n gcds at other incidences:',{math.gcd(nrm(V[i]),blocks[e]) for e in E for i in e if (e,i) not in iso},
        '\n |V|^2/n:',[float(Fraction(nrm(V[i]),n[i])) for i in range(5)],
        '\n log10(Delta*N^0.3)=',mp.nstr(mp.log10(mp.sqrt(mp.mpf(Dl2.numerator)/Dl2.denominator))+mp.mpf(3)/10*mp.log10(mp.mpf(N)),10),
        ' max cofactor digits',max(len(str(abs(int(c)))) for c in icof.values()),
        '\n defect exps (units P=N^(1/10)):',[mp.nstr(mp.log(abs(nrm(V[i])-n[i]))/(mp.log(mp.mpf(N))/10),5) if nrm(V[i])!=n[i] else 'exact' for i in range(5)])
