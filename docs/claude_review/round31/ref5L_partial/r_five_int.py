import sys; sys.path.insert(0,'.')
from glib import *
import pickle, sympy as sp, math
from itertools import combinations
F,h,E,cof=pickle.load(open('five_F.pkl','rb'))
T=sp.symbols('T',real=True)
coef=[sp.Poly(f,T).all_coeffs() for f in F]
dens=set()
for cs in coef:
    for c in cs: dens|={sp.re(c).q, sp.im(c).q}
A={};Bb={}
for e in E:
    A[e]=-1/(1+h[e]**2); Bb[e]=h[e]/(1+h[e]**2); dens|={A[e].q,Bb[e].q}
D=1
for d in dens: D=sp.ilcm(D,d)
# need D^k * coeff_k integral: D=lcm works since D^k multiple of D for k>=1
D=int(D); print('D digits',len(str(D)), 'D factor', sp.factorint(D))
cre=[[(int(sp.re(c)*D**k),int(sp.im(c)*D**k)) for k,c in enumerate(cs)] for cs in coef]
for cs in coef:
    for k,c in enumerate(cs): assert (sp.re(c)*D**k).is_Integer and (sp.im(c)*D**k).is_Integer
def col(i,m): return (sum(a*m**(4-k) for k,(a,b) in enumerate(cre[i])), sum(b*m**(4-k) for k,(a,b) in enumerate(cre[i])))
icof={e:sp.Integer(D)**6*cof[e] for e in E}
print('integral cofactors?',all(c.is_Integer for c in icof.values()))
bad=set(sp.primerange(2,24))|set(sp.factorint(D))
for c in icof.values(): bad|=set(sp.factorint(abs(int(c)), limit=10**6))
# include pairwise block resultants
Jp={e:sp.expand((T-A[e])**2+Bb[e]**2) for e in E}
for a,b in combinations(E,2):
    r=sp.resultant(Jp[a],Jp[b],T); bad|=set(sp.factorint(abs(int(sp.numer(r))),limit=10**6))
L=1
for p in bad: L*=int(p)
print('bad primes',len(bad))
D=D*L  # second dilation as in the lemma
cre=[[(int(sp.re(c)*D**k),int(sp.im(c)*D**k)) for k,c in enumerate(cs)] for cs in coef]
icof={e:sp.Integer(D)**6*cof[e] for e in E}
from collections import Counter; why=Counter()
import random; random.seed(7); SCALE=int(sys.argv[1])
found=None
for trial in range(2000):
    m=1+L*(D//L*SCALE+random.randrange(1,10**6))
    V=[col(i,m) for i in range(5)]
    if not all(math.gcd(*v)==1 for v in V): why['prim']+=1; continue
    blocks={e:(m-int(D*A[e]))**2+int(D*Bb[e])**2 for e in E}
    if any(math.gcd(blocks[a],blocks[b])!=1 for a,b in combinations(E,2)): why['cop']+=1; continue
    if any(math.gcd(blocks[e],abs(int(icof[f])))!=1 for e in E for f in E): why['cof']+=1; continue
    # trial squarefree up to 10^5
    sqf=all(blocks[e]%(p*p) for e in E for p in sp.primerange(2,10**5))
    if not sqf: why['sqf']+=1; continue
    found=m; break
print('m =',found,'trial',trial,why)
m=found; V=[col(i,m) for i in range(5)]
for e in E: assert blocks[e]==D*D*((sp.Rational(m,D)-A[e])**2+Bb[e]**2)
N=math.prod(blocks.values()); print('parent digits',len(str(N)))
n=[math.prod(blocks[e] for e in E if i in e) for i in range(5)]
for (i,j) in E:
    assert math.gcd(n[i],n[j])==blocks[i,j]
    assert det(V[i],V[j])==int(icof[i,j])*blocks[i,j]
iso=[(e,i) for e in E for i in e if nrm(V[i])%blocks[e]==0]
print('isotropic (edge,mark) 0-based:',iso)
nonisogcd={(e,i):math.gcd(nrm(V[i]),blocks[e]) for e in E for i in e if nrm(V[i])%blocks[e]}
print('gcds at non-isotropic incidences:',set(nonisogcd.values()))
print('|V|^2/n:',[float(Fraction(nrm(V[i]),n[i])) for i in range(5)])
print('blocks mod 4:',{b%4 for b in blocks.values()})
import mpmath as mp; mp.mp.dps=60
Dl2=max(chord2(V[i],V[j]) for i,j in E)
lg=mp.log(mp.sqrt(mp.mpf(Dl2.numerator)/Dl2.denominator))/mp.log(mp.mpf(N))
print('log Delta/log N =',mp.nstr(lg,8),'log10(Delta*N^0.3)=',mp.nstr(mp.log10(mp.sqrt(mp.mpf(Dl2.numerator)/Dl2.denominator))+mp.mpf(3)/10*mp.log10(mp.mpf(N)),10),' cofactor max digits',max(len(str(abs(int(c)))) for c in icof.values()))
