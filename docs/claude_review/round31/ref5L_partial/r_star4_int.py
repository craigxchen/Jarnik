# Independent integer specialisation of the four-block star at a NEW shape, plus Gaussian flips.
import sys; sys.path.insert(0,'.')
from glib import *
from itertools import combinations
import sympy as sp, math, random
T=sp.symbols('T')
shape=sys.argv[1] if len(sys.argv)>1 else 'new'
if shape=="astra": w=[sp.Integer(1)]+[sp.Integer(a) for a in (-9,-7,2,4)]
else: w=[sp.Integer(1)]+[sp.Integer(a) for a in (-9,-5,7,8)]
W=w[1:]
S=[sum(sp.prod(c) for c in combinations(W,m)) for m in (1,2,3,4)]
assert S[2]==S[0]*S[1] and S[0]*S[3]!=0
u=[a*a for a in w]
E=list(combinations(range(5),2))
assert len({abs(w[i]*w[j]) for i,j in E})==10
G=sp.expand(sp.prod([T-sp.I*a for a in W])); A=-T/S[0]+sp.I
F=[sp.expand(G+(ui-1)/ui*(ui+T**2)*A) for ui in u]
coef=[sp.Poly(f,T).all_coeffs() for f in F]
dens=set()
for cs in coef:
    for c in cs:
        dens.add(sp.re(c).q); dens.add(sp.im(c).q)
for i,j in E: dens.add((w[i]*w[j]).q)
D=1
for dd in dens: D=sp.ilcm(D,dd)
D=int(D)
if shape=='astra': D=6350400
print('dilation D =',D)
cre=[[ (int(sp.re(c)*D**k), int(sp.im(c)*D**k)) for k,c in enumerate(cs)] for cs in coef]
for cs in coef:
    for k,c in enumerate(cs): assert (sp.re(c)*D**k).is_Integer and (sp.im(c)*D**k).is_Integer
def col(i,m):
    return (sum(a*m**(4-k) for k,(a,b) in enumerate(cre[i])), sum(b*m**(4-k) for k,(a,b) in enumerate(cre[i])))
cof={}
for i,j in E: cof[i,j]=sp.Integer(D)**6*S[3]*(1/u[i]-1/u[j])
print('cofactors', [cof[e] for e in E])
assert all(c.is_Integer for c in cof.values())
bad=set(sp.primerange(2,24))
for c in cof.values(): bad|=set(sp.factorint(abs(int(c))))
bad|=set(sp.factorint(D))
L=prod(bad)
random.seed(1)
tries=0
while True:
    tries+=1
    m=1+L*random.randrange(1,40) if shape=='new' else 3175200023
    V=[col(i,m) for i in range(5)]
    blocks={ (i,j): m*m+int(D*w[i]*w[j])**2 for i,j in E}
    ok=all(gcd(*v)==1 for v in V)
    if not ok: continue
    if any(gcd(blocks[a],blocks[b])!=1 for a,b in combinations(E,2)) and shape=='new': continue
    badc=[(e,f,gcd(blocks[e],int(cof[f]))) for e in E for f in E if gcd(blocks[e],int(cof[f]))!=1]
    if shape!='new': print('block/cofactor common factors:',badc, flush=True)
    if badc and shape=='new': continue
    if shape!='new': print('pairwise coprime blocks:', all(gcd(blocks[a],blocks[b])==1 for a,b in combinations(E,2)))
    print('candidate',m,flush=True)
    fac={e:sp.factorint(blocks[e]) for e in E}
    if any(ex>1 for e in E for ex in fac[e].values()): 
        if shape!='new': print('not squarefree'); break
        continue
    break
print('tries',tries,'m digits',len(str(m)))
assert all(p%4==1 for e in E for p in fac[e])
N=prod(blocks.values()); print('parent digits',len(str(N)),'#primes',sum(len(fac[e]) for e in E))
n=[prod(blocks[e] for e in E if i in e) for i in range(5)]
for (i,j) in E:
    g=math.gcd(n[i],n[j]); assert g==blocks[i,j]
    assert det(V[i],V[j])==cof[i,j]*blocks[i,j], (i,j)
iso=[(i,e) for e in E for i in e if nrm(V[i])%blocks[e]==0]
print('whole-block isotropic incidences:',iso)
plev=[(i,e,p) for e in E for i in e for p in fac[e] if nrm(V[i])%p==0]
print('prime-level isotropic incidences count:',len(plev),' on blocks', sorted({x[1] for x in plev}))
print('|V|^2/n -1 :',[float(Fraction(nrm(V[i]),n[i])-1) for i in range(5)])
P=m  # block scale ~ m^2 ; P_block ~ N^(1/10)
Pb=N**(1/10) if False else None
import mpmath as mp; mp.mp.dps=50
Nm=mp.mpf(N)
print('defect exponents log|delta_i|/log(N^(1/10)):',[ (mp.log(abs(nrm(V[i])-n[i]))/(mp.log(Nm)/10)) if nrm(V[i])!=n[i] else 'exact' for i in range(5)])
Delta2=max(chord2(V[i],V[j]) for i,j in E)
print('Delta * N^(3/10) =', mp.sqrt(mp.mpf(Delta2.numerator)/Delta2.denominator)*Nm**(mp.mpf(3)/10))
# leading direction check: (y_j-1) ~ const * D_j  ; here y_j-1 ~ 2i(theta_j-theta_1), theta diff = asin(det/(|V1||Vj|))
th=[mp.asin(mp.mpf(det(V[0],V[j]))/mp.sqrt(mp.mpf(nrm(V[0])*nrm(V[j])))) for j in range(1,5)]
Dd=[w[j]**-2-1 for j in range(1,5)]
print('ratios theta_j/D_j (should be ~equal):',[mp.nstr(th[a]/mp.mpf(sp.N(Dd[a],40)),12) for a in range(4)])
# Gaussian flips on central edges (0,1),(0,2)
marked={e:set(e) for e in E}
Vf=list(V)
for e in [(0,1),(0,2)]:
    J=blocks[e]
    gam=ggcd(Vf[0],(J,0))
    assert nrm(gam)==J
    for i in e: assert gexactdiv(Vf[i],gam)  # divisible
    Vf=[ gexactdiv(Vf[i],gam) if i in e else gmul(Vf[i],gconj(gam)) for i in range(5)]
    marked[e]=set(range(5))-set(e)
nf=[prod(blocks[e] for e in E if i in marked[e]) for i in range(5)]
print('flipped n_i digits:',[len(str(x)) for x in nf])
assert all(gcd(*v)==1 for v in Vf)
inc=0; tot=0
for e in E:
    for i in marked[e]:
        tot+=1; inc+= nrm(Vf[i])%blocks[e]==0
for i,j in E:
    sh=prod(blocks[e] for e in E if i in marked[e] and j in marked[e])
    assert sh==math.gcd(nf[i],nf[j])
    assert det(Vf[i],Vf[j])==cof[i,j]*sh
    assert chord2(Vf[i],Vf[j])==chord2(V[i],V[j])
print('mixed orientation: isotropic incidences',inc,'of',tot,'marked; ratios, cofactors, |V|^2/n preserved:',
      all(Fraction(nrm(Vf[i]),nf[i])==Fraction(nrm(V[i]),n[i]) for i in range(5)))
