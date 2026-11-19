"""Post-process enum_chart output: L_0, per-decade best K5/exponent, cut profile of the best."""
import sys, re, math, itertools
from fractions import Fraction
from math import gcd, lcm
from cutprofile import tuple_of, profile
rows=[]
for fn in sys.argv[1:]:
    for line in open(fn):
        m=re.match(r'K5=([\d.]+) H=(\d+) L=(\d+) X=\(([\d, ]+)\) e=(-?[\d.]+)',line)
        if m:
            K5=float(m[1]);H=int(m[2]);L=int(m[3]);X=tuple(int(t) for t in m[4].split(','));e=float(m[5])
            rows.append((K5,H,L,X,e))
def L0(L,X):
    dens=[Fraction(x,L).denominator for x in X]
    for a,b in itertools.combinations(X,2):
        q=(a*b+L*L)//(a-b) if a>b else (a*b+L*L)//(b-a)
        dens.append(Fraction(q,L).denominator)
    r=1
    for d in dens: r=lcm(r,d)
    return r
# dedupe by ratio class
seen={};uniq=[]
for r in sorted(rows):
    K5,H,L,X,e=r; key=(H,frozenset(Fraction(x,L) for x in X))
    if key in seen: continue
    seen[key]=1; uniq.append(r)
print('unique clusters',len(uniq))
l0s=[L0(r[2],r[3]) for r in uniq]
print('L_0 values not divisible by 6:',sorted(set(l for l in l0s if l%6)))
narrow=[r for r in uniq if r[3][0]>=r[2]]   # arc angle <= pi/2
print('with A>=L (arc <= pi/2):',len(narrow),' min K5 %.4f'%min(r[0] for r in narrow))
print('largest H:',max(r[1] for r in uniq))
print('\nper decade of H (A>=L): count, min K5, max exponent e=-log(Delta)/log H (target 0.3)')
dec={}
for r in narrow: dec.setdefault(int(math.log10(r[1])),[]).append(r)
for d in sorted(dec):
    b=min(dec[d]); mx=max(dec[d],key=lambda r:r[4])
    pw,ws,fac=profile(*tuple_of(mx[2],mx[3]))
    dev=max(abs(v-0.1) for v in pw.values()); npair=sum(1 for v in pw.values() if v>0)
    print('1e%-2d n=%-4d minK5=%.3f  best e=%.4f (K5=%.3f H=%d L=%d X=%s) nonzero pair cuts=%d maxdev=%.3f singleton/W=%.3f'%(
        d,len(dec[d]),b[0],mx[4],mx[0],mx[1],mx[2],mx[3],npair,dev,ws))
print('\nclosest to the ten-pair profile (all ten pair cuts present), A>=L:')
best=[]
for r in narrow:
    pw,ws,fac=profile(*tuple_of(r[2],r[3]))
    if all(v>0 for v in pw.values()):
        best.append((max(abs(v-0.1) for v in pw.values()),ws,r))
best.sort()
for dev,ws,r in best[:15]:
    print('maxdev=%.3f singleton/W=%.3f K5=%.3f e=%.4f H=%d L=%d X=%s'%(dev,ws,r[0],r[4],r[1],r[2],r[3]))
print('number with all ten pair cuts nonzero:',len(best),'of',len(narrow))
