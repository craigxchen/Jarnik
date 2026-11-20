# (a) two-class -> square roots mod 4q; (b) root sets -> two classes mod Q2 after valuation division; exact, small q.
import random, math
from sympy import factorint
from sympy.ntheory.modular import crt
random.seed(3); bad=0; n=0
for q in range(2,3000):
    f=factorint(q); pe=[p**e for p,e in f.items()]
    for _ in range(3):
        r1=[random.randrange(m) for m in pe]; r2=[random.choice([random.randrange(m),r1[k]]) for k,m in enumerate(pe)]
        a=int(crt(pe,r1)[0]); b=int(crt(pe,r2)[0])
        for y in range(-q,2*q):
            if all(y%m in (r1[k],r2[k]) for k,m in enumerate(pe)):
                n+=1
                if ((2*y-a-b)**2-(a-b)**2)%(4*q): bad+=1
print("(a) admissible y checked",n,"failures",bad)
# (b) for b<=1500, all a: roots in [0,b); u divides every root, roots/u lie in <=2 classes mod each p^f || Q2, length sqrt(b)/u <= sqrt(2 Q2)
fails=0; cnt=0
for B in range(2,1501):
    f=factorint(B)
    sq={}
    for x in range(B): sq.setdefault(x*x%B,[]).append(x)
    for a,roots in sq.items():
        u=1; Q1=1
        for p,e in f.items():
            am=a%p**e
            v=e if am==0 else min(e,factorint(am).get(p,0))
            if v==e: u*=p**((e+1)//2)
            else: t=v//2; u*=p**t; Q1*=p**(e-2*t)
        Q2=Q1//math.gcd(Q1,2)
        cnt+=1
        if any(x%u for x in roots): fails+=1; continue
        for p,e in factorint(Q2).items():
            cl={ (x//u)%p**e for x in roots}
            if len(cl)>2: fails+=1
        if B> u*u*Q1 or (Q2 and math.sqrt(B)/u > math.sqrt(2*Q2)+1e-9): fails+=1
print("(b) classes checked",cnt,"failures",fails)
