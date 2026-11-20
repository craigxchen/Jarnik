# Judge spot check of the attack referee's Prop 2.3: for four distinct integers, q2 <= sqrt2 * L^(5/2).
# q2 = lcm over bipartitions of gcd of within-side differences (attack Lemma 1.1, re-implemented here).
import random, math, itertools
from math import gcd
random.seed(11)
def q2(y):
    B=len(y); res=1
    for mask in range(1<<(B-1)):
        g=0
        for side in (0,1):
            S=[y[i] for i in range(B) if (mask>>i&1)==side]
            for v in S: g=gcd(g,v-S[0])
        if g: res=res*g//gcd(res,g)
    return res
worst=0; n=0
for L in range(3,400):
    for _ in range(60):
        mid=random.sample(range(1,L),2); y=[0]+sorted(mid)+[L]
        r=math.log(q2(y))/math.log(L); worst=max(worst,r); n+=1
        assert q2(y)**2 <= 2*L**5, y
# exhaustive small L
for L in range(3,120):
    for a,b in itertools.combinations(range(1,L),2):
        y=[0,a,b,L]; assert q2(y)**2<=2*L**5,y; n+=1
print('Prop 2.3 checked on',n,'four-point sets; max log q2/log L = %.3f'%worst)
