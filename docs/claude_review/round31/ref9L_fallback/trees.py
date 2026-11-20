# Contact-integer splits at p: T (subset of [7], 2<=|T|<=6) is visible at level t if T is exactly a residue class of x mod p^t.
import random, itertools, math
from collections import Counter
def splits(x,p,maxlev=60):
    S={}
    for t in range(1,maxlev):
        m=p**t; cl={}
        for i,v in enumerate(x): cl.setdefault(v%m,[]).append(i)
        for c in cl.values():
            if 2<=len(c)<=6: S[tuple(c)]=S.get(tuple(c),0)+1
            if len(c)==7: S[('top',)]=1
        if all(len(c)==1 for c in cl.values()): break
    return S
random.seed(1)
stat={2:Counter(),3:Counter(),5:Counter()}
big={2:Counter(),3:Counter()}
for trial in range(20000):
    while True:
        x=[random.randrange(-10**6,10**6) for _ in range(7)]
        if len(set(x))==7:
            g=0
            for v in x: g=math.gcd(g,v-x[0])
            if g==1: break
    for p in (2,3,5):
        S=splits(x,p); n=len([T for T in S if T!=('top',)])
        stat[p][n]+=1
        if p in big: big[p][len([T for T in S if T!=('top',) and len(T)>=3])]+=1
print("number of distinct visible splits divisible by p, over 20000 random primitive 7-tuples:")
for p in stat: print(p, sorted(stat[p].items()))
print("number of those with |T|>=3:"); 
for p in big: print(p, sorted(big[p].items()))
