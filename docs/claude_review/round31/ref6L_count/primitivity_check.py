"""Astra's primitivity claim in the inverse bound: in the reduced word of a pair (i,j), all
separating layers at a prime point the same Gaussian way, so the word and every sub-product
of its cut blocks is a primitive plane vector.  Random primitive 5-tuples, exact arithmetic."""
import random, math, itertools
from math import gcd as igcd
from glib import *
from sympy import primerange
random.seed(11)
P1=[p for p in primerange(5,200) if p%4==1]
bad=0; tested=0
for trial in range(400):
    ps=random.sample(P1,random.randint(2,5)); pis=[split_prime(p) for p in ps]
    es=[random.randint(1,4) for _ in ps]
    a=[[random.randint(0,e) for e in es] for k in range(5)]
    for t in range(len(ps)):  # force primitivity at each prime
        a[0][t]=0; a[1][t]=es[t]
    for i,j in itertools.combinations(range(5),2):
        # blocks of the reduced word, indexed by cut (high side), with orientation
        blocks={}
        for t,(pi,e) in enumerate(zip(pis,es)):
            for tau in range(1,e+1):
                hi=frozenset(k for k in range(5) if a[k][t]>=tau)
                if (i in hi)==(j in hi): continue
                f=pi if j in hi else conj(pi)   # z_j/z_i gains pi where j is high
                blocks.setdefault(hi,[]).append(f)
        keys=list(blocks)
        for r in range(1,len(keys)+1):
            for sub in itertools.combinations(keys,r):
                U=(1,0)
                for kk in sub:
                    for f in blocks[kk]: U=mul(U,f)
                tested+=1
                if igcd(abs(U[0]),abs(U[1]))!=1: bad+=1
print('sub-products tested',tested,'non-primitive',bad)
