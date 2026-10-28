import random,sys
from rlib import *
rng=random.Random(12345)
nc=0; mism=0
for trial in range(150):
    k=rng.choice([5,6,7,8])
    x=rng.sample(range(-60,61),k)
    P=[(xi,1) for xi in x]
    A,_=contacts_minformula(P)
    B=contacts_classes(x)
    for S in A:
        nc+=1
        if A[S]!=B[S]: mism+=1
    # anchored at label 0 -> infinity: Moebius y=1/(x-x0); P_i=(1, x_i-x0) ~ point 1/(x_i-x0)
    # integral chart: X_i = L/(x_i-x0)
    from math import lcm
    L=1
    for i in range(1,k): L=lcm(L,x[i]-x[0])
    X=[L//(x[i]-x[0]) for i in range(1,k)]
    d,delta=contacts_gcd_anchor(X)
    for S in A:
        nc+=1
        if d[S]!=A[S]: mism+=1
print("comparisons",nc,"mismatches",mism)
