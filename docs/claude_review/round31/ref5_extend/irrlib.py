# Check, on random hyperplane-side families for random V (containing 1) and random T (small and large):
#  (a) F_i - F_j = P_J(1) P_J(Y+a)   (b) delta_ij := deg(F_i-F_j) - deg gcd(F_i,F_j) >= iota(P_J(1)) (odd-degree-prime part)
#  (c) the margin identity sum_{i<j} (D'/2 - deg gcd) <= s D'/4 for the reduced family (holds for ANY family with common norm,
#      cluster or not, as an identity on gcds: sum_{i<j} deg gcd >= D' s(s-2)/4)
#  (d) v_P(P_J(1)) >= |J cap K_P| unless one of the two forms vanishes on K_P = ker(ev_P|V).
import random, itertools, sys
from f2lib import *
def irreducibles(maxdeg):
    irr=[]
    for p in range(2,1<<(maxdeg+1)):
        if all(divmod2(p,q)[1] for q in irr if 2*deg(q)<=deg(p)): irr.append(p)
    return irr
IRR=irreducibles(10)
def vP(x,P):
    k=0
    while True:
        q,r=divmod2(x,P)
        if r: return k
        x=q;k+=1
def iota(x):  # degree of the part of x composed of odd-degree irreducibles (all factors here have degree <=10)
    s=0
    for P in IRR:
        if deg(P)%2==1:
            s+=deg(P)*vP(x,P)
    return s
