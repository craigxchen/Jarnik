from math import log
from sympy import primerange
P=[l for l in primerange(3,400000) if l%4==3]
def Kinit(T):
    if T<0: return 0
    s=0;K=0
    for l in P:
        s+=log(l)/(l+1);K+=l
        if s>T: return K
