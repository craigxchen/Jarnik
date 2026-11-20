"""K_C of (VIII.2): min sum of l over sets S of primes 3 mod 4 with sum log l/(l+1) > log C.
Exact optimum by branch and bound over small primes (initial-segment greedy also reported)."""
from math import log, sqrt
from sympy import primerange
P = [l for l in primerange(3, 400000) if l % 4 == 3]
w = [log(l)/(l+1) for l in P]
def initial(T):
    s=0; K=0
    if T < 0: return 0, 0
    for i,l in enumerate(P):
        s+=w[i]; K+=l
        if s > T: return i+1, K
def optimal(T, nmax=60):
    # branch and bound on the first nmax primes: minimise cost subject to weight > T
    if T < 0: return 0
    best=[initial(T)[1]]
    suf=[0]*(nmax+1)
    for i in range(nmax-1,-1,-1): suf[i]=suf[i+1]+w[i]
    def rec(i, s, c):
        if s > T:
            best[0]=min(best[0], c); return
        if i>=nmax or c+P[i] >= best[0] or s+suf[i] <= T: return
        rec(i+1, s+w[i], c+P[i]); rec(i+1, s, c)
    rec(0,0.0,0)
    return best[0]
print("C  log C  |S|init  K_init  K_opt  K_init(with free (1+i): T=logC-log2/2)  K with mod-2 block (T=logC-3log2/4, +1)")
for C in [0.99, 1.0, 1.2, 1.3, sqrt(2), 1.5, 1.7, 2, 2.5, 3, 4, 5, 5.9, 6, 8, 8.24, 10]:
    T=log(C)
    n,K=initial(T); Ko=optimal(T) if T<1.0 else None
    n2,K2=initial(T-0.5*log(2)); n3,K3=initial(T-0.75*log(2))
    print(f"{C:.4f} {T:.4f} {n} {K} {Ko} {K2} {K3+1 if T-0.75*log(2)>=0 else 1}")
# Mertens constant in the class 3 mod 4: sum_{l<=x} log l/(l+1) - (1/2)log x
import itertools
s=0
for i,l in enumerate(P):
    s+=w[i]
    if l in (P[999], P[9999], P[-1]): print("x=",l," sum - 0.5 log x =", s-0.5*log(l))
