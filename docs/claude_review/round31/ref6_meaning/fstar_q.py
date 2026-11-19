"""Own implementation: transfer of two_thirds/proof.md Thm 5.2 to F_p[t] with the norm of the
constant extension F_{p^2}(t)/F_p(t) (p=2: A^2+AB+B^2; p=3: A^2+B^2).
Prime P of degree e, Q=p^e: e odd inert (fibre Q+1, lifts *Q), e even split (fibre Q-1, lifts *Q).
No ramified place (disc a unit, infinity inert of degree 1). log -> deg.
Bound: F*(M) <= M D/4 + C(M,2) c, W = 2D  =>  W >= 8F*/M - 4(M-1)c.
All quantities kept as exact integers (8F*)."""
import sys
from fractions import Fraction
from math import comb, log

def Emin(n, Qc):            # least number of colliding pairs, n objects in Qc classes
    b, r = divmod(n, Qc)
    return Qc*comb(b,2) + r*b

def mobius(n):
    res=1; d=2
    while d*d<=n:
        if n%d==0:
            n//=d
            if n%d==0: return 0
            res=-res
        d+=1
    return -res if n>1 else res

def nirr(p,e):
    return sum(mobius(d)*p**(e//d) for d in range(1,e+1) if e%d==0)//e

def tower(n, Q0, Q):        # sum over a>=1 of E(n, Q0*Q^(a-1))
    s=0; c=Q0
    while c < n:
        s += Emin(n,c); c*=Q
    return s

def psi8(M, Q0, Q):         # 8*Psi: min over level profiles (>=2 nonempty levels)
    w=[8*tower(n,Q0,Q) for n in range(M+1)]
    INF=float('inf'); best=[INF]*(M+1)
    for S in range(1,M):
        cut=(2*S-M)**2      # 8*(1/2)(S-M/2)^2
        b=w[S]
        for S0 in range(1,S):
            v=best[S0]+w[S-S0]
            if v<b: b=v
        best[S]=b+cut
    return min(best[S]+w[M-S] for S in range(1,M))

def F8(p, M):
    tot=0; e=1
    while True:
        Q=p**e
        if Q-1 >= M: break
        if e%2==1:
            tot += 8*e*nirr(p,e)*tower(M,Q+1,Q)
        else:
            tot += e*nirr(p,e)*min(8*tower(M,Q-1,Q), psi8(M,Q-1,Q))
        e+=1
    return tot

if __name__=='__main__':
    for p in (2,3):
        print(f"p={p}:  M  8F*  W_min=8F*/M  3M log_p M  ratio")
        for M in [4,5,6,8,16,32,64,128,256]:
            f8=F8(p,M); Wmin=Fraction(f8,M)
            print(f"  {M:4d} {f8:10d} {float(Wmin):10.1f} {3*M*log(M)/log(p):10.1f} {float(Wmin)/(3*M*log(M)/log(p)):.3f}")
