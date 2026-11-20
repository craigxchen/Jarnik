# Own helpers: exact Gaussian integer arithmetic, exact rational rank.
from fractions import Fraction
from math import gcd, isqrt
def gmul(z,w): return (z[0]*w[0]-z[1]*w[1], z[0]*w[1]+z[1]*w[0])
def gconj(z): return (z[0],-z[1])
def gnorm(z): return z[0]*z[0]+z[1]*z[1]
def gpow(z,k):
    r=(1,0)
    for _ in range(k): r=gmul(r,z)
    return r
def gdivexact(z,w):
    n=gnorm(w); a,b=gmul(z,gconj(w)); assert a%n==0 and b%n==0; return (a//n,b//n)
def two_squares(p):
    for a in range(1,isqrt(p)+1):
        b2=p-a*a; b=isqrt(b2)
        if b*b==b2 and b>0: return (a,b)
    raise ValueError
def rank_Q(rows):
    """exact rank over Q of an integer/fraction matrix (list of lists)."""
    A=[[Fraction(x) for x in r] for r in rows]
    if not A: return 0
    m=len(A); n=len(A[0]); r=0
    for c in range(n):
        piv=None
        for i in range(r,m):
            if A[i][c]!=0: piv=i;break
        if piv is None: continue
        A[r],A[piv]=A[piv],A[r]
        for i in range(m):
            if i!=r and A[i][c]!=0:
                f=A[i][c]/A[r][c]
                A[i]=[x-f*y for x,y in zip(A[i],A[r])]
        r+=1
        if r==m: break
    return r
def kernel_Q(rows):
    """basis of {c : sum_i c_i rows[i] = 0} over Q (left kernel)."""
    import sympy
    M=sympy.Matrix(rows).T
    return [list(v) for v in M.nullspace()]
def aff_rank(rows):
    if len(rows)<=1: return 0
    r0=rows[0]
    return rank_Q([[x-y for x,y in zip(r,r0)] for r in rows[1:]])
