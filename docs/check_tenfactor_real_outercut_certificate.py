"""Rational contraction certificate for a real ten-factor outer-cut witness.

All acceptance checks use Fraction arithmetic. No numerical solver or
approximate root is trusted, and no rational/Gaussian realization is asserted.
"""

from fractions import Fraction as Q
from math import prod

S = (
    (-1,-1,-1,1,-1,-1,1,-1,1,1),
    (-1,1,-1,-1,1,-1,1,1,1,-1),
    (1,-1,-1,-1,1,-1,-1,1,-1,1),
    (1,1,1,-1,-1,-1,1,-1,-1,1),
    (1,1,-1,-1,-1,1,-1,-1,1,-1),
)
A = tuple(tuple(Q(s-t,2) for s,t in zip(row,S[0])) for row in S[1:])
# a=c+E*x; x=(a_4,a_6,a_7,a_8), with zero-based subscripts.
C = tuple(map(Q,(-1,Q(-7,5),Q(7,5),-1,0,1,0,0,0,Q(-2,5))))
E = tuple(tuple(map(Q,row)) for row in (
    (1,1,1,0), (1,0,1,-1), (0,-1,0,1), (2,0,2,-1),
    (1,0,0,0), (0,0,0,0), (0,1,0,0), (0,0,1,0),
    (0,0,0,1), (0,0,0,0),
))
X0 = tuple(Q(v,10**12) for v in
           (-48557368900,927791797667,805289163286,174670899854))
R = Q(1,10**8)


def mv(matrix,vector):
    return tuple(sum(a*b for a,b in zip(row,vector)) for row in matrix)


def inverse(matrix):
    n=len(matrix)
    aug=[list(row)+[Q(i==j) for j in range(n)] for i,row in enumerate(matrix)]
    for j in range(n):
        p=next(i for i in range(j,n) if aug[i][j])
        aug[p],aug[j]=aug[j],aug[p]
        z=aug[j][j]; aug[j]=[a/z for a in aug[j]]
        for i in range(n):
            if i!=j:
                z=aug[i][j]
                aug[i]=[a-z*b for a,b in zip(aug[i],aug[j])]
    return tuple(tuple(row[n:]) for row in aug)


def main():
    # These identities eliminate all four linear moments exactly.
    assert mv(A,C)==(0,)*4
    assert all(mv(A,tuple(row[j] for row in E))==(0,)*4 for j in range(4))
    a0=tuple(c+z for c,z in zip(C,mv(E,X0)))
    F0=mv(A,tuple(z**3 for z in a0))
    J0=tuple(tuple(sum(A[i][j]*3*a0[j]**2*E[j][k] for j in range(10))
                   for k in range(4)) for i in range(4))
    B=inverse(J0)
    assert all(mv(B,tuple(row[j] for row in J0))==tuple(Q(i==j) for i in range(4))
               for j in range(4))
    eta=max(map(abs,mv(B,F0)))
    bnorm=max(sum(map(abs,row)) for row in B)
    lips=tuple(sum(map(abs,row)) for row in E)
    # For ||x-X0||_infty<=R, |a_j-a0_j|<=lips_j R.
    # Bound the row-sum norm of DF(x)-DF(X0), retaining all coefficients.
    jvariation=max(sum(abs(A[i][j])*3*
                       (2*abs(a0[j])*lips[j]*R+lips[j]**2*R**2)*lips[j]
                       for j in range(10)) for i in range(4))
    contraction=bnorm*jvariation
    assert eta<Q(1,10**10)
    assert contraction<Q(1,1000)
    assert eta+contraction*R<R
    # T(x)=x-BF(x) is therefore a strict contraction of this closed box
    # into itself. Banach's theorem supplies a unique exact zero in it.
    amin=tuple(abs(z)-l*R for z,l in zip(a0,lips))
    amax=tuple(abs(z)+l*R for z,l in zip(a0,lips))
    assert min(amin)>Q(48,1000)
    assert all(amax[i]<amin[j] or amax[j]<amin[i]
               for i in range(10) for j in range(i))
    assert all(max(amin[i]-amax[j],amin[j]-amax[i])>Q(63,5000)
               for i in range(10) for j in range(i))
    fifth0=mv(S,tuple(z**5 for z in a0))
    fifth_error=sum(5*(abs(z)+l*R)**4*l*R for z,l in zip(a0,lips))
    order=(2,3,1,0,4)
    assert all(fifth0[j]-fifth0[i]>2*fifth_error for i,j in zip(order,order[1:]))
    assert all(S[i][6]==-1 for i in (2,4))
    assert all(S[i][6]==1 for i in (0,1,3))
    # The scale-invariant formal endpoint statistic requested in the note.
    span0=fifth0[4]-fifth0[2]
    span_low,span_high=span0-2*fifth_error,span0+2*fifth_error
    p_low,p_high=prod(amin),prod(amax)
    assert span_low**2>25*Q(693,100)**2*p_high
    assert span_high**2<25*Q(694,100)**2*p_low
    print("PASS: exact rational contraction maps the radius-1e-8 box strictly inside itself")
    print("PASS: residual < 1e-10; derivative contraction < 1/1000")
    print("PASS: ten nonzero distinct absolute values; fifth order (2,3,1,0,4)")
    print("PASS: column 6 isolates the two extreme rows")
    print("PASS: fifth-span/(5 sqrt(abs(product a))) lies strictly between 6.93 and 6.94")


if __name__=='__main__':
    main()
