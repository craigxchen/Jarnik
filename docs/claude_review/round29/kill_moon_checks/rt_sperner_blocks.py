"""kill_moon: 'Sperner blocks' over R(t). D=8, c=0 (s=4, m=3). Four disjoint pairs B_j={a_j,b_j} in the upper
half-plane with equal (Im a^r + Im b^r)_{r=1,2,3}. Then for the 6 sign patterns that are + on two blocks and - on
two, sum_k eps_k Im(a_k^r)=0 for r<=3, so by Lemma 4.5 the 6 points lie in one ball deg(z-z')<=4.
Per block: y2=t1-y1, x2=(t2/2-x1*y1)/y2, third equation is a polynomial in y1 (exact, sympy)."""
import sympy as sp, mpmath as mp, itertools
mp.mp.dps=60
y1=sp.symbols('y1')
def im_pow(x,y,r): return sp.im(sp.expand((x+sp.I*y)**r))
a0=(sp.Integer(0),sp.Integer(1)); b0=(sp.Integer(1),sp.Integer(2))
tau=[im_pow(*a0,r)+im_pow(*b0,r) for r in (1,2,3)]
print("tau =",tau)
blocks=[(a0,b0)]
for x1 in [sp.Rational(-1),sp.Rational(2),sp.Rational(-2),sp.Rational(3),sp.Rational(1,2)]:
    y2=tau[0]-y1; x2=(tau[1]/2-x1*y1)/y2
    e3=sp.together(3*x1**2*y1-y1**3+3*x2**2*y2-y2**3-tau[2])
    num=sp.numer(e3); P=sp.Poly(num,y1)
    roots=[r for r in sp.real_roots(P) if 0<r<tau[0]]
    for r in roots:
        yy2=tau[0]-r; xx2=x2.subs(y1,r)
        blk=((x1,r),(sp.nsimplify(xx2) if False else xx2,yy2))
        blocks.append(blk); break
    if len(blocks)==4: break
A=[]
for blk in blocks:
    for (x,y) in blk: A.append(mp.mpc(mp.mpf(sp.N(x,70)),mp.mpf(sp.N(y,70))))
print("roots a_k:",[mp.nstr(a,12) for a in A])
assert all(a.imag>0 for a in A)
assert min(abs(A[i]-A[j]) for i in range(8) for j in range(i)) > 1e-6
def polyfromroots(rs):
    c=[mp.mpc(1)]
    for r in rs:
        c=[mp.mpc(0)]+c; 
        for k in range(len(c)-1): c[k]-=r*c[k+1]
    return c  # low->high
pts=[]
for S in itertools.combinations(range(4),2):
    rs=[]
    for j,(blkA) in enumerate(range(4)):
        a,b=A[2*j],A[2*j+1]
        if j in S: rs+= [a,b]
        else: rs+= [mp.conj(a),mp.conj(b)]
    pts.append((S,polyfromroots(rs)))
for (S,z),(T,w) in itertools.combinations(pts,2):
    d=[abs(z[k]-w[k]) for k in range(9)]
    assert max(d[5:])<mp.mpf(10)**-50, (S,T,d)
    assert d[4]>1e-6
print("6 distinct points; all pairwise differences have zero coefficients in degrees 5,6,7 (to 1e-50) and nonzero degree-4 coefficient")
# also: N = prod (t-a)(t-abar) real, z*zbar=N by construction.  Plotkin bound for D=8: 2D=16.
