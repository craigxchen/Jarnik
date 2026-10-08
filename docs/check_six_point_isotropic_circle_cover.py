"""Exact pencil, genus-five witness, and circle-lift certificates."""
from fractions import Fraction as F
from itertools import permutations, combinations
from functools import reduce
from math import gcd, isqrt, lcm
from check_integer_cotangent_spectral import mul,div,power,add,neg,primitive_vector
from check_integer_cotangent_normalization import edge_norm
from check_least_radius_formula import gcd_gaussian, exact_div, norm, gcd_all
from check_segre_gradient_arithmetic import rref
from check_joint_affine_relation_lattice import polynomial_gcd,qtrim

def padd(f,g):
    return qtrim([(f[i] if i<len(f) else 0)+(g[i] if i<len(g) else 0) for i in range(max(len(f),len(g)))])
def pmul(f,g):
    a=[F(0)]*(len(f)+len(g)-1)
    for i,x in enumerate(f):
        for j,y in enumerate(g): a[i+j]+=x*y
    return qtrim(a)
def pscale(f,c): return qtrim([c*x for x in f])
def pquot(f,g):
    f=qtrim(f); ans=[F(0)]*max(1,len(f)-len(g)+1)
    while f!=[0] and len(f)>=len(g):
        d=len(f)-len(g); c=f[-1]/g[-1]; ans[d]=c
        f=padd(f,[F(0)]*d+pscale(g,-c))
    assert f==[0]
    return qtrim(ans)
def determinant(m):
    ans=[0]; n=len(m)
    for p in permutations(range(n)):
        s=(-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
        ans=padd(ans,pscale(reduce(pmul,[m[i][p[i]] for i in range(n)],[1]),s))
    return ans

def peval(f,t):
    return sum(c*t**i for i,c in enumerate(f))

def square_root(x):
    x=F(x)
    a,b=isqrt(x.numerator),isqrt(x.denominator)
    assert a*a==x.numerator and b*b==x.denominator
    return F(a,b)

def probe(params):
    points=[(F(1-t*t,1+t*t),F(2*t,1+t*t)) for t in params]
    raw=[]
    for z in points:
        d=(1,0)
        for w in points:
            if z!=w:d=mul(d,add(z,neg(w)))
        raw.append(div(power(z,2),d))
    central,_=primitive_vector(raw)
    u=[z[0] or z[1] for z in central]
    U=[[F(1)]*6,[z[0] for z in points],[z[1] for z in points]]
    dot=lambda a,b:sum(x*y*c for x,y,c in zip(a,b,u))
    assert all(dot(a,b)==0 for a in U for b in U)
    for inds in combinations(range(6),3):
        m=[[U[i][j]*u[j] for j in inds]+[F(i==k) for k in range(3)] for i in range(3)]
        reduced,pivots=rref(m,3)
        if len(pivots)==3:break
    f0=[[F(0)]*6 for _ in range(3)]
    for i in range(3):
        for k,j in enumerate(inds): f0[i][j]=reduced[k][3+i]
    b=[[dot(a,c) for c in f0] for a in f0]
    fs=[[f0[i][j]-sum(b[i][k]*U[k][j]/2 for k in range(3)) for j in range(6)] for i in range(3)]
    assert all(dot(a,c)==0 for a in fs for c in fs)
    assert all(dot(U[i],fs[j])==int(i==j) for i in range(3) for j in range(3))
    x=[[U[1][j],fs[2][j]] for j in range(6)]
    y=[[U[2][j],-fs[1][j]] for j in range(6)]
    frame_denominator=lcm(*(c.denominator for row in x+y for c in row))
    x=[pscale(row,frame_denominator) for row in x]
    y=[pscale(row,frame_denominator) for row in y]
    vectors=[[[F(1)] for _ in range(6)],x,y]
    for a in vectors:
        for b in vectors:
            assert reduce(padd,[pscale(pmul(v,w),c) for v,w,c in zip(a,b,u)],[0])==[0]
    rows=[[[1],a,b,pmul(a,a),pmul(a,b),pmul(b,b)] for a,b in zip(x,y)]
    cof=[pscale(determinant([[row[k] for k in range(6) if k!=j] for row in rows[:5]]),(-1)**j) for j in range(6)]
    common=reduce(polynomial_gcd,cof)
    co=[pquot(p,common) for p in cof]
    for row in rows:
        assert reduce(padd,[pmul(a,b) for a,b in zip(co,row)],[0])==[0]
    disc=padd(pscale(pmul(co[3],co[5]),4),pscale(pmul(co[4],co[4]),-1))
    assert [len(p)-1 for p in cof]==[8,7,7,6,6,6]
    assert len(common)==1
    assert len(disc)==13
    assert len(polynomial_gcd(disc,[i*disc[i] for i in range(1,len(disc))]))==1

    c0,c1,c2,A,B,C=co
    K=padd(padd(pmul(C,pmul(c1,c1)),pscale(pmul(B,pmul(c1,c2)),-1)),
           padd(pmul(A,pmul(c2,c2)),pscale(pmul(c0,disc),-1)))
    triangles=[determinant([[[1],x[j],y[j]] for j in inds])
               for inds in combinations(range(5),3)]
    assert K==pscale(reduce(pmul,triangles,[1]),-1)
    assert len(K)==21 and len(polynomial_gcd(K,disc))==1
    assert all(sum(u[i] for i in inds)!=0
               for count in (2,3) for inds in combinations(range(6),count))

    # The initial configuration provides a rational point on the cover.
    c0,c1,c2,A,B,C=[p[0] for p in co]
    d=square_root(disc[0])
    affine_points=[(x[i][0],y[i][0]) for i in range(6)]
    lifted=[(d*(2*A*a+B*b+c1),d*d*b-B*c1+2*A*c2) for a,b in affine_points]
    common_norm=4*A*(C*c1*c1-B*c1*c2+A*c2*c2-c0*d*d)
    assert common_norm>0
    assert all(a*a+b*b==common_norm for a,b in lifted)
    for exponent in range(3):
        assert reduce(add,[mul((c,0),power(z,exponent)) for c,z in zip(u,lifted)],(0,0))==(0,0)
    primitive,_=primitive_vector(lifted)
    radius_squared=sum(c*c for c in primitive[0])
    ratios=[div(z,lifted[0]) for z in lifted]
    original=[div(z,points[0]) for z in points]
    assert ratios==original or ratios==[(z[0],-z[1]) for z in original]
    cotangents=[]
    for z,w in combinations(lifted,2):
        product=mul(z,(w[0],-w[1]))
        cotangents.append(product[1]/(common_norm-product[0]))
    scale=lcm(*(c.denominator for c in cotangents))
    assert lcm(*(edge_norm(int(c*scale),scale) for c in cotangents))==radius_squared
    h=gcd(gcd(int(A),int(B)),int(C))
    metric=gcd_gaussian((int(2*A),0),(int(B),int(d)))
    assert norm(metric)==4*abs(A)*h
    extra_numerator=exact_div((int(2*A*c2-B*c1),int(-d*c1)),metric)
    extra=gcd_gaussian((int(d),0),extra_numerator)
    coefficient_gcd=gcd_all([(int(2*A*d),0),(int(B*d),int(d*d)),
                             (int(d*c1),int(2*A*c2-B*c1))])
    assert norm(coefficient_gcd)==norm(metric)*norm(extra)
    row_gcd=gcd_all([tuple(int(c) for c in z) for z in lifted])
    factor=exact_div(row_gcd,coefficient_gcd)
    triangle_values=[int(determinant([[[1],x[j],y[j]] for j in inds])[0])
                     for inds in combinations(range(6),3)]
    minor_gcd=reduce(gcd,triangle_values,0)
    exact_div((abs(minor_gcd),0),factor)
    assert radius_squared==abs(K[0])/(h*norm(extra)*norm(factor))
    content_denominator=h*norm(extra)*norm(factor)
    for i,j in combinations(range(6),2):
        dx=affine_points[i][0]-affine_points[j][0]
        dy=affine_points[i][1]-affine_points[j][1]
        chord_form=A*dx*dx+B*dx*dy+C*dy*dy
        lifted_difference=add(lifted[i],neg(lifted[j]))
        primitive_difference=add(primitive[i],neg(primitive[j]))
        assert norm(lifted_difference)==4*A*d*d*chord_form
        assert F(norm(primitive_difference),radius_squared)==F(d*d*chord_form,K[0])
        assert norm(primitive_difference)==d*d*abs(chord_form)/content_denominator
    print('PASS: u=',u,'; degree-12 squarefree cover, degree-20 factorization, '
          'and exact radius/content/chord formulas.')
    return u,U,fs,co,disc
if __name__=='__main__':
    probe([0,1,2,3,4,5])
    probe([0,1,2,4,7,11])
    cases=0
    for b0 in range(-8,9):
        for d0 in range(1,9):
            value=b0*b0+d0*d0
            for A in range(1,value+1):
                if value%A:
                    continue
                C,B,d=value//A,2*b0,2*d0
                h=gcd(gcd(A,B),C)
                metric=gcd_gaussian((2*A,0),(B,d))
                assert norm(metric)==4*A*h
                for c1,c2 in ((1,2),(2*d,3*d),(A-C,B+1)):
                    extra=gcd_gaussian((d,0),exact_div((2*A*c2-B*c1,-d*c1),metric))
                    coefficient_gcd=gcd_all([(2*A*d,0),(B*d,d*d),(d*c1,2*A*c2-B*c1)])
                    assert norm(coefficient_gcd)==norm(metric)*norm(extra)
                    cases+=1
    print('PASS:',cases,'all-prime metric/coefficient content checks, including ramification at 2.')
