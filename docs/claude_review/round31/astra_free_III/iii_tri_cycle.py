"""Generate the equal-cofactor four-cycle and certify it entirely in integers."""
from fractions import Fraction
from math import gcd, isqrt, prod


def add(z, w):
    return (z[0]+w[0], z[1]+w[1])


def mul(z, w):
    return (z[0]*w[0]-z[1]*w[1], z[0]*w[1]+z[1]*w[0])


def conj(z):
    return (z[0], -z[1])


def norm(z):
    return z[0]**2+z[1]**2


def exact_div(z, w):
    h = norm(w)
    a, b = mul(z,conj(w))
    assert a % h == b % h == 0
    return (a//h,b//h)


def gaussian_gcd(z, w):
    while norm(w):
        h = norm(w)
        a,b = mul(z,conj(w))
        quotient = ((2*a+h)//(2*h),(2*b+h)//(2*h))
        m,n = mul(quotient,w)
        z,w = w,(z[0]-m,z[1]-n)
    return z


def crt(congruences):
    x, modulus = 0,1
    for b,n in congruences:
        assert gcd(n,modulus)==1
        x += modulus*(((b-x)*pow(modulus,-1,n)) % n)
        modulus *= n
    return x+modulus, modulus


def four_cycle(M):
    assert M % 4 == 1 and all(M % j for j in range(2,isqrt(M)+1))
    assert M not in (5,37,61)
    nonresidue = next(j for j in range(2,M) if pow(j,(M-1)//2,M)==M-1)
    i = pow(nonresidue,(M-1)//4,M)
    assert (i*i+1) % M == 0
    v = 10*M
    u = 20*M
    p,q = u*u+1,v*v+1
    x,L = crt([(1,2),(1,5),(2+2*u,p),(-2+2*v,q),(-4*pow(3,-1,M)+2*i,M)])
    assert L == 10*p*q*M and L<=x<2*L
    s = x*x
    d = s*s+64
    plus, minus = (x+2)**2+4,(x-2)**2+4
    assert plus*minus == d and d%M
    factors = [(1,0),(s,8)]
    C = []
    real_a = []
    for r in [Fraction(5),Fraction(5*p,q)]:
        k,l=r*plus,minus/r
        assert k.denominator == l.denominator == 1
        k,l = k.numerator,l.numerator
        a,b,c=(k+s)//2,(l+s)//2,(k+l+2*s)//4
        assert a%2 == b%2 == c%2 == 1 and c%M==0
        assert a*a+16==k*c and b*b+16 == l*c
        assert a+b==2*c and a*b-16 == s*c
        factors.append((a,4))
        real_a.append(a)
        C.append(c)
    assert (real_a[0]-real_a[1]) % M == 0
    # The least common denominator of F_j / conjugate(F_j) gives the integer circle.
    den = (1,0)
    for f in factors:
        assert norm(gaussian_gcd(f,conj(f))) == 1
        c = conj(f)
        den = mul(den,exact_div(c,gaussian_gcd(den,c)))
    points = [exact_div(mul(den,f),conj(f)) for f in factors]
    N = norm(den)
    assert all(norm(z) == N for z in points) and len(set(points)) == 4
    assert N*M <= d*prod(C)
    # Compare every reduced Gaussian star ratio to the predicted labels, up to sign and inversion.
    output = []
    for i,j in [(0,2),(0,3),(1,2),(1,3)]:
        z,w = points[i],points[j]
        g = gaussian_gcd(z,w)
        fz = exact_div(z,g)
        fw = exact_div(w,g)
        assert mult_is_real_unit(fz,fw)
        output.append(min(abs(fz[0]),abs(fz[1])))
    assert output == [4]*4
    # Angular width <=16/s. Integer check verifies a shortened critical arc.
    assert (32**4)*N < s**4
    print(f'M={M}, x={x}, N={N}')
    print('parent digits',len(str(N)),'cross offsets',output)
    print('floor(s / (16*N^(1/4))) at least',s//(16*(isqrt(isqrt(N))+1)))
    print('exact multiplier gain d*C1*C2/N',d*prod(C)//N)
    return N


def mult_is_real_unit(z,w):
    return w in (conj(z),(-z[0],z[1]))


if __name__=='__main__':
    four_cycle(1000000009)
