"""Exact identities and full-circle checks for the Pell neighborhood theorem."""
from math import gcd, isqrt, prod

from check_direct_primitive_axis_lift import family
from check_gaussian_reflection_replacement import gmul, gconj, gnorm
from check_mobius_joint_anchor_residues import gaussian_prime, gpow


def factor_small(n):
    result = {}
    p = 2
    while p*p <= n:
        while n % p == 0:
            result[p] = result.get(p,0)+1
            n //= p
        p = 3 if p == 2 else p+2
    if n > 1:
        result[n] = result.get(n,0)+1
    return result


def complete_circle(factorization):
    points = [(1,0)]
    for p,e in factorization.items():
        assert p % 4 == 1
        pi = gaussian_prime(p)
        choices = [gmul(gpow(pi,j),gpow(gconj(pi),e-j)) for j in range(e+1)]
        points = [gmul(z,w) for z in points for w in choices]
    answer = set()
    for z in points:
        for unit in [(1,0),(-1,0),(0,1),(0,-1)]:
            answer.add(gmul(unit,z))
    assert len(answer) == 4*prod(e+1 for e in factorization.values())
    return answer


def matvec(T,z):
    p,q,r,s = T
    x,y = z
    return p*x+q*y,r*x+s*y


def forms(s,t):
    return s*(s-2*t-16), t*(2*s+t-2)


def check_identities(index):
    U,V,P,A,B,C,zs = family(index)
    a,b,c = gnorm(A),gnorm(B),gnorm(C)
    N = 5*a*b*c
    T = (2*V,V-U,-U-V,2*V)
    p,q,r,s = T
    assert p*s-q*r == -1
    inverse = (-s,q,r,-p)
    assert (p*p+r*r,p*q+r*s,q*q+s*s) == (a,b-a,b)
    assert a*b-(a-b)**2 == 1 and gcd(a,b) == 1
    assert c == 16*a-3*b
    assert (p*zs[0][0]+r*zs[0][1],q*zs[0][0]+s*zs[0][1]) == (-8*a,-b)
    expected = [(0,0),(16,0),(0,2),(4,-6)]
    for z,w in zip(zs,expected):
        diff = (z[0]-zs[0][0],z[1]-zs[0][1])
        assert matvec(inverse,diff) == w
        assert forms(*w) == (0,0)
    y0 = matvec(inverse,zs[0])
    v0 = matvec(inverse,(-zs[0][1],zs[0][0]))
    assert y0 == (b*b-9*a*b,7*a*b-8*a*a)
    assert v0 == (-b,8*a)
    assert all(gnorm(z) == N for z in zs)
    assert 14*V*V < a < 15*V*V
    assert 5*V*V < b < 6*V*V
    assert c > 206*V*V and N > 250**2*V**6
    assert max(map(abs,y0)) < 1800*V**4
    assert max(map(abs,v0)) < 120*V*V
    if V >= 16:
        # C=1/4: M <= (9/40)V+2 sqrt(V) <= (29/40)V < V.
        assert 4*V <= V*V//4
        assert V*(3*V+16) < b
        assert V*(3*V+2) < a
    return U,V,a,b,c,N,T,inverse,zs


def main():
    for index in range(1,31):
        check_identities(index)
    total = 0
    for index in range(1,5):
        U,V,a,b,c,N,T,inverse,zs = check_identities(index)
        ff = {}
        for integer in (5,a,b,c):
            for p,e in factor_small(integer).items():
                ff[p] = ff.get(p,0)+e
        assert prod(p**e for p,e in ff.items()) == N
        points = complete_circle(ff)
        inside = []
        other_d2 = []
        for w in points:
            assert gnorm(w) == N
            delta = (w[0]-zs[0][0],w[1]-zs[0][1])
            st = matvec(inverse,delta)
            F,H = forms(*st)
            assert a*F+b*H == 0
            assert F % b == 0 and H % a == 0
            d2 = gnorm(delta)
            # The fourth power removes both radius square roots exactly.
            if 256*d2*d2 <= N:
                inside.append(w)
                assert w in zs
            if w not in zs:
                other_d2.append(d2)
            total += 1
        print(f'index={index}, V={V}, full circle={len(points)}, chord ball={len(inside)}, minimum other squared chord={min(other_d2)}')
    print('PASS: unimodular basis, pencil forms, centers, and absolute constants on 30 Pell instances.')
    print(f'PASS: {total} full-circle points, including all partial Gaussian factor switches, on four circles.')
    print('The finite cases supplement the general divisibility and coordinate-bound proof.')


if __name__ == '__main__':
    main()
