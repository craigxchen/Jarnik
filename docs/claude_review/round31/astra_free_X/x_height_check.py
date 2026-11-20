"""Run python3 x_height_check.py for the exact fixed-parent ray certificates.

Labels are all split primes in [P,2P] or the 5 mod 8 subclass.
The relative-twist bound is Norm(gamma)<=2P with X>=2P.
Sieve units and conjugates together, which permutes the four tested units.
"""
from collections import Counter
from math import gcd, isqrt


def primes_to(n):
    sieve = bytearray([1]) * (n + 1)
    sieve[:2] = bytes(2)
    for d in range(2, isqrt(n) + 1):
        if sieve[d]:
            sieve[d*d:n+1:d] = bytes((n-d*d)//d+1)
    return [k for k in range(2, n+1) if sieve[k]]


def gauss(p):
    for x in range(1, isqrt(p)+1):
        y = isqrt(p-x*x)
        if x*x+y*y == p:
            return x, y
    raise ValueError(p)


def mult(z, w):
    x, y = z
    a, b = w
    return x*a-y*b, x*b+y*a


def external_lcm(X, labels):
    powers = []
    for t in primes_to(X):
        if t == 2 or t in labels:
            continue
        u = t
        while u <= X//t:
            u *= t
        powers.append(u)
    while len(powers) > 1:
        new = [a*b for a,b in zip(powers[::2], powers[1::2])]
        if len(powers) % 2:
            new.append(powers[-1])
        powers = new
    return powers[0] if powers else 1


def one_ray(alpha, p, multiplier, nmax):
    """First exponent where the four slabs have no common assigned angle.

    The slabs cut out a rational rectangle for (cos(2n phi),sin(2n phi))
    scaled to radius p^n. Compute its nearest and farthest squared radii.
    Returning zero leaves joint feasibility over different n unresolved.
    """
    z = (1,0)
    radius = 1
    for n in range(1,nmax+1):
        z = mult(z,alpha)
        radius *= p
        x, y = z
        d0,d1,d2,d3 = [gcd(multiplier,v) for v in (y,x-y,x,x+y)]
        left,right = 2*d2*d2-radius,radius-2*d0*d0
        low,high = d3*d3-radius,radius-d1*d1
        if left>right or low>high:
            return n
        near = max(0,left,-right)**2 + max(0,low,-high)**2
        far = max(abs(left),abs(right))**2 + max(abs(low),abs(high))**2
        if not near <= radius*radius <= far:
            return n
    return 0


def all_ray_certificate(alpha,p,multiplier,w):
    """Return last exponent checked. Every larger one passes by height.

    w is a signed Gaussian prime over 5. Its powers have four nonzero
    axis/diagonal coordinates for every n, by unique factorization.
    """
    assert sum(t*t for t in w)==5<p
    z=u=(1,0)
    pn=hn=1
    n=0
    while pn < hn*multiplier*multiplier:
        n += 1
        z,u = mult(z,alpha),mult(u,w)
        pn,hn = pn*p,hn*5
        sx,sy = z
        tx,ty = u
        source = (sy,sx-sy,sx,sx+sy)
        target = (ty,tx-ty,tx,tx+ty)
        assert all(target)
        if any(pn*t*t < hn*gcd(v,multiplier)**2 for v,t in zip(source,target)):
            return None
    return n


def parent_labels(P, only_five):
    return [p for p in primes_to(2*P) if p>=P and
            (p%8==5 if only_five else p%4==1)]


def candidates(labels):
    coords = {p: gauss(p) for p in labels}
    for p in labels:
        for ell in labels:
            a, b = coords[ell]
            for sign in (1,-1):
                direction = (a,sign*b)
                alpha = mult(coords[p],mult(direction,direction))
                if gcd(*alpha)!=1:
                    assert ell==p and sign==-1
                    continue
                assert sum(t*t for t in alpha)==p*ell*ell
                yield p,ell,alpha


def census(P,only_five,A,nmax=12):
    labels=parent_labels(P,only_five)
    s=len(labels)
    X=3*s**A
    assert X>=2*P
    multiplier=external_lcm(X,set(labels))
    counts=Counter()
    passes=[]
    for p in labels:
        assert one_ray(gauss(p),p,multiplier,nmax)==0
    for p,ell,alpha in candidates(labels):
        n=one_ray(alpha,p,multiplier,nmax)
        counts[n]+=1
        if n==0:
            for w in ((2,1),(1,-2),(1,2),(2,-1)):
                n=all_ray_certificate(alpha,p,multiplier,w)
                if n is not None:
                    passes.append((p,ell,alpha,w,n))
                    break
    assert sum(counts.values())==s*(2*s-1)
    return s,X,tuple(sorted(counts.items())),passes


def certify_least_cutoff(P,only_five,X,expected_last):
    """Exhaust all lifts just at X and verify the infinite survivor at X-1."""
    labels=parent_labels(P,only_five)
    assert X>2*P
    full=external_lcm(X,set(labels))
    prior=external_lcm(X-1,set(labels))
    assert full>prior
    examined=0
    last_candidates=[]
    for p,ell,alpha in candidates(labels):
        assert one_ray(alpha,p,full,24)>0
        if one_ray(alpha,p,prior,24)==0:
            last_candidates.append((p,ell,alpha))
        examined+=1
    assert last_candidates==[expected_last]
    assert examined==len(labels)*(2*len(labels)-1)
    p,ell,alpha=expected_last
    assert expected_last in candidates(labels)
    tail=all_ray_certificate(alpha,p,prior,(2,1))
    assert tail is not None
    print("certified least cutoff",P,"5mod8" if only_five else "all",
          X,"exclusions",examined,"prior survivor",p,ell,alpha,"through",tail)


if __name__=="__main__":
    print("P, pool, s, X, (first excluded ray exponent, lifts); 0 = unresolved in 12 tests")
    for P, only_five,A in ((50,False,3),(100,False,2),(100,True,3),
                         (200,False,2),(200,True,3),
                         (500,False,2),(500,True,3),
                         (1000,False,2),(1000,True,2)):
        s,X,counts,passes=census(P,only_five,A)
        print(P,"5mod8" if only_five else "all",s,X,counts)
        for p,ell,alpha,w,n in passes:
            print("  all rays pass:",p,ell,alpha,w,"tested through",n)

    for P,half,X,last in (
        (100,False,439,(181,113,(985,-1158))),
        (500,False,3881,(997,557,(-17482,-1923))),
        (1000,False,8681,(1601,1381,(39869,-38260))),
        (1000,True,5783,(1621,1493,(-34698,-49085))),
    ):
        certify_least_cutoff(P,half,X,last)
