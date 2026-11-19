"""Finite exact check of the rounded-cofactor construction, not an all-character test."""
from fractions import Fraction
from itertools import combinations
from math import gcd, isqrt, prod


def prime(n):
    if n % 2 == 0:
        return n == 2
    return all(n % d for d in range(3, isqrt(n) + 1, 2))


def skew_core(q):
    # Border index 0. At positive columns, flipping border+diagonal gives core 2.
    square = {x*x % q for x in range(1, q)}
    c = []
    for y in range(q):
        s = [-1] + [1 if x == y or (y-x) % q in square else -1 for x in range(q)]
        flipped = [(-v if i in (0,y+1) else v) for i,v in enumerate(s)]
        c.extend((s,flipped))
    return list(zip(*c))


def check_residue_tree(b,q,depth,N):
    classes = {0:(0,0)}
    power = 1
    for k in range(1,depth+1):
        modulus = q*power
        new = {}
        for u,(x,y) in classes.items():
            keys = sorted({v % modulus for v in b if v % power == u})
            lifts = [(xx,yy) for xx in range(x,modulus,power)
                     for yy in range(y,modulus,power) if (xx*xx+yy*yy-N) % modulus == 0]
            assert len(lifts) >= len(keys)
            new.update(zip(keys,lifts))
        classes, power = new,modulus
        for i,j in combinations(range(len(b)),2):
            assert (classes[b[i]%power] == classes[b[j]%power]) == ((b[i]-b[j])%power == 0)


def main():
    m = 8
    base = 2**36 + 1
    p = []
    while len(p) < 2*(m-1):
        if prime(base):
            p.append(base)
        base += 4
    E = p[0]**(m//2)
    N = prod(p)
    b = [4*(2*m*m+i)**2 for i in range(m)]
    rows = skew_core(m-1)
    columns = list(zip(*rows))
    first, second = columns[::2],columns[1::2]
    H = [[sum(x*y for x,y in zip(a,b))**2 for b in second] for a in first]
    assert all(sum(row)==m*m for row in H) and all(sum(col)==m*m for col in zip(*H))
    assert all(H[i][i]==(m-4)**2 for i in range(m-1))
    assert all({H[i][j],H[j][i]}=={0,16} for i,j in combinations(range(m-1),2))
    labels = set()
    eps = Fraction(1,m**2)
    max_relative = Fraction(0)
    for i,j in combinations(range(m),2):
        sep = [prime for prime,x,y in zip(p,rows[i],rows[j]) if x!=y]
        assert len(sep)==m
        d=prod(sep)
        T=abs(b[i]-b[j])
        assert T not in labels and T>=2
        labels.add(T)
        # For h = T/E, h - h^3/6 <= sin h <= h.
        assert d >= E*E
        assert d*T*T < E*E*(T+eps)**2
        assert Fraction(T**3,6*E*E) < eps
        max_relative = max(max_relative,Fraction(d,E*E)-1)
    # The angular width is <= 2 range(b)/E, and parent norm N.
    assert (4*(max(b)-min(b)))**4 * N < E**4
    first,second=1,2
    d=prod(pv for pv,x,y in zip(p,rows[first],rows[second]) if x!=y)
    assert b[second]-b[first]==1036 and 1036%37==0
    assert d%37==22 and pow(22,18,37)==36
    print("Retained marks 1,2 fail necessary squareclass modulo 37: T=1036, d mod 37=22 is nonsquare.")
    squareclass_failures=0
    for i,j in combinations(range(1,m),2):
        t=b[j]-b[i]
        d=prod(pv for pv,x,y in zip(p,rows[i],rows[j]) if x!=y)
        squareclass_failures += any(pow(d%q,(q-1)//2,q)!=1 for q in range(3,t+1,2)
                                    if t%q==0 and prime(q))
    assert squareclass_failures==17
    print("At least one necessary odd squareclass fails on 17 of all 21 retained edges.")
    for q,L in ((3,4),(5,3),(7,2),(13,2),(17,2)):
        check_residue_tree(b,q,L,N)
    print("Paley M =",m,"prime range",p[0],p[-1],"margins squared < 1 + 1/",
          max_relative.denominator//max_relative.numerator)
    print(len(labels),"distinct determinants; max",max(labels),
          "; all raw cofactors within 1/64 by rational inequalities")
    print("coherent norm-N collision trees checked through 3^4, 5^3, 7^2, 13^2, 17^2")


if __name__ == '__main__':
    main()
