"""Exact cofactor checks on sheared rational C10 and K2,5."""
from itertools import combinations
from math import gcd, prod


def mul(z,w):
    return (z[0]*w[0]-z[1]*w[1], z[0]*w[1]+z[1]*w[0])


def power(z,n):
    v=(1,0)
    for _ in range(n):
        v=mul(v,z)
    return v


def primitive(z):
    d=gcd(*z)
    return (z[0]//d,z[1]//d)


def edge(z,w):
    X=z[0]*w[0]+z[1]*w[1]
    Y=z[0]*w[1]-z[1]*w[0]
    g=gcd(X,Y)
    return abs(Y)//g,(X*X+Y*Y)//(g*g),abs(Y)


def check(name,points,marked):
    points=list(map(primitive,points))
    assert len(set(points))==len(points)
    pairs=list(combinations(range(len(points)),2))
    assert all(edge(points[i],points[j])[2] for i,j in pairs)
    L=2*prod(edge(points[i],points[j])[2]**2 for i,j in pairs)
    pts=[(x, y+L*x) for x,y in points]
    assert all(y for x,y in pts)
    pts=[(x+L*y,y) for x,y in pts]
    ds=set()
    for i,j in pairs:
        t,d,_=edge(pts[i],pts[j])
        assert t==edge(points[i],points[j])[0]
        if (i,j) in marked:
            assert t==1
        ds.add(d)
    assert len(ds)==len(pairs)
    print(name, len(points),"marks,",len(marked),"unit edges and",len(ds),"distinct divisor labels")


def main():
    vertices=[(a,0) for a in range(5)]+[(a,1) for a in range(4,-1,-1)]
    stars=[mul(power((2,1),a),power((4,1),b)) for a,b in vertices]
    c10={(i,i+1) for i in range(9)}|{(0,9)}
    check('C10',stars,c10)
    # n=211 has n²+4=44525=5²*13*137. Five distinct choices below sqrt.
    n=211
    assert n%4==3 and n*n+4==5**2*13*137
    div=[1,5,13,25,65]
    a=[(n+c)//2 for c in div]
    bs=[(n+(n*n+4)//c)//2 for c in div]
    assert all(u%2==v%2==0 and n*(u+v)==2*(u*v-1) for u,v in zip(a,bs))
    points=[(1,0),(n,2)]+[(u,1) for u in a]
    check('K2,5',points,{(j,i) for j in (0,1) for i in range(2,7)})


if __name__=='__main__':
    main()
