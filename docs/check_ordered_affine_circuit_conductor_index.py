"""Exact checks for banded affine indices and their nested conductor divisor."""

from collections import defaultdict
from fractions import Fraction as F
from functools import cmp_to_key, reduce
from itertools import combinations, product
from math import gcd, isqrt, prod
from random import Random

from check_joint_affine_relation_lattice import (
    determinant, gconj, gdiv_exact, ggcd, gmul, gnorm, gsub,
    minor_polynomials, minors, polynomial_gcd, primitive_points, qtrim,
)
from check_paley_compatible_phase_coordinate_gap import profile


def det(a):
    a = [list(map(F,row)) for row in a]
    value = F(1)
    for j in range(len(a)):
        pivot = next((i for i in range(j,len(a)) if a[i][j]),None)
        if pivot is None:
            return 0
        if pivot != j:
            a[j],a[pivot] = a[pivot],a[j]
            value = -value
        value *= a[j][j]
        scale = a[j][j]
        a[j] = [x/scale for x in a[j]]
        for i in range(j+1,len(a)):
            scale = a[i][j]
            a[i] = [x-scale*y for x,y in zip(a[i],a[j])]
    assert value.denominator == 1
    return value.numerator


def band_index(points, all_minors=False):
    m = len(points)
    ds = minors(points)
    assert all(ds.values())
    g = reduce(gcd,map(abs,ds.values()))
    gs,columns = [],[]
    for j in range(m-3):
        window = list(range(j,j+4))
        local = reduce(gcd,[abs(ds[t]) for t in combinations(window,3)])
        gs.append(local)
        c = [0]*m
        for r,i in enumerate(window):
            c[i] = (-1)**r*ds[tuple(k for k in window if k!=i)]//local
        assert sum(c)==0
        assert all(sum(c[i]*points[i][a] for i in range(m))==0 for a in (0,1))
        columns.append(c)
    numerator = g*prod(abs(ds[(j,j+1,j+2)]) for j in range(1,m-3))
    denominator = prod(gs)
    assert numerator%denominator == 0
    index = numerator//denominator
    assert prod(abs(c[j]) for j,c in enumerate(columns)) == index*abs(ds[(m-3,m-2,m-1)])//g
    assert prod(abs(c[j+3]) for j,c in enumerate(columns)) == index*abs(ds[(0,1,2)])//g
    if all_minors:
        for rows in combinations(range(m),m-3):
            complement = tuple(i for i in range(m) if i not in rows)
            actual = det([[c[i] for c in columns] for i in rows])
            assert abs(actual) == index*abs(ds[complement])//g
    return index,g,gs


def valuation(n,p):
    assert n
    n,answer = abs(n),0
    while n%p == 0:
        n//=p
        answer+=1
    return answer


def gvaluation(z,pi):
    assert z != (0,0)
    p,answer = gnorm(pi),0
    while True:
        quotient = gmul(z,gconj(pi))
        if any(x%p for x in quotient):
            return answer
        z = tuple(x//p for x in quotient)
        answer += 1


def gpower(z,n):
    return reduce(gmul,[z]*n,(1,0))


def nested_formula():
    for e in range(5):
        for t in product(range(e+1),repeat=5):
            u,v = min(t[1:4]),max(t[1:4])
            charge = max(0,min(t[0],t[4])-v)+max(0,u-max(t[0],t[4]))
            ranges = lambda x:max(x)-min(x)
            assert charge == ranges(t[:4])+ranges(t[1:])-ranges(t)-ranges(t[1:4])
            count = sum({i for i in range(5) if t[i]>=j} in ({0,4},{1,2,3})
                        for j in range(1,e+1))
            assert charge == count
    rng = Random(39577)
    tested = charged = excess = 0
    for pi in ((2,1),(3,2)):
        p = gnorm(pi)
        other = ((4,1),(5,2),(6,1))
        options = []
        for signs in product((0,1),repeat=3):
            base = reduce(gmul,[gconj(z) if s else z for z,s in zip(other,signs)],(1,0))
            options.extend(gmul(base,u) for u in ((1,0),(-1,0),(0,1),(0,-1)))
        for e in range(1,5):
            for trial in range(100):
                t = [rng.randrange(e+1) for _ in range(5)]
                if trial%2 == 0:
                    t = [e,0,0,0,e]
                points = [gmul(gmul(gpower(pi,a),gpower(gconj(pi),e-a)),z)
                          for a,z in zip(t,rng.sample(options,5))]
                if len(set(points))<5:
                    continue
                n = gnorm(points[0])
                assert all(gnorm(z)==n for z in points)
                index,_,_ = band_index(points)
                actual_t = [gvaluation(z,pi) for z in points]
                assert actual_t == t
                u,v = min(t[1:4]),max(t[1:4])
                h = max(0,min(t[0],t[4])-v)+max(0,u-max(t[0],t[4]))
                value = valuation(index,p)
                assert value>=h
                common = reduce(ggcd,points)
                primitive_n = n//gnorm(common)
                assert gcd(index,primitive_n)%(p**h)==0
                if h:
                    pairs = {(i,j):gvaluation(gsub(points[i],points[j]),pi)-t[i]
                             for i,j in combinations(range(5),2) if t[i]==t[j]}
                    if u!=v:
                        correction = sum(x for (i,j),x in pairs.items() if 1<=i<j<=3)
                    else:
                        triple = [pairs[key] for key in ((1,2),(1,3),(2,3))]
                        a = min(triple)
                        delta = min(a,pairs[(0,4)]) if t[0]==t[4] else 0
                        correction = sum(triple)-2*a+delta
                    assert value == h+correction
                    charged += 1
                    excess += correction>0
                tested += 1
    assert charged and excess
    print(f"nested valuations: {tested} actual tuples, {charged} charged, {excess} with excess")


def split_generators(count):
    answer = []
    p = 5
    while len(answer)<count:
        if p%4==1 and all(p%d for d in range(2,isqrt(p)+1)):
            for x in range(1,isqrt(p)+1):
                y = isqrt(p-x*x)
                if y*y == p-x*x:
                    answer.append((x,y))
                    break
        p += 1
    return answer


def singleton_and_banded():
    h,s,labels,_ = profile(11)
    generators = split_generators(len(s[0]))
    points = [reduce(gmul,[z if sign==1 else gconj(z) for z,sign in zip(generators,row)],(1,0))
              for row in s]
    ds = minors(points)
    for subset in combinations(range(12),4):
        local = reduce(gcd,[abs(ds[t]) for t in combinations(subset,3)])
        for pos,i in enumerate(subset):
            coefficient = (-1)**pos*ds[tuple(j for j in subset if j!=i)]//local
            original = sum(all(h[j][a]==-h[i][a] for j in subset if j!=i) for a in range(1,12))
            isolated = [a for a in range(len(generators))
                        if all(s[j][a]==-s[i][a] for j in subset if j!=i)]
            assert len(isolated)>=5*original-4
            divisor = prod(gnorm(generators[a]) for a in isolated)
            assert coefficient%divisor == 0
    for m in range(4,9):
        band_index(points[:m],all_minors=True)
    print("PASS: singleton rational-prime charges, bn_i-h count, and all banded minors")


def endpoint_families():
    polynomials = {key:qtrim(value) for key,value in minor_polynomials().items()}
    degrees = []
    for omitted in range(5):
        rows = [i for i in range(5) if i!=omitted]
        local = reduce(polynomial_gcd,[polynomials[t] for t in combinations(rows,3)])
        degrees.append(len(local)-1)
    assert degrees == [0,2,2,2,2]
    assert len(reduce(polynomial_gcd,polynomials.values()))==1
    for n,expected in ((12,(2937,89)),(80,(261393,89)),(100,(2937,1))):
        points = [primitive_points(n)[j] for j in (3,0,4,1,2)]
        index,_,_ = band_index(points,all_minors=True)
        assert (index,gcd(index,gnorm(points[0]))) == expected
    c,e,d,f = [404,428,464,484,524],[720,-720,-720,-720,720],[85,61,25,5,-35],[1164,396,324,-276,-1236]
    for t in (600,1200,1800,6000):
        points = [((-t**6-45*t**4-c[i]*t*t+e[i])//720,
                   (-t**5-d[i]*t**3-f[i]*t)//720) for i in range(5)]
        assert all(points[i][0]>points[i+1][0] and points[i][1]<points[i+1][1] for i in range(4))
        assert all(x<0 and y<0 for x,y in points)
        assert gnorm(reduce(ggcd,points))==1
        index,g,gs = band_index(points,all_minors=True)
        a = 1 if (t//6)%3==0 else 3
        assert g==t//15 and gs==[a*t//15,t//15]
        assert index == (21 if a==1 else 7)*(t//6)**2
        assert gcd(index,gnorm(points[0]))==1
    print("PASS: ordered quartic local-gcd degrees and cubic external-index growth")


def polar_compare(z,w):
    half = lambda a:0 if a[1]>0 or (a[1]==0 and a[0]>=0) else 1
    if half(z)!=half(w):
        return -1 if half(z)<half(w) else 1
    cross = determinant(z,w)
    return -1 if cross>0 else (1 if cross<0 else 0)


def finite_scan():
    circles = defaultdict(list)
    limit = 20000
    for x in range(-isqrt(limit),isqrt(limit)+1):
        for y in range(-isqrt(limit),isqrt(limit)+1):
            n = x*x+y*y
            if 0<n<=limit:
                circles[n].append((x,y))
    windows = overlap = maximum = 0
    for n,points in circles.items():
        if len(points)<5:
            continue
        points.sort(key=cmp_to_key(polar_compare))
        for i in range(len(points)):
            p = [points[(i+j)%len(points)] for j in range(5)]
            if determinant(p[0],p[-1])<=0 or gnorm(reduce(ggcd,p))!=1:
                continue
            index,_,_ = band_index(p)
            conductor = gcd(index,n)
            windows += 1
            overlap += conductor>1
            maximum = max(maximum,conductor)
    fixtures = [([(4,-33),(9,-32),(12,-31),(23,-24),(24,-23)],1105,5,F(4),(13,21,34)),
                ([(-110,-35),(-109,-38),(-98,-61),(-94,-67),(-86,-77)],13325,13,F(23,5),(1,9,82))]
    for p,n,expected,bound,metric in fixtures:
        assert all(gnorm(z)==n for z in p) and gnorm(reduce(ggcd,p))==1
        index,_,_ = band_index(p)
        assert index==expected==gcd(index,n)
        chord2 = gnorm(gsub(p[-1],p[0]))
        v = F(chord2,4*n)
        assert 0<v<1
        # arcsin(u)/u <= 1+u²/[6(1-u²)²], by the derivative bound.
        arc_fourth_upper = F(chord2*chord2,n)*(1+v/(6*(1-v)**2))**4
        assert arc_fourth_upper < bound**4
        e1,e2 = gsub(p[1],p[0]),gsub(p[2],p[0])
        gram = (gnorm(e1),sum(x*y for x,y in zip(e1,e2)),gnorm(e2))
        content = reduce(gcd,gram)
        a,b,c = tuple(x//content for x in gram)
        assert (a,b,c)==metric and a*c-b*b==1 and a>0
        ds = minors(p)
        local = reduce(gcd,(abs(ds[t]) for t in combinations(range(4),3)))
        circuit = [(-1)**i*ds[tuple(j for j in range(4) if j!=i)]//local
                   for i in range(4)]
        c0,c1,c2,c3 = circuit
        assert sum(circuit)==0 and c0==expected
        assert gcd(c1*c2,expected)==1
        assert (a+c-2*b)*c1*c2 == -c0*(a*c1+c*c2)
        assert valuation(a+c-2*b,expected)==1
    print("PASS: positive conic metric identity; conductor absorbed once in a primitive chord")
    print(f"primitive circle scan N<=20000: {windows} windows, {overlap} with overlap, max gcd={maximum}")


if __name__ == "__main__":
    nested_formula()
    singleton_and_banded()
    endpoint_families()
    finite_scan()
    print("PASS: conductor divisor and exact local-index checks; no uniform upper bound claimed")
