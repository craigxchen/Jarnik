"""Exact D=8 function-field and integer near-square certificates."""
from collections import Counter, defaultdict
from itertools import combinations, product
from math import gcd, isqrt, lcm, prod

ROOTS = ((-12,4),(-8,6),(-3,16),(0,10),(0,12),(3,16),(8,6),(12,4))
SELECT = (0, 2, 4, 5, 6)
LOCAL_SAMPLES = {2:2, 3:0, 5:10, 13:5, 37:12, 41:17, 53:0, 509:238}


def moments(root):
    x, y = root
    return (y, x*y, 3*x*x*y-y**3)


def directions(roots):
    points = list(map(moments,roots))
    pair_buckets = defaultdict(list)
    for i,j in combinations(range(len(points)),2):
        key=tuple(a+b for a,b in zip(points[i],points[j]))
        pair_buckets[key].append((i,j))
    return [tuple((k,1) for k in ij)+tuple((k,-1) for k in rs)
            for entries in pair_buckets.values() for ij,rs in combinations(entries,2)]


def midpoint_count(roots):
    points = list(map(moments,roots))
    return sum(all(2*a == b+c for a,b,c in zip(center,left,right))
               for center in points for left,right in combinations(points,2))


def cluster_census(roots, signs, selected=SELECT):
    moments_by_row = defaultdict(list)
    moments_list = list(map(moments,roots))
    for mask in range(2**len(roots)):
        key = tuple(sum(v[j] for k,v in enumerate(moments_list) if mask >> k & 1)
                    for j in range(3))
        moments_by_row[key].append(mask)
    support = [sum(1 << k for k,_ in word) for word in signs]
    assert len(set(support)) == len(support)
    max_selected=0
    for rows in moments_by_row.values():
        edge_labels = {i for i in selected
                       if any(a^b == support[i] for a,b in combinations(rows,2))}
        max_selected=max(max_selected,len(edge_labels))
    return Counter(map(len,moments_by_row.values())), max_selected


def evaluate(roots, word, n):
    a,b=1,0
    for k,sign in word:
        x,y = roots[k]
        r,s = n-x,-sign*y
        a,b=a*r-b*s,a*s+b*r
    return a,b


def specialized(roots, words, n):
    parent = prod((n-x)**2 + y*y for x,y in roots)
    ds, ts = [], []
    for w in words:
        u,v = evaluate(roots,w,n)
        g = gcd(u,v)
        ds.append((u*u+v*v)//g**2)
        ts.append(abs(v)//g)
    return parent,ds,ts


def valuation(n,p):
    e=0
    while n%p==0:
        e+=1
        n//=p
    return e


def crt_step(residue, modulus, other, prime_power):
    lift = (other-residue)*pow(modulus,-1,prime_power) % prime_power
    return residue+modulus*lift,modulus*prime_power


def second_example():
    roots=((-23,4),(-11,10),(-8,16),(2,16),(5,10),(5,12),(13,6),(17,4))
    indices=(0,1,2,3,5,7,8)
    words=directions(roots)
    hist,covered=cluster_census(roots,words,indices)
    assert len(words)==9 and midpoint_count(roots)==2 and hist=={1:104,2:32,3:16,5:8} and covered==6
    words=[words[i] for i in indices]
    b=[abs(evaluate(roots,w,0)[1]) for w in words]
    samples={2:3,5:314,41:12,53:16,349:131,541:242,953:544}
    L,residue,period=1,0,1
    table=[]
    for p,n in samples.items():
        parent,ds,_=specialized(roots,words,n)
        k=1+valuation(parent,p)
        e=valuation(parent//lcm(*ds),p)
        if p==2:
            assert all(d%2 for d in ds)
        L*=p**e
        table.append((p,n,k,e))
        residue,period=crt_step(residue,period,n,p**k)
    assert L==125119585606720000000000 and L>256*max(b)**4
    parent,ds,ts=specialized(roots,words,residue)
    N=lcm(*ds)
    assert N%2 and parent%(L*N)==0 and len(set(ds))==7
    lower_log=(N.bit_length()-1)//2
    for d,t in zip(ds,ts):
        s=isqrt(d)
        assert d-s*s==t*t and gcd(s,t)==1 and d*d>=256*t**4*N
        assert d*d <= N*lower_log**10
    print('Second example: 9 ternary directions, fibre census',sorted(hist.items()),
          'seven selected, max covered',covered)
    print('Second local table (p,n,k,e)',table,'L',L)
    print('Second parent has',len(str(N)),'digits; offsets',ts)


def main():
    words = directions(ROOTS)
    assert len(words)==7 and midpoint_count(ROOTS)==2
    hist,covered=cluster_census(ROOTS,words)
    assert hist == {1:120,2:32,3:8,4:12} and covered==3
    print('Pair-collision directions:',words)
    print('Cluster-size histogram:',sorted(hist.items()),'maximum selected labels:',covered)
    vectors = [tuple(dict(word).get(i,0) for i in range(8)) for word in words]
    triangles = [(i,j,k,s,t) for i,j,k in combinations(SELECT,3)
                 for s,t in product((-1,1),repeat=2)
                 if all(a+s*b+t*c == 0 for a,b,c in zip(vectors[i],vectors[j],vectors[k]))]
    assert triangles == [(0,2,6,-1,1)]
    print('Signed triangles on selected edges:', triangles)
    selected = [words[j] for j in SELECT]
    b = [abs(evaluate(ROOTS,w,0)[1]) for w in selected]
    assert b==[4800,9600,13200,13200,4800]
    L,residue,period=1,0,1
    for p,n0 in LOCAL_SAMPLES.items():
        parent,ds,_ = specialized(ROOTS,selected,n0)
        lc = lcm(*ds)
        loss = valuation(parent//lc,p)
        mod = p**(valuation(parent,p)+1)
        L*=p**loss
        if p==2:
            assert all(d%2 for d in ds)
        residue,period=crt_step(residue,period,n0,mod)
        print('prime, sample, gain exponent, modulus exponent:',p,n0,loss,valuation(parent,p)+1)
    print('Certified contraction:',L,'residue:',residue,'period:',period)
    assert L==306439728192000000000 and L>256*max(b)**4
    n = residue
    parent,ds,ts=specialized(ROOTS,selected,n)
    N = lcm(*ds)
    assert parent%(N*L)==0 and N%2 and len(set(ds))==5
    for d,t in zip(ds,ts):
        s=isqrt(d)
        assert d-s*s==t*t and t>0 and gcd(s,t)==1 and d*d >= 256*t**4*N
    # An explicit upper window at B=5, using floor(log N) >= floor(log2 N)/2.
    lower_log = (N.bit_length()-1)//2
    assert all(d*d <= N*lower_log**10 for d in ds)
    print('Example parent:',N)
    print('Its five (d,t) pairs:',list(zip(ds,ts)))
    print('Sizes in digits:',len(str(n)),len(str(N)))


if __name__=='__main__':
    main()
    second_example()
