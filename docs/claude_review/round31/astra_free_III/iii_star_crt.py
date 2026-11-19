"""Exact maxima for the squarefree horizontal integral-star model, for fixed parents."""
from itertools import combinations
from math import gcd, isqrt, prod


def fixed_parent(primes):
    n=prod(primes)
    bound=isqrt(isqrt(n))  # floor(n**(1/4)) by nested square roots.
    # n is squarefree; L>=sqrt(n) iff it contains some divisor d with d*d>=n.
    minimal=[]
    for k in range(len(primes)+1):
        for subset in combinations(primes,k):
            d=prod(subset)
            if d*d>=n and all((d//p)**2<n for p in subset):
                minimal.append(d)
    good=sorted({i for d in minimal for i in range(0,n,d)})
    residues=good+[x+n for x in good if x<=bound]
    end=best=0
    for start in range(len(good)):
        while end<len(residues) and residues[end]-residues[start]<=bound:
            end+=1
        if end-start>best:
            best=end-start
            witness=residues[start:end]
    assert all(gcd(x,n)**2>=n for x in witness)
    assert witness[-1]-witness[0]<=bound
    return n,bound,best,witness


if __name__=='__main__':
    primes=[5,13,17,29,37,41,53,61]
    for k in range(3,9):
        n,h,best,w=fixed_parent(primes[:k]); pos=[x-w[0] for x in w]
        print(n,h,best,w[0],pos)
