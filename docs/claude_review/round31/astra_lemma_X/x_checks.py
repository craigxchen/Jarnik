"""Exact checks for Part X; standard library and the Part IX family generators."""
from bisect import bisect_left
from collections import Counter
from itertools import combinations, combinations_with_replacement
from math import gcd, isqrt, prod
from unittest.mock import patch

import ix_checks as ix


def heron(a, b, c):
    return (a+b+c) * (-a+b+c) * (a-b+c) * (a+b-c)


def triangle_test(moduli, edges):
    """Return g if all local three-row square tests pass."""
    g = 1
    for p in moduli:
        residues = sorted((d*d) % p for d in edges)
        if residues == [0, 0, 0]:
            g *= p
        elif residues[0] != 0 or residues[1] == 0 or residues[1] != residues[2]:
            return None
    return g


def exhaustive_triangles():
    counts = Counter()
    for q in range(2, 301):
        pp = list(ix.prime_powers(q))
        if any(p != power for p, power in pp):
            continue
        counts["moduli"] += 1
        primes = [p for p, _ in pp]
        for edges in combinations_with_replacement(range(1, isqrt(4*q) + 1), 3):
            g = triangle_test(primes, edges)
            if g is None:
                continue
            counts["local"] += 1
            H = heron(*edges)
            divisor = (q*g)**2
            assert H % divisor == 0
            quotient = H // divisor
            assert quotient % 4 in (0, 3)
            if not (q//g) % 2:
                assert quotient % 4 == 0
            if H:
                counts["nonflat"] += 1
                assert max(edges)**2 >= q*g
            if max(edges)**2 <= q:
                counts["unit cap"] += 1
                assert H == 0
    print("small squarefree triangle census", dict(counts))


def pte_tests():
    # Feed the unmodified Part IX integer examples into the Part X local tests.
    def check(label, marks, factorization):
        primes = list(factorization)
        q = prod(primes)
        for p in primes:
            local = {}
            for i,j in combinations(range(len(marks)),2):
                d = marks[j]-marks[i]
                local.setdefault(bool(d % p),set()).add(d*d % p)
            assert all(len(residues) == 1 for residues in local.values())
        positions = sorted(marks)
        count = 0
        for a,b,c in combinations(positions,3):
            edges = (b-a, c-b, c-a)
            g = triangle_test(primes,edges)
            assert g == gcd(q,gcd(*edges))
            assert heron(*edges) == 0
            count += 1
        print(label, len(marks), "vertices", count, "triples")

    with patch.object(ix, "check_integer_family", check):
        ix.explicit_families()
        for k in (4,5,8):
            weights = [1]*((k+1)//2) + [-1]*(k//2)
            ix.pte(f"{k} PTE nodes", [(w,w*j) for j,w in enumerate(weights)],1,1009)
        vectors = [(-3,-3),(-3,-2),(0,-1),(1,2),(1,3),(2,0),(2,1)]
        ix.pte("two-moment PTE family",[(a+7*b,b) for a,b in vectors],2,1009)


def paley_test(m=43, lower=40000):
    primes = [p for p in ix.primes_upto(2*lower) if p >= lower and p % 4 == 1]
    target = (3*lower//2)**2
    used = set()
    factors = []
    for p in primes:
        if p in used:
            continue
        start = bisect_left(primes,(target+p-1)//p)
        for r in primes[start:]:
            if (p*r-target)*m*8 > target:
                break
            if r != p and r not in used:
                used.update((p,r))
                factors.append((p,r))
                break
        if len(factors) == m:
            break
    assert len(factors) == m and len(used) == 2*m
    Q = [p*r for p,r in factors]
    q = prod(Q)
    residues = {i*i % m for i in range(1,m)}
    rows = [sum(1 << v for v in range(m) if (v-i) % m in residues) for i in range(m)]
    labels = [(lower//(4*m))*i + (i == m-1) for i in range(m)]
    edges = {}
    agreements = {}
    separated_squares = dict((p,set()) for p in used)
    for i,j in combinations(range(m),2):
        separated = rows[i] ^ rows[j]
        assert separated.bit_count() == (m+1)//2
        h = prod(Q[v] for v in range(m) if not ((separated >> v) & 1))
        edge = h*(labels[j]-labels[i])
        assert edge**2 <= q and gcd(edge,q) == h
        edges[i,j], agreements[i,j] = edge,h
        for v,(p,r) in enumerate(factors):
            if (separated >> v) & 1:
                separated_squares[p].add(edge**2 % p)
                separated_squares[r].add(edge**2 % r)
    counts = Counter()
    for i,j,k in combinations(range(m),3):
        g = gcd(agreements[i,j],agreements[j,k])
        common = (~((rows[i]^rows[j]) | (rows[i]^rows[k]))) & ((1 << m)-1)
        assert common.bit_count() == (m-3)//4
        assert g == prod(Q[v] for v in range(m) if (common >> v) & 1)
        H = heron(edges[i,j],edges[j,k],edges[i,k])
        counts["triples"] += 1
        counts["Heron nonzero"] += bool(H)
        counts["divisibility failures"] += bool(H % (q*g)**2)
    local_counts = [len(values) for values in separated_squares.values()]
    print("IX.4 Paley on",m,"vertices",dict(counts),"distinct squares per prime",
          min(local_counts),max(local_counts))
    assert counts["triples"] == counts["divisibility failures"]
    assert min(local_counts) > 1


def small_det(matrix):
    if not matrix:
        return 1
    return sum((-1)**j * value * small_det([row[:j]+row[j+1:] for row in matrix[1:]])
               for j,value in enumerate(matrix[0]))


def erasure_test(m=11):
    residues = {x*x % m for x in range(1,m)}
    S = [[1 if (v-i) % m in residues else -1 for v in range(m)] for i in range(m)]
    sets = circuits = 0
    for k in range(3):
        for omitted in combinations(range(m),k):
            T = [[1]*m] + [S[i] for i in omitted]
            r = len(T)
            def minor(cols):
                return small_det([[row[j] for j in cols] for row in T])
            pivots = next(cols for cols in combinations(range(m),r) if minor(cols))
            sets += 1
            for j in range(m):
                if j in pivots:
                    continue
                coeff = [0]*m
                coeff[j] = minor(pivots)
                for i,v in enumerate(pivots):
                    coeff[v] = -minor(pivots[:i]+(j,)+pivots[i+1:])
                assert all(sum(b*x for b,x in zip(row,coeff))==0 for row in T)
                image = [sum(b*x for b,x in zip(row,coeff)) for row in S]
                assert any(image) and sum(image)==0 and all(image[i]==0 for i in omitted)
                for v in range(m):
                    if coeff[v]==0:
                        assert sum(image[i]*S[i][v] for i in range(m))==0
                circuits += 1
    print("Paley erasure check at",m,":",sets,"row sets,",circuits,"annihilators")


def excluded_small_cases():
    # Cycle edges alternate 17,29; opposite vertices are at distance 46.
    q = 17*23*29
    edges = {(0,1):17,(0,2):46,(0,3):29,(1,2):29,(1,3):46,(2,3):17}
    assert max(edges.values())**2 < q
    for p in (17,23,29):
        inside = [(i,j) for (i,j),d in edges.items() if d % p == 0]
        assert len(inside) == 2 and set(inside[0]).isdisjoint(inside[1])
        assert len({d*d % p for d in edges.values() if d % p}) == 1
    for i,j,k in combinations(range(4),3):
        assert heron(edges[i,j], edges[j,k], edges[i,k]) == 0
    # On the line the diameter between 0 and 2 would put 1 at 17 and 3 at 29.
    assert abs(edges[0,1]-edges[0,3]) != edges[1,3]
    # A new assignment on five marks of the referee's quadratic PTE template.
    a,b = 83,59
    u,v = a*(a+b),b*(a+2*b)
    positions = [0,u,v,u+v,(a+b)*(a+2*b)]
    primes = [59,67,71,83]
    q5 = prod(primes)
    original = {(i,j):abs(positions[j]-positions[i]) for i,j in combinations(range(5),2)}
    altered = dict(original)
    altered[0,3] = abs(u-v)
    assert len({y % q5 for y in positions}) == 5
    assert max(altered.values())**2 < 36*q5
    for p in primes:
        occupied = {x % p for x in positions}
        assert len(occupied) == 2
        r,s = occupied
        for (i,j),edge in altered.items():
            assert (edge % p == 0) == (positions[i] % p == positions[j] % p)
            if edge % p:
                assert edge**2 % p == (r-s)**2 % p
    H = heron(altered[0,3],altered[0,4],altered[3,4])
    assert H == -1152*q5*q5
    print("modified five-mark PTE: modulus",q5,"Heron quotient",H//(q5*q5))
    marks = [(4,7),(8,-1),(8,1)]
    assert all(x*x+y*y == 65 for x,y in marks)
    squares = [(x-X)**2+(y-Y)**2 for (x,y),(X,Y) in combinations(marks,2)]
    assert squares == [80,52,4]
    print("four-cycle on modulus",q,"Gaussian chord squares",squares)


if __name__ == "__main__":
    exhaustive_triangles()
    pte_tests()
    paley_test()
    erasure_test()
    excluded_small_cases()
