"""Exact horizontal integer stars with partially prescribed Gaussian factors."""
from itertools import combinations
from math import isqrt, prod

from iii_graph_rounding import prime, skew_core


def root_of_minus_one(p):
    nonresidue = 2
    while pow(nonresidue, (p - 1) // 2, p) == 1:
        nonresidue += 1
    root = pow(nonresidue, (p - 1) // 4, p)
    assert (root * root + 1) % p == 0
    return root


def main():
    m = 19
    primes = []
    candidate = 2**40 + 1
    while len(primes) < 2*m:
        if prime(candidate):
            primes.append(candidate)
        candidate += 4
    rows = [[s == 1 for s in row] for row in skew_core(m)[1:]]
    kept = [2*j for j in (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 15)]
    roots = {j:root_of_minus_one(primes[j]) for j in kept}
    h = isqrt(primes[0]**m)
    g0 = primes[0]**((m-1)//2)
    labels = [2*m*m*i + i*i for i in range(m)]
    xs = []
    n = prod(primes)
    for row, label in zip(rows, labels):
        residue, period = 0, 2
        for j in kept:
            if row[j]:
                p = primes[j]
                residue += period * (((roots[j]-residue) * pow(period, -1, p)) % p)
                period *= p
        target = h + g0*label
        start = (target-residue)//period
        for k in range(2*len(kept)+2):
            x = residue + (start+k)*period
            norm = x*x + 1
            if all((x-roots[j])%primes[j] != 0 or row[j] for j in kept) and all(
                (x+roots[j])%primes[j] and norm%(primes[j]**2) for j in kept
            ):
                break
        else:
            raise AssertionError('No CRT representative')
        for j in kept:
            assert ((x-roots[j])%primes[j]==0) == row[j]
        assert x%2 == 0
        intended = prod(p for bit,p in zip(row,primes) if bit)
        assert 100*abs(norm-intended)<3*intended
        xs.append(x)
    width = max(xs)-min(xs)
    assert 4**4 * n * width**4 < min(xs)**8
    differences, counts, defects, worst = set(), [], 0, 0
    for i,j in combinations(range(m), 2):
        shared = [k for k in range(2*m) if rows[i][k] and rows[j][k]]
        selected = [k for k in shared if k in kept]
        assert len(shared)==9
        g = prod(primes[k] for k in shared)
        det = xs[i]-xs[j]
        delta = labels[i]-labels[j]
        assert det % prod(primes[k] for k in selected) == 0
        defects += (det%g != 0)
        err = abs(det-g*delta)
        assert 1000*err < g
        worst = max(worst, (10**9*err+g-1)//g)
        counts.append(len(selected))
        differences.add(delta)
    assert len(differences)==m*(m-1)//2 and min(counts)>0 and defects>0
    print('marks, pairs, exact global atoms, selected shared range:', m,len(differences),len(kept),min(counts),max(counts))
    print('uniform normalized error <=',worst,'/ 10^9; incomplete full divisibilities:',defects)
    print('all squared norm errors < 3%; doubled angular diameter < N^(-1/4)/2')


if __name__=='__main__':
    main()
