"""Exact checks of the primitive near-critical parity-code allocation.

This checks allocation norms, not the angular condition of a short arc.
"""
from itertools import combinations
from math import log, prod


def main():
    primes = [5,13,17,29,37,41,53,61]
    for r,t in [(2,1000), (3,10000)]:
        m = 2**r
        base_primes, q = primes[:m-1], primes[m-1]
        powers = [2*max(0,round((t/log(p)-1)/2))+1 for p in base_primes]
        rows = [[s*(bin(x & ell).count('1') % 2)
                 for ell,s in zip(range(1,m),powers)] + [2*(2**x-1)]
                for x in range(m)]
        caps = [max(row[j] for row in rows) for j in range(m)]
        assert all(min(row[j] for row in rows) == 0 for j in range(m))
        assert caps == powers+[2*(2**(m-1)-1)]
        all_primes = base_primes+[q]
        N = prod(p**e for p,e in zip(all_primes,caps))
        pairs = {}
        classes = {}
        for i,j in combinations(range(m),2):
            differences = [abs(a-b) for a,b in zip(rows[i],rows[j])]
            norm = prod(p**e for p,e in zip(all_primes,differences))
            label = tuple(e % 2 for e in differences)
            assert any(label) and label[-1] == 0 and N % norm == 0
            assert norm**2 > N and norm**m < N**(m//2+1)
            assert norm not in pairs.values()
            pairs[i,j] = norm
            classes.setdefault(label, []).append((i,j))
        assert len(classes) == m-1
        for matching in classes.values():
            assert len(matching) == m//2
            assert len({i for pair in matching for i in pair}) == m
        print(f'PASS: M={m}, {len(pairs)} distinct nonsquare norms in the exact near-critical interval; {len(classes)} repeated matching classes.')
    print('SCOPE: Gaussian arguments are uncontrolled; this is not an endpoint construction.')


if __name__ == '__main__':
    main()
