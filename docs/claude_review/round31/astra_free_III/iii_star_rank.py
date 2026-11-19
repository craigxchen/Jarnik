"""Exact tests of shared-atom rank and parity signatures on the Paley double core."""
from itertools import combinations, permutations
from math import comb
from sympy import Matrix
from ii_paley import double_core


def field_basis(rows, prime):
    """Return the number of pivot rows over a prime field."""
    pivots = {}
    for v in rows:
        work=[x % prime for x in v]
        for i, x in enumerate(work):
            if x == 0:
                continue
            if i in pivots:
                w=pivots[i]
                work[i:] = [(a-x*b) % prime for a,b in zip(work[i:],w[i:])]
            else:
                u=pow(x,-1,prime)
                pivots[i]=[y*u % prime for y in work]
                break
    return len(pivots)


def verify(q):
    cols = sum(double_core(q), [])
    rows = list(zip(*cols))
    masks = [sum(x<<j for j,x in enumerate(row)) for row in rows]
    m = len(masks)
    labels = [i^j for i,j in combinations(masks,2)]
    assert 0 not in labels and len(set(labels)) == comb(m,2)
    def syndrome(i,j,k):
        # At each split prime the singleton valuations of a,b,c in II.1.
        a = (i^j) & (i^k)
        b = (i^j) & (j^k)
        c = (i^k) & (j^k)
        sig = (a^b,b^c,c^a)
        assert sig == (i^j,j^k,k^i)
        return sig
    signatures = {syndrome(i,j,k) for i,j,k in permutations(masks,3)}
    assert len(signatures) == m*(m-1)*(m-2)
    rank = Matrix([tuple(x-y for x,y in zip(v,rows[0])) for v in rows[1:]]).rank()
    assert rank == m-1
    # Degree-one vector signatures of five oriented marks, as in the proof.
    # One column per unprojected prime exponent for each of the four ratios.
    def phase_signature(t):
        return tuple(x-y for a in t[1:] for x,y in zip(rows[a],rows[t[0]]))
    base = phase_signature(range(5))
    def affine_differences():
        for t in permutations(range(m),5):
            yield (x-y for x,y in zip(phase_signature(t),base))
    # Theoretical Q-rank is at most 4(m-1); in the smaller example,
    # rank mod 101 attains it, certifying equality in characteristic 0.
    if m == 8:
        measured = field_basis(affine_differences(),101)
        assert measured == 4*(m-1)
    else:
        measured = None
    print('m',m,'edges',len(set(labels)),'ordered triangle signatures',len(signatures),
          'affine anchored-star rank',measured)


if __name__=='__main__':
    verify(7)
    verify(11)
