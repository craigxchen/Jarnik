"""Formal exponent tables e_j(v) = 1 + chi(eps_j lambda_j(v)) over V = i*F_q^k; common minimum g."""
import itertools, sys
from collections import Counter
from f9lib import Fq

def value_table(F, k, forms):
    """Counter over v in F_q^k of the tuple (chi(lambda_j(v)))_j"""
    q = F.q
    tab = Counter()
    for x in itertools.product(range(q), repeat=k):
        row = []
        for lam in forms:
            s = 0
            for c, xb in zip(lam, x):
                s = F.addt[s][F.smul(c, xb)]
            row.append(F.chi[s])
        tab[tuple(row)] += 1
    return tab

def best_g(tab, m):
    best = []
    for eps in itertools.product((1, -1), repeat=m):
        g = sum(f*min(1 + e*c for e, c in zip(eps, row)) for row, f in tab.items())
        best.append((g, eps))
    best.sort(reverse=True)
    return best

if __name__ == '__main__':
    F = Fq(3)
    tab = value_table(F, 2, [(1,0),(0,1),(1,1),(1,2)])
    print('q=27 k=2 four forms: n=729 h=27, best g, eps:', best_g(tab, 4)[:3])
    print('all-positive g:', [b for b in best_g(tab,4) if b[1]==(1,1,1,1)])
    tab5 = value_table(F, 3, [(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,0,1)])
    b5 = best_g(tab5, 5)
    print('q=27 k=3 five forms: n=19683 h=729, top g:', b5[:4])
    print('all-positive g:', [b for b in b5 if b[1]==(1,)*5])
    # each F_j has degree n: check
    for j in range(5):
        assert sum(f*(1+row[j]) for row, f in tab5.items()) == 27**3
