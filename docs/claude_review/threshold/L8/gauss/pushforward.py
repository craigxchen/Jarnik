from fractions import Fraction as Fr
from itertools import combinations
from gammahat import gamma, splits
for k in range(4, 10):
    labels, g = gamma(k)
    circ = list(range(1, k + 1)); full = frozenset(labels)
    def gD(S):
        S = frozenset(S); return g[S] if S in g else g[full - S]
    # pushforward to Mbar_{0,k} (forget P,Q): pi^* D_S = D_S + D_{S+P} + D_{S+Q} + D_{S+P+Q}, S a side of a split of [k]
    vals = set()
    for r in range(2, k - 1):
        for S in combinations(circ, r):
            S = set(S)
            vals.add(gD(S) + gD(S | {'P'}) + gD(S | {'Q'}) + gD(S | {'P', 'Q'}))
    # canonical class K = sum_S (s(n-s)/(n-1) - 2) D_S, n = k+2
    n = k + 2
    K = sum((Fr(len(S) * (n - len(S)), n - 1) - 2) * v for S, v in g.items())
    print('k=%d: pi_*gamma-hat . D_S over all boundary S of Mbar_{0,k}: %s (2^{1-k} = %s);  K.gamma-hat = %s' % (k, sorted(vals), Fr(1, 2 ** (k - 1)), K))
