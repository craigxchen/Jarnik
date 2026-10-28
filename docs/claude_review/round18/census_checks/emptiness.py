"""Theorem E (emptiness below an explicit scale) -- exact numbers.

Input facts (proved in census.md §2):
 (a) at one prime p, the visible blocks T with p | n_T lie in two nested chains on disjoint row sets,
     hence there are at most k-3 of them (also true for the contact blocks of any Phi: a tree on k
     leaves has at most k-3 internal edges);
 (b) hence the j-th smallest visible block satisfies n_(j) >= p_{ceil(j/(k-3))} whenever all n_T >= 2.
Outputs, for k = 5..8:
  V = #visible blocks, r = ceil(V/(k-3)) (min #primes), P = p_r (some block is >= P),
  S = min possible sum_T log n_T, and
  D = min possible max_{i<j} log G_ij, G_ij = prod_{T contains i,j} n_T (a lower bound for
      log|det(P_i,P_j)| for some pair), via rearrangement over the weights C(|T|,2)."""
from math import comb, log, ceil
from sympy import prime
for k in range(5, 9):
    m = k - 1
    V = 2 ** m - m - 2   # visible splits (ptolemy.md §1 prints 2^m-m-1: typo)
    r = ceil(V / (k - 3))
    lows = [log(prime(ceil(j / (k - 3)))) for j in range(1, V + 1)]  # lower bounds for log n_(j)
    S = sum(lows)
    # weights: block T of size s contributes log n_T to C(s,2) pairs; total over pairs
    weights = sorted([comb(s, 2) for s in range(2, m) for _ in range(comb(m, s))], reverse=True)
    tot_pairs = sum(w * l for w, l in zip(weights, lows))  # smallest blocks on heaviest weights
    D = tot_pairs / comb(m, 2)
    b = 2 ** (m - 2) - 1   # visible blocks containing a given pair
    print(f"k={k}: V={V}, max blocks/prime={k-3}, min #primes r={r}, p_r={prime(r)}, "
          f"min sum log n_T={S:.2f} (avg {S/V:.3f}), "
          f"some pair has log G_ij >= {D:.2f}; balanced: log G_ij >= {b}*(1-eta)*w for all pairs")
    for eta in (0.0, 0.1, 0.25, 0.5):
        wmin = log(prime(r)) / (1 + eta)
        wavg = (S / V) / (1 + eta)
        w0 = max(wmin, wavg)
        print(f"    eta={eta}: hypothesis class empty for w < {w0:.3f} (e^w < {2.718281828**w0:.1f});"
              f" any realisation has log|det| >= {b*(1-eta)*w0:.1f} for every pair (|det|>=10^{b*(1-eta)*w0/log(10):.1f})")
