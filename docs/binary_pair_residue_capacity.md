# Binary pair-residue capacity

This note checks whether the binary valuation

    k_ij = v_2(b_ij)

can turn the good-pair count into a stronger height estimate. Assume the
primitive odd-N Gaussian setup, so the pair parameters satisfy the binary
ultrametric rule:

* if two adjacent k's are unequal, the third is their minimum;
* if two adjacent k's are equal, the third is strictly larger.

This is the same valuation pattern as differences of distinct 2-adic
residues.

## Exact layer-cake minimum when all pairs are included

For M distinct residues, let S_2(M) be the minimum of

    sum_{i<j} v_2(b_ij).

At level r, pairs with v_2(b_ij) >= r lie in a common residue class modulo
2^r. Distributing M objects as evenly as possible among the 2^r classes
minimizes that level's contribution. Writing

    q_r = floor(M / 2^r),    s_r = M - 2^r q_r,

the exact minimum is

    S_2(M) = sum_{r>=1} [
      s_r * choose(q_r+1,2) + (2^r-s_r) * choose(q_r,2)
    ].                                                    (1)

The residues `0,1,...,M-1` attain balanced occupancy at every depth
simultaneously, proving that this lower bound is the exact minimum in the
abstract binary-tree model. No simultaneous endpoint-circle realization
is asserted. The sum is finite because its summand vanishes once 2^r > M.
For M=2^t,

    S_2(2^t) = (M^2-M)/2 - M*t/2.

Thus the minimum product of the binary parts is exponential in M^2, but its
minimal average binary valuation is bounded by a constant. It does not force an average
logarithm of M.

## The good-pair quota does not access this minimum

The argument supplies only

    Q >= M(M-2)/4

good pairs with b_ij <= H, where H=N^(1/(2M)). This quota is below the
maximum number of cross edges in a parity bipartition:

    floor(M^2/4) >= M(M-2)/4.                            (2)

Every cross-parity edge has v_2(b_ij)=0. Therefore the binary ultrametric
rules are compatible, at the valuation level, with all required good edges
having v_2(b)=0: put the vertices into two parity classes and take the cross
edges as the selected set. The same-side edges carry positive valuation and
are irrelevant to the quota.

Consequently no positive lower bound on the sum of v_2(b_ij) over good edges,
and no growing average binary valuation, follows from the stated number of
good pairs. The exact full-pair formula (1) cannot be applied after deleting
the same-parity edges.

## What the binary layer does prove

The usual occupancy argument gives a capacity statement: among a family
whose relevant b's all lie in [1,H], the binary tree has at most

    M <= 2H

leaves/vertices. The exact bound is `2^(floor(log_2 H)+1)<=2H`. This is useful when
all pair parameters are bounded. It cannot be applied when only Q of the
binom(M,2) pairs are bounded.

At level r, the strongest generic statement for a selected edge set would
involve the number of selected edges internal to residue classes modulo 2^r.
The quota (2) permits selecting only cross-class edges at level one, so those
internal counts can all be zero for every selected edge. No layer-cake
accumulation is forced.

## Consequence for extensions

To obtain a stronger average height, one needs additional information that
forces good edges inside parity classes, or uses odd-prime valuations and
controls the odd parts of b_ij. Binary valuations alone do not do so. The
proposed factorial/product route therefore yields no new bound at the current
Plotkin quota. This is an exact capacity obstruction, rather than a
counterexample to a theorem that uses further Gaussian or odd-prime
structure.
