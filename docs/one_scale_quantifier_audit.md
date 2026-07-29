# Fixed-scale parametric theorem: quantifier audit

## Exact statement used

For a fixed projective variety `Y`, a fixed twisted-weight system `c`, and a fixed
error margin `delta`, the Evertse--Ferretti twisted-height theorem produces a
finite family of bounded-degree hypersurfaces `F_1(c),...,F_t(c)` such that for
every sufficiently large `Q`, one member of this family contains the entire
small twisted-height locus at that `Q`.

The degree bound and the lower bound for `Q` are uniform in `c` once the
normalizations on `c` are imposed.  The hypersurfaces themselves are not stated
to be uniform in `c`.

## Consequence for one conductor

A common Gaussian conductor gives one common scale `Q`.  Thus there is no need
to cover many height intervals.  This removes the scale-decomposition issue.

It does not by itself remove the weight-system issue.  In the Ru--Vojta
reduction the local weight tuple is selected from the local section evaluations
of the point.  Hence different ordered divisor pairs may produce different
weight systems `c` even at the same conductor and the same `Q`.

The remaining implication is therefore

```text
for every pair P there exists c(P) with P in Z(c(P),Q)
```

to

```text
all pairs P lie in one bounded-degree Z(Q).
```

This implication is not a formal consequence of the fixed-`c` parametric
theorem: a continuously varying family of bounded-degree hypersurfaces can have
Zariski-dense union.

## Precise remaining Gaussian statement

For the sextic section space `H^0(P^2,I_[1:1:1](6))`, prove that all twisted
weight systems realized by divisor pairs of one Gaussian conductor have their
fixed-`Q` exceptional hypersurface in a finite family independent of the
conductor.

It would suffice to prove either of the following.

1. Every realized system has the same canonical destabilizing flat after
   conjugate-place pairing.
2. The canonical flat depends only on one of finitely many global signatures
   whose number is independent of the number of prime factors of the conductor.
3. A single fixed twisted-height system dominates every realized system with a
   loss smaller than the available numerical margin below `1/2`.

The conjugate-place identity controls the total Hilbert weight but has not yet
been shown to imply any of these three stronger statements.
