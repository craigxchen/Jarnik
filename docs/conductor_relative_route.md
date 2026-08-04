# Conductor-relative generalized-GCD route

## Purpose

The unrestricted uniform-in-`S` exceptional-degree statement is too broad. The
lattice-point problem supplies substantially more structure: every ratio comes
from one common Gaussian norm shell. This note records the corrected target,
the finite conductor-weight reductions that are rigorous, and the remaining
Diophantine gap.

## 1. Common conductor model

Fix a Gaussian integer `G`. A norm-one divisor ratio has the form

```text
u_A = A / conjugate(A),       A | G.
```

For a family arising from one circle, all ratios can be put in this form after
removing common unit and conjugation factors. The common conductor height is

```text
H_G = log |G|.
```

The square-root arc condition gives

```text
-log |u_A - 1| >= H_G / 2 - O_C(1).
```

The important distinction from arbitrary `S`-units is that the approximation
scale is measured against the total height of one common conductor, not against
the individual height of each pair.

## 2. Correct exceptional theorem

The route would be completed by the following conductor-relative statement.

### CRE(eta)

There are constants

```text
eta < 1/2,    D < infinity,    H0 < infinity
```

such that for every Gaussian integer `G` there is a nonzero polynomial

```text
F_G(X,Y) in Q(i)[X,Y],      deg F_G <= D,
```

and for all divisors `A,B | G`, outside `F_G = 0`,

```text
log GCD^+(u_A - 1, u_B - 1) < eta * log |G| + O(1),
```

provided `log |G| > H0`.

The polynomial may depend on `G`; its degree may not. This is enough for the
same all-pairs grid argument already formalized in
`UniformExceptionalReduction.lean`.

## 3. Why the unrestricted counterexample does not refute CRE

For

```text
x_a = (a+i)/(a-i),
y_b = (b+i)/(b-i),
```

the proximity is about one half of the sum of the two individual denominator
heights. A finite collection of unrelated pairs can be made into one `S`-unit
system by enlarging `S`, but the height of a common conductor containing every
denominator is approximately the sum of all those heights. Each individual
pair is then far below the `H_G/2` conductor-relative scale.

Thus the failure of uniformity for arbitrary `S`-units does not refute CRE.

## 4. Aggregating local conductor weights

At fixed finite Ru--Vojta level there is a finite set `T` of filtration types.
Every prime place `v` of the common conductor has a nonnegative normalized
conductor weight `w_v`, with total mass one. Grouping places of the same type,

```text
omega_t = sum_{v : tau(v)=t} w_v,
```

gives

```text
omega_t >= 0,       sum_t omega_t = 1.
```

For every candidate subspace `W`, its *slope score* factors through these
aggregate weights:

```text
score(W) = sum_t omega_t * c(t,W).
```

This statement is exact and is formalized in
`GaussianChain/ConductorWeights.lean`. It shows that the Harder--Narasimhan
slope optimization itself depends on a point of one fixed finite-dimensional
simplex, not on `|S|`.

## 5. Finite candidate signatures and near-maximal flats

The fixed adapted-section family gives a finite arrangement. The abstract
facts currently formalized are:

- finite subsets of a fixed family have only finitely many signatures
  (`FiniteFlatReduction.lean`);
- maximal-slope objects are closed under meet and join under the usual modular
  rank and supermodular degree hypotheses (`SlopeClosure.lean`);
- the set of flats whose score is within a fixed tolerance of maximal has only
  finitely many possible signatures (`NearMaximalReduction.lean`).

These facts remove the elementary chamber-wall problem: one can join all
`epsilon`-near-maximal flats and obtain a uniform slope gap outside that join.

## 6. Critical audit: slope aggregation is not the full theorem

The preceding finite reductions do **not** prove `CRE(eta)`.

The parametric Subspace Theorem contains two logically different pieces:

1. a canonical maximal-slope (Harder--Narasimhan) subspace;
2. exceptional high-parameter intervals outside that subspace.

Aggregation by conductor type controls the first piece. It does not control the
locations or the arithmetic spans of the interval-exception solutions. Those
subspaces can depend on the actual rational points and need not be flats of the
fixed section arrangement.

Equivalently, the equality

```text
sum_v w_v * c(tau(v),W) = sum_t omega_t * c(t,W)
```

aggregates the *slope score of a fixed subspace*. It does not aggregate the
pointwise twisted height

```text
sum_v c_{i,v} log |L_{i,v}(P)|_v,
```

because the values of the fixed forms at different prime places cannot be
recovered from the total weight carried by each filtration type.

Therefore compactness of the weight simplex, by itself, does not produce a
uniform height threshold and does not produce a bounded-degree exceptional
curve.

This is the precise gap in the proposed rescue proof.

## 7. What would actually complete the route

A valid completion needs an additional statement using the common conductor
beyond its aggregate slope weights. Two possible forms are:

### (a) Conductor-relative interval elimination

Show that the high-parameter interval exceptions in the fixed parametric
Subspace-Theorem system cannot occur when every finite valuation lies in one
common conductor box and the proximity is at least `(1/2-o(1)) H_G`.

### (b) A multipoint determinant theorem

For a fixed finite Ru--Vojta section space, show that any sufficiently large
collection of conductor-critical points has linearly dependent evaluation
vectors. This would put the whole collection on one fixed-degree hyperplane
section without enumerating the place-by-place adapted-basis assignments.

The obstruction in (b) is a mixed-basis assignment problem: a different point
may prefer a different adapted basis at the same prime. Any successful
multipoint determinant argument must recover the finite beta gain despite that
incompatibility.

Neither statement is presently proved in this branch.

## 8. Additional exact global constraint

For norm-one elements `x,y`, put

```text
r = (x-1)/(y-1).
```

Then exactly

```text
r / conjugate(r) = x / y.
```

This is formalized in `GaussianChain/ChordRatio.lean`. It implies that a
chord-ratio relation cannot have small arithmetic height when `x/y` has large
height. It is a useful constraint on any low-height exceptional curve, but it
does not by itself eliminate the critical Roth exponent.

The global affine-relation identity for an entire same-circle configuration is
formalized in `GaussianChain/AffineRelationGap.lean`.

## 9. Current formalization status

Compiler-checked before the latest two modules:

- explicit finite beta levels;
- numerical `1/2-eta` endgame;
- grid/fiber cardinality bound;
- finite signature count;
- maximal-slope meet/join closure;
- aggregation of arbitrary local weights into finitely many conductor types;
- near-maximal finite-signature reduction.

The latest branch also contains the chord-ratio and affine-relation modules;
CI should be consulted for their current compiler status.

Not formalized, and not yet proved:

- the conductor-relative interval-elimination theorem;
- a mixed-basis multipoint determinant theorem;
- `CRE(eta)` itself.

The conductor-relative route remains viable, but the missing input is genuinely
Diophantine. It is not a consequence of finite polyhedral chamber analysis
alone.
