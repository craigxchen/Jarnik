# Conductor-relative generalized-GCD route

## Purpose

The unrestricted uniform-in-`S` exceptional-degree statement is too broad.  The
lattice-point problem supplies substantially more structure: every ratio comes
from one common Gaussian norm shell.  This note records the corrected target
and the part of the Ru--Vojta optimization that can be made independent of the
number of places without assuming any new Diophantine theorem.

## 1. Common conductor model

Fix a Gaussian integer `G`.  A norm-one divisor ratio has the form

```text
u_A = A / conjugate(A),       A | G.
```

For a family arising from one circle, all ratios can be put in this form after
removing common unit and conjugation factors.  The common conductor height is

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

The polynomial may depend on `G`; its degree may not.

This is enough for the same all-pairs grid argument already formalized in
`UniformExceptionalReduction.lean`.

## 3. Why the unrestricted counterexample does not refute CRE

For

```text
x_a = (a+i)/(a-i),
y_b = (b+i)/(b-i),
```

the proximity is about one half of the sum of the two individual denominator
heights.  A finite collection of unrelated pairs can be made into one `S`-unit
system by enlarging `S`, but the height of a common conductor containing every
denominator is approximately the sum of all those heights.  Each individual
pair is then far below the `H_G/2` conductor-relative scale.

Thus the failure of uniformity for arbitrary `S`-units does not address CRE.

## 4. Aggregating local conductor weights

At fixed finite Ru--Vojta level there is a finite set `T` of filtration types.
Every prime place `v` of the common conductor has:

```text
- a type tau(v) in T;
- a nonnegative normalized conductor weight w_v;
- sum_v w_v = 1.
```

For any candidate subspace `W`, its local filtration score has the form

```text
score(W) = sum_v w_v * c(tau(v), W),
```

where `c(t,W)` depends only on the fixed section family.

Define aggregate weights

```text
omega_t = sum_{v : tau(v)=t} w_v.
```

Then exactly

```text
score(W) = sum_{t in T} omega_t * c(t,W),
```

with

```text
omega_t >= 0,       sum_t omega_t = 1.
```

This is formalized in `GaussianChain/ConductorWeights.lean`.

Consequently, the Harder--Narasimhan optimization does not depend on `|S|`.
It depends only on a point of the fixed simplex of aggregate type weights.

## 5. Finite candidate flats

The adapted section family is fixed and finite.  Candidate exceptional flats
are intersections of kernels of subsets of those sections, so only finitely
many flats can occur.  The abstract finite-signature count is formalized in
`FiniteFlatReduction.lean`.

For each aggregate weight vector, the maximal-slope flats are closed under meet
and join when rank is modular and filtration degree is supermodular.  Hence
there is a canonical largest maximal-slope flat.  This closure mechanism is
formalized in `SlopeClosure.lean`.

Therefore the possible canonical destabilizing flats belong to one finite
family independent of the number of prime places and independent of the
conductor.

## 6. What remains

The remaining input is narrower than a uniform quantitative Subspace Theorem:

> For the fixed finite section family and every aggregate conductor-weight
> vector, points satisfying the strict Ru--Vojta inequality at sufficiently
> large conductor height lie in the associated canonical maximal-slope flat,
> with a height threshold uniform over the compact weight simplex.

There are two distinct tasks.

1. **Pointwise parametric statement.**  Prove the large-height containment for
   one fixed aggregate weight vector.
2. **Uniformity over weights.**  Show the height threshold can be chosen
   uniformly over the simplex.  Since all slope functions are linear and only
   finitely many candidate flats occur, this should reduce to finitely many
   rational polyhedral chambers plus a positive margin from the critical
   slope.

The first task is the genuinely Diophantine part.  The second is finite
polyhedral analysis and should be formalizable once the first statement is
available with explicit dependence on the slope gap.

## 7. Current formalization status

Compiler-checked:

- explicit finite beta level;
- numerical `1/2 - eta` endgame;
- grid/fiber cardinality bound;
- finite signature count;
- maximal-slope meet/join closure;
- aggregation of arbitrary local weights into finitely many conductor types.

Not yet formalized:

- Gaussian conductor normalization for an actual circle;
- local Weil functions and `GCD^+`;
- the pointwise parametric Subspace-Theorem containment;
- uniform dependence of its height threshold on the aggregate weight vector.

The route remains viable, but the remaining theorem must be stated and proved
in conductor-relative form rather than as an unrestricted uniform-in-`S`
exceptional-degree theorem.
