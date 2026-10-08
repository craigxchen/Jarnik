# Finite adapted bases versus the quantitative exceptional set

## Status

The finite-type reduction for the Autissier--Ru--Vojta filtrations is valid: after fixing the approximation level, one may choose a finite family of representative adapted bases independent of the place set `S`.

This note audits whether that fact alone closes the uniform exceptional-degree argument. It does not.

## 1. The repaired finite local maximum

Let

```text
B_1, ..., B_M
```

be the fixed finite family of representative bases in the section space `V`. The local expression obtained from the filtration argument is of the form

```text
max_{1 <= j <= M} sum_{s in B_j} lambda_{s,v}(P).
```

Equivalently, after putting

```text
P_j = product_{s in B_j} s,
```

all `P_j` being sections of the same fixed tensor power, this becomes

```text
max_j lambda_{P_j,v}(P)
```

up to the fixed normalization of local Weil functions.

Thus the infinite adapted-basis issue is removed completely.

## 2. Why choosing the maximizing basis separately does not finish the proof

For each point `P` and each place `v`, choose an index `j(P,v)` attaining the maximum. This produces a placewise system

```text
v |-> B_{j(P,v)}.
```

For one fixed system, the quantitative Subspace Theorem has constants controlled by the ambient dimension, approximation margin, heights, and number of distinct forms. In the Evertse--Ferretti parametric formulation the normalization and large-parameter threshold do not depend directly on the number of places.

However, the exceptional subspaces supplied by the quantitative theorem depend on the chosen placewise system. As `P` varies, the map `v |-> j(P,v)` may vary. A naive union ranges over as many as

```text
M^|S|
```

systems. Uniformly bounding the number of subspaces for each individual system does not bound the total union over all systems.

## 3. Product sections do not remove the quantifier by themselves

The product-section identity

```text
max_j sum_{s in B_j} lambda_{s,v}(P)
  = max_j lambda_{P_j,v}(P)
```

is useful, but ordinary linear-form versions of the quantitative theorem treat one fixed independent system at each place. Replacing the local maximum by a choice of one `P_j` still creates the same point-dependent assignment.

A theorem stated directly for

```text
sum_v max_j lambda_{P_j,v}(P)
```

with a complete exceptional set of uniformly bounded total degree would close the argument. The qualitative multidivisor Evertse--Ferretti/Ru--Vojta theorem has precisely this max-type shape, but it does not give the needed uniform quantitative exceptional degree and height threshold.

## 4. What the non-absolute quantitative theorem does and does not supply

For `K`-rational points, Evertse--Ferretti's improved quantitative theorem gives bounds for a fixed system which depend on the number of distinct forms and their heights, rather than explicitly on `|S|`. This is important and salvages the threshold for any fixed assignment.

It does not identify one common exceptional union for all point-dependent assignments from the finite family. The auxiliary subspaces in the quantitative covering theorem need not be canonical and are not proved to range over a fixed finite list as the assignment varies.

The canonical Harder--Narasimhan subspace does range over a fixed finite family, but the interval theorem permits solutions outside that canonical subspace in finitely many height intervals. A common conductor supplies one height scale, which may lie in one of those intervals. Therefore the interval theorem alone does not remove the auxiliary subspaces.

## 5. Exact remaining theorem

The repaired proof is complete conditional on the following statement.

> **Uniform finite-maximum theorem.** Fix a number field, a projective variety, a finite family of equal-degree sections, and a positive approximation margin. Then the sufficiently large rational points satisfying a global inequality whose local terms are the maxima over that fixed section family lie in a proper closed set whose total degree and height threshold depend only on the fixed family and the margin, not on the finite place set.

This statement is strictly stronger than applying the ordinary quantitative theorem separately to each maximizing assignment. It is exactly the quantitative multidivisor/max-over-bases theorem required by the common-conductor argument.

## 6. Consequence for the endpoint proof

Everything else in the endpoint proof remains valid:

1. the common-`S` reduction;
2. the `1/2` generalized-GCD lower bound;
3. the joint-height upper bound;
4. the finite beta gain greater than one;
5. the finite representative basis family;
6. the coefficient contradiction;
7. the grid/fiber bound after obtaining a uniformly bounded-degree exceptional set.

The only unresolved implication is the uniform finite-maximum theorem above. Finite filtration type repairs the infinite-basis problem, but does not by itself repair the point-dependent global assignment problem.
