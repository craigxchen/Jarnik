# Projection, junta, and weighted sunflower route

## Status

This note audits the projection/junta idea and extracts a stronger surviving mechanism based on weighted sunflowers of difference supports.

The uniform endpoint theorem is not proved here.

## 1. Projection fibers and the basic recursion

Write the exponent box as

```text
B = product_j {0,...,e_j}
```

and the angular map as

```text
Phi(a) = sum_j (2 a_j - e_j) theta_j mod 2 pi.
```

For an endpoint interval I of width

```text
delta = C exp(-H/4),
H = sum_j e_j log p_j,
```

let

```text
A = Phi^{-1}(I).
```

Fix a coordinate j and slice by a_j=t. Each slice is again a fiber of the lower-dimensional map obtained by deleting coordinate j, but with a translated target interval.

This gives an exact recursion. However, large fiber size alone does not imply one large slice unless e_j+1 is bounded, and iterating the largest slice can lose a factor product_j(e_j+1), which is precisely the full divisor count.

## 2. Why a pure weighted-junta theorem is false

A large fiber of a product map need not be controlled by a bounded coordinate set. Abstract examples include parity fibers and linear-code fibers in binary cubes. Every proper coordinate projection may remain large while the set depends essentially on all coordinates.

The same obstruction can be approximated on the circle whenever many angle coordinates lie near a common rational grid. Thus any valid junta theorem must include a quantitative anti-resonance hypothesis on the actual Gaussian-prime angles.

Projection alone therefore does not escape the global arithmetic obstruction; it merely exposes it as rational-grid resonance.

## 3. Collision dichotomy under projection

For the projection deleting coordinate j, a collision means two cluster points agree in every coordinate except j. If their j-th exponents differ by d != 0, then endpoint concentration gives

```text
|| 2 d theta_j / (2 pi) || <= O(delta).
```

Equivalently,

```text
(pi_j / conjugate(pi_j))^d
```

is O(delta)-close to 1.

The elementary Gaussian separation bound then forces the block height

```text
d log p_j
```

to be at least H/2-O_C(1).

Hence only boundedly many coordinate blocks can support genuine one-coordinate collisions in a minimal endpoint counterexample.

So after deleting a bounded heavy set, almost every coordinate projection is injective on A. The remaining set is therefore a high-dimensional graph-like code rather than a branching product set.

This is a useful exact reduction, but injective projections do not preserve angular concentration after deleting the coordinate, so they do not give direct induction.

## 4. Difference supports

Fix a root point a_0 in A. After binary threshold-layer expansion, associate to each a in A the weighted support

```text
S(a) = {layers on which a differs from a_0}.
```

Each layer has weight log p_j. The conductor distance between two points is the weighted symmetric difference

```text
D(a,b) = weight(S(a) symmetric_difference S(b)).
```

Endpoint separation gives

```text
D(a,b) >= H/2-O_C(1)
```

for every distinct pair.

## 5. Exact weighted sunflower lemma for endpoint clusters

Suppose t supports form an exact sunflower:

```text
S_i = C union P_i,
```

where the petals P_i are pairwise disjoint.

The common core cancels in pairwise symmetric differences, so

```text
D(a_i,a_j) = weight(P_i)+weight(P_j).
```

Thus every pair satisfies

```text
weight(P_i)+weight(P_j) >= H/2-O_C(1).
```

At most one petal can have weight below H/4-O_C(1). All remaining petals have weight at least H/4-O_C(1).

Because the petals are disjoint and their total weight is at most H, there can be at most four such large petals, plus at most one exceptional small petal.

Therefore, for sufficiently large H:

> An endpoint cluster cannot contain a weighted sunflower of six difference supports.

This is a genuine global consequence of the circle geometry. It is not a determinant or pairwise-cut bound: the contradiction comes from simultaneous cancellation of one common core across many points and disjointness of all petals.

## 6. Approximate sunflowers

The same argument is stable. Suppose

```text
S_i = C union P_i union E_i
```

with petals P_i pairwise disjoint and the total weighted error in each pairwise symmetric difference at most eta H.

Then

```text
weight(P_i)+weight(P_j) >= (1/2-O(eta))H-O_C(1).
```

For eta sufficiently small, five or six disjoint petals remain impossible by summing their weights.

Thus a weighted approximate sunflower theorem with a fixed number of petals would prove the endpoint bound.

## 7. Why ordinary sunflower theory is insufficient

Classical sunflower bounds depend on the uniform set size or ambient combinatorial complexity. Here the difference supports can contain arbitrarily many threshold layers, and an unbounded cluster may grow much more slowly than the classical sunflower threshold.

Balanced Hadamard-type codes also show that large pairwise weighted distance does not force a fixed-size exact or approximate sunflower in an arbitrary weighted set system.

Therefore the missing input cannot be a purely combinatorial sunflower theorem. It must use the fact that all support sets arise from actual Gaussian divisor angles lying in one exponentially small interval.

## 8. New target: arithmetic sunflower extraction

The concrete theorem suggested by this route is:

> **Arithmetic weighted sunflower theorem.** For every fixed endpoint constant C, there exists M(C) such that any endpoint cluster with more than M(C) points contains six points whose threshold-layer difference supports relative to one root form an o(H)-approximate weighted sunflower.

The exact weighted sunflower lemma would then give a contradiction.

A plausible extraction mechanism is to combine:

1. projection collisions, which identify and remove boundedly many heavy resonant blocks;
2. injective-projection structure on the diffuse remainder;
3. Ramsey or regularity on the coordinate patterns of a fixed number of selected points;
4. the common angular constraint, used to rule out persistent non-sunflower overlap patterns by Gaussian algebraic separation.

The fourth step is essential. Without it, abstract balanced codes are counterexamples.

## 9. Verdict

The pure projection/junta route fails.

The surviving theorem is sharper:

```text
six-petal weighted sunflower => contradiction.
```

The unresolved problem is arithmetic extraction of such an approximate sunflower from a hypothetical unbounded endpoint cluster.

This is not obviously equivalent to the original conjecture: it identifies one finite global pattern whose presence is impossible and asks for a finite-pattern extraction theorem from the actual Gaussian realization.
