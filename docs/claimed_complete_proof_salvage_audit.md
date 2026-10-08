# Salvage audit of the earlier claimed complete proof

## Purpose

This note reconstructs the earlier claimed proof from the committed uniform-GCD route and separates the arguments that remain valid from the assertions that were never justified.

The conclusion is that the proof did not fail at the circle reduction, the generalized-GCD inequality, the finite beta calculation, or the final grid argument. It failed in the passage from the finite Ru--Vojta inequality to a uniform exceptional curve.

## 1. The proof skeleton

The claimed argument had the following form.

1. Normalize a cluster of Gaussian lattice points by one base point, obtaining S-units `u_j = z_j/z_0`.
2. Use the short arc to prove a generalized-GCD lower bound

   `log GCD^+(u_i-1,u_j-1) >= (1/2) log R - O_C(1)`.

3. Use the common norm to prove

   `h([1:u_i:u_j]) <= log R`.

4. Choose a finite Ru--Vojta level whose beta ratio is greater than one.
5. Invoke a Subspace-Theorem conclusion to place every violating pair `(u_i,u_j)` on an exceptional curve of degree bounded independently of `S`.
6. Apply the grid lemma to bound the size of the original cluster.

Only Step 5 was not proved.

## 2. Fully salvageable pieces

### 2.1 One common finite place set

For a fixed common norm, every quotient `u_j=z_j/z_0` is an S-unit for the same set S consisting of the archimedean place and the finite places above the rational primes dividing the norm.

This is exact and does not introduce a dependence problem by itself.

### 2.2 Archimedean generalized-GCD lower bound

For points on an arc of length at most `C sqrt R`,

`|u_j-1| <= C R^(-1/2)`.

The archimedean local contribution at `[1:1:1]` is therefore at least

`(1/2) log R - log C`.

The nonarchimedean positive-GCD contributions are nonnegative. Hence

`log GCD^+(u_i-1,u_j-1) >= (1/2) log R - log C`.

This part is valid.

### 2.3 Joint-height upper bound

At a split prime `p=pi bar(pi)`, if the common norm exponent is E and the three Gaussian integers have pi-adic exponents `a_0,a_i,a_j`, the combined contribution of the two conjugate places is

`(1/2)(max{a_0,a_i,a_j}-min{a_0,a_i,a_j}) log p`,

which is at most `(E/2) log p`.

Summing gives

`h([1:u_i:u_j]) <= log R`.

This is the crucial common-conductor estimate and remains valid.

### 2.4 The coefficient gap

Any uniform generalized-GCD theorem of the form

`log GCD^+(x-1,y-1) < eta h([1:x:y])`

outside a uniformly bounded-degree exceptional curve, with `eta<1/2`, contradicts the two previous estimates for sufficiently large R.

The numerical endpoint gap is therefore real and survives unchanged.

### 2.5 The grid argument

If all ordered off-diagonal pairs from `A x A` lie on one nonzero polynomial of total degree at most D, then

`|A|(|A|-1) <= D|A|`,

and hence `|A| <= D+1`.

This argument is elementary, correct, and already formalized conditionally in Lean.

### 2.6 The finite beta calculation

At finite approximation level `(N,k)=(100,247)`, the committed calculation gives

`beta = 45201/44900 > 1`.

Thus there is a genuine positive approximation margin at one fixed finite level. No limiting beta argument is required.

This calculation remains valid.

## 3. Abstract lemmas that remain true but do not finish the application

### 3.1 Supermodularity

For a genuinely fixed finite family of weighted flags, the degree function

`deg(U)=c dim U + sum_F a_F dim(U cap F)`, with `a_F>=0`,

is supermodular. Consequently maximal-slope subspaces are closed under sums and nonzero intersections, and there is a unique largest maximal-slope subspace.

This is a correct abstract linear-algebra theorem.

### 3.2 Finiteness of canonical Evertse--Ferretti subspaces

The naive claim that every canonical subspace is an intersection of kernels from the original finite form arrangement is false. However, Evertse--Ferretti give a uniform height bound for the canonical subspace. Since the field is fixed, Northcott finiteness on the Grassmannian implies that only finitely many canonical subspaces occur as the normalized weights vary.

Thus canonical-subspace finiteness survives, but from height boundedness rather than arrangement closure.

## 4. The exact failed chain

### 4.1 Unsupported finiteness of adapted flags

The earlier note said that the finite Ru--Vojta construction uses only finitely many adapted bases, hence finitely many flags.

That is not established and is generally false as stated: a fixed filtration may admit infinitely many adapted bases. One can sometimes replace bases by discrete filtration data, but that replacement has to be proved for the exact local maximum used in the inequality.

### 4.2 Pointwise maxima were treated as one global system

The finite Ru--Vojta inequality contains local maxima over adapted bases. The maximizing basis can depend on both the place and the point.

The claimed proof silently replaced these pointwise choices by one normalized twisted-height system. That is the same incompatibility later exposed by the determinant counterexample: optimization performed separately for each point does not automatically assemble into one global linear system.

### 4.3 Canonical subspace versus all solutions

Even when a canonical maximal-slope subspace is uniformly controlled, a quantitative Subspace Theorem may cover solutions using additional auxiliary subspaces, or a qualitative theorem may leave an S-dependent finite exceptional set outside the canonical closed set.

Neither remainder can be ignored for a fixed conductor, because the conductor grid is itself finite.

### 4.4 Uniform threshold was not checked

The final claim required a height threshold independent of S, the support of the conductor, and the normalized local weights. No line-by-line derivation of such a threshold was supplied.

Published quantitative bounds in the directly applicable formulation may contain `|S|`, while qualitative uniformity of the closed exceptional set does not control the finite S-dependent remainder.

## 5. What the earlier proof actually establishes

The valid result is the following conditional theorem.

> If, for some `eta<1/2`, the generalized-GCD inequality for the fixed divisor configuration has an exceptional plane curve whose total degree and height threshold are bounded independently of S, then the square-root arc conjecture follows.

In addition, the finite beta calculation proves that the fixed divisor configuration has more than enough local approximation weight to make such a theorem applicable if the uniform exceptional-set passage can be justified.

## 6. Best route forward from the salvaged proof

The most economical remaining task is not to invent a new angle, determinant, or box argument. It is to prove one precise bridge:

> Convert the finite-level Ru--Vojta max-over-adapted-bases inequality into a fixed-family Subspace-Theorem statement whose complete exceptional set, including any finite remainder or auxiliary covering spaces, has complexity and threshold independent of S.

There are two plausible ways this could succeed.

1. Show that the local maximum depends only on a finite collection of filtration or matroid types, despite the infinitude of adapted bases, and apply one uniform parametric theorem to each type.
2. Use the Ru--Vojta theorem's uniform closed exceptional set, then separately prove that the S-dependent finite remainder is absent above a common-conductor height threshold.

Everything else in the earlier proof can then be reused without modification.

## 7. Status summary

### Retain

- common-S reduction;
- short-arc generalized-GCD lower bound;
- common-conductor joint-height upper bound;
- endpoint coefficient gap;
- finite beta calculation;
- grid/fiber bound;
- abstract supermodularity;
- canonical-subspace finiteness via height bounds.

### Withdraw

- finite adapted-basis family as stated;
- arrangement-kernel description of the canonical subspace;
- conversion of pointwise maxima into one global twisted-height system;
- assertion that the canonical subspace alone captures every sufficiently large solution;
- uniform height threshold without a dependency audit.

The earlier proof therefore contains a nearly complete reduction, not a complete proof. Its missing content is one uniform exceptional-set bridge, concentrated in a single stage.