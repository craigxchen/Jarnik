# Deep second pass on global methods for the endpoint arc problem

## Status

This note revisits several previously compressed directions and pushes each to its strongest plausible formulation. The uniform bound is not proved.

The common setup is the divisor-angle model

```text
N = product_j p_j^(e_j),
H = log N,
alpha(a) = sum_j (2a_j-e_j) theta_j mod 2pi,
0 <= a_j <= e_j.
```

A bad family would give sets of exponent vectors with cardinality tending to infinity whose angles lie in intervals of width `delta = O(exp(-H/4))`.

---

## 1. Adelic capacity and potential theory

### Strongest plausible formulation

For selected norm-one ratios `u_1,...,u_M`, form the adelic logarithmic energy

```text
E = - sum_(i<j) sum_v log |u_i-u_j|_v.
```

The product formula couples the archimedean clustering to finite-place divisibility exactly. At the complex place, an endpoint arc gives approximately

```text
- log |u_i-u_j|_infinity >= H/4 - O_C(1).
```

The finite places record common Gaussian divisors of pair differences.

### Exact endpoint computation

After expanding valuations into binary threshold layers, the total finite energy is a nonnegative weighted sum of cut sizes

```text
k(M-k).
```

Its maximum is `floor(M^2/4)` on balanced cuts. Summing the archimedean lower bounds over pairs gives the same main coefficient. Hence ordinary adelic logarithmic capacity is exactly critical and admits the same balanced-code extremizers as the determinant method.

### Higher capacity

One could replace pair energy by sectional capacity or multivariate transfinite diameter. This amounts to using many monomials and Fekete determinants. It is not automatically dead, but its first nontrivial specializations reproduce interpolation determinants already found to be critical. A viable capacity argument would need a local set at finite places smaller than the full conductor box, so that the global capacity product becomes strictly less than one. No such uniformly smaller finite adelic set has yet been identified.

### Verdict

Pairwise adelic capacity fails exactly. Higher sectional capacity remains a language for a possible nonlinear auxiliary-polynomial proof, but currently offers no independent gain.

---

## 2. Delsarte linear programming and list decoding

### Strongest plausible formulation

After threshold expansion, selected points form a weighted binary code. Pairwise conductor distances satisfy

```text
D(i,j) >= H/2 - O_C(1).
```

Classical Delsarte linear programming studies positive-definite kernels on the Hamming scheme. Here one can enlarge the scheme by adjoining the angular syndrome

```text
Phi(sigma) = sum_l sigma_l theta_l mod 2pi.
```

The desired object is a positive-semidefinite kernel depending jointly on weighted Hamming overlap and angular difference. If it were positive on the diagonal and uniformly negative for every distinct pair inside the tiny arc, the usual Gram-matrix argument would bound the code size.

### Obstruction from continuity

Any translation-invariant positive-definite kernel on the circle has nonnegative Fourier coefficients and is continuous. Since distinct point angles can approach each other at the endpoint scale, no fixed kernel can be positive at zero and uniformly negative on all nonzero differences in the arc. A kernel of bandwidth about `delta^(-1)` can oscillate fast enough, but its coefficient mass and arithmetic complexity become conductor-sized; this returns to the short-Fourier-sum problem.

### Higher semidefinite hierarchies

Schrijver/Terwilliger-style semidefinite bounds use triple distributions rather than pair distances. They could potentially exclude abstract Plotkin extremizers when combined with consistency of nested prime layers. However, arbitrary large balanced codes have admissible triple distributions at relative distance one half. The angular syndrome must therefore enter essentially. This produces a new, concrete optimization problem rather than an immediate theorem:

> classify positive-semidefinite kernels on the joint orientation-angle association scheme whose degree may grow like `N^(1/4)` but whose trace bound remains absolute.

No current construction achieves this.

### Verdict

Ordinary LP bounds fail at Plotkin equality. A joint angle-code SDP is genuinely global and not fully exhausted, but presently it is equivalent to finding a phase-sensitive high-bandwidth kernel with controlled arithmetic cost.

---

## 3. Geometry of numbers on the lifted difference lattice

### Setup

Choose one cluster point as origin. For every other point, choose a winding integer and form

```text
b_i = (a_i-a_0, n_i) in Z^(r+1)
```

with

```text
|2 theta dot (a_i-a_0) - 2pi n_i| <= O(delta).
```

Thus all lifted differences lie in a very thin slab around one real hyperplane.

### Successive-minima strategy

If the affine rank is `d`, take `d` independent lifted differences. Their exterior product is a nonzero integral `d`-vector, so its Euclidean norm is at least one. On the other hand, the thin-slab condition makes the component involving the normal direction smaller by a factor `delta^d`.

A contradiction would follow if the tangential size of the wedge could be bounded by `exp((d/4-epsilon)H)`. The natural conductor-box bound is instead roughly `exp(dH/2)`, and even using pairwise endpoint separation only pushes the support width down to the critical `H/2` per independent direction. The resulting inequality again lands at equality.

### Possible unused gain

The full wedge has many Pluecker coordinates sharing the same prime blocks. A compound-body estimate might exploit overlap of coordinate supports and show that `d` independent differences cannot each consume disjoint half-conductors. This would be genuinely global. The precise target is:

> bound the weighted volume of the lattice generated by clustered differences strictly below the product of their individual weighted widths.

Balanced random sign systems show this is false for abstract weighted boxes. Any proof must use Gaussian-prime angles in the covolume, not only exponent supports.

### Verdict

Naive Minkowski/successive-minima estimates are critical. A compound-lattice inequality coupling covolume to the actual angle coefficients remains a possible but unproved global mechanism.

---

## 4. Sum-product, incidences, and rational cross-ratios

### Exact circle structure

For four points on one circle, their complex cross-ratio is real. Since the points lie in `Q(i)`, the cross-ratio is therefore rational. This is exact higher-order information:

```text
cr(u_1,u_2,u_3,u_4) in Q.
```

Under tangent parametrization of the rational circle, this is equivalent to the usual rational cross-ratio of four rational parameters.

### Potential approach

A large short arc produces many rational cross-ratios formed from very small Gaussian differences. One could hope that shared conductor factors make these rationals have much smaller height than generic cross-ratios, forcing coincidences. Repeated coincidences would put the point set on a low-complexity projective configuration.

### Height audit

The common small angular scale cancels from the cross-ratio, but so do the common denominators. The remaining numerator and denominator are products of primitive chord factors. Their heights can still be of full order `H`; no uniform improvement follows from the arc condition alone. Rational sets in a tiny interval can also have many distinct cross-ratios of comparable height. Standard sum-product and Szemeredi-Trotter estimates regard such a sparse set as generic rather than contradictory.

### Surviving target

The only plausible version is an ideal-theoretic cross-ratio theorem:

> prove that cross-ratios of four common-norm Gaussian points have conductor content that forces height below the generic product-of-chords bound.

Current gcd identities do not provide that gain.

### Verdict

Generic sum-product/incidence theory is too coarse. Rational cross-ratios expose exact four-point structure but no strict height saving has been found.

---

## 5. Concentration compactness and profile decomposition

### Strongest formulation

Assume a minimal bad sequence. Order prime blocks by normalized conductor mass

```text
lambda_(n,j) = e_(n,j) log p_(n,j) / H_n.
```

After passing to a subsequence, one can extract:

1. finitely many heavy profiles with positive limiting mass;
2. a diffuse remainder whose maximal block mass tends to zero;
3. possible resonant angular profiles after rescaling modulo rational grids.

This mirrors profile decompositions in critical PDE: compact profiles plus a dispersive remainder.

### What it buys

Heavy profiles are compatible with existing descent. Therefore a genuine counterexample must survive in the diffuse remainder. A global theorem would then need to show that a diffuse product of prime-block measures cannot retain an unbounded rare fiber at scale `exp(-H/4)` unless a rational-grid profile emerges.

### Exact obstruction

Compactness alone admits nonstandard limits consisting of balanced independent sign coordinates and an infinitesimal linear functional. Such limits model all abstract endpoint extremizers. Arithmetic separation is lost unless the limiting profile retains quantitative rates relative to the conductor.

### Verdict

Concentration compactness is a useful organizational framework, not the missing estimate. It cleanly separates heavy descent profiles from a diffuse arithmetic resonance problem.

---

## 6. Approximate groups and stabilizers of rare fibers

### Strongest formulation

For partial convolutions `sigma`, define the concentration function

```text
Q_delta(sigma) = sup_x sigma([x-delta,x+delta]).
```

If convolution by a packet `nu` barely decreases `Q_delta`, then many translations in the support of `nu` must almost stabilize a near-maximizing interval or near-maximizing level set of `sigma`.

Iterating failure of flattening should produce an approximate stabilizer `K` on the circle. Approximate subgroups of the circle are controlled by finite cyclic subgroups or short arcs. Hence repeated non-flattening should force many packet angles near one common rational grid.

### Why this is stronger than Fourier resonance

The grid must stabilize the same rare fiber through many convolution stages. It is not merely a frequency at which several blocks have large Fourier coefficients. Persistence could amplify endpoint closeness into supercritical closeness.

### Main technical gap

A short interval has trivial exact stabilizer. Near-stabilization of a probability level set need not imply that the translating measure lies near a finite subgroup unless one has quantitative additive-combinatorial control of the level set. For a needle event of mass `M/R_N -> 0`, standard approximate-group theorems do not give uniform constants.

### Verdict

This remains one of the strongest live directions. The needed theorem is a worst-fiber analogue of inverse convolution flattening with constants independent of the rare-event mass fraction.

---

## 7. Inverse Littlewood-Offord for rare fibers

### Existing mechanism

Sharp inverse Littlewood-Offord theorems characterize coefficient multisets whose random signed sums have unusually large concentration probability: most coefficients lie in a low-rank generalized arithmetic progression of controlled volume.

### Translation to this problem

After threshold expansion, the full divisor measure is a signed sum of Gaussian-prime angles. An arc with `M` atoms has probability

```text
rho = M / R_N.
```

Applying an inverse theorem at radius `delta` should put most angular coefficients in a low-rank progression whose volume depends polynomially on `rho^(-1)`.

### Quantitative mismatch

Here `rho` may be much smaller than any fixed power of the number of coordinates, while `M` still tends to infinity. Existing inverse theorems then allow a progression of complexity comparable to the whole coefficient set. More importantly, they count coordinates uniformly, whereas one enormous prime carries much more conductor height than many small primes.

### Required weighted theorem

The exact target is:

> if `rho R_N = M` is large at radius `exp(-H/4)`, then either finitely many coordinates carry positive conductor mass, or a positive fraction of conductor mass—not merely coordinate count—lies near a bounded-rank progression whose complexity depends only on `M`.

This is substantially stronger than current inverse Littlewood-Offord statements, but it is the cleanest formulation of the rare-fiber problem.

### Verdict

This is the most direct surviving external paradigm. It does not yet solve the weighting or vanishing-density issues.

---

## 8. Algebraic unlikely intersections

### Strongest formulation

The exponent vectors define points in a high-dimensional torus, while endpoint concentration places their character values near one archimedean target. Repeated approximate character relations could indicate proximity to a proper subtorus. Quantitative unlikely-intersection or Bogomolov-type results might then force exact containment.

### Obstruction

The ambient dimension and defining primes vary with `N`, and the points are not small-height points in the usual normalized sense. Known theorems are strongest for fixed tori, fixed subvarieties, or exact intersections. Uniform quantitative control in growing dimension would be at least as difficult as the original problem.

### Verdict

Useful only after another method reduces to bounded rank or bounded dimension.

---

## Comparative conclusion

The deeper pass changes the ranking but not the theorem status.

### Still live

1. **Rare-fiber inverse Littlewood-Offord with conductor weights.**
2. **Approximate stabilizers from repeated failure of concentration-function flattening.**
3. **Compound geometry of numbers using actual Gaussian angle coefficients.**
4. **Joint angle-code semidefinite programming.**

### Organizational or bounded-rank endgames

- concentration compactness;
- unlikely intersections;
- rational cross-ratios.

### Precisely critical in their standard forms

- pairwise adelic capacity;
- ordinary Delsarte/Plotkin bounds;
- naive Minkowski and wedge estimates;
- generic sum-product/incidence theory.

The most promising synthesis is now:

```text
minimal bad family
 -> profile decomposition into heavy and diffuse blocks
 -> heavy blocks descend
 -> diffuse blocks analyzed by rare-fiber inverse theory
 -> repeated non-flattening creates one persistent rational-grid stabilizer
 -> Gaussian arithmetic rules out the persistent grid.
```

The genuinely new theorem remains the middle implication. Unlike earlier compressed audits, the exact requirements are now explicit: worst-cell sensitivity, conductor weighting, and constants independent of the rare-event mass fraction.