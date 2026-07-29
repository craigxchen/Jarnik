# Codex handoff: quantitative uniformity audit for the common-conductor route

## Repository and branch

Repository: `craigxchen/Jarnik`

Branch: `agent/uniform-gcd-formalization`

Draft PR: #1

## Goal

Determine whether the existing Ru--Vojta / Evertse--Ferretti machinery actually yields a uniform bound for endpoint-scale lattice-point arcs on circles when all divisor ratios come from one common Gaussian conductor.

Do not assume the desired uniformity. The task is falsification-oriented: either produce a complete theorem chain with every dependency explicit, or identify the first parameter that still depends on the number of places or on point-dependent local data.

## Mathematical setup

Let `K = Q(i)`. For one Gaussian conductor `G`, the relevant norm-one ratios are

```text
u_A = A / conjugate(A),   A | G.
```

For a cluster on one circle, every ordered pair `(u_A,u_B)` satisfies:

1. a common conductor-height bound

```text
h([1:u_A:u_B]) <= log |G|;
```

2. endpoint proximity

```text
log GCD^+(u_A-1,u_B-1) >= (1/2) log |G| - O_C(1).
```

At the fixed Ru--Vojta level `(N,d)=(1,6)`, the branch proves

```text
h^0(6H-E) = 27,
beta = 50/27,
6 - 3 beta = 4/9 < 1/2.
```

Thus a theorem giving, outside a uniformly bounded exceptional family,

```text
log GCD^+(x-1,y-1) <= (4/9 + epsilon) log |G| + O(1)
```

with `epsilon < 1/18` would prove the desired uniform arc bound via the formalized grid argument.

## Current consensus

The branch explored and audited several approaches:

- unrestricted uniform-in-`S` exceptional degree: false in that generality;
- conductor-relative exceptional curve as a standalone theorem: essentially equivalent to the endpoint statement;
- canonical Harder--Narasimhan flats: finite for fixed form data, but do not automatically control all fixed-scale auxiliary exceptional spaces;
- arrangement-closure identity: false;
- product sections from adapted bases: attractive but the set of adapted bases may be infinite;
- pair-dependent twisted weights: probably unnecessary if the intrinsic Ru--Vojta beta theorem can be used directly.

The current favored route is:

```text
fixed blowup / fixed divisors
+ intrinsic finite-level beta inequality
+ quantitative higher-degree Subspace Theorem
+ uniform exceptional component/degree bounds
+ grid lemma.
```

## Exact task

Audit the strongest usable quantitative theorem of Ru--Vojta, Evertse--Ferretti, or a direct higher-degree Subspace Theorem and answer all of the following.

### 1. Fixed geometric data

Can the application be formulated with a projective variety and divisor/form family that are genuinely fixed as the finite place set `S` varies?

In particular, distinguish:

- repeated use of the same divisors/forms at different places;
- the number of distinct algebraic forms;
- the number of places carrying nonzero local data.

### 2. Exceptional complexity

Extract exact bounds or dependency statements for:

```text
T      = number of exceptional subvarieties,
Delta  = maximal degree of each exceptional subvariety,
H0     = height or parameter threshold.
```

Determine whether each depends on:

- ambient dimension;
- degree and height of the fixed variety;
- number and degree of distinct forms/divisors;
- approximation margin;
- field degree;
- `|S|`;
- the individual local weights;
- heights of place-dependent coefficients.

Do not infer independence from an omitted parameter. Trace definitions through every referenced theorem or proposition.

### 3. Repeated-form issue

If the theorem lists forms `L_i^(v)` separately at every place, determine whether repeating one fixed form at many places enlarges the theorem's complexity parameter. Check whether the relevant quantity is:

```text
number of indexed form occurrences
```

or

```text
number of distinct forms up to proportionality.
```

This is likely the decisive point.

### 4. Threshold uniformity

Even if component number and degree are uniform, verify that the lower threshold for height or the parameter `Q` is uniform as `S` varies.

Watch for hidden dependence through:

- a product of local coefficient heights;
- denominators introduced by normalization;
- the support size of the weight vector;
- a constant defined as a sum over places;
- conversion between Weil functions and chosen local equations.

### 5. Ru--Vojta interface

Write the exact finite-level inequality on `Bl_[1:1:1] P^2` at `(N,d)=(1,6)` and show how it matches the hypotheses of the quantitative theorem.

The proof must make clear whether:

- the adapted-basis maximum is already internal to the Ru--Vojta theorem;
- any auxiliary forms depend on the point `(u_A,u_B)`;
- any local coefficient heights depend on the conductor primes;
- one theorem application covers the entire divisor-pair grid for a fixed `G`.

### 6. Pushdown and grid endgame

If the theorem produces subvarieties in a blowup, embedding, or auxiliary projective space, bound the total degree of their images or pullbacks in `P^2` uniformly.

Then connect the result to:

- `GaussianChain/FiniteExceptionalFamily.lean`, or
- `GaussianChain/UniformExceptionalReduction.lean`.

State the final numerical coefficient and verify that all error terms fit inside the margin

```text
1/2 - 4/9 = 1/18.
```

## Required output

Create `docs/codex_quantitative_audit_result.md` with one of two outcomes.

### Outcome A: proof chain succeeds

Provide:

1. exact theorem statements and citations;
2. a dependency table for every constant;
3. the full specialization to `Q(i)` and `(N,d)=(1,6)`;
4. a complete uniform exceptional-degree theorem;
5. the final arc-cardinality deduction;
6. a list of remaining formalization tasks only, not mathematical gaps.

### Outcome B: proof chain fails

Provide:

1. the earliest failed implication;
2. the exact constant or object that still depends on `|S|`, conductor prime identities, or point-dependent data;
3. a quotation or precise theorem-definition reference supporting that dependence;
4. why the common-conductor condition does or does not remove it;
5. the narrowest genuinely new theorem that would close the gap.

## Repository files to read first

- `docs/three_track_consensus_audit.md`
- `docs/conductor_relative_route.md`
- `docs/conductor_relative_equivalence.md`
- `docs/one_scale_quantifier_audit.md`
- `docs/product_section_audit.md`
- `docs/product_section_pairing_audit.md`
- `docs/arrangement_closure_audit.md`
- `GaussianChain/FiniteBeta.lean`
- `GaussianChain/UniformExceptionalReduction.lean`
- `GaussianChain/FiniteExceptionalFamily.lean`

Treat later audit files as corrections to earlier optimistic claims.

## Standards

- Do not claim a proof from theorem names or abstracts.
- Track quantifiers and parameter dependencies line by line.
- Prefer primary papers and exact theorem numbering.
- Clearly separate proved statements, imported theorems, and new conjectural lemmas.
- Do not modify the PR description to claim completion unless the entire chain is verified.
