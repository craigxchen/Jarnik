# Pattern-rank cyclotomic rigidity

## Status

This note proves a restricted global theorem for the squarefree endpoint model.
It strictly extends the affine-cube/Walsh argument: no group structure on the
selected points is assumed.

The result uses Roth's theorem for approximation to one of finitely many fixed
algebraic slopes.  It does not prove the full uniform endpoint bound.  It shows
that every hypothetical unbounded family must have substantial linear
dependence among its active prime-orientation patterns.

---

## 1. Prime-pattern model

Fix `M` selected circle points from an odd squarefree conductor.  Group the
conductor primes according to their orientation pattern on these points.  For
each active nonconstant pattern `sigma_j in {+1,-1}^M`, let

```text
Gamma_j in Z[i]
```

be the product of the selected oriented Gaussian primes in that pattern class,
and put

```text
h_j = log Norm(Gamma_j),
H = sum_j h_j.
```

Patterns differing by one global sign are identified by conjugating `Gamma_j`.
The point angles have the form

```text
alpha_i = alpha_common + sum_j sigma_j(i) beta_j mod 2 pi,
beta_j = arg Gamma_j.
```

Assume all `alpha_i` have real lifts in one interval of width

```text
delta = C exp(-H/4).
```

Subtract the first row and define the integer pattern-difference matrix

```text
B_(i,j) = (sigma_j(i+1)-sigma_j(1))/2 in {-1,0,1},
1 <= i <= M-1.
```

Then for one integer vector `n` and an error vector `epsilon` with
`||epsilon||_infinity <= delta`, one has

```text
B beta = pi n + epsilon.                              (1.1)
```

---

## 2. Full column rank isolates every active block

Suppose `B` has full column rank `K`, where `K` is the number of active pattern
classes.  Choose `K` rows giving an invertible square integer matrix `A`, and
write

```text
D = |det A| > 0.
```

Multiplying the corresponding subsystem of (1.1) by the adjugate matrix gives,
for every `j`,

```text
D beta_j = pi m_j + O_A(delta)                       (2.1)
```

for some integer `m_j`.  Equivalently,

```text
dist(beta_j, (pi/D) Z) <= C_A delta.                 (2.2)
```

Thus every Gaussian pattern block lies exponentially close to one of the fixed
rays of angle `pi m/D`.

This uses only linear independence of the active patterns.  It does not require
the selected point set to be an affine cube or the pattern matrix to be
Hadamard.

---

## 3. Arithmetic separation from fixed cyclotomic rays

Write a primitive associate of `Gamma_j` as

```text
a_j + i b_j.
```

For a fixed target ray `phi = pi m/D`, there are two cases.

### Irrational slope

If `tan phi` is irrational, it is a fixed algebraic number.  Roth's theorem
implies that for every fixed `eta>0`,

```text
|b_j/a_j - tan phi|
  >= c_(phi,eta) |a_j|^(-2-eta)
```

unless equality holds.  Equality is impossible because the left side is
rational and `tan phi` is irrational.

Angular distance and slope distance are comparable in a fixed neighborhood of
the ray.  Since

```text
|a_j| <= Norm(Gamma_j)^(1/2),
```

we obtain

```text
dist(beta_j, phi)
  >= c'_(phi,eta) exp(-(1+eta/2) h_j).               (3.1)
```

Comparing (2.2), (3.1), and `delta=C exp(-H/4)` gives

```text
h_j >= H/(4+2 eta) - O_(C,A,eta)(1).                 (3.2)
```

### Rational Gaussian rays

For a rational multiple of `pi`, a rational value of the tangent can only give
the Gaussian rational rays with slopes `0`, `infinity`, or `+/-1`.  Elementary
integer separation gives a stronger lower bound of order

```text
Norm(Gamma_j)^(-1/2)
```

unless the ray is exact.

An exact odd-squarefree Gaussian block on one of these rays is impossible after
primitive normalization: a real or imaginary Gaussian integer has rational
prime support in conjugate pairs, and a slope `+/-1` introduces the ramified
factor `1+i`.  Neither is an active product of one selected factor above each
odd split prime.

Thus (3.2) also holds, with room to spare, for the rational target rays.

---

## 4. Consequence

Fix a small positive `eta`.  For sufficiently large `H`, every active pattern
block satisfies, for example,

```text
h_j > H/5.
```

Since the block heights are disjoint and sum to `H`, there can be at most four
active nonconstant pattern classes.

Four binary patterns distinguish at most `2^4=16` point rows.  Therefore:

### Pattern-rank rigidity theorem

For every fixed endpoint constant `C`, a sufficiently large squarefree endpoint
cluster whose active nonconstant orientation-pattern columns are linearly
independent has at most sixteen points.

Equivalently, every hypothetical unbounded endpoint family must eventually
have

```text
number of active pattern classes > rank_Q(B).        (4.1)
```

The hard case therefore contains a nontrivial and growing pattern kernel.  It is
not merely diffuse in prime support; many distinct conductor blocks act through
linearly dependent sign patterns on every fixed large selected subcluster.

---

## 5. Why the proof stops at dependent patterns

When `K > rank B`, equation (1.1) determines only the image of the block-angle
vector modulo `ker B`.  Large angle contributions can hide in this kernel, and
no individual `beta_j` is forced near a fixed cyclotomic ray.

A rational row-space combination supported on a small collection of columns
would still produce a near-root-of-unity Gaussian product.  The arithmetic
separation above would then force that collection to carry a positive fraction
of `H`.  The missing combinatorial statement is therefore a bounded-coefficient
cocircuit theorem for the active hypercube-pattern matrix.

Arbitrary integer matrices have cocircuits with very large coefficients, and
that coefficient growth destroys the conductor saving.  The next question is
whether matrices whose columns are actual binary orientation patterns admit a
stronger alternative:

```text
small bounded-coefficient cocircuit,
  or
large structured kernel generated by low-complexity binary relations.
```

The first branch would extend the present arithmetic argument.  The second
branch would expose a concrete global sign symmetry, potentially allowing
factorization or descent.
