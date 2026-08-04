# Affine-cube route: Walsh rigidity and the sparse-set obstruction

## Status

This note pushes the affine-cube program to a genuine arithmetic theorem in the
structured regime and identifies the precise remaining obstruction.

The uniform endpoint bound is **not** proved here.

What is proved, modulo the standard imported form of Roth's theorem, is:

> A sufficiently large endpoint cluster cannot contain a five-dimensional
> Boolean affine subcube on which every Gaussian prime-power layer restricts to
> an affine character.

Consequently, after the usual binary threshold-layer expansion, any hypothetical
unbounded counterexample must be asymptotically additive-sparse. In particular,
a bad family cannot have normalized additive energy bounded below by a positive
constant, once standard Balog--Szemeredi--Gowers/Freiman machinery is imported.

The remaining case is a Sidon-like family with almost all pair differences
distinct. Mere affine-cube extraction cannot handle that case; there are
arbitrarily large abstract cube-free sets.

---

## 1. Binary threshold-layer expansion

Write the split part of the common norm as

```text
N = product_j p_j^(e_j),
```

and choose `pi_j` above `p_j`. A point with exponent `a_j` may be encoded by the
`e_j` binary threshold signs

```text
sigma_(j,t)(a) = +1  if t <= a_j,
                 -1  if a_j < t,

1 <= t <= e_j.
```

The corresponding Gaussian factor in layer `(j,t)` is `pi_j` in the `+1`
state and `conjugate(pi_j)` in the `-1` state. Thus every point is a binary
orientation word in a cube

```text
{+1,-1}^L,       L = sum_j e_j,
```

subject to the monotonicity constraints inside each prime block.

Each layer has conductor weight

```text
w_(j,t) = log p_j
```

and angular coefficient

```text
theta_(j,t) = arg(pi_j).
```

The total conductor height is

```text
H = log N = sum_(j,t) w_(j,t).
```

The endpoint angular width is

```text
delta = C exp(-H/4).
```

---

## 2. Character-compatible Boolean cubes

Let

```text
V = F_2^d,
M = 2^d.
```

Suppose the cluster contains distinct points `z_x`, indexed by `x in V`, with
the following property.

For every binary Gaussian layer, its orientation pattern on these `M` points is
an affine character. After absorbing a constant sign by conjugating the chosen
Gaussian factor, the nonconstant pattern is

```text
chi_xi(x) = (-1)^(xi dot x)
```

for some nonzero `xi in V^*`.

Group every layer having the same character. Let `Gamma_xi` be the product of
its chosen Gaussian factors, let

```text
Q_xi = Norm(Gamma_xi),
h_xi = log Q_xi,
beta_xi = arg(Gamma_xi).
```

Constant character classes contribute the same oriented factor to all points
and can be removed from the discussion.

The points have the form

```text
z_x = common * product_xi
        Gamma_xi^((1+chi_xi(x))/2)
        conjugate(Gamma_xi)^((1-chi_xi(x))/2).
```

Hence

```text
arg z_x = constant + sum_xi chi_xi(x) beta_xi  mod 2 pi.
```

---

## 3. Exact Walsh isolation

Fix a nonzero character `eta`. Orthogonality gives

```text
sum_x chi_eta(x) = 0,

sum_x chi_eta(x) chi_xi(x)
  = M  if xi = eta,
    0  otherwise.
```

Therefore the multiplicative Walsh product is exactly

```text
product_x z_x^(chi_eta(x))
  = (Gamma_eta / conjugate(Gamma_eta))^(M/2).       (3.1)
```

Negative exponents are interpreted in `Q(i)^*`; all common factors cancel.

Choose real lifts `alpha_x` of the point arguments inside one interval of length
`delta`. Since exactly half the values of `chi_eta` are `+1` and half are `-1`,

```text
|sum_x chi_eta(x) alpha_x| <= (M/2) delta.
```

Combining this with (3.1), there is an integer `k_eta` such that

```text
|M beta_eta - 2 pi k_eta| <= (M/2) delta,
```

and consequently

```text
|beta_eta - 2 pi k_eta/M| <= delta/2.              (3.2)
```

Thus every active Gaussian block is forced into a `delta/2`-neighborhood of one
of the fixed `M` rays of angle `2 pi k/M`.

This is the higher-order gain supplied by the complete affine cube. Pairwise
separation sees only that ratios are near `1`; Walsh isolation recovers each
character block separately.

---

## 4. Arithmetic separation from the fixed rays

Fix `M=32` from now on. The target rays in (3.2) are fixed algebraic directions.

Write a primitive associate of `Gamma_eta` as

```text
a + i b,
```

without changing its angle. If the target ray has finite irrational slope

```text
tau = tan(2 pi k/32),
```

Roth's theorem gives, for every fixed `epsilon>0`,

```text
|b/a - tau| >= c_(k,epsilon) |a|^(-2-epsilon)
```

unless `b/a=tau`, which is impossible for irrational `tau`.

Near a fixed nonvertical ray, slope distance and angular distance are comparable.
Since

```text
|a| <= sqrt(Q_eta),
```

we obtain

```text
|beta_eta - 2 pi k/32|
  >= c'_(k,epsilon) Q_eta^(-1-epsilon/2).          (4.1)
```

For a vertical target, apply the same argument to `a/b`. For the rational rays
of slopes `0`, `infinity`, and `+/-1`, the elementary determinant estimate gives
an even stronger lower bound of order `Q_eta^(-1/2)`, unless the ray is exact.

If a nonconstant character block lies exactly on one of those rational Gaussian
rays, changing its sign rotates points by a nontrivial Gaussian unit. Such a
block cannot take both signs inside an interval of width tending to zero. The
only harmless exact case gives the same oriented factor on both signs and hence
cannot distinguish points.

Take, for example, `epsilon=1/4`. Equations (3.2) and (4.1), with

```text
delta = C exp(-H/4),
```

imply, for all sufficiently large `H`,

```text
h_eta >= (2/9) H - O_C(1).                         (4.2)
```

The numerical constant is not optimized. What matters is that it is strictly
larger than `H/5` asymptotically.

Since the total of the character-block heights is at most `H`, at most four
nonconstant character classes can be active.

But four binary characters distinguish at most

```text
2^4 = 16
```

points. They cannot distinguish all `32` points of `F_2^5`.

We have therefore proved:

### Walsh rigidity theorem

For every fixed endpoint constant `C`, a character-compatible Boolean affine
cube of dimension five cannot occur in a `C N^(1/4)` arc once `N` is sufficiently
large.

Equivalently, every such cube contains at most `16` distinct points for large
conductor height.

The threshold is ineffective because the proof imports Roth's theorem, but it
is uniform in the identities and number of Gaussian primes inside each
character block.

---

## 5. Why dimension five is enough globally

If a character-compatible cube has dimension `d>=5`, restrict it to any
five-dimensional affine subcube. Every original affine character restricts to
an affine character (possibly constant) on that subcube. Grouping layers by the
restricted character gives exactly the preceding `32`-point model.

Therefore **no character-compatible Boolean affine cube of dimension at least
five** can occur in a sufficiently large endpoint cluster.

---

## 6. Consequence in the high-energy regime

Let `B` denote the binary threshold-layer code of a cluster, viewed as a subset
of an ambient `F_2`-cube. If `B` contains a five-dimensional affine subspace,
then every coordinate restriction is an affine `F_2`-linear function, so its
sign is an affine character. The Walsh rigidity theorem applies.

Thus a hypothetical bad family must avoid five-dimensional affine subspaces in
its binary orientation code.

There is a standard imported additive-combinatorial consequence:

1. If the normalized additive energy satisfies

   ```text
   E(B) >= |B|^3 / K
   ```

   for one fixed `K`, Balog--Szemeredi--Gowers gives a large subset with bounded
   doubling.

2. Freiman-type theory in `F_2^n` places that subset densely inside a subspace of
   size `O_K(|B|)`.

3. A fixed positive-density subset of an `F_2`-space of sufficiently large
   dimension contains a five-dimensional affine subspace.

Combining these imported results with Walsh rigidity shows:

> In every hypothetical unbounded endpoint family,
>
> ```text
> E(B) / |B|^3 -> 0.
> ```

So the affine-cube program completely rules out the high-energy/small-doubling
regime.

This is a genuine structural reduction, not a restatement of the original
problem.

---

## 7. The concrete sparse-set obstruction

Mere unboundedness does not force an affine cube of any fixed dimension.

In an ambient `F_2^n`, a random set of size `2^(c n)` with sufficiently small
fixed `c>0` has, after deleting one point from every offending configuration,
an unbounded subset containing no nontrivial two-dimensional affine cube. A
simple first-moment count already gives such constructions for small `c`.

Thus there are arbitrarily large abstract Sidon-like subsets of binary cubes.
They have

```text
E(B) = (2+o(1)) |B|^2,
```

rather than energy comparable with `|B|^3`.

The endpoint hypothesis does not presently give a lower bound for additive
energy in the exponent lattice. Therefore the standard affine-cube extraction
machinery stops here.

This is the exact failure of the naive statement

```text
unbounded angular concentration -> many affine cubes.
```

It is false without additional Gaussian arithmetic.

---

## 8. The remaining arithmetic target

After the preceding reduction, a hypothetical counterexample must be
simultaneously:

- diffuse in conductor support;
- asymptotically balanced at almost every prime layer;
- free of large Boolean affine subcubes;
- additive-Sidon-like, with almost all pair differences distinct;
- nevertheless mapped by the Gaussian-prime angle form into an interval of
  width `C exp(-H/4)`.

The remaining theorem is therefore no longer a general inverse theorem. It is a
**sparse arithmetic inverse theorem**:

> Bound a Sidon-like family of Gaussian divisor exponent vectors whose pairwise
> ratios are all within `C exp(-H/4)` of `1`.

The abstract equal-angle and arbitrary-real-angle models show that combinatorics
alone cannot prove this. The next input must use arithmetic separation among
many distinct Gaussian-prime angle combinations simultaneously.

Possible next tools are:

- a simultaneous Roth/Subspace-Theorem statement tailored to a Sidon family of
  exponent differences;
- a higher common-divisor theorem for the many distinct Gaussian numerators
  produced by those differences;
- a sector-packing theorem for a sparse family of Gaussian divisors, stronger
  than pairwise slope separation because all divisors come from one conductor.

---

## 9. Verdict on the affine-cube approach

The approach does not finish the uniform bound, but it does more than merely hit
an obstruction:

```text
high additive energy
  -> affine subcube
  -> exact Walsh isolation
  -> approximation to fixed algebraic rays
  -> Roth contradiction / at most four active characters
  -> at most sixteen points.
```

The only surviving regime is low additive energy. An arbitrary large set can be
cube-free, so no further purely affine-cube argument can force progress there.
Gaussian arithmetic has to enter directly in the sparse regime.
