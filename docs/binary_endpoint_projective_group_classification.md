# Retained binary cuts force a free projective two-group

Let `Z={z_1,...,z_m}` be a Gaussian-primitive tuple of distinct lattice
points on the circle of squared radius `N`, in an arc of length at most
`2N^(1/4)`. Assume its allocations at every source prime are binary;
equivalently, every physical point has ordinary coordinate content one.
If

```text
m>=32,       N>2^36,
```

then its rational projective automorphism group is one of

```text
{1},       C_2,       C_2 x C_2.                    (1)
```

The group acts freely on the sampled points. Thus its order divides `m`,
and every odd-cardinality tuple under these hypotheses has trivial
projective automorphism group. This is a structural obstruction to the
self-permutation route to a nonconformal endpoint map, not a construction
of arbitrarily large endpoint tuples.

The earlier [fixed-anchor content obstruction](binary_endpoint_projective_selfmap_obstruction.md)
used a fixed point and a restricted cardinality class. The argument here
instead uses retained-cut compatibility to constrain every permutation
and every cardinality; it does not require a fixed point initially.

## 1. Retained cuts can only be preserved or complemented

Fix a nonconformal projective automorphism, with induced permutation
`sigma`, and write `Q=Nclip`, `E=N/Q`. At an odd split prime `p`, let
the binary physical allocations be `0,e`. They partition the labels into
a nonempty proper cut `S` and its complement. The image tuple is the
same physical tuple, permuted, so it too is binary with these extrema.

If `v_p(Q)>0`, the clipped values at the two source extrema are distinct,
and the same is true for the inverse clipped target extrema. The exact
forward/inverse clipping identity is an affine isometry between these
two clipped intervals. Therefore a source fiber cannot have images at
both target extrema. Bijectivity now forces

```text
sigma(S)=S or sigma(S)=S^c.                         (2)
```

This conclusion permits partial retained depth and all endpoint
cancellations. It only requires positive clipped width. In particular
`sigma^2` preserves every retained cut.

## 2. A permutation invisible on retained cuts cannot move a point

For distinct labels `i,j`, let the odd primitive edge conductor be

```text
n_ij=product p^e over source primes where the binary allocations differ.
```

If a permutation `tau` preserves every retained cut, every prime
contributing to `n_(i,tau(i))` has clipped width zero. Consequently

```text
n_(i,tau(i)) | E whenever tau(i)!=i.                (3)
```

The full primitive relative half-angle norm is `f_ij=epsilon_ij n_ij`,
where `epsilon_ij` is one or two. The endpoint chord inequality and
the nonzero integer primitive pair residue give

```text
f_ij>sqrt(N),       n_ij>sqrt(N)/2.                 (4)
```

The factor two is necessary for arbitrary Gaussian unit classes.
Strictness follows from the strict chord-versus-positive-arc inequality.

For a nonconformal self-map with both endpoint constants at most two,
the erased-conductor theorem gives

```text
E <= 2^(m q_m) N^q_m,
q_m<=8/m,
E <= 2^8 N^(1/4) < sqrt(N)/2                       (5)
```

under the stated hypotheses. Equations (3)--(5) imply `tau(i)=i` for
every label. Apply this first to `tau=sigma^2`, using (2), to obtain
`sigma^2=1`.

If `sigma` fixes any sampled label, the complement option in (2) is
impossible for every retained cut: a fixed vertex cannot cross the cut.
Thus `sigma` itself preserves every retained cut, and the same argument
gives `sigma=1`. Every nonidentity nonconformal automorphism is therefore
a fixed-point-free involution on the sample.

For odd `m`, one can also argue directly: complementation in (2)
requires equally sized cut fibers, so it never occurs. The same edge
argument immediately makes every nonconformal permutation trivial.

## 3. Conformal symmetries are excluded separately

A common rotation preserving the primitive physical tuple has multiplier
`rho` with `rho z_i` Gaussian integral for every `i`. Gaussian Bezout
forces `rho` to be Gaussian integral; `|rho|=1` makes it a Gaussian unit.
A nonidentity unit rotation cannot preserve a set in a proper arc of
width less than `pi`.

The same argument applied to the conjugate tuple shows that a preserving
reflection is `z -> rho conjugate(z)` for a Gaussian unit `rho`. Its axis
is a coordinate axis or a diagonal. The primitive squared radius is odd.
Distinct radial projection levels have spacing at least one on a
coordinate axis and at least `sqrt(2)` on a diagonal (the diagonal
numerators are odd).

The shortest containing arc is reflection-invariant and centered on the
relevant axis ray. Its radial sagitta is at most `C^2/8<=1/2`. Hence it
contains at most one projection level and at most two circle points,
contrary to `m>=32`. This is the existing
[primitive reflection-union radial bound](reflection_union_radial_count_obstruction.md),
applied to a tuple already invariant under the reflection.

Thus no nonidentity conformal automorphism remains.

## 4. Finite group classification and scope

A projective map fixing three distinct directions is the identity, so
the automorphism group embeds faithfully in the finite permutation group
of the sample. Its orientation-preserving subgroup is cyclic, because
it acts by cyclic shifts of the circle order. Every nonidentity group
element is now an involution. The cyclic orientation-preserving subgroup
therefore has order at most two, and the full group has order at most
four. Since every nonidentity element has order two, the possibilities
are exactly (1). Their actions are free by Section 2, proving the
divisibility of `m` by the group order.

No classification of distinct same-radius configurations in one projective
orbit is asserted. In particular, the result does not rule out a new
image configuration with suitable clipping and radius bounds, nor does
it supply such an image. The four-seed consequence without a binary-target
hypothesis is recorded separately in
[the fixed secant conductor note](fixed_secant_conductor_bridge.md).

The [checker](check_binary_endpoint_projective_group.py) verifies finite
permutation/cut compatibility, the erased-edge divisor implication,
the primitive-edge parity factor, and the uniform exponent certificates.
