# Audit of the rational-lift four-torsion obstruction

Let `C/Q` be the smooth genus-five square-norm curve, let `J=Jac(C)`, and
let `X->C` be a normalized degree-1024 square-descent twist.  Assume that
`X` has a rational point `x` above `P in C(Q)`.  The proposed obstruction is
valid: if the branch polynomial has Galois group `S_12` and the affine
translation kernel on a branch fiber is zero, then the branch splitting
field contains `J[4]` and hence contains `i`.  Consequently zero kernel is
possible only when the discriminant has rational squareclass `-1`.

## 1. Why the rational lift identifies the pointed cover

Over an algebraic closure, constants in the normalized diagonal equations
are squares.  Ratios of their beta coordinates adjoin the square roots of

```text
(x-alpha_i)/(x-alpha_j).
```

The divisors are `2(W_i-W_j)`, and ten independent such classes span the
geometric two-torsion of `J`.  Thus the degree-1024 cover is the maximal
elementary-abelian unramified two-cover of `C`.  Geometrically it is the
pullback of `[2]:J->J` along an Abel--Jacobi embedding.

For the rational base point `P`, write

```text
X_P = {(R,Q) in C x J : 2Q=[R-P]}.
```

This is the same maximal geometric cover and it has the distinguished
rational point `(P,0)`.  A maximal elementary-abelian cover equipped with a
point above `P` is the pointed mod-two homology cover, hence is unique.  More
explicitly, the geometric subgroup in the etale fundamental group is the
characteristic kernel of the maximal exponent-two abelian quotient.  The
rational point supplies the decomposition-group splitting at `P`, which
selects the same arithmetic subgroup as `X_P`.  Therefore there is a unique
pointed isomorphism

```text
(X,x)  ~=  (X_P,(P,0))
```

over `Q`.  Uniqueness also gives descent: two conjugate isomorphisms differ
by a deck transformation fixing the marked point, hence by the identity.

An initially chosen identification of the geometric deck group with
`J[2]` may differ by an automorphism.  This does not affect the pointed
cover, its residue fields, or the identification of its fibers with the
halves occurring in `X_P`.

## 2. Zero branch kernel puts all relevant halves in the splitting field

Let `f` be the degree-twelve branch polynomial, assume
`Gal(f/Q)=S_12`, and let `K` be its splitting field.  The action on
`J[2]` is the faithful even-subset representation, so `K` is also the
field over which the linear action on `J[2]` becomes trivial.

Fix a Weierstrass point `W_j`.  The fiber of `X_P` above it is exactly

```text
[2]^-1([W_j-P]).
```

Over `Q(W_j)`, the linear image on this 1024-point affine space is `S_11`.
If its translation kernel is zero, every element acting trivially on the
linear module acts trivially on the whole affine fiber.  In particular,
`Gal(Qbar/K)` fixes every point of this fiber.  Thus all halves of
`[W_j-P]` belong to `J(K)`.

Normality of `K/Q` and transitivity of `S_12` on the roots give the same
conclusion for every `j`.  Equivalently, zero kernel for one root fiber is
the zero-kernel alternative globally; the only other alternative is the
full translation module.

## 3. Branch differences generate all four-torsion

Choose, for every `j`, a point `Q_j in J(K)` satisfying

```text
2Q_j=[W_j-P].
```

Then

```text
2(Q_j-Q_k)=[W_j-W_k].
```

The pair classes `[W_j-W_k]` span `J[2]`: in the even-subset description
they are the two-element subsets, and fixing `k` leaves ten independent
classes among the eleven displayed pairs.  Hence the `K`-rational
differences `Q_j-Q_k` provide lifts under `[2]` of a spanning set of
`J[2]`.

Since `J[2]` itself is rational over `K`, these lifts generate all of
`J[4]`.  Indeed, for any `R in J[4]`, express `2R` in the spanning pair
classes and subtract the same combination of the corresponding
`Q_j-Q_k`; the difference lies in `J[2](K)`.  Therefore

```text
J[4] subset J(K).
```

The canonical principal polarization is defined over `Q`.  Its Weil
pairing on `J[4]` is perfect and assumes a primitive fourth root of unity.
If every four-torsion point is `K`-rational, Galois equivariance of the
pairing forces

```text
mu_4 subset K,       hence i in K.
```

## 4. The discriminant obstruction and its scope

An `S_12` splitting field has exactly one quadratic subfield, because
`S_12` has the unique index-two subgroup `A_12`.  That subfield is

```text
Q(sqrt(disc(f))).
```

Thus `i in K` forces

```text
disc(f) = -1 in Q*/Q*2.
```

Taking the contrapositive, if the actual branch polynomial has full
`S_12` and discriminant squareclass different from `-1`, a normalized
twist possessing a rational point cannot have zero branch translation
kernel.  Simplicity of the two-torsion module forces the full kernel, so
every branch point has degree `12*1024=12288`.  By equality of the branch
and osculating-hyperplane fields, none of the degree-12 or degree-132
osculation cases occurs for that twist.

This is a restriction on normalized twists with a rational lift.  It does
not apply without the `S_12` hypothesis, does not determine the kernel when
the discriminant squareclass is `-1`, and does not prove a statement for
all six-point or all full-fair weights.  A single certified polynomial can
show that the obstruction is nonvacuous but cannot supply that blanket
conclusion.
