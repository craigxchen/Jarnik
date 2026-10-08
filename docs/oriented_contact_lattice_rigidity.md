# Integral existence and rigidity of the oriented contact lattice

The contact congruences in
[the five-point norm-lift note](five_point_cm_norm_lift.md) always admit
a determinant-one integral frame. This remains true when the congruence
slopes come from the same five primitive vectors, so their barycentric
relations are retained. What the congruences do not control is the
height of that frame or of the additional row factors.

There is nevertheless an exact height restriction: two distinct short
frames cannot realize the same oriented data. More generally, two short
frames realizing different orientations must disagree on a substantial
amount of full prime-power weight. Neither statement proves that the
first short frame is impossible.

## 1. Data and integral existence

Put `O=Z[i]`. Let `G_e` be odd Gaussian integers, pairwise coprime and
coprime to every conjugate, including their own. Let
`lambda_e=(r_e,s_e)` be primitive integer vectors, and set

```text
Delta=product_e G_e,       D=Norm(Delta),
L={x=(z,w) in O^2: r_e z+s_e w=0 mod G_e for every e},
q(z,w)=Im(bar(z)w).
```

Then there exists `x in L` with `q(x)=1`.

To prove this, write `z=a+ic`, `w=b+id`, so `q(x)=det U` for
`U=[[a,b],[c,d]]`. At each rational prime power `p^e` contributing
to `D`, exactly one Gaussian orientation `pi^e` occurs in `Delta`.
The ring `O/(pi^e)` is `Z/(p^e)`, with `i` represented by a root
`iota` of `-1`. The congruence is

```text
(1,iota) U lambda_e=0 mod p^e.
```

Over `Z/(p^e)`, complete the primitive vector `lambda_e` to an
`SL_2` matrix. Complete the primitive target `(-iota,1)` to another
such matrix. Their quotient is an `SL_2` matrix satisfying the
required congruence. The Chinese remainder theorem combines these
local matrices into an element of `SL_2(Z/DZ)`.

The reduction map `SL_2(Z)->SL_2(Z/DZ)` is surjective. An explicit
proof is to make the first entry of a primitive column a unit by an
elementary shear, choosing its parameter independently at each prime
dividing `D`. Elementary row operations then reduce to a diagonal
matrix `diag(u,u^(-1))`; that matrix is a product of elementary
matrices as well. Lifting the finitely many elementary parameters
to integers gives the required lift `U in SL_2(Z)`.

If all slopes are selected from one collection of primitive vectors
`v_i`, and vectors in contact cluster `T_e` have determinants divisible
by `Norm(G_e)`, each is a unit multiple of `lambda_e` modulo `G_e`.
The same `U` therefore satisfies all cluster congruences. Every
projective cross ratio and barycentric identity of the original
vectors is preserved, since the transformation is common to them.

One may also impose any compatible condition `U mod 2` in the CRT.
In particular, if the input vectors admit an `SL_2(F_2)` image with
all row norms odd, that parity can be retained. This parity hypothesis
is necessary: all three nonzero directions of `F_2^2` cannot all be
sent to the two directions of odd norm. A five-vector frame already
arising from odd Gaussian rows has the required compatible parity.

This is an existence theorem for the congruences and determinant one.
It supplies no `exp(4w+o(w))` height bound. Extra prime factors of
the resulting rows can be large, so it does not supply the small
corrections required by a full profile.

## 2. Exact clean contact assignments

There is a stronger local assertion under an explicit cleanliness
hypothesis. Suppose that at a split `p^e`, a contact cluster has at
most three vectors, all primitive and proportional to `lambda` modulo
`p^e`, while every outside vector has a different projective reduction
from `lambda` modulo `p`. There are at most five vectors in total.
Then an `SL_2` frame can be chosen so that:

```text
v_pi(Q_i)=e,  v_bar(pi)(Q_i)=0       for i in the cluster,
v_pi(Q_i)=v_bar(pi)(Q_i)=0          for every outside i.       (1)
```

Indeed `p>=5`. Choose a projective direction distinct from the cluster
and all the outside directions as the preimage of the conjugate
isotropic line. There is room: `P^1(F_p)` has at least six points,
and at most four are excluded. An `SL_2(F_p)` matrix can send the
cluster and this second direction to the two respective isotropic
lines. Thus outsiders avoid both orientations modulo `p`.

Lift this matrix modulo `p^(e+1)`, making the cluster line isotropic
modulo `p^e`. In a primitive integral basis beginning with `lambda`,
each cluster vector has coordinates `(a_i,p^e b_i)` with `a_i` a
unit modulo `p`. The transverse coordinate of the image of `lambda`
at depth `e` is a free parameter `c in F_p`. The depth-`e` leading
coefficient of the image of the `i`th cluster vector is
`a_i c+b_i d`, where `d` is a fixed unit. At most three choices of
`c` make one of these coefficients zero. Choose another value.
This gives exact depth `e` for every cluster member; the outside
unit conditions persist from the reduction modulo `p`.

CRT modulo the resulting `p^(e+1)` powers, followed by the same
`SL_2(Z)` lifting, realizes all clean contacts simultaneously. Their
vectors still come from the same five-vector configuration.

The outside-distinctness hypothesis must not be silently imposed on
an arbitrary full profile. An outside vector can agree with a cluster
modulo `p` to a small residual depth while the core `p^e` is large.
Discarding that entire prime would then lose unbounded core weight.
Without cleanliness, Section 1 still proves the stated congruence
existence, but (1) is not being asserted.

## 3. Duality and an integral normal form

Use the Hermitian form

```text
h(x,y)=bar(x)^t [[0,-i],[i,0]] y,       h(x,x)=2q(x).
```

The real congruence slopes and conjugate-coprimality give

```text
L+bar(L)=O^2,
L^dual=(bar(Delta))^(-1) bar(L).                       (2)
```

For the second identity, `h(x,y)=-i det(bar(x),y)` and the
alternating dual of a rank-two module with determinant ideal
`(bar(Delta))` is that module divided by `bar(Delta)`. The first
identity follows locally: at each Gaussian prime, at least one
of `L` and `bar(L)` is the full ambient module.

Consequently the discriminant module has the exact form

```text
L^dual/L = O/(Delta) direct_sum O/(bar(Delta)).         (3)
```

Its subgroup `O^2/L` is automatically maximal isotropic for the
finite discriminant pairing: the ambient form is integral and
unimodular. Thus a local discriminant obstruction is not hidden
in the oriented congruences.

There is an integral normal form, not just a statement away from two.
Choose a generator of `(Delta)` with `Delta=1 mod 2`; this is possible
by a Gaussian unit because `Delta` is odd. For any `x in L` with
`q(x)=1`, put

```text
y=(x+Delta bar(x))/2.                                 (4)
```

The vector is integral, since `bar(x)=x mod 2`. Also `2y in L`,
and two is invertible modulo every contact, so `y in L`. Directly,

```text
h(x,bar(x))=0,       q(bar(x))=-1,
det(x,bar(x))=-2i,
det(x,y)=-i Delta.
```

Thus `(x,y)` is an `O`-basis of `L` and its Hermitian matrix is

```text
[[2,1],[1,(1-D)/2]].                                  (5)
```

Equivalently, in that basis,

```text
q(a x+b y)=Norm(a)+Re(bar(a)b)+(1-D)Norm(b)/4,
q(a x+b y)=1  iff  Norm(2a+b)-D Norm(b)=4.              (6)
```

The congruence `2a+b=b mod 2` is retained when the latter equation
is parametrized by `2a+b` and `b`. Section 1 guarantees existence
of the initial vector `x`; therefore all these contact lattices
have the integral Hermitian class (5). Orientations and barycentric
data remain encoded in the Euclidean position of the basis, not
in a different abstract class of this Hermitian form. A basis supplied
by the CRT can be very large.

## 4. Uniqueness of a short determinant-one lift

For any `x,x' in L`,

```text
det(x,x') is divisible by Delta.                      (7)
```

This is the determinant-ideal identity, or follows contact by
contact because both vectors lie in one rank-one kernel modulo
each `G_e`. If they are linearly independent over `Q(i)`, the
complex determinant inequality gives

```text
||x||_2 ||x'||_2 >= |det(x,x')| >= |Delta|.             (8)
```

If instead they are collinear and both represent `q=1`, then
`x'=u x` for a Gaussian unit `u`. Indeed `q=1` makes each vector
primitive over `O`; integrality forces the proportionality factor
to lie in `O`, and its norm must be one.

For the five-row height budget,

```text
log|Delta|=10w+o(w),       log||x||_2=4w+o(w).
```

It follows that there is at most one such short lift up to a common
Gaussian unit. Every nonunit alternative has

```text
log||x'||_2 >= 6w-o(w).                               (9)
```

This statement retains the complete contact modulus, including its
prime powers. It differs from reconstruction of an individual small
correction ratio: it is rigidity of the common determinant-one frame.

## 5. A lower bound on orientation disagreement

Consider two orientation allocations `G_e,G'_e` for the same contact
labels and the same coefficient vectors. Define their common oriented
modulus with exact exponent overlap by

```text
C=product_e gcd_G(G_e,G'_e).                          (10)
```

If `x` and `x'` solve their respective contact systems, their
determinant is divisible by `C`. Unless they differ by a common
Gaussian unit, (8) therefore holds with `C` in place of `Delta`.
In particular, if both have logarithmic height `4w+o(w)`,

```text
log Norm(C) <= 16w+o(w).                              (11)
```

When each allocation has total contact norm weight `20w+o(w)`,
each must have at least `4w-o(w)` of oriented prime-power weight
outside their common modulus. If the rational contact factors are
identical, this is a lower bound on the weight whose orientation
changes. Formula (10) uses minimum exponents, not radical support.

For a genuine orientation flip at a prime, unit-equivalent frames
cannot solve both allocations. The determinant-one identity makes
`lambda_e bar(x)` a unit at a prime where `lambda_e x` vanishes;
vanishing in both orientations is impossible. Thus two different
orientation allocations on fixed rational contacts, if both possess
short lifts, obey the disagreement bound (11).

The results give unrestricted existence, a canonical Hermitian form,
and quantitative rigidity of any exceptional short realization.
They do not decide whether one orientation allocation admits the
required short frame with small additional row factors.

The [exact checker](check_oriented_contact_lattice_rigidity.py) verifies
60 CRT lifts with arbitrary primitive slopes, prime-power moduli, and
parity. It also constructs one actual common five-vector system
realizing all ten clean median contact assignments, including a
depth-two contact at five, and checks its cross ratios, barycentric
identities, and canonical integral Hermitian basis. That example
illustrates unrestricted integral existence; its heights are not
claimed to satisfy the five-row full-profile budget.

## 6. The same rigidity for any number of rows

Fix `m>=3`, retain three median anchors, and put `A_m=2^(m-3)`.
After pairing complementary original cuts, a surviving nonconstant
pattern is a nonempty set `T` containing at most one median anchor.
There are exactly

```text
4 * 2^(m-3)-1 = 4 A_m-1
```

such patterns. Every paired block has logarithmic norm `2w+o(w)`.
The `m` singleton patterns are private row factors. Each other pattern
has at least two rows, so it gives a shared contact with a primitive
coefficient vector: its full norm divides the determinants of all
rows in that cluster. The correction trimming from the norm-lift note
applies to these patterns as well, with `m` fixed.

For the lattice made from the shared contacts alone, therefore,

```text
log|Delta|=(4 A_m-m-1)w+o(w),
log||x||_2=A_m w+o(w).
```

The determinant argument gives, for every nonunit alternative,

```text
log||x'||_2 >= (3 A_m-m-1)w-o(w).                     (12)
```

For two allocations on the same rational contacts with short frames,
their common oriented modulus has norm weight at most `4 A_m w+o(w)`.
Each allocation must consequently place at least

```text
(4 A_m-2m-2)w-o(w)                                   (13)
```

of its norm weight outside that overlap. These become nontrivial
uniqueness and positive-disagreement assertions for `m>=5`. The
relative disagreement threshold tends to one half as `m` grows;
it remains below one half. The fixed-`m` error notation here is not
a uniform error estimate when `m` varies with `w`.

If all private block orientations are prescribed as well, they too
may be included as congruences, using the corresponding individual
row vector. Then `log|Delta_all|=(4 A_m-1)w+o(w)`, and (12)--(13)
improve respectively to `(3 A_m-1)w-o(w)` and
`(4 A_m-2)w-o(w)`. This stronger conclusion requires the private
factors as input; the rational moduli contacts do not prescribe them.
Neither version produces a second frame or excludes the first one.
