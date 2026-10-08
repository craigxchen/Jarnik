# A Chow-form basepoint criterion at the first smaller cuts

At a first smaller cut, sufficiently short exact numerical-zero relations
have restriction hyperplanes with a common point on a fixed projective
variety. The point may lie on the boundary or on the image of an exceptional
divisor. The assertion does not make the restrictions vanish and gives no
uniform bound for lattice points on endpoint arcs.

The arithmetic input is the prime-power congruence and complement reflection
in [the all-cut hierarchy](all_cut_invariant_relation_hierarchy.md), Section 3.
We use the no-loss quadratic basis from
[the first smaller rank note](smaller_cut_relation_rank.md), Section 1.
The new step uses the Chow form of the **image** of the matching quadrics.
Taking a raw resultant of those quadrics would instead retain their
universal frame basepoints.

## 1. The matching image and its degree

Let `m=2q>=8`, and set

```text
n=q+1,       d=n-2,       t=n(n-3)/2,
U=Q^n/Q(1,...,1),       P(U)=P^d.
```

The first smaller restriction space `W_(n,2)` consists of squarefree,
translation-invariant quadrics in `z_1,...,z_n`. It is the complete space
of quadrics on `P(U)` vanishing at the `n` frame points `[e_i]`.
Indeed the full quadratic space has dimension `n(n-1)/2`, and the `n`
frame evaluations are independent, leaving dimension `t`. Every matching
product

```text
(z_a-z_b)(z_c-z_e),       a,b,c,e distinct,                 (1)
```

belongs to this space. The integral basis `q_ij` of the first smaller
rank note spans it and has no coefficient conversion loss; each basis
element is an integer sum of at most two matching products. Fix this
basis and call its coordinates `w_B`.

The rational matching map

```text
phi:P(U) --> P^(t-1),       z |--> [w_B(z)]                  (2)
```

has precisely the `n` frame points as basepoints. To see this, partition
the labels into classes of equal `z` values. If no two disjoint unequal
pairs exist, the complete multipartite graph of unequal pairs has matching
number at most one. For `n>=5`, this means that one class has `n-1`
labels and the last class has one; the all-equal vector is zero in `U`.
These are exactly the frame points. At `[e_i]`, products with one edge
incident to `i` have linear parts `z_a-z_b` on the other labels. They
span the tangent cotangent space, so one blowup at each frame point
resolves the base ideal. On the blowup, (2) is given by the basepoint-free
line bundle

```text
L=2H-sum_i E_i.                                           (3)
```

The map is birational onto its image `X_n`. On the open set of distinct
values, choose `a,b` outside any triple `c,d,e`; then the ratio of two
matching products with common first edge is

```text
w_(ab|cd)/w_(ab|ce)=(z_c-z_d)/(z_c-z_e).                  (4)
```

Each matching product is a linear form in the chosen `w_B` coordinates.
These ratios recover every `z_j` after fixing two values to zero and
one, which is exactly translation and projective scale. In particular,
`X_n` is integral of dimension `d`. Mixed top intersections of `H` and
the exceptional divisors vanish, while `H^d=1` and
`E_i^d=(-1)^(d-1)`. The degree of the birational image is therefore

```text
delta=deg X_n=L^d=2^(n-2)-n.                            (5)
```

For `n=5` this is the familiar cubic degree `3`. The corresponding
Chow evaluation is the cubic in four-row cofactors from
[the eight-row Segre congruence](eight_row_segre_cofactor_congruence.md).
Only (1)--(5) are used below.

## 2. A fixed integral incidence certificate

Let `Ch_n` be a primitive integral Chow form of `X_n`. It is a polynomial
in `d+1=n-1` hyperplane rows `u_0,...,u_d`, each of length `t`, of
degree `delta` in **each** row. It vanishes exactly when those
hyperplanes have a common point on `X_n` over `Qbar`.

These standard existence and multidegree statements are also recorded
in [the Chow-form construction, Propositions 8.8 and 8.11](https://people.kth.se/~dary/Chow.pdf).
The fixed integral incidence identity needed here is justified below.

There are a nonzero integer `c_n` and an integer `K_n>=0`, depending
only on `n` and the fixed coordinate basis, with the following
incidence identity for every coordinate `v_j` of the affine cone over
`X_n`:

```text
c_n v_j^K_n Ch_n(u_0,...,u_d)
    belongs to ( I_(X_n)(v), u_0 dot v,...,u_d dot v )
    in Z[u_0,...,u_d,v].                                  (6)
```

Here `I_(X_n)` means fixed integral homogeneous equations for the
rational vanishing ideal. To justify (6), first work over `Q` and
localize at `v_j`. The incidence equations solve for the `j`th
coefficient of each independent hyperplane row, since `v_j` is a unit.
The quotient is a polynomial algebra over the integral chart of
`X_n`, hence is a domain. Its projection is dense in the irreducible
Chow divisor: a generic incident tuple meets `X_n` at a point of this
chart. Consequently `Ch_n` lies in the localized incidence ideal,
not merely its radical. Clearing the finitely many chart denominators
and rational coefficients gives the common `c_n,K_n` in (6). No
effective value for these fixed integers is needed here.

The chart powers matter at primes where some matching coordinates
vanish. They cannot simply be dropped from (6).

## 3. Actual congruences and coordinate content

Use the centrally truncated actual Gaussian system of the all-cut
hierarchy. Its independent blocks `H_T` have odd split-prime support,
are coprime to each other and to all their conjugates, and obey

```text
log N(H_T)>=(1-eta)w,
P_i=K_i product_(T containing i) H_T,
log|K_i|<=sigma,       0<|t_ij|<=T.
```

Fix a first smaller inside set `S` of size `q-1`; the `n` outside
labels are temporarily `1,...,n`. For `Q` in the balanced kernel
`K_1`, let

```text
F_(S,Q)=sum_B c_(B,Q) w_B,       c_(B,Q) in Z,
sum_B |c_(B,Q)|<=C_Q.
```

The first smaller integral basis has no coefficient conversion loss;
`C_Q` is the coefficient sum of the original integer invariant.
Put `z_j=Y_j/P_j` on the outside rows and

```text
v_B=(product_outside P_j) w_B(z) in Z[i],
G_v=gcd_G(v_B:B),       v'_B=v_B/G_v,
kappa=product_outside K_j.                               (7)
```

Every matching product is nonzero for distinct actual directions.
In the chosen `q_ij` basis, `q_13=(z_1-z_n)(z_3-z_2)` is itself one
such product, so `v_13!=0` even though other basis values may vanish.
Thus `G_v` is nonzero. The vector `v'` has Gaussian
coordinate gcd one and
still lies on the affine cone over `X_n`, because its defining equations
are homogeneous. An actual numerical zero `Q` satisfies the previously
proved Gaussian divisibility

```text
H_S | kappa sum_B c_(B,Q) v_B.                           (8)
```

Remove **only** the common coordinate content and the correction
factor. Define the Gaussian divisor and its positive rational norm

```text
mathfrak_m_S=H_S/gcd_G(H_S,kappa G_v),
M_S=N(mathfrak_m_S).                                      (9)
```

Equation (8) implies `mathfrak_m_S | sum_B c_(B,Q) v'_B` for every
such relation. For any `d+1` relations, use their coefficient rows as
the `u_l` in (6). At each Gaussian prime power dividing
`mathfrak_m_S`, some coordinate `v'_j` is a unit, since the coordinates
have Gaussian gcd one. Applying that chart's identity gives

```text
mathfrak_m_S | c_n Ch_n(u_0,...,u_d).                    (10)
```

The Chow value is an ordinary integer. The divisor `mathfrak_m_S`
is conjugate-primitive, because it divides `H_S`; therefore the exact
real norm-divisibility rule turns (10) into

```text
M_S | c_n Ch_n(u_0,...,u_d)       in Z.                   (11)
```

The fixed denominator `c_n` may contain bad primes, but its logarithm
is a constant depending only on `m`. No coordinate has been inverted
globally modulo the block.

Since `G_v` divides each `v_B`, the modulus in (9) is at least the
fixed-coordinate modulus in the all-cut hierarchy. Its audited
prime-power estimate gives

```text
log M_S >= log N(H_S)-2m sigma-2log T.                  (12)
```

Apply exact complement reflection to the same invariant polynomials.
It preserves their labelled coefficient rows, while the corresponding
modulus `M'_S` is supported on the complementary block `H_(S^c)`.
The two rational norm moduli are coprime and both divide the **same**
integer on the right of (11). Keeping the reflected correction height
and its block trimming gives

```text
M_S M'_S | c_n Ch_n(u_0,...,u_d),
log(M_S M'_S)
 >= log N(H_S)+log N(H_(S^c))
       -2m(2m+1)sigma-4log T.                            (13)
```

The factor `2m(2m+1)` includes the first correction loss, the reflected
loss with height `(2m-1)sigma`, and the complementary-block trimming.

## 4. The basepoint theorem

Let `B_n` be the sum of the absolute coefficients of the fixed
integral Chow form, and put

```text
D=(n-1)delta,
C_m=|c_n| B_n (t delta)^D.                            (14)
```

Suppose every relation under consideration is an integer numerical
zero in `K_1` with `log C_Q<=h`, and assume

```text
D h+2m(2m+1)sigma+4log T+log C_m
       < 2(1-eta)w.                                     (15)
```

Then, **at each first smaller cut**, the restriction hyperplanes of
all these relations have a common point on `X_n(Qbar)`.

For the proof, fix one cut and suppose the restriction span has no
common point on `X_n`. It has a basis represented by at most `t`
of the short relations. A basepoint-free linear system on a
`d`-dimensional projective variety admits `d+1` members with empty
intersection: choose each member to avoid every component left by the
previous choices. The Chow form restricted to `(d+1)` copies of this
span is therefore a nonzero polynomial in the combination
coefficients. Its degree in each combination variable is at most
`delta`. The elementary polynomial grid lemma supplies a nonzero
value when all coefficients are chosen from `{0,...,delta}`. Each of
the resulting `d+1` integer relations has coefficient sum at most
`t delta exp(h)`. Thus the nonzero ordinary integer Chow value obeys

```text
|c_n Ch_n| <= C_m exp(D h).                            (16)
```

Equations (13), the full-profile lower bounds for both blocks, and
(15) make the positive divisor `M_S M'_S` larger than the nonzero
integer it divides, a contradiction.

The common point in this conclusion can be in the image of an
exceptional divisor or at another boundary point of `X_n`; it need
not represent distinct finite outside directions. If all restrictions
vanish, every point of `X_n` is a common point. The theorem controls
the existence of a projective basepoint, not the dimension of the
restriction space and not membership in the next kernel `K_2`.

## 5. A fixed large-row application and its limit

The existing primitive evaluation-height calculation gives, for
fixed even `m`, a full numerical relation lattice of rank `r-1` and
logarithmic covolume `A_m w+o(w)`, where

```text
r=[x^m](1+x+x^2)^m-[x^(m-1)](1+x+x^2)^m,
A_m=m 2^(m-2)-(m/2) binomial(m,m/2).                  (17)
```

Put `L=floor((r-1)/2)`. Minkowski's second theorem, with the
remaining `r-L` minima at least the `L`th minimum, supplies `L`
independent numerical-zero relations satisfying

```text
log C_Q <= A_m w/(r-L)+o(w).                           (18)
```

For `m=400`, the exact integer inequality is

```text
40 delta(n-1) A_m < r-L,       n=201.                  (19)
```

It puts the leading coefficient on the left of (15) below `1/40`.
It also puts each relation's height below the existing balanced-cut
threshold `2w-o(w)`, so these `L` relations lie in `K_1`. In a
hypothetical full-profile sequence with this fixed `m`, one has
`eta=o(1)` and `sigma,log T=o(w)`; as `w` tends to infinity the
fixed `log C_m` is negligible. Hence their first smaller restriction
hyperplanes have a common point on `X_201` at every cut.

This conclusion is compatible with all `L` restrictions vanishing,
and the basepoints can differ from cut to cut. Neither (15) nor (19)
produces a contradiction or improves the general growth bound
for lattice points on endpoint arcs.

The later [coherent boundary model](coherent_boundary_basepoints_of_numeric_relations.md)
shows that actual numerical vanishing, large dimension and simultaneous
compatibility alone can allow boundary points. The
[quartet height descent](quartet_boundary_height_descent.md) excludes
that model's nonzero quartet-supported extension at the required short
height. It does not classify every possible short relation space.

## Verification

Run `python3 -B docs/check_first_smaller_cut_chow_basepoint.py`. The
checker verifies the finite matching-space dimension and basepoint
examples over a small field, the birational coordinate-ratio identity
on distinct integer samples, the blowup intersection degree formula
for small `n`, and inequality (19) by exact integer arithmetic. It
does not compute the Chow form, the fixed incidence certificate, or
the actual Gaussian prime-power moduli; those are proved above and in
the cited hierarchy.
