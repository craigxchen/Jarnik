# Odd supports, multiplication, and five-row relation rigidity

The smaller-cut hierarchy extends to odd row counts, starting with
linear restrictions rather than balanced scalars. On five rows this
gives a new actual arithmetic consequence: two independent sufficiently
small numerical invariant relations cannot coexist. The proof uses
the ten boundary lines of the degree-five del Pezzo surface, together
with a coefficient bound excluding proper numerical-zero factors.

This is a rank-one conclusion, not a rank-zero conclusion. A single
small irreducible five-row relation is still not excluded. Its ideal
of multiples retains a nonconstant factor in every surviving
restriction, so multiplying it by fixed bracket or triangle invariants
does not supply the missing rank-to-zero implication.

The exact linear algebra and multiplication checks are in
[check_odd_support_relation_rigidity.py](check_odd_support_relation_rigidity.py).
That script checks only the displayed finite rank, incidence, and capacity
certificates; it does not prove the arithmetic divisibility estimates or
the surface intersection argument.
The geometric dictionary is established in
[five_row_del_pezzo_arithmetic.md](five_row_del_pezzo_arithmetic.md).

## 1. The odd hierarchy and its actual integer moduli

Let `m=2q+1>=3`, and let `V_m` be the simultaneous `SL_2` invariant
space of degree two in each row. Set

```text
K_k={Q: every restriction with |S|>=q-k+1 vanishes}.
```

As before, rows in `S` are specialized to `(0,1)` and the outside
rows to `(1,z_j)`. Weight balance makes `K_0=V_m`. If `Q in K_k`
and `|S|=q-k`, its restriction is squarefree, translation invariant,
and homogeneous of degree

```text
d=2k+1,             n=q+k+1 outside variables.
```

Indeed its total second-coordinate weight is `m`, the inside rows
account for `2|S|`, and a coefficient of `z_j^2` is a restriction
at the next larger cut. The same squarefree raising/lowering argument
as in the even hierarchy gives target dimension

```text
t=binomial(n,d)-binomial(n,d-1)   if 2d<=n,
t=0                             if 2d>n.              (1)
```

Thus a nonzero level requires `3k+1<=q`. Downward propagation to
the empty cut proves

```text
K_(floor((m-3)/6)+1)=0.                                 (2)
```

This termination is sharp: write `m=6ell+3+2r`, `r in {0,1,2}`,
and take `2ell+1` disjoint triangles and `r` disjoint squared
brackets. At most one row can be collapsed in each component without
vanishing. Collapsing exactly one in each component gives a nonzero
restriction at size `q-ell`, so the product lies in `K_ell` but not
`K_(ell+1)`. It is nonzero at every distinct actual configuration.

Now use the actual conjugate-primitive Gaussian rows and independent
odd split blocks of the earlier hierarchy:

```text
P_i=K_i product_(T containing i)H_T,
log |K_i|<=sigma,       0<|t_ij|<=T,
log N(H_T)>=(1-eta)w.
```

The integral ballot basis at degree `d` consists of products of `d`
differences on `2d` distinct outside rows. Let its fixed coordinate
conversion bound be `L=2^(t-1)`. Clearing the common outside
denominator leaves each basis value as a product of `d` nonzero
brackets and `n-2d` unused rows. The exact adjugate argument therefore
gives a common integer modulus for every maximal coefficient minor,
with one-pass loss

```text
2n sigma + 2(n-2d)sigma + 2d sigma + d log T
       = 2m sigma + d log T.
```

For clarity, the odd-row congruence used here is the same calculation,
with no scalar balanced-evaluation assumption. After the common
determinant-one change to `(P_i,Y_i)`, reduce at a cut `S` modulo
`H_S`. Terms with an internal edge in `S` vanish, while every surviving
term is
`(product_(i in S)Y_i^2)(product_(j outside S)P_j^2)F_(S,Q)(Y/P)`.
Writing `F_(S,Q)=sum_B c_B w_B` and
`v_B=(product_outside P_j)w_B(Y/P)`, cancellation of the inside
`Y_i` factors and of the outside core factors (which are units modulo
`H_S`) gives
`H_S | kappa sum_B c_B v_B`, with
`kappa=product_(j outside S)K_j`.
The distinguished `v_B` is a product of `d` outside brackets and the
remaining outside rows. The local pair-residue bound gives the displayed
loss, and an integer adjugate gives the divisibility for every maximal
minor. The reflected rows give the complementary modulus by the same
calculation. Thus (3) is an odd-row lemma, rather than a consequence of
the even-row scalar theorem by formal analogy.

Complement reflection preserves the same labelled coefficient matrix.
Its complementary modulus is coprime to the first, and hence

```text
M_S M_S' | every maximal coefficient minor,
log(M_S M_S') >=2(1-eta)w
                    -2m(2m+1)sigma-2d log T.           (3)
```

No new primitivity assertion is used here: the reflection lemma
already divides ordinary contents and charges the common core trim.
Zero imaginary coordinates are allowed by that lemma. Equations
(1)--(3) extend the actual integer restriction theorem, but (3)
still says only that the image is proper when all maximal minors are
too small.

For `m=5`, the first restrictions have `|S|=2`, `n=3`, `d=1`,
and `t=2`. Use the integral basis `z_a-z_c, z_b-z_c`. Its two
coefficients have sum of absolute values at most `C_Q`, the original
polynomial coefficient sum. Thus two numerical zeros satisfy

```text
M_S M_S' | delta_S,
|delta_S|<=C_Q C_R,
log(M_S M_S') >=2(1-eta)w-110sigma-2log T.              (4)
```

Here `delta_S` is the determinant of their two linear restriction
coefficient rows. These are restrictions of anticanonical sections
to the ten boundary lines, in a common fixed trivialization.

## 2. Every sufficiently small five-row zero is irreducible

Put `F_5=(3^5-1)log 2=242log 2`. The factor theorem gives a
primitive irreducible numerical-zero factor `G` of `Q` with

```text
log C_G<=log C_Q+F_5.                                  (5)
```

Proper factors with at most four supported rows are excluded below
`4(1-eta)w-8sigma` by support projection and the four-row theorem.
A proper factor involving all five rows has, up to permutation, one
of the two multidegrees

```text
(2,1,1,1,1),              (2,2,2,1,1).                (6)
```

For the first pattern, the invariant space and its weighted-balanced
map both have rank three. For the second, their ranks are four and
three; the scalar kernel is exactly

```text
Delta_12 Delta_23 Delta_31 Delta_45.                    (7)
```

The certificate enumerates the bracket monomials, verifies their
polynomial ranks using exact rational elimination, and computes the
weighted-balanced ranks. The standard weight multiplicity formula
gives dimensions three and four, so the checked spanning families
are complete. Expression (7) is nonzero at all distinct directions.

For these mixed degrees the complementary weighted-balanced scalar
bound is

```text
log C_G>=2(1-eta)w-2 sum_i d_i log|K_i|
       >=2(1-eta)w-16sigma                             (8)
```

unless every weighted-balanced scalar vanishes. In the first pattern
that would make `G=0`; in the second it makes `G` a multiple of
(7), which cannot be a numerical zero. Thus (8) holds for every
nonzero numerical zero of either proper full-support pattern.

It follows that any nonzero `Q in V_5` satisfying

```text
Q(P)=0,       log C_Q+F_5<2(1-eta)w-16sigma             (9)
```

is irreducible over `Q` up to its ordinary constant content. Its zero
divisor on the del Pezzo surface is an anticanonical curve irreducible
over `Q`; geometric irreducibility is not needed below. In particular
it has no boundary component. This uses actual
numerical nonvanishing of brackets and (7), rather than presuming
that arbitrary invariant factors are nonzero.

The passage from polynomial to surface irreducibility uses the standard
quotient dictionary implicit in the del Pezzo model. On the open locus
of distinct directions, the projective `PGL_2` quotient (equivalently the
even-degree `SL_2` invariant quotient) has the usual principal-orbit
dictionary. A proper rational divisor factor of the anticanonical
section pulls back to a proper multihomogeneous invariant factor; after
clearing boundary denominators, any extra factor is a boundary bracket.
Since the section has no boundary component, polynomial irreducibility
rules this out. The same argument applied to the Galois orbit of a
geometric common component shows that two rationally irreducible,
independent sections have no common component after base change. This
is the quotient-theoretic input needed in the intersection argument.

## 3. Two small actual zeros cannot share an interior point

Let `L=-K_Y` on the degree-five del Pezzo surface `Y`. The needed
intersection facts are

```text
L^2=5,       L.D_e=1,
there are ten boundary lines D_e,
at most two boundary lines meet at any point.          (10)
```

They follow either from the blowup of `P^1 x P^1` at the three
diagonal points or from the equivalent four-point plane blowup.
The boundary graph is the Petersen graph; there are no triple
boundary intersections.

Suppose `Q,R` are independent anticanonical sections irreducible over
`Q` and vanishing at the same interior rational configuration. Their
curves have no common component even after extending the field: a
shared geometric factor brings all its Galois conjugates, yielding
a nonconstant common rational factor. Independence and rational
irreducibility exclude this. No boundary component can be hidden in
the invariant-to-surface passage, since each such component would
give a bracket factor. If every `delta_S`
in (4) vanished, their restrictions on each boundary line would be
proportional nonzero sections of a degree-one line bundle. They
would therefore have a common zero on each of the ten lines.
At most two lines can share that zero, so the two curves have at
least five distinct common boundary points. They also have their
given common interior point. Since they have no common component,
this contradicts their total intersection number five.
Consequently some `delta_S` is nonzero.

Combining this argument with (4)--(9) gives the finite arithmetic
statement

```text
Q,R independent nonzero integer numerical zeros in V_5
  => log C_Q+log C_R
        >=2(1-eta)w-110sigma-2log T-F_5.                (11)
```

To check the constants, argue by strict contradiction to (11).
Since each coefficient sum is at least one, each polynomial
separately satisfies (9). Thus both are irreducible, and (4) makes
every boundary determinant zero, contradicting the preceding
intersection argument. The factor constant is deliberately not
optimized.

In particular, all five-row numerical zeros with

```text
log C_Q<(1-eta)w-55sigma-log T-F_5/2                    (12)
```

span a space of dimension at most one. This conclusion uses the
actual common numerical zero in the interior; an arbitrary pair of
sections whose boundary restrictions are proportional need not
satisfy it. The familiar two-pentagon pencil has five boundary
basepoints and no common interior basepoint, exactly the allowed
intersection budget.

## 4. Projection strengthens the result but does not force zero

If these five rows are retained from an `m`-row full-profile system,
their block weight is `w_5=2^(m-5)w`, their original corrections
remain bounded by `sigma`, and their primitive pair residues are
unchanged. Thus the leading term in (11) is

```text
2w_5=2^(m-4)w.                                        (13)
```

For a fixed triangle or bracket monomial `B` on disjoint rows,
`B(P)!=0`. Removing it from an actual relation `Q=B G` preserves
the numerical zero, and the factor bound costs only `O_m(1)` in
logarithmic coefficient height. Therefore two independent short
residual five-row factors in such a common-factor subspace are
excluded at the amplified threshold (13). This sharpens the bound
on its short dimension to one. It does not eliminate a single
residual factor.

There is an exact reason that multiplication alone does not make
its next restrictions vanish. For any nonzero invariant `F` on a
labelled row set, with the ambient set taken to be exactly its supported
rows, define its collision capacity by

```text
c(F)=max{|S|: the polynomial restriction F_S is nonzero}.
```

For disjoint row supports,

```text
c(FG)=c(F)+c(G).                                       (14)
```

The restriction factors, and the polynomial rings on the disjoint
outside variables have no zero divisors. This proves both directions
of (14). A triangle and a squared bracket each have capacity one.
Every nonzero five-row degree-two invariant has capacity two:
larger cuts vanish by weight balance and `K_1=0` in (2) excludes
vanishing at all two-row cuts.

Thus, if a nonzero small actual five-row relation `G` exists, its
product with a fixed disjoint triangle has capacity three. On eight
rows this product lies in the balanced kernel `K_1` but does not lie
in `K_2`. Its size differs from `G` by only a fixed multiplicative
constant. The improved five-row theorem permits one such `G`; it
does not prove its existence or exclude it. Hence promoting all
short common-factor relations to zero next restrictions would
already require a theorem excluding this single irreducible
five-row possibility.

The successive-minima budget is consistent with this residual gap.
The primitive five-row evaluation vector has height `10w_5+o(w_5)`
and its numerical relation lattice has rank five. Formula (11)
requires only

```text
log lambda_1+log lambda_2>=2w_5-o(w_5).
```

One bounded first minimum and four minima of exponent `5/2` have
total exponent ten and respect this inequality. This is a compatible
height budget, not a claimed actual elliptic or lattice construction.
The first minimum estimate supplied by the determinant is `2w_5`,
so it neither contradicts (11) nor supplies a second vector below
the threshold (12).

The odd hierarchy and the five-row pair theorem are therefore new
constraints on actual integer relations. A global rank-to-zero step,
or a positive lower coefficient exponent for every individual
irreducible five-row relation, remains unproved.
