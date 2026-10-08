# A five-swap cross-ratio forces one nearby root, without a height gain

The five Vieta swaps from **one** anticanonical section satisfy an
exact projective identity.  It forces at least one swapped direction
to stay within a constant multiple of the old five-point angular
span.  On the actual full-cut source, that proximity gives a primitive
residue allowance of about `exp(9w/4)`, weaker than the independent
swap bound `exp(2w)`.  The formal cut exponents of the identity also
leave multiple terms without a forced core factor at every cut, so degree counting alone
cannot force one exceptional content `E_i` to grow.

The universal identity and its converse are proved in
[five_row_simultaneous_vieta_identity.md](five_row_simultaneous_vieta_identity.md)
and certified symbolically in
[check_five_row_simultaneous_vieta_identity.py](check_five_row_simultaneous_vieta_identity.py).
The cut and rational-fixture audit for the consequences here is
[check_five_row_vieta_crossratio_tropical.py](check_five_row_vieta_crossratio_tropical.py).
Neither checker constructs an endpoint family.

## The same-section identity

Let `Q` be one nonzero five-row integer simultaneous `SL_2` invariant
of degree two in every row.  At a zero on five distinct real
projective directions, assume every row-gradient multiplier is
nonzero, so the five Vieta swaps are distinct from their old roots.
Write the old and new homogeneous binary vectors as `v_i,v_i'`.
Then the **projectively invariant** identity is

```text
sum_(i<j) det(v_i,v_j) det(v_i',v_j')
          /[det(v_i,v_i') det(v_j,v_j')] = 6.           (1)
```

It remains valid when a new root is at infinity and under independent
nonzero rescaling of every old or new vector.  This is a relationship
among all five swaps of the same `Q`, not five separate estimates.

For a proof in a finite affine chart, put `v_i=(1,x_i)`, let
`lambda_i=Q_(x_i)` and `a_i=Q_(x_i x_i)/2` at the actual zero, and put
`f(t)=product_j(t-x_j)`.  The exact `SL_2` Ward identities imply

```text
lambda_i f'(x_i)=A x_i+B                              (2)
```

for some nonzero pair `(A,B)`.  Define
`u_i=a_i/lambda_i` and `U_r=sum_i x_i^r u_i`.  The fixed cubic Hessian
certificate, checked on all six integral graph-basis polynomials in
the companion checker, implies at `Q=0`

```text
A U_1+B U_0=3A,
A U_2+B U_1=(sum_i x_i)A+2B.                     (3)
```

Eliminating `(A,B)` gives

```text
(U_1-3)(U_1-2)=U_0(U_2-sum_i x_i).               (4)
```

The affine Vieta root is `x_i'=x_i-1/u_i` when `a_i!=0`; the
homogeneous version covers `a_i=0`.  Expanding the left side of (1)
gives `U_0U_2-U_1^2+5U_1-(sum x_i)U_0`, which equals six by (4).
The checker verifies the Hessian certificate as thirty exact
polynomial identities and also exercises the projective infinity
case on rational fixtures.

Two immediate consequences require no sign assumption.  Some pair
term in (1) has absolute value at least `3/5`, so for one pair

```text
|det(v_i',v_j')|
 >=(3/5)|det(v_i,v_i')det(v_j,v_j')|
          /|det(v_i,v_j)|.                              (5)
```

This is a new-new separation relative to the two self swaps.  It
does not bound their primitive residues, since the new-new determinant
has no inherited small-residue upper bound.

The identity is also exactly blind to primitive normalization.  With
the actual integer exceptional contents of the one-row Vieta audit,
put `w_i=E_i v_i'`.  Each term of (1) can be rewritten as

```text
det(v_i,v_j)det(w_i,w_j)
 /[det(v_i,w_i)det(v_j,w_j)],                    (5a)
```

because its numerator and denominator both gain `E_iE_j`.  Clearing
the self-bracket denominators gives

```text
sum_(i<j) det(v_i,v_j)det(v_i',v_j')
           product_(k!=i,j) det(v_k,v_k')
 =6 product_k det(v_k,v_k').                          (5b)
```

The same cancellation persists there after replacing each `v_i'`
by `w_i`.  Thus the projective identity and all its real signs can
select a new **direction**, but by themselves cannot select an index
with large `E_i`; they carry no primitive coefficient-content data.
Using the common `Q` to force `E_i` requires additional integer
coefficient or Gaussian-anchor arithmetic beyond (1).

## A nearby swap, even if another root is at infinity

Choose a common real affine chart in which the five old roots are
finite, and let

```text
D=max_(i<j)|x_i-x_j|,       u_i=a_i/lambda_i.
```

The affine summand of (1) is
`(x_i-x_j)^2u_i u_j+(x_j-x_i)(u_j-u_i)`, valid also when
`u_i=0` and the corresponding new root is at infinity.  If every
`|u_i|<=1/(4D)`, its absolute value is at most
`1/16+1/4+1/4=9/16`.  Ten such terms sum in absolute value to at
most `90/16<6`, a contradiction.  Hence one `|u_i|>1/(4D)`; that
new root is finite and satisfies

```text
min_i |x_i-x_i'|<4D.                                   (6)
```

For an old circle arc with phase span `Delta<pi/4`, choose the
stereographic coordinate
`x=tan((theta-theta_left)/2)` with pole opposite the arc, where
`theta` is the phase of `P/bar(P)` in a consistent branch.  The old
coordinates lie in `[0,tan(Delta/2)]`, so (6) puts one new coordinate
in `[-4tan(Delta/2),5tan(Delta/2)]`.  Its exact phase lies within

```text
2 arctan(5 tan(Delta/2))=O(Delta)                    (7)
```

of the old arc's left endpoint.  If the chosen pole coincides with a
new direction, the `u_i=0` version above still yields a finite nearby
root.  No assertion is made that
the near swap lies inside the old arc or that its Gaussian denominator
divides the old source anchor.

The single-near conclusion cannot be increased to two near swaps using
the projective identity.  For an exact rational family, take old
coordinates `x=(0,1,2,3,4)`, set four shifts
`k_2=k_3=k_4=k_5=T`, and put

```text
k_1=t(T)=(-10+30/T)/(6-10/T-20/T^2).                  (6a)
```

The ten-term left side of (1) is exactly
`(30/T-10)/t+10/T+20/T^2=6`.  As rational `T` tends to infinity,
`t(T)` tends to `-5/3`, while the other four swaps escape every fixed
affine neighborhood of the old cluster.  The moment-kernel generator
in the [converse reconstruction](five_row_simultaneous_vieta_identity.md)
tends to `(A,B)=(1,-5)` and has `A x_i+B!=0` for all five old rows
when `T` is sufficiently large.  Therefore each such rational tuple
comes from one rational invariant section, scalable to a primitive
integer section.  The checker reconstructs the same section at
`T=100,1000,10000` and verifies every other root exactly.  This is a
sharp geometric barrier for (1).  At `T=100` the primitive graph-basis
coefficient vector is

```text
(-19948,3804,-2018,-1843,4038,-59).
```

Fixing the first three projective rows at `0,1,2`, the normalized
two-variable section is `F(x,y)=A(x)y^2+B(x)y+C(x)`, with

```text
A=35388-5871x-1961x^2,
B=-71712-18906x+15558x^2,
C=79320x-31816x^2.
```

Its discriminant modulo eleven is
`B^2-4AC=9+8x+2x^2+5x^3+x^4`, whose gcd with its derivative is one.
The exact quartic is therefore squarefree over the rationals and is
not a square in `C(x)`; `F` is geometrically irreducible.  The checker
also evaluates all ten labelled boundary restrictions and finds all
three outside coefficients nonzero in each, so this particular
section has no boundary component or boundary-node contact.  The
reconstructed coefficients, normalized discriminant, and boundary
restriction coefficients are rational functions of `T`: the old-source
matrix in the converse is fixed, while `t(T)` and the proposed new
directions are rational in `T`.  Since the discriminant squarefreeness
resultant and every boundary coefficient are nonzero at `T=100`, none
is the zero rational function.  They can fail at only finitely many
values of `T`.  Consequently **arbitrarily large rational `T`** in
this family still give geometrically irreducible, boundary-nonnodal
integer sections with exactly one bounded shift and four escaping
shifts.  It is
not asserted to have small coefficients relative to inherited full-cut
heights, actual full-cut source factors, or an endpoint arc.

In the actual five-row full-cut source, the Gaussian anchor `z_0`
contains one conjugate copy of all thirty-one nonempty odd core
blocks.  Thus `R=|z_0|>=exp(31(1-eta)w/2)`, and an endpoint arc has
`Delta<=C/sqrt(R)`, with leading exponent `-31w/4`.  The one-row
[Vieta audit](five_row_vieta_swap_content_and_residue_audit.md)
allows `|P_i'|` of order `C_Q exp(13w)/E_i`; a retained row has
modulus of order `exp(8w)`.  Equation (7) consequently gives a
new-old bracket allowance of order

```text
|det(P_i',P_j)| <= O(C C_Q) exp((53/4+O(eta))w
                                  +O(beta+sigma))/E_i.  (8)
```

After the transferred eleven-cut Gaussian gcd, the primitive
new-old residue allowance is about `exp(9w/4)`.  It is weaker in its
main exponent than the exact `exp(2w)` bound from evaluating
`q_i(P_j)`; near-angle control therefore does not reduce the residue
cost.  The anchor multiplier bound of order `exp(5w)/sqrt(E_i)` is
also unchanged by (7).

## An original-anchor imaginary-coordinate upper bound on content

Keep the original Gaussian anchor `z_0` and its Cartesian rows
`P_j=X_j+iY_j`; an arbitrary common rotation would change the
integer coefficients in the next identity.  Put `S=max_j|Y_j|`.
For the row quadratic

```text
q_i(X,Y)=a_iX^2+b_iXY+c_iY^2,
```

the canonical primitive second root in the
[one-row Vieta audit](five_row_vieta_swap_content_and_residue_audit.md)
satisfies the **signed exact** factorization
`q_i(W)=g_i det(P_i,W)det(P_i',W)`.  Comparing the `X^2`
coefficient gives

```text
a_i=g_iY_iY_i'.                                        (8a)
```

In every integral graph-basis monomial, row `i` has two incident
brackets.  Evaluating that monomial at row `i=(1,0)` replaces those
brackets by `±Y_j` and `±Y_k`, while its three remaining retained-row
brackets have absolute value at most
`exp(8(1+eta)w+beta)` each.  The integral basis conversion
`H_G<=2C_Q` therefore yields

```text
|a_i|<=2C_Q S^2 exp(24(1+eta)w+3beta).                 (8b)
```

Actual full-cut rows have `Y_i!=0`: a primitive real-axis row would
be a unit, whereas its sixteen inherited nonunit blocks make
`|P_i|>1`.  If the new root also has `Y_i'!=0`, then its ordinary
integer coordinate has `|Y_i'|>=1`.  From
`E_i=g_i Norm(K_i)/F_i`, `F_i>=exp(19(1-eta)w)`, and
`Norm(K_i)<=exp(2sigma)`, equations (8a)--(8b) give

```text
E_i <= (2C_Q S^2/|Y_iY_i'|)
       exp((5+43eta)w+3beta+2sigma)
    <= (2C_Q S^2/|Y_i|)
       exp((5+43eta)w+3beta+2sigma).                   (8c)
```

If `Y_i'=0`, Cartesian primitivity makes `P_i'=(±1,0)`.
Its exact phase `P_i'/bar(P_i')=1` is the original anchor point
`z_0`, so this is an anchor collision and (8c) does not apply.
If `log S+log C_Q+beta+sigma+eta w=o(w)`, (8c) limits
`E_i` to `exp((5+o(1))w)`.
The endpoint angle bound alone gives only
`log S<=(1/4+63eta/4)w+sigma+log(C/2)` from
`|Y_i|<=|P_i|Delta/2`, `|P_i|<=exp(8(1+eta)w+sigma)`, and
`R>=exp(31(1-eta)w/2)`.  Thus it does not make `log S=o(w)`
without an extra arithmetic input.
It is an **upper** bound and supplies neither the large exceptional
content needed to shrink the anchor multiplier nor a radius descent.
The checker verifies (8a) coefficient by coefficient on integer
same-section fixtures and includes a real-axis second-root fixture;
it does not test the full-cut height hypotheses of (8b)--(8c).

## Why a cut valuation cannot force exceptional content

The formal saturated cut pattern of
[five_row_simultaneous_vieta_cut_patterns.md](five_row_simultaneous_vieta_cut_patterns.md)
has old incidence bit `b_i=1(i in U)` and new selected bit

```text
b_i'=1[(i outside U and |U|>=2)
       or (i inside U and |U|>=4)].
```

Delete all exceptional-content, reduced-multiplier, row-correction,
pair-residue, and accidental new-new determinant valuations from this
formal bookkeeping.  The listed core blocks alone give the
self-bracket exponent
`d_i=1(i in U and |U|>=4)`.  The core exponent of a term of (1) is

```text
e_ij=b_i b_j+b_i'b_j'-d_i-d_j.                        (9)
```

The exact enumeration of all thirty-one nonempty cuts is

| Cut size | Terms with exponent zero | Terms with exponent one |
| ---: | ---: | ---: |
| 1 | 10 | 0 |
| 2 or 3 | 6 | 4 |
| 4 or 5 | 10 | 0 |

Every cut has at least six terms with **no formally forced core
factor**.  There is no positive common power arising from these block
degrees on the left of (1) that could contradict `6` or compel extra
valuation in some `E_i`.  The table is deliberately a formal count:
actual `E_i`, row corrections, unselected conjugate support, and
new-new determinant residues can change individual valuations, and
even all zero-exponent terms could be divisible by the cut prime for
value-level reasons.  It rules out only a forced-degree contradiction
from (1); sign/order information and simultaneous congruences remain
open inputs.
