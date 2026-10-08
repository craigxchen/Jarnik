# An exact Vieta swap of one five-row invariant root

A small five-row numerical-zero invariant gives an exact rational
second root when one primitive row varies.  The swap carries substantial
core divisibility, but the resulting height and primitive pair-residue
budgets are too large to prove an endpoint descent. In particular,
these estimates do not establish preservation of the source radius,
the short arc, or the small-residue regime. This note keeps the ordinary coefficient
content, Gaussian core factors, and source circle anchor separate.

The identities are checked on literal full-cut Gaussian tuples by
[check_five_row_vieta_swap_content_and_residue.py](check_five_row_vieta_swap_content_and_residue.py)
and, independently, by
[check_five_row_vieta_literal_core_transfer.py](check_five_row_vieta_literal_core_transfer.py).
Those fixtures use invariants of unrestricted height and are not
endpoint examples.

## 1. Primitive second root and all pair brackets

Use the actual five-row source hypotheses (2) of
[five_row_gradient_discriminant_core_count.md](five_row_gradient_discriminant_core_count.md).
Let `Q` be a nonzero integer simultaneous `SL_2` invariant of degree
two in each row, with `Q(P_1,...,P_5)=0`, and let `C_Q` be its polynomial
coefficient sum.  Fix a row `i`; with the other four rows fixed, write

```text
q_i(X,Y)=a_i X^2+b_i XY+c_i Y^2,       q_i(P_i)=0.
```

Assume `q_i` is not identically zero.  This is automatic for an
irreducible anticanonical `Q`: a whole forgetful fiber in its zero curve
would be a component.  All displayed swap claims are scoped to such a
row.

Set `g_i=gcd(a_i,b_i,c_i)>0`, the **ordinary integer** content.
Choose an integer frame `u_i` with `det(P_i,u_i)=1`, possible because
`gcd(X_i,Y_i)=1`.  In the coordinates `W=tP_i+s u_i`, Euler's identity
and [the integral-gradient identity](five_row_gradient_discriminant_core_count.md)
give

```text
q_i(tP_i+s u_i)=s(lambda_i t+C_i s),
C_i=q_i(u_i),       g_i=gcd(lambda_i,C_i).              (1)
```

The second primitive root is the Cartesian integer vector

```text
P_i'=(C_i P_i-lambda_i u_i)/g_i,
gcd(X_i',Y_i')=1,
det(P_i,P_i')=-lambda_i/g_i.                            (2)
```

If `lambda_i=0`, the two roots coincide projectively.  When `Q` lies
below the ramification threshold in
[five_row_singular_point_height.md](five_row_singular_point_height.md),
every `lambda_i` is nonzero and (2) is a genuine swap.  A
different determinant-one frame `u_i+kP_i` changes `C_i` to
`C_i+k lambda_i` and leaves the numerator in (2) unchanged.  Swapping
the primitive roots again returns `P_i` up to sign.  This is a fiberwise
involution and does not assert that `Q` is globally reducible.

For every retained row `j`, let `Delta_ij=det(P_i,P_j)` and
`Delta'_ij=det(P_i',P_j)`.  The exact binary-form factorization in the
frame is

```text
q_i(W)=g_i det(P_i,W)det(P_i',W),
q_i(P_j)=g_i Delta_ij Delta'_ij.                         (3)
```

Thus `Delta'_ij=q_i(P_j)/(g_i Delta_ij)` when the original directions
are distinct.  The primitive new residue is the **integer**

```text
t'_ij=Delta'_ij/Norm(gcd_G(P_i',P_j)).                   (4)
```

No old `t_ij` is silently reused.  The new residue can be zero if the
swapped row collides projectively with `j`; a small-height collision
exclusion is proved below.  Cartesian primitivity of `P_i'`
makes its odd Gaussian support conjugate-primitive; a common ramified
`1+i` can remain and must be tracked when reducing its phase.

## 2. The on-zero coefficient content has weight nineteen

For a labelled cut `U`, put `a=|U\{i}|`, `n_U=Norm H_U`, and
`m=v_pi(H_U)` at one of its odd split Gaussian primes.  Use the
determinant-one Gaussian chart `(P=X+iY,Y)`, so that

```text
q_i=A_i P^2+B_i P Y+C_i^* Y^2.                         (5)
```

Here `A_i=a_i`, `B_i=b_i-2ia_i`, and
`C_i^*=c_i-ib_i-a_i`.  The triangular change is unimodular over
`Z[i]`, so these three coefficients generate the same Gaussian ideal
as the ordinary coefficients and its valuation is `v_pi(g_i)`.
The simultaneous invariant has first-coordinate weight five in this
chart.  Among the other four rows, each of the `4-a` rows outside `U`
can contribute at most two first coordinates.  Consequently

```text
v_pi(A_i)>=max(0,2a-5)m,
v_pi(B_i)>=max(0,2a-4)m,
v_pi(C_i^*)>=max(0,2a-3)m.                              (6)
```

If `i` is inside `U`, the first bound is the content bound.  If `i`
is outside `U` and `a>=3`, the actual root improves it.  Put
`t=v_pi(K_i)=v_pi(P_i)`.  From
`A_i P_i^2+B_i P_iY_i+C_i^*Y_i^2=0`, the lowest term must cancel.
Using only `v_pi(Y_i)>=0`, (6), and the baseline bound gives

```text
v_pi(A_i)>=(2a-4)m-min(t,m).                           (7)
```

For `a<=2`, the baseline already equals the required nonnegative
bound.  The correction prime may share the core orientation, which is
why `min(t,m)` is retained.  Define

```text
h_i(U)=max(0,2a-5)       if i in U,
       =max(0,2a-4)       if i outside U,
F_i=product_U n_U^h_i(U).                              (8)
```

The disjoint odd split supports make the primewise conclusions
multiplicative.  Because `g_i` is an ordinary integer and
`Norm(K_i)` charges every lost same-orientation correction valuation,

```text
F_i divides g_i Norm(K_i),
E_i=g_i Norm(K_i)/F_i in Z_>0,
sum_U h_i(U)=19,
g_i>=exp(19(1-eta)w-2sigma).                            (9)
```

The nineteen consists of twelve from the eight cuts with `a=3`
and seven from the two cuts with `a=4`.  This improves the three
nonincident-bracket baseline weight fourteen by five, but the
`Norm(K_i)` charge is necessary.

## 3. Transferred core and the surviving residue cost

In the chart (5), the exact factorization (3) gives

```text
C_i^*=g_i P_i P_i'.                                    (10)
```

Select a cut when `(i outside U, a>=2)` or `(i inside U, a>=3)`.
For each selected `U`, subtracting the content valuation in (9) and
the old `P_i` valuation from the last bound in (6) proves

```text
H_U divides E_i P_i' in Z[i].                          (11)
```

There are exactly sixteen selected cuts, and exactly eleven contain
each retained row `j`.  Let `Gamma'_ij` be the product of those eleven
blocks.  It divides both `E_i P_i'` and `P_j`.  If
`h=gcd_G(Gamma'_ij,P_i')`, then
`Gamma'_ij/h` divides the ordinary integer `E_i`.  It uses only one
Gaussian orientation at each split norm-prime, so
`Norm(Gamma'_ij/h)` divides `E_i`, not `E_i^2`.  Hence

```text
Norm(gcd_G(P_i',P_j))
       >=Norm(Gamma'_ij)/E_i
       >=exp(11(1-eta)w)/E_i.                           (12)
```

The determinant `Delta_ij` contains eight inherited cut norms and
has `|Delta_ij|>=exp(8(1-eta)w)`.  Evaluating `Q` with row `i`
replaced by `P_j` evaluates five pair brackets among the original
five rows; every graph value is at most
`exp(40(1+eta)w+5beta)`.  The integral graph coefficient conversion
is at most `2C_Q`.  Equations (3), (9), and (12) therefore yield

```text
|Delta'_ij|<=2C_Q exp((13+67eta)w+5beta+2sigma)/E_i,
|t'_ij|    <=2C_Q exp((2+78eta)w+5beta+2sigma).        (13)
```

The unknown extra content `E_i` cancels in the primitive residue
bound.  The surviving `2w` coefficient is much larger than the
subexponential original residue budget.  The exact transfer (11)
therefore does not establish that the new row belongs to the same
small-residue source class. It also gives no bound on the new angle
relative to the old short arc. The generic binary majority pattern has sixteen selected
blocks and eleven shared with each retained row, but (11) only
transfers them after multiplying by `E_i`; it does not assert exact
new prime allocations.

For **irreducible** `Q`, a collision with a retained row is separately
excluded below a `4w` coefficient threshold.  If `P_i'` has the same
direction as `P_j`, Cartesian primitivity makes it `±P_j`, and (3)
gives `q_i(P_j)=0`.  Put the duplicated direction at infinity and the
other three retained directions at `0,1,t`.  The collision restriction
of `Q` is a linear polynomial `R(t)`; it is nonzero, since an
identically zero restriction would make the bracket `Delta_ij` a
factor of `Q`. After permuting the duplicated labels into the first
two positions, the integral six-graph basis has restrictions
`0,0,0,0,t-1,t`. Permutation preserves `C_Q`. Thus the coefficient sum is at most
`2H_G<=4C_Q`.

The actual four retained rows have three matching coordinates, each a
product of two pair brackets.  Each has raw core norm-logarithm at
least `16(1-eta)w`.  Their common core exponent is exactly twelve;
any extra ordinary gcd divides the product of their six bracket
corrections, of size at most `exp(6beta)`.  Their primitive projective
height therefore satisfies

```text
h_matching>=(4-28eta)w-6beta.
```

In the `infinity,0,1,t` chart, the three matching coordinates are
projectively `1-t,-t,-1`; their projective logarithmic height is at most
`h(t)+log2`. A nonzero linear `R(t)` vanishing at `t`
gives `h(t)<=log(4C_Q)`.  Consequently a collision forces

```text
log C_Q>=(4-28eta)w-6beta-log 8.                    (13a)
```

Below this threshold all four new pair residues are nonzero.  The
upper bound (13) still permits them to be about `exp(2w)`, so this
collision exclusion does not complete the descent.

There is a corresponding norm allowance.  A graph coefficient of
`q_i` is a product of three other-row determinants and two coefficients
from incident brackets, so its absolute value is at most
`4C_Q exp(40(1+eta)w+3beta+2sigma)`.  For two linear forms whose
product is `q_i/g_i`, the coefficient vector has Euclidean norm at
least `|P_i||P_i'|/sqrt(2)`.  Since
`|P_i|>=exp(8(1-eta)w)`, (9) gives

```text
|P_i'|<=4sqrt(6) C_Q
         exp((13+67eta)w+3beta+4sigma)/E_i.             (14)
```

The old row has modulus about `exp(8w)`; (14) permits the new row to
be about `exp(13w)`.  It does not force a smaller row or a fixed norm.

## 4. The source anchor and the exact radius cost

The original integer circle is represented by a **Gaussian anchor**
`z_0` of radius `R=sqrt(Norm z_0)`, with

```text
z_j/z_0=P_j/bar(P_j),      bar(P_j) divides z_0.       (15)
```

Reduce `P_i'/bar(P_i')=A_i'/B_i'` to coprime Gaussian numerator and
denominator.  At odd split primes the reduction is conjugate-primitive;
if `1+i` is removed, keep the unit in the exact ratio `A_i'/B_i'`
because it may rotate the phase by a quarter-turn.  With the **same**
anchor, the swapped point is integral exactly when

```text
B_i' divides z_0.                                     (16)
```

If one instead scales the whole old configuration, the minimal
Gaussian multiplier needed for that anchor is, up to a unit,

```text
D_i=B_i'/gcd_G(B_i',z_0),
z_0'=D_i z_0,       R'=|D_i| R.                          (17)
```

The transferred core also bounds this multiplier.  Let
`Gamma'_i=product_(selected U) H_U`, of norm-logarithm at least
`16(1-eta)w`, and put
`R_i=Gamma'_i/gcd_G(Gamma'_i,P_i')`.  By (11), `R_i` divides the
ordinary integer `E_i`, and its one-orientation support gives
`Norm(R_i)|E_i`.  Every selected cut contains at least two of the
retained rows, so `bar(Gamma'_i)` divides the original `z_0`.
At odd primes the reduced denominator `B_i'` contains
`bar(Gamma'_i/R_i)`; the possible ramified phase reduction does not
touch those primes.  Thus

```text
Norm(gcd_G(B_i',z_0))>=Norm(Gamma'_i)/E_i,
|D_i|<=|P_i'| sqrt(E_i)/exp(8(1-eta)w)
     <=4sqrt(6) C_Q
       exp((5+75eta)w+3beta+4sigma)/sqrt(E_i).         (18)
```

This charges the new conductor through the actual old anchor.  The
`5w` allowance remains large enough for radius inflation; extra
content `E_i` can reduce it but is uncontrolled.

All retained points scale by `D_i`; their angular differences are
unchanged, but the new row may lie outside the original arc.  Thus
conductor or radius changes must be charged through `z_0` itself.
The least radius is not obtained by replacing `R` with an `lcm` of
the Cartesian row norms.

As an elementary algebraic obstruction, the checker also uses five
small primitive integer rows and an invariant with `C_Q=440`.
At one row its Vieta replacement changes squared norm `5` to `149`;
at another it collides projectively with a retained row.  That fixture
does **not** satisfy the full-cut small-residue endpoint hypotheses.
The literal full-cut fixtures verify (1)--(13) with unrestricted
invariant height; they do not turn (14)--(17) into a descent.
