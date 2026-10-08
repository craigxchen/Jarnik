# Small exact invariant relations vanish on every balanced cut

This is an exact arithmetic consequence of the centrally truncated integer
system. It concerns small-coefficient **zero relations**, which are outside
the nonzero-value auxiliary-polynomial ceiling. It does not prove that these
relations vanish identically or that their values cannot be zero at the actual
configuration.

## 1. Integer data and the invariant relation

Use an even number `m` of nonanchor rows, with

```text
P_i=X_i+iY_i=K_i A_i,
A_i=product_(S containing i) H_S,
X_i,Y_i in Z,       Y_i!=0.
```

The Gaussian blocks are pairwise coprime and conjugate-coprime, as in
[endpoint_central_truncation.md](endpoint_central_truncation.md). In particular
`gcd(H_S,A_j)=1` when `j` is outside `S`. Corrections `K_j` may share primes
with `H_S`; their contribution is retained below.

Let `Q(v_1,...,v_m)` be a nonzero integer polynomial in binary vectors
`v_i=(x_i,y_i)`, homogeneous of degree two in each vector and invariant under
simultaneous `SL_2` transformations. Let `C_Q` be the sum of the absolute
values of its coefficients. Suppose

```text
Q((X_1,Y_1),...,(X_m,Y_m))=0.                              (1)
```

For each subset `S` of size `m/2`, define the Gaussian integer

```text
c_S=Q(v_i=u for i in S, v_j=e for j outside S),
u=(-i,1),       e=(1,0).
```

Since every coordinate of `u,e` has modulus at most one,

```text
|c_S|<=C_Q.                                                (2)
```

In fact `c_S` is an ordinary integer. An invariant of total degree
`2m` transforms under `GL_2` by `det^m`, by rescaling a matrix to
determinant one over `C`. The matrix with columns `u,e` has determinant
`-1`, so `c_S=(-1)^m Q(e_1 on S,e_2 outside S)`, an integer.

## 2. Exact balanced-cut factorization

For arbitrary outside vectors `v_j=(x_j,y_j)`, invariant weight balance gives

```text
Q(v_i=Y_i u for i in S, v_j for j outside S)
 =c_S product_(i in S)Y_i^2
       product_(j outside S)(x_j+i y_j)^2.                 (3)
```

Here is an explicit justification. The vectors `u,e` form a basis with
determinant `-1`, and

```text
(x_j,y_j)=y_j u+(x_j+i y_j)e.
```

Conjugating the diagonal map with eigenvalues `lambda,lambda^(-1)` in
this basis gives an `SL_2` transformation. The fixed `m/2` inside vectors
have total positive weight `m`. The outside vectors have total degree
`m`, so invariance permits only their component of weight `-m`: every
outside factor must use the `e` coordinate twice. Its coefficient is
exactly `c_S`. No division by two, coefficient denominator, or choice of
Gaussian unit occurs in (3).

At the actual integer configuration, `H_S | P_i` for `i in S`, hence
`X_i=-iY_i modulo H_S`. Equations (1) and (3) imply

```text
H_S divides c_S product_(i in S)Y_i^2
                    product_(j outside S)P_j^2.
```

Removing only the outside factors `A_j`, which are coprime to `H_S`, gives
the exact divisibility

```text
H_S divides c_S product_(i in S)Y_i^2
                    product_(j outside S)K_j^2.           (4)
```

Thus, if

```text
log |H_S| > log C_Q
               +2 sum_(i in S)log |Y_i|
               +2 sum_(j outside S)log |K_j|,             (5)
```

then necessarily `c_S=0`. This implication is valid for arbitrary prime
powers and does not presume the correcting factors or residues are coprime
to the conductor. It uses the modulus of an actual nonzero Gaussian
integer divisor.

For example, if `log N(H_S)>=(1-eta)w`, `log|Y_i|<=tau w`, and
`log|K_i|<=kappa w`, a sufficient simultaneous threshold is

```text
log C_Q/w +m(tau+kappa) < (1-eta)/2.                       (6)
```

The fixed-`m`, increasingly accurate central extraction can satisfy (6)
for any relation with `log C_Q=o(w)`. For growing `m`, the factors `m`
in (6) must be retained.

If the `P_i` are the actual primitive numerators `h_i` supplied by central
truncation, there is a sharper version. For `i in S`,
`gcd(H_S,Y_i)=1`: a Gaussian prime dividing both would divide `P_i`
and `bar(P_i)=P_i-2iY_i`, contrary to primitivity; the conductor primes
are odd. Thus (4) sharpens to

```text
H_S divides c_S product_(j outside S)K_j^2.                (6a)
```

The sufficient threshold becomes

```text
log C_Q/w +m kappa < (1-eta)/2.                           (6b)
```

The more general bounds (4)--(6) remain valid if nonprimitive correcting
numerators are used instead.

The fact that `c_S` is a real integer gives a stronger exact version.
Define

```text
J_S=gcd_G(H_S, product_(j outside S)K_j^2).
```

Equation (6a) implies `H_S/J_S | c_S`. This quotient is coprime to
its conjugate, so conjugating and multiplying the two coprime
divisors yields

```text
N(H_S/J_S) divides c_S in Z.                              (6c)
```

In particular, a sufficient vanishing threshold is

```text
log C_Q < log N(H_S)-4 sum_(j outside S)log|K_j|,
```

or, using the uniform bounds,

```text
log C_Q/w+2m kappa < 1-eta.                               (6d)
```

There is no coprimality assumption between `H_S/J_S` and the
correcting product in this argument. Primewise subtraction of its
valuation from (6a) gives the asserted divisibility directly.

For even `m`, complementary balanced evaluations agree:
`c_S=c_(S^c)`. The two norm divisors in (6c) are coprime, so their
product divides this same integer. The resulting stronger threshold is

```text
log C_Q < 2(1-eta)w-4 sum_i log|K_i|.
```

Every exact relation below this threshold vanishes on all balanced
cuts. The exact balanced rank, both primitive evaluation heights, and
the remaining successive-minima limitation are proved in
[the full relation-lattice audit](invariant_relation_lattice_threshold.md).
In particular, the existence of many small relations still does not
force a relation outside the balanced kernel.

## 3. What the resulting vanishing means

Under (5) for all balanced cuts, `Q` vanishes on every configuration
whose binary vectors lie on at most two projective directions.

For two independent directions, an invertible linear change of coordinates
reduces to the pair `u,e` above. An `SL_2`-invariant polynomial of total
degree `2m` transforms under `GL_2` by the scalar `det^m`; this follows by
rescaling a matrix to determinant one over `C`. Row scalings contribute
their squares. If exactly `m/2` rows use each direction, the remaining
coefficient is `c_S=0`. If the two counts are unequal, diagonal weight
balance forces zero without using (5). The case of only one direction
follows either by the same weight argument or by polynomial continuity.

This is a genuine restriction on a small numerical relation, beyond the
statement that its monomials have tied minimum valuations. It is still a
restriction on an ideal of polynomials, not a factorization of each
polynomial into individual pair differences.

## 4. Why this restriction alone does not finish the argument

Write `Delta_ij=x_i y_j-x_j y_i`. For `m=6`, the nonzero invariant

```text
Q_0=Delta_12 Delta_23 Delta_31
              Delta_45 Delta_56 Delta_64                  (7)
```

has degree two in every row and vanishes on every configuration with at
most two projective directions. Each triple necessarily contains two
rows on the same direction. For any larger even `m`, multiply (7) by
`Delta_78^2 Delta_9,10^2 ...` to obtain another example in the same
kernel. At a configuration of distinct projective directions these
particular products are nonzero. However, a linear combination of such
polynomials can vanish at a particular distinct configuration. Nonzero
values of generators cannot be used to rule out a zero value of every
polynomial in the generated ideal.

There is also an elementary dimension warning for a Siegel construction.
Suppose an integral invariant basis has `r` elements, each coefficient
sum at most `B`, and each numerical value has modulus at most `Z`.
Impose both the numerical zero relation and all balanced conditions
`c_S=0`. There are at most `s=binomial(m,m/2)` balanced conditions,
and their entries have modulus at most `B`. If `r>s+1`, box
pigeonholing still supplies a nonzero integer relation: comparing
coefficient vectors in `{0,...,H}^r`, the number of images is at most

```text
(2rZH+1)(2rBH+1)^s.
```

Whenever this is less than `(H+1)^r`, two coefficient vectors have
the same image. Their difference gives a numerical zero relation
already lying in the balanced-cut kernel. In particular the leading
large-height term in `log H` is only

```text
log Z/(r-s-1),                                             (8)
```

with the displayed bounded-entry costs also retained. Thus an invariant
space whose dimension is much larger than the number of balanced cuts
still permits very small exact relations after this new restriction.
The argument is conditional only on the supplied basis dimension and
bounds; it does not assume a function-field parametrization of the
integer data.

## 5. The first smaller cuts leave four-cycle quadratics

The smaller-cut specializations can be described explicitly. Put
`a=|S|<=m/2` and again write each outside vector as

```text
v_j=Y_j u+P_j e,       z_j=Y_j/P_j.
```

There is a polynomial `F_S` with ordinary integer coefficients such that

```text
Q(v_i=Y_i u inside S, v_j outside S)
 =product_(i in S)Y_i^2 product_(j outside S)P_j^2
                         F_S((z_j)_(j outside S)).        (9)
```

This is a polynomial identity after the denominators on the right are
multiplied out. Each variable of `F_S` has degree at most two. More explicitly,

```text
F_S(z)=Q(e_1 inside S, (z_j,1) outside S),
coefficient sum of F_S <= C_Q.                            (9a)
```

Indeed the change of basis with columns `u,e` has determinant `-1`,
and its `GL_2` multiplier is `(-1)^m=1`. Consequently there is no
coefficient growth from the complex translations in this invariant
subspace. The same diagonal weight argument as before gives

```text
F_S is homogeneous of degree m-2a.                        (10)
```

Furthermore it is invariant under simultaneous translation of the `z_j`.
Indeed the unipotent `SL_2` map fixing `u` and sending `e` to `e+c u`
leaves the left side invariant and sends each `z_j` to `z_j+c`.
The bound in (9a) is independent of the numerical outside coordinates.

Now take `a=m/2-1`, and assume the balanced scalars of Section 2 vanish.
Then `F_S` has degree two, and the coefficient of `z_j^2` is exactly
`c_(S union {j})`: set that outside vector to `u` and all others to
`e`. Consequently

```text
F_S(z)=sum_(j<l outside S) b_jl z_j z_l,
sum_(l!=j) b_jl=0       for every outside row j.           (11)
```

The second assertion follows by comparing the coefficients linear in
the translation parameter. Thus the first smaller-cut restriction
lands in exactly the zero-row-degree edge space that is integrally
generated by four-cycles in
[coupled_crossratio_height_norm.md](coupled_crossratio_height_norm.md).
For a four-cycle its quadratic is a product of two differences of
the `z` coordinates. No new diagonal coefficient remains to isolate.

This identifies the failure of an immediate scalar induction. The
actual relation (1) makes the polynomial on the right of (9) divisible
by `H_S`, but (11) can still depend nontrivially on the outside
directions. Its coefficient height is small; its evaluated numerator
need not be small. Dividing by the outside `P_j` without retaining
their denominators would incorrectly turn this into a bound on the
coefficients `b_jl` themselves.

The example (7) makes the obstruction explicit. Take `m=6` and
`S={1,4}`. After fixing those two vectors to `u`, its remaining
quadratic, up to a fixed sign, is

```text
F_S(z)=(z_2-z_3)(z_5-z_6).                                (12)
```

It is nonzero and has exactly the form (11), while all balanced
scalars vanish. Its differences at the actual configuration are

```text
z_j-z_l=(Y_j X_l-Y_l X_j)/(P_j P_l).
```

They contain the actual primitive pair residues together with their
gcd factors and the displayed denominators. For primitive `P_j`, this
gives more than a formal description. Put

```text
Delta_jl=X_j Y_l-X_l Y_j
        =N(gcd_G(P_j,P_l)) t_jl.
```

The residue orientation in this display is chosen consistently; its
sign has no effect on any of the following divisibilities or bounds.

For `j,l` outside `S`, both `A_j` and `A_l`, including their conjugates,
are coprime to `H_S`. Primewise comparison therefore gives

```text
gcd_G(H_S,Delta_jl)
 divides N(gcd_G(K_j,K_l)) t_jl.                           (13)
```

Thus an individual matching product such as (12), after the displayed
outside denominators have been cleared, can acquire an `H_S` divisor
only through its small coefficient, primitive residues, and correcting
factors. In the negligible-height regime it cannot supply the full
nonunit block. Consequently, if an exact small relation has such a
single nonzero matching product as one of its first smaller-cut
restrictions, that relation is excluded by (13). The example (7) is
therefore an obstruction to scalar-only classification, rather than
an actual possible zero relation at distinct integer directions.

More quantitatively, if
`F_S=a(z_j-z_l)(z_r-z_s)` with four distinct outside indices and
nonzero integer `a`, then a necessary condition for (1) is

```text
log|H_S| <= log|a|+2 sum_(v outside S)log|K_v|
                        +log|t_jl|+log|t_rs|.             (14)
```

To check the correction constant, clearing the four displayed
denominators leaves the factor `K_v` to exponent one at the four
matched rows and exponent two at every other outside row. Each
bracket contributes `N(gcd(K_j,K_l))`, whose modulus is at most
`|K_j K_l|`. Hence the total correcting modulus is at most
`product_(v outside S)|K_v|^2`, as used in (14). The inside `Y_i`
are units modulo `H_S` by primitivity.

The remaining problem is a **sum** of matching products in (11).
After one matching product is divided out, (13) removes its expected
local contents. With exactly four outside rows (`m=6`), the remaining
sum is a small-coefficient linear combination of projective
cross-ratios. Its numerator can have an additional divisor `H_S` at
a prime where the individual cross-ratio units have no such divisor.
With more than four outside rows, different matchings may use
different row sets; their ratios also retain the corresponding
moving `P_j` factors and need not be projective cross-ratios. The
general expression must therefore be kept as an evaluated matching
combination with its denominators. The simultaneous multiplicative
height norm does not bound this additional cancellation. Ruling it
out for all these coupled smaller-cut combinations would be a
further arithmetic result; it has not been proved here.

The companion [near_balanced_invariant_congruences.md](near_balanced_invariant_congruences.md)
gives the sharper real modulus for six rows and an explicit invariant
whose fifteen smaller restrictions all avoid the single-matching cases.

Thus the scalar restriction at balanced cuts is proved, and its
first recursive stage is now explicit. A further arithmetic estimate
on the evaluated quadratics (11), or a higher-order restriction
beyond the determinantal vanishing already exhibited by triangles,
is still missing.

## Audit

The root agent independently checked Sections 1--5, including the
exact correction factors in (4), the threshold (6), invariant weight
balance, the two-direction conclusion, and the small-relation
pigeonhole limitation, together with the first smaller-cut formulas
and the finite single-matching exclusion (14). The fresh-algebraic agent independently checked
the same normalizations and finite bounds, and the first smaller-cut
formula (9)--(11). It also supplied the primitive improvement (6a)
and the sharpened bracket-content comparison (13). Exact Gaussian
arithmetic checked 1,200 balanced and unbalanced invariant
specializations on four through ten rows. No uniform circle-point
bound is claimed.
