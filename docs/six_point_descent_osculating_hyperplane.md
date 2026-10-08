# An explicit osculating hyperplane in normalized square descent

The two-plane condition gives a single-point restriction stronger than the
ordinary small-embedding comparison in a specified low-degree case. At a
branch point of the projectivized descent curve, there is a unique
hyperplane with contact order twenty. Its value at an endpoint point of
projective beta-height `T` is `O(U^C T^-199)`. If the **actual field of
definition of this hyperplane** has degree `D<200`, this implies

```text
H^((200-D)/2) <= C(C_end) U^C.                         (1)
```

The exponent C includes the controlled frame and descent costs. It is
not evaluated here, so (1) does not establish the fair-scale cotangent
target `X<=U^(c+o(1))` for `c<8/7`. The degree-12288 case allowed by the
arithmetic classification remains outside the stated criterion.

## 1. Setup and the projective height conversion

Use the degree-twelve etale algebra and normalized representatives from
[the polynomial representative theorem](six_point_polynomial_squareclass_representatives.md):

```text
a0 s-theta t = eta beta^2,
gcd(s,t)=1,                 H=max(|s|,|t|),
N beta integral,           N<=U^C.
```

All coefficients, discriminants, the embeddings of eta and their
inverses have polynomial cost in U. Choose a rational integral basis of
the etale algebra with polynomially bounded conjugates and inverse
embedding matrix. Here a rational integral basis means a Q-basis of
integral elements; it need not generate the full maximal order. Such a
basis follows from the small-vector argument in the representative
theorem. Its index in the maximal order is polynomially bounded by its
embedding determinant and the field discriminant. Absorb this index
into N so that

```text
N beta = sum_j x_j e_j,       x_j in Z.
```

Write `x=g z`, where z is primitive integral, and `T=max|z_j|`.
Extracting s and t from the rational coefficients of
`eta (sum z_j e_j)^2` gives two rational quadratic forms with a common
denominator `Q<=U^C`. Thus

```text
Q N^2 s = g^2 S(z),       Q N^2 t = g^2 T0(z),
S,T0 in Z[z_1,...,z_12].
```

Consequently `g^2` divides `Q N^2`, by primitivity of `(s,t)`. The
embedding bounds and their inverse then give

```text
U^-C sqrt(H) <= T <= U^C sqrt(H).                       (2)
```

The lower bound uses at least one ordinary embedding of beta of size
`U^-C sqrt(H)` near a simple branch root. The upper bound and the
bounded content g follow as above. This is why passing to projective
beta coordinates does not silently discard an H-sized common factor.

Over the algebraic closure, diagonal coordinates `b_i=sigma_i(beta)`
put the projectivized two-plane locus in the form

```text
eta_i b_i^2 = a0(s-alpha_i t),      0<=i<=11.
```

Eliminating s,t gives ten independent diagonal quadrics in P^11. For
distinct alpha_i, their intersection is a smooth connected curve of
degree `2^10=1024` and genus 4097. It maps to P^1 with degree 2048 and
to the genus-five square-norm curve with degree 1024. These three
degrees must not be interchanged. The local calculation below only
needs this diagonal presentation and distinctness of the roots.

## 2. The complete local vanishing sequence

Fix the approached real root alpha_0 and one local sign sheet. Work in
the chart `r=s/t`, put `delta=r-alpha_0`, and choose a local square
root of t over C. Define

```text
v_i^2 = a0(alpha_0-alpha_i)/eta_i,
q_i = 1/(alpha_0-alpha_i),                         1<=i<=11.
```

The signs of the nonzero v_i specify the branch point B. On its local
sheet,

```text
b_0 = sqrt(t) sqrt(a0/eta_0) sqrt(delta),
b_i = sqrt(t) v_i (1+q_i delta)^(1/2),             1<=i<=11.  (3)
```

One can use `sqrt(|t|)` and absorb the sign into the constants when
bounding absolute values. All constants in (3), their inverses, and
the separations of the distinct nonzero q_i have polynomial cost in U.

Put `zeta=sqrt(delta)`. The first coordinate has order one in zeta;
the other eleven have expansions in even powers. Their coefficient
matrix through degree twenty is, up to nonzero column and row factors,
the eleven-by-eleven Vandermonde matrix `(q_i^k)_(0<=k<=10,1<=i<=11)`.
Every binomial coefficient `(1/2 choose k)` in this range is nonzero.
The vanishing sequence of ambient linear forms at B is therefore

```text
0, 1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20.             (4)
```

This proves both uniqueness of the order-twenty hyperplane and the
absence of a linear form of still higher contact. It also shows that
the relevant linear forms span all twelve ambient coordinates.

Here is an explicit form attaining order twenty. Set

```text
lambda_i = 1 / product_(j!=i,1<=j<=11)(q_i-q_j),
L_B(b) = sum_(i=1)^11 (lambda_i/v_i) b_i.                (5)
```

Lagrange interpolation gives `sum lambda_i q_i^k=0` for `0<=k<10`
and `sum lambda_i q_i^10=1`. Hence

```text
L_B(b) = sqrt(t) [(1/2 choose 10) delta^10
                         + O(U^C delta^11)].           (6)
```

The leading coefficient is a fixed nonzero rational number. Thus, for
`0<|delta|<U^-C` with an adjusted constant, `L_B(b)` is nonzero. This
nonvanishing is essential: a norm estimate alone would say nothing
about a rational point lying exactly on the auxiliary hyperplane.

## 3. Endpoint size and a one-point Liouville inequality

The endpoint branch estimate is

```text
0<|delta| <= C(C_end) U^C H^-10.
```

Equations (2) and (6), including the factor `N/g` relating b to z,
give, after increasing the polynomial threshold in U,

```text
0<|L_B(z)| <= C(C_end) U^C T^-199.                     (7)
```

Here and below L_B is expressed in the rational coordinate basis e_j.
Its coefficients have bounded algebraic degree and absolute
multiplicative height `U^C`: formula (5), root separation, and the
controlled inverse embedding matrix verify this directly.

Let E be the smallest field over which the hyperplane is defined, and
put `D=[E:Q]`. Divide by one nonzero rational-coordinate coefficient,
so that all coefficients belong to E. This changes their heights,
their relevant absolute values, and the estimate (7) by polynomial
factors only. Multiply by a positive integer `M<=U^C` to make all
coefficients integral. Such an M exists from their bounded degrees
and polynomial heights. At every embedding of E, the coefficients of
`M L_B` are at most `U^C`, so

```text
|(M L_B)^sigma(z)| <= U^C T.
```

The nonzero integral norm of `M L_B(z)` now gives

```text
1 <= |Norm_E/Q(M L_B(z))|
  <= C(C_end) U^C T^(D-200).                            (8)
```

For `D<200`, this proves (1), using (2). For `D>=200`, (8) supplies
no upper bound on H. In particular, small field discriminant alone
does not change this degree threshold.

The hyperplane is intrinsically the unique member of the rational
ambient linear series vanishing to order twenty at B. In fact its
field E is **equal** to the residue field of B, as proved next. An
arbitrary larger splitting field chosen to display (5) must not be
substituted for E when deciding applicability of (8).

## 4. The hyperplane field equals the branch field

The assignment `B |-> {L_B=0}` is Galois-equivariant, by the uniqueness
of the order-twenty hyperplane. It is also injective on the complete
branch locus. To see this, fix one diagonal coordinate system over
the algebraic closure. For a point B above alpha_j, the coefficient
of b_j in L_B is zero and every other coefficient is nonzero. Thus
the hyperplane determines j. For fixed j, the numbers lambda_i in
(5) are fixed. If the hyperplanes of B and B' are equal, their
coefficient vectors are proportional, giving

```text
lambda_i/v'_i = c lambda_i/v_i        (i!=j).
```

Hence `v'_i=c^-1 v_i` for every `i!=j`, and B'=B projectively. The
common invertible embedding matrix makes equality in these diagonal
coordinates equivalent to equality in the rational coordinate basis.
Galois-equivariance and injectivity give identical stabilizers, so

```text
E = Q(B),             D = [Q(B):Q].                     (9)
```

Under the **additional hypothesis** that the branch polynomial has
full Galois group S_12 over Q, the possible D are particularly
explicit. Fix one root alpha_j. Its field has degree twelve, and
over that field the linear Galois action on the 1024-point branch
fiber is S_11 on

```text
V = {even subsets of the remaining eleven labels},
dim_F2(V)=10.
```

The module V is simple. If a nonzero invariant subspace contains the
nonempty even subset S, choose a label in S and one outside S; adding
S to its transposed image gives that pair. All pairs are conjugate
and generate V. The translation kernel of the affine fiber action
is therefore either zero or V.

There is no nontrivial affine cocycle class for this S_11 action:

```text
H^1(S_11,V)=0.                                          (10)
```

For a short proof, the eleven-dimensional permutation module splits
as `F2^11 = 1 direct-sum V`, because eleven is odd. Shapiro's
identification gives
`H^1(S_11,F2^11)=H^1(S_10,F2)=F2`. The map from the constant summand
`H^1(S_11,F2)=F2` is restriction of the sign character to S_10, hence
an isomorphism. The remaining summand in H^1 is zero, proving (10).

If the translation kernel is zero, (10) makes the fiber action
translation-conjugate to the linear action. Its orbit sizes, from
even subset weights `0,2,4,6,8,10`, are

```text
1, 55, 330, 462, 165, 11.
```

If the kernel is V, its translations act transitively and the fiber
orbit has size 1024. The full S_12 action is transitive on the twelve
roots, so the resulting absolute branch degrees, and by (9) the
actual hyperplane degrees, are precisely the following alternatives:

```text
zero translation kernel:  12, 660, 3960, 5544, 1980, 132;
full translation kernel:  12288.                        (11)
```

The two degrees below 200 are 12 and 132. They give respectively

```text
H^94 <= C(C_end) U^C,       H^34 <= C(C_end) U^C.         (12)
```

This classification assumes full S_12 and does not identify which
translation kernel or branch orbit a particular twist realizes.
For twists carrying a rational point there is now an additional
restriction: if the discriminant squareclass is not -1, all branch
fibers have full translation kernel, so only D=12288 occurs. The
[four-torsion argument](six_point_pointed_cover_four_torsion.md)
proves this and verifies the discriminant condition for the first
existing S_12 witness. In particular the cases in (12) do not apply
to pointed twists of that witness.
The independent exact checker
[check_s12_even_subset_h1.py](check_s12_even_subset_h1.py) also checks
(10): its S_11 cocycle system has rank 90 on 100 variables, while
the coboundary space has dimension ten. It enumerates all six
displayed S_11 orbit sizes.

## 5. What this adds and what remains

The first two levels of (4) imply the basic tangent estimate: in
primitive beta coordinates there are eleven independent forms of
size `O(U^C T^-9)`. Ten forms annihilating the full tangent line
improve this to `O(U^C T^-19)`. The entire flag gives successive
bounds

```text
T, T^-9, T^-19, T^-39, ..., T^-199,
```

up to polynomial U factors. Thus the two-plane equations really do
contain stronger simultaneous linear approximation than the single
small conjugate. Their coefficients, however, need not live in the
original degree-twelve algebra. That field cost cannot be omitted.

When the branch fiber has full translation kernel, its residue
degree over Q is `12*1024=12288`. The order-twenty hyperplane
criterion does not cover this case. Taking a product of all conjugate
hyperplanes does not retain the exponent -199 without paying for the
other D-1 factors in (8).

A concrete next test is to determine the actual hyperplane field (or
a sufficiently small branch residue field) for the normalized twists
arising in the six-point family. Cases with degree below 200 obey the
proved inequality (1). Even there, its coefficient exponent must be
made explicit and compared with the fair moving-weight scale before
claiming any circle-arc consequence. No general uniform endpoint
bound or improvement in the radius-dependent growth estimate is
claimed here.
