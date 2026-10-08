# Two small exact invariant relations are impossible on six rows

For an actual full-profile central configuration on six nonanchor rows,
two independent exact invariant relations of degree two in each row
cannot both have coefficient height below `w-o(w)`. More precisely,
their coefficient sums have product at least `exp(2w-o(w))` once both
relations lie in the balanced kernel. The one-pass proof below first
gives `w/2` and `w`; the exact complement reflection in
[complement_reflection_invariant_relations.md](complement_reflection_invariant_relations.md)
doubles those thresholds. This uses the simultaneous
first-smaller-cut congruences; it is stronger than the balanced linear
congruences and the primitive evaluation heights alone.

It does not prove a uniform count. The available six-row
successive-minima estimate does not supply two relations this small.
The exact coefficient certificate is
[check_two_small_invariant_relations.py](check_two_small_invariant_relations.py).

## 1. Statement with finite losses

Use the actual integer hypotheses of
[near_balanced_invariant_congruences.md](near_balanced_invariant_congruences.md):
there are primitive Gaussian numerators `P_i=K_i A_i`, six nonanchor
rows, independent odd split-prime blocks `H_S`, and nonzero primitive
pair residues satisfying

```text
log|K_i| <= sigma,       0<|t_ij|<=T,       T>=1,
log N(H_S) >= (1-eta)w.
```

Let `V` be the rational simultaneous `SL_2` invariant space of degree
two in every row. It has dimension `15`. Let `K` be the kernel of its
balanced-cut evaluations; it has dimension `5`. For an integer
polynomial `Q`, let `C_Q` be the sum of the absolute values of its
coefficients in the original twelve binary coordinates.

**One-pass theorem.** If `Q,R` are independent integer invariants in `K` and
both vanish at the actual rows, then

```text
log C_Q + log C_R >= (1-eta)w - 12 sigma - 2 log T.        (1)
```

There is also a version without presupposing membership in `K`:
any two nonzero exact numerical relations in `V` satisfying

```text
log C_Q + log C_R + 12 sigma + 2 log T < (1-eta)w         (2)
```

are proportional over `Q`. In particular all exact relations with

```text
log C_Q < ((1-eta)w - 12 sigma - 2 log T)/2              (3)
```

span a vector space of dimension at most one. The span assertion
does not claim that every vector in the span obeys the same norm bound.

**Stronger two-pass theorem.** Complement reflection preserves every
numerical invariant zero and its coefficients, preserves the primitive
pair residues up to sign, and supplies a second coprime modulus. Hence
independent numerical zero relations in `K` satisfy

```text
log C_Q+log C_R >= 2(1-eta)w-156sigma-4log T.              (3a)
```

Without presupposing membership in `K`, all exact numerical relations
whose coefficient heights are strictly below

```text
(1-eta)w-78sigma-2log T                                   (3b)
```

span a space of dimension at most one. The reflection note gives the
full finite proof, including content removal at Gaussian primes and
zero imaginary coordinates. The one-pass proof remains useful as the
local ingredient in this stronger conclusion.

The proof retains every prime-power exponent. It needs neither a
fixed set of places nor an approximation of the original rows by a
function-field configuration.

## 2. Two equations remove the varying cross-ratio

Fix an inside pair `S`, and order the outside rows as `a<b<c<d`.
For `Q in K`, its first-smaller-cut restriction is

```text
F_(S,Q)(z)=alpha_Q (z_a-z_b)(z_c-z_d)
          +beta_Q (z_a-z_c)(z_b-z_d),
alpha_Q,beta_Q in Z,       |alpha_Q|+|beta_Q|<=C_Q.       (4)
```

Write `U=D_ab D_cd`, `V_0=D_ac D_bd`, and divide both by their
ordinary positive gcd to obtain `A,B` with `gcd(A,B)=1`.
Equations (9)--(12) of the preceding congruence note give an integer
`M_(1,S)`, depending only on the actual rows and `S`, for which

```text
M_(1,S) | alpha_Q A + beta_Q B,
M_(1,S) >= N(H_S) exp(-12 sigma) T^(-2).                 (5)
```

The same modulus and the same `A,B` work for every numerical zero
relation. Therefore, for two such relations in `K`, eliminating in
the two equations gives

```text
M_(1,S) | delta_S A,       M_(1,S) | delta_S B,
delta_S=alpha_Q beta_R-beta_Q alpha_R.
```

Bezout for the primitive pair `(A,B)` now gives the exact divisibility

```text
M_(1,S) | delta_S.                                      (6)
```

No individual matching product is presumed coprime to the modulus.
In particular, this step does not incorrectly pass coprimality from
summands to their sum. Also

```text
|delta_S| <= (|alpha_Q|+|beta_Q|)(|alpha_R|+|beta_R|)
          <= C_Q C_R.                                  (7)
```

Under the strict reverse of (1), equations (5)--(7) force
`delta_S=0` for every one of the fifteen inside pairs. Unlike the
bound for one relation, the cross-ratio height does not remain on
the right side of (7): the two simultaneous equations have eliminated
its numerical value.

## 3. The fifteen determinants determine the exterior product

For a sorted triple put `T_abc=Delta_ab Delta_bc Delta_ca`. Use
the following ordered basis of `K`:

```text
B_1=T_123 T_456,    B_2=T_124 T_356,    B_3=T_125 T_346,
B_4=T_134 T_256,    B_5=T_135 T_246.                     (8)
```

The dimension and integral basis conversion are proved in
[segre_gradient_arithmetic.md](segre_gradient_arithmetic.md).
Only its rational basis assertion is needed here.
Taking the determinant of the two restrictions defines a linear map

```text
L: exterior^2 K -> Q^15,
L(Q wedge R)=(delta_S)_(|S|=2).                          (9)
```

**This map is injective.** Here is a small exact rank certificate.
Order its ten columns by
`12,13,14,15,23,24,25,34,35,45` in the basis indices of (8).
Take its ten rows with exactly those labels, now referring to the
six actual row indices. The resulting matrix is

```text
 0  0  0  0  0  0  0  0  0  1
 0  0  0  0  1  0  0  0  0  0
 0 -1  0 -1  0  0  0  0  1  0
 1  0  1  0  0  1  0  0  0  0
 0  0  0  0  1  0  1 -1  0  1
 0 -1  0  0  0  0  0 -1  0  0
 1  0  0  0  0  0 -1  0  0  0
 0  0  0  1  0  0  1  0  0  0
 0  0 -1  0  0  0  0  1  0  0
 0  0  0  0  0  1  1  1  1  0
```

Its determinant is `-2`. The checker constructs all restrictions
directly from the binary determinant polynomials, verifies their
expansions in (4), and calculates this determinant with exact rational
arithmetic. Thus the assertion includes rank-degenerate choices of
`Q,R`: if all fifteen determinants vanish, `Q wedge R=0`, so they
are proportional. There is no generic-position qualification.

This proves (1). To deduce (2), the complementary balanced-divisor
theorem from
[invariant_relation_lattice_threshold.md](invariant_relation_lattice_threshold.md)
first places both relations in `K` whenever each coefficient height
is below

```text
2(1-eta)w-24 sigma.                                     (10)
```

Condition (2) implies this. Indeed, put
`D=(1-eta)w-12 sigma-2log T`. If (2) holds then `D>0`,
each coefficient height is less than `D`, and the bound in (10) is
`2D+4log T>=2D`. The just-proved wedge assertion applies.

## 4. What this rules out, and what it does not

For six rows the earlier CRT linear model has a primitive numerical
functional on `K` proportional to `(N,1,0,0,0)`. It has three
independent constant-coefficient numerical zero relations. Fixed
changes between integral bases only change their coefficient sums
by a fixed factor. Consequently, as `w` grows, this model violates
(2). The additional actual first-smaller-cut congruences exclude it;
matching the balanced congruences and the two primitive heights is
not enough to model the actual relation ideal.

Let `lambda_1<=...<=lambda_4` be successive minima of the actual
relation lattice in `K`, using any fixed coefficient norm. The
theorem gives

```text
log lambda_1 + log lambda_2 >= 2w-o(w).                  (11)
```

The primitive evaluation vector on `K` has height `16w+o(w)`,
so the lattice covolume has that logarithmic size. Minkowski alone
only bounds the second minimum by `exp(16w/3+o(w))`, because its
rank is four and its first nonzero integral vector has norm at
least a fixed positive constant. It does not guarantee two vectors
at the stronger threshold in (3b). A spectrum with all four logarithmic
minima equal to `4w+o(w)` is compatible with both the covolume and
(11). The full fifteen-dimensional evaluation lattice does not
remove this gap merely by projecting to `K`.

The four-outside-row feature is essential to this proof: the
restriction space in (4) has dimension two and clearing its
denominator leaves ordinary real matching products. For more than
six rows the first smaller restrictions have larger dimension,
and two numerical zero relations do not directly eliminate all
their varying values. No larger-row analogue of (1), no
classification of common factors, and no uniform arc bound is
asserted here.

Independent audit: the one-pass proof and exact determinant certificate
were checked by `symmetry_descent`. The common-reflection construction
and doubled finite constants were independently checked by the root
agent; both exact checkers pass.
