# Real congruences at the first smaller cuts

This note sharpens the local restriction in
[small_invariant_relations_at_balanced_cuts.md](small_invariant_relations_at_balanced_cuts.md).
For six nonanchor rows, a single matching specialization of a sufficiently
small exact invariant relation is impossible. A general specialization is
instead a cross-ratio congruence. Its available height is larger than the
new modulus, so this does not prove a residue-growth or uniform-count theorem.

A later simultaneous argument does eliminate the varying cross-ratio
between two relations. See
[two_small_invariant_relations.md](two_small_invariant_relations.md)
and its exact complement-reflection strengthening: all six-row exact
relations of coefficient height below `w-o(w)` span at most one
dimension. The one-relation height comparison below remains valid.

The accompanying exact certificate is
[check_near_balanced_invariant_congruences.py](check_near_balanced_invariant_congruences.py).
It checks a fixed invariant whose smaller-cut restrictions all avoid the
single-matching alternatives. It does **not** construct a numerical zero
relation at an actual lattice configuration.

For two simultaneous numerical zero relations, the varying cross-ratio
can now be eliminated. The exact
[two-small-relations theorem](two_small_invariant_relations.md) proves
that independent relations in the six-row balanced kernel have
`log C_Q+log C_R >= (1-eta)w-12s-2log T`. A full-rank exterior-product
certificate makes this stronger than the single-relation restriction
proved below, while the available lattice bounds still do not force
two relations under that threshold.

## 1. Primitive integer data and two local facts

Use the actual centrally truncated system, with its fixed common unit class:

```text
P_i=X_i+iY_i=K_i A_i,     A_i=product_(T containing i) H_T,
gcd_G(P_i,bar(P_i))=1,    0<|Y_i|,
log|K_i|<=s.
```

All core blocks and all their conjugates have disjoint prime support. They
use odd split primes. For each pair put

```text
g_jl=gcd_G(P_j,P_l),
D_jl=Im(P_j bar(P_l))=Y_j X_l-Y_l X_j,
t_jl=Im(P_j bar(P_l)/N(g_jl)),     0<|t_jl|<=T,  T>=1.
```

The quotient is a conjugate-primitive Gaussian integer, and

```text
D_jl=N(g_jl)t_jl.                                          (1)
```

Fix a block `H=H_S`, and write `N=N(H)`. For an inside row `i`,

```text
gcd_G(H,Y_i)=1.                                            (2)
```

Indeed `H|P_i`, whereas `H` is coprime to `bar(P_i)`; use
`2iY_i=P_i-bar(P_i)` and the oddness of `N`.

For two outside rows `j,l`, the valuations of both `P_j,P_l` at
every prime over a rational prime dividing `N` come entirely from
`K_j,K_l`. Consequently, as an ordinary positive integer,

```text
gcd(N,D_jl) divides N(gcd_G(K_j,K_l)) |t_jl|.                (3)
```

This retains arbitrary prime powers. At a rational prime `p|N`,
the contribution of `N(g_jl)` is the sum of the common valuations
in the two Gaussian orientations; each agrees with the corresponding
valuation of `gcd_G(K_j,K_l)`. In particular

```text
gcd(N,D_jl)<=exp(2s)T,
gcd(N,D_ab D_cd)<=exp(4s)T^2.                              (4)
```

The second inequality follows prime by prime from
`gcd(N,uv)|gcd(N,u)gcd(N,v)`.

We also use an exact real divisibility rule. For a Gaussian integer
`kappa` and an ordinary integer `b`,

```text
H | kappa b
  if and only if
N/N(gcd_G(H,kappa)) | b.                                  (5)
```

At an oriented Gaussian prime `pi^e||H`, an ordinary integer has the
same valuation at `pi` and `bar(pi)`. The remaining required exponent
is therefore `max(0,e-v_pi(kappa))`, exactly the exponent of `p`
in the displayed real modulus. No inversion modulo a prime power is
being assumed.

For example, (2) and (5) strengthen the balanced scalar conclusion of
the earlier note: its ordinary integer `c_S` satisfies

```text
N/N(gcd_G(H, product_(j outside S) K_j^2)) | c_S.
```

Thus a nonzero balanced scalar has

```text
log|c_S| >= log N - 4 sum_(j outside S)log|K_j|.             (6)
```

The imaginary-residue losses in the earlier sufficient bound are
unnecessary under the explicit primitive-row hypothesis used here.

## 2. Six rows give a real matching congruence

Let `m=6`. Let `Q` be an integer simultaneous `SL_2` invariant of
degree two in each row, with coefficient sum `C_Q`, which vanishes
at the actual six rows. Assume all its balanced scalars vanish.
For `|S|=2`, label the four outside rows `a<b<c<d`. The smaller
restriction is uniquely

```text
F_S(z)=alpha (z_a-z_b)(z_c-z_d)
       +beta (z_a-z_c)(z_b-z_d),
alpha,beta in Z,       |alpha|+|beta|<=C_Q.                 (7)
```

This follows from the zero diagonal and zero row sums of the quadratic.
Its coefficient of `z_a z_c` is `alpha`, and its coefficient of
`z_a z_b` is `beta`. At the actual points, `z_j=Y_j/P_j`, so

```text
F_S(z)=(alpha U+beta V)/(P_a P_b P_c P_d),
U=D_ab D_cd,       V=D_ac D_bd.                            (8)
```

The specialized invariant is
`product_inside Y_i^2 product_outside P_j^2 F_S`.
Its divisibility by `H`, together with (2) and coprimality of the
outside core factors, gives

```text
H | kappa (alpha U+beta V),
kappa=K_a K_b K_c K_d.
```

The quantity in parentheses is an ordinary integer. By (5),

```text
M | alpha U+beta V,
M=N/N(gcd_G(H,kappa)) >= N exp(-8s).                       (9)
```

This real normalization is specific to four outside rows. With more
outside rows, clearing a four-cycle denominator leaves products of
the other, large Gaussian coordinates. Equation (9) is not asserted
for those cases.

## 3. Single matchings are excluded, not all quadratics

The three single-matching alternatives in (7) are

```text
alpha=0,       beta=0,       alpha+beta=0.
```

For the last case use the exact identity
`D_ab D_cd-D_ac D_bd+D_ad D_bc=0`. If `F_S` is a nonzero
single matching with coefficient `c`, (9) and (4) imply

```text
|c| >= M/gcd(M, matching product)
    >= N exp(-12s) T^(-2).                                (10)
```

Therefore the sufficient threshold

```text
log C_Q+12s+2log T < log N                                (11)
```

forces every `F_S` either to be zero or to satisfy
`alpha beta(alpha+beta)!=0`. In the latter case it is irreducible
even over `C`: after setting `z_d=0`, its symmetric coefficient
matrix has determinant `-alpha beta(alpha+beta)/4`, so its rank is
three. A product of two linear forms has rank at most two.

This is a further genuine restriction on a small numerical zero
relation, but irreducible restrictions are still possible.

## 4. The remaining cross-ratio modulus is smaller than its height

Let `g=gcd(|U|,|V|)`, and write `U/g=A`, `V/g=B`. Then
`gcd(A,B)=1`; change both signs if a positive denominator is wanted.
Equation (9) becomes the exact cleared congruence

```text
M_1 | alpha A+beta B,
M_1=M/gcd(M,g) >= N exp(-12s)T^(-2).                       (12)
```

This does not presume that `A` or `B` is a unit modulo `M_1`.
For example the congruence itself implies
`gcd(M_1,A)|beta` and `gcd(M_1,B)|alpha`.
It remains essential that a sum of matching products can have a
large divisor even though neither matching does.

The rational number `A/B` is an ordinary projective cross-ratio of
the same circle points. For seven selected rows (one anchor and six
nonanchor rows), there are `64` full-profile blocks. If their common
log norm is `w`, the norm in
[coupled_crossratio_height_norm.md](coupled_crossratio_height_norm.md)
has value `ell=8` on this quadruple cross-ratio. Thus

```text
h(A/B)=8w+o(w),       log M_1>=w-o(w).                    (13)
```

If `log C_Q=o(w)`, then `alpha A+beta B` cannot be zero for a
nonzero restriction: an exact zero would give
`h(A/B)<=log max(|alpha|,|beta|)`. Its nonzero value yields only

```text
log(|alpha|+|beta|) >= log M_1-h(A/B)
                   >= -7w-o(w).                          (14)
```

There is no positive coefficient lower bound in (14). A simultaneous
argument among the different restrictions would need additional
information; the multiplicative cross-ratio height identity alone
does not supply it.

## 5. A fixed invariant avoiding every single-matching restriction

Number the rows `1,...,6`. List the ten triples containing row `1`
in lexicographic order:

```text
123,124,125,126,134,135,136,145,146,156.
```

For a sorted triple `abc`, put
`T_abc=Delta_ab Delta_bc Delta_ca`, where
`Delta_ij=x_i y_j-x_j y_i`. If `J_j` is the `j`th listed triple,
with indices `j=0,...,9`, define the fixed integer invariant

```text
Q_* = sum_(j=0..9) 2^j T_(J_j) T_(J_j complement).         (15)
```

Each summand has degree two in each row. Every balanced restriction
vanishes: at least one of its triangles has two inside vertices,
making one determinant zero. Its coefficient sum is at most
`64(1+2+...+512)=65472`.

The exact certificate substitutes the binary vectors directly, without
using numerical approximations. All fifteen smaller restrictions have
the following coefficients in (7):

| Inside pair | alpha | beta | alpha+beta |
|---|---:|---:|---:|
| 12 | 336 | 480 | 816 |
| 13 | -390 | -372 | -762 |
| 14 | 553 | -108 | 445 |
| 15 | 263 | 54 | 317 |
| 16 | 109 | 54 | 163 |
| 23 | 54 | -108 | -54 |
| 24 | 263 | -372 | -109 |
| 25 | 553 | -390 | 163 |
| 26 | -445 | 762 | 317 |
| 34 | -317 | 480 | 163 |
| 35 | -445 | 336 | -109 |
| 36 | 553 | -816 | -263 |
| 45 | 762 | -816 | -54 |
| 46 | -390 | 336 | -54 |
| 56 | -372 | 480 | 108 |

Thus `Q_*` is nonzero and every smaller restriction has rank three.
This is an example in the invariant polynomial space, not an example
of `Q_*(P_1,...,P_6)=0` for the arithmetic data. It shows precisely
why ruling out every single matching does not rule out the entire
balanced kernel.

## 6. Exact outside arithmetic can create a new divisor in the sum

The outside-bracket estimates do not, by themselves, exclude new large
divisors in a sum of two matchings. This can be seen with four actual,
equal-degree, constant-imaginary numerators. The following example uses
the rational five-point family in
[four_row_rational_quartic_construction.md](four_row_rational_quartic_construction.md).
Its additional exact certificate is
[check_outside_matching_cancellation.py](check_outside_matching_cancellation.py).

Use a real polynomial variable `X`. In that construction, write the four
unit-circle roots as `z_1,...,z_4`, the private roots as `w_1,...,w_4`, and

```text
h_j=(X-w_j) product_(l!=j)(X-z_l),
tau_j=Im h_j != 0,
q_j=N(X-z_j)=X^2-2 Re(z_j) X+1.
```

The `tau_j` are distinct nonzero rational constants. The first two roots
have real parts `-3/5` and `-2117/3085`. Consequently the exact real
polynomial identity

```text
2117 q_1-1851 q_2=266(X^2+1)                              (16)
```

holds. Neither `i` nor `-i` occurs among any of the eight roots or their
conjugates.

Take original point `4` as the new anchor, and retain original outside
labels `0,1,2,3`. Their canonical polynomial numerators are

```text
P'_0=bar(h_4),
P'_j=h_j bar(h_4)/N(g_j4),     j=1,2,3,
g_j4=product_(l notin {j,4})(X-z_l).
```

Every `P'_j` is a monic quartic with nonzero constant imaginary part:
these constants are `-tau_4,tau_1-tau_4,tau_2-tau_4,tau_3-tau_4`.
Put

```text
D'_jl=Im(P'_j bar(P'_l)),
U'=D'_01 D'_23,       V'=D'_02 D'_13,
a_1=-tau_1(tau_2-tau_3),
a_2=-tau_2(tau_1-tau_3),
W_4=N(X-w_4)^2.
```

Reanchoring multiplies every matching product on these four labels by
the same real polynomial `W_4`. More explicitly,

```text
U'=a_1 W_4 q_1 q_4,
V'=a_2 W_4 q_2 q_4.
```

Equation (16) therefore gives the exact identity

```text
(2117/a_1)U'-(1851/a_2)V'
    =266 W_4 q_4 (X^2+1).                                (17)
```

The coefficients on the left are fixed nonzero rationals. All six
individual brackets `D'_jl` are coprime as polynomials to `X^2+1`.
Thus the factor on the right is a new factor of the sum, not a shared
factor of the rows or individual brackets. The certificate checks both
the polynomial identity and these six nonvanishings at `i` and `-i`
using exact rational arithmetic.

This identity also has its literal integer interpretation. Specialize
at even integers `n`. The integer block normalization in the construction
has odd norms. Clearing the fixed rational coefficient denominators of
each `P'_j(n)` gives an odd-norm Gaussian integer. Dividing its rational
coordinate gcd gives the actual conjugate-primitive numerator without
a ramified-prime rotation. These coordinate gcds are bounded: each divides
the fixed nonzero imaginary coordinate of the corresponding cleared row.
There are therefore only finitely many possible positive real scaling
factors between the polynomial and primitive integer numerators.

All matching products acquire the same product of these four row scaling
factors. A single further fixed integer multiplier clears the rational
coefficients and every possible denominator on the right of (17).
Consequently there are fixed integers `alpha,beta` such that

```text
n^2+1 divides alpha U'_primitive(n)+beta V'_primitive(n)
```

for every such specialization. The primitive pair residues stay bounded
by the same fixed-factor specialization argument used in the original
construction. The new conjugate-primitive Gaussian factor `n-i` has
unbounded modulus, while its gcd with every outside row and individual
bracket is bounded. The latter assertion follows from the nonzero fixed
resultants with `X-i`; clearing denominators only changes the fixed bound.

The four outside rows have size comparable to `n^4`. They are part of
the actual five-point circle construction of radius comparable to `n^8`
and an arc of length `O(n^4)`. This example therefore retains actual
integer residues and endpoint-scale angular geometry. It does **not**
realize `n-i` as a nonempty conductor block of an enlarged seven-point
central system, does not have its full cut profile, and does not produce
a numerical zero of the invariant in Section 5. Those extra simultaneous
conditions could still obstruct the congruences (12). The example isolates
why an outside-only gcd estimate for a matching sum cannot provide that
obstruction.

## Independent audit

The uniformity-audit agent independently verified (1)--(14), including
the exact identity

```text
N(gcd_G(H,D_jl))=gcd(N(H),D_jl),
```

the correction and residue costs in (3), the real normalization (5),
and the constants `8s` and `12s` in (9)--(12). No additional ramification
factor is needed under the stated odd-block and primitive-row hypotheses.
The same agent supplied and exactly checked the construction in Section 6.
