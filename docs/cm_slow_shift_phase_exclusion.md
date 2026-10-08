# CM phase exclusion with slowly varying shifts

Fix `m>=3`, put `B=2^m`, and keep
`E:y^2=x^3-2x`, `P=(2,2)`, and `q=hhat(P)>0`. The fixed-shift
residue theorem in
[the collective phase audit](cm_denominator_phase_lift_audit.md)
extends to distinct integer shifts that vary with `N`, provided that

```text
K_N=max_b |k_b(N)| <= N^(1/2-epsilon)                  (1)
```

for some fixed `epsilon>0`. The number of rows is fixed throughout.

Use the Gaussian denominator blocks of
`[N-k_b(N)+i]P`, with the bounded-common-factor trimming of the CM
construction, and form conjugate-primitive rows

```text
P_(N,i)=K_(N,i) product_(b containing i) H_(N,b),
s_N=max_i log max(2,|K_(N,i)|)=o(N^2).
```

Assume the resulting directions are distinct. Then their primitive
residues satisfy

```text
liminf_(N->infinity) log max_(i<j)|t_(N,ij)|/N^2
 >= (2^(m-1)-2)q.                                    (2)
```

Thus this larger family still cannot supply subquadratic logarithmic
residues. The corresponding actual-arc corollaries remain valid, with
the same leading exponents as in Section 5 of the phase audit. This is
not uniform in the number of rows, the curve, the generator, or
unrestricted shift growth.

## 1. Uniform data bounds

Write `K=K_N` and `H=1+K^2`. All constants below can depend on `m` and
the fixed curve and point, but not on the particular shifts satisfying
(1). The point `[k]P` has projective coordinate height `O(1+k^2)`.
The chord functions defining the phase words consequently have bounded
degrees and coefficient heights `O_m(H)`. Products defining the
`m-1` collective phase functions preserve these bounds, since their
number is fixed. Thus each `F_(N,i)` is a rational function on the
fixed Weierstrass curve, of degree `O_m(1)` and coefficient height
`O_m(H)`.

The exact phase constants have height `O_m(H)` as well. The Gaussian
denominator generator of `[k+i]P` has logarithmic modulus `O(1+k^2)`.
The common factors removed between two blocks divide denominators of
`[k_b-k_c]P` or `[k_b-k_c+2i]P`, whose logarithmic norms are `O(H)`.
There are only `O_m(1)` pairs. Ramified factors remain bounded. Hence
the entire trimming and phase correction costs `O_m(H)`, even though
it is no longer bounded independently of `N`.

In particular, uniformly in the shifts,

```text
log|H_(N,b)|=qN^2+O_m(NK+H),
log|h_(N,ij)|=(B/2)qN^2+O_m(NK+H+s_N).               (3)
```

The phase targets obtained from the row corrections and the exact
finite/unit/trimming factors have height `O_m(s_N+H)`.

## 2. A bounded-degree elimination lemma

The following elementary effective-algebra observation records the
specialization issue explicitly. Work over the fixed field `Q(i)` and
the fixed Weierstrass curve. Suppose a fixed number of rational
functions have degrees bounded by a fixed number and coefficient
heights at most `H`. Then:

* The ideal of their image in the algebraic torus admits generators of
  bounded degree and coefficient height `O(H+1)`.
* If a further bounded-degree function `u`, of coefficient height
  `O(H+1)`, belongs to their generated function field, one can write
  `u=A(F)/D(F)` with `A,D` of bounded degree and coefficient height
  `O(H+1)`, and `D(F)` nonzero as a function.
* The zeros of `D(F)` on the curve have bounded degree over the base
  field and height `O(H+1)`.

Here every degree bound and implied constant depends only on the
fixed input count and degree bounds. No chosen denominator or pivot
minor is asserted to remain nonzero under all specializations.

To justify the first assertion, describe the graph by clearing the
function denominators and adjoining one variable with equation
`z product(denominators)-1=0`; include inverses of torus coordinates
in the same way. This explicitly saturates away the unwanted
denominator components. Elimination uses a bounded number of variables
and equations of bounded degree. Homogenize with one further variable
and saturate by that homogenizing coordinate, equivalently by adjoining
one more inverse variable, before forming coefficient matrices. This
prevents components supported at infinity from entering the affine
elimination ideal. A uniform Groebner-basis degree
bound, such as [Dubé, Theorem 8.2](https://epubs.siam.org/doi/10.1137/0219053),
then reduces this computation to Macaulay coefficient matrices of
bounded size; homogenization keeps the required homogeneous ideal
membership computations within bounded degrees as well. After each
specialization, choose a nonzero pivot minor
appropriate to the actual ranks. Row-reduction coefficients are ratios
of minors of those bounded matrices, so their heights are `O(H+1)`.
This proves the generator bound even when a previously chosen pivot
vanishes. The degree and height bounds for implicitization are also
given by [D'Andrea--Krick--Sombra, Theorem 3.24](https://arxiv.org/abs/1103.4561).

For the second assertion, use the graph including `u`, and work over
the function field of the image. The generic fiber algebra has bounded
dimension. Its multiplication matrices are obtained by the same
bounded-degree elimination and row reduction. Since `u` belongs to
the base function field, its multiplication matrix is scalar; its
trace divided by the fiber degree expresses that scalar as a rational
function in the image coordinates. The matrix entries and this trace
are ratios of polynomial minors of bounded size. Clearing their
denominators gives `A,D` with the stated bounds. This argument is
applied separately to each actual specialization; only the bounds are
uniform. It does not assume that an inverse formula from a generic
member specializes successfully.

Finally, substitution in `D(F)` leaves a nonzero rational function of
bounded degree and coefficient height `O(H+1)` on the fixed curve.
Elimination with the fixed Weierstrass equation, followed by the
polynomial root-height bound, gives the assertion about its zeros.
These particular zeros suffice for the argument below; no uniform
claim about the entire singular locus of every image is needed.

## 3. The rational coordinate and its height

The divisor-vector classification in the phase audit holds for every
tuple of distinct shifts. If the nontrivial complement sums are not
all equal, the phase field is `Q(i)(E)` and take `u_N=x`. If those
sums equal `C_N`, put `T_N=C_NP`. For `T_N!=O` take

```text
u_N(Q)=(y(Q)+y(T_N))/(x(Q)-x(T_N));                   (4)
```

for `T_N=O` take `u_N=x`. This function is invariant under
`Q -> T_N-Q`, has degree two, and generates the quotient field.
It is defined over `Q`, and its coefficient height is `O(H)`.
The elimination lemma therefore supplies
`u_N=R_N(F_N)=A_N(F_N)/D_N(F_N)` with uniformly bounded degrees and
coefficient heights `O_m(H)`.

The height estimate needed here must retain the moving translation.
On the fixed surface `E x E`, the rational slope function

```text
u(Q,T)=(y(Q)+y(T))/(x(Q)-x(T))
```

has polar divisor
`{O} x E + E x {O} + Delta`, where `Delta` is the diagonal. It is a
symmetric ample divisor, with canonical height
`hhat(Q)+hhat(T)+hhat(Q-T)`. The two sections defining the rational
map form a subsystem of this fixed line bundle. Passing to a fixed
very ample multiple bounds the height of that subsystem by the full
projective height. Consequently, wherever the slope is defined,

```text
h(u(Q,T)) <= hhat(Q)+hhat(T)+hhat(Q-T)+O(1),
h(u_N(NP)) <= 2qN^2+O(NK+H).                         (5)
```

The constant in the first inequality is uniform because its surface,
divisor, and rational map are fixed. The case `u_N=x` has the usual
bound `2qN^2+O(1)`.

Equality of `u_N(NP)` with a target `z` of height `o(N^2)` is impossible
under (1). In the slope case the line through `-T_N` with slope `z`
has coefficient height `O(h(z)+H)`. Substituting this line into the
fixed cubic gives a cubic polynomial in `x` with coefficient height
`O(h(z)+H)`. Every root, including `x(NP)`, would have height
`O(h(z)+H)=o(N^2)`, contradicting
`h(x(NP))=2qN^2+O(1)`. In the case `u_N=x`, that contradiction is
immediate. This avoids importing a fixed-translation height error as
an unjustified uniform constant.

## 4. Generic targets

Suppose (2) fails along a subsequence. Then for some fixed `eta>0`
all anchor residues on that subsequence satisfy

```text
log max_(i<m)|t_(N,i,m)|
 <= (B/2-2)qN^2-eta N^2.
```

Equations (3) and the exact phase identity give targets `a_N` in
`Q(i)^(m-1)` with `h(a_N)=O_m(s_N+H)=o(N^2)` and

```text
max_i |F_(N,i)(NP)-a_(N,i)|
 <= exp(-(2q+eta)N^2+o(N^2)).                         (6)
```

The functions and targets have bounded archimedean coordinates: the
phases have modulus one. The defining equations furnished by the
elimination lemma have bounded degree and coefficients of height
`O_m(H)`. Their Lipschitz constants on this compact region are at most
`exp(O_m(H))`. Thus every defining equation evaluated at `a_N` is
exponentially small with a fixed positive quadratic exponent. A
nonzero value would have height `O_m(s_N+H)=o(N^2)` and hence modulus
at least `exp(-o(N^2))`. All defining equations therefore vanish, so
`a_N` is exactly on its corresponding image curve.

If `D_N(a_N)!=0`, the product formula bounds its modulus below by
`exp(-O_m(s_N+H))`. Equation (6) gives the same lower bound, up to a
factor two, at `F_N(NP)`. Applying `R_N=A_N/D_N` consequently yields

```text
z_N=R_N(a_N) in Q(i),       h(z_N)=o(N^2),
|u_N(NP)-z_N| <= exp(-(2q+eta)N^2+o(N^2)).             (7)
```

Write `u_N(NP)=a/b` in reduced rational form with `b>0`, and
`z_N=c/d` in primitive Gaussian form. If they differ,

```text
|u_N(NP)-z_N|=|ad-bc|/(b|d|)
 >= exp(-2qN^2-o(N^2)),                              (8)
```

by (5), because a nonzero Gaussian integer has modulus at least one
and `log|d|<=h(z_N)`. Equality was excluded in Section 3.
This contradicts (7).

## 5. Moving denominator exceptions

It remains to treat `D_N(a_N)=0`. Set
`G_N=D_N(F_N)`, a nonzero rational function on the fixed elliptic
curve of bounded degree and coefficient height `O_m(H)`. On this
subsequence (6) implies

```text
|G_N(NP)| <= exp(-(2q+eta)N^2+o(N^2)).                 (9)
```

It is not necessary to list singular image values. The zeros of
`G_N` already have degree `O_m(1)` and canonical height `O_m(H)`.
If `G_N(NP)=0`, the same bounded-degree root-height estimate would
give `h(x(NP))=O_m(H)=o(N^2)`, which is impossible.

For a nonzero evaluation, apply the moving-target elliptic-logarithm
bound audited in Section 3 of
[the phase audit](cm_denominator_phase_lift_audit.md). Its explicit
height dependence, with target height `O_m(H)`, gives an exponent
bounded by

```text
O_m((H+1) max(H,log(N+2))
        (log(H+log(N+3)+2))^5).                       (10)
```

To pass from distance to the zero set to `G_N(NP)`, factor the
bounded-degree numerator on the fixed curve. Its zero multiplicities
are bounded, all zero heights are `O_m(H)`, and the nonzero
normalizing coefficients have height `O_m(H)`. Algebraic root
separation and the local parameters of the fixed curve therefore
cost at most `exp(O_m(H))`; the denominator has the same bound,
apart from the fixed-orbit archimedean coordinate growth
`exp(O(polylog N))`. The latter follows from the fixed-target
elliptic-logarithm bound recorded in Section 3 of the phase audit,
applied to the fixed pole set of the Weierstrass coordinates. Thus
(10), with an additional `O_m(H+polylog N)` term, also bounds
`-log|G_N(NP)|`.

Under (1), `H=O(N^(1-2epsilon))`, and (10) is `o(N^2)`.
This contradicts (9). Both the generic and exceptional cases are now
excluded, proving (2).

The proof more generally works when `K=o(N)` and the expression in
(10), with `H=1+K^2`, is `o(N^2)`. It does not cover arbitrary shifts
of size `o(N)` without that additional condition. Nor does it make
the constants uniform when `m` grows. The actual-arc radius calculation
from the phase audit changes only by `O_m(NK+H+s_N)=o(N^2)`, so its
positive normalized-arc growth exponents persist in the range (1).
