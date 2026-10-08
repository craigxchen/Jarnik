# A seven-row Christoffel positivity constraint

This is a magnitude-sensitive necessary condition for the punctured
orthogonal moment template. It uses arbitrary real combinations of all
eight rows, and is not a condition on pairwise root orders alone. It does
not exclude every template, and it gives no transfer to arbitrary circles.

The follow-up [relaxation theorem](critical_moment_christoffel_relaxation_limit.md)
characterizes the combined scalar and `s`-subset determinant conditions
for `s>=100`: they are feasible exactly when the centered row space has
no nonzero vector supported on at most `s` columns. The existing admissible
word family passes these conditions with explicit Chebyshev magnitudes.
The full matrix inequality, equivalent to the original moments, remains
unresolved on that family.

Let `s>=2`, `n=4s+2`, and suppose an eight-row sign matrix `T` satisfies
`TT^t=(4s+4)I`. Delete a constant column and a balanced column `b`, obtaining
`S`; normalize the constant column to `1`. Thus `b^t1=0`, `b^tb=8`, and

```
SS^t = (4s+4)I - 11^t - bb^t.
```

Suppose nonzero `a_j` with distinct absolute values satisfy

```
S(a^(2r-1)) = alpha_r 1,       1<=r<=s.
```

Put `x_j=a_j^2`, absorb the signs into `V_ij=S_ij sign(a_j)`, and let
`A,B` be the two four-row groups of `b`. Write

```
h_Aj = sum_(i in A) V_ij,       h_Bj = sum_(i in B) V_ij,
h_j = h_Aj+h_Bj,               d_j = h_Aj-h_Bj.
```

Define the positive definite `s` by `s` Hankel matrix and its evaluation
kernel by

```
H_kl = sum_j x_j^(k+l+1),                  0<=k,l<s,
K_jk = sqrt(x_j x_k) v(x_j)^t H^(-1) v(x_k),
v(x) = (1,x,...,x^(s-1))^t.
```

Define another `n` by `n` matrix explicitly from the signed columns:

```
Q_jk = [V_j^t V_k - h_j h_k/8]/[4(s+1)]
       + d_j d_k/[16(s^2-1)].                         (1)
```

Then the necessary matrix positivity condition is

```
                         K + Q <= I_n.               (2)
```

Here the inequality is in positive semidefinite order. In particular,

```
x_j v(x_j)^t H^(-1) v(x_j) <= 1-q_j,                 (3)

q_j = [8-(h_Aj^2+h_Bj^2)/4]/[4(s+1)]
      + (h_Aj-h_Bj)^2/[32(s-1)].
```

The full matrix condition (2), for fixed `V,x`, is in fact equivalent to
the original common-row odd-moment equations. The scalar and determinant
consequences below are relaxations of it. No sufficiency is asserted for
those weaker consequences or for a discrete path that passes them.

## Proof

Let `E=I-11^t/8`, `C=EV`, and let `U_jk=sqrt(x_j)x_j^k`. The equations are
precisely `CU=0`. Their two Gram matrices are

```
CC^t = (4s+4)E-bb^t,       U^tU=H.
```

On `1^perp`, the first matrix has eigenvalue `4s-4` on the line spanned
by `b`, and `4s+4` on its six-dimensional orthogonal complement. Therefore
`C` has rank seven. Its row-space orthogonal projector is
`Q=C^t(CC^t)^+C`; decomposing these two eigenspaces gives (1).
The column-space projector of `U` is `K=UH^(-1)U^t`. The moment equations
make these two subspaces orthogonal, proving (2). Taking diagonal entries
gives (3), after separating the `b` component from the other six components.

Conversely, if two orthogonal projectors satisfy `K+Q<=I`, then for
`z` in the range of `K`, positivity gives
`0<=z^t(I-K-Q)z=-z^tQz`. Hence `Qz=0`; thus `CU=0`.

Equivalently, a direct general-row-combination proof of (3) starts with
any `mu` satisfying `mu^t1=0` and any polynomial `p` of degree at most
`s-1`. Set `c=C^t mu`, `f_j=sqrt(x_j)p(x_j)`. Since `c^tf=0`,

```
c_j^2 f_j^2 <= (||c||^2-c_j^2)(||f||^2-f_j^2).
```

The maximum of `c_j^2/||c||^2` over all such real row combinations is
`q_j`; the maximum of `f_j^2/||f||^2` over all such polynomials is `K_jj`.
They can be optimized independently, giving `q_j+K_jj<=1`.

## A concrete separation inequality

For every `j` and every real polynomial `p` of degree at most `s-1`,

```
q_j x_j p(x_j)^2
  <= (1-q_j) sum_(k!=j) x_k p(x_k)^2.                (4)
```

In particular, order `|a_1|<...<|a_n|`. If the last signed column is
nonconstant, then `q_n>0`. Every solution has `q_n<1`, and choosing
`p(x)=x^(s-1)` gives

```
(|a_n|/|a_(n-1)|)^(4s-2)
  <= (1-q_n)(n-1)/q_n.                              (5)
```

Thus magnitude gaps cannot be freely assigned to a passing discrete
path. A sharper version of (4) can set `p` equal to the product of
`x-x_k` at any selected `s-1` other nodes, eliminating those terms.

For comparison, restricting `mu` to one row difference gives the
single-pair quantity `1/D_ik` at a column separating that pair, where
`D_ik` is `2s+1` or `2s+2`. If a signed column has two plus entries in
each of `A,B`, then the full multirow value is

```
q_j=2/(s+1),
```

strictly larger than every single-pair value `1/(2s+1)` or `1/(2s+2)`.
In fact this strict comparison holds for every nonconstant signed column:

```
q_j >= (7s-5)/[8(s^2-1)] > 1/(2s+1) > 1/(2s+2).     (7)
```

To prove the first inequality, put `a=h_Aj`, `b=h_Bj`. After multiplication
by `32(s^2-1)`, its difference is

```
[36-(a+b)^2]s -44+3a^2+3b^2-2ab.
```

For the 23 nonconstant types `a,b in {-4,-2,0,2,4}`, the slope is
nonnegative and the value at `s=2` is nonnegative. Equality holds, for
example, at `(a,b)=(4,2)`. The next strict inequality has positive
numerator `6s^2-3s+3` after clearing its positive denominators.
This comparison concerns the Bessel inequalities; it does not claim
that (3) logically strengthens every nonlinear pairwise moment identity.

## A coupled Vandermonde inequality

For every `s`-element column set `J`, restriction of (2) and determinant
monotonicity give

```
 [product_(j in J) x_j]
 [product_(j<k in J) (x_k-x_j)^2]
 ---------------------------------  <= det(I_J-Q_J).       (6)
            det H
```

Indeed the numerator divided by `det H` is `det K_J`. By Cauchy--Binet,

```
det H = sum_(|L|=s) [product_(j in L) x_j]
                    [product_(j<k in L) (x_k-x_j)^2].
```

Thus (6) bounds the fraction of the entire weighted Vandermonde sum
that can concentrate on any one `s`-subset, using the complete seven-row
Gram data. This is a precise all-degree positivity constraint retaining
the unknown magnitudes, rather than another finite-state exclusion.

## Exact algebra audit

[The standard-library Fraction checker](check_critical_moment_multirow_christoffel.py)
verifies the projector formulas on an eight-row partial Walsh fixture
with `s=3`, including symmetry, idempotence, trace seven, and `CQ=C`.
Independent distinct positive square nodes verify the Hankel projector,
all 364 three-subset identities `det K_J` on the left of (6), and the
Cauchy--Binet sum. It also checks the affine numerator proof of (7) for
all 23 nonconstant group-sum types uniformly in `s>=2`.

The chosen nodes do not solve the moment equations. The checker does not
assert `K+Q<=I`, inequality (6), or existence of a distinct-node solution
on that arbitrary-node fixture.
