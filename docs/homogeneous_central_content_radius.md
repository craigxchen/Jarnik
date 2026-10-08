# Homogeneous central content and the intrinsic radius

This note gives a coordinate-free arithmetic formula for the ordinary
central vector of six rational circle nodes.  It separates its exact
primewise primitive normalization from its archimedean size and compares
both pieces with the intrinsic least circle radius.  The formulas are
invariant under integral unimodular changes of frame.

They do not prove the unproved uniform short-arc bound.  The comparison
shows the remaining obstruction explicitly: central primitivization uses
one minimum over six determinant-star valuations, whereas the least radius
uses a positive-part maximum over fifteen edge valuations.

## 1. Homogeneous setup and covariance

Write

```text
Q=[[a,b],[b,c]],       q(v)=v^T Q v,       det(Q)=r^2>0,
```

where `Q` is primitive positive integral and `r` is the positive integer
square root of its determinant.  Let `v_1,...,v_6` be primitive integer
vectors representing six distinct points of `P^1(Q)`, and put

```text
q_i=q(v_i),
Delta_ij=det(v_i,v_j),
D_i=product_(j!=i) Delta_ij.
```

All `q_i` are positive and all `Delta_ij` are nonzero.  The rational
central weights are, up to one common sign,

```text
w_i=q_i^2/D_i.                                         (1)
```

In a chart `v_i=y_i(t_i,1)`, formula (1) differs by the common factor
`1/product_i y_i` from

```text
q((t_i,1))^2/product_(j!=i)(t_i-t_j),
```

which is the usual weighted Vandermonde formula.  This also proves the
formula when one of the nodes is at infinity, either in the other chart or
by homogeneity.

For `S in SL_2(Z)`, make the simultaneous change

```text
Q -> S^T Q S,       v_i -> S^(-1)v_i.                (2)
```

Then every `q_i`, `Delta_ij`, `D_i`, and `w_i` is unchanged.  The new
vectors remain primitive.  Thus (1) is an actual invariant of (2), rather
than a chart-dependent denominator formula.

There is also the expected projective covariance.  Replacing
`v_i` by `lambda_i v_i`, for arbitrary nonzero rational `lambda_i`, gives

```text
w_i -> (product_k lambda_k)^(-1) w_i                 (3)
```

for every `i`.  Indeed, the numerator gains `lambda_i^4`, while `D_i`
gains `lambda_i^5 product_(j!=i)lambda_j`.  Hence the primitive integer
vector represented by `(w_i)` is independent of all projective choices.

## 2. Exact all-prime primitive normalization

For every rational prime `p`, define

```text
alpha_i(p)=2 v_p(q_i)-sum_(j!=i) v_p(Delta_ij),
m_p=min_i alpha_i(p).                                 (4)
```

Only finitely many of these integers are nonzero.  Put

```text
kappa=product_p p^(-m_p).                             (5)
```

Then the primitive central vector and its height are exactly

```text
u_i=kappa w_i,
U=max_i |u_i|
 =kappa max_i q_i^2/|D_i|.                            (6)
```

The signs in (6) can all be reversed.  The proof is just valuation by
valuation:

```text
v_p(u_i)=alpha_i(p)-m_p>=0,
min_i v_p(u_i)=0.                                     (7)
```

Thus all coordinates are integers and their gcd is one.  Conversely these
conditions uniquely specify the positive rational multiplier `kappa`.

This identifies exactly what the specialized affine weight content was
measuring.  Let `ell_W` be the least positive integer that clears the
denominators of the six weights, and put

```text
C_W=gcd_i(ell_W w_i).
```

Then

```text
kappa=ell_W/C_W,
v_p(ell_W)=max(0,-m_p),
v_p(C_W)=v_p(ell_W)+m_p.                             (8a)
```

The evaluation-matrix clearer `ell_E` in the moving-base notes can be
strictly larger: it clears every entry `t_i^k/D_i`, not only the six
specialized weights.  Nevertheless, if

```text
C_E=gcd_i(ell_E w_i),
```

then, for the same affine representatives,

```text
kappa=ell_E/C_E,
v_p(C_E)=v_p(ell_E)+m_p.                             (8b)
```

Only (8a) characterizes the minimal weight clearer.  For example, the
balanced `S={3,4,5}` cut in the fixed-Gram family has `m_p=0` and hence
`v_p(ell_W)=0`, but `v_p(ell_E)=v_p(C_E)=2e`.

The two large chart-dependent evaluation integers can therefore cancel
completely.  Their quotient is the normalizer for the fixed affine
representatives.  The multiplier `kappa` itself is covariant: changing the
representatives as in (3) changes `kappa` by the reciprocal common scalar
needed to keep `u` fixed.  The primitive vector and the valuation
differences (7), rather than `kappa` alone, are intrinsic.  Formulas
(8a)--(8b) include primes in the node denominators, the determinant factors,
accidental evaluation content, and the prime two.

Under (3), every `alpha_i(p)` and `m_p` is shifted by the same number
`-sum_k v_p(lambda_k)`.  Hence the differences in (7), and therefore the
primitive vector, are projectively invariant.  This is the primewise form
of the common-scalar cancellation.

## 3. Exact archimedean chord identity

There is a canonical rational unit-circle realization.  For
`v=(x,y)`, define

```text
h_v=a*x+(b+i*r)*y,       zeta_v=h_v/conjugate(h_v).
```

The identity

```text
|h_v|^2=a q(v)
```

shows that `zeta_v` has modulus one.  For two nodes, their normalized
chord length is

```text
delta_ij=|zeta_i-zeta_j|
        =2r |Delta_ij|/sqrt(q_i q_j).                 (9)
```

Multiplying (9) over the five edges incident to `i` and using (1) gives
the exact star identity

```text
|w_i| product_(j!=i) delta_ij
   =(2r)^5/sqrt(product_k q_k).                       (10)
```

The right side is independent of `i`.  Consequently

```text
U=kappa (2r)^5 /
  (sqrt(product_k q_k) min_i product_(j!=i)delta_ij). (11)
```

Thus the archimedean part of the central height is exactly a reciprocal
five-chord star product.  Formula (11) is invariant under scaling `Q` or
the individual node representatives: the changes in `kappa` cancel the
changes in the displayed common factor.

## 4. The intrinsic least-radius lcm in the same data

Let

```text
B_ij=v_i^T Q v_j,
g_ij=gcd(|B_ij|,r|Delta_ij|),
A_ij=B_ij/g_ij,       T_ij=r Delta_ij/g_ij,
epsilon_ij=2 if A_ij and T_ij are both odd, and 1 otherwise.
```

The binary Gram determinant identity is

```text
B_ij^2+r^2 Delta_ij^2=q_i q_j.                       (12)
```

It removes the apparent dependence on the polar value.  Indeed,

```text
g_ij^2=gcd(B_ij^2,r^2 Delta_ij^2)
      =gcd(q_i q_j,r^2 Delta_ij^2).                 (13)
```

The first equality uses `gcd(x,y)^2=gcd(x^2,y^2)`, and the second uses
`gcd(x,y)=gcd(x+y,y)`.  Thus `g_ij`, despite its definition through
`B_ij`, is determined by `q_i,q_j,r`, and `Delta_ij` alone.

Hence

```text
n_ij=(A_ij^2+T_ij^2)/epsilon_ij
    =q_i q_j/(g_ij^2 epsilon_ij)                    (14)
```

is a positive integer.  It is the norm of the reduced Gaussian
denominator of `zeta_j/zeta_i`: that quotient has the form

```text
(A_ij+i T_ij)/(A_ij-i T_ij),
```

and coprime real and imaginary parts have a common Gaussian factor with
norm two exactly when they are both odd.

Let `odd(m)` denote the largest odd divisor of a positive integer `m`.
Since coprime `A_ij,T_ij` have a sum of squares that is odd unless both are
odd, in which case it is exactly twice an odd integer, (13)--(14) give the
polar-free formula

```text
n_ij=odd(q_i q_j/gcd(q_i q_j,r^2 Delta_ij^2)).       (15)
```

This integer is also exactly the reduced denominator of the squared
normalized chord

```text
delta_ij^2=4r^2 Delta_ij^2/(q_i q_j).
```

After the gcd in (15) is removed, the denominator has `2`-valuation at
most one by the coprime sum-of-two-squares description.  The numerator
factor `4` removes that entire `2`-part, while every odd part remains
uncancelled.  Thus (16) can equivalently be read as

```text
N=lcm of the reduced denominators of all delta_ij^2.
```

The exact all-edge least-radius theorem now gives

```text
N=lcm_(i<j) n_ij,       R_min=sqrt(N).               (16)
```

Here `R_min` is the radius of the primitive Gaussian integral realization
of the six relative phases.  In primewise form,

```text
v_p(N)=max_(i<j) max(0, v_p(q_i)+v_p(q_j)
                         -2v_p(r)-2v_p(Delta_ij))     (p odd),
v_2(N)=0.                                            (17)
```

The nonnegative part in (17) comes from the gcd in (15).  The second line
is the complete ramified correction: every reduced edge-denominator norm
is odd.

Equations (9), (14), and (16) also give the exact squared chord in a
least-radius realization:

```text
N delta_ij^2=(4/epsilon_ij) T_ij^2 (N/n_ij).         (18)
```

In particular every nonzero chord gives `delta_ij>=N^(-1/2)`, while
`delta_ij<=2`.  Combining this elementary range with (11) yields

```text
kappa (2r)^5/(32 sqrt(product q_k))
 <= U <=
kappa (2r)^5 N^(5/2)/sqrt(product q_k).              (19)
```

This is an exact-data comparison, not a useful uniform lower bound, because
the homogeneous content multiplier `kappa` remains present.

## 5. Full-profile exhaustion makes the fair-table totals global

There is a useful scoped upgrade of the core calculation in
[`moving_base_fair_content_table.md`](moving_base_fair_content_table.md).
Use its five finite integer cotangents `X_i/L`, and suppose all of the
following hold along a hypothetical family:

```text
log|X_i|=8w+o(w)                    for every i,
log L=o(w),
log N(H_S)=w+o(w)                  for all nonempty S subset [5]. (20)
```

Here the `H_S` are the independent exact good-cut core blocks: `H_S`
divides `X_i+iL` when `i in S`, and its rational norm contributes to
`X_i^2+L^2` and to `X_i-X_j` when `{i,j} subset S`, with the orientations
used in the fair-content table.  The exact blocks have independent prime
support, so their indicated products divide the corresponding integers.

For fixed `i`, exactly 16 of the 31 nonempty subsets contain `i`.  Therefore
the core contributes

```text
16w+o(w)
```

to `log(X_i^2+L^2)`.  The first two lines of (20) give
`log(X_i^2+L^2)=16w+o(w)`, so the total valuation mass outside those core
divisors is `o(w)`.  Similarly, exactly eight subsets contain a fixed pair
`{i,j}`.  Their independent product divides the nonzero integer
`X_i-X_j`, giving

```text
log|X_i-X_j|>=8w+o(w).
```

On the other hand,
`|X_i-X_j|<=|X_i|+|X_j|=exp(8w+o(w))`.  Hence

```text
log|X_i-X_j|=8w+o(w),
```

and the total remaining difference valuation mass is `o(w)`.  Notice that
equal individual heights alone would not give this conclusion; the lower
bound uses the full exact core divisibility.  These exhaustion statements
use only exact divisibility and nonnegativity of all prime valuations.

The formulas for `v_p(ell_E)` and `v_p(C_E)` are fixed minima and maxima
of finitely many integer linear forms in these norm and difference
valuations.  They are Lipschitz with an absolute constant.  At primes
dividing `2L`, the deviation from the ideal cut formulas has total size
`O(log(2L))=o(w)`, and exhaustion makes the remaining noncore valuation mass
`o(w)` as well.  This does not assert that the total valuation at those
primes is small if some exact core mass is supported there.  Consequently
the core totals from the fair table exhaust the global quantities:

```text
log ell_E=29w+o(w),       log C_E=38w+o(w),
log kappa_aff=-9w+o(w).                              (21)
```

Here `kappa_aff=ell_E/C_E` is the covariant normalizer for this specified
affine frame and its specified representatives.

This also determines the previously conditional evaluation scale.  Select
the first two rows, so `delta=X_2-X_1` and `d|L`.  Every raw normalized
weight has logarithmic size

```text
2 log|delta/d|+2 log(X_i^2+L^2)
  -sum_(j!=i)log|X_i-X_j|=16w+o(w),                 (22)
```

with the infinity row giving the same value.  Since
`Nhat*f=ell_E*(W_1,...,W_5,W_infinity)`, while
`log H(Q0)^2=16w+o(w)`, one gets

```text
log ||Nhat*f||_infinity=45w+o(w),
log R_eval=29w+o(w),
log U=7w+o(w).                                      (23)
```

Thus the full fair profile does make the `29w` and `38w` core counts global
under (20), but it does not make `R_eval` subexponential.  Rather, the
evaluation scale carries exactly the same `29w` denominator-clearer growth.

## 6. Why the two local formulas do not close

The central formula (4) and radius formula (17) use exactly the same local
input.  Put

```text
x_i=v_p(q_i),       d_ij=v_p(Delta_ij),       rho=v_p(r).
lambda_ij=x_i+x_j-2rho-2d_ij.                         (24)
```

Let `s_i=sum_(j!=i)lambda_ij` be the signed degree of vertex `i` in this
edge-labelled complete graph.  If `X=sum_i x_i`, then direct summation gives

```text
s_i=4x_i+X-10rho-2sum_(j!=i)d_ij
   =2alpha_i+X-10rho.                                (25)
```

Consequently there is an exact primewise cancellation identity:

```text
v_p(u_i)=alpha_i-m_p=(s_i-min_k s_k)/2,              (26)
v_p(N)=max_(i<j) max(0,lambda_ij)                    (p odd),
v_2(N)=0.
```

The edge defects `lambda_ij` are unchanged by scaling `Q` or by individually
rescaling the node representatives.  Formula (26) is therefore a fully
projective local description: primitive central content records the centered
signed vertex degrees, while the radius records the largest positive edge
defect.

This is sharper than treating the polar gcds as independent data.  It also
pinpoints why the radius alone does not determine the central height.  It
controls every positive `lambda_ij`, but it does not control the negative
edge defects or their signed row sums.  A further bound would have to use the
global compatibility of all fifteen determinants of six primitive vectors,
together with the positive binary form and the archimedean short-arc
condition.

This is a real obstruction rather than an omitted chart cancellation.
The fixed-moment elliptic family in
[`six_point_fixed_moment_elliptic_obstruction.md`](six_point_fixed_moment_elliptic_obstruction.md)
has the constant primitive central vector

```text
(5,-5,5,-5,-1,1),       U=5,
```

while its angular span tends to zero and `N` tends to infinity.  Thus no
bound `U >= N^c` with fixed `c>0`, or any other proper lower growth of `U`
in `N` alone, can follow from the homogeneous formulas.  That family is not
a short-arc counterexample because its radius growth compensates for the
shrinking angle.

Likewise, the fixed-Gram family in
[`moving_base_fair_content_table.md`](moving_base_fair_content_table.md)
has unbounded affine `C_E` with fixed primitive `Q`.  Formula (8b) explains
why this does not contradict invariance: the evaluation-matrix denominator
clearer grows at the same primes and cancels in its quotient with `C_E` for
those representatives.

Finally, a small determinant does not bound rational-node height.  The
primitive vectors `(T,1)` and `(T+1,1)` have determinant `-1` for every
integer `T`, while their heights tend to infinity.  Reducing `Q` by (2)
can therefore transfer height into the node vectors without changing any
quantity in (4), (10), or (17).

The usable conclusion is scoped.  Formulas (6), (11), (17), and (26) give exact,
simultaneously invariant arithmetic and archimedean descriptions of `U`
and `N`.  A proof of the uniform short-arc statement still needs a new
compatibility principle linking the two extrema in (26) with the small
chord-star products in (11).  Covariance alone does not supply that link.

The companion checker
[`check_homogeneous_central_content_radius.py`](check_homogeneous_central_content_radius.py)
verifies the valuation normalization, unimodular covariance, projective
rescaling, Gram and chord identities, edge denominator formula, and the
all-edge Gaussian least-radius reconstruction on exact examples.  It also
checks the `p=13` balanced-cut regression distinguishing `ell_E` from
`ell_W`.
