# Exact degree-twenty radius and all-prime content of the circle cover

For the integer isotropic pencils in
[the fixed-weight cover construction](six_point_isotropic_circle_cover.md),
the primitive squared radius is an exact degree-twenty expression divided
by explicitly specified contents. No exceptional prime valuation is discarded.
For a fixed frame whose weights have no vanishing proper subsum, these
contents are bounded by actual homogeneous Bezout constants. This proves
`N asymp H^20` on its real circle cover. Near angular coalescence, the
normalized endpoint scale is `asymp |d|/H`, where `d^2=Delta(s,t)` is the
degree-twelve cover equation.

Every comparison constant in this note depends on the fixed cleared frame.
The bounds are not uniform when the central weights or frame vary with the
extracted fair core.

The [prime-cut separation](central_subsum_cut_separation.md) subsequently
verifies the nonvanishing-subsum condition for a fixed six-label full fair
core with subpower errors. The [local branch calculation](six_point_simple_branch_genus.md)
then proves that `Delta` is squarefree and both covers have genus five.
These results do not remove the dependence of the comparison constants.

## 1. Raw cofactor conventions and the universal factorization

Fix an integer three-by-six pencil matrix

```text
V(s,t)=(1,x(s,t),y(s,t))^T,
```

where `x_i,y_i` are homogeneous linear forms. Its rowspace is maximal
isotropic for the fixed diagonal form `diag(u)`, and contains `1`.
Choose the first five feature rows

```text
(1,x_i,y_i,x_i^2,x_i y_i,y_i^2),       i=0,...,4,
```

and use their RAW signed five-by-five cofactors, with sign `(-1)^j` for
omitted column `j=0,...,5`, to define

```text
(c0,c1,c2,A,B,C),
Delta=4AC-B^2,
K=C c1^2-B c1 c2+A c2^2-c0 Delta.                    (1)
```

Put `m_I=det(V_I)` in increasing column order. Then the following polynomial
identity is exact:

```text
K = - product_(I subset {0,...,4}, |I|=3) m_I.        (2)
```

This identity holds for arbitrary five points, before imposing isotropy.
To prove it, homogenize their coordinates in `P^2`. The determinant of the
cofactor conic has degree six in each point. It vanishes whenever any chosen
triple is collinear: generically the unique conic is that line times the
line through the other two points. Each collinearity determinant is an
irreducible polynomial, so their product divides the conic determinant.
The ten determinants have the same multidegree, proving proportionality.
For the cofactor convention above, the points
`(0,0),(1,0),(0,1),(1,1),(2,3)` give `K=48` and the product `-48`, fixing
the sign. Equivalently, `K=-4 det(Q)` for the symmetric conic matrix `Q`.

Each `m_I` is homogeneous quadratic in `(s,t)`. There are ten complementary
`3|3` partitions, and exactly one member of each partition avoids label five.
Self-duality of the isotropic rowspace makes complementary minors proportional
by fixed nonzero rational factors. In squared form their relation is

```text
m_I^2/m_(I^c)^2 = -product_(j in I^c)u_j/product_(i in I)u_i.
```

Thus (2) is also the product over all ten balanced partitions, with exact
representatives fixed as above. If the cofactor vector is divided by a
polynomial content `h0`, its `K` is divided by `h0^3`; this factor may not
be silently removed from (2).

## 2. The affine coefficient ideal determines the Gaussian row content

Evaluate at coprime integers `(s,t)`, and suppose

```text
d^2=Delta(s,t)>0,
```

with a nonsingular conic and six distinct columns. Since `Delta` is an
integer, rational `d` is an integer. The exact Gaussian lift is

```text
Z_i=alpha x_i+beta y_i+gamma,
alpha=2Ad,
beta=d(B+i d),
gamma=d c1+i(2Ac2-Bc1).                             (3)
```

All rows have norm `4AK>0`. Define

```text
h=gcd(A,B,C)>0,
g0=gcd_G(2A,B+i d),
eta=[2Ac2-(B+i d)c1]/g0,
e=gcd_G(d,eta).                                     (4)
```

The quotient defining `eta` is Gaussian integral. The coefficient ideal has
Gaussian gcd exactly `g0 e`, up to a unit. Indeed, put
`a'=2A/g0`, `b'=(B+i d)/g0`. These are coprime, and

```text
gamma=i[2Ac2-(B+i d)c1].
```

After dividing the three coefficients by `g0`, their gcd is consequently
`gcd_G(d a',d b',i eta)=gcd_G(d,eta)`.

Let `g_V` be the positive ordinary gcd of all twenty maximal minors of the
evaluated integer matrix `V`, and let `G=gcd_G(Z_0,...,Z_5)`. There is a
Gaussian integer `f` such that

```text
G=g0 e f,       f divides g_V in Z[i].               (5)
```

For proof, every row is a linear combination of the three coefficients,
so their gcd divides `G`. Conversely, the adjugate of each three-column
submatrix shows that its determinant times every coefficient lies in the
row ideal. An ordinary Bezout combination of all the minors then gives
`g_V` times every coefficient in that ideal. This proves (5), including all
primes dividing the frame or parameter denominators.

## 3. The metric norm identity includes the prime two

The exact identity is

```text
N_G(g0)=4|A|h.                                     (6)
```

Here `B` and `d` must both be even: their squares sum to `4AC`, and two odd
squares sum to two modulo four. Put `b=B/2`, `c=d/2`, so `AC=b^2+c^2`.
The Gaussian ideal `(A,b+i c)` is the rank-two ordinary lattice generated
by `(A,0),(0,A),(b,c),(-c,b)`. Its index, computed from all two-by-two minors,
is

```text
|A| gcd(A,b,c,C).
```

The norm relation gives `gcd(A,b,c,C)=gcd(A,2b,C)=h`. At odd primes this
follows directly by comparing valuations in `AC=b^2+c^2`. At two, if both
`A,C` are divisible by `2^r`, the same relation forces `b,c` divisible by
`2^r`: a smaller common valuation would leave either one odd square or
two odd squares, of valuation zero or one, instead of valuation at least
two. This also proves equality at two. Finally `g0` is twice a Gaussian
gcd of `(A,b+i c)`, proving (6).

Divide (3) by its full Gaussian gcd. Its least primitive squared radius is
therefore exactly

```text
N=R^2=|K|/[h N_G(e) N_G(f)].                         (7)
```

In particular `N` divides `|K|`. The denominator in (7) is an integer
divisor of `|K|`, rather than an unspecified factor supported on a finite
set of primes. Moreover

```text
h divides d,       N_G(e) divides d^2=Delta,
N_G(f) divides g_V^2.                               (8)
```

Thus at a prime not dividing `d g_V`, the valuation of the primitive squared
radius equals the valuation of `K`. Equations (4)--(8), rather than a
finite-support assertion, describe every remaining valuation.

## 4. Actual fixed-frame Bezout bounds

The twenty homogeneous minor quadratics span `s^2,st,t^2` over `Q`.
Indeed, in hyperbolic-frame coordinates the exterior product of the pencil
is a quadratic Veronese with three independent coefficients; an invertible
rational change of ambient basis preserves their independence. Choose a
positive integer `E_V` clearing a rational left inverse of their coefficient
matrix. Then

```text
E_V s^2, E_V st, E_V t^2
```

are integer linear combinations of the minors. At every primitive integer
parameter pair, it follows that `g_V` divides `E_V`.

Now impose the fixed-weight condition

```text
sum_(i in S)u_i != 0
for every nonempty proper subset S of {0,...,5}.      (9)
```

Since the total sum is zero, it suffices to check singletons, pairs and
triples. No pair summing to zero means that the columns never coincide,
even over the algebraic closure; no four can be collinear. The argument is
given in Section 2 of the cover note. Consequently every five columns
determine a unique conic and the raw cofactor vector has no common zero.

Under (9), the binary forms `Delta` and `K` have no common projective zero.
For if `K=0`, the conic consists of two lines, each containing exactly three
of the six points. If also `Delta=0`, the lines are parallel. Their affine
linear equations then give a vector in `W` which is zero on one triple and
one on the other. Its squared isotropy forces that triple's weight sum to
vanish, contradicting (9). A double line or a line at infinity would force
all the affine columns to be collinear and is already excluded.

Choose an explicit positive integer `E_R` and exponent `r` such that

```text
E_R s^r, E_R t^r belong to the integer ideal (Delta,K).
```

Such certificates are obtained by the homogeneous resultant or rational
linear algebra followed by clearing denominators. Evaluation at a primitive
parameter pair proves

```text
gcd(|Delta(s,t)|,|K(s,t)|) divides E_R.               (10)
```

This bounds valuations, not merely prime support. By (7)--(8), at each prime
the exponent in `h N_G(e)` is at most both `v_p(K)` and `3v_p(Delta)/2`.
It follows that

```text
1 <= h N_G(e) <= E_R^(3/2),
|K|/(E_V^2 E_R^(3/2)) <= N <= |K|.                  (11)
```

Both constants come from the actual fixed frame. There is no claim that
they are bounded or subpower in the height of a moving central vector.

## 5. Degree-twenty growth and the remaining square-value problem

On the real parameter set `Delta>=0`, the homogeneous form `K` does not
vanish under (9). At a real zero of `K`, the two component lines above are
real. Nonparallel real lines give `Delta<0`; parallel lines are excluded
by (9). Therefore on the compact set

```text
max(|s|,|t|)=1,       Delta(s,t)>=0,
```

the absolute value of `K` has a positive minimum `k_min`, provided this
set is nonempty. Let `k_max` be its maximum. For integer parameters of height
`H=max(|s|,|t|)`, equations (2) and (11) give

```text
[k_min/(E_V^2 E_R^(3/2))] H^20 <= N <= k_max H^20.   (12)
```

The angular statement retains the branch coordinate as well. Near a real
parameter with `Delta=0`, use the `A` chart or interchange `x,y` so that
`A!=0` there. At the branch, `K!=0` and

```text
delta=2Ac2-Bc1 !=0,       delta^2=4AK.
```

After normalizing `(s,t)` to height one, (3) has common limiting value
`i delta`, while

```text
Z_i-Z_j=d[2A(x_i-x_j)+B(y_i-y_j)]
        +i d^2(y_i-y_j).
```

At least one of the first-order real brackets is nonzero; otherwise all
six affine columns lie on one line. Hence their angular span `theta`
is bounded above and below by positive fixed multiples of the normalized
branch coordinate. Restoring homogeneity yields

```text
theta asymp |d|/H^6,
theta N^(1/4) asymp |d|/H.                          (13)
```

All constants depend on the fixed frame and branch neighborhoods. Outside
those neighborhoods the six distinct circle nodes have angular diameter
bounded away from zero. Thus any bounded-endpoint-scale sequence with this
fixed frame and `H->infinity` must solve the exact square-value problem

```text
d^2=Delta(s,t),       |d|=O(H),       gcd(s,t)=1.      (14)
```

The degree-twelve square equation and degree-twenty radius have both survived
primitive normalization. What remains for the original fair-six task is
control uniform in the moving weights, their isotropic frame, and the
resultant constants; (12)--(14) alone do not provide that control.

## 6. All twenty conic-degeneration parameters are simple

Under (9), `K` is a squarefree homogeneous binary form of degree twenty.
Its zeros are disjoint from those of `Delta`, and each of the ten balanced
triangle quadratics has two simple roots, distinct from the roots belonging
to every other balanced partition.

Here is a local proof, over an algebraically closed field of characteristic
zero. At a zero of `K`, the unique conic is a union of two nonparallel lines,
with three distinct marked points on each. No point is their intersection:
that would force four marked points on one of the two components. Make a
fixed affine coordinate change so the lines are `x=0` and `y=0`. Let `S`
index the three points on `y=0`, so `x_i!=0` there; on `S^c`, `x_i=0`
and `y_i!=0`.

Choose a hyperbolic dual frame for `E=(1,x,y)`. A local parameter on its
isotropic ruling has tangent

```text
x'=f2,       y'=-f1,
beta(x,y')=-1,       beta(y,x')=1.
```

The rank-two matrix `(1,x,0)_S` has its unique column relation proportional
to `(u_i x_i)_(i in S)`, since isotropy gives
`sum_S u_i x_i=sum_S u_i x_i^2=0`. All these relation coefficients are
nonzero. The derivative of its triangle determinant is

```text
m_S'=det(1,x,y')_S.
```

The cofactor vector of the first two rows is a nonzero multiple of that
column relation, and

```text
sum_(i in S)u_i x_i y_i'=beta(x,y')=-1.
```

Therefore `m_S'!=0`. The same argument applies to the complementary
minor, or follows from its fixed nonzero proportionality to `m_S`.
The skew-matrix parameter is an actual local coordinate on the ruling,
so this proves simplicity in any parameter chart, including infinity.

No different balanced partition can vanish at this parameter. A triple
using both components contains two points on one line; if its third point
were collinear with them, it would be the intersection point, already
excluded. Thus exactly one of the ten factors in (2) vanishes, and its
zero is simple. This proves the assertion without a genericity assumption.

The independent checker `check_six_point_simple_conic_degenerations.py`
constructs rational two-line configurations and differentiates their
triangle determinants exactly. It checks both complementary derivatives
and that all other triangle determinants are nonzero.

## 7. Exact verification

The extended `check_six_point_isotropic_circle_cover.py` verifies (2) as a
rational polynomial identity on both degree-twelve squarefree witnesses,
checks `gcd(Delta,K)=1`, and checks the exact metric, coefficient, row-content,
and least-radius formula at their seed configurations. It additionally tests
1,770 metric/coefficient-ideal cases, including ramification at two. These
checks supplement the identities and Bezout arguments above.
