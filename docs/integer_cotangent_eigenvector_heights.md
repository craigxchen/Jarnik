# Primitive cotangent eigenvector heights and the fair-core scale

The spectral determinant audit removes every conductor prime outside
`2L` after column normalization. This note computes the individual
primitive column heights, including the exact bounded correction at
the remaining primes. The result retains the actual least circle
radius, but gives no new lower bound for that radius.

Use the notation and consistently oriented cotangent matrix of
[integer_cotangent_spectral_audit.md](integer_cotangent_spectral_audit.md).
There are `n` labels, including infinity. Write the least primitive
Gaussian integral circle tuple as `z_0,...,z_(n-1)`, with common norm
`N=R^2`, and put

```text
r_i=conjugate(z_i)/conjugate(z_0),
w_i=1/product_(j!=i)(r_i-r_j),
u_r=c_r (w_i r_i^r)_i,             0<=r<n,
```

where `u_r` is Gaussian integral and primitive. Its multiplier `c_r`
is uniquely determined up to a Gaussian unit.

## 1. An intrinsic integer for each column

For each split rational prime `p=pi conjugate(pi)` dividing `N`, let
`e_1<=...<=e_n` be the sorted valuations `v_pi(z_i)`. If `v_p(N)=a`,
primitivity gives `e_1=0,e_n=a`. For `0<=m<=n`, define

```text
H_m=sum_(p|N) [sum largest m e_i - sum smallest m e_i] log p,
K_m=exp(H_m).
```

The definition is independent of the chosen orientation over `p`,
since conjugation replaces each exponent by `a-e_i`. It gives ordinary
positive integers with

```text
K_0=K_n=1,       K_1=K_(n-1)=N,       K_m=K_(n-m).
```

There is also a factorization-free definition. Form every product of
`m` distinct coordinates `z_i`, and take their Gaussian lcm `M_m` and
gcd `G_m`. Then

```text
K_m=|M_m/G_m|.                                        (1)
```

The quotient has equal valuations at the two orientations over each
rational prime, so it is a Gaussian unit times an ordinary positive
integer. This proves (1), including arbitrary prime powers.

Put

```text
S^2=sum_i |w_i|^2
   =sum_i product_(j!=i) (Q_ij^2+L^2)/(4L^2).
```

Every raw eigenvector has squared Euclidean norm `S^2`, because all
nodes have modulus one.

**Exact column-height formula.** For column `r`, put `m=n-r-1`. There
is a positive ordinary integer `d_r` such that

```text
d_r=K_m |c_r|^2,        d_r | (2L)^(2(n-1)),
||u_r||_2^2=S^2 d_r/K_m.                              (2)
```

In particular,

```text
S exp(-H_m/2) <= ||u_r||_2
             <= (2L)^(n-1) S exp(-H_m/2).             (3)
```

## 2. The error at every prime is bounded explicitly

For any Gaussian prime, let `v` be its valuation and set `b=v(2L)`.
Cotangent integrality implies the two-sided estimate

```text
min(v(r_i),v(r_j)) <= v(r_i-r_j)
                  <= min(v(r_i),v(r_j))+b.            (4)
```

Only equal node valuations require proof. If the difference has
valuation more than `v(r_i)+v(2)`, the sum has valuation exactly
`v(r_i)+v(2)`. Then integrality of
`Q_ij=iL(r_i+r_j)/(r_i-r_j)` gives the asserted upper bound.
If the extra valuation is at most `v(2)`, the bound is immediate.

Let the ordered node valuations be `f_1<=...<=f_n`. Replacing all
pair-difference valuations by their minima gives the smallest raw
column valuation

```text
-sum_(j=1)^m f_j,
```

as proved in the spectral audit. Each coordinate's actual weight
valuation differs from this model by a number in `[-(n-1)b,0]`.
Consequently the primitive multiplier satisfies

```text
v(c_r)=sum_(j=1)^m f_j + eta_pi,
0<=eta_pi<=(n-1)v_pi(2L).                             (5)
```

For a prime over `p|N`, the node valuations are `e_0-e_i`; at the
conjugate prime they are `e_i-e_0`. Adding the two baseline sums in
(5) gives minus the spread of the largest and smallest `m` allocation
exponents. Primes outside `N` have baseline zero. Taking Gaussian
norms therefore proves that `K_m |c_r|^2` is an integer, and (5)
shows that it divides `(2L)^(2(n-1))`. This proves (2) at split,
inert, and ramified primes alike.

## 3. The central columns have exact conjugacy and moment identities

Let `U=product_i r_i`. Since the nodes lie on the unit circle,

```text
conjugate(w_i r_i^r)
    =(-1)^(n-1) U w_i r_i^(n-2-r).                   (6)
```

If `n=2h`, column `r=h-1` is therefore collinear with its conjugate.
After Gaussian primitive normalization it is, up to a Gaussian unit,
an ordinary primitive integer vector. A non-axis rational Gaussian
line fixed by conjugation up to a unit would have slope `1` or `-1`;
every Gaussian integer on either of those lines is divisible by
`1+i`, contradicting primitivity. Thus only the real or imaginary
axis can occur for this column.

The Lagrange identities `sum_i w_i r_i^j=0` for `0<=j<=n-2` give,
after choosing the central column real,

```text
sum_i (u_(h-1))_i z_i^j=0,          0<=j<=h-1.         (7)
```

The negative moments give the conjugate equations. This is an exact
integer relation among the first `h-1` circle powers, not merely a
small-angle approximation.

If `n=2h+1`, the two columns `h-1,h` are conjugate up to a Gaussian
unit after primitive normalization. Their Euclidean norms are exactly
equal. All columns have proportional coordinate absolute-value
profiles, since `|r_i|=1`.

There is additional exact support information for the even central
column. At a good prime `p` not dividing `2L`, let `e_i` be the
physical point's allocation exponent, and retain `H_h(p)` for the
unweighted largest-minus-smallest sum at this prime. The ordinary
integer central coefficient has valuation

```text
v_p((u_(h-1))_i)
    = [sum_j |e_i-e_j|-H_h(p)]/2.                     (8)
```

Indeed the norm of the raw coordinate contributes
`sum_j |e_i-e_j|`, while its primitive multiplier contributes
`-H_h(p)` by (5). The column is an ordinary integer vector, so the
Gaussian norm doubles its ordinary prime valuation. The number
`H_h(p)` is also the minimum of `sum_j |e_j-t|` over real `t`,
attained at every median between the middle two exponents.

For a two-level cut of size `s`, (8) means that the coefficients on
the smaller side contain `p^(|h-s|a)`, and the others contain no
factor `p`. A perfectly balanced cut contributes nothing to any
central coefficient. In the eight-row full profile, outside primes
dividing `2L`, coefficient `i` is therefore a product of the norms
of the minority singleton cuts containing `i` cubed, of the minority
pair cuts containing `i` squared, and of the minority triple cuts
containing `i` to the first power. The balanced four-cuts are absent.
This statement retains arbitrary prime-power layers by (8).

## 4. Explicit scale on a hypothetical full fair profile

Suppose an actual sequence has one cut for each nonempty subset of
the `n-1` nonanchor labels, each of log norm `w+o(w)`, and all the
residue-normalization errors have `log L=o(w)`. Suppose also its
actual pair cotangents have the extracted common height

```text
log |Q_ij|=2^(n-3)w+o(w).
```

These are hypotheses about an extracted integer sequence, not an
assertion that such short-arc sequences exist. Then

```text
log N=(2^(n-1)-1)w+o(w),
log S=(n-1)2^(n-3)w+o(w),
H_m=w sum_(s=1)^(n-1) binom(n-1,s) min(m,n-m,s,n-s)+o(w).
```

Formula (2) consequently computes every primitive column height.
For `n=8`, in column order `r=0,...,7`, the coefficients of `w` in
`log ||u_r||_2` are exactly

```text
(321/2, 101, 111/2, 38, 111/2, 101, 321/2, 224).
```

The smallest column has height `exp(38w+o(w))`, whereas the pair
half-angle height is `exp(32w+o(w))`. Integer nonvanishing therefore
provides no contradiction, even for the central real column (7).

More generally, if `n=2h`, the minimum coefficient is

```text
(h/4) binom(2h,h)-2^(2h-3).
```

If `n=2h+1`, it is

```text
((2h+1)/4) binom(2h,h)-2^(2h-2).
```

For `n=4,5,6,7,8,9,10`, these coefficients are respectively

```text
1, 7/2, 7, 19, 38, 187/2, 187.
```

Relative to the pair-cotangent exponent `2^(n-3)`, the minimum height
already exceeds that exponent at `n=7,8`, by the factor `19/16`, and
increases further as the number of labels grows. The primitive spectral
columns do not become short integer vectors on the growing fair core.

## 5. Fixing the central coefficients gives a curve, at their actual height

Let `n=2h` and fix the nonzero ordinary integer central coefficients
`u_i`. On the locus of distinct nonzero complex nodes, the equations

```text
sum_i u_i q_i^j=0,          j=-(h-1),...,-1,1,...,h-1,
q_0=1                                                    (9)
```

define a smooth complex curve whenever they have a solution. Before
fixing `q_0`, their Jacobian has rank `n-2`. A dependence among its
rows, after multiplying each column by its nonzero `q_i/u_i`, would
give a Laurent polynomial with exponents `-(h-1),...,-1,1,...,h-1`
vanishing at all `n` distinct nodes. Multiplication by `q^(h-1)`
produces a polynomial of degree at most `n-2`; it must be zero, so
the row dependence is trivial. The free common-scaling tangent
`(q_i)_i` lies in the Jacobian kernel and has nonzero zeroth coordinate.
Fixing `q_0=1` therefore removes one dimension transversely.

This is a genuine one-dimensional algebraic reduction, unlike merely
fixing the pair imaginary residues. Its coefficients have the height
of the central integer vector computed above. On the eight-row fair
profile that logarithmic height is `38w+o(w)`, not `o(w)`. Thus it
does not meet a small-coefficient curve-extraction hypothesis.

Equivalently, the rational function

```text
F(z)=product_i (z-q_i)^(u_i),
P(z)=product_i (z-q_i)
```

satisfies

```text
F'(z)/F(z)=c z^(h-1)/P(z),             c!=0.            (10)
```

The nonnegative moment vanishings give a numerator of degree at most
`h-1`, and the negative moments give a zero of order at least `h-1`
at zero. Distinct poles and nonzero `u_i` make the numerator nonzero.
Its degree as a rational map depends on the sizes of the integers
`u_i`; no bounded-degree claim about `F` follows from a fixed number
of labels. This identity supplies no additional height theorem for
the moving curves (9).

## 6. Scope compared with earlier moment arguments

Equation (7) retains weighted integer moments, but its coefficients
grow at the heights in (2). It neither supplies a product equality nor
bounded multiplicities for
[angular moment product rigidity](angular_moment_product_rigidity.md).
Likewise, the
[conic feature audit](affine_conic_feature_content_audit.md) explains
why vanishing higher feature determinants must retain their primitive
contents. Formula (2) performs that normalization for these particular
spectral columns; it does not introduce an independent lower bound on
their heights.

No estimate here proves `N>=c(min X_i/L)^4`. A useful continuation
would need additional restrictions on the actual large integer
moment vector (7), or on several columns together, beyond their
norms, conjugacy, and already audited determinant.

The subsequent [six-point elliptic construction](six_point_fixed_moment_elliptic_obstruction.md)
shows why the six-label coefficient comparison alone is insufficient:
the primitive central vector can even stay fixed while all angles
coalesce. Its product collisions force compensating least-radius growth,
and its unbalanced cut mass excludes the full fair profile. Thus that
construction does not settle the six-label full-fair case either.

The expanded
[spectral checker](check_integer_cotangent_spectral.py) verifies (2)
on twelve actual cliques by recovering each `K_m` from Gaussian
subset-product lcms and gcds, without prime factorization. It also
checks every resulting bounded integer defect and the central real
moment identities or conjugate-column relation, as appropriate.
