# Universal circle coordinates and the exact-content obstruction

This note combines the fixed-factor content criterion with the universal
cotangent family. There is no fixed-family extraction gap for these
coordinates: every fixed-size rational circle tuple occurs. What fails
is the proposed certificate. No nonzero homogeneous polynomial can
simultaneously satisfy its exact common-factor condition and the needed
endpoint chord inequality. A separate
[mixed-ideal argument](mixed_cotangent_content_ideal.md) extends this
obstruction to certificates using the common factor and its conjugate.
Neither argument rules out a bound on the absolute value of that factor.

## 1. A universal homogeneous family

Fix `k>=1`, put `H_i=X_i+iL`, and define

```text
W_0=product_i conjugate(H_i),
W_i=H_i product_(j!=i) conjugate(H_j),       1<=i<=k,
rho=product_i |H_i|.
```

The rows have degree `k`, and each has modulus `rho` at real parameters.
Their ratios are

```text
W_i/W_0=(X_i+iL)/(X_i-iL).
```

Given distinct rational Gaussian points `z_0,...,z_k` of equal nonzero
modulus, put `q_i=z_i/z_0`. Each `q_i` is a rational unit-circle point
other than one. The rational number `i(q_i+1)/(q_i-1)` is its half-angle
cotangent. Taking these as `X_i/L` recovers the tuple up to a common
rational Gaussian multiplier. Clearing rational parameters to integers
and dividing all `W_i` by their actual Gaussian gcd `G` therefore
recovers the least primitive integral circle tuple up to a Gaussian unit.
Its radius is exactly

```text
R=rho/|G|.                                           (1)
```

Every rational parameter tuple with distinct `X_i` can moreover be
scaled to an integer cotangent clique. Indeed the finitely many numbers

```text
X_i, L, (X_i X_j+L^2)/(X_i-X_j)
```

are rational, and all scale linearly under common parameter scaling.
A common denominator makes all of them integers. Neither the primitive
circle nor the homogeneous norm ratios change. Thus imposing integer
pair cotangents by this common scaling does not restrict the rational
projective parameter domain.

For an arc anchored at one endpoint with actual angular differences
`0<theta_i<=pi/2`, the ratios are `X_i/L=cot(theta_i/2)>=1`. The
corresponding parameter cone is `X_i>=L>0`.

## 2. The universal row ideal has no extra degree-balanced sections

Set `Y_i=X_i-iL`. Since

```text
W_i-W_0=2iL product_(j!=i)Y_j,
```

the row ideal over `Q(i)` is

```text
I=(product_i Y_i, L product_(j!=i)Y_j : 1<=i<=k).
```

It consists of all squarefree products omitting one of the `k+1`
variables `L,Y_1,...,Y_k`. The elementary proof in
[universal_cotangent_monomial_ideal.md](universal_cotangent_monomial_ideal.md)
shows, for every integer `q>=1`, that

```text
overline(I^q)=I^q.
```

More precisely, a homogeneous polynomial `P` of degree `qk` satisfies
the chartwise integral-closure condition for `I^q` if and only if each
monomial in its `L,Y` expansion has exponent at most `q` in every
variable. In particular,

```text
P!=0 and chartwise IC  =>  ord_(L=0)(P)<=q.           (2)
```

The invertible substitution `Y_i=X_i-iL` preserves divisibility by
powers of `L`, so the same order bound holds in the original variables.

The fixed-factor criterion applies to the products of `q` rows: their
common Gaussian gcd is `G^q`, and their common degree is `qk`. Therefore
chartwise IC is equivalent to the existence of a fixed nonzero Gaussian
integer `C` with

```text
G(x)^q divides C P(x)                                (3)
```

for every ordinary integer specialization at which all relevant values
are nonzero. See
[fixed_factor_primitive_content_criterion.md](fixed_factor_primitive_content_criterion.md).
The criterion concerns exact divisibility, not just absolute values.

It is enough to impose (3) on integer cotangent cliques in the positive
open cone `X_i>L>0` with distinct `X_i`. To see the necessity, take the
unbounded local divisibility defect supplied by the fixed-factor proof.
Approximate its integer vector to sufficiently high precision at the
fixed split prime, then add multiples of that prime power to the
coordinates to enter the positive cone and avoid equality hyperplanes.
Such vectors exist because any open real cone contains arbitrarily
large boxes. Polynomial continuity preserves the local defect. Remove
any common integer factor; the defect is unchanged because both sides
have degree `qk`. Finally scale to clear all pair-cotangent denominators,
which again leaves the divisibility defect unchanged. Thus cone or
integrality restrictions do not evade (2) within this homogeneous class.

## 3. No balanced chord certificate in any power

Fix `t>=1` selected pairs of distinct rows; pairs may be repeated.
Suppose a nonzero homogeneous degree-`qk` polynomial `P` were to satisfy
both (3) and the endpoint certificate

```text
product_(ell=1)^t |W_(a_ell)-W_(b_ell)|
    >= c |P|^(t/(2q)) rho^(t/2),        c>0,          (4)
```

throughout that cone. Dividing by `(rho |G|)^(t/2)` would imply that
some primitive chord divided by `sqrt(R)` is at least

```text
c^(1/t) / |C|^(1/(2q)),
```

which is precisely the desired uniform separation for this tuple size.

However, the exact difference formulas are

```text
W_i-W_0=2iL product_(j!=i) conjugate(H_j),
W_i-W_j=2iL (X_j-X_i) product_(ell!=i,j) conjugate(H_ell).
```

Fix distinct positive rational `X_i` generically and let rational
`L -> 0+`. Every selected chord vanishes to exactly first order in `L`,
while `rho` tends to a positive number. If `v=ord_(L=0)(P)`, choose the
`X_i` so that the leading coefficient of `P` is nonzero. Inequality (4)
would force

```text
L^t / L^(v t/(2q)) cannot tend to zero,
so t-v t/(2q)<=0, equivalently v>=2q.                 (5)
```

This contradicts `v<=q` from (2). Hence no such nonzero `P` exists.
The same argument rules out a finite menu of these certificates, even
with different powers and chord selections: choose positive rational
`X_i` avoiding the finitely many leading-coefficient zero sets. Every
member fails for all sufficiently small `L`.

The rational limiting points may all be cleared to integer cotangent
cliques as in section 1. Homogeneity preserves the inequalities, so this
is a contradiction on the integer domain as well.

## 4. Proper birational changes do not repair this certificate class

Let `pi:B -> P^k` be a proper birational morphism with `B` normal.
Let the transformed rows `m pi^*W_j` be global sections of one line
bundle `M`, where `m` is a nonzero rational section of
`M tensor pi^*O(-k)`. Require the transformed certificate `P'` to be
a section of `M^q`; this is the degree-balancing condition on the new
model. Suppose that, in local trivializations everywhere on `B`,
`P'` is integral over the `q`-th power of the transformed row ideal.
Dividing a monic dependence equation by the appropriate power of `m`
shows that the rational section of `pi^*O(qk)`

```text
f=P'/m^q
```

is integral over the `q`-th power of the original pulled-back ideal.
Normality of `B` makes `f` a regular section of `pi^*O(qk)`.

This section descends to a degree-`qk` section `P` on `P^k`. One can see
the descent directly on each factorial affine chart: properness lifts
the spectrum of every codimension-one valuation ring to `B`, so the
rational function represented by `f` has no pole at any irreducible
polynomial downstairs. A rational function with that property is a
polynomial. The chart sections glue with the `O(qk)` transition rules.
The same proper lifting, now for every valuation ring, descends the
integral-closure inequalities, so `P` satisfies the original chartwise
IC condition. The valuative criterion is the same one used in the
fixed-factor note.

Every chord and the common row modulus acquire a factor `|m|`.
The two sides of (4) acquire exactly `|m|^t`, which cancels. Thus the
transformed certificate would give the impossible original one.
This statement requires a global certificate on a proper model.
An estimate on a single affine chart with omitted boundary denominators
is not covered by it.

## 5. A remaining arithmetic norm target

For positive `X_i`, put

```text
P_ij=L^2 product_(ell!=i,j) X_ell,        i!=j.
```

The universal chord formulas give

```text
|W_i-W_0| |W_j-W_0|
   =4 L^2 rho^2/(|H_i| |H_j|)
   >=4 P_ij rho.                                    (6)
```

Consequently an absolute-value bound

```text
|G| <= C P_ij                                       (7)
```

on the actual integer cliques, for a fixed sufficiently large `k` and
an appropriately selected pair, would force normalized arc span at
least `2/sqrt(C)`. Taking the two smallest `X_i` is one definite target.
This is sufficient for a uniform count after partitioning longer arcs.
For the same fixed number of labels it asks for a stronger norm
comparison than the exact lcm target (H), which only uses the smallest
cotangent. But existence of such a bound for some fixed number of
labels is equivalent to (H), and hence to uniformity.

Indeed, suppose (H) holds for `k_0` finite cotangents with constant `c`.
Take `k=k_0+1` finite cotangents ordered as
`L<=X_1=A<X_2<...<X_k`. Delete `X_1`, retaining the anchor. Its all-edge
lcm divides that of the full clique, so (H) applied to the remaining
`k_0` finite cotangents gives

```text
R_full >= R_subset >= sqrt(c) (X_2/L)^2.
```

Since `rho<=2^(k/2) product_i X_i`, formula (1) now yields

```text
|G| <= [2^(k/2)/sqrt(c)] L^2 product_i X_i/X_2^2
    <= [2^(k/2)/sqrt(c)] L^2 product_(i>=3) X_i,
```

which is (7) for the two smallest cotangents. Conversely, (7) for
that pair implies, directly from (1),

```text
R >= rho/(C P_12) >= X_1 X_2/(C L^2)
  >= (A/L)^2/C.
```

Squaring gives (H) with `c=C^(-2)` for the same `k`. Thus the
absolute-content route is a second exact formulation of the goal,
with at most one extra finite cotangent in the implication from (H).

For `k=3`, the Pell four-point family disproves (7): its normalized
span tends to zero, contradicting the consequence of (6)--(7).
See [integer_cotangent_lcm_height_target.md](integer_cotangent_lcm_height_target.md).
For a larger fixed `k`, (7) remains unproved.

Although `P_ij` cannot satisfy the exact divisibility criterion, that
fact does not disprove (7). The fixed-factor note contains an explicit
example where a uniform norm bound holds while every fixed exact
Gaussian divisibility bound fails. The argument above treats `G^q`;
the [mixed-ideal theorem](mixed_cotangent_content_ideal.md) separately
proves that a degree-`k(a+b)` chartwise certificate for
`G^a conjugate(G)^b` has contact order strictly less than `2(a+b)`.
Thus that exact mixed divisibility class also misses the chord order.
The absolute-value target remains open, and the general bound on arc
counts is unchanged.
