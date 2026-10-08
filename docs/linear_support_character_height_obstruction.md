# A growth bound for repeated Walsh endpoint profiles

For repeated Walsh sign profiles with `b>=5` copies, canonical or
relabelled single-entry flips, and nearly equal actual prime logarithms,
any actual endpoint realization satisfies
`M=O_C((log N/loglog N)^(2/3))`. Higher-dimensional affine-plane
certificates force `b` to grow at least proportionally to `sqrt(M)`;
the number of distinct prime factors then supplies the growth bound.
The theorem permits arbitrary row units and common Gaussian content.
It is restricted to this profile family, not arbitrary endpoint tuples.

The later [weighted parity argument](walsh_ordinary_dot_weighted_phase_growth.md)
removes the nearly equal prime-log hypothesis from the count estimate
for both assignments. The copy-count conclusion here still uses that
hypothesis.

At five copies, an explicit **four-row integral character** already
excludes the canonical assignment for `M>=32`, and a natural relabelled
assignment for `M>=16`, once the prime scale is large. The actual
normalized arc constant grows at least like `M^3` in the archived
nearby-prime family. Both results use a fixed `pi/4` target and an
elementary integer-coordinate gap, with no transcendence estimate.

There is a complementary exact obstruction to the proposed generic
short-character shortcut. After grouping whole cut classes, these
matrices have only `2M-1` classes, but their saturated integer character
lattice has minimum at least `M/2`. With the actual nearby-prime weights,
every nonzero character has logarithmic Gaussian height above `3W_0/10`,
where `W_0` is the total physical-prime log weight before common content.
This blocks a sublinear-character or generic `W/4`-height deduction from
obtuse geometry alone. It does **not** block the sharper `W/2` threshold
for a `t=1` rational-grid certificate, which gives the positive result.

The construction and its primitive Gaussian realization are from
[strict_obtuse_prime_box_countermodels.md](strict_obtuse_prime_box_countermodels.md).
Its arguments were not asserted to form endpoint arcs. The exclusion
here supplies a specific missing phase obstruction for its canonical
assignment; it does not contradict that note or cover arbitrary flip
assignments. The full-fair exponential-support minimum already appears
in [the pattern-lattice note](corank_one_local_bessel_pair_filter.md).
The uniform endpoint count remains unproved.

## Exact construction and saturation

Let `M=2^t>=8`, with rows indexed by `x in F_2^t`. Write
`h_a(x)=(-1)^(a dot x)` for the `M-1` nonconstant Walsh columns.
Take five physical copies of each `h_a`. At row `x`, flip the entry in
one assigned physical column; all `M` assigned columns are distinct,
and each Walsh class receives at most two assignments. For example,
list the physical columns in five successive copies of `a=1,...,M-1`
and assign row `x` to physical column `x` (zero-based).

Each class retains three, four or five unflipped copies. After grouping
identical whole classes, the sign matrix is

```text
S=[H_0,F],    F_x=h_(a_x)-2 h_(a_x)(x)e_x,
T=(M-1)+M=2M-1.                                      (1)
```

These columns represent distinct nonconstant unoriented cuts. Different
Walsh columns have Hamming distance `M/2` even after one is negated;
two flips cannot erase that distance when `M>=8`. Columns arising from
the same Walsh column have different flipped positions.

For a rational row vector `lambda` with `sum lambda=0`, put
`u=lambda^t H_0` and `v=lambda^t S`. Then

```text
v_(F_x)=u_(a_x)-2 lambda_x h_(a_x)(x).                 (2)
```

If `v` is integral, (2) forces `c=2lambda` to be integral. Conversely,
if `c` is integral with zero sum, then
`u_a=(1/2)sum_x c_x h_a(x)` is integral: subtracting the zero sum
expresses it as a signed integer subset sum. Equation (2) makes the
remaining coordinates integral as well. Hence the saturated lattice is
exactly

```text
{lambda^t S: lambda rational, sum lambda=0} intersect Z^T
   = {(c^t S)/2: c in Z^M, sum c=0}.                  (3)
```

Walsh inversion on the zero-mean subspace gives

```text
c_x=(2/M)sum_(a!=0) u_a h_a(x).
```

For nonzero `c`, its maximum absolute coordinate is at least one.
Consequently every nonzero lattice vector satisfies

```text
||v||_1 >= ||u||_1 >= (M/2)||c||_infinity >= M/2.     (4)
```

This is a saturated-lattice statement, not merely a bound for integral
row coefficients. Its order is sharp: take `c=e_i-e_j`. Exactly `M/2`
old coordinates are nonzero, each of absolute value one. Since each old
class receives at most two assignments, the uncorrected new coordinates
have total absolute value at most `M`. The two corrections in (2) add
at most two. Thus

```text
M/2 <= min_(0!=v) ||v||_1 <= 3M/2+2.                 (5)
```

For the canonical and relabelled assignments below, each old class is
assigned once and one class twice. Their upper bound improves to `M+3`.

## Strict obtuseness and actual Gaussian heights

Before the flips, the physical-column row Gram has off-diagonal entry
`-5`. Each of the two affected entries for a row pair changes this by
at most two, so after the flips every off-diagonal entry is at most
`-1`. Grouping columns retains this Gram by assigning each old class
its remaining multiplicity three, four or five and each new class weight one.

Use the actual nearby distinct split primes supplied in the earlier
countermodel, one per physical column. For `r=5(M-1)` they can satisfy

```text
ell <= log p_j < ell+1/r,    ell>1,
W=sum_j log p_j < r ell+1,  W=O(M log M).              (6)
```

The weighted off-diagonal Gram is at most `-ell+1<0`. The literal
Gaussian realization is primitive, has width one at every prime, and
has distinct pair prime-norm supports, as proved in the earlier note.
No prescribed real prime weights are substituted for the actual logs.
The archived checker uses the opposite global sign convention for the
Walsh columns. This conjugates all literal Gaussian source rows and the
composite directions; their norms and containing-arc lengths are unchanged.

Each whole cut block `Gamma_j` uses disjoint rational primes, and each
old block contains at least three. For an integral character `v`, form
`Beta_v` by multiplying `Gamma_j^v_j` for positive entries and
`bar(Gamma_j)^(-v_j)` for negative entries. There is no common conjugate
content to remove: each rational prime belongs to only one block and
appears in only one orientation in this product. Thus the exact height is

```text
V(Beta_v)=sum_j |v_j| log Norm(Gamma_j)
          >=3 ell ||u||_1 >=(3M/2)ell >3W/10.         (7)
```

The last strict inequality follows from (6) and `ell>1/5`.
The direction is nonunit for every nonzero `v`, without an additional
endpoint assumption.

The certificate itself can have tiny row coefficients: for
`c=e_i-e_j`, equation (3) gives `c^t S=2v`, with `||c||_1=2`.
The certificate has `t=1`, so its composite phase has a rational
Gaussian-unit half-angle target, a multiple of `pi/4`. The elementary
integer-coordinate gap at this grid requires height only about `W/2`,
which is stronger than the general Roth threshold. Accordingly (7)
stays above the asymptotic `W/4` scale of the generic fixed-target
argument but does not prevent a useful `t=1` certificate below `W/2`. It also prevents using the Bessel-plus-
Matveev sufficient condition `||v||_1 log B=o(M)` from the
[effective short-character theorem](effective_hadamard_short_character_phase_gap.md).
This does not rule out a sharper arithmetic estimate for some special
composite, simultaneous phase information, or another construction.

The [checker](check_linear_support_character_height_obstruction.py)
verifies the Walsh identities, distinct cut classes, strict Gram,
saturation equivalence on exact rational fixtures, the character bounds,
and a literal prime/Gaussian height fixture. Its checks do not assert
that the configurations lie on endpoint arcs.


## An explicit four-row phase exclusion

Here are two assignments with an explicit certificate. In both, row
`x` flips a distinct physical column of class `a_x`, with no class
receiving more than two flips.

* **Canonical assignment:** `a_x=1+(x mod (M-1))`, as in the archived
  checker. For `M>=32`, take `x0=2`, `p=6`, `q=24`.
* **Relabelled assignment:** `a_x=x` for `x!=0`, and `a_0=1`, using a
  second physical copy at row zero. For `M>=16`, take `x0=1`, `p=3`,
  `q=12`.

Interpret indices as binary vectors. In each case the plane
`V=span(p,q)` is totally isotropic for the binary dot product:
`p dot p=q dot q=p dot q=0`. The restriction
`l(v)=x0 dot v` on `V` is nonzero. On the four rows `x=x0+v`, put
`c_x=(-1)^l(v)`, and put `c=0` elsewhere. Thus `sum c=0` and
`||c||_1=4`. The canonical rows are `[2,4,26,28]`, the relabelled rows
`[1,2,13,14]`, and in these orders the coefficients are `[1,-1,1,-1]`.

For the canonical rows, all indices are even and avoid the exceptional
last row, so `a_x=x+1=x XOR 1`; the extra bit is orthogonal to `V`.
For the relabelled rows, `a_x=x`. In either case every assigned class
has restriction `a_x|V=l`. For every old Walsh class,

```text
u_a=(1/2)sum_x c_x h_a(x)
   =2 h_a(x0) if a|V=l, and 0 otherwise.              (8)
```

Exactly `M/4` old classes qualify, all nonconstant. At each of the four
flipped rows, `c_x h_(a_x)(x)=h_(a_x)(x0)`. Thus (2) changes the
absolute value in that physical column from two to one. The other
flips do not change its character coordinate. The full physical-prime
character has exact absolute coefficient sum

```text
sum_(physical columns j)|v_j| = 5M/2-4,  |v_j|<=2.    (9)
```

For actual prime logs in (6), its conjugate-primitive composite obeys

```text
V(Beta)-W/2 < -(3/2)ell+2.                            (10)
```

Indeed its main equal-weight contribution is `(5M/2-4)ell`; the
weighted error is less than two, while `W/2>=5(M-1)ell/2`.

If the actual source rows lay on an arc of angular width
`Delta<=C exp(-W/4)`, arbitrary row units would give

```text
product_x z_x^c_x = u Beta/bar(Beta),  u in {1,-1,i,-i}.
```

The zero row sum cancels the arc's central direction, and
`|sum c_x theta_x|<=2Delta`. Therefore

```text
dist(arg Beta, (pi/4)Z) <= Delta.                      (11)
```

A conjugate-primitive nonunit Gaussian integer cannot lie on any of
these axis or diagonal lines. Its distance to such a line is at least
`1/sqrt(2)`, since the relevant integer linear form is nonzero. Hence
`dist(arg Beta,(pi/4)Z)>=1/sqrt(2 Norm(Beta))`. Combining with (11)
proves the necessary endpoint height bound

```text
V(Beta) >= W/2-log(2C^2).                              (12)
```

More directly, if `C_source` is the normalized constant of any containing
arc, these same inequalities give the quantitative anticoncentration
bound

```text
C_source > exp(3ell/4)/(sqrt(2) e).                    (13a)
```

Thus when `ell=4 log M+O(1)`, as in the archived nearby-prime family,
its actual normalized arc constant grows at least a fixed positive
multiple of `M^3`. No minor-arc assumption is needed: (11) uses lifted
arguments on any containing arc. If all source rows are multiplied by
a common nonzero Gaussian factor `d`, the lower bound (13a) acquires
the additional factor `Norm(d)^(1/4)`; the common factor cancels from
the signed product while remaining in the original radius.

Equations (10)--(12) in particular contradict one another whenever

```text
(3/2)ell > 2+log(2C^2).                               (13)
```

This explicit sufficient threshold applies to either assignment and
any fixed `C>0`. It is eventually satisfied by the nearby-prime family
of the archived construction. No claim is made for an arbitrary
injection of the flipped entries: the four-row alignment is a genuine
additional combinatorial condition. The gain uses a long pattern
character with a short row certificate and its exact rational phase
normalization, rather than producing a sublinear pattern character.

## Varying copy counts: a growth bound in the repeated-Walsh family

The same certificate works in higher dimension and yields a growth
estimate when the number of copies is allowed to increase. Fix `C>0`.
Take `M=2^t` rows and `b>=5` physical copies of every nonconstant Walsh
column, with one entry flipped per row according to either the canonical
or relabelled assignment above. Thus `r=b(M-1)` distinct physical split
primes occur, each to the first power. Assume their actual logarithms
satisfy

```text
ell <= log p_j < ell+1/r,       ell>=log 5.             (14)
```

Independent row units and a complete common Gaussian factor `d` are
allowed. Write `W_0=sum_j log p_j`, `D=log Norm(d)>=0`, and
`W=log N=W_0+D` for the actual squared radius. If the actual rows lie
on an arc of length `C N^(1/4)`, then within this profile family

```text
M = O_C((W/log W)^(2/3)).                              (15)
```

The claim is restricted to these assignments and the nearly equal
prime logs in (14). It neither constructs such endpoint families nor
bounds arbitrary endpoint tuples.

**General affine-coset certificate.** The height argument below applies
to any one-entry-per-row assignment, not just the two specified ones.
Let `P=x0+V` be an affine coset of an F_2-linear subspace, `h=|V|>=2`.
Suppose a fixed nonzero linear form `l` on `V` satisfies

```text
a_x restricted to V = l       for every x in P.       (16a)
```

The physical column flipped at each row is distinct; there are `b`
physical copies of every nonconstant Walsh class. No isotropy of `V`
is needed for this lemma. Define `c_(x0+v)=(-1)^l(v)` and `c=0` off
`P`. Then `sum c=0`, `||c||_1=h`, and the old Walsh character equals
`(h/2)h_a(x0)` on exactly the `M/h` nonconstant classes restricting to
`l`, and vanishes elsewhere. On each row of `P`, condition (16a) gives
`c_x h_(a_x)(x)=h_(a_x)(x0)`. Consequently all `h` flips on `P` lower
their coefficient's absolute value by one; the other flips have no
coefficient change. Thus

```text
sum_(physical j)|v_j| = bM/2-h,    |v_j|<=h/2.         (16)
```

The total prime-log error is less than `(bM/2-h)/r<1`, using the sum
in (16), rather than the maximum single coefficient. The composite
is conjugate-primitive and nonunit because its nonzero exponents are
on disjoint physical primes. Its exact height therefore satisfies

```text
V(Beta)-W/2 < -(h-b/2)ell+1-D/2.                      (17)
```

The exact row certificate has `t=1`:
`product z_x^c_x=u Beta/bar(Beta)`. The lifted arc estimate gives
`dist(arg Beta,(pi/4)Z)<=h Delta/4`. Applying the integer line gap
from (11)--(12) yields

```text
V(Beta) >= W/2+log 8-2log(h C).
```

Thus an actual endpoint realization admitting (16a) necessarily obeys

```text
(h-b/2)ell+D/2 < 1+2log(h C)-log 8.                   (18)
```

In particular `b>=2h-O_C(log h/ell)`. Equations (16)--(18) are the
general certificate; the remaining step finds a large aligned coset
for the two assignments in the theorem.

For the relabelled assignment, take

```text
dim V=floor(t/2),  V=span(3,12,48,...),  x0=1.
```

For the canonical assignment, take

```text
dim V=floor((t-1)/2),  V=span(6,24,96,...),  x0=2.
```

The listed binary vectors have disjoint pairs of nonzero bits and span
a totally isotropic subspace. The restriction `l(v)=x0 dot v` is
nonzero, and `x0+V` excludes zero. In the canonical case every row in
the coset is even, so `a_x=x XOR 1` there, with no exceptional last
row. Since bit zero is orthogonal to `V`, both assignments satisfy
(16a). This gives a lower bound on their required copy count rather
than just a fixed-copy exclusion.

Here is an explicit coarse threshold for converting it to (15). Put

```text
H_C=max(16, ceil(8 log^+(C)/log 5)).
```

For `h>=16`, `2log h<=(h/4)log 5`; check at sixteen and differentiate.
For `h>=H_C`, also `2log^+C<=(h/4)log 5`. Since `1-log 8<0`, (18)
therefore rules out `b<=h`. In either assignment `h>=sqrt(M)/2`.
Hence for `M>=4H_C^2`,

```text
b>h>=sqrt(M)/2,
r=b(M-1)>=M^(3/2)/4.                                  (19)
```

Order the `r` distinct prime norms. The j-th is at least `j+1`, so

```text
W>=W_0>=log((r+1)!) >= (r/2)log(r/2)   (r>=2).        (20)
```

Equations (19)--(20) give `W>=c M^(3/2)log M` for an absolute positive
`c` and all sufficiently large `M`. Inverting this increasing function
proves (15); bounded `M<4H_C^2` is absorbed in its C-dependent constant.
This improves the general logarithmic rank scale within the stated
varying-support profile family. Common source content remains in `W`
and strengthens (18); no primitive renormalization loses its radius.
