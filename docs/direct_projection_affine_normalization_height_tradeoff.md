# Affine recentering cannot shrink the direct projection roots on a large short-arc cohort

Affine recentering of the direct monic square model can make its root
**differences** look small.  The parameter remains an integer only up to
an exact effective gcd.  A rational parameter can be cleared, but that
clearing restores the same effective root scale.  On an actual primitive
short endpoint arc, every large set of distinct projection roots has a
radial diameter of order at least `N^(1/2−epsilon)`, while the effective
gcd is `N^o(1)` after choosing a small-content actual anchor.  This
excludes the small-root premise for this **direct affine monic scheme**.
It does not exclude other identities or prove the general uniform
point-count goal.

Use [the direct projection model](direct_projection_monic_height_bridge.md).
Let `z_0,...,z_(m−1)` be a Gaussian-primitive physical tuple on squared
radius `N>1`.  Fix an actual anchor `z_0` and select `q>=2` other rows
`I`, with distinct radial deficits.  For every `i∈I` choose primitive
ordinary integers `a_i,b_i`, `b_i≠0`, such that

```text
z_i/z_0=(a_i+i b_i)/(a_i−i b_i),
f_i=a_i²+b_i² | 2N,       d_i=2N/f_i,
X_i=Re(conj(z_0)z_i),    c_i=N−X_i=d_i b_i²>0,
2N−c_i=d_i a_i².
```

The antipodal row, if present, has `(a_i,b_i)=(0,1)`; with distinct
rows there is at most one, so `gcd_I |a_i|` is positive.

## 1. Exact affine gcds, including residue cancellations

Choose one reference row `r∈I` and put

```text
G_aff=gcd_(i∈I)|c_i−c_r|,
V=2N−c_r,
G_eff=gcd(G_aff,V).
```

The gcd includes the zero difference at `i=r`.  Since
`V−(c_i−c_r)=2N−c_i`, the maximal **integral** affine divisor is

\[
 \boxed{G_{\rm eff}=\gcd_{i\in I}(2N-c_i)
       =\gcd_{i\in I}(d_i a_i^2).}                    \tag{1}
\]

It is independent of the chosen reference.

More generally, let `x -> alpha x+beta` be any nonconstant rational
affine change that makes both `2N` and every selected `c_i` integers,
with `alpha=u/v` in lowest terms and `v>0`.  Subtracting transformed
roots from the transformed parameter gives
`(u/v)(2N−c_i) in Z` for every `i`, so `v|G_eff`.
Thus `|alpha|>=1/G_eff` and the transformed root diameter is at
least `diam(c_I)/G_eff`, regardless of the translation or integer
numerator.  For an even square product, its transformed monic value
is the old square multiplied by `alpha^(2k)`, an integer rational
square and therefore an integer square.  Equation (1) is the exact
denominator limit for **every** rational affine monic normalization
of this selected set.

There is an exact primewise formula for `G_aff` itself, which shows
why allocation widths alone do not determine it.  At a rational prime
`p`, set

```text
alpha_i=v_p(d_i)=v_p(2N)−v_p(f_i),
s_i=alpha_i+2v_p(b_i)=v_p(c_i),
u_i=c_i/p^(s_i),            a p-unit.
```

For `i≠r`,

\[
 v_p(c_i-c_r)=
 \begin{cases}
  \min(s_i,s_r),&s_i\ne s_r,\\
  s_r+v_p(u_i-u_r),&s_i=s_r.
 \end{cases}                                         \tag{2}
\]

Taking the minimum over `i≠r` gives `v_p G_aff`.  In (2), the first
case is forced by the cut and primitive imaginary-residue valuations;
the second is an exact unit congruence and can produce extra primes or
powers.  At an odd split prime of `N`,
`alpha_i=e_p−|e_i−e_0|`; at two it is `1−eta_i`; outside `2N` it is
zero.  Thus (2) retains split, inert, ramified, and residue-only
primes without pretending that a high affine gcd is a cut width.

The effective gcd has a much simpler **factorization**.  Write

\[
 D_I=\gcd_{i\in I}d_i,\qquad
 A_I=\gcd_{i\in I}|a_i|.
\]

Then

\[
 \boxed{G_{\rm eff}=D_I A_I^2.}                    \tag{3}
\]

For the primewise proof, put `E=v_p(2N)`, `F_i=v_p(f_i)`, and
`alpha_i=E−F_i`.  Ordinary primitivity gives
`gcd(a_i,f_i)=1`.  If some `F_i>0`, a row maximizing `F_i`
minimizes `alpha_i` and has `v_p(a_i)=0`; the minimum of all
`v_p(a_i)` is also zero.  If every `F_i=0`, all `alpha_i=E`.
In either case

```text
min_i(alpha_i+2v_p(a_i))
  =min_i alpha_i+2 min_i v_p(a_i),
```

including `p=2`.  This proves (3) with all prime powers and actual
half-angle parity intact.

If one instead divides root differences by the larger `G_aff`, put
`g=G_eff`, `a=V/g`, `b=G_aff/g`, so `gcd(a,b)=1`, and
`u_i=(c_i−c_r)/G_aff`.  The apparent parameter is `a/b`.
For any even subset whose original monic product is an integer square,

```text
b^(2k) product (a/b−u_i)
 =product (a−b u_i)
 =product ((2N−c_i)/g),
```

which is an integer rational square, hence an integer square.  Its
integer roots are `b u_i=(c_i−c_r)/g`, **exactly** those produced by
integral normalization with `G_eff`.  Clearing the rational denominator
therefore undoes any extra shrinkage from `G_aff` in this monic Runge
workflow.  No unproved divisibility of an individual square root by
`g^k` is used.

An actual arithmetic example is `N=13`, `z_0=−3−2i`, and selected
rows `3−2i, 2−3i, −3+2i`.  Their distinct deficits are `18,13,8`.
Here `G_aff=5` but `V=8`, so `G_eff=1`: the rational roots
`0,−1,−2` acquire denominator five, and clearing it restores
integer roots `0,−5,−10`.  This is a gcd illustration, not a short
endpoint family.

## 2. The effective divisor is small on a large actual endpoint set

Assume now that the full primitive tuple lies in an arc of length at
most `C N^(1/4)`, with `0<C<=sqrt(2)`.  All Gaussian units remain
allowed.  For `S={0}∪I`, let `gamma_S` be the exact threshold-layer
weight constant on every row of `S`.  The allocation formula for
`d_i` gives, at an odd split prime `p|N`,

```text
v_p(D_I)=e_p−max_(i∈I)|e_i−e_0|,
v_p(g_0)=min(e_0,e_p−e_0),
gamma_(S,p)/log p=e_p−(max_S e−min_S e),
```

where `g_0=gcd_Z(Re z_0,Im z_0)`.  The difference between the first
and third exponents is at most the second: it is the smaller of the
largest excursions to opposite sides of `e_0`.  At two, `D_I` costs at
most one extra factor two.  Hence

\[
 D_I\le 2g_0\exp(\gamma_S).                       \tag{4}
\]

For `C<=sqrt(2)`, the pair allocation-sign Gram is obtuse even across
different literal Gaussian units.  The usual positive-semidefinite
sign-moment inequality bounds constant layers by

\[
 \gamma_S\le W/(q+1),\qquad W=\log N.             \tag{5}
\]

There is also a direct upper bound for the common half-angle real
content `A_I`.  Put `h_i=a_i+i b_i` and

```text
H_ij=h_j bar(h_i)/Norm(gcd_G(h_i,h_j)),
t_ij=Im H_ij=(a_i b_j−a_j b_i)/Norm(gcd_G(h_i,h_j)).
```

This is an exact Gaussian integer representing `z_j/z_i` as
`H_ij/bar(H_ij)`; its norm is `epsilon_ij n_ij`, with
`epsilon_ij∈{1,2}` the actual pair half-angle parity.  Since
`A_I|a_i,a_j` and `gcd(a_i,b_i)=1`, every rational prime of `A_I`
is absent from `f_i,f_j` and hence from the Gaussian gcd norm.
If `2|A_I`, all `a_i` are even, all `b_i` odd, and all `f_i` odd,
so there is no ramified loss either.  Therefore

\[
 A_I\mid t_{ij}\ne0\qquad(i,j\in I,\ i\ne j).     \tag{6}
\]

The literal chord bound with `C<=sqrt(2)` gives
`log|t_ij|<=−W mu_ij/4`, where
`mu_ij=1−2log n_ij/W<=0`.  Positivity of the sign Gram gives
`sum_(i<j)(−mu_ij)<=q/2`.  Multiplying (6) over the
`binom(q,2)` pairs yields

\[
 \boxed{A_I\le N^{1/[4(q-1)]}.}                 \tag{7}
\]

Combining (3)--(5) and (7),

\[
 \boxed{G_{\rm eff}\le
  2g_0 N^{1/(q+1)+1/[2(q-1)]}.}                    \tag{8}
\]

This inherited-height estimate uses the selected rows and their
actual prime allocations.  It does not remove their newly common
Gaussian factor or reset `N`.

The same interval calculation gives a complementary **lower** bound
for the ordinary lcm of the labels.  Put `L_I=lcm_(i∈I)d_i`.
At a split prime,

```text
v_p(L_I)=e_p−min_(i∈I)|e_i−e_0|,
gamma_(I,p)/log p=e_p−(max_I e−min_I e).
```

The selected allocations may leave a discrete hole around `e_0`;
the nearest selected allocation need not lie at an interval endpoint.
Nevertheless choosing `l=min_I e` when `e_0<=e_p/2`, or
`u=max_I e` when `e_0>=e_p/2`, proves the exact interval inequality:
in the first case `|l−e_0|<=l+e_0<=e_0+gamma_(I,p)/log p`;
in the second case
`|u−e_0|<=(e_p−u)+(e_p−e_0)<=(e_p−e_0)+gamma_(I,p)/log p`.
Hence

```text
min_(i∈I)|e_i−e_0|
 <= min(e_0,e_p−e_0)+e_p−(max_I e−min_I e).
```

Summing this primewise inequality, with no negative contribution at
two, gives

\[
 \boxed{L_I\ge\frac{N}{g_0\exp(\gamma_I)}
             \ge\frac{N^{1-1/q}}{g_0}.}            \tag{8a}
\]

The second inequality uses `gamma_I<=W/q` from the same obtuse
sign-moment argument for `I`.  With the small-content anchor below,
`log L_I>=(1−1/(2sqrt(m))−1/q)W`.  This is an ordinary lcm height
statement; it does not by itself lower-bound the rank of the labels'
squareclasses.

The total ordinary anchor-content estimate in
[endpoint allocation rounding](endpoint_allocation_rounding.md) is
`sum_j log gcd_Z(Re z_j,Im z_j)<=W sqrt(m)/2` in this same
`C<=sqrt(2)` range.  Choose an **actual** anchor with

\[
 g_0\le N^{1/(2\sqrt m)}.                         \tag{9}
\]

At most two physical rows share one radial deficit, so one can choose
`I` with distinct deficits and `q>=ceil((m−1)/2)`; thus `q` is
proportional to the full row count.  Equation (8) then saves only
`N^{O(1/sqrt(m))}` for this small-content choice of actual anchor.

## 3. Recentring cannot remove the actual root diameter

The Gram bound from
[the direct projection note, Section 4](direct_projection_monic_height_bridge.md)
applies at every actual anchor and every selected row: for any
`epsilon>0`, fewer than `B=1/(4epsilon²)` rows have
`n_ij>N^(1/2+epsilon)` for any fixed row `i`.  Suppose

\[
 q>3B+4.                                         \tag{10}
\]

Discard the fewer than `B` selected rows whose anchor pair norm is
bad.  On one angular side of the anchor, more than `B+2` good rows
remain.  Choose one of them; by the same rowwise Gram bound it has
a second good row on that side.  Order these two as
`0<theta_i<theta_j<pi/2`, reflecting angles if necessary.  The
source arc has angular span at most `C N^(−1/4)<pi/2` for primitive
`N>1` and `C<=sqrt(2)`.  The exact half-angle and root identities
give

```text
sin(theta_i/2)>=1/sqrt(2 n_0i),
sin((theta_j−theta_i)/2)>=1/sqrt(2 n_ij),
|c_j−c_i|=2N sin((theta_i+theta_j)/2)
                sin((theta_j−theta_i)/2)
          >=N/sqrt(n_0i n_ij)>=N^(1/2−epsilon).
```

Therefore the actual selected root diameter obeys

\[
 \operatorname{diam}(c_I)\ge N^{1/2-\epsilon}.   \tag{11}
\]

The maximal integer-affine normalization gives roots
`(c_i−c_r)/G_eff`; their diameter is at least (11) divided by
(8).  With the small-content actual anchor (9),

\[
 \boxed{\operatorname{diam}\left(
    \frac{c_i-c_r}{G_{\rm eff}}:i\in I\right)
  \ge\frac12 N^{1/2-\epsilon-1/(2\sqrt m)
                  -1/(q+1)-1/[2(q-1)]}.}         \tag{12}
\]

For a fixed small `epsilon` and `q` proportional to large `m`, this
is a fixed positive power of `N`; it is far above `N^{O(1/m)}`.
The rational-affine calculation in Section 1 shows that a larger
`G_aff` does not change this integer Runge root scale.  These are
conditions on a common affine normalization across **many** actual
roots.  A short subset may have a different gcd, and a new nonlinear
or multivariate relation is not excluded.

The [checker](check_direct_projection_affine_normalization_height_tradeoff.py)
tests (1)--(3), the exact relative Gaussian residues, and rational
denominator clearing on 48,688 literal physical subsets, including the
`G_aff>G_eff` example.  It checks (4) and (8a) primewise on the
small subsets and (7)--(8a) on an explicit short-arc Pell tuple.
The large-`q` diameter proof is the elementary Gram and sine argument
above; no synthetic large endpoint family is claimed.
