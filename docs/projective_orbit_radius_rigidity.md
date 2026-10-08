# Endpoint projective orbits determine the radius to relative logarithmic error O(1/M)

## Status

This note treats rational projective transformations of arbitrary height.
If two labelled configurations of `M` lattice points both lie in arcs
of length at most `2 sqrt(radius)` and are projectively equivalent in
their rational half-angle coordinates, their primitive radii satisfy

```text
log R_new / log R_old = 1+O(1/M).                (1)
```

The implied constant is absolute; explicit finite bounds appear below.
The same assertion holds for any two fixed normalized arc constants,
with a constant depending on those two constants. There is no bound
on the matrix height or determinant.

The new step is to average the exact cross-ratio height identity already
proved in Section 7 of
[uniformity_audit_parallel.md](uniformity_audit_parallel.md). The average
is an invariant of the labelled projective configuration. Pairwise
obtuseness controls its error by `O(W/M)`, improving the `O(W/sqrt(M))`
error supplied by selecting one jointly near-uniform quadruple.

Consequently a projective map preserving the endpoint condition cannot
lower the radius by any fixed positive power in a sequence with
`M -> infinity`. This is not a proof of uniformity. In particular,
it neither constructs a descent nor excludes reductions by a factor
`exp(O(log R/M))`.

## 1. Normalization and a common-unit cohort

First divide the common Gaussian gcd out of each entire configuration.
This decreases its normalized arc constant and leaves its cross-ratios
unchanged. The radius of each resulting primitive tuple is its least
radius as an integral realization of its relative phases. Write its
Gaussian factorization as

```text
z_i=epsilon_i product_p pi_p^(a_ip) bar(pi_p)^(e_p-a_ip),
N=R^2,       W=log N=sum_p e_p log p.
```

Only odd split primes remain in `N`. A prime inert in `Z[i]`, or the
ramified prime `1+i`, would otherwise divide every point.

For the moment fix a cohort of `L>=8` rows with the same literal
Gaussian unit `epsilon_i`. They remain on the inherited parent circle
of radius `R`. **Do not divide out a new cohort gcd.** Every height and
threshold weight below continues to refer to the parent `W`.

Suppose its parent arc has length at most `C sqrt(R)`, where `0<C<=2`.
For distinct rows define

```text
d_ij=sum_p |a_ip-a_jp| log p,
h_ij=product_p pi_p^max(a_jp-a_ip,0)
                  bar(pi_p)^max(a_ip-a_jp,0),
t_ij=Im(h_ij),       t_ij is a nonzero integer.
```

The exact chord formula is

```text
|z_j-z_i|=2R |t_ij|/|h_ij|,
log|h_ij|=d_ij/2.
```

For normalized threshold signs and their inner products put

```text
f_i(p,t)=2 indicator(a_ip>=t)-1,
weight(p,t)=log p/W,
mu_ij=E_layers f_i f_j=1-2d_ij/W.
```

The chord bound and integrality of `t_ij` imply

```text
mu_ij<=0,
0<=log|t_ij|<=-W mu_ij/4+log(C/2)
                         <=-W mu_ij/4.          (2)
```

Let a pair average always mean an average over ordered distinct rows.
The norm square of `sum_i f_i` is nonnegative. Thus

```text
0 <= E_layers (sum_i f_i)^2
   = L+sum_(i!=j) mu_ij <= L,
mean_pair(-mu_ij)<=1/(L-1),
mean_pair log|t_ij|<=W/[4(L-1)].                 (3)
```

This argument takes place in the inherited height. For example, if a
proportion `gamma` of the threshold measure is constant on the whole
cohort, the first line of (3) gives `gamma L^2<=L`, hence
`gamma<=1/L`. A newly common cohort factor is therefore controlled
without changing the normalization.

## 2. Exact conductor and residue parts of the cross-ratio

For four distinct labelled rows `0,1,2,3`, use

```text
lambda=(z_0-z_1)(z_2-z_3)/((z_0-z_2)(z_1-z_3)).
```

This is also the cross-ratio of their rational half-angle parameters:
the complex radial coordinate is a fractional-linear function of the
half-angle parameter. It is consequently rational and invariant under
`PGL_2(Q)` changes of those parameters. We use the projective cross-ratio,
not a quotient of the squared Euclidean chord lengths.

The exact primitive-chord factorization gives

```text
lambda=(A/B) (t_01 t_23)/(t_02 t_13),
A=product_(b_p>0) p^b_p,
B=product_(b_p<0) p^(-b_p),       gcd(A,B)=1,
b_p=min(a_0p,a_1p)+min(a_2p,a_3p)
       -min(a_0p,a_2p)-min(a_1p,a_3p).           (4)
```

Equivalently,

```text
b_p=(|a_0p-a_2p|+|a_1p-a_3p|
       -|a_0p-a_1p|-|a_2p-a_3p|)/2.
```

At a threshold, its contribution is `+1` on the cut `01|23`, `-1`
on the cut `02|13`, and zero on every other cut. At a fixed prime
the threshold subsets are nested, so there is at most one unordered
two-against-two cut. Positive and negative contributions cannot cancel
inside that prime. Therefore `log A` and `log B` are exactly the
total threshold weights of these respective cut types.

There is also no cancellation against the opposing residues. If
`b_p>0`, both pairs `02` and `13` have distinct allocations at `p`;
their primitive residues are `p`-units. If `b_p<0`, the same assertion
holds for `01` and `23`. Thus for the reduced rational form
`lambda=n/q`, `q>0`,

```text
A divides n,       B divides q.
```

Define the logarithmic rational height
`h(lambda)=log max(|n|,q)`. It satisfies the exact bounds

```text
max(log A,log B) <= h(lambda)
 <= max(log A+log|t_01 t_23|,
        log B+log|t_02 t_13|)
 <= max(log A,log B)+sum_four_pairs log|t_ij|.   (5)
```

These statements retain arbitrary prime powers and all residue signs.
Their earlier independent proof is in the audit note cited above.

## 3. Averaging the conductor height

All quadruple averages below range over ordered distinct rows of the
fixed cohort. At a threshold write

```text
S=sum_(i=1..L) f_i,
n_+=(L+S)/2,       n_-=(L-S)/2.
```

Four of the ordered sign patterns contribute to `log A+log B`.
Their total sampling probability is

```text
4 n_+(n_+-1)n_-(n_--1)/[L(L-1)(L-2)(L-3)]
 = (L^2-S^2)((L-2)^2-S^2)
     /[4L(L-1)(L-2)(L-3)].                     (6)
```

Since a prime has no cancellation between cut signs, averaging (6)
over layers gives the exact identity

```text
mean_quad(log A+log B)/W
 = E_layers[(L^2-S^2)((L-2)^2-S^2)]
     /[4L(L-1)(L-2)(L-3)].                     (7)
```

Set `v=E S^2` and `u=E S^4`. Equation (3) and elementary moment
bounds give

```text
0<=v<=L,       v^2<=u<=L^2 v.
```

The numerator of (7) is

```text
L^2(L-2)^2-[L^2+(L-2)^2]v+u.
```

The upper bound follows by using `u<=L^2 v` and then `v>=0`.
For the lower bound use `u>=v^2`; the resulting quadratic is decreasing
on `0<=v<=L`, so substitute `v=L`. After division by the denominator,

```text
1/4-1/[2(L-2)(L-3)]
 <= mean_quad(log A+log B)/W
 <= 1/4+(2L-3)/[4(L-1)(L-3)].                  (8)
```

To pass from the sum to the maximum, the imbalance has a stronger
pairwise bound. Equation (4) gives

```text
log A-log B=(d_02+d_13-d_01-d_23)/2
 = (W/4)(mu_01+mu_23-mu_02-mu_13).
```

Using `mu_ij<=0` and (3),

```text
mean_quad |log A-log B| <= W/(L-1).             (9)
```

Combine (8)--(9) with
`max(a,b)=(a+b+|a-b|)/2`. Finally, the four-pair residue sum in
(5) has mean at most `W/(L-1)` by (3). This proves

```text
a_L W <= I_L <= b_L W,
I_L=mean_quad h(lambda),
a_L=1/8-1/[4(L-2)(L-3)],
b_L=1/8+(2L-3)/[8(L-1)(L-3)]+3/[2(L-1)].       (10)
```

In particular `a_L>0` for `L>=8`, and

```text
I_L=W/8+O(W/L),
b_L/a_L=1+14/L+O(1/L^2).                       (11)
```

The lower error in (10) is actually `O(W/L^2)`. The larger upper
error charges the absolute imbalance and the small primitive residues.

## 4. Radius rigidity at every projective matrix height

Consider two labelled primitive endpoint configurations of `M` points,
with heights `W=log(R^2)` and `W'=log((R')^2)` and constants at most
two. Suppose their rational half-angle configurations are related by
an invertible rational projective map.

Partition the labels by the pair consisting of their input and output
Gaussian units. There are at most 16 such classes. Choose a common
cohort with `L>=M/16`, and suppose `L>=8`. Its units are constant
separately in both configurations. No new common divisor is removed
from either cohort.

The map preserves each labelled quadruple's `lambda`, so the same
average `I_L` appears in (10) for both inherited heights. Therefore

```text
a_L/b_L <= W'/W <= b_L/a_L.                     (12)
```

Equivalently,

```text
R^(a_L/b_L) <= R' <= R^(b_L/a_L),
log R'/log R=1+O(1/M).                         (13)
```

This proof has no matrix-height, determinant, content, or prime-support
restriction. Any such effects have already been included in the actual
integer configurations and in their exact invariant cross-ratios.

For fixed constants `C,C'` greater than two, partition each arc into
`ceil(C/2)` and `ceil(C'/2)` pieces of the required length. Intersect
these arc-piece labels and the two unit labels. The selected cohort
still has size at least

```text
M/[16 ceil(C/2) ceil(C'/2)].
```

Applying (12) on this inherited cohort proves the same asymptotic
statement, with a constant depending on `C,C'`.

## 5. Consequence for orbit minimizers and the remaining gap

Fix an `M`-point endpoint configuration and a bound on the normalized
arc constant. Among all its projective images satisfying that bound,
the least primitive radius exists: the set is nonempty and squared
radii are positive integers. For `M` large, (13) places that minimum
between

```text
R^(1-O(1/M)) and R.
```

Thus a high-height transformation cannot furnish a fixed-power radius
descent while preserving all `M` endpoint points. This complements
[projective_radius_inflation.md](projective_radius_inflation.md), which
excludes small matrices by the destruction of shared Gaussian divisors.
Here large matrices are permitted, and cross-ratio invariance controls
the radius independently of how their new shared factors are arranged.

The conclusion does not classify all minimizing matrices, prove that
the initial radius is exactly minimal, or exclude a smaller descent
inside the factor `exp(O(log R/M))`. There is also no reason established
here that a nontrivial endpoint-preserving map must exist at all.
Configurations in different projective orbits can have unrelated
cross-ratio data and arbitrarily large heights within these necessary
conditions. A radius-independent point count therefore remains open.

## 6. Stability after deleting points

The same allocation model controls the common factor exposed by
discarding rows. Start with a primitive `k`-point configuration of radius
`R` in an arc of length at most `C sqrt(R)`, and retain an arbitrary set
`S` of `m>=1` rows. The radius-one case has at most four points and
needs no logarithmic argument; assume `R>1`. In the normalization of
Section 1 put

```text
gamma_S=sum_p [min_(i in S) a_ip+e_p-max_(i in S) a_ip] log p.
```

At `p`, the common Gaussian divisor of the retained rows has exponents
`min_S a_ip` at `pi_p` and `e_p-max_S a_ip` at `bar(pi_p)`. After this
divisor is removed, the least primitive squared radius of the retained
directions is therefore

```text
R_S^2=product_p p^(max_S a_ip-min_S a_ip),
gamma_S=log(R^2/R_S^2)=2 log(R/R_S).           (14)
```

Equivalently, `gamma_S` is the total threshold-layer weight on which
all retained sign rows have the same value. This identity uses the
primitive parent normalization: inert and ramified common factors have
already been removed, and literal Gaussian units do not affect it.

Set `J=ceil(C/2)`. Partition the original arc into `J` subarcs of length
at most `2 sqrt(R)`, assign each retained row to one subarc, and then
partition by its four possible literal Gaussian units. One cohort has
size

```text
L>=ceil(m/(4J)).                               (15)
```

Every threshold layer constant on all of `S` is constant on this cohort.
For its sign rows `f_1,...,f_L`, the chord calculation in Section 1 gives
`mu_ij<=0` for every distinct pair, while the squared-norm identity gives

```text
E_layers (sum_(i=1..L) f_i)^2
 = L+sum_(i!=j) mu_ij <= L.
```

The constant layers alone contribute `L^2 gamma_S/W` to the left side.
Consequently

```text
L^2 gamma_S/W<=L,       gamma_S<=W/L,
log R_S>= (1-1/L) log R.                       (16)
```

The division by the retained common factor preserves all directions and
can only improve the normalized arc constant: if their angular span is
`Theta`, then

```text
Theta sqrt(R_S)<=Theta sqrt(R)<=C.
```

Now suppose a rational projective change of half-angle coordinates sends
these `m` retained directions to another primitive endpoint configuration
of radius `R'` and fixed normalized constant `C'`. Applying Section 4 to
the primitive source of radius `R_S` and combining it with (16) yields,
for `m` large,

```text
log R'/log R_S=1+O_(C,C')(1/m),
log R'>=(1-O_(C,C')(1/m)) log R.                (17)
```

Thus deleting `O(1)` points before one such transformation still cannot
produce a fixed-power descent `R'<=R^(1-epsilon)` as `k` tends to
infinity. This applies to rational rotations and changes of
stereographic reference, and more generally whenever a proposed
recircularization induces a `PGL_2(Q)` map on the rational circle
directions. An off-center inversion followed by a translation requires
that condition to be checked; rational integrality does not follow merely
from writing down the inversion.

Equation (17) is only a one-step obstruction. It permits a factor
`R^O(1/m)`, and an iteration that repeatedly deletes points can accumulate
the successive `O(1/m)` losses. It therefore does not exclude a new
multi-step induction or prove the desired radius-independent bound.

## Verification

The uniformity-audit and fresh-algebraic agents independently read the
written proof. Both checked the conductor noncancellation, the quartic
sampling factor, the moment bounds and finite constants in (10), and
the inherited-height convention after intersecting unit classes.
The fresh-algebraic audit also checked the equality of the half-angle
and complex radial cross-ratios.

Exact finite checks covered 89,270 four-row allocation tuples with
prime exponents between one and twelve. They verified the conductor
valuation identity, absence of opposite-sign threshold contributions
at one prime, and the allocation inequalities used to prevent residue
cancellation. Enumeration of all ordered quadruples for 56 sign
populations with four through ten rows verified the sampling polynomial
(6). These checks supplement the algebraic and averaging proofs.
