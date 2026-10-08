# An effective phase gap for a short Hadamard character

An exact Hadamard core with only a small number of additional sign
classes admits a short integral row character.  For endpoint arcs in
the unit regimes below, that character gives an **effective**
large-row exclusion by Matveev's two-logarithm bound, without Roth's
ineffective fixed-target constants.  In particular a fixed number of
extra classes cannot accompany Hadamard cores of unbounded order.
More quantitatively, such a profile at large order must have at least
on the order of `M/(log M)^2` extra classes, with a very small explicit
coefficient.  This is a result for the specified exact sign profile;
it does not imply a general bound on endpoint clusters.

Use either a literal common row unit with `0<C<=2`, or arbitrary
row units with `0<C<=sqrt(2)`.  Let `M>=8` be a Hadamard order and
`H_0` the `M x (M-1)` matrix obtained by deleting the all-ones
column of a normalized Hadamard matrix.  Let `F` be `M x k`, with
`k>=0` additional sign columns.  The columns represent **distinct
whole unoriented nonconstant cut classes**; orient them once and
group all nested prime-power layers of each class into a Gaussian
block `Gamma_j`.  Assume the `M-1` Hadamard-core classes are genuinely
present, so each core block is nonunit.  The source rows have the
literal factorization

```text
z_i=d epsilon_i product_j Gamma_j^((1+s_ij)/2)
                         bar(Gamma_j)^((1-s_ij)/2),
S=[H_0,F],  epsilon_i in {1,-1,i,-i},
W=log R^2,     D=log Norm(d)>=0,
V_j=log Norm(Gamma_j),
Delta<=C exp(-W/4).                                    (1)
```

The factor `d` is the complete common Gaussian gcd, including
fixed split allocations and ramified or inert factors.  Blocks from
different classes may share a rational prime through nested
threshold layers; the exact sign factorization and their
conjugate-primitivity are the requirements.  The core nonunits give
the source lower bound

```text
W=D+sum_j V_j >=(M-1)log 5.                            (2)
```

## 1. A general short-character theorem

Suppose there is a nonzero integral vector
`v in {-1,0,1}^(M-1)` with

```text
v^t H_0^t F=0,      ell=||v||_1>=1,
lambda=H_0 v in Z^M.                                  (3)
```

Hadamard orthogonality gives the exact certificate

```text
sum_i lambda_i=0,
lambda^t S=M(v^t,0),
L=||lambda||_1<=M sqrt(ell).                            (4)
```

The `L` bound follows from
`||lambda||_2^2=M||v||_2^2=M ell` and Cauchy--Schwarz.

Let `Beta_v` be the conjugate-primitive reduction of the oriented
Gaussian product of the core blocks prescribed by `v`, exactly as in
[the integer pattern-lattice theorem](corank_one_local_bessel_pair_filter.md).
Primewise signed-row cancellation, including nested layers and the
common factor, gives

```text
Q=product_i z_i^lambda_i
 =u (Beta_v/bar(Beta_v))^(M/2),
u=product_i epsilon_i^lambda_i in mu_4.                (5)
```

Because `lambda!=0`, endpoint affine allocation rigidity makes `Q`
nonunit for `R>1`.  Thus `Beta_v` is nonunit.  Its conjugate-primitivity
implies that it has no ramified common factor, has norm at least five,
and that `alpha=Beta_v/bar(Beta_v)` is not a root of unity.  The
whole-pattern Bessel estimate gives

```text
V=log Norm(Beta_v)
 <=sum_j |v_j|V_j<=2 ell W/M.                        (6)
```

Take the principal logarithm `log alpha=i psi`, `|psi|<=pi`.
For the row-unit product write `u=i^s`, `s in Z`.  Lift all source
arguments on their containing arc.  Equation (5), the zero row-sum
in (4), and the angular width in (1) yield an integer `n` with

```text
|(M/2)psi+(pi/2)s-2pi n|<=L Delta/2.
```

After **doubling** this relation, put `m=4n-s`.  The actual
two-logarithm form and its integer coefficient bound are

```text
Lambda=M log alpha-m log(-1) !=0,
|Lambda|<=L Delta<=M sqrt(ell) C exp(-W/4),
B=max(M,|m|)<=M+L Delta/pi
                  <=M(1+C sqrt(ell)/pi).               (7)
```

The nonzero assertion follows because `alpha^M=(-1)^m` would make
`alpha` a root of unity.  Independent units merely change `m` by
the recorded integer `s`; they do not introduce a moving target.

For conjugate-primitive `Beta_v`, the primitive degree-two polynomial
of `alpha` has leading coefficient `Norm(Beta_v)`, so
`h(alpha)=V/2`.  The degree-two specialization of Matveev's
[Corollary 2.3](minimal_prime_rank_width_bound.md) therefore applies
to this **composite** Gaussian block with the same `K=2^40`:

```text
log|Lambda|>-K V log(eB).                               (8)
```

Combining (6)--(8) gives the fully effective gap inequality

```text
gamma(M,ell,C) W < log(M sqrt(ell) C),
gamma=1/4-(2K ell/M)
        log(eM(1+C sqrt(ell)/pi)).                      (9)
```

If `gamma>=1/8`, then
`W<8 log^+(M sqrt(ell) C)`; together with (2)
this excludes every such source whenever

```text
(M-1)log 5 > 8 log^+(M sqrt(ell) C).                   (10)
```

One convenient sufficient condition for `gamma>=1/8` in both
unit regimes, since `C<=2`, is

```text
M>=16K ell log(2eM sqrt(ell)).                          (11)
```

The constants are huge but all thresholds in (9)--(11) are
computable from `M,ell,C`.

## 2. Bounded augmentation and a quantitative redundancy cost

The exact combinatorial construction in
[hadamard_augmentation_short_character.md](hadamard_augmentation_short_character.md)
applies to any `M x k` additional sign matrix `F` when `k>=1`.
With

```text
a=ceil(log_2(2k)),       q=4ka,                         (12)
```

it gives, for `M>=16q`, a nonzero vector in (3) with entries
`-1,0,1` and `ell<=2q`.  It uses Parseval to find many old columns
whose extra-column Fourier vectors are short, then a collision of
two equal-size `q`-subset sums.  Its checker verifies the exact
Hadamard certificate and the finite counting threshold.

Equations (9)--(11) now give a direct effective exclusion.  There is
no endpoint profile (1) if

```text
M>=max(128,16q)
and
M>=32K q log(2eM sqrt(2q)).                           (13)
```

Indeed, (11) holds for `ell<=2q`.  Also `q<=M/16` makes
`8 log(2M sqrt(2q))<(M-1)log 5` for `M>=128`, so (10) holds.
For each fixed `k`, (13) eventually holds as `M` grows, giving a
uniform **effective** row bound independent of source primes, powers,
complete common content, and independent row units.

The same inequalities permit `k` to grow with `M`.  For `k>=1` and
`M>=128`, any such actual endpoint profile must satisfy the coarse explicit
necessary condition

```text
k > M/[2048 K (log M)^2],       K=2^40,                (14)
```

where logarithms are natural.  To check it, suppose the opposite and
`1<=k<=M`.  Then `a<=2 log_2 M`, so
`q<=8k log_2 M`.  The hypothesized upper bound makes
`q<=M/16`.  Moreover
`log(2eM sqrt(2q))<=2log M`, hence the second lower bound in (13)
is at most
`(512K/log 2)k(log M)^2<M`.  Thus (13) would exclude the
profile.  If `k>M`, (14) is immediate.  The coefficient in (14) is
not meant to be sharp; the meaningful feature is the required
`M/(log M)^2` redundancy in this structured family.

For `k=0`, a single old coordinate has `v=e_j`, `ell=1`, and (9)--(11)
give an even simpler effective large-order exclusion.  No
pigeonhole lemma is needed.

## 3. What this resolves

The ineffective Roth constant in the fixed-Hadamard argument is
replaced by a degree-two linear-form estimate for the composite
`Beta_v`.  Its logarithmic coefficient is only `O(M sqrt(ell))`
because the exact Hadamard certificate is explicit; a general
moving allocation inverse can have a factorial determinant and does
not meet (11).  Arbitrary endpoint profiles need not contain any
Hadamard core or a short vector in the integral row-difference
lattice.  In particular this theorem does not improve the general
`M(R)=O(log R/loglog R)` rank bound.

There is also a fully effective elementary substitute for the
fixed-grid Roth lemma, with a sharp degree cost.  For a Gaussian
prime `pi` of norm `p` and a rational-`pi` target `beta`, put
`zeta=exp(2i beta)` and `d_beta=[Q(i,zeta):Q(i)]`.  The nonzero
Gaussian integer
`Norm_(Q(i,zeta)/Q(i))(pi-zeta bar(pi))` has modulus at least one.
Its `d_beta-1` other conjugate factors have modulus at most
`2sqrt(p)`, while the near factor has modulus at most
`2sqrt(p) dist(arg(pi),beta)`.  Hence

```text
dist(arg(pi),beta)>=2^(-d_beta) p^(-d_beta/2).           (15)
```

In the full-rank width argument, let `d_j` be the degree of the
torsion target isolated for prime `p_j`, and let `K_0 C R^(-1/2)`
be its phase-error bound.  Equation (15) gives
`p_j>=R^(1/d_j)/[4(K_0 C)^(2/d_j)]` and an effective radius
power gap precisely when `sum_j e_j/d_j>2`.  Quadratic targets with
total width at least five qualify.  The available uniform determinant
bound is `|det B|<=r! H_width^r`, where `H_width` is the explicit
allocation-width cap in
[the full-rank theorem](minimal_prime_rank_width_bound.md).  A primitive torsion target
of order `4|det B|` can have degree comparable to `|det B|`.
Thus the reciprocal-degree criterion gives no general rate as
`r=M-1` grows.  Likewise the available Matveev coefficient bound in
that moving inverse has `log B=O(r log r)` after the effective width
bound; substituting it into a Bessel height `O(W/M)` would demand
`M` much larger than `K r log r`, which is not supplied by `M=r+1`.

The companion [checker](check_effective_hadamard_short_character_phase_gap.py)
verifies exact sign certificates, a literal Gaussian signed product
with independent units and common content, and the finite threshold
algebra.  Matveev and the obtuse Bessel lemma are cited inputs, not
finite-search assertions.
