# A uniform endpoint bound for bounded deletions and augmentations of Hadamard profiles

For each fixed number of deleted rows and extra cut columns, the family
below has an effective point bound independent of the circle radius.
The argument uses an ambient Hadamard symmetry, integer subset collisions,
and the effective phase estimate in
[effective_hadamard_short_character_phase_gap.md](effective_hadamard_short_character_phase_gap.md).
It does not prove that arbitrary endpoint clusters contain such a family.

This extends the complete-core combinatorics in
[hadamard_augmentation_short_character.md](hadamard_augmentation_short_character.md).
In particular it includes the profiles obtained by deleting two rows
from a Hadamard matrix: their remaining number of rows is congruent to
two modulo four, so they need not have a square Hadamard core of their
own order.

## 1. Exact family and endpoint regime

Let `H` be a Hadamard matrix of even order `N`, with constant first
column, and write `H0` for its other `N-1` columns. Delete a set `D` of
`d` rows and retain `I`, of size `M=N-d`. The exact sign-class matrix of
the Gaussian rows is, after reordering and orienting cut columns,

```text
S = [ H0|I , F ],              F in {+1,-1}^{M x k}.       (1)
```

Here all columns are distinct nonconstant unoriented whole cut classes,
and every associated block is a conjugate-primitive nonunit. Duplicate
extra cuts are first combined with their existing whole class. The
factorization is literal:

```text
z_i = g epsilon_i product_j Gamma_j^((1+s_ij)/2)
                              bar(Gamma_j)^((1-s_ij)/2),
W = log R^2.                                               (2)
```

Moving prime factors, shared primes in distinct nested threshold layers,
and the complete common Gaussian factor `g` are allowed. This is the
same factorization as in
[corank_one_local_bessel_pair_filter.md](corank_one_local_bessel_pair_filter.md).
The endpoint assumptions are either

- arbitrary literal row units and `0<C<=sqrt(2)`; or
- one common literal row unit and `0<C<=2`.

The points are distinct and lie in one circle arc of length at most
`C sqrt(R)`. We use `M>=5`, hence `R>1`. Under these hypotheses, the
whole-pattern Bessel bound and affine allocation rigidity give

```text
log Norm(Gamma_j) <= 2W/M,
W >= (M-1) log 5.                                        (3)
```

For the second assertion, affine rigidity requires at least `M-1`
distinct varying split rational primes. Their total contribution to
`R^2` is at least `5^(M-1)`. Neither the number of threshold columns
nor independence of the block prime supports is substituted for that
prime count.

If `d<N/2`, restriction preserves every old cut as a nonconstant cut,
and preserves distinctness even modulo sign. Each original column has
`N/2` entries of either sign, and each pair of original columns agrees
in `N/2` rows and disagrees in `N/2` rows. Deleting fewer than `N/2`
rows cannot remove any of these four required distinctions.

## 2. A certificate that vanishes on the deleted rows

Extend each extra column by zero on `D`, and define integer Fourier
numerators for the old columns:

```text
w_jl = sum_(i in I) H_ij F_il,       1<=j<=N-1, 1<=l<=k.
```

Hadamard orthogonality, including the unused constant coefficient,
gives

```text
sum_(j=1)^(N-1) w_jl^2 <= N M <= N^2.                   (4)
```

Suppose first that `k>=1`, and put

```text
B = floor sqrt(2k N^2/(N-1)),
n = floor((N-1)/2).
```

At least `n` old columns satisfy `|w_jl|<=B` for every `l`.
Indeed, every bad column contributes strictly more than
`2kN^2/(N-1)` to the sum of squares over all `j,l`, whose total is at
most `kN^2`.

Choose any `n` good columns. Record for each `q`-element subset the
sum of its `k` Fourier vectors and the sum of its entries in each of
the `d` deleted Hadamard rows. There are at most

```text
(2qB+1)^k (2q+1)^d                                     (5)
```

possible records. Thus the exact sufficient condition

```text
binom(n,q) > (2qB+1)^k (2q+1)^d                         (6)
```

gives two distinct subsets with the same record. Their difference is
a nonzero vector `v in {0,+1,-1}^{N-1}` with

```text
ell = ||v||_1 <= 2q,
w^t v = 0,
H0|D v = 0.                                             (7)
```

Define `lambda_full=H0 v` and retain its coordinates in `I` as
`lambda`. The deleted coordinates vanish by construction. Therefore

```text
sum_(i in I) lambda_i = 0,
lambda^t S = N(v,0),
||lambda||_2^2 = N ell,
L^2 := ||lambda||_1^2 <= M N ell <= N^2 ell.              (8)
```

All these are integer identities. The extra coordinates in the second
identity vanish by `w^t v=0`; the old coordinates follow by using
`H0^t H0=N Id` and then removing only zero coordinates of
`lambda_full`. In particular `lambda` is nonzero.

For `k=0`, use all old columns and only the deleted-row records. With
`q=1`, these are sign vectors in `{+1,-1}^d`. Consequently

```text
N-1 > 2^d                                               (9)
```

already gives (7)--(8), with `ell=2`. This also covers `d=0` for
`N>=4`; when no rows are deleted one could instead use a single
Hadamard column and `ell=1`.

## 3. Effective phase exclusion with the ambient denominator

The signed source product from (8) is

```text
product_i z_i^lambda_i
   = u (Beta/bar(Beta))^(N/2),          u in mu_4,         (10)
V := log Norm(Beta) <= 2 ell W/M.
```

Here `Beta` is obtained from the signed old-block product by dividing
out all rational Gaussian content. The exact valuations include all
shared-prime cancellations. Affine allocation rigidity makes `Beta`
a nonunit: otherwise (10) would be a forbidden unit-valued nonzero
affine multiplicative relation among the source points. Thus
`alpha=Beta/bar(Beta)` is not a root of unity.

Lift the point arguments in their containing arc. Squaring (10) removes
the quarter-turn unit obstruction and gives a nonzero two-logarithm form

```text
Lambda = N log(alpha) - m log(-1),       m in Z,
0 < |Lambda| <= L C exp(-W/4).                         (11)
```

Take the principal logarithm of `alpha` and `log(-1)=pi i`. Then

```text
|m| <= N + LC/pi,
B_star := max(2,N,|m|) <= 2N sqrt(ell).                 (12)
```

The last inequality uses `C<=2` and (8). This is the denominator benefit
of the symmetry: the coefficient is the ambient order `N`, even though
some rows have been removed.

The degree-two Matveev specialization proved in the effective companion,
with `K=2^40`, gives

```text
log |Lambda| > -K V log(e B_star).
```

Combining it with (3), (10)--(12) shows that

```text
M >= 16 K ell log(2eN sqrt(ell))                        (13)
```

would force

```text
W < 8 log(2N sqrt(ell)).                               (14)
```

In particular, there is no profile satisfying both (13) and

```text
(M-1) log 5 >= 8 log(2N sqrt(ell)).                     (15)
```

The analytic estimate is effective and uniform in the moving Gaussian
blocks, their widths, the common factor, and the row units in the stated
regimes. It does not use an ineffective fixed-target Roth cutoff.

## 4. A radius-independent bound for fixed `k,d`

For `k>=1`, choose

```text
q = floor(k/2)+1.
```

Then (6) holds for every sufficiently large `N`, uniformly over the
Hadamard matrix, the deleted rows and the extra columns. Its left side
has order `N^q`; its right side is `O_(k,d)(N^(k/2))`; and `q>k/2`.
Since `ell<=2q` is bounded while `M=N-d` grows linearly with `N`,
(13) and (15) then hold as well. For `k=0`, (9) provides the same
conclusion with `q=1` and `ell=2`.

For completeness, the following deliberately large integer cutoff makes
the effectivity explicit. Set `q=floor(k/2)+1` for all `k>=0`, and put

```text
A(k,d) = 2^d+2                                      if k=0,
A(k,d) = (4q)^(2q) (25q^2 k)^k (2q+1)^(2d) + 1     if k>=1,

U(k,d) = max(2d+1, 4q, A(k,d), (256 K q)^2+1),
K = 2^40.                                             (16)
```

There are no such endpoint profiles with `N>=U(k,d)`. Thus, in
particular, their point count satisfies `M<U(k,d)`, independently of
`R` and uniformly throughout the specified range of `C`.

To verify the combinatorial part of this cutoff, for `N>=4q` use
`n>=N/4`, `binom(n,q)>=(n/q)^q`, and `B<=2sqrt(kN)`. The left side
of (6) is at least `(N/(4q))^q`, whereas its right side is at most
`(5q sqrt(kN))^k (2q+1)^d`. Squaring, it suffices that

```text
N^(2q-k) > (4q)^(2q) (25q^2 k)^k (2q+1)^(2d),
```

which follows from `N>=A(k,d)` since `2q-k>=1`.
For the analytic part, `N>=2d+1` gives `M>N/2`. The last term of
(16) implies `N>(256Kq)^2`, as well as `N>=16q^2` and `N>1024`.
Consequently

```text
32Kq log(2eN sqrt(2q)) <= 64Kq log N
                      <= 64Kq sqrt(N) < N/4 < M,

8 log(2N sqrt(2q)) <= 16 log N <= 16sqrt(N)
                  < N/2 <= (N/2-1)log 5 < (M-1)log 5.
```

These are stronger than (13) and (15) for every `ell<=2q`.

## 5. Scope and verification

The theorem is a uniform point-count result for this symmetry family;
its bound is not a new estimate for general endpoint arcs. A general
cut matrix need not be an augmentation of a bounded-row deletion of
an ambient Hadamard matrix. Partitioning an arbitrary `C sqrt(R)` arc
into smaller arcs also does not preserve hypothesis (1), so the
restricted unit/arc regimes cannot be silently removed that way.

The [exact checker](check_deleted_hadamard_uniform_family_bound.py)
constructs integer certificates on finite Sylvester fixtures and checks
the deleted-row conditions, the full source exponent identity, the
nonzero certificate and norm bounds, and the distinctness of old cuts.
Finite checks validate these algebraic calculations; the prime-count
and logarithmic-form statements are proved inputs, not conclusions of
those computations.
