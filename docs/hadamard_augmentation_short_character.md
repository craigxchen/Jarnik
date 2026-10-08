# Short integer characters for Hadamard augmentations

A complete Hadamard core plus a bounded number of arbitrary sign columns
has a short integer character with an explicit row certificate. The bounds
below are uniform in the Hadamard matrix and in the extra columns. Combined
with [the effective phase estimate](effective_hadamard_short_character_phase_gap.md),
they give an effective bound on the row count in this family, even when
the row count and matrix vary. This improves the fixed-matrix application
of Roth, but does not prove that general low-redundancy profiles have such
a core or that their ratio `T/M` must diverge.

## 1. Exact lattice and certificate

Let `H=[1,H_0]` be a normalized Hadamard matrix of order `M>=4`:
its entries are signs and `H^t H=M I`. Let `F` consist of `k>=1`
additional sign columns, and put

```text
S=[H_0,F],              T=M-1+k,
w_j=((H_0)^t F)_j in Z^k,             1<=j<=M-1.
```

For the whole-pattern interpretation, all columns of `S` represent distinct
nonconstant unoriented cuts. The algebra below remains valid without this
restriction, but such repetitions are not separate whole-pattern blocks.

Suppose two distinct `q`-subsets of the old column indices have the same
sum of their `w_j`. The difference of their indicators is a nonzero vector
`v in {-1,0,1}^{M-1}` satisfying

```text
(H_0^t F)^t v=0,          ell:=||v||_1<=2q.
```

Extend `v` by zeros in the extra coordinates. Its literal integer row
certificate is

```text
lambda=H_0 v,
sum_i lambda_i=0,
lambda^t S=M(v,0),
||lambda||_2^2=M ell,       ||lambda||_1<=M sqrt(ell).       (1)
```

Consequently `(v,0)` belongs to
`L_S=Z^T intersect {mu^t S:mu in Q^M, sum mu_i=0}`.
In the notation of [the integer pattern-lattice theorem](corank_one_local_bessel_pair_filter.md),
the explicit denominator is `2t=M`. No determinant, unknown integral
saturation, or bound on the prime exponents enters (1).

The `k`-dimensional column kernel is, explicitly,

```text
ker(S_i-S_0)={(-(H_0^t F)y/M,y):y in Q^k}.               (2)
```

This is a relation among columns. It is not a row barycentric dependence,
and none of its entries is an actual logarithmic block norm.

## 2. Parseval and a finite pigeonhole inequality

For each extra column `f`, Parseval gives

```text
(1^t f)^2+sum_j (H_0^t f)_j^2=M^2.
```

Hence `sum_j ||w_j||_2^2<=k M^2`. At least `M/2` old column indices obey

```text
||w_j||_2^2<=2k M^2/(M-1).
```

Indeed, more than half of the `M-1` indices cannot exceed twice their
upper average. Put

```text
B=floor(sqrt(2k M^2/(M-1))),        N=M/2.
```

Choose any `N` of these good indices. Every sum over a `q`-subset lies
in the integer box `[-qB,qB]^k`. Therefore

```text
binom(N,q)>(2qB+1)^k                                  (3)
```

guarantees the nonzero character (1). This is a finite, constructive
pigeonhole argument: store the vector sum and the first subset producing
it, and stop at the first repetition. Overlap of the two subsets only
decreases `ell`.

For fixed `k`, taking `q=floor(k/2)+1` makes the left side grow as
`M^q`, while the right side is `O_k(M^(k/2))`. Thus (3) holds for all
sufficiently large `M`, with an effectively computable threshold. This
improves the elementary box estimate `|w_j,l|<=M`, which instead requires
`q>k` in the same argument.

There is also a simple bound when `k` grows. Set

```text
a=ceil(log_2(2k)),             q=4ka.
```

Then

```text
M>=64ka  implies  lambda_1(L_S,ell_1)<=8ka.              (4)
```

To prove this without asymptotic notation, use `B<=2sqrt(kM)` and
`binom(M/2,q)>=(M/(2q))^q`. Since `2qB+1<=5q sqrt(kM)`, it suffices
that

```text
(M/(2q))^q > (5q sqrt(kM))^k.
```

The ratio increases with `M`, as `q>k/2`. At `M=16q=64ka` this reduces to

```text
2^(12a)>160 k^2 a^(3/2).
```

Use `k<=2^(a-1)` and `2^(10a)>40a^(3/2)`, valid for every integer
`a>=1`. The latter follows by squaring and induction from
`2^20>1600`; its successive left-side ratio is `2^20`, while the
right-side ratio is at most eight. This proves (3), hence (4).

## 3. One extra column has an even shorter character

When `k=1` and `M>=8`, one always has `lambda_1(L_S,ell_1)<=2`.
If any old Fourier numerator is zero, its standard coordinate works.
Otherwise, if two numerators have the same absolute value, their signed
coordinate sum or difference works. If neither occurs, the `M-1`
distinct positive absolute values have squared sum at least

```text
1^2+...+(M-1)^2=(M-1)M(2M-1)/6>M^2,
```

contradicting Parseval. The strict numerical inequality already holds for
`M>=5`; an admissible Hadamard order in that range is at least eight.
These certificates also have (1), even if the two signs add rather than
subtract. The corank-one pair test therefore excludes a fixed such
profile at unbounded endpoint radius once `M>16`. The companion effective
phase estimate removes the fixed-`M` qualification for a sufficiently
large, explicit `M`.

## 4. Arithmetic use and scope of the reduction

On an actual endpoint cluster in the unit ranges of the
[allocation rigidity theorem](linear_allocation_affine_rigidity.md),
the exact certificate (1) supplies

```text
product_i z_i^lambda_i=u (Beta/bar(Beta))^(M/2),
log Norm(Beta)<=2ell W/M,             W=log(R^2).         (5)
```

Here `Beta` is the conjugate-primitive reduction of the signed old-block
product, and `u` is the literal Gaussian row-unit product. Nested prime
layers cancel exactly against the full matrix, including the extra
columns. The common Gaussian divisor cancels in the signed product but
remains in the original `W` and angular error. Affine allocation rigidity
makes `Beta` a nonunit. Equation (1) controls both the integer exponent
and the phase error, permitting a two-logarithm estimate over `Q(i)`
uniformly in `M`; see the [effective proof](effective_hadamard_short_character_phase_gap.md).
It must not be replaced by an assertion that Roth's fixed-target constant
is uniform as `M` grows.

The present norm constraints alone do not extract this complete core.
For any extras `F`, assigning weight one to each old column and total
weight less than one to the extras makes every off-diagonal row Gram
entry negative: its old contribution is exactly `-1`. Thus the family
is a substantial class of weighted obtuse profiles, but the hypothesis
is additional structure.

In fact a same-row-count Hadamard extraction statement is false, even at
redundancy two. Start with a normalized Hadamard matrix of order
`N=2^d>=8`; delete two rows and its constant column. Then

```text
M=N-2,          T=N-1=M+1,
<S_i,S_j>=-1 for i!=j.
```

All columns remain distinct nonconstant unoriented cuts: originally a
column has `N/2` signs of each kind and two different columns agree and
disagree at `N/2` positions, so deleting two rows cannot make either
count zero. But `M=2 mod 4` admits no Hadamard matrix when `M>2`.
The elementary obstruction follows by normalizing one row and splitting
the equal plus/minus positions of a second row against a third orthogonal
row; their four sign classes must each have `M/4` entries. These truncated
profiles still have short characters; the example refutes only this
literal structural extraction, not a possible general bound on `L_S`.
The [bounded-deletion extension](deleted_hadamard_uniform_family_bound.md)
handles them by imposing zero row-certificate coordinates at deleted
rows of the larger matrix.

Rational independence of distinct prime logarithms does not repair a
missing extraction. Strict weighted Gram feasibility is open. Assign one
distinct split prime to each column and choose its positive integer power
so that its logarithmic norm approximates a large multiple of any strictly
feasible weight vector. The approximation error per coordinate is bounded
by the fixed prime logarithm, so all strict inequalities eventually hold.
The resulting logarithmic block norms are rationally independent, while no phase
control follows. Neither these norms nor (2) can be substituted for row
central weights.

## Verification

[check_hadamard_augmentation_short_character.py](check_hadamard_augmentation_short_character.py)
checks exact Parseval identities and integer certificates, exhausts the
one-extra-column case at Walsh order eight, includes order-twelve Paley
fixtures, and verifies the finite pigeonhole bounds at integer parameters.
It also checks the truncated-Hadamard obstruction. Its sign examples have
no asserted small angular span; the arithmetic hypotheses and effective
linear-form theorem are supplied by the companion note.
