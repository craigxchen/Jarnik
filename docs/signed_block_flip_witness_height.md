# Height cost of flipping signed Gaussian blocks

This note sharpens the deterministic part of the composite-block transfer in [the Gaussian moat note](gaussian_moat_composite_block_transfer.md). The source of the common-divisor determinant step is [*Bounded-Step Walks on Gaussian Primes*, Lemma 4.2](https://github.com/openai/math/blob/adc7f1241b42e322a6451854ab7e4b4c146bf78a/preprints/Bounded-Step-Walks-on-Gaussian-Primes-September-26-2026/paper.pdf). The witness-height bound below is an additional deduction. It does not bound the number of lattice points on a short circle arc.

## A forbidden interval for two signed witnesses

Let `H_1,...,H_r` be nonunit Gaussian integers with `gcd_G(H_j,bar(H_j))=1` and pairwise disjoint rational-prime norm supports. Put `n_j=N(H_j)` and `P=prod_j n_j`. For each sign vector `sigma`, let `B_sigma` choose `H_j` or `bar(H_j)` independently in the product. For two sign vectors `sigma,tau`, put

```text
q = prod_{j: sigma_j != tau_j} n_j,
N(gcd_G(B_sigma,B_tau)) = P/q.
```

Let `Q` be an origin-centered rectangle in any orientation, with half-lengths `R_1 >= W_1 >= 1` and `A=R_1 W_1`. Suppose `Q` contains nonzero witnesses `x in B_sigma Z[i]` and `y in B_tau Z[i]`. Then necessarily

```text
q <= 2 R_1^2/P       or       q >= P/(2A).       (1)
```

Equivalently, the open interval

```text
2 R_1^2/P < q < P/(2A)                         (2)
```

contains no sign-change norm for which both signs have witnesses in `Q`. The two thresholds are exact outputs of the argument; no claim of necessity outside this interval is made.

Independently, every witness `x=B_sigma u in Q` has `|u|<=sqrt(2)R_1/sqrt(P)`. At endpoint height this leaves only a multiplier of size `exp(o(w))`; a proposed ghost construction must meet that budget as well as its divisibility conditions.

**Proof.** If `q<P/(2A)`, the common Gaussian factor `c=gcd_G(B_sigma,B_tau)` has `N(c)=P/q>2A`. Hence `N(c)` divides the integer `det(x,y)`, while `|det(x,y)|<=2A`; thus `x,y` are collinear. Write `x=a v`, `y=b v`, with `v` a primitive direction in `Z^2` and nonzero `a,b in Z`. Set `G_sigma=N(gcd_G(v,B_sigma))` and likewise `G_tau`. Primewise divisibility gives the exact least positive integer multipliers

```text
a_sigma = P/G_sigma,       a_tau = P/G_tau.
```

For nonnegative Gaussian-prime valuations `t,u,u'`,

```text
min(t,u)+min(t,u') <= t+min(u,u').
```

Multiplying this inequality over all Gaussian primes yields

```text
G_sigma G_tau
  <= N(v) N(gcd_G(B_sigma,B_tau))
  = |v|^2 P/q.
```

Since `|a|>=a_sigma`, `|b|>=a_tau`, the larger witness has length at least

```text
max(|x|,|y|)
  >= P|v|/min(G_sigma,G_tau)
  >= sqrt(Pq).                                    (3)
```

Every point of `Q` has modulus at most `sqrt(2)R_1`; thus `Pq<=2R_1^2`. This proves (1).

The collinear height estimate (3) is sharp: for `H=2+i`, `P=q=5`, both conjugate sign lattices contain the real point `5`, whose length is `sqrt(Pq)=5`. The rectangle determinant condition is the separate ingredient that forces collinearity for nearby signs.

If `B_sigma K` is itself an actual witness with primitive integer coordinates, the collinearity conclusion is sharper. Every collinear Gaussian integer is then an integer multiple `t(B_sigma K)`. On every flipped block, the orientation opposite to that selected in `B_sigma` is coprime to `B_sigma K` (provided the actual witness is conjugate-primitive), so `B_tau | t(B_sigma K)` forces `q|t`. Any opposite-sign witness on that same line consequently has length at least `q |B_sigma K|`, rather than merely `sqrt(Pq)`.

## Consequence for an extracted near-real row

Fix `m>=3` and write `r=2^(m-1)`. An actual endpoint row has `r` incident balanced blocks and the exact form

```text
P_i = K_i prod_{T containing i} H_T,
log N(H_T)=w+o(w),
log max(1,|K_i|,|Im P_i|)=o(w),
gcd_G(P_i,bar(P_i))=1.
```

Throughout this application the entire block family `{H_T}` is
conjugate-primitive and has pairwise disjoint rational-prime norm supports,
as guaranteed by [the extraction](endpoint_full_profile_quantifiers.md).
The correcting factors may overlap those supports.

Let `P` denote the norm product of its `r` incident blocks. Then `log P=rw+o(w)` and `|P_i|=exp(rw/2+o(w))`. A near-real origin-centered rectangle containing `P_i` can be chosen with

```text
R_1=exp(rw/2+o(w)),       W_1=exp(o(w)),
A=exp(rw/2+o(w)).                             (4)
```

If another sign vector flips `h` of these incident blocks, then `log q=hw+o(w)`. For every fixed `1<=h<r/2`,

```text
log(P/q)-log(2A)=(r/2-h)w+o(w)>0,
log(Pq)-log(2R_1^2)=hw+o(w)>0.
```

Equation (2) therefore forbids **every** near-real ghost witness at those sign distances. Since `P_i` is primitive as a vector in `Z^2`, the sharper line argument says a hypothetical collinear flipped witness would have to be at least `q|P_i|=exp((r/2+h)w+o(w))` long. Flipping exactly `r/2` blocks lies at the determinant threshold and is not excluded here; conjugating the whole row flips all `r` blocks and gives the obvious second near-real witness. For `m=2`, `r=2` and the interval `1<=h<r/2` is empty, consistent with the explicit balanced two-row family.

This rules out filling a Hamming neighborhood around an actual row by any Bezout or linear-combination construction that still produces a witness in the same endpoint rectangle. In the real-axis rectangle used here, conjugation strengthens the conclusion to every proper nonempty distance except `r/2`, as proved in the next section. Neither statement bounds the number of original rows, whose supports differ.

## Conjugation leaves only the half-distance boundary

Assume now that `Q` is parallel to the real and imaginary axes, so
`bar(Q)=Q`. Put `L=2R_1^2/P` and `U=P/(2R_1 W_1)`. If signs
`sigma,tau` have witnesses in `Q`, their conjugates do too. The pair
`sigma,-tau` has flip norm `P/q`. Applying (1) to both pairs gives
the exact simultaneous restrictions

```text
q notin (L,U),       P/q notin (L,U).                 (11)
```

This additional restriction does not apply to an arbitrarily tilted
rectangle. It is available for the actual near-real rows in (4).

Fix the number `r` of blocks, let `w` tend to infinity, and assume
`log n_j=w+o(w)` for every block and the height bounds (4). Then
`log L=o(w)` and `log U=rw/2+o(w)`. If `h` coordinates differ,
`log q=hw+o(w)` and `log(P/q)=(r-h)w+o(w)`. The first restriction
excludes `0<h<r/2`; the second excludes `r/2<h<r`. Thus, for all
sufficiently large `w`, every pair of signs admitting witnesses has
distance in

```text
{0,r/2,r}.                                          (12)
```

All such signs form an antipodally closed set: conjugate a witness
to reverse every sign. Regard signs as vectors in `{+1,-1}^r` and
choose one representative from each antipodal pair. Distinct chosen
vectors have inner product `r-2(r/2)=0`. There are at most `r`
nonzero mutually orthogonal vectors in `R^r`, so the number of
admissible signs is at most `2r`. If `r` is odd, the half distance
is not integral and there are at most two signs. The thresholds may
depend on the particular asymptotic family; no uniform finite-height
threshold is asserted from unspecified `o(w)` errors.

Equality `2r` requires the chosen signs to be the rows of a Hadamard
matrix. This is only a necessary combinatorial condition. The Pell
family (8) does give arithmetic equality when `r=2`: its four sign
lattices contain `U,y,bar(y),bar(U)`, with multipliers
`1,2+i,2-i,1` and imaginary parts `1,2,-2,-1`.

For the four-pattern boundary, allow a common signed group `D` and
consider

```text
DABC,   DA bar(B) bar(C),
D bar(A) B bar(C),   D bar(A) bar(B) C.                (13)
```

Let `a,b,c,d` be their respective group log norms and `S=a+b+c+d`.
Suppose each of `a,b,c` is bounded below by a positive proportion of
`S`, the groups have disjoint conjugate-primitive supports, and all
four products admit nonzero witnesses with log modulus at most
`S/2+o(S)` and `log max(1,|Im|)=o(S)`. Their three distinct
pair distances have log norms `a+b,a+c,b+c`. Each is bounded away
from both zero and `S`, so (11) forces each to be `S/2+o(S)`.
Consequently

```text
a=b=c=d=S/4+o(S).                                  (14)
```

In particular, four such witnesses are impossible with `D` a unit
and all three changing groups of positive asymptotic mass. The
positive-mass hypothesis matters: the two-group Pell family has
`C=D=1` and is not excluded. The three-group algebraic circuit in
[the four-term note](four_term_signed_circuit_height.md) likewise
does not provide four endpoint witnesses.

The surviving four-group pattern is now excluded arithmetically in
[the matching-norm proof](four_signed_witness_matching_norm_gap.md).
In fact, without any balance assumption it gives
`min(N(A),N(B),N(C))<=8 M^2 max(1,Y)^4`, where `M` and `Y` are the
largest correction modulus and imaginary coordinate of the four
witnesses. The proof uses two independent integer relations among
the repeated opposite-edge norm divisors. Thus (14) describes the
boundary of the rectangle argument, not an arithmetic construction
of four thin witnesses.

There is still no common-support collection of this kind containing
all the original endpoint rows. If `P_i=K_i prod_{T containing i}H_T`,
then the guaranteed thin product

```text
Z_ij=P_i bar(P_j)/G_ij,
G_ij=prod_{T containing i,j} N(H_T)
```

has support `{T: exactly one of i,j belongs to T}` and imaginary
part `-t_ij`. The residual correction is `K_i bar(K_j)`. Different
unordered pairs have different supports: their memberships on the
singleton subsets already distinguish them. Including the anchor
gives the analogous single-row supports. The bound `2r` therefore
does not bound the number of these rows or pair products.

The [conjugation checker](check_signed_witness_conjugation_code.py)
checks the complementary flip norms, orthogonality constraints,
four-group balance, and the exact Pell witnesses. The argument above
is deterministic; it needs no distributional assumption on the signs.

## Why another row cannot be multiplied into the missing blocks cheaply

For distinct endpoint rows `i,j`, among the `r=2^(m-1)` blocks incident to `i`, exactly `r/2` do not contain `j`. Let `B` be **any** signed product choosing one orientation on all blocks incident to `i`. For the blocks `T` containing `i` but not `j`, the core factorization of `P_j` contains neither orientation of `H_T`. Any overlap with `B` on those supports comes from its correcting factor `K_j`. Therefore

```text
N(gcd_G(P_j,B))
  <= N(K_j) prod_{T: i in T, j in T} N(H_T),
N(B/gcd_G(P_j,B))
  >= [prod_{T: i in T, j notin T} N(H_T)]/N(K_j)
  = exp((r/2)w-o(w)).                            (5)
```

The least Gaussian multiplier `L` with `B | L P_j` is, up to a unit, `B/gcd_G(P_j,B)`. Equation (5) gives

```text
|L| >= exp(rw/4-o(w)),
|L P_j| >= exp(3rw/4-o(w)),                     (6)
```

which exceeds the near-real length budget `R_1=exp(rw/2+o(w))`. Thus directly multiplying a different row cannot furnish the flipped witness at endpoint height. A linear combination of several rows can cancel before divisibility is tested; (5)--(6) do not bound that possibility. The general short-basis coefficient obstruction in [the row-lattice transference note](small_imaginary_row_lattice_transference.md) is consistent with this cost.

The [exact checker](check_signed_block_flip_witness_height.py) verifies
12,288 instances of the collinear height and content inequalities, and
128 primitive-row multiplier identities. Its correcting factors overlap
the original incident blocks. These finite checks supplement the proof;
the full interval and endpoint assertions do not depend on a search.

## Raw flips with the original multiplier

There is a stronger restriction if the multiplier is kept fixed. Let
`U=X+iY=FA` be conjugate-primitive, where `F=a+ib` is a nonunit Gaussian
factor and `q=N(F)`. The raw flipped product is

```text
V=bar(F) A=U bar(F)/F,
Im V=[(a^2-b^2)Y-2abX]/q.
```

Conjugate primitiveness makes both `a,b` nonzero: a nonunit factor on
either coordinate axis would share a nonunit factor with its conjugate.
Since `(a^2-1)(b^2-1)>=0`, we have `a^2 b^2>=q-1`. Consequently

```text
|Im V| >= 2 sqrt(q-1)|X|/q-|Y|.                 (7)
```

For a balanced `r`-block row, flipping exactly `h` blocks while retaining
its original `K_i` gives `q=exp(hw+o(w))`. If `1<=h<r`, (7) yields
`|Im V|>=exp((r-h)w/2-o(w))`. Thus every proper nonempty raw flip fails
the required subpower imaginary height. Unlike (1), this conclusion
does not allow a new Gaussian multiplier on the flipped product.

That distinction is necessary even for primitive witnesses and disjoint
blocks. Let

```text
u_j+b_j sqrt(5)=(2+sqrt(5))(9+4sqrt(5))^j,       j>=1,
u_j^2-5b_j^2=-1.
```

Suppress the index and set

```text
F=(u+2b)+ib,            A=b+i(2b-u),
q=N(F)=10b^2+4ub-1,     n=N(A)=10b^2-4ub-1.
```

Direct multiplication gives the exact identities

```text
U=FA=2ub+i,
V=bar(F)A=4b^2+i(1-2b^2),
(2+i)V=(10b^2-1)+2i.                            (8)
```

The recurrence `(u,b)->(9u+20b,4u+9b)` preserves `u` even, `b` odd,
and `gcd(u,b)=1`. Both blocks have coprime coordinates of opposite parity,
and so do both thin witnesses in (8). All four are therefore
conjugate-primitive. The two blocks are nonunits for `j>=1`.

Their norm supports are disjoint. A common prime divisor of `q,n` is
odd and divides `q-n=8ub`. It cannot divide `b`, since then `q=-1`
modulo that prime. If it divides `u`, the Pell equation gives `5b^2=1`
and hence `q=1`, another contradiction. Thus `gcd(q,n)=1`. Moreover
`q/n -> 9+4sqrt(5)`, so `log q=log n+O(1)`.

This is an infinite balanced two-block example: a half flip has raw
imaginary part `1-2b^2`, but the fixed multiplier `2+i` makes its
imaginary part exactly two at only a constant cost in modulus.
The correction may overlap a core prime, as allowed throughout this
note. It is not a full many-row profile and does not contradict (1).
It rules out extending (7) to arbitrary subpower Gaussian multipliers.

## The first signed relation not excluded by singleton divisibility

Consider an exact relation `sum_j c_j B_j=0` among distinct signed
products on the same disjoint blocks, with nonzero Gaussian integer
coefficients. If, on a particular block, only term `j` chooses one
orientation, the opposite orientation divides every other term and is
coprime to `B_j`. It therefore divides `c_j`. Its modulus is
`exp(w/2+o(w))` for a balanced block.

Every nonconstant coordinate of two or three signs has a singleton
side. Thus no relation among at most three distinct signed products
can have all coefficients of size `exp(o(w))`. More generally, if a
target signed product divides a combination of two distinct inputs with
both coefficients nonzero and subpower, and the resulting nonzero
witness has endpoint modulus,
the quotient is also subpower; moving it to the other side gives the
same forbidden three-term relation. If the target equals one input,
the resulting two-term relation gives the exclusion directly.

For four terms, the first remaining possibility has every changing
block split two against two. Orient the first term positively and
group blocks by the three partitions of the four terms. After removing
the common signed factor, the remaining equation has the form

```text
c0 ABC+c1 A bar(B) bar(C)
  +c2 bar(A) B bar(C)+c3 bar(A) bar(B) C=0.       (9)
```

Some of `A,B,C` may be units. Indeed the Pell family (8) realizes a
four-term circuit when only two of the groups are nonunits:

```text
2FA+(2-i)F bar(A)-(2+i)bar(F)A-2bar(F)bar(A)=0.   (10)
```

This is just `2(U-bar(U))-(y-bar(y))=0`, since `Im U=1` and
`Im y=2`. Hence the threshold of three terms is sharp; a universal
coefficient-height gap for every instance of (9) is false. An explicit
family with all three groups large is now proved in
[the four-term circuit note](four_term_signed_circuit_height.md),
which also proves a height-product bound for two independent relations.
There is still no reduction of the full many-row problem to two
sufficiently small signed circuits on the same support.
For combinations of actual rows with differing core supports, the
singleton argument must be applied to their actual factors; it does
not automatically give this same signed-product model.

The checker additionally verifies (7) without floating-point square
roots, thirty-two primitive Pell half-flip and four-term fixtures in
(8) and (10), and the
two-, three- and four-term sign-incidence assertions on a three-block
cube. The separate checker linked from the four-term note covers
the three-large-group family.
