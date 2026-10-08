# Sparse reconstruction from Gaussian ideal windows

## Status

The exact small-value equations give an injective reconstruction
principle once enough oriented core factors are fixed. The key
fact is elementary: a Gaussian ideal cannot contain two primitive
directions in a sufficiently narrow rectangle. It applies across
the entire allowed coefficient and residue window, rather than
fixing a target residue first.

For eight extracted rows, 49 of the 127 nonconstant cut blocks
suffice to determine the labelled **primitive circle configuration**
up to a common Gaussian unit, within an explicit common height
window. This 49-block family is optimal for the stated unweighted
strict whole-cut anchor-majority criterion. A 63-block Fano family also determines
every primitive pair numerator directly.

This is an injectivity theorem, not a bound on the number of
possible fixed blocks. It does not determine every unfixed block
and correcting factor separately, and it does not prove endpoint
uniformity.

## 1. Exact oriented pair numerators

Use common-unit selected rows and the central truncation from
[endpoint_central_truncation.md](endpoint_central_truncation.md).
Rows are indexed by `0,...,k-1`, the anchor is zero, and
`B=2^(k-1)`. Retain the empty block when defining the total weight.
Write

```text
w'_S=log N(H'_S),
W_core=sum_S w'_S,
w=W_core/B,
s=log D.
```

For this note orient the pair numerator by

```text
h_ij=product_p pi_p^((a_jp-a_ip)_+)
                   bar(pi_p)^((a_ip-a_jp)_+),
h_ij/bar(h_ij)=z_j/z_i,
gcd(h_ij,bar(h_ij))=1.                           (1)
```

The core pair product is

```text
C_ij=product_(j in S, i not in S) H'_S
      product_(i in S, j not in S) bar(H'_S),
```

where row zero is in none of the indexing sets `S`. Central
truncation gives the exact Gaussian factorization

```text
h_ij=K_ij C_ij,
log|K_ij|<=s,
N(K_ij) divides D^2.                            (2)
```

The proof is the same for every pair: at a prime separating their
rounded endpoints the remaining exponent is `2r_p-u_ip-u_jp`,
and at a prime with equal rounded endpoints it is `|u_ip-u_jp|`.
Both lie between zero and `2r_p`.

Let `F` be a collection of frozen cut blocks. Define the fixed
oriented Gaussian divisor

```text
Q_ij=product_(S in F, j in S, i not in S) H'_S
      product_(S in F, i in S, j not in S) bar(H'_S).
```

Then

```text
h_ij=Q_ij V_ij,
V_ij=K_ij product_(unfrozen separating blocks, oriented as above).
```

All `Q_ij` are coprime to their conjugates. Their orientation,
including their actual Gaussian values, is part of the fixed
data; a list of rational norms alone is insufficient.

## 2. One primitive direction in a Gaussian ideal window

Let `P` be a fixed nonzero Gaussian integer, and let `U,T>=1`
be fixed common bounds satisfying

```text
2UT<|P|.                                         (3)
```

There is at most one Gaussian integer `h`, up to sign, satisfying

```text
P divides h,
gcd(h,bar(h))=1,
|h/P|<=U,
0<|Im h|<=T.                                    (4)
```

Indeed, two candidates `h=x+it`, `h'=x'+it'` have determinant

```text
xt'-x't in N(P) Z,
|xt'-x't| <= 2|P|UT < N(P).
```

The determinant is therefore zero. Conjugate primitivity implies
that each integer coordinate pair is primitive, and two collinear
primitive integer vectors differ by sign. The argument even
applies to the larger rectangle `|x|<=|P|U`, `|t|<=T`.

The shared bounds are essential. Fixing `P` while permitting
arbitrarily large multipliers or residues does not give uniqueness.

There is an exact inverse-residue algorithm for the possible
candidate. If `P=a+ib` is conjugate primitive, `Q=N(P)>1`, and
`r=-b a^(-1) mod Q` is chosen in `[0,Q)`, then

```text
P Z[i]={x+it : x=r t mod Q},
r^2=-1 mod Q.
```

Choose the candidate sign with `t>0`, and put `b=(rt-x)/Q`.
The reduced fraction `b/t` is a convergent of `r/Q`; in fact it is
the **last** convergent whose denominator is at most `T`.
Its denominator is `t` itself. The next denominator satisfies

```text
t_next > |P|/U-T,
a_next > |P|/(UT)-2.                            (5)
```

Thus one continued-fraction index, followed by the norm and
primitivity checks, decides whether (4) has a solution. The full
proof, including finite continued fractions and their initial
denominator-one convention, is in
[sparse_convergent_parametrization.md](sparse_convergent_parametrization.md).
The standard Legendre criterion used there is also recorded in
[Waldschmidt's continued-fraction lectures, slides 77--79](https://webusers.imj-prg.fr/~michel.waldschmidt/articles/pdf/WAMS2017Erbil2VI.pdf).

## 3. Reconstruction along a connected graph

Fix the Gaussian blocks in `F`. Choose a connected graph on the
labelled rows. For each of its oriented edges, fix common bounds
`U_ij,T_ij>=1` on the candidate class. Suppose every candidate
completion satisfies

```text
|h_ij/Q_ij|<=U_ij,
0<|Im h_ij|<=T_ij,
2U_ij T_ij<|Q_ij|.                              (6)
```

Section 2 determines each edge numerator up to sign and hence
determines `h_ij/bar(h_ij)=z_j/z_i` exactly. Connectedness then
determines every pair ratio. A star from row zero already suffices;
there is no need to impose (6) for all pairs.

Any two actual labelled circle configurations in this candidate
class differ by a common multiplier `lambda in Q(i)`. Divide each
tuple by its own common Gaussian divisor. For the resulting
primitive tuples, integrality of `lambda z_i` for every `i`,
together with a Gaussian Bezout identity for the `z_i`, implies
`lambda in Z[i]`. The reverse argument makes its inverse integral.
Therefore the primitive tuples differ by a Gaussian unit.

This conclusion concerns the primitive configuration. It does
not uniquely recover all factors `H'_S,K_ij`: a small divisor can
sometimes be transferred from an unfrozen block to incident
corrections without changing any `h_ij`. The common divisor of
the original selected tuple is also not determined by its ratios.

Here is a convenient finite common window. Require

```text
(1-eta)w <= w'_S <= (1+eta)w,
log|K_ij|<=s,
0<|Im h_ij|<=T,                                (7)
```

with the same specified `w>0`, `0<=eta<1`, `s>=0`, `T>=1` for
all candidate completions being compared. Every pair is separated
by `n=B/2` core cuts. If `f_ij` of them are frozen, valid common
cofactor bounds can be chosen as follows, together with the
displayed divisor lower bounds:

```text
log U_ij = s+(n-f_ij)(1+eta)w/2,
log|Q_ij| >= f_ij(1-eta)w/2.
```

Consequently a sufficient condition for (6) is

```text
[2f_ij-n-eta n]w/2 > s+log(2T).                (8)
```

More generally the actual fixed divisor weights may replace the
counting lower bound. The content of (8) is that the frozen
oriented divisor contains strictly more than half of the pair's
logarithmic norm, with room for the residue and correction bounds.
The theorem is always used with an explicit common window such
as (6) or (7), not with unrestricted remaining heights.

## 4. Removing private blocks

The original exact row equations can first be viewed one row at
a time. Put

```text
P_i=product_(S contains i, S!={i}) H'_S,
U_i=K_0i H'_{ {i} },
h_0i=P_i U_i.
```

Holding the other blocks fixed gives the sparse window of
Section 2. Since `P_i` has `B/2-1` block factors while `U_i`
has one block factor and a small correction, the strict exponent
gap holds for `B>4`, or `k>=4`.

One can remove **all** singleton-cut blocks at once. Orient the
singleton block toward its unique row by setting

```text
D_i=H'_{ {i} }  for i>0,
D_0=bar(H'_{ {1,...,k-1} }).
```

Freeze only nonconstant cuts whose two sides both have size at
least two. Then every pair numerator factors as

```text
h_ij=Q_ij (K_ij D_j bar(D_i)),                 (9)
```

where the fixed `Q_ij` has `B/2-2` block factors and the remaining
cofactor has two. The condition in (8) becomes

```text
[(B/2-4)-eta B/2]w/2 > s+log(2T).              (10)
```

For any fixed `k>=5`, or in the growing central reduction, this
holds eventually because `s,log T=o(w)` and `eta=o(1)`.
Thus the nonprivate blocks determine the entire primitive
configuration in the stated window, while every singleton block
is allowed to vary.

## 5. An optimal 49-block anchor design for eight rows

Now take `k=8`, so `B=128` and every pair is separated by 64
blocks. Identify the seven nonanchor rows with `Z/7Z`.
Freeze all 29 subsets of size five, six, or seven. We will add
20 of the 35 subsets of size four.

Consider the two disjoint cyclic Fano systems of triples

```text
L_a={a,a+1,a+3},
L'_a={a,a+1,a+5},        a in Z/7Z.
```

Each is a Steiner triple system: every vertex occurs in three
lines and every pair occurs in one. Their 14 lines are distinct.
Omit the 14 complementary four-subsets, and also omit
`{0,1,2,3}`. This additional four-subset is different from those
14, because its complementary triple is consecutive. Freeze the
remaining 20 four-subsets.

Each vertex occurs in 20 of all 35 four-subsets. It occurs in
eight of the 14 omitted complements, and in the additional omitted
four-subset exactly when its label is `0,1,2,3`. Its retained
four-subset degree is therefore 11 or 12. The large subsets
contribute degree

```text
binom(6,4)+binom(6,5)+binom(6,6)=22.
```

Thus each anchor pair has 33 or 34 frozen separating blocks.
There are 49 frozen blocks in total. Condition (8) holds on the
anchor star whenever

```text
(1-32eta)w > s+log(2T).                         (11)
```

This is automatic for sufficiently large extracted eight-row
families, where `eta=o(1)` and `s,log T=o(w)`.

The count 49 is optimal for this unweighted strict whole-cut anchor-majority
criterion. The largest possible total anchor incidence of any
48 subsets is obtained by taking all 29 subsets of size at
least five and 19 of size four. It is

```text
21*5+7*6+1*7+19*4=230.
```

Strict majority at all seven anchor edges requires at least
`7*33=231` incidences. This rules out 48; it is not a claim
that 49 is optimal for every possible arithmetic reconstruction
method, partial-block freezing, or every nonuniform weight profile.

## 6. A 63-block design covering every pair

For an optional all-pair version, freeze all 35 four-subsets,
all 21 five-subsets, and the seven triples of one Fano system.
There are 63 frozen blocks, fewer than half of the 128 inherited
cut blocks. The separating counts are

| Pair | Four-subsets | Five-subsets | Fano triples | Total |
| --- | ---: | ---: | ---: | ---: |
| Anchor and a nonanchor | 20 | 15 | 3 | 38 |
| Two nonanchors | 20 | 10 | 4 | 34 |

The Fano contribution four in the second row is `3+3-2`, since
each vertex has degree three and the pair occurs in one line.
Every primitive pair numerator is therefore determined directly
whenever

```text
(2-32eta)w > s+log(2T).                         (12)
```

The 49-block construction is sufficient for configuration
injectivity; the 63-block construction additionally treats all
pair numerators by the same frozen-divisor test.

## 7. Roughly half the blocks for growing numbers of rows

For general `k`, distinguish nonanchor row one. Freeze

```text
F={S : 1 in S} union {{j}:2<=j<=k-1}.
```

This uses `B/2+k-2` blocks. The pair `(0,1)` has all `B/2`
separating cuts frozen. A pair incident to exactly one of `0,1`
has `B/4+1`; a pair between two other rows has `B/4+2`.
Condition (8) is therefore implied by

```text
(1-eta B/4)w > s+log(2T).                       (13)
```

With only the original full-pattern error estimate, this already
works for `k=floor(c log_2 M)` with any fixed `c<1/4`.
The reason for this temporary threshold is explicit: the
accumulated error `eta B` must tend to zero when the combinatorial
surplus is only one block.

The full range `c<1/2` is available by controlling a few additional
moments during extraction. Let the original sign functions be
`f_0,...,f_(k-1)`, with Fourier moments `mu_A`. In addition to all
pair moments, control the nonempty moments indexed by

```text
{0,1} symmetric_difference {i,j},   i<j, {i,j}!={0,1}.
```

There are at most `k(k-1)` distinct targeted sets, each of size
two or four. Assume `L>=k>=4`, as holds eventually in this
growing selection. If `L` good common-unit rows are available, the
same two-objective averaging used in central truncation gives
the full-pattern bound there together with

```text
|mu_A|<=zeta,
zeta=2sqrt(|target_sets|/(L-3))=O(k/sqrt(M)).     (14)
```

For clarity, sum the squares of these targeted moments and use
the Bessel expectation bound `2/(L-3)` for each. Adding its
normalized energy to the normalized full-pattern energy selects
a tuple with each energy at most twice its expectation. No
moment of order growing with `k` is added to this second budget.

For a pair other than `(0,1)`, the original signed difference
between its frozen-base weight and its unfrozen-base weight is

```text
-(W/2)(mu_{0,1}-mu_({0,1} symmetric_difference {i,j})).
```

This has magnitude at most `zeta W`. Removing the outer central
layers changes it by at most `2s`. The extra singleton cuts add
at least `2(1-eta)w` to twice the frozen core weight minus the
full core pair weight. Finally, the pair correction in (2) has
logarithmic norm at most `2s`. Therefore, writing
`q_ij=log N(Q_ij)` and `d_ij=log N(h_ij)`,

```text
2q_ij-d_ij
 >= 2(1-eta)w-zeta W-4s.                        (15)
```

The pair `(0,1)` has all its core factors frozen and satisfies
the required inequality with still more room. All pair moments
are included in (14), so the original arc also gives the common
bound `T=exp(zeta W/4)`, since `C/2<1`. Thus a sufficient common
slack condition is

```text
(1-eta)w > (3/4)zeta W+2s+log 2.               (16)
```

In the central growing reduction,
`zeta W/w=O(k2^k/sqrt(M))=o(1)` for every fixed `c<1/2`,
as are `s/w` and `eta`. Equation (16) recovers the full range.
As throughout, configurations are compared using common bounds
on these quantities; the frozen values are not being asserted
to control arbitrarily large remaining heights.

## 8. What remains after this sparsification

The fixed divisors on different edges are not independent.
Whenever they share an oriented Gaussian divisor, their modular
roots agree on its rational norm. Their uniquely selected
continued-fraction candidates must also give compatible ratios
around every cycle and admit the prescribed Gaussian block
factorization. Those are exact cross-row arithmetic conditions.
No theorem presently proves that one of them must fail for all
choices of the frozen blocks.

For a fixed 49-block collection the procedure is concrete: form
the seven anchor ideals, use the common `T` and their cofactor
bounds to compute seven prescribed last convergents, and perform
the primitive and norm checks. Failure of any check excludes
every completion in the window. If all pass, the ratios reconstruct
one primitive configuration up to a Gaussian unit; its remaining
cut profile and correction constraints must still be checked.
This decides each fixed input collection, rather than bounding
the number of input collections.

There is also no direct high-genus conclusion from the original
equations alone. A private singleton block `u_i+i v_i` occurs
only in its own row equation. Once the other blocks are fixed,
that equation has the form

```text
b_i u_i+a_i v_i=t_i.
```

On a nonempty open chart one coordinate can be solved rationally.
Hence the unconstrained algebraic system is rational with many
free shared-block coordinates; eliminating the private variables
alone cannot impose an algebraic equation on all those free
coordinates. The new restriction is the exceptionally small
**integral** representative in each ideal, made explicit by the
last-convergent test. Any Thue, Siegel, or resultant argument
must retain that height condition and the shared divisors.

The available result is consequently a smaller, explicit search
space and a rigorous injectivity theorem. The number and height
of admissible frozen-block choices remain unbounded in the
present argument, so endpoint uniformity is still unproved.

## Verification

The ideal-window and last-convergent arguments were independently
audited in
[sparse_convergent_parametrization.md](sparse_convergent_parametrization.md),
including 2,072 exact finite windows, of which 600 were nonempty.
The 49-block and 63-block designs, their incidence counts, finite
window constants, and primitive-normalization conclusion were
independently checked in
[fano_freezing_audit.md](fano_freezing_audit.md).
These are prose proofs, not new Lean theorems.
