# Attack on a strict finite local-test margin

## Status

The uniform endpoint bound and the proposed balanced-layer bonus remain
unproved. This continuation gives two rigorous refinements:

1. a capped local test that needs an inequality only for separately
   normalized eight-point clusters, eliminating the strongest inherited
   common-factor requirement of the earlier test;
2. an exact obstruction to obtaining a strict binomial-mean margin from
   nonnegative combinations of the existing subset determinant inequalities.

No genuine endpoint counterexample to either proposed arithmetic bonus
has been found here. The scaling calculations below identify differences
between the proposals; they do not claim that a counterexample exists.

Throughout, the endpoint constant is \(C_0=1/2\). Prime-exponent layers
have weight \(\log p\), so a primitive circle of radius \(R\) has total
weight \(W=\log(R^2)\).

Subsequent exact tests are recorded in
[the Pfaffian content analysis](eight_point_pfaffian_content.md),
[the reciprocal-chord analysis](reciprocal_chord_content.md), and
[the sum-kernel denominator analysis](sum_kernel_denominator_budget.md).
In particular, the first proves `B_X | T` and shows that the still
unproved endpoint estimate `2log B_X>=W_4-B` would suffice for (1).
Its near-equality audit does not establish that estimate. The positive
sum-kernel remainders also retain a large termwise denominator; a
smaller reduced denominator of the final remainder has not been proved.

## 1. A weaker target with automatic lifting

For a primitive eight-point cluster, let
\[
D=\sum_a w_a(r_a-4)^2,\qquad
W_4=\sum_{a:r_a=4}w_a.
\]
Primitivity means division by the common Gaussian gcd has already been
performed. Consequently all its retained layer sizes satisfy
\(1\leq r_a\leq7\).

Consider the following **unproved** intrinsic arithmetic statement:

> There exist \(B\geq0\) and \(W_0\geq0\), independent of all conductor
> primes and exponents, such that every primitive eight-point cluster
> on an arc of length at most \(\frac12\sqrt R\), with \(W\geq W_0\),
> satisfies
> \[
> D+W_4\leq2W+B.
> \tag{1}
> \]

It is enough to establish (1) when the eight points belong to one
Gaussian-unit class.

Use the symmetric test
\[
\begin{array}{c|rrrrrrrrr}
r&0&1&2&3&4&5&6&7&8\\ \hline
f(r)&2&9&4&1&1&1&4&9&2.
\end{array}
\tag{2}
\]
Thus \(f(r)=(r-4)^2+\mathbf1_{\{r=4\}}\) for nonconstant layers,
but \(f(0)=f(8)=2\). Its oscillation is \(8\), and
\[
\mu_f
=\frac{2\cdot2+16\cdot9+56\cdot4+112+70}{256}
=\frac{277}{128}
=2+\frac{21}{128}.
\tag{3}
\]

Now take eight points from a larger primitive endpoint cluster. Let \(W\)
be the inherited ambient conductor weight, let \(A\) be the total weight
of the inherited layers constant on these eight rows, and let
\(W'=W-A\) be their conductor weight after dividing their new common
Gaussian gcd. That division preserves cardinality and decreases the
endpoint constant. Its remaining layer sizes are unchanged.

If (1) applies to the normalized eight points, then
\[
\sum_a w_a f(r_a)
=D'+W'_4+2A
\leq2W'+B+2A=2W+B.
\tag{4}
\]
Thus the capped endpoint values in (2) give exactly the inherited-layer
inequality required by the sampling criterion in Appendix D of
docs/multipoint_continuation.md.

There is no hidden height-threshold problem in this lifting. The ordinary
eight-point determinant on the original circle gives
\[
16A\leq D'+16A\leq2W-56\log2.
\]
Consequently
\[
W'\geq\frac78W+\frac72\log2.
\tag{5}
\]
In particular \(W\geq8W_0/7\) is sufficient to ensure \(W'\geq W_0\)
for every eight-point subset.

For an even ambient cluster of \(M\geq8\) points, the proved sampling
estimate is
\[
\left|\mathbb E_T\sum_a w_a f(r_{a,T})-\mu_fW\right|
\leq\frac{42\,\operatorname{osc}(f)}M W.
\]
Combining (3)--(5) with that estimate and the bound
\(W\geq4(M-1)\log2\) gives the conditional cardinality bound
\[
\boxed{M\leq2048+\frac{64B}{21\log2}.}
\tag{6}
\]
If \(W<8W_0/7\), the latter determinant bound instead gives
\[
M\leq1+\frac{2W_0}{7\log2}.
\tag{7}
\]
Odd cardinalities cost at most one point. A restriction to a largest
Gaussian-unit class costs at most a factor four. A general fixed arc
constant \(C\) costs at most a further factor \(\lceil2C\rceil\), by
partitioning the original arc.

This proves that (1) is a sufficient arithmetic target. It does not
prove (1).

### The exact weighted content

For a primitive eight-point cluster, denote by \(V_j\) the total weight
of layers with \(\min(r,8-r)=j\), for \(j=1,2,3,4\). Then (1) is
precisely
\[
7V_1+2V_2\leq V_3+V_4+B.
\tag{8}
\]
The existing determinant only gives
\[
7V_1+2V_2\leq V_3+2V_4-56\log2.
\tag{9}
\]
The missing improvement is therefore explicit: remove one copy of the
balanced-layer mass \(V_4\) from the right side, with a uniform additive
cost.

### A source-gcd formulation

There is an equivalent formulation using only literal source gcds. Put
`N=R^2`, `H_i=gcd_G(z_k : k!=i)`, and
`L_ij=gcd_G(z_k : k not in {i,j})/(H_i H_j)`, with Gaussian gcds
defined up to units. Primitivity makes every `L_ij` a Gaussian integer:
at a Gaussian prime at most one `H_i` has positive valuation, and its
valuation is no larger than that of the numerator whenever its index
is deleted. If `A=Norm(product_i H_i)` and
`F=Norm(product_(i<j) L_ij)`, threshold-layer counting gives exactly
`log A=V_1` and `log F=V_2`. A layer whose zero side has one row
contributes to its `H_i`; a layer whose zero side has two rows
contributes to their `L_ij`; all larger zero sides contribute neither.
Counting both conjugate primes includes cuts of sizes `1,7` and `2,6`,
respectively, with their full multiplicities. The nested zero sides
also show that the `H_i` are pairwise Gaussian-coprime, the `L_ij`
are pairwise Gaussian-coprime over distinct unordered pairs, and
`gcd_G(H_i,L_jk)=1` when `i` is not in `{j,k}`; overlap between
`H_i` and `L_ij` is allowed. More strongly, the edges `ij` for which
a fixed rational prime divides `Norm(L_ij)` form a matching of at
most two edges. Each Gaussian orientation contributes to at most
one edge, and the two orientations select the two lowest and the
two highest allocation rows, respectively; those pairs are disjoint
whenever their defining gaps are positive, including tied allocations.
In particular adjacent edge norms are coprime. Thus (8) is precisely
`A^8 F^3 <= exp(B) N`. Equivalently, if
`G_2=Norm(product_(i<j) gcd_G(z_k : k not in {i,j}))`, then
`G_2=A^7 F` and the target is `G_2^3 <= exp(B) N A^13`.
These equalities do not supply that inequality: its unresolved input
is still the endpoint arc hypothesis. The existing global matching
and difference-content bounds concern different gcds and cannot be
substituted for these source deletion contents.

The [global content checker](check_global_matching_triangle_content_normalization.py)
also checks this source-gcd dictionary on actual primitive eight-point
tuples, including repeated prime layers and independent row units.

## 2. Why the raw inherited bonus is stronger

Suppose eight points are obtained by multiplying a primitive eight-point
configuration by a Gaussian integer \(d\) whose norm has only split
prime factors. Put
\[
A=2\log|d|.
\]
The new layers contributed by this common factor are constant, so
\[
W_{\rm new}=W_{\rm old}+A,\qquad
D_{\rm new}=D_{\rm old}+16A,\qquad
(W_4)_{\rm new}=(W_4)_{\rm old}.
\tag{10}
\]
These identities remain valid when the primes of \(d\) overlap the
original conductor: the additional initial and final threshold layers
are exactly the constant layers.

The excess in the uncapped proposal transforms by
\[
(D+\delta W_4-2W)_{\rm new}
=(D+\delta W_4-2W)_{\rm old}+14A.
\tag{11}
\]
Thus a primitive inequality with an additive constant independent of
the improved arc scale does not automatically prove the uncapped
inherited inequality.

The geometric scale is important. Multiplication by \(d\) multiplies
the endpoint arc constant by \(\sqrt{|d|}\). Therefore a fixed
configuration cannot simply be scaled by arbitrarily large \(d\)
while remaining on a \(C_0=1/2\) arc. Moreover, an arbitrary commonly
scaled eight-point set need not occur inside a larger primitive
endpoint cluster. Equation (11) is not by itself a counterexample
to the inherited proposal.

For the capped test (2), the corresponding excess is exactly invariant:
the constant layers add \(2A\) to both sides of the lifted inequality.
This is why (1) is a more natural intrinsic target than the uncapped
inherited version.

## 3. The full cone of symmetrized subset determinant tests

For a \(j\)-point subset, \(2\leq j\leq8\), define its cut defect
\[
q_j(u)=\left(u-\frac j2\right)^2-\frac{\eta_j}{4},
\qquad
\eta_j=\begin{cases}0,&j\text{ even},\\1,&j\text{ odd}.\end{cases}
\]
The corresponding endpoint coefficient is
\[
a_j=\frac{j-\eta_j}{4}.
\tag{12}
\]
This includes the odd formula
\((u-s)(u-s-1)\) when \(j=2s+1\).

Fix an eight-row cut with \(r\) marked rows. Average \(q_j\) over all
\(j\)-element subsets of those eight rows, and call the resulting
function \(F_j(r)\). The hypergeometric mean and variance give exactly
\[
\boxed{
F_j(r)=
\frac{j(8-j)}{28}-\frac{\eta_j}{4}
+\frac{j(j-1)}{56}(r-4)^2.
}
\tag{13}
\]
Indeed, if \(p=r/8\) and \(U\sim\operatorname{Hyp}(8,r,j)\), then
\[
\mathbb E(U-j/2)^2
=jp(1-p)\frac{8-j}{7}+j^2(p-\tfrac12)^2,
\]
which simplifies to (13).

Every symmetrized subset determinant test is therefore affine in
\((r-4)^2\). Under the binomial law \(r\sim\operatorname{Bin}(8,1/2)\),
its mean is
\[
\boxed{\mathbb E F_j(r)=a_j.}
\tag{14}
\]
The averaged \(j\)-point determinant inequality is
\[
\sum_a w_a F_j(r_a)
\leq a_jW+j(j-1)\log C_0.
\tag{15}
\]
Thus each such test has precisely zero strict binomial-mean margin.

More generally, suppose one forms nonnegative combinations of these
inequalities and pointwise weakens the left side. In particular, suppose
\[
f(r)\leq c+\sum_{j=2}^8\lambda_jF_j(r),
\qquad \lambda_j\geq0.
\]
The resulting coefficient is
\(a=c+\sum_j\lambda_j a_j\). Taking binomial means and using (14)
gives
\[
\mu_f\leq c+\sum_j\lambda_j\mathbb E F_j=a.
\tag{16}
\]
No such construction can prove the strict condition \(\mu_f>a\).
Favorable additive constants do not change this coefficient statement.

This is a precise limitation of this cone of inequalities. It does
not rule out a new multipoint identity involving phases, residues,
or arithmetic information absent from the cut-defect polynomials.

## 4. Exact balance is not the critical sampled distribution

Normalize the cut weights to a probability distribution \(\nu_r\) on
\(\{0,\ldots,8\}\). The existing determinant controls
\[
\sum_r\nu_r(r-4)^2\leq2+O(1/W).
\]
The binomial distribution lies on the limiting boundary, since its
quadratic mean is two. In contrast, the exactly balanced distribution
\(\nu=\delta_4\) has quadratic mean zero.

The new rational-ray argument excludes this latter, very special
arithmetic profile in the relevant unit model. Removing \(\delta_4\)
from the cut-distribution feasible set does not remove the binomial
boundary point. In fact the binomial mass at four is only \(35/128\),
and the total variation distance from \(\delta_4\) is \(93/128\).
Excluding just a small neighborhood of exact balance would still not
reach that boundary point.

This explains why the exact-balance phase theorem does not turn into
(1) by a direct convexity or averaging argument. The weighted arithmetic
statement must control a broader family of mixed layer profiles.

## 5. Keeping actual chord lengths: the remaining coupling requirement

Writing \(\ell_{ij}=|z_i-z_j|\), the exact eight-point determinant
retains an archimedean upper bound of the form
\[
D_T\leq2W+
2\sum_{i<j}\log\left(\frac{\ell_{ij}}{\sqrt R}\right).
\tag{17}
\]
Replacing all normalized chords by \(C_0\) gives the familiar
\(56\log C_0\). Thus actual chord lengths can improve the additive
term, especially in a very short consecutive subcluster.

However, choosing a subcluster by its geometric span changes the row
sampling law. The hypergeometric/binomial comparison in the finite-test
criterion applies to uniform row subsets; it does not automatically
apply to consecutive subsets or to subsets weighted by their chord
lengths. A successful version of this idea needs a proved joint
estimate coupling conductor cuts to the spatial selection rule.

For a fixed normalized eight-point shape, the sum in (17) is a fixed
constant; it supplies no coefficient proportional to \(W\). To obtain
such a coefficient from (17) alone requires normalized chords shrinking
by a power of the radius, together with control of the cut statistics
on the selected points. Neither requirement has been established here.

The current sharper finite target is therefore (1), or equivalently
(8), with the exact conditional consequence (6). It remains an open
arithmetic step, not a proved bound.

## 6. A sharp cut-and-geometry relaxation for every subset at once

The sampling issue has an explicit obstruction stronger than the cone
calculation in Section 3. The following construction supplies balanced
weighted cuts and actual ordered points on a real circle satisfying
**all** subset determinant inequalities, including their actual chord
lengths, for arbitrarily large cardinalities. It does not supply Gaussian
integer points with those cuts.

Fix an even \(M\geq8\), \(0<C\leq1/2\), and a total weight satisfying
\[
W\geq4(M-1)\log\frac{4(M-1)}C.
\tag{18}
\]
Give equal total weight \(W\) to the complete family of \(M/2\)-subsets
of the \(M\) row labels. One may equivalently use one representative
from each complementary pair, since every statistic below is symmetric
under complementation.

For a fixed \(k\)-row subset \(K\), put
\[
D_K=\sum_a w_a
\left((|K\cap S_a|-k/2)^2-\frac{\eta_k}{4}\right),
\qquad
a_k=\frac{k-\eta_k}{4},
\]
where \(\eta_k\) is zero or one according as \(k\) is even or odd.
Symmetry and the hypergeometric variance give the exact formula
\[
\boxed{
D_K=a_kW-\frac{k(k-1)}{4(M-1)}W.
}
\tag{19}
\]
It is the same for every \(K\) of that cardinality, so selecting
consecutive rows, or selecting by a different geometric rule, cannot
change its value in this family.

Set \(R=e^{W/2}\), and place the rows at the equally spaced angles
\[
\theta_i=\frac{iC}{(M-1)\sqrt R},
\qquad z_i=Re^{i\theta_i},
\qquad 0\leq i<M.
\tag{20}
\]
These are distinct points on an actual real circle and their containing
arc has length \(C\sqrt R\). For any distinct \(i,j\), the elementary
bound \(\sin x\geq x/2\) in the present range gives
\[
\frac{|z_i-z_j|}{\sqrt R}
\geq\frac{C|i-j|}{2(M-1)}
\geq\frac{C}{2(M-1)}.
\tag{21}
\]

Consequently every \(k\)-subset simultaneously satisfies the exact
chord-retaining determinant inequality
\[
\boxed{
D_K\leq a_kW+
2\sum_{\{i,j\}\subset K}
\log\frac{|z_i-z_j|}{\sqrt R}.
}
\tag{22}
\]
Indeed, by (19), the gap from \(a_kW\) is
\(k(k-1)W/[4(M-1)]\). The lower bound (21) makes the sum on the
right at least \(k(k-1)\log(C/[2(M-1)])\). Inequality (18) is
stronger than the sufficient bound obtained by comparing these two
quantities.

There is also a matching positive real cofactor model. The weighted
cut distance of every row pair is
\[
d_{ij}=\frac{M}{2(M-1)}W.
\]
Define
\[
\tau_{ij}
=\frac{|z_i-z_j|}{2R}e^{d_{ij}/2}.
\tag{23}
\]
Then the exact chord-size identity
\[
|z_i-z_j|^2
=4\tau_{ij}^2R^2e^{-d_{ij}}
\]
holds, and (18)--(21) give
\[
\tau_{ij}
\geq\frac{C}{4(M-1)}
e^{W/[4(M-1)]}\geq1.
\tag{24}
\]
Thus the individual nonzero-integer lower bound, considered only as the
inequality \(\tau_{ij}\geq1\), does not break this relaxation.

For any fixed eight-row subset \(T\), let
\(\mathcal T_T=\prod_{\{i,j\}\subset T}\tau_{ij}\). As \(W\to\infty\)
with \(M,C,T\) fixed, the exact sine expression in (20) gives
\[
2\log\mathcal T_T
=\frac{14}{M-1}W+O_{C,M,T}(1).
\tag{25}
\]
Its balanced-layer mass is exactly
\[
W_{T,4}=b_MW,\qquad
b_M=\frac{\binom{M/2}{4}^2}{\binom M8}
\longrightarrow\frac{35}{128}.
\tag{26}
\]
Therefore
\[
\lim_{W\to\infty}
\frac{2\log\mathcal T_T}{W_{T,4}}
=\frac{14}{(M-1)b_M}\longrightarrow0
\quad(M\to\infty).
\tag{27}
\]
For the concrete value \(M=48\), the inner limit is already
\(1763/1771<1\).

This proves a sharp limitation of a specified relaxation: balanced
conductor cuts, true ordered circle geometry, all chord-retaining
subset determinant inequalities, and the individual cofactor lower
bounds can coexist with arbitrarily small residue-growth exponent.
The \(\tau_{ij}\) have not been proved to be integers, and the chosen
angles have not been realized by the Gaussian blocks producing the
cuts. Those simultaneous arithmetic requirements are precisely what
the relaxation omits. It is not a lattice-point counterexample.

## 7. The intrinsic target after normalization: an explicit 64-row audit

The 48-row value above concerns the stronger residue-only bonus. It
does not by itself violate the separately normalized target (1).
The constant-layer subtraction must be retained.

For eight selected rows in the complete balanced `M`-row model, write

```text
a_M = 2 binom(M/2,8)/binom(M,8),
b_M = binom(M/2,4)^2/binom(M,8).
```

Here `a_M W` is the weight constant on the selected rows. Their
inherited defect is `(2-14/(M-1))W`. After removing constant layers,

```text
W' = (1-a_M)W,
D' = (2-14/(M-1)-16 a_M)W,
W'_4 = b_M W.
```

Thus the intrinsic excess is exactly

```text
D' + W'_4 - 2W' = (b_M-14/(M-1)-14 a_M)W.        (28)
```

At `M=48` this coefficient is `-4410/82861`, so it is negative.
At `M=64`, exact binomial arithmetic gives

```text
a_64 = 325/68381,
b_64 = 179800/615429,
D' + W'_4 - 2W' = (232/68381)W = (29/8507)W'.    (29)
```

Consequently the 64-row relaxation exceeds every fixed additive
constant in (1) as its total weight tends to infinity. All eight-row
subsets have these same intrinsic statistics.

### Compatibility with the later relation and phase tests

One can use the compatible short phases from
[the all-integer-combinations construction](integer_combination_separation.md)
instead of the equally spaced phases in Section 6. Choose its total
weight sufficiently large and `C<=1/2`. For every pair, the
integer-combination imaginary-part bound then implies

```text
|z_i-z_j|/sqrt(R) >= 2 exp(-W/(4(M-1))).           (30)
```

Indeed its pair weight is `M W/(2(M-1))`, and the chord is twice
the radius times the corresponding sine. For any `k`-row subset,
(30) gives

```text
2 sum_(i<j) log(|z_i-z_j|/sqrt(R))
 >= k(k-1) log 2 - k(k-1)W/(4(M-1)).
```

Combined with (19), this verifies every chord-retaining subset
determinant inequality again. It also gives all individual imaginary
cofactor lower bounds. Thus compatible factor phases, all integer
combination separation tests, and all subset determinant tests can
hold in this same relaxation, while (29) still violates the desired
intrinsic inequality.

The [full-rank calculation](relation_free_balanced_patterns.md) shows
that every nonzero zero-sum integer combination has nonzero combined
exponents. The all-combinations imaginary-part bound then makes its
reduced monomial nonreal, so its quotient by its conjugate cannot
equal one. Together these facts prove every exact product-rigidity
conclusion for the formal complex points, including `B_32`.
Full exponent rank alone would not exclude accidental relations among
arbitrary complex factors. Common-factor normalization
preserves their distinctness and their affine independence, while
decreasing the endpoint arc constant.

The qualification is unchanged and essential: these formal factors
are not proved to be Gaussian integers. Equation (29) is an explicit
obstruction to this combined relaxation, not a counterexample to (1)
for lattice points. Any completion of the finite-test route must use
arithmetic information beyond the simultaneously verified tests above.
