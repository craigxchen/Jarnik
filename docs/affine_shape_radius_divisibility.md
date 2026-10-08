# Primitive radius from affine triangle data

This is an exact radius bound in terms of affine shape, not a uniform
short-arc point-count theorem. It uses the original Gaussian coordinates
and their integer triangle determinants.

Let `m >= 5` distinct Gaussian integers `z_i` have common norm `N=R^2`
and Gaussian gcd a unit. For each unordered triple `I={i,j,k}`, set

\[
D_I=\det(z_j-z_i,z_k-z_i),\qquad
g=\gcd_I|D_I|,\qquad p_I=D_I/g,
\qquad H=\max_I|p_I|.
\]

All `D_I` are nonzero. Put

\[
K=\binom m3,\qquad
b_m=\binom{\lfloor m/2\rfloor}3+
\binom{\lceil m/2\rceil}3.
\]

Then

\[
\boxed{N^{b_m}\mid\prod_{|I|=3}|p_I|},\qquad
\boxed{R\le H^{K/(2b_m)}}.
\tag{1}
\]

In particular, five or six points satisfy `R <= H^5` with constant
one. Primitive triangle ratios are invariant under nonsingular real
affine transformations whenever both configurations are integral:
every triangle determinant changes by the same nonzero factor. Thus
(1) gives a quantitative fixed-affine-shape radius bound. If the
original Gaussian gcd is `G`, the same statement applies to
`R/|G|`; dividing out `G` does not change the primitive triangle data.

## 1. The local valuation formula

Gaussian primitivity implies that `N` is odd and has no inert prime
factor. Indeed, a prime above two, or an inert rational prime dividing
`N`, would divide every `z_i`.

Fix a split rational prime `p`, with `p^e || N`, and one Gaussian
prime `pi` above it. Put `t_i=v_pi(z_i)`. Primitivity at both `pi` and
its conjugate gives

\[
\min_i t_i=0,\qquad \max_i t_i=e.
\]

The exact identity, up to the orientation sign, is

\[
2iD_{ijk}=
N\frac{(z_j-z_i)(z_k-z_i)(z_j-z_k)}{z_i z_j z_k}.
\]

Since `p` is odd and `D_I` is a rational integer,

\[
v_p(D_{ijk})=
e+\sum_{\{a,b\}\subset I}v_\pi(z_a-z_b)-\sum_{a\in I}t_a.
\tag{2}
\]

For unequal `t_a,t_b`, the difference valuation is their minimum.
For equal levels, write its nonnegative excess above that minimum
as `r_ab`. Consequently

\[
v_p(D_I)=e-\bigl(\max_{i\in I}t_i-\min_{i\in I}t_i\bigr)
+\sum_{\substack{\{a,b\}\subset I\\t_a=t_b}}r_{ab}.
\tag{3}
\]

## 2. Three or more valuation levels

Choose a triple containing a minimum-level point, a maximum-level
point, and a point at a strictly intermediate level. Equation (3)
gives valuation zero. Thus `v_p(g)=0`.

For each integer `1 <= h <= e`, split the rows according to
`t_i >= h`. The first term in (3) counts exactly the levels at which
the triple is contained entirely on one side of this split. A split
of sizes `s,m-s` has

\[
\binom s3+\binom{m-s}3\ge b_m
\]

such triples. Summing (3) over all triples gives

\[
\sum_I v_p(D_I)\ge e b_m.
\tag{4}
\]

## 3. Exactly two valuation levels

The two levels must be zero and `e`; let their sizes be `s,m-s`.
Let `r` be the minimum `r_ab` over all pairs within either level.
Such a pair exists. Every mixed triple has one repeated-level pair,
so its valuation is that pair's excess. Every pure triple has
valuation `e` plus three excesses, hence exceeds `r`. A pair attaining
`r`, together with any point on the other side, gives a mixed triple
of valuation `r`. Therefore `v_p(g)=r` exactly.

Each within-level pair occurs in `m-2` triples. Write

\[
E_s=\binom s2+\binom{m-s}2,\qquad
b_s=\binom s3+\binom{m-s}3.
\]

Equation (3) gives

\[
\begin{aligned}
\sum_I v_p(p_I)
&=e b_s+(m-2)\sum_{t_a=t_b}r_{ab}-K r\\
&\ge e b_s+\bigl((m-2)E_s-K\bigr)r\\
&\ge e b_m.
\end{aligned}
\tag{5}
\]

For completeness, `E_s >= m(m-2)/4`, so

\[
(m-2)E_s-K\ge \frac{m(m-2)(m-4)}{12}\ge0.
\]

For five rows, the `2+3` split gives the sharper expression
`e+2r`; a `1+4` split gives `4e+8r`. These calculations retain the
possible common triangle content instead of assuming it is coprime
to the radius.

Equations (4) and (5), for every prime dividing `N`, prove (1).

## 4. Scope and exponent barrier

The exponents `K/(2b_m)` are five for `m=5,6`, `7/2` for `m=7,8`,
and three for `m=9,10`. They approach two strictly from above.
For four rows the balanced `2+2` split has no pure triple; the same
argument has zero radius exponent. This is consistent with the
fixed affine relation and bounded determinants in the four-point
Pell family of
[the counterexample note](four_point_bonus_counterexample.md).

On a sufficiently short minor arc of length `L`, each actual triangle
determinant satisfies `|D_I| <= L^3/(8R)`: its three chord lengths
are at most `a,b,a+b`, where the successive arc gaps have
`a+b <= L`, and the determinant is their product divided by `2R`.
Since `g>=1`, the same upper bound applies to `H`.
At `L=C sqrt(R)`, this gives only `H <= (C^3/8) sqrt(R)`.
Substitution in (1) leaves a positive power of `R` for every fixed
`m`, because `K/(2b_m)>2`. Thus the reconstruction estimate does not
close the endpoint problem, and bounded affine shape must not be
assumed for arbitrary endpoint clusters.

Retaining `g` gives a further exact restriction. Equations (1) and
the triangle bound imply

\[
R^{2b_m}\le
\left(\frac{C^3\sqrt R}{8g}\right)^K,
\qquad
\boxed{g\le\frac{C^3}{8}R^{\gamma_m}},
\]

where

\[
\gamma_m=\frac12-\frac{2b_m}{K}
=\begin{cases}
\displaystyle\frac{3}{2(m-1)},&m\text{ even},\\[4pt]
\displaystyle\frac{3}{2m},&m\text{ odd}.
\end{cases}
\]

The integer `g` is also the index of the rank-two sublattice generated
by all differences `z_i-z_0`, since its maximal minors are precisely
the anchored triangle determinants and these generate the same
integer ideal as all triangle determinants. Thus any hypothetical
sequence with growing point count has difference-lattice index
`R^{o(1)}`. This restriction still leaves the determinant-method
endpoint exponent unchanged.

The proof is algebraic and does not invoke an external theorem.
An exact integer audit checked 39,840 five-point samples with
`N <= 2000`, including 541 with nontrivial overlap between triangle
content and the primitive norm. The proof, rather than these finite
checks, establishes the stated divisibility. No Lean formalization
is claimed.

The local valuation argument and all parity-dependent exponents passed
an independent review. The [joint relation-lattice formulation](joint_affine_relation_lattice.md)
records the exact complementary-minor identity. A separate
[five-point endpoint family](five_point_affine_shape_cubic_family.md)
has `R` of order `H^3`, so replacing the universal reconstruction
exponent by two is false, even at a fixed endpoint constant.
