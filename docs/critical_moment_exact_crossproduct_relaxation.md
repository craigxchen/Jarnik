# First moments and all exact cross products can coexist

This note treats a relaxation of the punctured orthogonal eight-row
moment problem. In the even-`s`, `q=2` compressed orientation, there are
exact-distance signed-column words with strictly increasing positive
real magnitudes that satisfy the eight common first-moment equations
and **all 36 exact scalar quartet cross-ratio equations** in
[the common-jet product note](critical_moment_common_jet_products.md),
equation (10), for one prescribed strict coefficient order. The
construction works for every `s=2+20,805,120 K` with integer `K>=33`.

The prescribed coefficients are free auxiliary real numbers. They are
not asserted to equal the actual middle coefficients of the row-root
polynomials. Higher odd moments and the full four-row polynomial
identities remain unresolved. Thus this is a feasibility result for
the first-moment plus scalar-cross-product relaxation, not an actual
full-moment solution or a circle construction.

The later [reciprocal-prefix obstruction](critical_moment_reciprocal_prefix_obstruction.md)
excludes every magnitude tuple produced below, even without imposing the
large-shift threshold. This does not invalidate the relaxation construction: it
proves that these first moments and exact scalar cross products do not
supply the full moment equations. Other magnitudes on the same words
remain unresolved.

## Exact 31-block word

Number the four rows in each group `0,...,3` and `4,...,7`. Fix the
prescribed strict order

```
gamma_6<gamma_7<gamma_0<gamma_1<gamma_2<gamma_3<gamma_4<gamma_5.
```

Let `D_ik` denote a row-pair Hamming distance and take `(0,4)` as
reference among the sixteen cross pairs. The 29 terminal-rooted active
closed walks in [the saved graph data](critical_moment_compressed_graph_even_q2_summary.json)
have full rational label span. More concretely, the
[exact checker](check_critical_moment_exact_crossproducts.py) solves
their 29-by-29 integer label matrix to obtain an **integer** walk
combination for each of the fifteen relative cross-label directions:

```
length=0;  within-pair labels=0;  D_04=0;
D_ik=4 at one other cross pair and 0 at the remaining fourteen. (1)
```

Every coefficient in these fifteen combinations has absolute value at
most six. Let `e_l in Z^15` be the cumulative relative cross-label
vector after block `l`. Prescribe

```
e_0=0,
e_(2j-1)=+4 unit_j,  e_(2j)=-4 unit_j,  1<=j<=15,
e_31=0.                                                     (2)
```

Block `l` uses the integer walk combination with relative label
`e_l-e_(l-1)`. Add the earlier integer `s=2` closed-walk correction
to the last block; that correction has zero cross-label increments.
Each block is a signed terminal-rooted circulation. Add one copy of
the positive active circulation to each block, except blocks 24 and
26, which need two copies. The checker replays every walk and verifies
that these 31 flows are strictly positive on every active edge of the
closed giant component. Their circulation multipliers sum to 33.
Here conservation, positive support, and the isotropic labels of that
active circulation are the previously proved facts in
[the compressed-graph construction](critical_moment_compressed_graph_even_q2.md).
For `K>33`, add the extra `K-33` copies to the last block.

Every block therefore has a terminal-rooted Euler circuit. Concatenate
the published eleven-column initial-to-terminal path and these 31
circuits. The signed labels in (2) telescope, leaving just the earlier
`s=2` correction plus `K` active circulations. The resulting word has
exactly `n=4s+2` columns, cross-pair distance `2s+1`, and within-group
distance `2s+2` for `s=2+20,805,120 K`. Every block ends at the
terminal compressed height

```
delta=(0,0,0,0,-1,-1,1,1).                               (3)
```

Because every block uses every active giant-component edge, the full
word visits all fourteen axial heights `delta±unit_r`, `1<=r<=7`,
at proper prefixes. [The barycentric first-moment construction](critical_moment_first_moment_barycentric.md)
therefore supplies strictly increasing positive integer baseline
magnitudes `a_j<=A=16n` with common first moments.

## Independent shifts preserve the first moments

Keep the eleven initial-path magnitudes fixed. Add a constant `H_l>0`
to every baseline magnitude in block `l`, with
`H_1<...<H_31`. The resulting magnitudes remain strictly increasing.
Each block starts and ends at height `delta`, so its signed column sum
is constant across the eight rows. Adding `H_l` to that block therefore
changes all eight first moments by the same amount. The common
first-moment equations remain exact.

For a cross pair `p=(i,k)`, put

```
M_p=product_(h:S_ih!=S_kh)(S_ih a_h),
F_p=log|M_p|-log|M_(0,4)|,  p!=(0,4).                   (4)
```

Let `v_l=e_l-e_(l-1)` be the relative cross-distance vector of block
`l`. Choose `H_l=R exp(t_l)`, with `t_1=0`, successive gaps
`g_l=t_(l+1)-t_l>0`, and `R>A`. The eleven-column path has equal cross
distances of five, so its contribution to (4) is a fixed vector `b`.
Since `sum_l v_l=e_31-e_0=0`, expanding
`log(a_h+H_l)=log R+t_l+log(1+a_h/H_l)` gives

```
F(R,g)=b-sum_(l=1)^30 e_l g_l+E(R,g),
||E(R,g)||_infinity <= 2nA/R             when all t_l>=0. (5)
```

The error bound counts at most `n` terms for each of the two pairs in
one relative coordinate and uses `0<=log(1+x)<=x` for `x>=0`.

Fix **any** eight real `gamma_i` in the prescribed strict order, and
set for the fifteen nonreference cross pairs

```
d_(i,k)=log|gamma_k-gamma_i|-log|gamma_4-gamma_0|.     (6)
```

In fact any vector `d in R^15` can be realized: there are positive gaps
`g_l` for which `F(R,g)=d` exactly. For coordinate `j`, put

```
L_j=1+max(0,(b_j-d_j)/4),
g_(2j-1)=L_j,
g_(2j)=L_j+(d_j-b_j)/4+u_j,
u_j in [-1/4,1/4].                                      (7)
```

Every gap in (7) is at least `3/4`. By (2), the part of (5) before
the error is `d+4u`. Choose `R>=8nA`; then
`||E(R,g(u))||_infinity<=1/4` uniformly over the cube. The continuous
map `u -> -E(R,g(u))/4` maps that cube into itself. Brouwer's fixed
point theorem gives `u` with `u=-E/4`, so (5) equals `d` exactly.
The same estimates make this constructive. Each `t_l` depends on
any one `u_k` with derivative either zero or one, and

```
|d/dt log(1+a/(R exp(t)))| = a/(R exp(t)+a) <= A/R.
```

Thus every entry of the derivative matrix of `E` has absolute value
at most `2nA/R<=1/4`. The map `u -> -E/4` has supremum-norm
Lipschitz constant at most `15/16<1` on the cube. Iterating it from
`u=0` stays in the cube and converges to its unique fixed point there.
This gives finite positive real shifts and retains strict magnitude
order. It does not claim rational or integer shifted magnitudes.

Finally, the signs of the sixteen `M_ik` match the signs of
`gamma_k-gamma_i`. In this even-`s` orientation, the final height (3)
implies that a cross pair has `2s+1` differing columns and that the
number of negative `S_ij` among them is `(2s+1+delta_k)/2`. It is even
for `k=4,5` and odd for `k=6,7`. The former `gamma_k` lie above all
four first-group values; the latter lie below them. Thus (4)--(6)
give one scalar `rho>0` with

```
M_ik=rho (gamma_k-gamma_i)    for all 16 cross pairs.  (8)
```

Setting *formal* coefficients `c_i=2rho gamma_i` makes
`M_ik=(c_k-c_i)/2`. Consequently the complete cross-product matrix
has rank two and satisfies every exact scalar quartet cross-ratio
equation. The actual coefficients `[X^(2s+1)]P_i(X)` of these row
polynomials need not equal the formal `c_i`. The construction does not
assert the within-group product formulas, three-row identities, or
four-row polynomial identities of the full moment problem.

## Exact audit and scope

The standard-library checker performs rational elimination on the
saved 29-walk labels, replays every walk through the compressed
automaton, checks the 31 cumulative labels, and proves positive
edge multiplicities with total active-circulation multiplier 33.
It also checks that the signed block labels telescope to the old
integer correction and that the resulting overall length has the
required degree. It does not materialize a word with billions of
columns, numerically approximate the fixed point, or claim higher odd moments.
