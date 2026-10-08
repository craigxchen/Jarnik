# Simultaneous product identities from the common jet

This note concerns the punctured orthogonal eight-row template, with
`n=4s+2`, `s>=1`, and positive, pairwise distinct magnitudes `a_j`.
It derives simultaneous algebraic conditions on a solution of the **full**
odd-moment equations. The common jet, pair-difference factorizations, and
distinct middle coefficients were already proved in
[the spectral orientation note](critical_moment_spectral_orientation.md).
The product matrix and three- and four-row identities below are consequences
of those facts. They do not, by themselves, exclude the template or imply a
bound for arbitrary circles.

Let `S` be obtained by deleting a constant column and a balanced column
`b` from eight orthogonal sign rows of length `4s+4`. Normalize the
constant column to `1`; the four entries of each sign in `b` define two
four-row groups. Assume

```
sum_j S_ij a_j^(2r+1) = alpha_r,
          0<=r<s,  0<=i<8.                              (1)
```

Define the row polynomials, their middle coefficients, and pairwise
agreement factors by

```
m=2s+1,       P_i(X)=product_j (X-S_ij a_j),
c_i=[X^m]P_i(X),
D_ik={j:S_ij!=S_kj},
G_ik(X)=product_(j not in D_ik) (X-S_ij a_j).             (2)
```

The cited note proves that the eight `c_i` are pairwise distinct and that
each pair difference factors by its agreement roots. Expanding the remaining
constant or linear factor gives, for rows in different groups,

```
P_i-P_k=(c_i-c_k)G_ik,
c_i-c_k=-2 product_(j in D_ik)(S_ij a_j).                 (3)
```

For rows in the same group,

```
P_i-P_k=(c_i-c_k) X G_ik,
c_i-c_k=-2 [product_(j in D_ik)(S_ij a_j)]
                  [sum_(j in D_ik) 1/(S_ij a_j)].        (4)
```

**Product-matrix condition.** If rows `i,j` are in the first group and rows `k,l` in
the second, the four signed products in (3) obey

```
M_ik-M_il-M_jk+M_jl=0,
M_uv=product_(h in D_uv)(S_uh a_h).                        (5)
```

The complete `4 x 4` cross-group matrix `M` has rank exactly two. All its
`3 x 3` minors vanish, while its `2 x 2` minors are nonzero because

```
det M_({i,j},{k,l})=(c_i-c_j)(c_k-c_l)/4.             (6)
```

**Three-row condition.** For a triple
`i,j,k`, let `T(X)` be the product of `(X-S_ih a_h)` over columns where
all three signs agree. Let `R_i(X)` be the product over columns where
row `i` is the unique dissenting sign, using the root sign shared by
the other two rows. Define `R_j,R_k` similarly. If `i,j` are in one
group and `k` is in the other, then `deg T=s`, `deg R_i=deg R_j=s+1`,
`deg R_k=s`, and

```
(c_i-c_k) R_j
  = (c_i-c_j) X R_k + (c_j-c_k) R_i.                     (7)
```

If all three rows lie in one group, then `deg T=s-1`, all three `R`
polynomials have degree `s+1`, and

```
(c_i-c_k) R_j
  = (c_i-c_j) R_k + (c_j-c_k) R_i.                       (8)
```

Thus the coefficient vectors of the three products in (7), counting
`X R_k` as one product, or in (8) have rank at most two. These are
simultaneous magnitude constraints beyond the separate pairwise
root-order rules.

**Four-row condition.** For `i,j` in one group and `k,l` in the other,
partition the retained columns with a two-versus-two split of their four
signs into three classes:

```
A: S_i=S_k=-S_j=-S_l,
B: S_i=S_l=-S_j=-S_k,
C: S_i=S_j=-S_k=-S_l.
```

Their sizes satisfy `|A|=|B|=|C|+1=:t`, so `t>=1`. With
`F_E(Y)=product_(h in E)(Y-a_h^2)`, put

```
u=(c_i-c_k)(c_j-c_l),
v=(c_i-c_l)(c_j-c_k),
w=(c_i-c_j)(c_k-c_l)=u-v.
```

Then the following polynomial identity holds:

```
u F_A(Y)-v F_B(Y)=w Y F_C(Y).                        (9)
```

In particular, evaluation at zero gives the positive cross-ratio identity

```
u/v = product_(h in B) a_h^2 / product_(h in A) a_h^2 > 0. (10)
```

Write `Q_E(Z)=product_(h in E)(1+Z/a_h^2)`. Dividing (9) by its
constant-term normalization and setting `Y=-Z` gives

```
Q_A(Z)-Q_B(Z)=lambda Z Q_C(Z),
lambda=(w/u) [product_(h in C) a_h^2]
             /[product_(h in A) a_h^2] != 0.          (11)
```

Thus, for every `1<=r<=t`,

```
e_r((a_h^-2)_(h in A))-e_r((a_h^-2)_(h in B))
  =lambda e_(r-1)((a_h^-2)_(h in C)).                 (12)
```

All these differences are nonzero and have the same sign. This gives
a magnitude-sensitive four-row test. No sufficiency is claimed for
these identities.

The identity also determines the entire `C` magnitude set from the
`A,B` magnitudes. Its first positive-degree coefficient fixes

```
lambda=sum_(h in A) a_h^-2 - sum_(h in B) a_h^-2,
Q_C(Z)=[Q_A(Z)-Q_B(Z)]/(lambda Z).                       (13)
```

The quotient must have exactly the `t-1` simple negative roots
`{-a_h^2:h in C}`. For `t=2`, where `C={h_C}`, this reduces to the
explicit equation

```
a_(h_C)^2 =
 [sum_(h in A) a_h^-2 - sum_(h in B) a_h^-2]
 /[product_(h in A) a_h^-2 - product_(h in B) a_h^-2]. (14)
```

For `t>=3`, Newton's inequalities give an additional quick numerical
check on the coefficient differences `d_r` on the left of (12): for
`2<=r<=t-1`,

```
 [d_r / binom(t-1,r-1)]^2
  > [d_(r-1) / binom(t-1,r-2)]
    [d_(r+1) / binom(t-1,r)].                           (15)
```

The inequality is strict because the `C` magnitudes are distinct.

There is also a condition using only the order of the magnitudes and
the order of the `c_i`. Write the column indices in `A,B` in increasing
magnitude order as `h^A_1<...<h^A_t` and `h^B_1<...<h^B_t`. If
`h^A_r<h^B_r` for every `r`, then (10) forces `w/v>0`; if
`h^B_r<h^A_r` for every `r`, it forces `w/v<0`. More generally (10)
gives the strict linear
inequality

```
sign(w/v) [sum_(h in B) log(a_h)-sum_(h in A) log(a_h)] > 0. (16)
```

Because `|A|=|B|`, this test is unchanged by rescaling every
magnitude. A contradictory collection of these linear inequalities
would exclude a proposed ordered sign word without solving the
nonlinear moment equations.

Evaluating (9) at the positive roots of `F_A,F_B,F_C` gives a second
order-only test. Set `epsilon=sign(w/v)` and
`N_E(h)=#{e in E:e>h}` for each class `E` and column position `h`.
Then necessarily

```
h in A: epsilon=(-1)^(1+N_B(h)+N_C(h)),
h in B: epsilon=(-1)^(N_A(h)+N_C(h)),
h in C: N_A(h)+N_B(h) is even.                         (17)
```

For example, the five-column order `A,B,A,B,C` with two `A` columns,
two `B` columns, and one `C` column is impossible: the two `A`
positions demand opposite values of `epsilon`. This conclusion needs
only the strict magnitude order and the four-row identity.
This parity rule is already implicit in the full pair-sign and zero-gap
order tracked by the existing compressed automaton. At an `A` column,
for example, the selected quartet pair is either `{i,k}` or `{j,l}`.
The interval rule places the selected rows together across the zero
gap, so both have the same comparison sign against either unselected
row. The automaton's pair-difference sign formula, with the exact
agreement-root counts, translates those comparison signs into the
first line of (17). The other matching classes give the remaining
lines in the same way. Thus (17) is a convenient extracted check, not
a stronger graph exclusion.

**A two-quartet bound on exact values.** In the even-`s`, `q=2`
orientation, the prescribed coefficient order is
`c_6<c_7<c_0<c_1<c_2<c_3<c_4<c_5`. Let `R_L` be the squared-magnitude
product ratio on the right of (10) for quartet `(1,2,6,7)`, and let
`R_H` be that ratio for `(0,1,4,5)`. Every full moment solution in
this orientation must satisfy

```
(R_L-1)(R_H-1)<1;     in particular min(R_L,R_H)<2.    (19)
```

This is an exact-value condition beyond the signs `R_L,R_H>1`.
It also separates the explicit first-moment and sign witness in
[the quartet log-sign note](critical_moment_quartet_log_sign_feasibility.md)
from the full common-jet equations. That witness shifts every suffix
magnitude by `H >= (2T)^n+T` after a prefix where every quartet has
`h_q>=2`, with original magnitudes at most `T`. Its ratio of products
of magnitudes is at least `H^2/(2^n T^n)>2` for both selected
quartets; the squared-magnitude ratios `R_L,R_H` are larger still.
Thus those particular magnitudes cannot satisfy (19). This does not
exclude other magnitudes on the same exact-distance Euler words.
Indeed, [the multiple-block construction](critical_moment_exact_crossproduct_relaxation.md)
subsequently realizes every scalar cross ratio, together with the first
moment, using different positive real magnitudes and auxiliary coefficients.
The actual middle coefficients and the full polynomial identities remain
additional requirements.

The reduced four-row identity alone does permit distinct positive nodes.
For example, with `A={1,4}`, `B={2,8/3}`, `C={6}` as sets of squared
magnitudes, the exact identity

```
4(Y-1)(Y-4)-3(Y-2)(Y-8/3)=Y(Y-6)
```

holds. The coefficients `u=4,v=3,w=1` can even be realized as the
four differences in (9) by `c_i=0,c_j=1,c_k=2,c_l=3`. This is only a
local algebraic countertest to a generic real-root or polynomial-Mason
exclusion; it is not an eight-row moment solution.

The order test is consistent with one actual admissible automaton path.
In the eleven-mask terminal path displayed in
[the even-`q=2` graph note](critical_moment_compressed_graph_even_q2.md),
the quartet `(i,j,k,l)=(0,2,4,6)` has `A={3,6}`, `B={10,11}`,
and `C={7}` in path order. Its prescribed coefficient order is
`6<7<0<1<2<3<4<5`; hence `w/v>0`, as required by the strict
coordinatewise order of `A` before `B`. This path has eleven columns
and incomplete within-group distance labels, so this check does not
assert a full punctured template or a moment solution. In particular,
the order test has not excluded the complete Euler-word family.

## Proof

The pair distances are `|D_ik|=m` across groups and `m+1` within a
group. The spectral orientation note proves the two factorizations in
(3)--(4), with their coefficients equal to twice the constant or linear
coefficient of the changed-root polynomial. For a monic odd-degree
polynomial with signed roots `z_h=S_ih a_h`, its constant coefficient is
`-product z_h`. For a monic even-degree polynomial, its linear
coefficient is `-[product z_h][sum 1/z_h]`. This gives the explicit
product formulas in (3)--(4).

Equation (3) says `M_ik=(c_k-c_i)/2`, proving (5), the rank bound, and
(6). The rank is exactly two because both groups contain four distinct
`c` values.

For three rows, every column is either an all-agree column or has
exactly one dissenting row. Write `r_i,r_j,r_k` for the counts of the
three dissent types. If `i,j` are together and `k` is in the other
group, the pair distances give

```
r_i+r_j=m+1,   r_i+r_k=r_j+r_k=m,
(r_i,r_j,r_k)=(s+1,s+1,s).
```

The remaining `s` columns are all-agree. The factors in (3)--(4) are
therefore `G_ij=T R_k`, `G_jk=T R_i`, and `G_ik=T R_j`.
Substitute these into `P_i-P_k=(P_i-P_j)+(P_j-P_k)` and cancel `T` to
obtain (7). If all three are together, each pair distance is `m+1`,
so `r_i=r_j=r_k=s+1` and there are `s-1` all-agree columns. The
same substitution and cancellation of `X T` gives (8).

For the four-row assertion, use the polynomial identity

```
(P_i-P_k)(P_j-P_l)-(P_i-P_l)(P_j-P_k)
  =(P_i-P_j)(P_k-P_l).                                  (18)
```

Substitute the pair factorizations (3)--(4). A column whose four signs
all agree contributes the same squared factor `(X-S_ih a_h)^2` to
each term. A column with one dissenting sign contributes the same
single factor `(X-epsilon_h a_h)`, where `epsilon_h` is the majority
sign. Cancel these common factors. In a two-versus-two column, only
the term corresponding to its matching class contains factors, and
they multiply to `X^2-a_h^2`. The right-hand term has the additional
factor `X^2`. Setting `Y=X^2` gives (9).

To count the three classes, let `q` be the number of one-dissent
columns. Summing the two cross distances for pairs `(i,k),(j,l)`
gives `2m=2(|B|+|C|)+q`; the other cross matching gives
`2m=2(|A|+|C|)+q`. Summing the two within-group distances gives
`2(m+1)=2(|A|+|B|)+q`. Solving yields
`|A|=|B|=|C|+1`. Evaluation of (9) at `Y=0` proves (10): the two
left-hand products have the same parity of degree, and all `a_h^2`
are positive. Normalizing those products at zero and substituting
`Y=-Z` proves the stated coefficient-sign property, because every
coefficient of `Z product_(C)(1+Z/a_h^2)` is strictly positive.
The coefficient of `Z` in (11) proves (13), after which the quotient
and its roots prove (14). The normalized Newton inequalities for the
positive numbers `a_h^-2`, `h in C`, give (15). Equality would require
all of them to be equal, contrary to distinct magnitudes.
At a positive root `a_h^2` with `h in A`, equation (9) says
`-v F_B(a_h^2)=w a_h^2 F_C(a_h^2)`. The sign of `F_E(a_h^2)` is
`(-1)^N_E(h)`, giving the first line of (17). Evaluation at a root
in `B` gives the second line, using `sign(w/u)=sign(w/v)` from (10).
At a root in `C`, the positive ratio `u/v` forces `F_A` and `F_B`
to have the same sign, giving the third line.

For (19), write `d_1=c_1-c_0>0`, `d_2=c_2-c_1>0`,
`e_L=c_7-c_6>0`, and `e_H=c_5-c_4>0`. By `R-1=w/v` in (9),

```
R_L-1 = d_2 e_L / [(c_1-c_7)(c_2-c_6)] < d_2/d_1,
R_H-1 = d_1 e_H / [(c_5-c_0)(c_4-c_1)] < d_1/d_2.
```

The first inequality uses `c_1-c_7>d_1` and `c_2-c_6>e_L`;
the second uses `c_5-c_0>e_H` and `c_4-c_1>d_2`.
Multiplication proves (19). The suffix-witness comparison follows
from the stated product-ratio lower bound; no claim about actual
middle coefficients of that first-moment witness is needed.

[The exact checker](check_critical_moment_common_jet_products.py) verifies
the cancellation and class sizes for all 36 cross-group quadruples of a
punctured Walsh-eight fixture using six arbitrary distinct rational
magnitudes. It also verifies the displayed local reduced identity. The
fixture magnitudes do not solve the moment equations, so the checker
does not assert (9) for that fixture. It separately verifies the
eleven-mask path example and its prescribed coefficient order.
It also checks (19) on an independent ordered rational coefficient
fixture, yielding `R_L=9/8`, `R_H=16/15`.
