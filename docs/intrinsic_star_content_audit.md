# Intrinsic chord products and small-residue elimination

The signed six-factor family has an exact two-chord identity with a
uniformly positive normalized value after its specified clearing.
The corresponding intrinsic identity retains a quotient of two
different conductor products. On a full eight-row cut profile this
quotient can tend to zero simultaneously for every star. This audit
does not show that the associated imaginary residues can be small.

## 1. Exact content of a two-chord star

For three distinct equal-norm Gaussian integers in one common-unit
factorization, separate the constant prime-allocation layers from the
three nonconstant cut types. This gives

\[
z_0=G\overline{A_0}A_1A_2,\qquad
z_1=GA_0\overline{A_1}A_2,\qquad
z_2=GA_0A_1\overline{A_2}.
\tag{1}
\]

Here `G` is their actual common Gaussian gcd, up to a unit. The
`A_i` are pairwise Gaussian-coprime; each collects the layers isolating
row `i`, with conjugation chosen as in (1). Different `A_i` may share
conjugate prime support. Arbitrary prime-power allocations are allowed:
for a sorted local allocation `(0,a,e)`, the minimum and maximum rows
receive opposite Gaussian orientations with exponents `a` and `e-a`.

Put `s_ij=Im(bar(A_i) A_j)`. These are the nonzero primitive
imaginary pair residues, with signs depending on orientation. Indeed
`bar(A_i) A_j` is coprime to its conjugate, while the Gaussian gcd of
the two corresponding points is `G A_k`, up to a unit. Thus

\[
|z_i-z_j|=2|G A_k|\,|s_{ij}|,
\qquad R=|G A_0A_1A_2|.
\]

In particular the star centered at zero satisfies exactly

\[
\boxed{\frac{|z_0-z_1|\,|z_0-z_2|}{R}
=4|s_{01}s_{02}|\frac{|G|}{|A_0|}.}
\tag{2}
\]

This does not discard the common factor of the three-row subset.
It also identifies the missing input in a proposed direct extension
of the signed-family cancellation: a lower bound on the actual
right-hand side must control both residues and the content quotient.

## 2. Every star has a deficit on the full eight-row profile

Use rows `0,...,7`, with zero as the reference row. Choose a Gaussian
prime `pi_S` above a different odd split rational prime for every nonempty subset
`S subset {1,...,7}`. For a real parameter `w` tending to infinity, put

\[
H_S=\pi_S^{\lfloor w/\log N(\pi_S)\rfloor},\qquad
\log N(H_S)=w+O(1).
\]

The constants here involve only the fixed 127 chosen primes. Set

\[
L=\prod_{S\ne\varnothing}H_S,\quad
B_i=\prod_{S\ni i}H_S,\quad
z_0=\overline L,\quad
z_i=\overline L\,B_i/\overline{B_i}.
\tag{3}
\]

For sufficiently large `w` these are eight distinct Gaussian integers,
are Gaussian-primitive, and have common radius `|L|`. These assertions
follow directly from the independent disjoint prime supports: every
prime orientation occurs at some row and is absent at another, and
the singleton cuts distinguish all rows. No short-arc assertion is
made about (3).

Fix any of the 168 stars, namely a center and two other rows. Among
all 128 subsets including the empty one, exactly 32 cuts are constant
on that triple, and exactly 32 isolate its center. Removing the empty
cut removes one constant cut and no center-isolating cut. Consequently
the exact factors in (1) satisfy

\[
\log|G|=\frac{31w}{2}+O(1),\qquad
\log|A_{\rm center}|=\frac{32w}{2}+O(1),
\]

and hence, simultaneously for every star,

\[
\boxed{\frac{|G|}{|A_{\rm center}|}
=\exp(-w/2+O(1)).}
\tag{4}
\]

The ambient conductor has `W=127w+O(1)` and balanced-layer mass
`W_4=35w+O(1)`. Thus (4) is a failure of a proposed uniform content
quotient bound even for actual primitive full-cut Gaussian data with
growing balanced mass. It is not an endpoint counterexample: the
actual residues in (2) have not been shown to be small.

## 3. Multiplying stars does not cancel this deficit

Take a multiset of `s` stars and let `n_ij` be its edge multiplicities,
so `sum n_ij=2s`. At a cut `S`, the exponent of its block modulus in
the product of the content quotients in (2) is

\[
s-\sum_{i<j}n_{ij}\mathbf1_{\{i,j\}\text{ crosses }S}.
\]

Each edge crosses exactly 64 of the 127 nonempty cuts. Summing this
expression therefore gives

\[
127s-64(2s)=-s.
\]

For the actual family (3), every fixed such product of star-content
quotients is consequently `exp(-sw/2+O(s))`. Selecting stars, or
multiplying several of them, cannot by itself produce the signed
family's radius-free cancellation on this profile.

The precise unresolved arithmetic is simultaneous control of the
actual imaginary residues in (2). Neither the construction (3) nor
the cut count supplies that control. The audit rules out only the
content-cancellation shortcut; it does not exclude a stronger identity
that uses those residues and their phases together.

## 4. Eliminating the contents gives a seven-dimensional variety

There is a related exact test of whether all 28 small pair residues
automatically give a low-degree curve with small coefficients. Fix row
zero and write its seven conjugate-primitive numerators as
`h_i=x_i+i t_i`. For each nonanchor pair define

\[
G_{ij}=N(\gcd_{\mathbb Z[i]}(h_i,h_j)),\qquad
h_{ij}=h_j\overline{h_i}/G_{ij}=x_{ij}+it_{ij}.
\]

The exact coordinate equations, with their full normalized content,
are established in
[the pair-defect note](maximal_order_pair_defects.md):

\[
G_{ij}x_{ij}=x_ix_j+t_it_j,\qquad
G_{ij}t_{ij}=x_it_j-x_jt_i.
\tag{5}
\]

Fix the 28 signed nonzero integers `t_i,t_ij`, of maximum absolute
value `T`. Eliminating the `G_ij` gives the 21 quadrics

\[
\boxed{x_{ij}(x_it_j-x_jt_i)=t_{ij}(x_ix_j+t_it_j).}
\tag{6}
\]

Their coefficient heights are `O(log(2T))`: their integer coefficients
have absolute value at most `T^3`. On the open set with all
`x_it_j-x_jt_i!=0`, however, this variety is exactly the graph

\[
x_{ij}=t_{ij}\frac{x_ix_j+t_it_j}{x_it_j-x_jt_i}.
\tag{7}
\]

Projection onto the seven algebraically free `x_i` is a birational
isomorphism to an open subset of affine seven-space. All triangle
phase equations hold identically after (7); for example

\[
t_{ij}x_{jk}x_{ik}+t_{jk}x_{ij}x_{ik}
-t_{ik}x_{ij}x_{jk}+t_{ij}t_{jk}t_{ik}=0.
\]

The quotients `x_ij/t_ij` in (7) are the cotangents of differences of
angles with rational tangents `t_i/x_i`. Thus adding these phase equations does
not turn (6) into a curve. This is only a statement about its open
rational compatibility variety, not about integer realizability of an
arbitrary fixed residue tuple.

For actual primitive data one must restore, simultaneously,

\[
(x_it_j-x_jt_i)\mid t_{ij}(x_ix_j+t_it_j),\qquad
G_{ij}=\frac{x_it_j-x_jt_i}{t_{ij}}\in\mathbb Z_{>0},
\]

as well as the exact condition
`G_ij=N(gcd_G(h_i,h_j))` and conjugate-primitivity of every numerator.
Arbitrary rational points of (6) satisfy none of these arithmetic
requirements automatically.

In the actual full-cut construction (3), the two anchor numerators
`B_i,B_j` share exactly 32 blocks. Therefore

\[
\log G_{ij}=32w+O(1),\qquad
\log R=127w/2+O(1).
\tag{8}
\]

These contents have size `R^(64/127+o(1))`. If the desired residues
satisfy `log T=o(w)`, a coefficient containing `G_ij` cannot be counted
as having height `O(log T)`. Keeping it as a variable retains the
arithmetic equations; fixing it as a coefficient imports the large
conductor height. As before, (3) does not establish small residues.

Consequently the small coefficient heights in (6) do not extract the
curve or fixed-formula hypotheses of the existing moving-height
exclusions. The missing height restriction must use the simultaneous
integer quotient and actual Gaussian-gcd conditions. Neither this
dimension calculation nor the star-content count rules out such a
stronger arithmetic step.
