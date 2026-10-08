# Ordered Ptolemy identities for primitive chord residues

## Status and normalization

This note proves an unconditional residue inequality that retains the
geometric cyclic order of the points. It does not prove the uniform
endpoint bound or the general mixed-layer bonus D11.

Fix Gaussian primes \(\pi_p\) of norm \(p\), one above each rational split
prime \(p\equiv1\pmod4\). Consider distinct points in one common Gaussian
unit class,

\[
z_i=\epsilon\prod_p\pi_p^{a_i(p)}
                  \overline{\pi_p}^{\,e_p-a_i(p)},\qquad
0\leq a_i(p)\leq e_p,\qquad
N=|z_i|^2=R^2=\prod_pp^{e_p}.
\tag{1}
\]

The choices of Gaussian prime representatives are fixed throughout.
The common-unit hypothesis refers to the literal expression (1). A common
Gaussian factor may first be divided out; inert and ramified factors of a
circle are common factors. Inherited split-prime layers, including layers
constant on a selected subset, may be retained in (1).

Write

\[
S_{p,\ell}=\{i:a_i(p)\geq\ell\},\qquad
1\leq\ell\leq e_p,\qquad w_{p,\ell}=\log p,\qquad
W=\log N=\sum_{p,\ell}w_{p,\ell}.
\tag{2}
\]

Thus a prime of exponent \(e_p\) contributes exactly \(e_p\) threshold
layers, with multiplicity. Reorienting \(\pi_p\) exchanges a cut with its
complement and reverses the layer order; all quantities below are invariant
under that operation.

For an unordered pair define

\[
g_{ij}=\prod_p\pi_p^{\min(a_i,a_j)}
                   \overline{\pi_p}^{\,e_p-\max(a_i,a_j)},\qquad
h_{ij}=\prod_{a_i>a_j}\pi_p^{a_i-a_j}
       \prod_{a_i<a_j}\overline{\pi_p}^{\,a_j-a_i}.
\]

Then \(z_i=\epsilon g_{ij}h_{ij}\) and
\(z_j=\epsilon g_{ij}\overline{h_{ij}}\). Put

\[
t_{ij}=|\operatorname{Im}h_{ij}|\in\mathbf Z_{\geq1},\qquad
d_{ij}=\sum_p|a_i(p)-a_j(p)|\log p.
\]

The integer is positive because the points are distinct and their units
agree. Its value is independent of the orientation of the pair. Exactly,

\[
|z_i-z_j|=2t_{ij}|g_{ij}|,\qquad
|g_{ij}|^2=N e^{-d_{ij}}.
\tag{3}
\]

All logarithms in this note are natural.

## 1. A four-point identity with prime-power coefficients

Take four of the points, labeled \(1,2,3,4\) in their geometric cyclic
order. The three perfect matchings are

\[
M_1=(12)(34),\qquad M_2=(13)(24),\qquad M_3=(14)(23).
\]

For a fixed prime \(p\), sort the four valuation values as
\(b_1\leq b_2\leq b_3\leq b_4\), and put \(m_p=b_3-b_2\).
If \(m_p>0\), exactly one matching joins the two smaller valuations
together and the two larger valuations together. Call it \(M_{\kappa(p)}\).
Define positive integers

\[
X_k=\prod_{\substack{p:m_p>0\\\kappa(p)=k}}p^{m_p},
\qquad k=1,2,3.
\tag{4}
\]

The three \(X_k\) are pairwise coprime. Ties within the lower pair or
within the upper pair do not affect the definition; when the middle two
values tie, the prime contributes nothing.

Equivalently, each threshold layer whose restriction to these four
vertices is a two-versus-two cut contributes one factor \(p\) to the
matching that joins the vertices on the same side. There are exactly
\(b_3-b_2\) such layers at \(p\), all giving the same matching. This
equivalence accounts for arbitrary prime powers rather than just binary
valuation rows.

There is a positive real number \(K\), depending on the four points, such
that the product of the two \(|g_{ij}|\) in matching \(M_k\) is \(KX_k\).
Indeed, the within-side matching has valuation-distance sum

\[
(b_2-b_1)+(b_4-b_3),
\]

whereas both other matchings have distance sum
\((b_3-b_1)+(b_4-b_2)\), larger by \(2m_p\).
Equation (3) therefore makes its product of common-factor lengths larger
by precisely \(p^{m_p}\).

Ptolemy's identity for four points on a circle, in this cyclic order, is

\[
|z_1-z_3||z_2-z_4|
=|z_1-z_2||z_3-z_4|+|z_1-z_4||z_2-z_3|.
\]

Substituting (3) and canceling \(4K\) gives the exact positive-integer
identity

\[
\boxed{t_{13}t_{24}X_2
       =t_{12}t_{34}X_1+t_{14}t_{23}X_3.}
\tag{5}
\]

The fact that \(M_2\) consists of the two intersecting diagonals is
essential. The three terms cannot be interchanged freely.

Put \(u_{ij}=\log t_{ij}\geq0\).
Bounding each right-hand summand in (5) by the left-hand side gives

\[
\log X_1+\log X_3-2\log X_2
\leq2(u_{13}+u_{24})
 -(u_{12}+u_{34}+u_{14}+u_{23}).
\tag{6}
\]

In fact, the arithmetic-geometric mean inequality improves the right-hand
side by \(-2\log2\). This improvement will only affect an additive constant
in the final result.

## 2. All seventy quadruples of eight ordered points

Now label eight points \(1,\ldots,8\) in geometric cyclic order and put
\(T=\prod_{i<j}t_{ij}\).
For a cut \(S\subseteq[8]\), write

\[
r(S)=|S|,\qquad
L(S)=\binom{r(S)}2\binom{8-r(S)}2.
\tag{7}
\]

There are exactly \(L(S)\) four-element vertex subsets on which the cut
has two vertices on each side.

Let \(a(S)\) count the alternating quadruples: increasing indices
\(i<j<k<l\) with

\[
S\cap\{i,j,k,l\}=\{i,k\}
\quad\hbox{or}\quad
S\cap\{i,j,k,l\}=\{j,l\}.
\tag{8}
\]

For these quadruples the matching joining equal cut values is the
diagonal matching \(M_2\). It contributes coefficient \(-2\) to the
left side of (6); each other two-versus-two quadruple contributes \(+1\).
Consequently summing the left sides of (6) over all
\(\binom84=70\) quadruples gives exactly

\[
\mathcal Q
:=
\sum_{p,\ell}w_{p,\ell}
  \bigl(L(S_{p,\ell})-3a(S_{p,\ell})\bigr).
\tag{9}
\]

A fixed pair \(i<j\), with \(s=j-i\), is a diagonal in precisely

\[
m_{ij}=(s-1)(7-s)\leq9
\tag{10}
\]

quadruples: choose one remaining vertex from each of the two open cyclic
intervals cut out by the pair. Their sizes are \(s-1\) and \(7-s\),
which sum to six. Dropping the negative residue terms in (6) and using
\(u_{ij}\geq0\) therefore gives

\[
\boxed{\mathcal Q\leq18\log T.}
\tag{11}
\]

There is a useful sharpening. Each pair belongs to
\(\binom62=15\) quadruples in total, so retaining all terms in (6)
gives coefficient \(2m_{ij}-(15-m_{ij})=3m_{ij}-15\leq12\).
Retaining the arithmetic-geometric mean constant as well gives

\[
\boxed{\mathcal Q
\leq\sum_{i<j}(3m_{ij}-15)\log t_{ij}-140\log2
\leq12\log T-140\log2.}
\tag{12}
\]

Both bounds allow constant inherited layers; their \(L\) and \(a\)
are zero.

## 3. A cyclic-interval residue bound

Suppose every nonconstant threshold cut is a cyclic interval in the
geometric ordering of the eight points. Its complement is also a cyclic
interval. Such a cut has no alternating quadruple, so \(a(S)=0\).

Define the balanced-layer weight by

\[
W_4=\sum_{\substack{p,\ell\\|S_{p,\ell}|=4}}\log p.
\]

Every balanced cut contributes

\[
L(S)=\binom42\binom42=36.
\]

All other contributions to (9) are nonnegative. Thus (11) proves

\[
36W_4\leq18\log T,\qquad
\boxed{2\log T\geq4W_4.}
\tag{13}
\]

The sharper inequality (12) yields

\[
\boxed{2\log T\geq6W_4+\frac{70}{3}\log2.}
\tag{14}
\]

These are unconditional arithmetic inequalities within the stated
cyclic-interval class: no small-arc assumption enters their proof.

At the endpoint scale there is already a simpler and stronger point-count
argument for this class. Suppose \(M\geq3\) points of the form (1) lie on
an arc of length at most \(C\sqrt R\). Equation (3) implies for every
distinct pair

\[
d_{ij}\geq W/2+2\log(2/C).
\tag{15}
\]

If every nonconstant cut is a cyclic interval among these \(M\) points,
each such cut separates exactly two cyclic neighboring pairs. Hence

\[
M\bigl(W/2+2\log(2/C)\bigr)
\leq\sum_{\text{cyclic neighbors}}d_{ij}\leq2W.
\tag{16}
\]

For \(C<2\), this forces \(M<4\). In particular at \(C=1/2\) the
inequality is \(M(W/2+4\log2)\leq2W\), and there are at most three
points. Thus (13) and (14) do not establish a previously unresolved
endpoint point count for interval cuts; their additional content is the
unconditioned residue inequality.

For the point-count conclusion alone the common-unit assumption can also
be removed at \(C<1\). For arbitrary units, the primitive Gaussian
numerator of a nonzero chord has modulus at least one, giving
\(d_{ij}\geq W/2-2\log C\) instead of (15). The same neighboring-pair
sum then forces \(M<4\). This does not change the unit hypothesis in the
exact residue identities above.

## 4. Alternating mass and the surviving obstruction

For arbitrary cuts define

\[
\mathcal L=\sum_{p,\ell}w_{p,\ell}L(S_{p,\ell}),\qquad
\mathcal A=\sum_{p,\ell}w_{p,\ell}a(S_{p,\ell}).
\]

The general result is exactly

\[
\boxed{\mathcal L-3\mathcal A\leq18\log T,}
\tag{17}
\]

or the sharper version with \(12\log T-140\log2\) on the right.
Since \(\mathcal L\geq36W_4\), these imply, respectively,

\[
2\log T\geq4W_4-\mathcal A/3,\qquad
2\log T\geq6W_4-\mathcal A/2+\frac{70}{3}\log2.
\tag{18}
\]

A large alternating mass can therefore erase this route to a positive
balanced-layer bonus.

This is not merely a loss caused by a loose multiplicity bound. If the
cut weights are uniform over all subsets of any fixed size \(r\), then

\[
\operatorname{average}a(S)=L(S)/3,
\]

because, conditional on a quadruple splitting two-versus-two, two of its
six assignments are alternating. The weighted expression
\(\mathcal L-3\mathcal A\) is exactly zero for any distribution obtained
by mixing such uniform fixed-size distributions.
In particular the fair independent-bit distribution on eight vertices
has

\[
\mathbf E(r-4)^2=2,\qquad
\Pr(r=4)=70/256,\qquad
\mathbf E\bigl(L(S)-3a(S)\bigr)=0.
\]

These are the critical defect and positive balanced mass in the sampling
obstruction. The present inequalities then give only a fixed additive
residue lower bound, not growth proportional to the conductor weight.
No theorem forcing arbitrary arithmetic threshold cuts to be cyclic
intervals, or forcing an adequate deficit of alternating mass, is proved
here.
