# Radius growth under rational Möbius maps

For an anchored short-arc source with at least five points, every
nonconformal rational Möbius map that does not increase its least primitive
radius must have height growing by a fixed power of the source radius.
The statement allows moving sources, arbitrary nested Gaussian prime
powers, and arbitrary nonzero matrix determinant. It requires no fairness
hypothesis. This is a constraint on a proposed descent, not a uniform bound
on the number of source points.

## Source and target normalization

Let \(z_0,\ldots,z_{m-1}\) be a primitive fixed-Gaussian-unit source
cluster of common norm \(N\), on an arc of length at most
\(C N^{1/4}\), where \(0<C\le2\) and \(m\ge3\). Thus the original
circle radius is \(\sqrt N\). The usual common Gaussian gcd and unit-class
reduction gives this setting; extracting one unit class may retain only a
quarter of a general cluster. All point counts below refer to this retained
primitive cluster.

Choose primitive, conjugate-coprime Gaussian half-angle factors
\[
 H_0=1,\qquad z_i/z_0=H_i/\overline H_i,\qquad f_i=|H_i|^2.
\]
Their norms are odd. The \(H_i\) are integer columns when acted on by a
matrix. They need not have equal norms: their *phases* have a primitive
circle realization of squared radius \(N\). For each pair put
\[
 e_{ij}=|\gcd_{\mathbb Z[i]}(H_i,H_j)|^2,\qquad
 b_{ij}=\frac{|\det(H_i,H_j)|}{e_{ij}},\quad
 B_* =\prod_{i<j}b_{ij},\quad E=\binom m2.
\]
These are positive integers. The endpoint chord bound and the arbitrary
prime-power layer-cake estimate give
\[
 f_i\ge(4/C^2)\sqrt N\quad(i>0),\qquad
 B_*\le(C/2)^E N^{m/8}.                                      \tag{1}
\]
See [the pair reduction](kurlberg_wigman_shrinking_arc_audit.md) and
[triangle compatibility](pair_norm_triangle_compatibility.md). In the
product estimate, \(b_{ij}\le(C/2)\sqrt{d_{ij}}N^{-1/4}\) and
\(\sum_{i<j}\log d_{ij}\le m^2\log N/4\); subtracting the angular
contribution leaves the exponent \(m/8\).

Represent a rational Möbius map by a primitive integral matrix
\[
 M=\begin{pmatrix}a&b\\c&d\end{pmatrix},\quad\delta=\det M\ne0,
 \qquad Q=M^TM=\begin{pmatrix}A&B\\B&D\end{pmatrix},
\]
and define
\[
 T=A+D,\qquad K=(A-D)^2+4B^2=T^2-4\delta^2.
\]
For each transformed column remove its ordinary coordinate gcd:
\[
 r_i=\gcd_{\mathbb Z}((MH_i)_1,(MH_i)_2),\qquad
 B_i=MH_i/r_i,\qquad F_i=|B_i|^2.
\]
Let \(N'\) be the least primitive squared radius realizing all phases
\(B_i/\overline B_i\), **including the transformed anchor**.

## Determinant-free pair transfer

For \(K>0\), the precise arithmetic input is
\[
 \gcd(f_i,|MH_i|^2)\mid K,\qquad
 |\gcd_{\mathbb Z[i]}(B_i,B_j)|^2\mid Kb_{ij}.                \tag{2}
\]
The second statement includes primes dividing \(\delta\); ordinary row
normalization removes the determinant loss that appeared in a preliminary
bound. A complete proof is in
[the determinant-free transfer audit](mobius_determinant_free_transfer_audit.md),
independently checked by Astra, Sol, and the root agent.

For the first statement, set \(n=\gcd(f_i,|MH_i|^2)\), write
\(H_i=(x,y)\), and invert \(y\) modulo \(n\). Primitivity ensures
that inverse exists. The quotient \(t=x/y\) obeys
\(t^2=-1\) and \(2Bt+(D-A)=0\). Squaring the latter relation gives
\(K=0\pmod n\). This proof includes all prime powers and does not divide
by two. Consequently the lcm of the row gcds divides \(K\); their product
need not divide \(K\).

The adjugate identity implies \(r_i\mid|\delta|\). The smaller
eigenvalue of \(Q\) is at least \(\delta^2/T\), so
\[
 F_i\ge\frac{\delta^2 f_i}{T r_i^2}\ge\frac{f_i}{T},
 \qquad F_0\le T.                                           \tag{3}
\]
Thus no determinant factor reappears in the real norm comparison.

## Compressing the common-factor loss

At an odd split prime, let \(t_i\) be the signed orientation valuation of
\(B_i\). Integer primitivity means it has positive valuation at at most
one of the two conjugate primes. The exponent of \(p\) in \(N'\) is
\(\max_i t_i-\min_i t_i\), while \(v_p(F_i)=|t_i|\).
Put \(k=v_p(K)\) and \(\beta_{ij}=v_p(b_{ij})\). The pair transfer
implies
\[
 t_i t_j>0\quad\Longrightarrow\quad
 \min(|t_i|,|t_j|)\le k+\beta_{ij}.
\]
The useful anchored inequality is
\[
 \sum_i|t_i|-(\max_i t_i-\min_i t_i)
 \le(m-2)k+\sum_{i<j}\beta_{ij}+2|t_0|.                     \tag{4}
\]
If both signs occur, retain one largest positive and one largest negative
magnitude. Each remaining nonzero entry is bounded by one comparison with
the extreme of its own sign, using at most \(m-2\) comparisons. If all
entries have one sign and none is zero, the loss is twice the minimum plus
the sum of the middle \(m-2\) magnitudes. Bound the minimum by
\(|t_0|\), and compare each middle entry with the maximum. If a zero
occurs, there are at most \(m-2\) nonmaximum nonzero entries, which are
handled directly. This proves (4).

Exponentiating over odd split primes gives
\[
 N'\ge\frac{\prod_i F_i}{2F_0^2 K^{m-2} B_*}.               \tag{5}
\]
The factor two covers the ramified prime. Indeed, a primitive Gaussian row
has norm valuation at two either zero or one. Let \(r\le m\) rows have
odd-odd coordinates. If \(K\) is odd, each pair of these rows forces an
even \(b_{ij}\) by (2), leaving exponent at most
\(r-\binom r2\le1\). If \(K\) is even, its displayed square form
shows \(4\mid K\); the denominator supplies at least
\(1+2(m-2)\ge m\) powers of two. Inert primes divide no primitive row
norm. Thus no hidden conductor factor remains in (5).

## The resulting radius and height bounds

Cancel one factor \(F_0\) in (5), then use (3) on the other rows and
\(F_0\le T\). Combining with (1) proves
\[
 \boxed{\displaystyle
 N'\ge\frac12\left(\frac2C\right)^{2(m-1)+E}
 \frac{N^{(3m-4)/8}}{T^m K^{m-2}}.}                         \tag{6}
\]
The source exponent is \((m-1)/2-m/8=(3m-4)/8\). This formula
supersedes the preliminary all-pair loss \(K^E\) and the determinant
factor \(|\delta|^{m-2}\).

For \(H=\max|M_{ab}|\), we have \(T\le4H^2\) and
\(K\le16H^4\). Therefore
\[
 N'\le N,\quad m\ge5
 \quad\Longrightarrow\quad
 H\gg_{m,C}N^{\,3(m-4)/(16(3m-4))}.                        \tag{7}
\]
The exponent tends to \(1/16\) as fixed \(m\) increases. Formula (7)
holds for every rational Möbius matrix with \(K>0\), not just determinant
one. It rules out bounded-height descent at unbounded source radius.
It does not rule out maps whose heights grow, and it does not improve the
current general point-count growth rate.

## A determinant gain from joint accounting

There is a further improvement. Let D_delta be the entire positive part
of det(M) supported on primes not dividing K, including two if K is odd.
Jointly tracking ordinary row contents and the overlap of normalized
Gaussian images proves that the right side of (6) can be multiplied by

`max(1, D_delta^(2m-3)/T^2)`.

At each odd determinant prime not dividing K, the combined row-content
and phase-overlap loss is at most three determinant valuations plus the
old pair-residue valuations. A separate ramified calculation gives the
same conclusion at two with the existing factor two. See
[the joint proof](mobius_joint_content_conductor.md).
For the family `M=diag(n+1,n)`, the resulting estimate is

`N' >>_(m,C) N^[(3m-4)/8]/n^6`.

The exponent six in this matrix-parameter loss is independent of m.
The result constrains this family of maps, without bounding the number
of points of an arbitrary source cluster.

## Removing artificial output-rotation height

The common-content charge has subsequently been sharpened. Define
\(J=\gcd(K,|\delta|N)\). The exact transfer
\[
 |\gcd(B_i,B_j)|^2\mid\gcd(K,|\delta|e_{ij})b_{ij}
\]
allows every occurrence of \(K^{m-2}\) in (5)--(6), including the joint
determinant refinement, to be replaced by \(J^{m-2}\). Thus the strongest
version of (6) currently proved is
\[
 N'\ge\frac12(2/C)^{2(m-1)+E}
 \frac{N^{(3m-4)/8}}{T^mJ^{m-2}}
 \max\{1,D_\delta^{2m-3}/T^2\}.
\]
The complete parity argument and a separate large-trace growth estimate
for unimodular maps are in
[the source-capped content proof](mobius_source_capped_content.md).
The output minimum below also optimizes this stronger cap; see
[the exact padding and anchor formulas](mobius_capped_anchor_clipping.md).

Rational output rotations preserve all relative phases and N'. They can
still inflate the height of a primitive integral matrix, so (7) is best
used after removing this freedom. Write the complex action as
`Mz=alpha z+beta conjugate(z)`, reduce `beta/alpha=v/u` in coprime
Gaussian integers, and put `U=Norm(u)`, `V=Norm(v)`,
`epsilon=4/Norm(gcd_G(2,u-v))`. The exact simultaneous minimum in the
output-conformal class is

`T_*=epsilon(U+V)/2`, `K_*=epsilon^2 UV`.

These minima are attained by one primitive integral matrix, so they can
both be substituted into (6), with the source anchor unchanged. In
particular, a radius-nonincreasing image with m>=5 requires

`U+V >>_(m,C) N^[3(m-4)/(8(3m-4))]`.

The exponent tends to 1/8 for this squared-height invariant. The exact
normalization, its proof, and the caveat about additionally rotating the
input are in [the conformal quotient audit](mobius_conformal_quotient_normalization.md).
The invariant U+V remains unbounded over rational maps.

## Conformal maps and validation

When \(K=0\), the matrix acts as multiplication by a nonzero rational
Gaussian number, possibly followed by conjugation. All phases change by a
common rational unit-circle factor, possibly with conjugation. Every
relative phase is therefore preserved or conjugated, so its least primitive
squared radius is unchanged: \(N'=N\). Removing or fixing the transformed
anchor would invalidate this conclusion and the preceding estimates.

The persistent [exact checker](check_mobius_conductor_transfer.py) verifies
row content, pair transfer, the compressed conductor inequality, and direct
all-pair least-radius reconstruction. Its configurations include an actual
four-point Pell endpoint cluster. The check also tests conformal invariance
with all rows retained and exhaustive finite signed-valuation cases. These
are independent arithmetic checks of the displayed proofs; they are not a
proof of the uniform endpoint conjecture.
