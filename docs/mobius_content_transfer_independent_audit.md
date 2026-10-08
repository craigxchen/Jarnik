# Integral Möbius content transfer: independent audit

This note proves a bounded-content statement for a common matrix in
\(SL_2(\mathbb Z)\). It does not establish the uniform short-circle-arc
bound. The later sections connect the pair bounds to the common circle
conductor with an explicit matrix-height loss. The determinant-free
extension to all rational matrices is proved separately in
[`mobius_determinant_free_transfer_audit.md`](mobius_determinant_free_transfer_audit.md).

Let
\[
 M=\begin{pmatrix}a&b\\c&d\end{pmatrix},\quad
 Q=M^TM=\begin{pmatrix}A&B\\B&D\end{pmatrix},\quad
 K=(A-D)^2+4B^2=(A+D)^2-4.
\]
Here vectors are columns. If the matrix entries have absolute value at most
\(H\), then \(0\le K\le16H^4\). Also \(K=0\) precisely when
\(Q=I\), that is, when \(M\) is an integral rotation. All divisibility
statements involving a numerical upper bound below require \(K>0\).

## One primitive row

For an integer primitive vector \(h=(x,y)^T\), put
\(f=x^2+y^2\) and \(F=h^TQh\). Then
\[
 \gcd(f,F)\mid K. \tag{1}
\]
Indeed, let \(n=\gcd(f,F)\). Each prime dividing \(n\) cannot divide
\(y\), since it would then divide \(x\). Thus \(t=xy^{-1}\) exists
modulo \(n\) and satisfies
\[
 t^2=-1,\qquad 2Bt+(D-A)=0\pmod n.
\]
Squaring the second relation gives
\((D-A)^2=-4B^2\pmod n\). This proves (1) for arbitrary prime powers,
including powers of two, without division by two.

Consequently, across any collection of primitive rows,
\[
 \operatorname{lcm}_i\gcd(f_i,F_i)\mid K. \tag{2}
\]
A product version of (2) is justified only for pairwise coprime selected
factors. Repeated appearances of the same prime across rows must not be
multiplied. For example, the matrix \(\left(\begin{smallmatrix}1&1\\0&1
\end{smallmatrix}\right)\) has \(K=5\); both primitive rows
\((2,1)\) and \((7,1)\) have old and new norm divisible by five.
Their two gcd factors have product 25, which does not divide \(K\).

## Two primitive rows

Identify a primitive vector with its Gaussian integer. Let \(h_i,h_j\)
be nonparallel primitive integer vectors, and write
\[
 g=\gcd_{\mathbb Z[i]}(h_i,h_j),\quad e=N(g),\qquad
 g'=\gcd_{\mathbb Z[i]}(Mh_i,Mh_j),\quad e'=N(g').
\]
Gaussian gcds are defined up to units, so their norms are unambiguous.
The determinant is divisible by \(e\): writing \(h_i=gu\) and
\(h_j=gv\) gives
\[
 \det(h_i,h_j)=e\det(u,v)=e b_{ij}.
\]
Because \(\det M=1\), the same determinant is also divisible by \(e'\).
Then
\[
 e'\mid K|b_{ij}|. \tag{3}
\]
To prove this at a prime \(p\), set
\(\alpha=v_p(e)\), \(\beta=v_p(e')\),
\(t=v_p(|b_{ij}|)\), and \(k=v_p(K)\).
Determinant divisibility gives \(\beta\le\alpha+t\). Since \(e\mid f_i\)
and \(e'\mid F_i\), (1) gives \(\min(\alpha,\beta)\le k\).
If \(\beta\le k\), the claimed bound is immediate; otherwise
\(\alpha\le k\) and again \(\beta\le k+t\). This proves (3).
No conjugate-coprimality hypothesis is needed.

There is a useful interpretation of (3) in terms of ordinary coordinate
content. For integer primitive Gaussian integers \(z,w\),
\[
 \gcd\bigl(\Re(z\bar w),\Im(z\bar w)\bigr)
 =N\bigl(\gcd_{\mathbb Z[i]}(z,w)\bigr). \tag{4}
\]
At a split prime, integer primitivity allows positive valuation on only
one of the two conjugate Gaussian primes in each row; the valuation of the
coordinate content is therefore the common-orientation minimum. An inert
prime divides neither row. At the ramified prime \(1+i\), each row has
valuation at most one, and coordinate content has a factor two exactly
when both rows have that valuation. These observations prove (4).

Thus if \(P'_{ij}=(Mh_i)\overline{(Mh_j)}/e'\) denotes the Gaussian
integer obtained by removing its integer coordinate content, then
\[
 N(P'_{ij})=\frac{F_iF_j}{(e')^2}
 \ge \frac{F_iF_j}{K^2 b_{ij}^2}. \tag{5}
\]
The norm in (5) belongs to this specified primitive Gaussian pair factor;
any subsequent circle parametrization and its possible factor of two must
be tracked separately. Multiplying (5) over edges is valid as an inequality,
but identifying that product with a common circle conductor requires its
own cancellation accounting.

## Finite verification

An independent exhaustive integer check used all 116 determinant-one
matrices with entries in \([-3,3]\), and all 80 primitive vectors in
\([-5,5]^2\). It checked (1), and checked (3) for 349,440 nonparallel
unordered row pairs with \(K>0\), computing Gaussian gcd norms through
(4). All checks passed. The arguments above, rather than this finite check,
cover arbitrary heights and prime powers.

## Independent check of the global extension

For \(m\ge3\) arbitrary nonzero integer-primitive Gaussian rows \(B_i\),
let \(N'\) denote the least squared radius of a Gaussian-integer realization
of the phases \(B_i/\bar B_i\), allowing one common rational Gaussian
multiplier. Then
\[
 N'\ge\frac{\prod_iN(B_i)}{2\prod_{i<j}N(\gcd(B_i,B_j))}. \tag{6}
\]
At a split odd prime, let \(t_i\) be the signed row valuation. The least
radius exponent is \(\max t_i-\min t_i\): the common multiplier requires
valuations at least \(-\min t_i\) and \(\max t_i\) at the two Gaussian
primes. Their sum is the displayed width, even if one multiplier valuation
is negative. Integer primitivity gives row norm exponent \(|t_i|\), and
the pair gcd exponent is the same-sign minimum (zero for opposite signs).
Thus the needed inequality is
\[
 \max t_i-\min t_i\ge
 \sum_i|t_i|-\sum_{i<j,\ t_it_j>0}\min(|t_i|,|t_j|).
\]
For a sign group with ordered magnitudes \(a_1\le\cdots\le a_s\), its
contribution on the right is
\(a_s+\sum_{j<s}(j-s+1)a_j\), which is at most \(a_s\).
This proves the mixed-sign case and the case with a zero. If every entry
has the same nonzero sign, \(s=m\ge3\) and the coefficient of \(a_1\)
is \(2-m\le-1\), proving the sharper bound \(a_m-a_1\).
At two, if \(r\) rows have odd coordinates, the exponent of the quotient
on the right before its factor two is \(r-\binom r2\le1\).
There are no inert-prime row factors. This proves (6).

Combining (6) with (3) and \(F_i\ge f_i/T\), where
\(T=\operatorname{tr}(M^TM)\), gives
\[
 N'\ge\frac{\prod_i f_i}{2T^m K^E\prod_{i<j}|b_{ij}|},
 \qquad E=\binom m2. \tag{7}
\]
For an anchored source with \(H_0=1\),
\(f_i\ge(4/C^2)\sqrt N\) for \(i\ne0\), and
\(\prod|b_{ij}|\le(C/2)^E N^{m/8}\), this yields
\[
 N'\ge\frac12(2/C)^{2(m-1)+E}
 \frac{N^{(3m-4)/8}}{T^mK^E}. \tag{8}
\]
Indeed the standard layer-cake bound gives
\(\prod d_{ij}\le N^{m^2/4}\); multiplying
\(|b_{ij}|\le(C/2)\sqrt{d_{ij}}N^{-1/4}\) gives exactly the exponent
\(m^2/8-E/4=m/8\). With \(H=\max|M_{ab}|\), the denominator loss obeys
\(T^mK^E\le4^m16^E H^{2m^2}\). Therefore \(N'\le N\), \(m\ge5\), and
\(0<C\le2\) imply
\[
 H\gg_m N^{3(m-4)/(16m^2)}.
\]
The exponent is largest at \(m=8\), with value \(3/256\).

**Source normalization qualification.** For these particular source
hypotheses, the usual construction \(H_i=A_{0i}\) of odd norm applies to
points in a single Gaussian-unit class. A general cluster is partitioned
into at most four such classes, as in the existing near-real divisor
reduction; hence \(m\) here counts the retained class. Residual factors
\(\pm i\) in general point ratios cannot be absorbed by multiplying an
odd-norm row by a Gaussian unit, since this changes its phase by only
\(\pm1\). Within a fixed unit class, the primitive pair factor of
\(H_i\bar H_j\) equals the source pair factor up to sign (or conjugation
if pair order is reversed), so the absolute \(b_{ij}\) in (7) is exactly
the one in the source chord and Plotkin estimates. The general cluster
still has the stated obstruction after this bounded extraction. The
estimate does not supply a height bound for unrestricted Möbius maps and
does not resolve the uniform endpoint question.

## Sharpening the matrix-content exponent to \(m-2\)

For an anchored source \(H_0=1\), the loss \(K^E\) above can be replaced
by \(K^{m-2}\), with an intermediate factor \(F_0^2\). This is a stronger
bound than the preceding all-pairs multiplication argument.

At an odd split prime let \(k=v_p(K)\),
\(\beta_{ij}=v_p(|b_{ij}|)\), and let \(t_i\) be the signed transformed
row valuations. Pair transfer implies
\[
 \min(|t_i|,|t_j|)\le k+\beta_{ij}\qquad(t_it_j>0).
\]
The following bound is valid for \(m\ge3\):
\[
 \sum_i|t_i|-(\max_i t_i-\min_i t_i)
 \le(m-2)k+\sum_{i<j}\beta_{ij}+2|t_0|. \tag{9}
\]
When both signs occur, retain an entry of maximum magnitude in each sign.
The left side is the sum of the other magnitudes. Each nonzero remaining
entry is bounded by \(k+\beta_{ij}\) on its edge to its own sign's retained
maximum. These are distinct edges and there are at most \(m-2\) of them.
If all entries have one weak sign, order their magnitudes
\(a_1\le\cdots\le a_m\). The left side equals
\(2a_1+\sum_{j=2}^{m-1}a_j\). Each nonzero middle term is bounded using its
edge to the maximum; zero middle terms need no edge. Finally
\(a_1\le|t_0|\). This proves (9), including zeros and repeated extrema.

At two, every transformed row norm has valuation zero or one. If \(K\)
is odd, pair transfer forces \(\beta_{ij}\ge1\) whenever both transformed
rows have odd coordinates. If there are \(r\) such rows, the remaining
exponent is at most \(r-\binom r2\le1\). If \(K\) is even, then
\(K=T^2-4\) implies \(v_2(K)\ge2\); thus
\(1+(m-2)v_2(K)\ge m\) covers every possible row-norm contribution.
In both cases the additional factor \(F_0^2\) can only help.
Consequently (9), multiplied over primes, yields
\[
 N'\ge\frac{\prod_iF_i}
 {2F_0^2K^{m-2}\prod_{i<j}|b_{ij}|}. \tag{10}
\]
Because \(F_0\le T\) and \(F_i\ge f_i/T\) for the other \(m-1\) rows,
\[
 \frac{\prod_iF_i}{F_0^2}
 =\frac{\prod_{i\ne0}F_i}{F_0}
 \ge T^{-m}\prod_{i\ne0}f_i.
\]
Thus in the same fixed-unit-class clustered source regime as above,
\[
 N'\ge\frac12(2/C)^{2(m-1)+E}
 \frac{N^{(3m-4)/8}}{T^mK^{m-2}}. \tag{11}
\]
The matrix-height loss is at most
\(4^m16^{m-2}H^{6m-8}\). For \(m\ge5\), \(C\le2\), and \(N'\le N\),
this gives the sharpened necessary condition
\[
 H\gg_m N^{3(m-4)/(16(3m-4))}. \tag{12}
\]
The exponent increases toward \(1/16\) as \(m\) increases. This remains a
scoped obstruction to small-height descent, not a bound for arbitrary-height
maps or a solution of the endpoint question.

An exhaustive independent check of (9) tested all
\(t_i\in\{-3,-2,\ldots,3\}\), \(k\in\{0,1,2,3\}\), and
\(m\in\{3,4,5,6\}\), taking each \(\beta_{ij}\) to be its least admissible
nonnegative value. All 548,800 checks passed. The proof above does not rely
on this finite range.
