# External anchor height is a conductor in every target pair

This note combines the exact actual-anchor formula with a condition on the
target configuration. It controls the common external part of anchor height,
including the explicit obstruction family from the previous note. It does
not control the remaining rational norm-ratio height and does not prove a
uniform endpoint bound.

Use the setup and notation of
[actual source reanchoring](mobius_actual_anchor_normalization.md). Thus
\(H_0=1\), the source half-angle columns are integer primitive of odd norm,
their least squared radius is \(N\), and the reduced complex coefficient
ratio is \(w=v/u\). For an odd split prime \(p=\pi\bar\pi\), put
\[
 d_p=v_\pi(w),\quad e_p=v_{\bar\pi}(w),\quad
 t_{i,p}=v_\pi(H_i)-v_{\bar\pi}(H_i),
\]
\[
 I_p=[\min(d_p,-e_p),\max(d_p,-e_p)],\qquad
 J_p=[\min_i t_{i,p},\max_i t_{i,p}],
\]
and define
\[
 h_p=\operatorname{dist}(I_p,J_p),\qquad
 G=\prod_{p\equiv1\pmod4}p^{h_p}.
 \tag{1}
\]
The distance is the gap between the two closed intervals, and is zero when
they meet. All quantities are integral, and only finitely many \(h_p\)
are nonzero.

## 1. Pair residues are invariant under a common phase rotation

For any two integer-primitive Gaussian columns \(B_i,B_j\), let
\[
 e_{ij}=N\gcd(B_i,B_j),\qquad
 R_{ij}=\bar B_iB_j/e_{ij}.
\]
The Gaussian integer \(R_{ij}\) is integer primitive. To check this at an
odd split prime, primitive rows each have at most one nonzero orientation
valuation; division by the pair gcd removes the common orientation, leaving
at most one orientation in \(R_{ij}\). Inert primes occur in neither
primitive row. At \(1+i\), the two row exponents are each zero or one,
and the exponent of \(R_{ij}\) is their absolute difference, hence at most
one. Thus no rational prime divides both coordinates of \(R_{ij}\).

Moreover, \(R_{ij}/\bar R_{ij}=q_j/q_i\), where
\(q_i=B_i/\bar B_i\). A primitive integer half-angle representative of a
fixed phase is unique up to sign. It follows that
\[
 \boxed{
 b_{ij}=|\det(B_i,B_j)|/e_{ij}=|\operatorname{Im}R_{ij}|}
 \tag{2}
\]
depends only on the relative phase \(q_j/q_i\). In particular it is exactly
invariant under any common rational rotation of all phases. This statement
does not require the row norms to be odd.

## 2. The common external conductor theorem

Let \(q_i'\) be the entire target configuration under \(M\), let \(N'\)
be its least squared realization radius, and let \(b_{ij}'\) be its primitive
pair residues as in (2). Then
\[
 \boxed{G\mid b_{ij}'\quad\text{for every }i\ne j,
 \qquad\gcd(G,N')=1.}
 \tag{3}
\]
The local statement requires no unit-class or odd-norm hypothesis on the
target. The source retains the odd-norm hypotheses of the anchor theorem.

Fix a prime with \(h_p>0\). Suppose first that \(I_p\) lies strictly to
the right of \(J_p\). Choose as actual anchor a source row attaining the
largest valuation in \(J_p\). After that reanchoring, every source valuation
is at most zero, while the two endpoints of the coefficient interval are
strictly positive. For the new reduced coefficients, write these endpoints
as \(d>0\) and \(-e>0\); their minimum is \(h_p\). In isotropic
coordinates their valuations are
\[
 (v_\pi(u),v_\pi(v))=(0,d),\qquad
 (v_{\bar\pi}(u),v_{\bar\pi}(v))=(-e,0).
\]
The minimal integral output-conformal representative has complex coefficients
\((\gamma/2)(u,v)\), where \(\gamma\) has Gaussian norm one, two, or
four. At the odd prime \(p\), these prefactors are units. Its isotropic
matrix has coefficient valuations
\[
 \begin{pmatrix}0&d\\0&-e\end{pmatrix}.
\]
Every reanchored source column has first isotropic coordinate a unit and
second coordinate of nonnegative valuation. Thus both isotropic output
coordinates are units: in each row of the displayed matrix the first term
is a unit and the second is divisible by \(p\). All ordinary output row
contents are therefore \(p\)-units, and every target row norm is also a
\(p\)-unit.

On the other hand, the determinant is a unit multiple of \(Nu-Nv\).
Those two norms have valuations \(-e\) and \(d\), so
\[
 v_p(\delta)\ge\min(d,-e)=h_p.
\]
The determinant identity, with the now-vanishing ordinary row contents and
target Gaussian gcd valuations, gives
\[
 v_p(b_{ij}')
 =v_p(\delta)+v_p\det(H_i^{\mathrm{new}},H_j^{\mathrm{new}})
 \ge h_p.
\]
If \(I_p\) lies to the left of \(J_p\), use the source row attaining the
minimum instead, or interchange the two Gaussian orientations. The same
argument applies.

Actual source reanchoring preserves the output phases, and output conformal
normalization preserves every \(b_{ij}'\) by (2). The prime-dependent
normalizations in this proof therefore all constrain the same target pair
residues. Multiplying the local divisibilities proves the first assertion
of (3). At each prime with \(h_p>0\), the proof also made every target
signed phase valuation zero. Their width is therefore zero, so \(p\nmid N'\).
This width is invariant under common output rotation, proving the asserted
coprimality.

## 3. Consequence when the target also satisfies the endpoint bound

Suppose the entire \(m\)-point target has least squared radius \(N'\)
and endpoint constant \(C'\). No common Gaussian unit class is required.
Let \(E=\binom m2\), and let \(r\) be the number of primitive target
columns with two odd coordinates. Their primitive pair norms satisfy
\[
 \prod_{i<j}f'_{ij}\le 2^{r(m-r)}(N')^{m^2/4}.
\]
At an odd split prime this is the signed-width Plotkin inequality. At two,
the pair norm has valuation one exactly when one of the two columns has
two odd coordinates. Multiplying the individual chord inequalities gives
\[
 \prod_{i<j}b_{ij}'\le(C'/2)^E2^{r(m-r)/2}(N')^{m/8}.
\]
Since \(G^E\) divides this product, (3) implies
\[
 \boxed{G\le A'(N')^{1/(4(m-1))},\qquad
 A'=\frac{C'}2\,2^{r(m-r)/(2E)}.}
 \tag{4}
\]
Here \(A'\le(C'/2)2^{m/(4(m-1))}\). If all target columns have odd
norm, the power of two disappears. The full configuration is retained,
including its transformed anchor; no extraction or radius change is used.

## 4. The internal anchor factor is charged to source radius

Write the exact anchor factors as before:
\[
 k_j=\prod_p p^{\operatorname{dist}(t_{j,p},I_p)}.
\]
Let \(W_p\) be the length of \(J_p\), so \(N=\prod_p p^{W_p}\).
For every source coordinate in \(J_p\),
\[
 h_p\le\operatorname{dist}(t_{j,p},I_p)\le h_p+W_p.
\]
Consequently the following divisibilities are exact:
\[
 \boxed{G\mid k_j\mid GN\quad\text{for every anchor }j.}
 \tag{5}
\]
There is also an unconditional improvement after taking a product over all
anchors. At the two attained extrema of the source interval, the sum of the
residual distances \(\operatorname{dist}(t,I_p)-h_p\) is at most \(W_p\).
This follows directly whether the intervals overlap or are disjoint. Each
of the remaining \(m-2\) residual distances is at most \(W_p\). Therefore
their total is at most \((m-1)W_p\), giving
\[
 \prod_j(k_j/G)\mid N^{m-1},\qquad
 \min_j k_j\le G N^{(m-1)/m}.
 \tag{6}
\]
Combining (4) with (5)--(6) yields
\[
 \boxed{\begin{aligned}
 k_j&\le A'N(N')^{1/(4(m-1))}\quad\text{for every }j,\\
 \min_j k_j&\le A'N^{(m-1)/m}(N')^{1/(4(m-1))}.
 \end{aligned}}
 \tag{7}
\]
The proof uses the attained source extrema and does not assume that an
interval intersection contains one of the discrete source valuations.

## 5. Effect on the previous obstruction family and remaining scope

For the family \(w=\pi^L/(2\bar\pi^L)\) with \(p\nmid N\), the map
interval at \(p\) is the single point \(L\), while the source interval is
the single point zero. Hence \(h_p=L\), and the theorem gives
\[
 p^L\mid b_{ij}'\quad\text{for every pair},\qquad
 p^L\le A'(N')^{1/(4(m-1))}
\]
whenever the target endpoint hypotheses hold. Thus the arbitrarily large
external factor in that family cannot coexist with an arbitrarily small
target endpoint radius.

The normalized trace remains
\(T_j=\epsilon(a+b)k_j/2\), and the discriminant remains
\(K_j=\epsilon^2ab\,k_j^2\). Equations (4)--(7) control the external and
internal factors in \(k_j\), but do not bound the invariant coprime integers
\(a,b\). Nor does their present size suffice to close the endpoint radius
inequality. New information is still needed to obtain a uniform bound.
