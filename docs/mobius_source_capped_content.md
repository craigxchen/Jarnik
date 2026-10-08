# Source-capped pair content and large-trace growth

This note sharpens the discriminant charge in the Möbius radius bound. The
charge can be capped by the source Gaussian pair contents and the determinant.
For unimodular maps, a separate expansion argument also turns sufficiently
large trace into radius growth. These give additional necessary conditions
for radius descent, not a uniform endpoint obstruction.

Let \(M\) be an integral real \(2\times2\) matrix, with
\[
 \delta=\det M\ne0,\qquad T=\operatorname{tr}(M^TM),\qquad
 K=T^2-4\delta^2>0.
\]
Let \(H_i\) be pairwise nonparallel integer-primitive Gaussian columns of
odd norms \(f_i\), and set
\[
 r_i=\gcd_{\mathbb Z}(MH_i),\qquad B_i=MH_i/r_i,\qquad F_i=NB_i,
\]
\[
 e_{ij}=N\gcd(H_i,H_j),\qquad
 b_{ij}=|\det(H_i,H_j)|/e_{ij},\qquad
 e_{ij}'=N\gcd(B_i,B_j).
\]

## 1. The exact pairwise cap

The determinant-free transfer theorem gives \(e_{ij}'\mid Kb_{ij}\).
Independently, a Gaussian pair gcd norm divides the integer determinant, so
\[
 e_{ij}'\mid |\det(B_i,B_j)|
 =\frac{|\delta|e_{ij}b_{ij}}{r_ir_j}.
\]
Taking the gcd of the two integral upper multiples yields
\[
 \boxed{e_{ij}'\mid\gcd(K,|\delta|e_{ij})\,b_{ij}.}
 \tag{1}
\]
The ordinary contents give a still sharper valuation statement. For every
rational prime, write \(k=v_pK\), \(\ell=v_p\delta\),
\(s=v_pe_{ij}\), \(\beta=v_pb_{ij}\), and \(c_r=v_pr_r\). Then
\[
 \boxed{v_p(e_{ij}')\le\beta+\min(k,\ell+s-c_i-c_j).}
 \tag{2}
\]
The second entry of this minimum can be negative; (2) is a valuation
inequality and does not assert that this entry defines an integer divisor.
The simpler integral divisor (1) is sufficient below.

Suppose the source tuple includes \(H_0=1\), and let \(N\) be its least
squared realization radius. Then \(e_{ij}\mid N\). Indeed, at every split
prime, the source signed valuations include zero, and their width is the
exponent of \(N\). A same-orientation pair minimum is at most that width.
Inert and ramified primes divide no source norm under the stated hypotheses.
Consequently the common cap
\[
 \boxed{J=\gcd(K,|\delta|N)}
 \tag{3}
\]
satisfies \(e_{ij}'\mid Jb_{ij}\) for every pair. For a unimodular map,
\(J=\gcd(K,N)\le N\). New discriminant primes outside the source radius
then incur no pair-content charge.

## 2. Radius bounds with the capped content

Let \(N'\) be the least squared realization radius of the entire target
phase configuration, and assume \(m\ge3\). The unanchored signed-width
inequality, together with (1), gives the pair-specific bound
\[
 \boxed{
 N'\ge\frac{\prod_iF_i}
 {2\prod_{i<j}\gcd(K,|\delta|e_{ij})\,b_{ij}}.}
 \tag{4}
\]
With the retained source anchor, the compressed signed-width inequality gives
\[
 \boxed{
 N'\ge\frac{\prod_iF_i}
 {2F_0^2J^{m-2}\prod_{i<j}b_{ij}}.}
 \tag{5}
\]
At odd split primes this is the usual anchored inequality with common charge
\(v_pJ\), now justified by (1)--(3). The parity audit remains valid even
when \(v_2J=1\), unlike the possible valuations of \(K\) itself:

* If \(v_2J=0\), every pair of odd-odd primitive outputs forces its
  \(b_{ij}\) to be even. For \(r\) such rows, the uncancelled exponent
  is at most \(r-\binom r2\le1\).
* If \(v_2J\ge2\), the factor \(2J^{m-2}\) supplies at least
  \(1+2(m-2)\ge m\) powers of two, covering every row norm parity.
* If \(v_2J=1\), that factor supplies \(m-1\). This suffices unless
  every output is odd-odd. In that remaining case, the \(m\ge3\) odd-norm
  source rows occupy only two nonzero parity directions. Two have the same
  direction, so their source determinant is even. Their old Gaussian gcd
  norm is odd, hence the corresponding \(b_{ij}\) supplies the last factor
  two.

The usual primitive row estimate \(F_i\ge f_i/T\), and \(F_0\le T\),
then imply the endpoint consequence below. Let \(E=\binom m2\), and suppose
the primitive fixed-unit source cluster has constant \(C\le2\), so
\[
 f_i\ge(4/C^2)\sqrt N\quad(i\ne0),\qquad
 \prod_{i<j}b_{ij}\le(C/2)^E N^{m/8}.
\]
Then
\[
 \boxed{
 N'\ge\frac12(2/C)^{2(m-1)+E}
 \frac{N^{(3m-4)/8}}{T^mJ^{m-2}}.}
 \tag{6}
\]

The joint determinant refinement remains available with the same cap. If
\(D\) is the full part of \(|\delta|\) supported on primes not dividing
\(K\), then the lower bound (6) can be multiplied by
\[
 \max\{1,D^{2m-3}/T^2\}.
 \tag{7}
\]
The proof at good determinant primes is unchanged, and all other compressed
charges use \(J\). For primes dividing \(\delta\), the conditions
\(p\nmid K\) and \(p\nmid J\) are equivalent, so this does not redefine
which determinant primes are usable.

## 3. All but one source row must expand when the trace is large

Now specialize to \(|\delta|=1\). Every output row is already integer
primitive. For a pair of source directions, its primitive relative
half-angle numerator has norm \(f_{ij}\mid N\), and imaginary coordinate
of absolute value \(b_{ij}\ge1\). Therefore
\[
 |\sin\angle(H_i,H_j)|=b_{ij}/\sqrt{f_{ij}}\ge N^{-1/2}.
 \tag{8}
\]
Let \(z\) be a unit eigenvector for the larger eigenvalue of \(M^TM\).
That eigenvalue is at least \(T/2\). There cannot be two source unit
directions whose absolute projections onto \(z\) are both less than
\(1/(2\sqrt N)\): expanding their determinant in the basis \((z,z^\perp)\)
would make its absolute value less than \(1/\sqrt N\), contrary to (8).

Thus there is at most one exceptional row, and every other row satisfies
\[
 \boxed{F_i\ge\frac{T}{8N}f_i.}
 \tag{9}
\]
The exceptional row, if present, still satisfies \(F_i\ge1\), because
it is a nonzero integer vector. This avoids paying an inverse-eigenvalue
factor for that row.

At least \(m-2\) nonanchor rows are nonexceptional. Since
\((4/C^2)\sqrt N\ge1\), the endpoint norm estimates give
\[
 \prod_iF_i\ge
 \left(\frac{T}{8N}\right)^{m-1}
 (4/C^2)^{m-2}N^{(m-2)/2}.
\]
Substitute this into (5), use \(F_0^2\le T^2\), and use the source
pair-product estimate. The resulting genuinely growing-trace lower bound is
\[
 \boxed{
 N'\ge c_{m,C}\frac{T^{m-3}}{J^{m-2}N^{5m/8}},\qquad
 c_{m,C}=\frac{(2/C)^{E+2(m-2)}}{2\,8^{m-1}}.}
 \tag{10}
\]
All constants are explicit and independent of the source configuration.

For \(m\ge4\), a radius-nonincreasing image must therefore satisfy
\[
 \boxed{
 T\le c_{m,C}^{-1/(m-3)}
 \left(N^{1+5m/8}J^{m-2}\right)^{1/(m-3)}.}
 \tag{11}
\]
Since \(J\le N\), this implies
\[
 T\ll_{m,C}N^{(13m-8)/(8(m-3))},\qquad
 \max|M_{ab}|\ll_{m,C}N^{(13m-8)/(16(m-3))}.
 \tag{12}
\]
The latter exponent tends to \(13/16\). This upper restriction on the
height of an actual unimodular radius-descent map complements the earlier
lower restriction; the two leave a substantial range and do not contradict
one another. Equation (11) is stronger when the capped content \(J\) is
smaller than \(N\). No corresponding large-trace conclusion is asserted
here for arbitrary determinants, where ordinary row contents can remove
the expanding factors.
