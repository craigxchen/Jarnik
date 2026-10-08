# Finite samples force the critical condition scale and a near-shear form

This note uses only the actual source and target points when estimating the
geometry of a rational Möbius map. It does not assume that the entire arc
between sampled source points has a small image. Two finite-sample lemmas
give new restrictions:

* the condition number is at most a product of explicit powers of the two
  radii, with total exponent \(1/2+O(1/m)\) when \(N'\le N\);
* angular pair products select an actual anchor with controlled normalized
  stretch. Combining this selection with the established conductor budgets
  gives an integral triangular representative with small diagonal entries.

Together with the previous lower condition bound, these results isolate a
critical near-shear regime. They do not exclude that regime or prove a
uniform endpoint count.

## 1. A finite-subset diameter bound

Suppose \(q\ge2\) distinct rational circle phases are realized on an
inherited integer circle of squared radius \(N\). Do not divide out a new
subset gcd. Let \(f_{ij}\) be their primitive relative half-angle norms,
and \(b_{ij}\) the corresponding nonzero integer pair residues. Put
\[
 h_q=\frac{\lfloor q^2/4\rfloor}{q(q-1)}
 =\begin{cases}q/[4(q-1)],&q\text{ even},\\
 (q+1)/(4q),&q\text{ odd}.
 \end{cases}
\]
The odd-prime layer bound and the binary parity cut give
\[
 \prod_{i<j}f_{ij}\le(2N)^{\lfloor q^2/4\rfloor}.
\]
Since
\(b_{ij}=\sqrt{f_{ij}}\,|\sin\theta_{ij}|\), where \(\theta_{ij}\) is
the angle between the half-angle directions, integrality yields
\[
 \boxed{\max_{i<j}|\sin\theta_{ij}|\ge(2N)^{-h_q}.}
 \tag{1}
\]
Indeed, multiply the bound \(b_{ij}\le D\sqrt{f_{ij}}\), with \(D\)
the maximum sine, and use \(\prod b_{ij}\ge1\). No common-unit extraction
or endpoint hypothesis is needed for this lemma.

## 2. A shared-pivot singular-vector argument

Let \(M\) map a configuration of \(m\) phases on squared radius \(N\)
to one on squared radius \(N'\). Let \(\kappa\ge1\) be its Euclidean
condition number. Choose unit singular coordinates in the source, and order
the source directions by the absolute value of their projection onto the
expanding singular vector.

Choose \(2\le n\le m-1\), put \(r=m-n+1\), and let \(a\) be the
\(n\)-th smallest such absolute projection. At most one projective direction
has zero projection, so \(a>0\). The \(n\) nearest directions all have
expanding projection at most \(a\). Expanding their pair determinants in
singular coordinates gives
\[
 |\sin\theta_{ij}|\le2a
\]
on that source subset.

The \(r\) farthest directions, including the same pivot, have expanding
projection at least \(a\). After applying \(M\), each corresponding unit
target direction has contracting projection at most \(1/(\kappa a)\).
Their pair sines are therefore at most \(2/(\kappa a)\). Apply (1) to
both subsets to obtain
\[
 2a\ge(2N)^{-h_n},\qquad
 \frac2{\kappa a}\ge(2N')^{-h_r}.
\]
Eliminating \(a\) proves the finite-sample condition bound
\[
 \boxed{\kappa\le4(2N)^{h_n}(2N')^{h_r},
 \qquad n+r=m+1.}
 \tag{2}
\]
This argument is in projective unit directions and crosses affine poles
without any change. Its two subsets need not be the same subset, and no
assertion is made about intermediate points on either arc.

For example, any nine sampled points allow \(n=r=5\), giving
\[
 \kappa\le4(2N)^{3/10}(2N')^{3/10}.
\]
The subset sizes can grow with \(m\). If \(N'\le N\), define
\[
 u_m=\min_{\substack{n,r\ge2\\n+r=m+1}}(h_n+h_r).
\]
Then
\[
 \boxed{\kappa\le4(2N)^{u_m},\qquad
 u_m=\tfrac12+\tfrac1m+O(m^{-2}).}
 \tag{3}
\]
There is an exact formula. Let \(L=2\lfloor(m+1)/2\rfloor\). For \(m\ge5\),
\[
 u_m=\begin{cases}
 1/2+1/L,&L\equiv2\pmod4,\\
 1/2+L/(L^2-4),&L\equiv0\pmod4.
 \end{cases}
\]
To verify the optimization, write
\(h_q=1/4+1/(4\ell_q)\), where \(\ell_q\) is the largest odd integer
at most \(q\). The two effective odd sizes should have maximal total \(L\)
and be as equal as possible. This gives the formula directly.

Combined with the arbitrary-unit lower condition theorem, (3) confines
nonconformal endpoint-to-endpoint maps with \(N'\le N\) and \(m\ge8\) to
\[
 \boxed{\frac1{512}N^{z_m}\le\kappa\le4(2N)^{u_m},
 \quad z_m=\tfrac12-\tfrac3m+O(m^{-2}).}
 \tag{4}
\]
Thus the condition scale is \(N^{1/2+O(1/m)}\), rather than an arbitrary
large condition number.

## 3. Angular products control actual anchor stretch

Now assume both full configurations obey endpoint constants at most two.
Write \(W=\log N\), \(W'=\log N'\). Let \(h_i\) be the source unit
half-angle vector and define the normalized stretch
\[
 s_i=\frac{\|Mh_i\|}{\sqrt{|\det M|}},\qquad
 x_i=\log s_i-\frac{W'-W}{8}.
\]
The exact pair determinant identity is
\[
 |\sin\theta_{ij}'|
 =\frac{|\sin\theta_{ij}|}{s_is_j}.
\]
Define nonnegative angular deficits
\[
 d_{ij}=-W/4-\log|\sin\theta_{ij}|,\qquad
 d'_{ij}=-W'/4-\log|\sin\theta'_{ij}|.
\]
Nonnegativity follows from the endpoint chord bound. The determinant identity
becomes the finite-sample product identity
\[
 \boxed{x_i+x_j=d'_{ij}-d_{ij}.}
 \tag{5}
\]
If \(\sigma\) source half-angle rows and \(r\) target rows have even norms,
integral pair residues and the full pair-norm product bounds imply
\[
 \sum_{i<j}d_{ij}\le mW/8+\sigma(m-\sigma)\log2/2,
\]
and the same inequality with primes and \(r\). Therefore
\[
 \sum_{i<j}|x_i+x_j|
 \le m(W+W')/8+m^2\log2/4.
\]

For arbitrary real numbers \(x_i\), let
\(F_m=\lfloor(m-1)^2/4\rfloor\). The exact useful inequality is
\[
 \sum_{i<j}|x_i+x_j|\ge\frac{2F_m}{m}\sum_i|x_i|.
\]
For a proof, integrate signed superlevel sets of \(|x_i|\). At one level
with \(p\) positive and \(n\) negative entries, and \(v=p+n\), the pair
sum is \(v(m-1)-2pn\). Dividing by \(v\), this is at least
\(m-1-2\lfloor m^2/4\rfloor/m=2F_m/m\); the minimum occurs when all
\(m\) entries are present and their signs are as balanced as possible.
Integration proves the inequality.

Consequently
\[
 \boxed{\frac1m\sum_i|x_i|
 \le\frac{m(W+W')}{16F_m}+\frac{m^2}{8F_m}\log2.}
 \tag{6}
\]
In particular an actual source anchor has normalized squared stretch
\(s_i^2\) within the explicit factors
\[
 2^{\pm m^2/(4F_m)}(NN')^{\pm m/(8F_m)}
\]
of \((N'/N)^{1/4}\). This is stronger than bounding only the condition
number, which concerns extreme singular directions that need not be sampled.

## 4. One anchor simultaneously controls content and stretch

Assume henceforth \(N'\le N\) and \(m\ge8\). Use \(A_m,c_m,\rho_m\)
and
\[
 q_m=\frac{m(1+2A_m)}{8c_m},\quad
 Q_m=\lfloor m^2/4\rfloor
\]
from the arbitrary-unit erased-conductor audit. Its explicit bound is
\[
 E_0\le q_m(W+m\log2).
\]
The odd-prime overlap identity gives
\[
 (m-1)\sum_j\log k_j
 \le2\log B'+\log P_{\rm odd}-\log P_{\rm clip}
 \le mW/4+(m^2/4)\log2+Q_mE_0.
\]
The last step charges each erased source layer by at most \(Q_m\).
Define
\[
 \ell_m=\frac1{4(m-1)}+\frac{Q_mq_m}{m(m-1)},\qquad
 b_m=\ell_m+\frac{m}{4F_m}.
\]
Then
\[
 \frac1m\sum_j\log k_j\le\ell_m(W+m\log2),
\]
while (6) bounds the average of \(2|x_j|\) by
\(m(W+m\log2)/(4F_m)\). Hence there is one actual anchor \(j\) with
\[
 \boxed{\log k_j+2|x_j|\le b_m(W+m\log2).}
 \tag{7}
\]
In particular \(b_m=9/(4m)+O(m^{-2})\). This simultaneous choice does
not incorrectly combine the best content anchor with a different best
stretch anchor.

Take the minimal integral representative at this anchor, with determinant
\(\delta_j\), first-column norm \(F_0\), and core \(c_j\). Then
\[
 |\delta_j|=c_jk_j,\qquad F_0=s_j^2|\delta_j|.
\]
An even-norm anchor can change \(\epsilon_j\), but \(1\le\epsilon_j\le4\)
implies \(c_j\le4c_0\) relative to the fixed anchor used in the full-core
theorem. Put
\[
 g_m=\frac{m}{4F_m},\qquad
 D_m=4\,2^{m(m-1)/F_m+mb_m}.
\]
The full-core bound and (7), using \((N'/N)^{1/4}\le1\), give
\[
 \boxed{1\le F_0\le D_mN^{g_m+b_m},\qquad
 1\le|\delta_j|\le D_mN^{g_m+b_m}.}
 \tag{8}
\]
The common exponent is \(g_m+b_m=13/(4m)+O(m^{-2})\).

## 5. An integral triangular representative with small diagonals

Write the two integer columns of the selected map as \((x,y)\) and
\((p,q)\). Output multiplication by \(x-iy\) produces the exact integral
matrix
\[
 \begin{pmatrix}x&y\\-y&x\end{pmatrix}M_j
 =\begin{pmatrix}F_0&xp+yq\\0&\delta_j\end{pmatrix}.
\]
Divide its ordinary common content if desired. Thus, after actual source
and target anchoring, there is a primitive integral triangular representative
\[
 \boxed{\widetilde M=\begin{pmatrix}a&b\\0&d\end{pmatrix},\qquad
 1\le a\le D_mN^{g_m+b_m},\quad
 1\le|d|\le D_mN^{g_m+b_m}.}
 \tag{9}
\]
The inequality also includes \(1\le|d|\); \(a>0\) fixes the projective
sign. Let \(R=T/|\delta|=\kappa+\kappa^{-1}\). For this representative,
\[
 R=\frac a{|d|}+\frac{|d|}a+\frac{b^2}{a|d|}.
 \tag{10}
\]
The finite-sample bound implies \(R\le5(2N)^{u_m}\), so
\[
 \boxed{|b|\le\sqrt5\,2^{u_m/2}D_m
 N^{u_m/2+g_m+b_m}
 =N^{1/4+O(1/m)}\text{ up to the displayed factor}.}
 \tag{11}
\]
In affine half-angle coordinates at the actual anchors, the map is exactly
\[
 t' =\frac{dt}{a+bt}.
\]
The denominator is nonzero at the sampled points. It need not remain
nonzero, or have a small image, throughout an interval containing them.

Equations (7) and the radius comparison
\(N'\ge N^{\rho_m}/16\) also bound the actual anchor ratio:
\[
 2^{-1-mb_m}N^{-b_m-(1-\rho_m)/4}
 \le a/|d|\le2^{mb_m}N^{b_m}.
 \tag{12}
\]
If \(z_m>b_m+(1-\rho_m)/4\), equations (10), (12), and the lower bound
\(R\ge N^{z_m}/256\) force the shear term to dominate for sufficiently
large \(N\). In particular then \(|b|\gg_m N^{z_m/2}\), since
\(a|d|\ge1\). This condition holds for all sufficiently large \(m\);
the exact displayed inequality is the hypothesis, rather than an assumed
uniform lower bound on every finite source angle.

There is an exact arithmetic near-square consequence for the triangular
matrix, whether or not its ordinary content is divided out. Its discriminant
is
\[
 K_{\rm tri}=(a^2+b^2+d^2)^2-4a^2d^2.
\]
Before dividing ordinary content it also equals \(F_0^2K_j\). More
intrinsically, for any integral representative write the Gaussian integer
coefficients \(U=2\alpha,V=2\beta\). At each odd prime,
\[
 |v_p(NV)-v_p(NU)|\le v_p(NU)+v_p(NV)=v_p(K).
\]
The left side is the invariant coefficient-interval length, so it bounds
every clipped width from above independently of matrix content. Therefore
\[
 \boxed{N_{\rm clip}\mid K_{\rm tri},\qquad
 4a^2d^2\le4D_m^4N^{4(g_m+b_m)}.}
 \tag{13}
\]
The divisor \(N_{\rm clip}\ge2^{-mq_m}N^{1-q_m}\) contains almost all
source logarithmic conductor, while the gap from the displayed square has
size \(N^{O(1/m)}\). The interval-length proof shows that the displayed
divisibility survives ordinary primitive reduction; one need not discard
source conductor merely because the numerical discriminant loses a fourth
power under that reduction.

## 6. Why the critical case is still compatible with finite geometry

The near-shear form is a reduction, not a contradiction. For example,
\[
 M_B=\begin{pmatrix}1&B\\0&1\end{pmatrix},\qquad
 t'=t/(1+Bt)
\]
has condition number of order \(B^2\) and determinant one. Choose any fixed
finite set of distinct rational \(a_i\in[1,2]\), together with zero, and
put \(t_i=a_i/B\). Both these source samples and their images have angular
diameter of order \(B^{-1}\). Relative to the formal scale \(N=B^4\),
this is exactly the endpoint angular scale and \(\kappa\) is of order
\(\sqrt N\). No power saving in either sampled window follows from the
singular-vector geometry alone.

This example does not assert that the least integer-circle radii of these
rational samples equal \(B^4\). For many points they generally do not.
The missing arithmetic compatibility between that critical geometric scale
and the actual least radii is precisely the unresolved issue.

Conversely, if \(\kappa\ge N^{1/2+\varepsilon}\) with \(N'\le N\),
the elementary threshold \(a=\kappa^{-1/2}\) partitions the points so that
at least half the points lie in a source or target chord window of size at
most four times its own endpoint length scale with power saving
\(\varepsilon/2\). Lemma (1) excludes this when the subset is sufficiently
large compared with \(1/\varepsilon\). This is the mechanism behind the
upper exponent in (3); it does not yield a saving at the critical scale.

The earlier fixed-six work studies a fixed source orbit, whereas (2), (6),
and (9) are uniform over moving source configurations. The integer triangular
form has diagonals of size \(N^{O(1/m)}\) and shear of size
\(N^{1/4+O(1/m)}\). Excluding large endpoint configurations for precisely
these maps would require a further arithmetic argument. No such exclusion
or uniform point-count theorem is claimed here.

The independent [finite-geometry checker](check_mobius_finite_geometry.py)
passes 11,404 integer profiles for the signed-sum inequality, 345 inherited
subset diameter checks, and 2,270 matrix/source cases including mixed
half-angle parities and shears up to \(10^6\). It verifies 28,148 exact
squared stretch identities and rational upper certificates for the
condition number for every permitted pivot count in those cases.
