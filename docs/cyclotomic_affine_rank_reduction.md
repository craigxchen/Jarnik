# Exact affine allocation bridge for one resonant Fibonacci class

This note connects the one-class cyclotomic fibre in
[resonant_fibonacci_cyclotomic_reduction.md](resonant_fibonacci_cyclotomic_reduction.md)
to the exact geometric theorem in
[linear_allocation_affine_rigidity.md](linear_allocation_affine_rigidity.md).
It applies to a **growing fixed finite template** of Gaussian Fibonacci
factors, with positive effective layer degree. Its conclusion is a count
for actual points in that template, not a uniform count over arbitrary
circle radii or stationary templates.

## 1. Hypotheses and the phase constraint

Fix odd multiples \(m_s\), multiplicities \(r_s\), row exponents
\(0\leq a_{is}\leq r_s\), and nonzero Gaussian prefactors \(u_i\)
of one common norm. Let \(d_0(n)=An+B>0\) eventually, with
\(A>0\) even and \(B\) odd, and put
\[
Z_i(n)=u_i\prod_s
G_{m_s d_0(n)}^{a_{is}}
\overline{G_{m_s d_0(n)}}^{\,r_s-a_{is}}.
                                                               \tag{1}
\]
All \(Z_i(n)\) have the same modulus. Divide them by their true
common Gaussian gcd \(D(n)\), choosing any associate, and write
\(P_i(n)=Z_i(n)/D(n)\) and \(|P_i(n)|=R_n\).
Assume at least two \(P_i(n)\) are distinct for all sufficiently
large \(n\) on a fixed parity subsequence and that these points
lie in an arc of length at most \(C\sqrt{R_n}\), with fixed
\(C>0\). The family and its prefactors remain fixed as \(n\to\infty\).
For the count below, keep only the distinct eventual points.

Let \(C_e(n)\) be the exact Gaussian cyclotomic layers in the
linked reduction, with total multiplicity \(M_e\) and row count
\[
A_{i,e}=\sum_{s:e\mid m_s}a_{is}.
                                                               \tag{2}
\]
Write \(W_e=\max_i A_{i,e}-\min_i A_{i,e}\) and
\(L=\sum_e\phi(e)W_e\). **Assume \(L>0\).** The effective
primitive-radius formula (12) in the linked reduction then gives
\(R_n\to\infty\). Distinct fixed prefactors alone do not imply
this: if \(L=0\), the primitive points can lie on a stationary
circle, which is outside the theorem. The arc's angular width is
\(O(R_n^{-1/2})\), so every ratio \(P_i(n)/P_0(n)\) tends to
\(1\). The common divisor cancels from this ratio.

Set \(F=1+2i\) and
\[
\zeta=\frac{\alpha+\beta}{\alpha-\beta}
     =\frac F{\sqrt5},\qquad
\zeta^2=\frac F{\overline F}.
                                                               \tag{3}
\]
The leading coefficient of \(C_1\) relative to
\(\overline{C_1}\) is \(\zeta\); all layers \(e>1\) have the
same positive leading coefficient in both orientations. Hence
\[
\lim_{n\to\infty}\frac{Z_i(n)}{Z_0(n)}
=\frac{u_i}{u_0}\,
  \zeta^{A_{i,1}-A_{0,1}}=1.                              \tag{4}
\]
This is an exact equality of the limiting complex numbers, not
merely an equality of absolute values. Because \(u_i/u_0\in
\mathbf Q(i)\), while odd powers of \(\zeta\) are outside
\(\mathbf Q(i)\), all differences \(A_{i,1}-A_{0,1}\) are even.
With \(k_i=(A_{i,1}-A_{0,1})/2\in\mathbf Z\), (4) forces
\[
\boxed{\displaystyle
\frac{u_i}{u_0}
=\left(\frac{\overline F}{F}\right)^{k_i}.}             \tag{5}
\]
Thus the apparently arbitrary row-dependent prefactors are
determined by the base-layer orientation, up to one common
prefactor. There is no extra row-dependent Gaussian unit:
such a unit would change the limiting ratio in (4).
The argument is unchanged by which associate of \(D(n)\) is used.

## 2. Exact affine map to prime allocations

For each rational split prime \(p\) dividing the common squared
modulus of the \(P_i(n)\), choose one Gaussian prime \(\pi_p\)
above it. Let \(v_{\pi_p}\) be its valuation. The exact layer
factorization gives
\[
Z_i=u_i\prod_e C_e^{A_{i,e}}\overline C_e^{\,M_e-A_{i,e}}.
                                                               \tag{6}
\]
Combining (5)--(6), for every such \(p\) and each sufficiently
large \(n\),
\[
\begin{aligned}
v_{\pi_p}(P_i)
={}&c_p(n)
+\sum_e A_{i,e}
       \bigl(v_{\pi_p}(C_e)-v_{\pi_p}(\overline C_e)\bigr)\\
&+\frac{A_{i,1}-A_{0,1}}2
       \bigl(v_{\pi_p}(\overline F)-v_{\pi_p}(F)\bigr),\\
c_p(n)
={}&v_{\pi_p}(u_0)
 +\sum_e M_e v_{\pi_p}(\overline C_e)-v_{\pi_p}(D(n)).
\end{aligned}                                                \tag{7}
\]
The last coefficient is zero unless \(p=5\), where it is
\(+1\) or \(-1\) according to the chosen orientation. The
half in (7) is integral on the actual rows by the parity
already proved.

Equation (7) is an **exact affine-linear map** from the layer
row \((A_{i,e})_e\) to the split-prime allocation row. A
rational prime may divide several \(C_e\), or both \(C_e\)
and a prefactor; valuations still add, so bounded layer
overlaps do not disturb this map. Dividing by the true common
gcd subtracts the same valuation from every row. Gaussian
units have zero valuations. Inert and ramified prime powers
are common to the equal-norm rows and are not allocation
coordinates. Consequently, any real affine relation among
the layer rows would induce the same relation among the
prime-allocation rows.

## 3. Low jets and the actual small-arc count

Normalize the layer counts by subtracting \(\min_i A_{i,e}\).
The resulting orientation polynomials \(\mathcal A_i(z)\)
have degree \(L\) and fixed norm product, as in equations
(16)--(20) of the linked reduction. Let
\[
h=\min_{i\ne j}\operatorname{ord}_{z=0}
       \bigl(\mathcal A_i-\mathcal A_j\bigr).
                                                               \tag{8}
\]
By (5), identical layer counts would give identical rows, so
the polynomials of distinct eventual points are distinct.
Thus \(h\) is finite and odd, and the endpoint
condition gives \(L\leq4h\). The shared coefficients below
\(h\) imply that every difference of layer rows lies in
the real kernel of the integer matrix
\[
\mathsf M_{d,e}
=\begin{cases}
\mu(e/d),&d\mid e,\\
0,&d\nmid e,
\end{cases}
\quad d<h\text{ odd},\quad e\text{ active odd layer}.
                                                               \tag{9}
\]
Indeed, the logarithmic odd moments vanish below \(h\);
their Möbius inversion is equation (20) of the reduction.
Write \(\rho=\dim_{\mathbf R}\ker\mathsf M\), equal to its
rational nullity.

The geometric affine-rigidity theorem cited above says that
distinct equal-modulus Gaussian points on an arc of length at
most \(\sqrt2\sqrt R\) have affinely independent split-prime
allocation rows. Its inequality is strict even at \(C=\sqrt2\).
By (7), the corresponding layer rows are therefore affinely
independent as well. They lie in one affine translate of
\(\ker\mathsf M\), so
\[
\boxed{\displaystyle
M\leq\rho+1\qquad(C\leq\sqrt2).}                         \tag{10}
\]
This conclusion is about **actual points** in the small arc.
The polynomial fibre itself can have more orientation vectors:
affine independence is supplied by the circle geometry, not
by the low-jet equations alone.

For arbitrary fixed \(C>0\), partition the original arc into
\(J=\max(1,\lceil C/\sqrt2\rceil)\) consecutive arcs of length
at most \(\sqrt2\sqrt{R_n}\). Each subfamily retains the
same affine map (7) and lies in the same affine translate of
\(\ker\mathsf M\). Applying (10) separately gives
\[
\boxed{\displaystyle M\leq J(\rho+1).}                \tag{11}
\]
Assign shared subdivision endpoints to one part. This
partition does not assert affine independence across parts.

## 4. Uniform model bounds and the remaining rank problem

Projection from \(\ker\mathsf M\) to its coordinates \(e\geq h\)
is injective: a nonzero kernel vector supported entirely below
\(h\) has a largest nonzero coordinate \(d<h\), and then its
\(d\)-th row of (9) equals that coordinate. Thus
\[
\rho\leq\#\{e\geq h:W_e>0\}.                            \tag{12}
\]
If all active indices have at most \(r\) distinct prime
divisors, define
\[
c_r=\prod_{j=1}^{r}\left(1-\frac1{p_j}\right),
\qquad p_1=3,p_2=5,\ldots,\qquad c_0=1.
\]
Every active \(e\geq h\) has
\(\phi(e)\geq c_r e\geq c_r h\), while
\(\sum_e\phi(e)W_e=L\leq4h\). Therefore
\[
\rho\leq\left\lfloor\frac4{c_r}\right\rfloor,\qquad
M\leq J\left(1+\left\lfloor\frac4{c_r}\right\rfloor\right).
                                                               \tag{13}
\]
For at most two distinct prime divisors per active index,
\(c_2=8/15\), and (13) gives \(M\leq8J\).
These constants are independent of the growing fixed template and
its layer indices. The point at which its asymptotic arc
condition and distinctness become valid may depend on that
template. Under these growing endpoint hypotheses, the
separate prime-power fibre theorem in Section 7 of the
reduction gives the stronger four-row bound even without
the small-arc argument.

There is also a degree-dependent bound without a restriction
on prime divisors. For odd \(e\geq3\), the elementary
estimate proved in Section 8 of the reduction gives
\[
\phi(e)\geq\frac{e}
 {\sqrt{\exp(1)\log_3 e}}.
\]
For \(h\geq3\), its right-hand side increases with \(e\),
so (12) and \(L\leq4h\) yield
\[
\rho\leq
\left\lfloor4\sqrt{\exp(1)\log_3h}\right\rfloor
\leq
\left\lfloor4\sqrt{\exp(1)\log_3(L/2)}\right\rfloor.
                                                               \tag{14}
\]
The second inequality uses reciprocity of a nonzero
degree-\(L\) difference, which gives \(h\leq L/2\).
Thus (11) is \(O(J\sqrt{\log L})\), or
\(O(C\sqrt{\log L})\) for \(C\geq\sqrt2\), **inside this
growing fixed-template model**. This is not a bound stated in terms
of a general circle radius. When \(h=1\), \(L\leq4\);
the only possible active supports have at most two indices,
so \(\rho\leq2\) and (11) gives \(M\leq3J\).

The proposed universal estimate \(\rho\leq3\) is false: the exact
rank-four support in
[the nullity counterexample note](cyclotomic_nullity_four_counterexample.md)
satisfies the same degree budget, even with squarefree divisor-closed
support. The support budget alone does not produce an endpoint family;
a separate [five-row construction](cyclotomic_five_row_subendpoint_family.md)
now verifies the actual widths, parity and Gaussian realization.
The bounded-prime-divisor theorem and the prime-power
four-row theorem remain valid in their stated scopes.

The dependency-free
[weighted-nullity checker](check_cyclotomic_weighted_nullity.py)
exhaustively enumerates every odd support within the minimum
degree budget \(4h\), charging the active base layer two,
for odd \(3\leq h\leq31\) with `--full`. It finds modular
nullity at most three (two at \(h=3\)). Rank modulo \(101\)
cannot exceed rational rank, so this certifies the rational bound in
this **finite range**; it does not control larger thresholds. The
default checker run stops at \(h=19\) for speed.

## 5. An affine-dimension bound for every whole polynomial fibre

The low-jet equations also control the number of polynomials in a
whole fibre, without assuming that the rows are realized by points on
one circle. Let \(\mathcal F\) be any family of distinct orientation
polynomials (16) with common widths \(W_e\), even base-layer counts,
and pairwise contact order at least \(h\), where \(L\leq4h\). Let
\(\rho=\dim_{\mathbb R}\ker\mathsf M\) as in (9). Then
\[
\boxed{|\mathcal F|\leq2\rho+1.}                         \tag{15}
\]

To prove this, omit coordinates with \(W_e=0\) and map each row to
\[
 f_i(e)=\sqrt{\phi(e)W_e}
       \left(\frac{2a_{i,e}}{W_e}-1\right).
                                                               \tag{16}
\]
This is an injective affine map on the active layer coordinates, so
the \(f_i\) are distinct. Their affine hull has dimension \(r\leq\rho\),
since every difference of layer rows is in \(\ker\mathsf M\).

For a pair \(i\ne j\), put
\[
\delta_{ij}=\sum_e\phi(e)|a_{i,e}-a_{j,e}|.
\]
Cancel the common cyclotomic factors from \(\mathcal A_i\) and
\(\mathcal A_j\). The remaining pair is \(U(z),U(-z)\), where both
polynomials are reciprocal of the same even degree
\(d=\delta_{ij}\), and the canceled factor has constant coefficient
one. Their difference therefore has the same contact order as
\(\mathcal A_i-\mathcal A_j\). A nonzero reciprocal difference with
first nonzero term at order \(h_{ij}\) also has a nonzero term at
degree \(d-h_{ij}\); hence \(d\geq2h_{ij}\geq2h\). The difference
is odd under \(z\mapsto-z\), so its first nonzero order is odd.
The elementary inequality
\[
|x-y|\leq x+y-2xy\qquad(0\leq x,y\leq1)
\]
applied in each active coordinate gives
\[
\langle f_i,f_j\rangle
\leq L-2\delta_{ij}\leq L-4h\leq0.
                                                               \tag{17}
\]

Let \(p\) be the orthogonal projection of the origin onto the affine
hull of the \(f_i\), and write \(f_i=p+v_i\) with \(v_i\) in its
parallel linear space of dimension \(r\). If \(p\ne0\), then
\(\langle v_i,v_j\rangle\leq-\|p\|^2<0\) for distinct rows. Such a
set is affinely independent: a nontrivial affine relation splits into
positive and negative coefficients, and the squared norm of the
common positive/negative sum would be a sum of strictly negative
cross inner products. Thus \(|\mathcal F|\leq r+1\).

If \(p=0\), the vectors themselves have pairwise nonpositive inner
products. There is at most one zero vector. For the \(m\) nonzero
vectors, normalize to unit vectors \(u_i\) and let \(G\) be their
Gram matrix. Its rank is at most \(r\), its diagonal entries are
one, and its off-diagonal entries \(g_{ij}\) lie in \([-1,0]\).
Since \(\|\sum_i u_i\|^2\geq0\),
\[
 -\sum_{i\ne j}g_{ij}\leq m.
\]
Using \(g_{ij}^2\leq-g_{ij}\), we obtain
\[
\operatorname{tr}(G^2)=m+\sum_{i\ne j}g_{ij}^2\leq2m.
\]
Eigenvalue Cauchy--Schwarz now yields
\[
m^2=(\operatorname{tr}G)^2
\leq\operatorname{rank}(G)\operatorname{tr}(G^2)
\leq2rm,
\]
so \(m\leq2r\). Including a possible zero vector gives
\(|\mathcal F|\leq2r+1\leq2\rho+1\), proving (15).

When \(L<4h\), equation (17) is strictly negative for every pair.
The affine-relation argument applies directly to the \(f_i\), giving
the stronger bound \(|\mathcal F|\leq\rho+1\).

Combining (15) with (13), when every active index has at most \(s\)
distinct prime divisors,
\[
|\mathcal F|\leq2\left\lfloor\frac4{c_s}\right\rfloor+1.
\]
For \(s=2\), this gives \(|\mathcal F|\leq15\), improving the
previous 128-row whole-fibre estimate. Without a bound on prime
divisors, (14) gives \(|\mathcal F|=O(\sqrt{\log L})\). These are
whole-fibre bounds in the cyclotomic template, with no subdivision
factor \(J\); they are still degree-dependent in the unrestricted
case and do not imply a uniform bound over arbitrary lattice circles.

There is now a degree-independent result for all active odd indices:
[prime-power compression and an intersection graph](cyclotomic_uniform_rank.md)
give \(\rho\leq1125\), hence at most 2251 polynomials, or 1126
when \(L<4h\). Arbitrary multiplicities remain allowed. This supersedes
the degree-dependent estimates above within the one-class model.
Nonproportional affine classes and arbitrary circles remain open.
