# Mixed affine classes: coupled phase equations and a rank limitation

This note records the exact low-frequency equations when several
nonproportional affine Fibonacci classes are present. Their equations
couple different classes through powers of the golden ratio. A
class-by-class Möbius rank bound therefore does not follow. The note
gives (i) an exact three-class cancellation, (ii) a large rational
kernel under the height budget, and (iii) a proof that this kernel
count does not itself produce endpoint rows in the width-one example.
No reduction of arbitrary lattice circles to affine Fibonacci
templates is asserted.

**Continuation:** bounded integer differences now give uniform counts
for [common-slope classes](common_slope_cyclotomic_uniform_count.md)
and [distinct slopes with distinct 2-adic valuations](frequency_separated_affine_uniform_count.md).
These results overcome the rational-rank obstruction within those
families. The unbounded rational kernel below remains a valid example;
the new proofs constrain only the admissible bounded integer differences.

## 1. Mixed frequency expansion

Fix finitely many nonproportional affine classes
\[
d_c(n)=A_c n+B_c>0,
\qquad A_c>0\text{ even},\quad B_c\text{ odd},
\]
and restrict to a parity subsequence of \(n\). Let
\(\varphi=(1+\sqrt5)/2\), \(\xi=\varphi^{-n}\), and
\[
\gamma_c=i^{B_c}(-1)^{A_cn/2}\varphi^{-B_c},
\qquad z_c=\gamma_c\xi^{A_c}
      =i^{d_c(n)}\varphi^{-d_c(n)}.                    \tag{1}
\]
On the chosen subsequence \(\gamma_c\) is constant. For odd layer
indices, put \(\Psi_1(z)=1-z\) and \(\Psi_e(z)=\Phi_e(1,z)\)
for \(e>1\). After removing the explicit common layer factors,
the orientation polynomial of row \(i\) is
\[
\mathcal P_i(\xi)=\prod_{c,e}
 \Psi_e(z_c)^{a_{i,c,e}}\Psi_e(-z_c)^{W_{c,e}-a_{i,c,e}},
\quad 0\leq a_{i,c,e}\leq W_{c,e}.                    \tag{2}
\]
Here \(W_{c,e}\) is the width in class \(c\), and the effective
height slope is
\[
\mathcal L=\sum_c A_c\sum_e\phi(e)W_{c,e}.             \tag{3}
\]
For the circle application assume fixed equal-norm Gaussian prefactors,
positive effective slope \(\mathcal L>0\), and distinct eventual rows
whose limiting phases coincide. Phase alignment gives common total
base parity, and cancels the fixed prefactor from each normalized
phase ratio. These hypotheses exclude stationary templates made from
arbitrary fixed prefactors.
The primitive radius satisfies
\(\log R_n=\tfrac12\mathcal L n\log\varphi+O(1)\).

For a pair of rows, set
\(b_{c,e}=a_{i,c,e}-a_{j,c,e}\) and
\[
S_{c,k}=\sum_e b_{c,e}c_e(k),
\]
where \(c_e(k)\) is the Ramanujan sum. The odd logarithmic phase
expansion is
\[
\log\frac{\mathcal P_i}{\mathcal P_j}
=-2\sum_c\sum_{k\geq1,\ k\,\mathrm{odd}}
       \frac{S_{c,k}}{k}\,\gamma_c^k\xi^{A_ck}.          \tag{4}
\]
Thus the coefficient at frequency \(t\) is
\[
-2\sum_{\substack{c:A_c\mid t\\t/A_c\text{ odd}}}
   \frac{S_{c,t/A_c}}{t/A_c}\,\gamma_c^{t/A_c}.       \tag{5}
\]
All terms in (5) lie in \(i\mathbb Q(\sqrt5)\), a two-dimensional
vector space over \(\mathbb Q\). Vanishing at one frequency is
therefore at most two rational linear equations, and it is generally
a coupled sum over classes. For distinct eventual rows whose limiting
phases are aligned, let \(\tau_{ij}\) be the first nonzero frequency
in (4). The fixed-endpoint arc condition of size \(C\sqrt{R_n}\) forces
\[
\mathcal L\leq4\tau_{ij}.                              \tag{6}
\]
The frequency \(\tau_{ij}\) is even because every \(A_c\) is even;
there is no general odd-contact assertion in the mixed setting.

## 2. A three-class cancellation with a genuine endpoint pair

Take three nonproportional classes
\[
(A_c,B_c)=(2,1),(2,3),(2,5),
\]
with base-layer multiplicities \(W_{c,1}=2,6,2\). Let row 0 use
the conjugate of each Gaussian Fibonacci factor and row 1 use the
unconjugated factor, with the indicated multiplicity. Their total
base-layer difference is \(2+6+2=10\), which is even. The Gaussian
prefactors
\[
u_0=F^5\overline F^{5},\qquad
u_1=\overline F^{10},\qquad F=1+2i,
\]
have equal norm and ratio
\(u_1/u_0=(\overline F/F)^5\). This cancels the limiting phase
\(\zeta^{10}\), where \(\zeta=F/\sqrt5\).

After this cancellation the exact phase ratio is
\[
\prod_{c=0}^2\left(\frac{1-z_c}{1+z_c}\right)^{r_c},
\qquad (r_0,r_1,r_2)=(2,6,2).
\]
On a fixed parity subsequence, its first odd logarithmic coefficient
is proportional to
\[
2(\varphi^{-1}-3\varphi^{-3}+\varphi^{-5})=0,
\]
since \(\varphi^4-3\varphi^2+1=0\). The next odd coefficient is
proportional to
\[
2(\varphi^{-3}-3\varphi^{-9}+\varphi^{-15})\ne0;
\]
indeed the parenthesis equals
\(\varphi^{-9}(\varphi^6-3)+\varphi^{-15}>0\). Hence the first
frequency is \(\tau=2\cdot3=6\). The height slope is
\[
\mathcal L=2(2+6+2)=20<24=4\tau.
\]
The two primitive points are distinct for all sufficiently large
indices and satisfy
\(\operatorname{diam}(P_0,P_1)/\sqrt{R_n}=O(\varphi^{-n})\to0\).
The cross-class determinants
\(A_cB_{c'}-A_{c'}B_c=2(B_{c'}-B_c)\) are nonzero constants, so these
are genuinely distinct affine classes. This exact pair shows that
low-frequency cancellation need not hold class by class.

## 3. A large coupled rational kernel is not a fiber construction

For \(N\) classes take
\[
A_c=2,\qquad B_c=2c+1\quad(0\leq c<N),
\qquad e=1,\quad W_{c,1}=1.
\]
The total height slope is \(\mathcal L=2N\). At odd harmonic
\(k=2\ell+1\), the frequency is \(t=2k\), and after removing a
common nonzero scalar its coefficient equation is
\[
\sum_{c=0}^{N-1} b_c(-1)^c
        \varphi^{-k(2c+1)}=0
        \quad\text{in }\mathbb Q(\sqrt5).              \tag{7}
\]
Each \(\ell\) supplies at most two rational equations. Requiring
vanishing for \(\ell=0,\ldots,K-1\) therefore leaves a rational
kernel of dimension at least \(N-2K\). For \(N=8K\), this is at
least \(6K\), while
\[
\mathcal L=16K<16K+8
 =4\bigl(2(2K+1)\bigr).                                \tag{8}
\]
Thus the coupled linear kernel can have unbounded dimension under the
formal endpoint height budget. The total base-parity condition
\(\sum_c b_c\) even is only a parity restriction, not an additional
rational rank equation. This dimension count says nothing by itself
about the number of admissible orientation rows: in a width-one box,
pair differences must have \(b_c\in\{-1,0,1\}\).

In fact no nonzero such bounded vector cancels even the first
frequency. Put \(q=-\varphi^{-2}\), so \(|q|<1/2\). Equation (7) at
\(k=1\) is a nonzero scalar times
\(\sum_c b_cq^c=0\). If \(c_0\) is the first nonzero coefficient,
the leading term has magnitude \(|q|^{c_0}\), whereas the remaining
tail has magnitude at most
\[
\sum_{c>c_0}|q|^c
\leq |q|^{c_0}\frac{|q|}{1-|q|}
=|q|^{c_0}\varphi^{-1}<|q|^{c_0}.
\]
So the sum cannot vanish. For \(N>4\), the first frequency has
slope 2 and \(\mathcal L=2N>8=4\cdot2\); a pair satisfying the
contact condition for that full width budget would have to cancel it.
Thus no growing endpoint family with at least two distinct rows can
have all \(N>4\) of these width-one coordinates active, despite the
large rational kernel. A smaller subfamily may have fewer active
coordinates after its own common factor is removed; this statement
does not assign it the original \(2N\) effective slope.

## 4. The coupled kernel still admits an obtuse-fiber bound

For fixed widths, let \(V_{<\mathcal L/4}\) be the rational vector
space of all layer-difference arrays \(b=(b_{c,e})\) for which the
coefficient (5) vanishes at every frequency \(t<\mathcal L/4\). Write
\(\rho=\dim_{\mathbb Q}V_{<\mathcal L/4}\). An endpoint-compatible
family has every pair difference in this kernel, by (6). This is the
full coupled kernel across all classes; it is not the sum of separate
classwise kernels.

There is a useful bound once the dimension \(\rho\) of the full
rational coupled low-frequency kernel is known. Let \(\mathcal F\)
be a family of distinct rows with common widths, common total
base-layer parity, positive \(\mathcal L\), and pairwise endpoint
contacts \(\tau_{ij}\) satisfying (6). For a pair, cancel common
factors class by class. The remaining polynomials \(U(\xi),V(\xi)\)
have the same degree
\[
\delta_{ij}=\sum_{c,e}A_c\phi(e)|b_{c,e}|.             \tag{9}
\]
Every \(\Psi_e(\pm z)\), \(e>1\), is reciprocal of degree
\(\phi(e)\). For \(e=1\), the two leading coefficients differ by
\(-1\); their total sign in the pair is
\((-1)^{\sum_c|b_{c,1}|}=1\), since the total base difference is
even. Thus \(U,V\) have the same leading coefficient. Let
\(K=\mathbb Q(i,\sqrt5)\) and let
\(\sigma(i)=i,\ \sigma(\varphi)=-1/\varphi\), acting on the
coefficients while \(X\) is a formal variable. Since \(B_c\) is odd,
\(\sigma(\gamma_c)=\gamma_c^{-1}\). For \(e>1\), reciprocity gives
\[
\sigma_{\rm coeff}\Psi_e(\pm\gamma_c X^{A_c})
= (\pm\gamma_c)^{-\phi(e)}X^{A_c\phi(e)}
  \Psi_e(\pm\gamma_c X^{-A_c}).
\]
The analogous identity for \(e=1\) has an extra minus sign for the
plus-argument orientation. The numbers of these factors in \(U,V\)
differ by \(\sum_c b_{c,1}\), which is even by total base parity.
Their reciprocal scalars therefore agree. Consequently there is a common nonzero scalar
\(\Lambda_{ij}\) such that
\[
\sigma_{\rm coeff}U(X)=\Lambda_{ij}X^{\delta_{ij}}U(X^{-1}),\qquad
\sigma_{\rm coeff}V(X)=\Lambda_{ij}X^{\delta_{ij}}V(X^{-1}).
\]
The residual multiplicities \(|b_{c,e}|\) agree in both products,
and the preceding parity calculation handles their signs. Therefore \(U-V\) has coefficient
symmetry about degree \(\delta_{ij}/2\).
If its first nonzero phase term is at \(\tau_{ij}\), it also has a
nonzero term at degree \(\delta_{ij}-\tau_{ij}\); hence
\[
\delta_{ij}\geq2\tau_{ij}\geq\mathcal L/2.           \tag{10}
\]

Map each row to active layer coordinates
\[
f_i(c,e)=\sqrt{A_c\phi(e)W_{c,e}}
       \left(2a_{i,c,e}/W_{c,e}-1\right).
\]
The coupled low-frequency equations place all these vectors in an
affine space of dimension at most \(\rho\). The same inequality
\(|x-y|\leq x+y-2xy\) used in the one-class obtuse embedding gives
\(\langle f_i,f_j\rangle\leq\mathcal L-2\delta_{ij}\leq0\).
The affine-projection and Gram-matrix argument in
[the one-class affine-rank note](cyclotomic_affine_rank_reduction.md#5-an-affine-dimension-bound-for-every-whole-polynomial-fibre)
therefore gives
\[
|\mathcal F|\leq2\rho+1.                              \tag{11}
\]
The mixed example (7) shows why (11) alone does not settle the
problem: under the endpoint budget, the full rational kernel dimension
need not be bounded independently of the number of affine classes.
The width-one example also shows that a large kernel need not contain
any admissible nonzero row difference. A useful mixed-class theorem
must control this integer-box intersection, not just rational rank.

All these statements concern fixed finite affine templates along a
growing parity subsequence. They imply no bound for arbitrary
origin-centered lattice circles.

The common-slope continuation separates all sufficiently large canceled
harmonics by conjugating the quadratic unit, then recovers the remaining
small moments modulo three primes. Its colored compression theorem bounds
the affine dimension of the actual rows, even though the full coupled
rational kernel is unbounded. Distinct slopes with distinct 2-adic
valuations have disjoint frequency sets, permitting a weighted extension.
Overlapping slope frequencies remain outside those results.
