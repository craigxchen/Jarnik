# Arbitrary-unit addendum to Möbius conductor retention

This note audits the extension of
[the erased-conductor theorem](mobius_erased_conductor_budget.md) when the
primitive source points are not first restricted to one odd Gaussian-unit
class. The extension is valid. For \(m\ge5\), if the full primitive source
and its full primitive image both have endpoint constants at most two, then
\[
 \boxed{N'\ge\frac1{16}N^{\,1-4/m}.}                     \tag{1}
\]
Neither side is subjected to a unit-class extraction, and neither least
squared radius is changed.

The only structural caveat is that an even-norm half-angle anchor can change
the two-adic normalization factor \(\epsilon_j\). Thus one should not import
the constant-\(\epsilon\) all-anchor formula verbatim. The proof below fixes
one rational phase anchor \(H_0=1\) and uses only odd split-prime coefficient
intervals, relative phases, and intrinsic pair residues. Those quantities
are unaffected by the source unit classes.

## 1. Primitive half-angle rows with arbitrary units

Choose any source point as phase anchor. Every relative phase is a rational
point of the unit circle and can be written as
\[
 \frac{z_i/z_0}{1}=\frac{H_i}{\overline{H_i}}
\]
for a coordinate-primitive Gaussian integer \(H_i\), unique up to sign.
Explicitly, for a relative phase \(x+iy\ne-1\), clear denominators in
\((1+x)+iy\) and remove the coordinate gcd; for the phase \(-1\), use
\(H_i=i\).
Take \(H_0=1\). A primitive \(H_i\) contains \(1+i\) to exponent zero or
one. Let
\[
 \sigma=\#\{i:N(H_i)\text{ is even}\}.
\]
Since \(H_0\) is odd, \(0\le\sigma\le m-1\).

For each pair let \(f_{ij}\) be the full primitive relative half-angle norm,
and put
\[
 P_{\rm full}=\prod_{i<j}f_{ij}.
\]
The two-adic exponent of \(f_{ij}\) is one exactly when one of \(H_i,H_j\)
has even norm and the other has odd norm. Therefore
\[
 P_{\rm full}=2^{\sigma(m-\sigma)}P_{\rm odd},            \tag{2}
\]
where \(P_{\rm odd}\) is the product of the odd parts of all pair norms.
This is the complete source-unit cost.

The chord identity is intrinsic to the pair phase and does not require a
common unit class. Distinct source points in an endpoint arc of constant
\(C\le2\) satisfy
\[
 f_{ij}\ge\frac4{C^2}\sqrt N\ge\sqrt N.
\]
Writing \(E=\binom m2\) and \(W=\log N\), multiplication gives
\[
 \log P_{\rm odd}
 \ge\frac{m(m-1)}4W-\sigma(m-\sigma)\log2.                \tag{3}
\]

## 2. Odd-prime deficit and erased layers

The least primitive Gaussian-circle squared radius \(N\) is odd: if it
were even, every point would be divisible by \(1+i\), and dividing all
points by that common factor would give a smaller integral circle with the
same relative phases. Inert prime factors are similarly excluded by common
Gaussian division. This remains true when the half-angle parities mix;
their ramified factors belong to pair numerators, not to the least circle
radius. Thus \(W=\log N\) is entirely supported on odd split primes.

At every odd split prime, use the signed valuations of the primitive
\(H_i\). The exact layer-cake deficit is
\[
 D_{\rm odd}
 =\frac{m^2}{4}W-\log P_{\rm odd}
 =\sum_{\text{odd source layers}}(n-m/2)^2\log p.         \tag{4}
\]
In particular, (3) gives the safe bound
\[
 D_{\rm odd}\le
 \frac m4W+\sigma(m-\sigma)\log2.                         \tag{5}
\]
Clip these odd signed valuations to the coefficient intervals exactly as
in the fixed-unit proof. Let \(W_{\rm clip}\) be the resulting clipped
width and
\[
 E_{\rm odd}=W-W_{\rm clip}.
\]
The clipped conductor divides \(N'\), so \(W_{\rm clip}\le W'\).
At each erased odd layer, if \(n\) rows lie beyond the relevant interval
endpoint, their \(\binom n2\) pairs contribute to the target primitive
residue product \(B'\). Thus, for
\[
 c(A)=\min_{1\le n\le m-1}
 \left\{\binom n2+A(n-m/2)^2\right\},
\]
the unchanged odd-prime argument gives
\[
 c(A)E_{\rm odd}\le\log B'+AD_{\rm odd}.                 \tag{6}
\]
Rows at an interval endpoint contribute zero; disjoint intervals and
primes outside \(N\) add only nonnegative target overlap.

Using \(E_{\rm odd}\ge W-W'\) in (6), followed by (5), gives directly
\[
 (c(A)-Am/4)W
 \le(c(A)+m/8)W'
 +\eta_{\rm target}
 +A\sigma(m-\sigma)\log2.                                \tag{7}
\]

## 3. Optimized source and target parity constants

Take
\[
 k=\lfloor m/4\rfloor,\qquad L=m-2k-1,\qquad
 A_m=\frac{k}{L},
\]
\[
 c_m=\frac{k[(m-1)L+1]}{4L},\qquad
 G_m=c_m+\frac m8.
\]
As proved and independently audited in the fixed-unit argument, this gives
\[
 \rho_m=\frac{c_m-A_mm/4}{G_m}
 =\frac{2k(m-2k-2)}{m+2k(m-2k-2)}
 \ge1-\frac4m.                                           \tag{8}
\]

Let \(r\) of the \(m\) primitive target half-angle columns have even norm.
The full target parity-cut estimate is
\[
 \log B'\le\frac m8W'
 +E\log(C'/2)+\frac{r(m-r)}2\log2.                       \tag{9}
\]
For \(C'\le2\), equations (7), (8)--(9) yield
\[
 W'\ge\rho_mW-
 \frac{\frac12r(m-r)+A_m\sigma(m-\sigma)}{G_m}\log2.     \tag{10}
\]

The additive loss is at most \(4\log2\). Indeed, put
\(q_m=\lfloor m^2/4\rfloor\). Then both
\(r(m-r)\) and \(\sigma(m-\sigma)\) are at most \(q_m\), while
\[
 q_m(A_m+1/2)\le4G_m.                                    \tag{11}
\]
For a direct verification, write \(m=4k+s\), \(0\le s\le3\), and set
\[
 H=m+2k(m-2k-2).
\]
The identity
\[
 4H-m^2=s(4-s)
\]
implies \(H\ge q_m\). Also
\[
 A_m+\frac12=\frac{m-1}{2L},\qquad
 4G_m=\frac{(m-1)H}{2L},
\]
which proves (11). Combining (8), (10), and (11) gives
\[
 W'\ge(1-4/m)W-4\log2.
\]
Exponentiation proves (1). With specified \(\sigma,r\), equation (10) is
the sharper exact result.

## 4. Full determinant-core parity with arbitrary source units

The same source parity must be retained in the two-adic determinant-core
lemma. Let
\[
 c=\frac{\epsilon}{4}|b-a|,\qquad
 F_m=\left\lfloor\frac{(m-1)^2}{4}\right\rfloor,
\]
using one fixed primitive integral normalization with anchor \(H_0=1\).
In this normalization the common coefficient scale
\(k_0=\gcd(Nu,Nv)\) is odd: if \(2\mid k_0\), then the coprime Gaussian
coefficients \(u,v\) would both be divisible by \(1+i\). Consequently the
determinant scalar has \(v_2(\delta)=v_2(c)=\ell\), with no hidden
coefficient contribution at two.

At two, let \(o_i\) and \(o_i'\) be the source and target even-norm
indicators. Then
\[
 v_2(e_{ij})=o_io_j,\qquad v_2(e'_{ij})=o_i'o_j'.
\]
If \(\rho_i\) are the raw ordinary row-content valuations and
\(\ell=v_2(c)\), a primitive matrix row gives
\[
 \beta_{ij}+v_2(e_{ij})\ge\min(\rho_i,\rho_j).
\]
Combining this with the determinant identity gives the conservative
pair inequality
\[
 \beta_{ij}+\beta'_{ij}
 \ge\ell-|\rho_i-\rho_j|-o_io_j-o_i'o_j'.                \tag{12}
\]
After summing pairs and applying the content layer-cake bound,
\[
 \boxed{\displaystyle
 c^{F_m}\mid
 2^{\binom{\sigma}{2}+\binom r2}B_*B'_* .}               \tag{13}
\]
Odd split and inert core primes have no parity loss.

The full source and target endpoint product estimates are
\[
\begin{aligned}
 B_*&\le(C/2)^E2^{\sigma(m-\sigma)/2}N^{m/8},\\
 B'_*&\le(C'/2)^E2^{r(m-r)/2}(N')^{m/8}.
\end{aligned}
\]
Since
\[
 \binom x2+\frac{x(m-x)}2=\frac{x(m-1)}2,
\]
equation (13) gives
\[
 \boxed{\displaystyle
 c^{F_m}\le
 \left(\frac{CC'}4\right)^E
 2^{(\sigma+r)(m-1)/2}(NN')^{m/8}.}                      \tag{14}
\]
For \(C,C'\le2\), the uniform version is
\[
 \boxed{\displaystyle
 c\le
 2^{m(m-1)/F_m}(NN')^{m/(8F_m)}.}                        \tag{15}
\]
Thus arbitrary source units change only an explicit \(m\)-dependent
factor and do not change the radius exponent in the determinant-core
bound.

## 5. Scope

The full arbitrary-unit retention theorem (1) is valid because every
integer-circle relative phase has a primitive Gaussian half-angle row and
because all uses of units reduce to the exact binary parity cuts (2) and
(9). No fairness or target extraction is hidden.

The theorem remains a near-unit comparison, not a radius-independent point
bound. Equation (15) controls the determinant difference \(a-b\), not the
coefficient sum \(a+b\) or every anchor scale \(k_j\). The possible change
of \(\epsilon_j\) at even-norm anchors is why no stronger arbitrary-unit
all-anchor invariant is asserted here.

The persistent [joint arithmetic checker](check_mobius_joint_anchor_residues.py)
now includes mixed-parity source rows. Its 2,187 matrix/configuration cases
check the full source-and-target parity charge, odd pair-product deficit,
clipped conductor, and the explicit cut constants. The geometric consequence
for the full point count is proved in
[the condition-number theorem, Section 7](mobius_endpoint_condition_number.md).
