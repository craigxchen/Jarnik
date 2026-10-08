# Uniform counts for divisibility-chain families of affine slopes

The uniform count extends to arbitrarily many affine slopes whose
distinct values in each fixed 2-adic valuation bucket form a divisibility
chain. Slopes in one bucket have overlapping odd-harmonic frequencies.
An integral packet
projection moves the smallest slope upwards while preserving each
original class's low phase coefficients and not increasing its width
cost. All losses of injectivity are bounded together in the original
coordinates, rather than charged once per projection.

This is a theorem about fixed growing Fibonacci templates. It does not
classify arbitrary circle configurations or give a uniform asymptotic
entry index for varying templates.

## 1. Statement and original geometry

Use the normalization in
[the mixed affine reduction](mixed_affine_cyclotomic_rank_scope.md).
Fix finitely many nonproportional affine class representatives

\[
d_c(n)=A_c n+B_c,\qquad A_c>0\text{ even},\quad B_c\text{ odd}.
                                                               \tag{1}
\]

For each value of \(v_2(A_c)\), assume its distinct slope values can
be ordered as \(A_1\mid A_2\mid\cdots\mid A_m\). The number of
buckets and the lengths and ratios of the chains are unrestricted.
This includes arbitrary odd-prime-power chains, arbitrary divisibility
chains of one 2-adic valuation, and the previously treated case of one
slope per 2-adic valuation.
Section 6 also proves the same bound when each bucket has at most two
distinct slopes, even when they are incomparable under divisibility.
Buckets satisfying either condition may be combined.

All data are fixed, and \(n\) runs through a fixed parity subsequence.
Retain the original class labels throughout the argument. Write the
integer orientation rows as \(a_{i,c,e}\), for odd layer indices
\(e\), normalized to minimum zero in each coordinate, with actual
widths \(W_{c,e}\). Omit zero-width coordinates. Put

\[
L=\sum_{c,e}A_c\phi(e)W_{c,e}>0,
\qquad D_0=2^{32}.                                     \tag{2}
\]

Assume fixed equal-norm Gaussian prefactors, aligned limiting phases,
and distinct eventual rows on a fixed \(C\sqrt{R_n}\) arc after
division by the entire Gaussian gcd. A family of at most one row is
trivial. Otherwise let \(\tau\) be the smallest first nonzero pair
frequency. The original endpoint geometry implies

\[
L\le4\tau.                                             \tag{3}
\]

It also gives common total base parity. Global sigma-reciprocity and
the weighted obtuse embedding therefore give at most \(2r+1\) rows
if their original affine dimension is at most \(r\).

**Theorem.** Under these hypotheses the original affine dimension is
at most \(D_0+1125\), and the number of rows is at most

\[
2(D_0+1125)+1=8589936843.                               \tag{4}
\]

The bound is independent of chain lengths and ratios, slopes, intercepts,
layer indices, multiplicities, and the fixed arc constant.

## 2. Arithmetic separation with one excluded prime

Fix any odd prime \(p\). We first extend the common-slope arithmetic
lemma. Suppose a group of
classes shares a slope \(A\), its intercepts are distinct odd integers,
and its pair differences have coefficient budget

\[
\sum_{c,e}\phi(e)|b_{c,e}|\le B,
\qquad D_0\le B\le4h,
\]

where \(B\) is an integer and \(h\) is odd. Suppose its coupled phase
moments vanish for every odd \(k<h\) with \(p\nmid k\). Then the
individual class moments

\[
S_{c,k}=\sum_e b_{c,e}c_e(k)
\]

vanish for those same harmonics.

Here is the explicit modification of
[the common-slope proof](common_slope_cyclotomic_uniform_count.md).
Let \(\varphi=(1+\sqrt5)/2\), \(q=-\varphi^{-2}\), and write the
intercepts as \(B_c=B_*+2s_c\). At a known vanishing harmonic,

\[
P_k(q^k)=0,\qquad P_k(X)=\sum_c S_{c,k}X^{s_c},
\qquad \|P_k\|_1\le B.                                 \tag{5}
\]

The quadratic conjugate has modulus \(\varphi^{2k}\). The integer
leading-term root bound makes \(P_k\) identically zero whenever
\(\varphi^{2k}>B\). Thus all permitted harmonics at least

\[
K=\left\lceil\frac{\log B}{2\log\varphi}\right\rceil+1
\]

separate class by class.

For a smaller odd \(d<K\), with \(p\nmid d\), put

\[
n=\left\lfloor\frac{h-1}{16d}\right\rfloor.
\]

For \(B\ge2^{32}\), writing \(x=\log B\), we have
\(K\le3x\) and \(\sqrt B>384x\). The latter holds at \(2^{32}\)
since \(2^{16}>384\cdot32\), and then follows by monotonicity.
Consequently

\[
n\ge\frac{B}{64K}-2
 \ge\frac{B}{192x}-2
 >2\sqrt B-2\ge\sqrt B>K.                              \tag{6}
\]

Choose four Bertrand primes, one in each of
\((n,2n),(2n,4n),(4n,8n),(8n,16n)\). Discard \(p\) if it occurs,
and retain any three of the remaining primes \(r_1,r_2,r_3\).
They satisfy

\[
K<dr_j<h,\qquad p\nmid dr_j,\qquad
r_1r_2r_3>8n^3>B.
\]

Large-harmonic separation and the exact congruence
\(c_e(dr_j)\equiv c_e(d)\pmod{r_j}\) show that their product divides
\(S_{c,d}\). Since \(|S_{c,d}|\le B\), it is zero. This proves the
excluded-prime lemma. It uses only bounded integer differences, not
arbitrary rational vectors in the coupled kernel.

## 3. Integral packet projection for one class

Fix a class of current slope \(A\) and intercept \(B_c\). Write its
coordinates as \(a_{i,j,s}=a_{i,p^j s}\), with \(p\nmid s\),
using zero for absent coordinates. Suppose every pair difference has
\(S_k=0\) for odd \(k<h\) with \(p\nmid k\). For a difference set

\[
u_s=b_{0,s}-b_{1,s}.
\]

At a harmonic not divisible by \(p\), multiplicativity of Ramanujan
sums gives

\[
S_k=\sum_s u_s c_s(k),
\]

because \(c_p(k)=-1\) and \(c_{p^j}(k)=0\) for \(j\ge2\).
Moreover \(c_s(p^v k)=c_s(k)\) for \(p\nmid s\). Removing powers
of \(p\) from any odd \(k<h\) therefore proves

\[
\sum_s u_s c_s(k)=0\quad\text{for every odd }k<h.       \tag{7}
\]

One may now equalize the entire level-0 and level-1 slices by either
of the following two choices:

\[
(a'_{i,0,s},a'_{i,1,s})=(a_{i,0,s},a_{i,0,s})
\quad\text{for every }s,
\]

or

\[
(a'_{i,0,s},a'_{i,1,s})=(a_{i,1,s},a_{i,1,s})
\quad\text{for every }s.                              \tag{8}
\]

Leave higher levels unchanged. The second choice changes each low
moment by minus the expression in (7); the first changes it by
\(c_p(k)\) times that expression. Hence either choice preserves all
class moments below \(h\). The same choice must be made for **all**
\(s\) in this class. Coordinatewise choices could truncate \(u\)
and destroy its cancellations.

Let

\[
C_0=\sum_s\phi(s)W_{0,s},\qquad
C_1=\sum_s\phi(s)W_{1,s}.
\]

Use level 0 throughout if \(C_0\le C_1\), and level 1 throughout
otherwise. The equalized prefix costs \(p\min(C_0,C_1)\), at most
its former cost \(C_0+(p-1)C_1\).

For \(\Psi_1(z)=1-z\) and \(\Psi_e(z)=\Phi_e(1,z)\) otherwise,
the exact identities are

\[
\Psi_s(z)\Psi_{ps}(z)=\Psi_s(z^p),\qquad
\Psi_{p^j s}(z)=\Psi_{p^{j-1}s}(z^p)\quad(j\ge2).      \tag{9}
\]

Thus rescale the class to

\[
(A,B_c)\longmapsto(pA,pB_c),
\]

replace the equalized prefix by its one common slice indexed by
\(s\), and shift all levels \(j\ge2\) to \(j-1\). These are integral
coordinate copies, with no averaging denominators. The higher-level
weighted costs are unchanged, since
\(p\phi(p^{j-1}s)=\phi(p^j s)\). The total weighted width cost of
the class therefore does not increase.

The phase variable changes from \(z_c\) to \(z_c^p\): its constant
satisfies \(\gamma_{pA,pB_c}=\gamma_{A,B_c}^p\). Identity (9) shows
that the class's logarithmic phase coefficients at all physical
frequencies \(t<Ah\) are preserved. In terms of cleared moments,
an equalized packet satisfies

\[
S_{\mathrm{packet},k}=0\ (p\nmid k),\qquad
S_{\mathrm{packet},pk}=pS_{\mathrm{rescaled},k},
\]

so the slope multiplier in \(A S_{c,k}\) is preserved as well.
This statement concerns pair differences; common row-independent
factors and coordinate normalizations do not affect it.

## 4. Raising a smallest slope and terminating

At every stage retain the original colors, even after rescaling a
class. Its ratio \(B_c/A_c\) is unchanged. Original nonproportionality
therefore ensures distinct intercepts among classes at the same current
slope. Empty projected supports can be discarded.

If all surviving slopes exceed \(L/D_0\), stop. Otherwise choose a
2-adic valuation bucket that still has a slope at most \(L/D_0\),
and let \(A\) be its smallest slope. If its next distinct slope is
\(A_2\), choose an odd prime \(p\mid A_2/A\). The ratio is odd
because the two slopes have the same 2-adic valuation. If there is
no next slope in this bucket, choose \(p=3\). Set

\[
h_A=\min\{h\ge1:h\text{ odd},\ Ah\ge\tau\},\qquad
B=\lceil L/A\rceil.
\]

Every other current slope in this bucket is divisible by \(pA\).
Other buckets have disjoint odd-harmonic frequency sets, because
\(v_2(A'k')=v_2(A')\) for odd \(k'\). At a frequency \(Ak\) with
\(p\nmid k\), only the chosen smallest-slope group can contribute.
Its coupled moments therefore vanish whenever \(Ak<\tau\),
equivalently for odd \(k<h_A\). Its current integer width degree is
at most \(L/A\), while

\[
D_0\le B\le4h_A.
\]

Section 2 separates these moments class by class. Apply Section 3 to
every class at this smallest slope. All move to \(pA\), their physical
phase coefficients below \(\tau\) are preserved individually, and
the total width cost remains at most the fixed original \(L\).

The moved slope \(pA\) still divides \(A_2\) when a next slope
exists, so the divisibility-chain hypothesis is preserved. Equal
slopes are treated together at subsequent steps, but original colors
are not identified or merged.

Repeat, choosing a prime anew when necessary. Each processed original
color has its slope multiplied by at least three. There are finitely
many original colors, and a color is processed only while its current
slope is at most \(L/D_0\). Thus only finitely many steps are possible.
It is legitimate to continue above a bucket's largest original slope:
the same argument applies with \(p=3\). At termination either no
coordinate remains, or every remaining
coordinate has slope greater than \(L/D_0\). Since a nonconstant
integer coordinate of slope \(A\) costs at least \(A\), the number
of final active coordinates is then less than \(D_0\).

Intermediate rows need not retain the original total base parity or
be realized by an actual Gaussian circle. Neither property is used
in the iteration. Integrality, width cost, and the individually
preserved low physical phase coefficients are the invariants needed.

## 5. Bounding all lost dimensions once

Let \(V\) be the linear span of differences of the original rows,
and let \(T\) be the composition of the packet projections on this
space. Each choice in (8) is fixed across the family, so these maps
are linear on differences. All low-moment assertions are linear and
therefore hold on the entire span, not just on individual pairs.

The image has dimension less than \(D_0\), because it uses fewer
than \(D_0\) final coordinates. If \(v\in\ker(T|_V)\), its final
phase coefficients vanish in every original color. Individual
preservation through all steps implies that its original class
moments vanish at every physical frequency below \(\tau\). Thus

\[
\ker(T|_V)\subseteq
\bigoplus_c\ker M_c,
\quad
M_c[d,e]=\mu(e/d)\mathbf1_{d\mid e},\quad
d<h_{A_c}\text{ odd},                                 \tag{10}
\]

where the slopes and active column supports on the right are the
**original** ones. Triangular Ramanujan inversion justifies (10).

The weighted colored compression argument in
[the frequency-separated theorem, Section 4](frequency_separated_affine_uniform_count.md)
applies to this direct sum without any frequency-separation assumption.
Its input is only the original cost
\(\sum_c A_c\sum_{e\in E_c}\phi(e)\le L\le4\tau\).
After compression a high order \(e\ge h_{A_c}\) has weighted root
group size \(A_ce\ge\tau\); six sufficiently weakly intersecting
groups exceed the total budget. Consequently

\[
\sum_c\dim\ker M_c\le1125.                             \tag{11}
\]

This bounds the kernel of the **composed** map. There is no sum of
separate defect bounds over the chain length. Rank-nullity now gives
\(\dim V\le D_0+1125\). Applying the original global sigma-reciprocal
obtuse embedding proves (4).

## 6. A more general promotion criterion, including two arbitrary slopes

The only isolation property used in a promotion of slope \(A\) is
the existence of an odd prime \(p\) such that

\[
v_p(A)<v_p(B)\quad\text{for every other current slope }B
\text{ in the same 2-adic bucket}.                       \tag{12}
\]

Indeed \(p\nmid k\) then gives \(v_p(Ak)=v_p(A)<v_p(B)\),
so no such \(B\) divides \(Ak\). No requirement that \(B/A\) be
an integer is needed for this argument. Whenever \(A\le L/D_0\),
Sections 2--3 still promote \(A\) to \(pA\), preserving all the
same per-color low coefficients and width costs. If a finite sequence
of such promotions takes every remaining slope above \(L/D_0\),
the proof in Section 5 applies unchanged. A bucket with a single
slope may always use \(p=3\).

In particular, suppose a bucket has exactly two distinct slopes
\(A<B\), with their common 2-adic valuation. If \(A\le L/D_0\),
there must be an odd prime \(p\) with \(v_p(A)<v_p(B)\): otherwise
\(B\mid A\), contradicting \(A<B\). Choose that prime and promote
the smaller slope. Condition (12) holds. While two distinct slopes
remain, their least common multiple \(H\) is unchanged, since the
promoted \(pA\) still divides \(H\), and its p-valuation does not
exceed that of \(B\). Both current slopes divide this fixed \(H\).
The positive integer product of the two distinct slope values strictly
increases at every step and is bounded by \(H^2\). Thus finitely many
steps either put both slopes above \(L/D_0\) or leave one slope value.
The latter case is handled by repeated promotion with \(p=3\).
If a projected group becomes empty, the one-slope case applies earlier.

This proves (4) for at most two arbitrary slopes in each bucket, with
arbitrary class counts and intercepts at those slopes. For example,
the incomparable slopes 6 and 10 can be promoted to their common
multiple 30. The argument also permits combining such buckets with
the divisibility-chain buckets of the main statement, since odd
promotions preserve their frequency separation.

Condition (12) need not be available for three incomparable slopes.
For example, at slopes 6, 10 and 14, each prime has its minimum
valuation at least twice, so no group can be promoted by this rule.
This is a limitation of the proved operation, not a large-row
counterexample. More general incomparable slope sets and the arbitrary
circle problem remain unresolved.

The [three-class extension](three_class_affine_uniform_count.md) now
bypasses this isolation obstruction when a bucket contains at most
three original affine classes. It uses exact two-term irrationality,
so in particular covers one class at each of the slopes 6, 10 and 14.
It does not cover arbitrary numbers of intercept classes at those slopes.
