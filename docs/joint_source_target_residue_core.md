# Joint source-target residue control of the Möbius determinant core

This note proves a two-ended local packing inequality for the determinant
core left by all-source-anchor normalization. At every odd prime of that
core, source and target primitive pair residues jointly absorb a fixed
quadratic multiple of the determinant valuation. Combining this with the
endpoint product bounds first controls the odd part of the norm-difference
parameter \(a-b\). A separate ramified calculation controls its two-part
with an explicit target-parity factor, giving exponent \(O(1/m)\) for the
complete determinant core.

The result does not control \(a+b\), the primes already occurring in the
common norm factor \(k_j\), or the remaining intersecting-interval loss.
It is therefore a rigorous component of a possible radius-minimality
argument, not a uniform endpoint theorem.

## 1. The determinant core and its good odd primes

Use the all-anchor notation from
[the all-anchor audit](mobius_all_anchor_descent_audit.md). For
\[
 w=\frac vu,\qquad (u,v)=1,\qquad
 \frac{N(v)}{N(u)}=\frac ab,\qquad (a,b)=1,\quad a\ne b,
\]
the output-conformally normalized invariants at source anchor \(j\) are
\[
 T_j=\frac{\epsilon}{2}(a+b)k_j,\qquad
 K_j=\epsilon^2ab\,k_j^2,\qquad
 |\delta_j|=\frac{\epsilon}{4}|b-a|k_j.                  \tag{1}
\]
Put
\[
 c=\frac{\epsilon}{4}|b-a|,\qquad K_{\rm core}=\epsilon^2ab. \tag{2}
\]
The integrality conditions in the exact conformal normalization make \(c\)
an integer. If an odd prime \(p\mid c\), then \(p\nmid ab\), since
\((a,b)=1\). Hence every odd prime of \(c\) is a good determinant prime
for the core matrix:
\[
 p\mid c,\qquad p\nmid K_{\rm core}.                     \tag{3}
\]

At such a prime write
\[
 d=v_\pi(w),\qquad e=v_{\bar\pi}(w)
\]
when \(p=\pi\bar\pi\) is split. Since
\(v_p(Nw)=v_p(a/b)=0\), one has \(d+e=0\). Thus the coefficient interval
in the all-anchor distance formula collapses to the single integer
\[
 I_p=\{q\},\qquad q=d=-e,\qquad
 v_p(k_j)=|t_j-q|.                                       \tag{4}
\]
Rows with \(t_j=q\) retain \(p\) in the good determinant part \(D_j\);
rows away from \(q\) move the prime into \(k_j\) and hence into \(K_j\).
The lemma below accounts for both regimes simultaneously.

If \(p\equiv3\pmod4\) is inert, it cannot divide the norm of a primitive
Gaussian source column. It also cannot occur in the common norm scale
\(k_j\): an inert Gaussian prime occurring in both coefficient norms would
divide both \(u\) and \(v\), contrary to \((u,v)=1\). Thus an odd inert
prime dividing \(c\) is already a good determinant prime at every anchor.

## 2. A pairwise source-target inequality

Let \(H_1,\ldots,H_m\) be primitive Gaussian source columns and let
\(B_1,\ldots,B_m\) be their primitive target columns under the same
nonsingular rational projective map. At the fixed odd prime \(p\) define
\[
 \beta_{ij}=v_p(b_{ij}),\qquad
 \beta'_{ij}=v_p(b'_{ij}),\qquad
 \ell=v_p(c)>0.                                          \tag{5}
\]
Here \(b_{ij}\) and \(b'_{ij}\) are the full primitive source and target
pair residues: the Gaussian gcd norm has been divided out on each side.

When \(p\) is split, locally rotate the input isotropic coordinates by
\(q\) from (4). This subtracts \(q\) from every signed source valuation
and turns the coefficient interval into \(\{0\}\). It preserves all
relative source phases, all source pair residues, and the target
configuration. After the corresponding output-conformal normalization,
the local matrix has
\[
 v_p(\delta)=\ell,\qquad p\nmid K,                        \tag{6}
\]
so all four of its isotropic entries are \(p\)-adic units.

For an inert \(p\), use any actual anchor normalization. Equation (6) still
holds. Primitive Gaussian source and target norms are \(p\)-adic units, so
both Gaussian pair gcd norms below have zero \(p\)-valuation.

Let
\[
 \rho_i=v_p(r_i),\qquad 0\le\rho_i\le\ell,                \tag{7}
\]
where \(r_i\) is the ordinary coordinate content of the raw image of row
\(i\). Let
\[
 e_{ij}=N(\gcd_{\mathbb Z[i]}(H_i,H_j)),\qquad
 e'_{ij}=N(\gcd_{\mathbb Z[i]}(B_i,B_j)).
\]
The good-determinant joint pair inequality from
[the joint-content theorem](mobius_joint_content_conductor.md) is exactly
\[
 \beta_{ij}\ge
 \min(\rho_i,\rho_j)+v_p(e'_{ij}).                        \tag{8}
\]
Indeed, the second term in that theorem is the minimum of the two target
signed valuations when they have the same nonzero orientation, which is
precisely \(v_p(e'_{ij})\).
At an inert prime the same statement is the ordinary primitive-functional
inequality \(\beta_{ij}\ge\min(\rho_i,\rho_j)\), with
\(v_p(e'_{ij})=0\).

On the other hand, the determinant identity before and after ordinary row
normalization gives the exact equality
\[
 \beta'_{ij}+v_p(e'_{ij})
 =\ell+\beta_{ij}+v_p(e_{ij})-\rho_i-\rho_j.              \tag{9}
\]
Add \(\beta_{ij}\) to (9), and use (8). This proves the local pair lemma
\[
 \boxed{\displaystyle
 \beta_{ij}+\beta'_{ij}
 \ge \ell-|\rho_i-\rho_j|.}                              \tag{10}
\]
More precisely, the proof leaves the two nonnegative bonuses
\(v_p(e_{ij})+v_p(e'_{ij})\) on the right after (8) is used. Thus (10)
loses no untracked Gaussian common content.

The extreme zero-cost case in (10) requires one row content to be zero and
the other to be \(\ell\). Consequently a prime can escape one pair only by
placing its two rows at opposite ends of the full content interval. It
cannot do this for most pairs.

## 3. Sharp packing over all pairs

Sum (10) over the \(E=\binom m2\) unordered pairs:
\[
 \sum_{i<j}(\beta_{ij}+\beta'_{ij})
 \ge E\ell-\sum_{i<j}|\rho_i-\rho_j|.                    \tag{11}
\]
For arbitrary integers \(0\le\rho_i\le\ell\), layer cake gives
\[
 \sum_{i<j}|\rho_i-\rho_j|
 =\sum_{h=1}^{\ell}n_h(m-n_h)
 \le \left\lfloor\frac{m^2}{4}\right\rfloor\ell,          \tag{12}
\]
where \(n_h=\#\{i:\rho_i\ge h\}\). Since
\[
 \binom m2-\left\lfloor\frac{m^2}{4}\right\rfloor
 =\left\lfloor\frac{(m-1)^2}{4}\right\rfloor,
\]
we obtain the local packing theorem:
\[
 \boxed{\displaystyle
 \sum_{i<j}(\beta_{ij}+\beta'_{ij})
 \ge
 F_m\ell,\qquad
 F_m=\left\lfloor\frac{(m-1)^2}{4}\right\rfloor.}         \tag{13}
\]
The coefficient is the sharp consequence of the pair lemma for unrestricted
content values: equality in (12) occurs when every \(\rho_i\) is either
zero or \(\ell\), with the two groups as balanced as possible. No fairness
of source prime allocations or row contents has been assumed.

For \(m=3\), equality is attained by the explicit good-prime example from
the joint-content audit. Put \(Q=p^\ell\), with \(p=17\), and take
\[
 M=\begin{pmatrix}1&Q\\1&0\end{pmatrix},\qquad
 H_0=(1,0),\quad H_1=(Q,2),\quad H_2=(2Q,1).
\]
The source residue valuations are \((0,0,\ell)\), the row contents are
\((0,\ell,\ell)\), and the primitive target rows
\[
 (1,1),\qquad(3,1),\qquad(3,2)
\]
have target residue valuations \((0,0,0)\). Thus the left side of (13) is
exactly \(\ell=F_3\ell\). This example makes no endpoint-clustering claim,
but it shows that a larger local coefficient cannot follow from content
and pair transfer alone.

Let
\[
 c_{\rm odd}=\prod_{\substack{p\mid c\\p\ {\rm odd}}}
 p^{v_p(c)},\qquad
 B_*=\prod_{i<j}b_{ij},\qquad B'_*=\prod_{i<j}b'_{ij}.
\]
Multiplying (13) over all odd core primes yields
\[
 \boxed{\displaystyle c_{\rm odd}^{\,F_m}\mid B_*B'_*,
 \qquad\hbox{hence}\qquad
 c_{\rm odd}^{\,F_m}\le B_*B'_* .}                       \tag{14}
\]
This is a joint source-target statement. Neither residue product alone
need contain a positive fraction of every core prime.

## 4. The ramified prime

It remains to account for \(2\mid c\). Coprime sums-of-two-squares
integers \(a,b\) have only the following relevant parity types.

If one of \(a,b\) is even, the other is odd, and the parity factor is
\(\epsilon=4\). The even norm can contain a higher even power (for
example, a coprime Gaussian coefficient can have norm four), so no
congruence \(2\pmod4\) is asserted. What is needed here remains true:
\(b-a\) is odd and \(c=|b-a|\) is odd.

If both are odd, then \(a\equiv b\equiv1\pmod4\). The parity factor is
either \(\epsilon=1\) or \(\epsilon=2\). In the first case
\(K_{\rm core}=ab\) is odd. When \(2\mid c\), the ramified
good-determinant inequality already proved in the joint-content theorem
applies. It is the exact analogue of (8), including the target parity gcd,
so (10)--(13) hold at two without a loss.

The remaining case has \(\epsilon=2\),
\[
 K_{\rm core}=4ab,\qquad c=\frac{|b-a|}{2},\qquad
 \ell=v_2(c)\ge1.                                       \tag{14a}
\]
Here every \(k_j\) is odd: otherwise the reduced coprime Gaussian
coefficient pair would have both norms even and both members would be
divisible by \(1+i\). Thus (14a) gives the exact two-adic determinant and
discriminant valuations of every anchor-normalized matrix.
All source half-angle columns have odd norm. For the core matrix write
\(\rho_i=v_2(r_i)\), as before, and let
\[
 o_i=\begin{cases}
 1,&B_i\text{ has two odd coordinates},\\
 0,&\text{otherwise}.
 \end{cases}
\]
A primitive Gaussian target column contains \(1+i\) to exact exponent
\(o_i\). Hence
\[
 v_2(e'_{ij})=o_io_j,\qquad v_2(e_{ij})=0.               \tag{14b}
\]

Choose a row of the primitive core matrix containing an odd entry. It is a
primitive linear functional over \(\mathbb Z_2\). If the two raw image
contents have valuations \(\rho_i,\rho_j\), completing this functional to
an invertible coordinate system proves
\[
 \beta_{ij}\ge\min(\rho_i,\rho_j).
\]
The determinant identity (9), now with (14b), gives
\[
 \beta_{ij}+\beta'_{ij}
 \ge\ell-|\rho_i-\rho_j|-o_io_j.
\]
Summing and applying (12), with
\(r=\sum_i o_i\), yields
\[
 \sum_{i<j}(\beta_{ij}+\beta'_{ij})
 \ge F_m\ell-\binom r2.                                 \tag{14c}
\]
Combining the odd-prime result (14) with (14c), and using the stronger
loss-free inequality in the \(\epsilon=1\) case, gives the uniform
full-core estimate
\[
 \boxed{\displaystyle
 c^{F_m}\mid 2^{\binom r2}B_*B'_*,
 \qquad\hbox{hence}\qquad
 c^{F_m}\le 2^{\binom r2}B_*B'_* .}                     \tag{14d}
\]
Thus the entire two-part is controlled. The only extra charge is explicit
and depends on the parity cut of the full target rows.

## 5. Endpoint consequence without a target unit extraction

Suppose the primitive source has least squared radius \(N\), lies in an arc
of length at most \(C N^{1/4}\), and is in the odd fixed-unit normalization
used by the source radius theorem. Its pair product satisfies
\[
 B_*\le\left(\frac C2\right)^E N^{m/8}.                  \tag{15}
\]

Suppose independently that the full primitive target configuration has
least squared radius \(N'\) and lies in an arc of length at most
\(C'(N')^{1/4}\). No target unit-class extraction is required. If \(r\)
of its \(m\) primitive rows have two odd coordinates, the exact parity-cut
version of the target pair-product bound is
\[
 B'_*\le
 \left(\frac{C'}2\right)^E
 2^{r(m-r)/2}(N')^{m/8}.                                 \tag{16}
\]
This follows by multiplying the individual chord bounds and using
\[
 \prod_{i<j}d'_{ij}
 \le2^{r(m-r)}(N')^{m^2/4}
\]
for the primitive target pair norms \(d'_{ij}\). At two, a pair norm has
valuation one exactly across the odd-odd parity cut.

Combining (14d)--(16) proves
\[
 \boxed{\displaystyle
 c
 \le
 \left[
 \left(\frac{CC'}4\right)^E
 2^{r(m-1)/2}
 (NN')^{m/8}
 \right]^{1/F_m}.}                                      \tag{17}
\]
In logarithmic form,
\[
\log c\le
\frac{E}{F_m}\log\frac{CC'}4+
\frac{r(m-1)}{2F_m}\log2+
\frac{m}{8F_m}(\log N+\log N').                         \tag{18}
\]
For \(C,C'\le2\), the first term is nonpositive and may be dropped. Since
\(r\le m\),
\[
 c\le
 2^{m(m-1)/(2F_m)}(NN')^{m/(8F_m)}.                     \tag{19}
\]
If also \(N'\le N\), this becomes
\[
 \boxed{\displaystyle
 c\le
 2^{m(m-1)/(2F_m)}N^{m/(4F_m)}
 =2^{O(1)}N^{1/m+O(1/m^2)}.}                            \tag{20}
\]

Thus the complete determinant core, equivalently the reduced norm
difference \(\epsilon|b-a|/4\), has the desired two-ended
\(O((\log N+\log N')/m)\) logarithmic bound. This closes the proposed
local packing step, including high-content rows and full Gaussian pair
contents.

## 6. Exact remaining scope

The theorem controls \(a-b\), not \(a+b\). Near-singular maps can have a
large coefficient sum even when their full difference is small.
The two-adic calculation is included in (14a)--(14d), with its exact target
parity charge. Primes in the anchor scales \(k_j\), including
interval-distance primes from (4), are not charged by (14d) unless they
also occur in the fixed core.

Accordingly, (20) is a genuine new restriction on a radius-nonincreasing
endpoint map, but it does not by itself bound the full normalized losses
\(T_j^mK_j^{m-2}\), prove radius minimality, or improve the established
general point-count growth rate.

The persistent
[joint-anchor checker](check_mobius_joint_anchor_residues.py) verifies the
exact divisibility in (14d), the odd split and inert local inequalities,
the external-conductor statement from the all-anchor audit, and the target
parity-cut bound in 1,752 matrix/configuration cases. These include powers
of \(2,3,5,7,17\) through exponent three and an actual Pell source. The
finite checks supplement the proofs above.
