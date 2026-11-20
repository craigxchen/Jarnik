## A fixed-parent model with all powers of all bounded unbalanced words

**PROVED.** For every fixed integer \(H\ge2\) and critical arc constant \(C>0\), there is a sequence of models with \(m\asymp_H s\to\infty\), \(W=O_{H,C}(m\log m)\), on each of which the ordinary sharp monomial gaps, P3-char at every cutoff, and the strengthened e-all inequalities at **all** odd exterior prime powers and every unit hold for every \(n v\), \(n\ge1\), whenever
\[
 v\in\mathbb Z^s,\qquad 0<\|v\|_1\le H,\qquad \textstyle\sum_j v_j\ne0 .
\]
This includes every coordinate ray and every two-label sum ray with the respective quantifiers of IX.7 on the same fixed parent. These local data cannot have all the rational lifts in IX.7, even allowing finite sets of exceptional exterior primes. This uses no random-parent statement.

*Proof.* Start with the replicated all-height P3-char construction of IX.5 with \(g\ge 256 H^2\) copies per coarse Paley column and coarse Paley prime number \(q\). Use distinct physical primes in a bin \([P_0,P_0e^\sigma]\subset[P,2P]\) with \(P=q^d\), fixed \(d>2\), and the small fixed angular sector, as in IX.3 and the elementary completion in IX.5. The coarse angle of each block and all windings may be fixed before the fine angular split. The physical scope here requires the unconditioned Haar split only. The norm-one ratios associated to the IX.3 witnesses are units at exterior prime powers. At each such power replace their individual values within a block, preserving the product and the parities prescribed by the physical primes. For two mark rows the exponent difference is constant on each block. Thus the new point character residues agree with the old ones up to the common point factor that cancels in each character. External e-char and all same-level pair clauses at parent primes continue to hold together: the witnesses at parent primes are left alone.

Here is the local replacement. For an exterior odd rational prime \(l\) and its power \(l^a\), the norm-one group is cyclic of order \(b=(l-\chi_4(l))l^{a-1}\), with the generator chosen compatibly down the tower. Reduction to lower levels gives the corresponding reductions of generator exponents. A ratio with discrete log \(L_j\) is obtained from a unit of norm \(p_j\) if \(L_j\) has the same parity as the initial ratio. Indeed starting from the latter factor one multiplies it by a norm-one element, squaring that element on the ratio. These basic local facts are also proved in round 6, sharp.md, Proposition 2.1. Denote by \(B\) the prescribed discrete log of the block product. Its parity equals the sum of the physical parities.

For each prime divisor \(r\) of \(b\) put \(k_r=\min\{k\ge1:r^k\ge 256H^2\}\). If \(r^k\Vert b\) with \(k<k_r\), impose only the prescribed block sum modulo \(r^k\) and the individual parity at 2; this is always possible. If \(k\ge k_r\), put \(D=r^{k_r}\) and assign each \(L_j\bmod D\) a representative in the integer interval
\[
 \Big[\frac D{2H}-\frac D{16H^2},\ \frac D{2H}+\frac D{16H^2}\Big].
\]
Choose the parity prescribed for this copy when \(r=2\). The block sums modulo \(D\) attain every possible residue, of the correct parity when \(r=2\). Indeed the allowed integers in one copy form an interval with step one or two. Even after rounding with the parity restriction their span is at least \(D/(8H^2)-4\ge D/(16H^2)\). The Minkowski sum is an interval of the same step and has span at least \(gD/(16H^2)\ge D\). Lift the entries to obtain the block sum \(B\) modulo \(r^k\), adjusting a multiple of \(D\) on one slot. CRT over the prime-power factors of \(b\) gives the local discrete logs. At every lower level the same nonvanishing and divisibility bounds hold by reduction. Finite branching and compactness (or König's lemma at fixed \(l\)) now give compatible norms and replacements for all powers at once. At finite stages choose any norm-\(p_j\) representatives of these ratios. Norm and block-ratio conditions descend and are closed, so the same compactness gives compatible representatives.

Write \(u_v=\sum v_j L_j\). At each \(D\) imposed above, its integer value using the chosen representatives differs from \(D\sum v_j/(2H)\) by at most \(D/(16H)\). Hence when \(0<|\sum v_j|\le H\) its distance from \(D\mathbb Z\) is positive. Thus \(\gcd(u_v,b)\) divides the fixed constant
\[
 K_H=\prod_{r<256H^2}r^{k_r-1},
\]
including when \(k<k_r\). If \(n v\) hits any fourth root of unity modulo \(l^a\), then \(b\mid4K_H n\). In particular \(l\le 4K_Hn+1\) and \(l^{a-1}\mid4K_H n\). The product of all maximal colliding powers, even if each prime is charged at the unit of its own choice, is at most
\[
 4K_H n\,\exp(\vartheta(4K_H n+1))\le \exp(c_H n)
\]
for a fixed \(c_H\), by Chebyshev's bound. This holds at every cutoff simultaneously.

Complete the fine angles by Haar independently on each of the specified sum fibres, matching actual real lifts and their windings. Every \(v\) under test is nonconstant on some block since \(\|v\|_1\le H<g\). Its dot product, and that of \(n v\), have uniform image on the circle conditional on the other blocks. Working modulo \(2\pi\) on the fine fibres, with each real block sum equal to the desired angle, handles all four units; the measure bound under any equivalent quarter-turn parametrization with fixed real lifts is also absolute by exposing a single coordinate with nonzero coefficient difference. For mass \(h\) there are at most \((2gq)^h\) words. The total new failure probability is
\[
 \ll\sum_{h=1}^H(2gq)^h\sum_{n\ge1} \exp(c_H n) P^{-hn/2}=o(1).
\]
The ordinary fine mixed-word gaps cost \(\ll\sum_{h\ge1}(2gq)^h P^{-h/2}=o(1)\). The coarse words were already protected in IX.5. The exact mark offsets need no change. Almost surely a fine angle differs from either genuine Gaussian factor's angle modulo \(\pi\); selecting among the successful positive-measure set therefore gives nongenuine prime angles. Any cofinite rational lift at every physical label would contradict IX.7 applied to the singleton rays, so in this construction at least one such lift fails (and in fact all can be made to fail). This proves the claimed model.

**PROVED (scope).** The same argument works for each nonzero word whose specified coefficient sum on a fixed finite partition is nonzero, by replacing the common center \(1/(2H)\) by sufficiently small fixed rational multiples encoding the finitely many colors; the width and replication factor then depend on that partition and \(H\). In particular, all positive one-sign words of bounded mass are already covered by the single-color version. Arbitrary balanced signed words are essential omissions. The smallest example is a difference between two distinct physical labels. Uniformly high local order for every such difference cannot hold on growing sets in general, due to pigeonhole mod divisors smaller than their cardinality. The preceding proof does not give a full e-all fake.

## A height barrier for all right inverses and all finite lifts of the owner design

**PROVED.** On one **unreplicated** Paley profile of IX.3, let \(q\) be its order, \(r=q-1\) the number of labels and \(k=m-1\) the number of independent anchored rows. Put \(\mathcal B_X=\prod_{l^a\le X,\,l\nmid2N}l\), where each admissible power contributes one copy of \(l\). Use physical prime labels \(p_j\in[P,2P]\). Suppose norm-correct rational lifts \(r_j\) of prime factors of those labels reproduce the IX.3 exterior *anchored mark ratios*, including prescribed windings, at all powers up to \(X\). For primitive one-sided directions write \(r_j/\bar r_j=u_j\alpha_j/\bar\alpha_j\), with \(u_j\) a unit, and \(J_j=|\alpha_j|^2\). Then
\[
 \sum_{j=1}^r\log J_j\ \ge \min\left(\log\mathcal B_X-\log(4P),\ (r-k)q^{k/(2(r-k))}\log P\right).
\]
In particular at \(X\ge m^{3/2+\epsilon}\) and along the high-density diagonal \(m/q\to1/2\), even the mean logarithmic height is at least \(q^{1/2-o(1)}\log P\). The bound for exact witnesses (including every right inverse) is the second term alone. All estimates apply to fixed primes and the specified mark-ratio design. They give no nonexistence result for other local designs.

*Proof.* With the anchor mark indexed by \(a\), the target ratio for \(i\ne a\) is \(\gamma_i/\gamma_a\), \(\gamma_i=\pi_i/\bar\pi_i\). If any of these exterior congruences fails as a rational identity, clear both true owner primes and the primitive denominators of \(\prod_j(r_j/\bar r_j)^{D_{ij}}\), including the phase units. Since all reduced denominators are invertible at tested primes, this yields a **nonzero** Gaussian integer divisible by the rational integer \(\mathcal B_X\). Its modulus is at most \(2|\pi_i\pi_a|\prod_j J_j^{|D_{ij}|/2}\le4P\prod_j J_j\), because \(|D_{ij}|\le2\). This proves the first alternative.

Otherwise the rational identities are exact. Form the integer physical valuation matrix
\[
 F_{jt}=v_{\pi_t}(r_j/\bar r_j).
\]
The norm equation gives \(F\equiv I\pmod2\), so \(\det F\) is odd, even if the lifts involve other primes. Each true anchored ratio then gives \(DF=H\), where \(H_i=e_i-e_a\). The truncated Paley difference Gram is \(DD^T=q(I+\mathbf1\mathbf1^T)\) on the \(k\) rows by the Paley frame (the deleted column was constant). Meanwhile \(HH^T=I+\mathbf1\mathbf1^T\). Consequently right multiplication by \(F\) decreases the \(k\)-volume on the span of the rows of \(D\) by exactly \(q^{-k/2}\). If \(t_1\ge\cdots\ge t_r>0\) are its singular values, the product of the last \(k\) is at most \(q^{-k/2}\) by the exterior-power minimum singular value. Therefore \(\prod_{j\le r-k} t_j\ge q^{k/2}\), as \(\prod_j t_j=|\det F|\ge1\). The nuclear norm is thus at least \((r-k)q^{k/(2(r-k))}\) by AM–GM. By the triangle inequality for that norm and the matrix units, \(\sum_{j,t}|F_{jt}|\ge\sum_jt_j\). Finally \(\log J_j\ge(\log P)\sum_t|F_{jt}|\), since \(\alpha_j\) contains one side of each physical factor with exactly the signed exponent's absolute value; all omitted factors have height at least one. The bound follows. For the asymptotic application \(\log\mathcal B_X\ge\psi(X)-(r+1)\log X\gg X\) for \(X\ge m^{3/2+\epsilon}\). \(\square\)
