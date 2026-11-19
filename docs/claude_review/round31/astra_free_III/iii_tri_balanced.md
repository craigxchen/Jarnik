### A balanced four-cycle with the smallest primary-class cofactor

**PROVED using the exact COMPUTED certificate below.** There are one-unit, Gaussian-gcd-primitive four-point clusters on arbitrarily large norms \(N\) with **all seven nonconstant four-mark cuts of weight \((\log N)/7+O(1)\)**, normalized angular width \(\asymp N^{-2/7}\), and the same primitive cofactor **two** on every edge of a \(K_{2,2}\). In the primary prime convention, two is the smallest positive cofactor on a one-unit edge: its primitive factor is a real unit times a product of primary factors, whose imaginary coordinate is even. This gives actual integer common-parent examples at the four-mark exponent obtained by restricting an asymptotically balanced growth-scale profile.

*Proof.* Use the parameters \(x=-5,m=3\) of Part I §6 and clear the root denominators by six. Here is the complete resulting construction. Index the following Gaussian roots starting at zero:
\[
 (0,6),\quad(-135,-36),\quad(-135,-51),\quad(-135,75),\quad
 (-145,81),\quad(-117,-45),\quad(-108,-30),
\]
where \((u,v)=u+iv\). Put \(H_j=n-\rho_j\). The four rows are
\[
 S_0=\varnothing,\quad S_1=\{0,1,2,4\},\quad
 S_2=\{0,1,3,5\},\quad S_3=\{0,2,3,6\}.
\]
Take \(z_i=\prod_{j\in S_i}H_j\prod_{j\notin S_i}\bar H_j\) and divide the four by their common Gaussian gcd. For \(i<j\) their ratio is \(B_{ij}/\bar B_{ij}\), where \(B_{ij}\) uses \(H_k\) on \(S_i\setminus S_j\) and \(\bar H_k\) on \(S_j\setminus S_i\). The two sides have four factors together. Each resulting monic quartic has constant nonzero imaginary part \(b_{ij}\). In lexicographic edge order \(01,02,03,12,13,23\), multiplying the four linear factors gives
\[
 (|b_{ij}|)=(19966500,16110900,14871600,3855600,5094900,1239300).
\]
In particular the angular width is \(\asymp n^{-4}\). On the arithmetic progression
\[
 n\equiv79658575767\pmod{86434978500}
\]
the four gcds \(\gcd(\Re B_{ij}(n),b_{ij})\) for \(ij=01,02,13,23\) are respectively \(9983250,8055450,2547450,619650\). Verification at one residue suffices since the period is the lcm of the four corresponding \(|b_{ij}|\). After dividing out the gcd, each primitive near-real representation has imaginary coordinate of absolute value two and odd real coordinate. Thus the four relevant norms are odd and every cross edge has cofactor two. Ratios of these near-real factors to their conjugates have primary unit one, so connectivity of the four-cycle puts all four points in one primary unit class.

For the cut weights, the seven root signatures run once through all seven nonempty subsets of marks \(1,2,3\), with mark zero as anchor. The seven norms \(\mathrm NH_k\) are \(\asymp n^2\). All fourteen roots including conjugates are distinct. Any common Gaussian factor of two of the corresponding fourteen linear values divides their fixed nonzero root difference. Hence, after discarding a bounded Gaussian factor from each linear value, the fourteen residual ideals are pairwise coprime, conjugation permutes the residual ideals, and every residual prime contributes only the prescribed cut. For example, take a fixed ordinary integer divisible by all the pairwise root differences in \(\mathbb Z[i]\), and remove from each prime valuation the smaller of the valuation itself and that fixed bound. Only bounded factors have been removed. Restoring them affects at most a bounded number of weighted layers in total and changes the common gcd by a bounded factor. This proves \(N\asymp n^{14}\) and gives each cut weight \(2\log n+O(1)\), as claimed.

**COMPUTED.** The stand-alone script multiplies the quartics over \(\mathbb Z[i]\), checks the constant offsets and the progression, then produces the four Gaussian integers at the displayed residue and computes their primitive half-angle tangents directly. Their six absolute offsets are \((2,2,7344,560,2,2)\) and the norm has 135 digits. Rational inequalities certify critical arc constant smaller than \(1/128\). The seven prime blocks can be composite, as in a growth profile; no claim of squarefree values is needed. [Balanced four-cycle certificate](iii_tri_balanced.py)

**PROVED — the numerical calibration.** In a large squarefree-profile model with sixteen unoriented cut weights for five marks (or with fair bits on four marks), a balanced four-subtuple has primitive log norm \((7/8+o(1))W_{\rm parent}\) after deletion of its constant cut. A parent arc of width \(\ll\exp(-W_{\rm parent}/4)\) has the four-mark exponent \(2/7+o(1)\). Our example realizes this scale with exact cofactors on a cycle. Thus a forbidden-*unordered*-four-cycle argument cannot be justified by those four-point consequences of multiplicity, even after imposing all seven cut weights and the integer origin of the star matrix. A restriction on how many such cycles can share vertices of a much larger parent would go beyond this example.
