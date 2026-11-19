### A rational configuration with unrelated minimum-support divisors, including specialization

**PROVED using the COMPUTED certificates below.** In degree eight there is a rational root configuration with seven signed ternary directions of weight four up to reversal and with maximum cluster size four. Five of its seven directions simultaneously give divisors counted by \(D_B(N)\) along one infinite sequence of odd split-supported integer parents, for *each* fixed \(B>0\). No four integer lattice points in a short arc account for all five selected divisor labels on that sequence once the specialization parameter is large. This is a fixed lower bound and gives no violation of the growth conjecture for \(D_B\).

**COMPUTED — configuration.** Number the roots from 0 through 7. Their pairs \((x,y)\), denoting \(x+iy\), are
\[
 (-12,4),(-8,6),(-3,16),(0,10),(0,12),(3,16),(8,6),(12,4).
\]
Put \(u_j=(y_j,x_jy_j,3x_j^2y_j-y_j^3)\). A weight-four ternary solution must have two positive and two negative entries by the mass bound. Thus it comes exactly from an equality \(u_i+u_j=u_k+u_l\) with distinct indices. All pair-sum coincidences, abbreviated by their indices, are
\[
 04=13,\qquad 05=27,\qquad 06=17,\qquad 15=26=34,\qquad 36=47.
\]
Count all three pairwise equalities on the length-three chain. For a cluster of sign rows, the corresponding subsets must have one and the same sum of the \(u_j\), and this is sufficient by Newton's identities. Enumerating the 256 subsets gives fibre-size frequencies \(120,32,8,12\) for sizes one through four. The elementary checks use only integer arithmetic. [Certificate](ii_template.py)

**PROVED — the specialization cost.** Choose the five equalities \(04=13,06=17,15=34,26=34,36=47\). For any one of them take \(B_v(n)\) to be the product of \(n-\rho_j\) for the two indices on the left and \(n-\bar\rho_j\) for the two on the right. The equality of the three moments makes its first three non-leading coefficients real, so \(b_v=|\Im B_v(n)|\) is constant. In the listed order the five constants are \(4800,9600,13200,13200,4800\), evaluated at \(n=0\). Put
\[
 f_j(n)=(n-x_j)^2+y_j^2,\quad N_0(n)=\prod_j f_j(n),\quad
 g_v=\gcd(\Re B_v(n),b_v),\quad d_v=\mathrm N B_v(n)/g_v^2,\quad t_v=b_v/g_v.
\]
For the natural monic parent \(N_0\), the ratio \(d_v/(16t_v^2\sqrt{N_0})\) tends to \(1/(16b_v^2)\), so passage to integer polynomials alone loses the desired cutoff. Take instead \(N=\operatorname{lcm}_v d_v\) and \(G=N_0/N\). Then the ratio is multiplied by \(\sqrt G\).

**COMPUTED — a simultaneous CRT certificate.** For each prime in the first row below impose \(n\equiv n_p\pmod{p^{k_p}}\) with \(k_p=1+v_p(N_0(n_p))\). The last row is the valuation of \(G(n_p)\).
\[
\begin{array}{c|rrrrrrrr}
p&2&3&5&13&37&41&53&509\\
n_p&2&0&10&5&12&17&0&238\\
k_p&16&3&12&5&3&3&3&3\\
v_p(G)&15&2&9&1&1&1&1&1
\end{array}
\]
All five \(d_v(n_2)\) are odd. The CRT representative printed by the script has 43 digits and gives a 666-digit odd parent with five distinct candidates. At that parent it verifies the upper cutoff even at \(B=5\) and obtains \(t_v=(2,12,11,11,2)\). [Certificate](ii_template.py)

**PROVED — passage to the progression.** Each displayed valuation persists on the progression: \(p^{k_p}\) fixes the valuations of \(N_0,\mathrm N B_v\) since the latter divides the former, and of the gcd \(g_v\) since \(g_v^2\mid\mathrm N B_v\). Thus \(G\) is divisible by \(L=306439728192000000000>256\max_v b_v^4\). Each \(d_v\) is an odd primitive sum of squares, hence \(N\) is odd and split-supported. Uniformly on the progression,
\[
 \frac{d_v}{16t_v^2\sqrt N}
 =\frac{\sqrt G}{16b_v^2}\frac{\prod_{j\in\mathrm{supp}(v)} f_j}{\sqrt{N_0}}>1
\]
eventually because the last quotient tends to one. The \(G\)'s are bounded above independently of \(n\). Indeed each of the eight \(f_j\)'s occurs on some selected support, and each such \(f_j(n)\) divides \((\operatorname{lcm}_v b_v)^2 N(n)\). The quadratics \(f_j\) are pairwise coprime over \(\mathbb Q[t]\), so pairwise gcds of their integer values divide fixed nonzero resultants. In particular \(\prod f_j/\operatorname{lcm}_j f_j\) is bounded (it divides the product of the pairwise gcds). This gives \(N\asymp n^{16}\), \(d_v/\sqrt N=O(1)\) and the upper cutoff for every fixed \(B>0\) eventually. Distinct supports yield five distinct \(d_v\)'s for large \(n\), since there are finitely many possibilities for \(g_v\mid b_v\) and the different products of the monic, distinct irreducible \(f_j\)'s cannot be proportional as polynomials. The floor condition follows either from Part I or because \(t_v\) stays bounded while the reduced real part tends to infinity.

**PROVED using a finite check — meaning of unrelated.** No polynomial subset-sum fibre supplies more than three of the selected labels. More intrinsically, of their ten distinct triples exactly one has a signed sum of direction vectors zero. In zero-based order of *all seven* equalities (break the chain into \(15=26,15=34,26=34\)), it is \(v_0-v_2+v_6=0\). Gaussian unique factorization of the monic linear factors shows that any of the other three-direction products \(\prod (B_v/\bar B_v)^{\pm1}\) is a nonconstant rational function. Hence these other triangles still fail for all sufficiently large integer \(n\). Five distinct edges among four points would contain two triangles. On a short arc each such triangle would impose a signed phase closure, since \(B_v/\bar B_v\) is fixed up to inversion from its near-square divisor. This proves the claim about the five selected labels. The construction does not bound the size of integer clusters involving further divisors.

**EVIDENCE — growth.** This template supplies a fixed finite number while \((\log N/\log\log N)^2\) tends to infinity. Even an increasing sequence of polynomial examples would need a lower bound on the number of distinct eligible directions at one specialization in terms of that specialization's actual \(\log N\), including denominators, CRT modulus and angle cutoff. The existence statement of the bridge theorem by itself provides none of those uniform estimates.

**COMPUTED with the same PROVED progression argument — larger count.** The eight roots \((-23,4),(-11,10),(-8,16),(2,16),(5,10),(5,12),(13,6),(17,4)\) give nine ternary directions and fibre-size frequencies \(104,32,16,8\) for sizes \(1,2,3,5\) respectively. Seven of their directions simultaneously specialize to \(D_B(N)\) divisors on an infinite progression. No polynomial fibre accounts for all seven. The check at a 752-digit parent gives reduced offsets \(123,24,3,84,12,11,2\). In the script's order choose directions \(0,1,2,3,5,7,8\). The local tuples \((p,n_p,k_p,v_p(G))\) are \((2,3,16,15),(5,314,13,10),(41,12,3,1),(53,16,3,1),(349,131,3,1),(541,242,3,1),(953,544,3,1)\), giving \(L=125119585606720000000000>256(98400)^4\). There are also two midpoint relations in each configuration, taking the full weight-four counts up to sign to nine and eleven respectively. They are verified by testing all distinguished midpoints and all pairs of distinct endpoints. This improves the exhibited finite divisor count. The previous example gives the stronger four-vertex incompatibility. [Second certificate](ii_template.py)

### Two obstructions to obtaining growing examples by symmetry or composition

**PROVED.** For any even \(D=2m\), suppose the \(D\) distinct upper roots lie on a common ellipse with real center and axes parallel to the real and imaginary axes. There are at most two ternary directions of weight \(m\), up to reversal. Indeed, parametrize points as \(x_0+Aw+B/w\), \(|w|=1\), with real \(A=(a+b)/2\), \(B=(a-b)/2\) for positive semiaxes \(a,b\). The imaginary parts of the powers through degree \(m-1\) are triangular real combinations of \(\Im w^j\) with nonzero diagonal \(A^j-B^j\). Thus the moment equations for any signed root selection are equivalent to the unit-circle equations. For a selection of \(m\) roots upstairs on that unit circle, let \(F(z)=z^m+\sum_{j=0}^{m-1}a_jz^j\) be the corresponding polynomial. Its interior coefficients are real by Newton's identities, and \(a_0\) is nonreal because a real polynomial would contain every selected conjugate. The unit-root identity \(a_j=a_0\bar a_{m-j}\) forces every interior coefficient to vanish. Consequently every selection is a whole fibre \(w^m=\eta\). Each pair of conjugate selections uses \(2m\) roots, distinct such pairs are disjoint, and the full upper and lower configuration has \(4m\) roots.

**PROVED.** Given a fixed rational Gaussian template with distinct nonreal roots \(\rho_j,\bar\rho_j\), there is one \(a\in\mathbb Z\) such that for *every* positive integer \(k\), substitution \(t\mapsto t^k+a\) introduces no new half-degree constant-offset divisors over \(\mathbb Q(i)[t]\). Choose distinct odd split rational primes \(p_j\) avoiding the denominators and a Gaussian prime \(\pi_j\mid p_j\) for each upper root. The map \(\mathbb Z/p_j^2\mathbb Z\to\mathbb Z[i]/\pi_j^2\) is an isomorphism by equal cardinalities and contraction; rational denominators prime to \(p_j\) are units there. Choose a residue \(a\pmod {p_j^2}\) for which \(v_{\pi_j}(a-\rho_j)=1\), and combine the conditions by CRT. At the conjugate root the analogous valuation holds at \(\bar\pi_j\). Every polynomial \(t^k+a-\rho_j\) and its conjugate is now Eisenstein in the relevant discrete valuation ring. The substituted polynomials are distinct monic irreducibles, so every monic divisor is exactly a union of whole fibres and equals \(F(t^k+a)\) for a divisor \(F\) downstairs. The degree is multiplied by \(k\), and \(\Im F(t^k+a)\) is constant precisely when \(\Im F\) is constant. This gives an obstruction for rational-polynomial divisors under these substitutions. The roots upstairs need not be rational Gaussian and the proof leaves split substitutions with partial fibres open.
