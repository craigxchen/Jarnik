### The sharp max-scale cross-gap against Gaussian primes

**PROVED.** There exist an absolute \(p_0\) and a probability measure \(\mu\) on quarter-turns with power Fourier decay, such that for each \(y\) in its support, each rational prime label \(p\ge p_0\), and every Gaussian prime \(\varpi\),
\[
 |\sin((\pi/2)y-\arg\varpi)|\ \ge\ 1/\max(p,\operatorname{Norm}\varpi).
\]
The same support satisfies the lattice max-scale inequality against every nonzero Gaussian integer with an absolute positive constant. These assertions include all associates and conjugates.

*Proof.* We first arrange that certain primitive integer vectors have composite Gaussian norm. For each nonzero \((a,b)\in[-5,5]^2\cap\mathbb Z^2\) choose a distinct split prime \(\ell_{a,b}>50\). Write \(u=a,v=7a+b\); then \(0<|v|\le40\). Choose an integer \(L\) by CRT so that \(u+Lv\equiv\sqrt{-1}\,v\pmod{\ell_{a,b}}\) for all the selected pairs. The unimodular integer matrix \(B:(a,b)\mapsto (u+Lv,v)\) thus satisfies \(\ell_{a,b}\mid\|B(a,b)\|^2\). Put \(Q=\prod\ell_{a,b}\).

Choose a fixed multiple \(D\) of \(Q\) large enough that \(D\ge102+\|Be_1\|/\|Be_2\|\), and let \(a_i=iD\) for \(i=1,2,3\). In the upper half plane let \(g_i\) be the reflection in the semicircle of radius one with center \(a_i\). On the real boundary \(g_i(x)=a_i+1/(x-a_i)\). The orientation-preserving subgroup \(\Gamma\) of the three-reflection group is a classical rank-two Schottky group. Indeed, the generators \(g_1g_2,g_1g_3\) pair the two disjoint half-disks \(D_2,D_3\) with their respective disjoint images \(g_1D_2,g_1D_3\) inside \(D_1\), mapping exterior to interior. Pairing these disjoint geodesic boundaries is the classical Schottky construction of a convex-cocompact group. Its limit points belong to the compact set of infinite reduced words \(i_1,i_2,\ldots\), with successive letters different. Indeed finite reduced words nest successive images of the three disjoint disks, and on the boundary interval of a different letter each map contracts by at least the factor \((D-1)^{-2}\) in the sense that its Lipschitz constant is at most that number. The resulting boundary points are the limits of \(g_{i_1}\cdots g_{i_{n-1}}(a_{i_n})\). One can use the same coding for the even subgroup since it has finite index.

Map a real slope \(x\) to the planar direction of \(B(1,x)\). For the limit of any such word let \(\phi\) be this direction in radians. Its convergents are represented by unimodular-neighbor vectors
\[
 P_0=Be_2,\quad P_1=B(e_1+a_{i_1}e_2),\quad
 P_{n+1}=(a_{i_{n+1}}-a_{i_n})P_n+P_{n-1}.
\]
This follows by expanding the boundary compositions as signed continued fractions with successive entries \(a_{i_1},a_{i_2}-a_{i_1},\ldots\). Put \(h_n=\|P_n\|\). By the choice of \(D\) and induction \(h_{n+1}\ge100h_n\) for all \(n\), and \(h_{n+1}\le(2D+1)h_n\) for \(n\ge1\). The adjacent determinants have absolute value one. Summing \(|\sin(\arg P_j-\arg P_{j+1})|=1/(h_j h_{j+1})\) gives
\[
 |\sin(\phi-\arg P_n)|\le 2/(h_n h_{n+1}).
\]
Here \(|\sin(s+t)|\le|\sin s|+|\sin t|\), which also handles sign reversals of the vectors.

For a primitive integer vector \(w\), put \(t=\|w\|\). If \(t^2|\sin(\arg w-\phi)|<1\), then
\[
 |\det(w,P_n)|\le h_n/t+2t/h_{n+1}.
\]
Take \(h_n\le t\le h_{n+1}\) with \(n\ge1\). The finitely many lengths below \(h_1\) will be handled below. For \(4h_n\le t\le h_{n+1}/4\) the determinant is zero, impossible since the parallel primitive vector would have length \(h_n\). For \(t\le4h_n\) write \(w=A P_n+B'P_{n-1}\). Its determinant gives \(|B'|\le1\), and the length bound gives \(|A|\le5\). For \(t\ge h_{n+1}/4\) instead use the basis \(P_{n+1},P_n\). The determinant bound at \(n+1\) gives \(|B'|\le4\), and \(|A|\le2\). In either case the basis modulo \(Q\) consists of \(Be_1,Be_2\) in alternating order by the recurrence. Therefore some chosen \(\ell_{a,b}\mid\operatorname{Norm}(w)\). A primitive Gaussian prime of norm exceeding the largest chosen \(\ell_{a,b}\) cannot have such a direction. Imprimitive Gaussian primes are inert rational primes times units and have axial directions.

Every limit direction is irrational, since a fixed rational direction would give a zero integer determinant with all sufficiently large \(P_n\) by the sine estimate. The directions form a fixed compact set. For the finitely many primitive Gaussian prime vectors of length at most \(\max(\max_i\|B(e_1+a_i e_2)\|,\sqrt{\max\ell_{a,b}})\), and for the axial directions, the sine distance to this set has a positive minimum (including associates). Choose \(p_0\) larger than its reciprocal. This proves the claimed sharp prime inequality at all norms. For the fixed-constant lattice version, given a sufficiently long primitive \(w\), choose \(n\) minimal with \(h_n\ge4t\). Its determinant with \(P_n\) is a nonzero integer, so
\[
 1\le h_n t\, |\sin(\arg w-\phi)|+2t/h_{n+1},\qquad h_n\le4(2D+1)t.
\]
The remaining finite directions use compactness, and imprimitive multiples only increase the norm. Thus a fixed \(c>0\) works for the lattice version.

Take the Patterson–Sullivan probability measure of \(\Gamma\) on its limit set and push it forward by \(x\mapsto (2/\pi)\arg B(1,x)\) modulo one. The phase has a smooth real lift on a compact interval containing the limit set, with derivative \((2/\pi)/\|B(1,x)\|^2\) bounded away from zero there. Bourgain–Dyatlov's oscillatory-integral theorem for convex-cocompact Fuchsian groups therefore gives \(|\widehat\mu(k)|\le C(1+|k|)^{-\varepsilon}\) for fixed \(C,\varepsilon>0\). The theorem applies with smooth amplitude one on the limit set and the displayed phase, allowing the fixed smooth density between Euclidean and spherical normalizations. [Bourgain and Dyatlov, Theorem 2](https://arxiv.org/abs/1704.02909) \(\square\)

**PROVED (conditional measures on every angular sum).** Let \(\mu\) be any probability on \(\mathbb T\) with such a power Fourier bound. There are a fixed integer \(g\), constants \(A<\infty,\beta>0\), and, for every \(\Theta\in\mathbb T\), a probability supported on \((\operatorname{supp}\mu)^g\cap\{\sum x_i=\Theta\}\), such that for every nonconstant integer coefficient list \(v_1,\ldots,v_g\), with \(h=\sum|v_i|\), and every interval \(J\),
\[
      \Pr\left(\sum v_i x_i\in J\right)\ \le\ A(1+h)^2 |J|^\beta .
\]
The constants are uniform in \(\Theta,v,J\).

*Proof.* Choose \(0<\beta<\min(1,\varepsilon)\). The Fejér kernel at order \(n\asymp1/|J|\), chosen small enough that its central peak bounds \(c n\) on \(J\) translated to the origin, gives \(\mu(J)\ll n^{-1}\sum_{|k|<n}|\widehat\mu(k)|\ll |J|^\beta\). Thus \(\mu\) is atomless and \(|\widehat\mu(k)|<1\) for every \(k\ne0\), since the equality case would give support on a finite coset. The numbers \(\sum_{k\ne0}|\widehat\mu(k)|^r\) tend to zero by dominated convergence. Choose \(g_0\) so every \(\mu^{*r}\) for \(r\ge g_0\) has continuous density in \([1/2,3/2]\). Choose fixed \(g>g_0+3/\varepsilon+2\).

Write \(\Sigma=\sum x_i\), \(V=\sum v_i x_i\) for independent samples under \(\mu^g\). If some coefficient \(a\) occurs at least \(g_0\) times, then conditional on the other samples the density of \(\Sigma\) is at most \(3/2\). If simultaneously \(\Sigma\in I\) and \(V\in J\), the sum of the other variables with coefficients \(v_i-a\) belongs to \(J-a I\). This arc has length at most \(\min(1,|J|+h|I|)\). At least one outside coefficient differs from \(a\), and testing it with the Frostman bound and the integer multiplication map gives probability \(O((1+h)^{1-\beta}(|J|+h|I|)^\beta)\). Hence the joint probability has this bound times \(3|I|/2\).

If no value occurs \(g_0\) times, the joint Fourier coefficients of \((\Sigma,V)\) are \(\prod_i\widehat\mu(k+v_i l)\). Let \(M=\max(|k|,|l|)>0\). If \(|l|<M/(4(h+1))\), then \(|k|=M\) and every argument has absolute value at least \(M/2\). Otherwise all but at most one distinct coefficient value have \(|k+v_i l|\ge|l|/3\ge M/(12(h+1))\). In both cases at least \(g-g_0\) factors have that lower bound on their arguments. Summing the resulting bound \(\min(1,O(((1+h)/M)^{(g-g_0)\varepsilon}))\) over integer frequency pairs yields \(O((1+h)^2)\), since \((g-g_0)\varepsilon>2\). Thus the joint law has bounded density \(O((1+h)^2)\), giving joint probability \(O((1+h)^2 |I||J|)\).

For fixed \(\Theta\) condition \(\mu^g\) on \(\Sigma\in I\) for shrinking intervals centered at \(\Theta\); their mass is at least \(|I|/2\). Compactness gives a weakly convergent subsequence with limit supported on the required fibre. For an open \(J\) the preimage under \(V\) is open, so the Portmanteau inequality carries both divided joint bounds to the stated estimate. Enlarge general intervals by an arbitrarily small amount. The same weak limit works for every integer word and interval, and no assertion for merely almost every \(\Theta\) is used. \(\square\)

**PROVED (application).** Replicate each Paley column \(g\) times with disjoint primes of size \(P=q^d\), where \(d>\max(2,2/\beta)\) is fixed. In each block fix the coarse angle as before and choose a sample having that prescribed angular sum. For a relevant nonconstant block the displayed estimate gives a risk \(O((1+h)^2P^{-\beta h/2})\) conditional on every other block. The number of signed words is at most \((2gq)^h\); thus every mixed gap passes with probability \(1-o(1)\). All block-constant monomials are handled by coarse weights \(\sum_a\log p_{ja}\), and the original realizations and residues are preserved. Real lifts of each block sum can be matched exactly by adding integers to any one coordinate. Take the auxiliary owners \(\prod_a\pi_{ia}\) in the all-height factored construction; their arguments lie in a short sector by taking the Gaussian-prime bin of width \(O(C/g)\). Parent log weight remains \(O(m\log m)\). Each selected angle now has the sharp constant one against every genuine Gaussian prime at all norms, and a fixed positive constant against every Gaussian integer.

**PROVED (the finite, sharp test against every lattice direction, on this same model).** Along an unbounded sequence of the above parents one can also require, for each label prime \(p\),
\[
 |\sin(\phi_p-\arg w)|\ge1/p\qquad(0\ne w\in\mathbb Z[i],\ \operatorname{Norm}w\le p).
\]
This is the finite test used in the Part VIII referee illustration. The prime-only, all-height statement and the fixed-constant lattice statement continue to hold.

*Proof.* For a limit angle in the Schottky construction, the first-step determinant and the tail sine estimate give
\[
 |\sin(\arg P_n-\phi)|\ge {1\over h_n h_{n+1}}-{2\over h_{n+1} h_{n+2}}
 \ge {1-2/100^2\over h_n h_{n+1}}.
\]
Every real number \(P\) in \([2h_n h_{n+1},\,h_{n+1}^2/32]\) is good for this angle, where good means sine distance at least \(1/P\) from *all* lattice directions of norm at most \(2P\). A direction parallel to \(P_n\) is covered by this lower estimate. At any other \(w\) with \(\|w\|\le\sqrt{2P}\le h_{n+1}/4\), integrality of its determinant with \(P_n\) gives
\[
 |\sin(\arg w-\phi)|\ge{1\over h_n\|w\|}-{2\over h_n h_{n+1}}
     \ge{1\over 2h_n\sqrt{2P}}\ge 1/P.
\]
At every step the good log-interval occupies at least the fraction \(b=(1-\log64/\log100)/2>0\) of \([2\log h_n,2\log h_{n+1}]\). Those step lengths are at most \(2\log(2D+1)\) for \(n\ge1\), and the finitely many initial \(h_1\) are fixed. Thus for each sufficiently large \(T\), each limit angle is good on a set of measure at least \(bT/2\) of the log-parameters in \([T,2T]\). By Fubini there is a \(P\) in the resulting range with \(\mu(G_P)\ge a=b/2\), where \(G_P\) is the closed set of good angles. Its definition is invariant under quarter-turns. Choose such a sequence \(P\to\infty\). Any physical label prime in \([P,2P]\) passes the advertised finite test on \(G_P\).

The preceding angular-sum construction applies uniformly to the normalized restrictions \(\nu=\mu|_{G_P}/\mu(G_P)\). We supply the needed uniform lower density since the locations of these sets vary. For any probability \(\nu\le\mu/a\), one has \(\nu^{*g_0}\le a^{-g_0}\mu^{*g_0}\), so Plancherel gives the uniform bound \(\sum_k|\widehat\nu(k)|^{2g_0}\ll_a1\). Also \(\sup_{k\ne0,\nu\le\mu/a}|\widehat\nu(k)|<1\). Indeed the images of \(\mu\) under multiplication by \(k\) converge weakly to Haar as \(|k|\to\infty\), by Fourier decay, and each of the finitely many remaining images is atomless. Consequently there is \(t=t(a)>0\) for which every arc of radius \(t\) has mass less than \(a/2\) in every such image: for large \(|k|\) use a fixed finite mesh of slightly larger arcs and Haar convergence, then decrease \(t\) for the finitely many other \(k\). In each image of \(\nu\) this gives arc mass at most \(1/2\). Projecting onto the direction of its Fourier average bounds the absolute Fourier coefficient by \((1+\cos(2\pi t))/2<1\). Increasing the fixed block size \(g\) as needed, \(\sum_{k\ne0}|\widehat\nu(k)|^g<1/2\) uniformly. Thus the sum of \(g\) independent samples from \(\nu\) has continuous density at least \(1/2\). The two joint upper bounds in the conditional-measure proof are inherited directly from \(\nu^g\le a^{-g}\mu^g\) before passing to the shrinking intervals. All constants are independent of \(P\). The weak limit is supported in \(G_P^g\), because this set is closed. The same mixed-gap union bound now works.

Finally, along the good \(P\)'s choose a Paley prime \(q\equiv3\pmod4\) in \([P^{1/d},2P^{1/d}]\) by the prime number theorem for progressions. Then \(P\asymp q^d\), which is all the construction requires. The finite sharp test holds simultaneously with all the earlier conditions. \(\square\)
