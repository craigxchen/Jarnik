**PROVED (a deeper low-height window).** There are nonarithmetic angle models with point residues through \(2m\), and also with the prime-angle spacing of VIII.4, such that
\[
 \limsup\frac{\log P_+(N')}{\log m}\le\frac{11}{6},\qquad
 \limsup\frac{W}{m\log m}\le\frac{22}{3}.
\]
Thus the earlier exponent \(23/12\) can be replaced by \(11/6\). This statement gives a fake for each fixed cutoff exponent strictly larger than \(11/6\).

*Proof.* Keep the Paley matrix, mark class and independent sampling with fixed probability \(\rho<1\) from VIII.3. Fix \(\tau<1/2\). For every primitive balanced integer word \(v\), put \(h=\|v\|_1\le q^\tau\), and let \(w=(Kv)_{S^c}\) be its coefficients on the retained complement. Simultaneously with the squared-norm conclusion of VIII.3, we may demand
\[
 \#\{i:w_i\not\equiv0\pmod n\}\ge 2c_\rho\sqrt q
 \qquad (n=4\ \text{or any odd prime}\ n\le h).
\]
Indeed, before sampling at least \(\sqrt q-h-1\) entries detect each modulus, by the square-root bound in VIII.3. Each is independently in the complement with probability at least \(1-\rho\). Choose for example \(c_\rho=(1-\rho)/10\). The lower-tail inequality gives failure probability \(\exp(-\Omega_\rho(\sqrt q))\). There are at most \((2q+1)^{q^\tau+1} q^\tau\) word/modulus pairs, so all tests and the previous norm tests pass for one choice. In particular \(\|w\|_2^2\ge c_\rho qh\), since \(\|v\|_2^2\ge h\).

Let \(A=q^a\) with any fixed \(a>0\). For independent discrete Gaussians \(X_i\) of width \(A\) as in VIII.3, we claim
\[
 \sup_b {\mathbb P}\Big(\sum_iw_iX_i=b\Big)\le q^{o(1)}(A\sqrt{qh})^{-1}
 \qquad(h\le q^\tau).
\]
Put \(F(t)=\sum_i\|w_i t\|_{\mathbb R/\mathbb Z}^2\), \(Z=(\log q)/A\), and use the representative \(0\le t<1\). Uniformly in \(u\), the one-variable characteristic function satisfies
\[
 |\psi_A(u)|\le(1+O(e^{-cA^2}))\exp(-c A^2\|u\|^2)
\]
for an absolute \(c>0\). This follows from the periodized Gaussian formula in VIII.3. On \(\|u\|\le1/4\) separate the nearest term, and on \(1/4<\|u\|\le1/2\) absorb the bounded prefactor by reducing \(c\). The factor to the power \(q\) is \(1+o(1)\). Outside \(F(t)\le Z^2\) the contribution to the integral for Fourier inversion is negligible. We bound the measure of that sublevel set.

There are \(\gg_\rho q\) nonzero \(w_i\) by the polynomial support modulo \(q\) in VIII.3. Call an index good at \(t\) if \(w_i\ne0\) and \(\|w_i t\|<1/(3h)\). At most \(9h^2 Z^2=o(q)\) nonzero indices are excluded. All good indices have a common nearest-integer ratio \(r=k_i/w_i\), because \(|k_iw_j-k_jw_i|<2/3\) follows from \(|w_i|\le h\). Write \(r\) in lowest terms with positive denominator \(d_0\). Every good \(w_i\) is a nonzero multiple of \(d_0\). It follows that
\[
 |t-r|\ll_\rho Z/(d_0\sqrt q).
\]
On every index with \(d_0\nmid w_i\), we have \(\|w_i r\|\ge1/d_0\) and \(|w_i(t-r)|\ll_\rho hZ/(d_0\sqrt q)=o(1/d_0)\). Hence \(F(r)\le2Z^2\) for all sufficiently large \(q\), since the other terms of \(F(r)\) vanish. For \(d_0=1,2\) this implies that every \(w_i\) is divisible by \(d_0\) once \(Z=o(1)\). Then \(F(t)=|t-r|^2\|w\|_2^2\), so these rational centers contribute total length \(O_\rho(Z/\sqrt{qh})\).

For each \(d_0>2\) with \(d_0\le h\), consider reduced numerators \(b_0\bmod d_0\) with \(F(b_0/d_0)\le2 Z^2\). At least \(N_0\ge c_\rho\sqrt q\) coefficients \(w_i\) are nonzero modulo \(d_0\), because \(d_0\) has either an odd prime factor or is divisible by four. For each such index, let \(d_i=d_0/\gcd(w_i,d_0)>1\). Under a uniform unit numerator \(b_0\) modulo \(d_0\), \(w_i b_0/d_0\) is uniform on the reduced fractions of denominator \(d_i\) modulo one, since reduction of units is onto and each fiber has equal size. For \(b_0\) satisfying the sublevel condition the product of \(\exp(-\|w_i b_0/d_0\|^2/Z^2)\) over those \(N_0\) indices is at least \(e^{-2}\). By Hölder its average is bounded by the geometric mean of the averages of the individual factors to the power \(N_0\). Each such average is at most
\[
 \frac2{\varphi(d_i)}\sum_{k=1}^{\lfloor d_i/2\rfloor}\exp[-N_0 k^2/(d_i Z)^2]
 \ll \frac{d_i}{\varphi(d_i)}\frac{Z}{\sqrt{N_0}}
 \le q^{o(1)} Z q^{-1/4}.
\]
For the integral comparison in the middle inequality it matters that \(k=0\) is absent. Thus there are at most \(d_0 q^{o(1)}Z q^{-1/4}\) possible reduced centers modulo one of that denominator.

Average over the choice of a good anchor index \(j\). Every \(t\) in this part of the sublevel set has \(\gg_\rho q\) such anchors. For fixed \(j\) its denominator \(d_0\) must divide \(|w_j|\le h\), which has \(q^{o(1)}\) divisors. At any divisor \(d_0>2\), the total length of the intervals of radius \(O_\rho(Z/(d_0\sqrt q))\) at the centers just counted is at most \(q^{o(1)} Z^2 q^{-3/4}\). The same bound holds for their union over divisors. Summing over at most \(q\) anchors and dividing by \(\gg_\rho q\) proves that the whole sublevel set has measure bounded by
\[
 O_\rho(Z/\sqrt{qh})+q^{o(1)}Z^2 q^{-3/4}.
\]
The standard divisor bound \(\#\{b:b\mid n\}=n^{o(1)}\) also bounds \(d_i/\varphi(d_i)\le2^{\omega(d_i)}\), uniformly here. Since \(h\le q^{\tau}<\sqrt q\), the second term is absorbed in \(q^{o(1)}/(A\sqrt{qh})\). Fourier inversion proves the stated atom bound, and specifying an integer moment digit in VIII.3 either has probability zero or costs one such atom. In particular this estimate needs no hypothesis \(a>\tau\).

For precision put \(d=2-\xi\), with \(0<\xi<1/6\). We may now choose \(\kappa>2\xi,\ \tau>\kappa+\xi,\ 0<a\) with \(\tau<1/2,\ 1.5+a+\kappa<d\), because at the lower endpoints both strict inequalities reduce to \(3\xi<1/2\). Use the \(\varepsilon,E,B,H,H_u\) of the proof of (VIII.10). At balanced primitive lengths \(q^\kappa<h\le q^\tau\) the improved atom is \(q^{-a-1/2+o(1)}h^{-1/2}\). Above \(q^\tau\) the same old Hölder estimate \(O_\rho(q^{-a-1/2})\) suffices. Writing \(\gamma=\log_q h\), the erasure count gives the same exponent per token \((\xi+\kappa+2\varepsilon)/2-3\gamma/4+o(1)\) in the first range and \((\xi+\kappa+2\varepsilon-\gamma)/2+o(1)\) in the second. Both are uniformly negative. The \(O(q^2)\) possible pairs of length and first surviving digit contribute \(o(1)\). Below \(H\), unbalanced words, nonprimitive words and the final \(h\ge q/2\) tail use no strong atom. Their inequalities in (VIII.10) hold unchanged since \(q^{1.5+a+\kappa+2\varepsilon}<P\). The weighted-bin estimate still has error \(O(q^{1.5-d}\log^{3/2}q)=o(1)\). For the residue character bound use \(\log L_0/\log P < .56\) instead of \(.53\). At \(8\le s\le q/2\), \(4/3-.56/2>1.05\); at \(s\ge q/2\), \((8/9)(2s+q/2)-.56s\ge(4/3-.56/2)q\); and always \(h-.56s\ge.72h\). Thus (VIII.9) still supplies the strict slack in the required sums. The spread variant requires the same rank estimate, \(|T|=O(q^{\xi/2})=o(\sqrt q)\), a fixed finite digit permutation and polynomial factors in \(q\), none affecting the preceding strict estimates. Lastly let \(\xi\uparrow1/6\) and \(\rho\uparrow1\) along a diagonal of sufficiently large Paley orders. The previous bound \(W/(m\log m)\le4d/\rho+o_\rho(1)\) gives the claimed coefficient. \(\square\)
