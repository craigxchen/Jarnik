# Round 29, lane spectral: the Li–Liu template and lattice points on short arcs

## 0. Verdict

* **Li–Liu (task 1).** The paper is a correct but formal reformulation. Its hypotheses S and R are, by the paper's own Thm 6.1 and (2.28)/(2.30), the operator-norm bounds $\|T_N\|\le\mathfrak S/(2h\sqrt{1-\eta})$ and $\|T_N\|\le\delta_N$, which are stronger than the conclusion. GRH is used only for holomorphy (contour shifts). The GRH bound $\|T_N\|\ll(\log N)^2$ is weaker than the trivial unconditional bound $\ll\log N$. **No step produces cancellation.** The coordinator's reading is confirmed. Numerically (evidence), the residual kernel does cancel to $\sim10^{-2}$ at $N\approx2\cdot10^5$, so R is not refuted; it is a uniform square-root Hardy–Littlewood statement far beyond GRH or pair correlation.
* **Transplant (task 2).** Three versions were built: H (Hecke angle characters), Q (auxiliary modulus on the strip), and A (Li–Liu verbatim on $N=a^2+b^2$). Each has an exact projection, a diagonal, a model and the quadratic-form identity.
* **Decisive negatives (task 3).**
  * Prop. 2.0: the arc problem has only a sub-polynomial margin, so every $Q^\varepsilon$-lossy analytic input is vacuous.
  * Prop. 2.1: H is fake-blind (F1).
  * Lemma 2.2 / Prop. 2.4: Q is Parseval of an elementary equal-norm pair count on neighbouring circles. Its uniform (operator-norm) version is **false** for $W\ge5R$: rational-direction resonances give $\sup E^W\ge c_C\sqrt R$ at non-trivial positions. Small windows reduce to the original problem.
  * Prop. 2.6: A is exactly diagonal. Off-diagonal families are empty, and the residual is the count itself.
  * Lemma 2.7: the arc localisation costs $\|v\|^2\asymp R^{1/2}/C$.
* **Outcome.** No hypothesis that is not a restatement (GRH, pair correlation, ratios, Gaussian Elliott–Halberstam, sup-norm, subconvexity) implies the uniform bound or a growth improvement through this template. The only new conditional statement (Remark 3.1, $M\ll C\sqrt{\log R}$ from a non-resonant uniform pair count) is a restatement. The growth rate and the uniform bound are unchanged.

## 1. The Li–Liu paper: what is proved, what is assumed

*Source.* J. Li, J. Liu, "A conditional proof of the Goldbach conjecture", 63 pp., the user's upload (read in full, pages 1–63).
Page and equation numbers below are the paper's.

### 1.1 What the paper does, step by step

Fix $\phi\in C_c^\infty$, $0\le\phi\le1$, supported in $[1/3,2/3]$, $\phi(u)=\phi(1-u)$; $\mathcal G(N)=N^{-1}\sum_{n+m=N,(nm,N)=1}\Lambda(n)\Lambda(m)\phi(n/N)\phi(m/N)$, $c_\phi=\int\phi^2$.

1. **Finite projection (Lemma 3.1; unconditional, elementary).** With $X_N^-$ the odd characters mod $N$,
   $\frac{2}{\varphi(N)}\sum_{\chi\in X_N^-}\chi(n)\bar\chi(m)=\mathbf 1_{(nm,N)=1}(\mathbf 1_{n\equiv m}-\mathbf 1_{n\equiv-m})$, and the support of $\phi$ turns the two congruences into $n=m$ and $n+m=N$. Hence $\mathcal G=\mathcal D-\mathcal P$ with $\mathcal D=N^{-1}\sum\Lambda^2\phi^2$ (diagonal) and $\mathcal P=\frac2{\varphi(N)}\sum_{\chi}|P_\chi|^2$, and $0\le\mathcal P\le\mathcal D\ll\log N$.
2. **Mellin form (Lemma 3.5, Prop. 3.2).** $\mathcal P(N)=(2\pi i)^{-2}\iint N^{s+z}\tilde\Phi(s)\tilde\Phi(z)\mathscr P_N(s,z)$ with $\mathscr P_N=\frac2{\varphi(N)}\sum_\chi\frac{L'}{L}(\tfrac12+s,\chi)\frac{L'}{L}(\tfrac12+z,\bar\chi)$. **GRH($N$) is used only here**: it makes $L'/L(\frac12+s,\chi)$ holomorphic on $\Re s>0$, so the contours move to $0<\sigma\le1/9$.
3. **Model (§4, Lemmas 4.1–4.3, Prop. 4.4; unconditional).** The off-diagonal families $m=n+kN$, $n=m+kN$, $n+m=kN$ are replaced by their Hardy–Littlewood densities $\mathfrak S(kN)$ times Beta integrals; $\sum_k\mathfrak S(kN)k^{-w}=\mathfrak S(N)\zeta(w)P_N(w)$, giving the continued model $\mathscr K_N(s,z)=\mathfrak S(N)(\pi/N)^{s+z}P_N(s+z)\zeta(1-s-z)\mathcal X_1(s)\mathcal X_1(z)$, whose weighted integral is exactly $-c_\phi\mathfrak S(N)$ (Lemma 7.2).
4. **Residual and operator (Props. 4.5, 5.3).** $\mathscr E_N=\mathscr P_N-\mathscr D_N-\mathscr K_N$ is holomorphic on $0<\Re s,\Re z<3/2$ under GRH; $T_N$ is the integral operator with kernel $\mathscr E_N(s,\bar z)/((1+s)(1+\bar z))$ on the Hardy space $H^2$ of the strip $1/\log N<\Re s<1/8$. It is selfadjoint and Hilbert–Schmidt.
5. **Quadratic-form identity (Thm 7.1).** With $v_N$ the Riesz vector of $\ell_N(h)=\frac1{2\pi i}\int N^s\tilde\Phi(s)(1+s)h(s)\,ds$, $\|v_N\|^2\ll_\phi1$ and $c_\phi\mathfrak S(N)-\mathcal G(N)=\langle T_Nv_N,v_N\rangle$; hence $|\mathcal G-c_\phi\mathfrak S|\ll_\phi\|T_N\|$.
6. **Spectral function (§6).** $\mathfrak X_N(q)=\prod_j(1-4\mu_j^2\delta_N^{-2}q(1-q))=\det(I-4\delta_N^{-2}q(1-q)T_N^2)$, $\delta_N^2=\mathfrak S(N)(\log N)^3/N$. Solving the quadratic factor by factor (Prop. 6.3) gives the "zero–norm dictionary" (Thm 6.1): $\mathfrak X_N\ne0$ on $\Re q>1-\eta$ $\iff\|T_N\|\le\delta_N/(2\sqrt{\eta(1-\eta)})$, and "all zeros on $\Re q=\frac12$" $\iff\|T_N\|\le\delta_N$.
7. **What GRH alone gives (Thm 6.2, Prop. 6.5).** $\|T_N\|\le\|T_N\|_{HS}\ll(\log N)^2$, by bounding $\mathcal P_N,\mathcal D_N,\mathcal K_N$ *separately* (the paper says so: "they do not measure cancellation between the components", p. 12).
8. **Hypotheses and theorems.** Hypothesis S: no zeros of $\mathfrak X_N$ in $\Re q>1-\eta_N$, $\eta_N=h^2(\log N)^3/(N\mathfrak S(N))$; by step 6 this is **verbatim** $\|T_N\|\le\mathfrak S(N)/(2h(N)\sqrt{1-\eta_N})$ ((2.28)). Hypothesis R: all zeros on $\Re q=\frac12$, i.e. **verbatim** $\|T_N\|\le\delta_N$ ((2.30)). Thm 1.1 (GRH + S): $\mathcal G=c_\phi\mathfrak S+O(\mathfrak S/h)$. Thm 1.2 (GRH + R): error $O(\delta_N)$.
9. **Appendices.** A: the left edge $\alpha=c/\log N$ is a free parameter; Prop. A.4 proves that for fixed $N$, $\|T_{N,c}\|\ge z_N/c\to\infty$ as $c\downarrow0$, so R$(N,c)$ is false for small $c$. B: at $N=30030$ the character integrals recover the exact count 905 (unconditional finite orthogonality), and a discretised operator has $\|C\|/\delta_N\approx0.83$; the authors state this does not certify R.

### 1.2 Assessment

**(a) The deduction is correct and formal.** Steps 1–6 are exact identities. I found no gap in them, and none of them estimates anything.

**(b) The hypotheses are restatements of the needed bound, not new inputs.** The zeros of $\mathfrak X_N$ are by construction the numbers $\frac12\pm\frac12\sqrt{1-\delta_N^2/\mu_j^2}$. Statements about their location are statements about $\max|\mu_j|=\|T_N\|$; the paper itself calls them equivalences ((2.28), (2.30), Thm 6.1). The quaternionic extension, slice regularity and Fueter regularity play no role in any estimate (the paper says the real coefficients extend the scalar function to other slices; nothing is deduced from it).

**(c) The hypothesis is strictly stronger than the conclusion.** Thm 1.1's conclusion is the value of the quadratic form at one vector $v_N$; the hypothesis bounds it at every unit vector. In particular (proof: the argument of Prop. B.4, which works for any weight in $C_c^6(0,1)$ with $\tilde H_f(1)=0$), S$(N;h)$ gives the same relative error simultaneously for every admissible weight, and, since $\|\ell_X\|^2\asymp(X/N)^{2/\log N}$, for the odd-character moment at every scale $X\le N^{O(1)}$, i.e. for Hardy–Littlewood prime-pair sums along all shifts $kN$, $k\le X/N$. So the hypothesis is a uniform square-root-type Hardy–Littlewood statement for the progression of shifts $N\mathbb Z$, of which binary Goldbach at $N$ is one coordinate.

**(d) GRH supplies no cancellation.** GRH enters only through holomorphy of $L'/L$ (contour shifts) and the local zero count. The GRH bound $\|T_N\|\ll(\log N)^2$ gives $|\mathcal G-c_\phi\mathfrak S|\ll(\log N)^2$, which is **weaker than the trivial unconditional bound** $|\mathcal G-c_\phi\mathfrak S|\le\mathcal D+c_\phi\mathfrak S\ll\log N$ from step 1. The averaged bound for $L'/L$ (Lemma 6.4) is the size of the diagonal of a Dirichlet polynomial (orthogonality), again no cancellation. I found **no step in which GRH, or anything else, produces cancellation** between $\mathscr P_N$, $\mathscr D_N$ and $\mathscr K_N$. The non-circular content is: the exact finite identity (classical orthogonality), the exact evaluation of the Hardy–Littlewood main term as a continued Beta/zeta model (a clean bookkeeping of the HL singular series average over $kN$, of the type $\sum_k\mathfrak S(kN)k^{-w}=\mathfrak S(N)\zeta(w)P_N(w)$), and the sensitivity result A.4.

**(e) How strong is R?** At a real point $s=\sigma$, $\mathscr K_N(\sigma,\sigma)=\mathfrak S(N)P_N(2\sigma)(\pi/N)^{2\sigma}\zeta(1-2\sigma)\mathcal X_1^2\approx-e^{-2a}/(4\sigma^2)$ with $a=\sigma\log N$, while $\mathscr D_N(2\sigma)\approx1/(4\sigma^2)$. This is the diagonal/off-diagonal split of the Goldston–Gonek–Montgomery mean value of $|L'/L|^2$ near the critical line, which for $\zeta$ is equivalent to Montgomery's pair-correlation conjecture (cited from memory; see §3.4). So $\mathscr E_N$ is the error term in a *family-averaged* pair-correlation law at **fixed height**, and R asks this error, uniformly over the Hardy space, to be $\le\delta_N\approx N^{-1/2}(\log N)^{3/2}$: square-root cancellation over the $\varphi(N)/2$ characters. Heuristically (independent characters, GUE statistics) this is the expected order, so R is not obviously false. But it lies far beyond GRH and pair correlation: it is a power-saving, uniform-in-test-vector Hardy–Littlewood statement. My own numerical check of the residual kernel is in §1.3.

**Verdict on the paper.** The coordinator's reading is correct. The conditional theorem is a reformulation: "if the residual operator, whose quadratic form at $v_N$ is the Goldbach error, has small norm, then the Goldbach error is small." Nothing in it is a mechanism by which GRH, or any standard conjecture, gives the needed cancellation. What survives as a tool is the **template**: exact projection, exact diagonal, exact continuous model, and the residual as a quadratic form. §2 transplants that template.

### 1.3 Numerical check of the residual kernel (evidence)

`spectral_checks/ll_kernel.py` computes $\mathscr E_N(\sigma,\sigma)=\mathscr P_N-\mathscr D_N-\mathscr K_N$ at real points for $N=2p$. $\mathscr P_N$ is computed from $L'/L(\frac12+\sigma,\chi)$ for all odd $\chi$ (Hurwitz $\zeta$ with python-flint, FFT over discrete logarithms). $\mathscr D_N$ comes from prime-zeta derivatives, and $\mathscr K_N$ from the paper's closed form (2.12). Output is in `ll_kernel_out.txt`. All values are double precision.

| $N$ | $\sigma$ | $\mathscr P$ | $\mathscr D$ | $\mathscr K$ | $\mathscr E$ | $\delta_N$ |
|---|---|---|---|---|---|---|
| 2018 | 1/9 | 12.60 | 19.11 | −6.75 | 0.237 | 0.537 |
| 20014 | 1/9 | 15.30 | 19.12 | −4.05 | 0.234 | 0.253 |
| 60022 | 1/9 | 15.89 | 19.12 | −3.17 | −0.059 | 0.171 |
| 200006 | 1/9 | 16.68 | 19.12 | −2.43 | −0.015 | 0.110 |
| 200006 | 0.05 | 59.52 | 98.47 | −38.99 | 0.037 | 0.110 |

The three components are each of order 10–100, and they cancel to $|\mathscr E|\approx10^{-2}$. So the Hardy–Littlewood model is numerically right, and the cancellation that R needs is visible. The pointwise lower bound $\|T_N\|\ge8d(s)|\mathscr E(s,\bar s)|/|1+s|^2$ (from (5.8)) gives values $\le0.015$, below $\delta_N$, so R is not contradicted in this range. The proved GRH bound $(\log N)^2\approx150$ is $10^4$ times too large. *Evidence only:* these are single real points, not operator norms.

## 2. Transplanting the template to short arcs

Notation. $N=R^2$; $\Theta=C\sqrt R/R=CR^{-1/2}$ is the angular width of an arc of length $C\sqrt R$; $M_N(I)$ is the number of lattice points of $|z|^2=N$ with $\arg z\in I$, $|I|=\Theta$. "Non-trivial position" means $I$ at angular distance $\ge2\Theta$ from the 8 axis/diagonal directions.

### 2.0 A scale obstruction that applies to every analytic transplant

**Proposition 2.0 (proved; elementary).** $M_N(I)\le r_2(N)\le4\tau(N)\le\exp\big((\log2+o(1))\log N/\log\log N\big)=R^{o(1)}$. The record's Theorem 1.1 gives $M\le(\frac23+\varepsilon)\log R/\log\log R$. Hence an upper bound for $M$ (or $M^2$, or $M(M-1)$) of the form $B(Q)$ is weaker than what is already proved whenever $\log B(Q)\gg\log\log R$. Here $Q$ is any analytic conductor with $\log Q\asymp\log R$, for example $B=Q^\varepsilon$ or $B=\exp(c\log Q/\log\log Q)$.

*Proof.* The first two inequalities are standard: $r_2(N)=4\sum_{d|N}\chi_4(d)$, together with Wigert's bound for $\tau$. The last claim is a comparison of sizes. $\square$

*Consequence.* Goldbach has a margin of a full power of $N$: the main term is $\asymp1$ after normalisation, while the trivial error is $\log N$ and R asks for $N^{-1/2+o(1)}$. So Li–Liu can afford $N^\varepsilon$ losses. The arc problem lives entirely inside a **sub-polynomial window**. Every analytic input whose error term carries a conductor-power loss $Q^\varepsilon$ (GRH bounds for non-prime-weighted sums, Lindelöf, subconvexity, CFKRS moments and ratios with $O(Q^{-1/2+\varepsilon})$, sup-norm bounds with $\lambda^\varepsilon$) is therefore vacuous for the growth question before any structure is examined. This covers GRH used through Lindelöf, which loses $\exp(c\log Q/\log\log Q)$ (memory-quote: Littlewood's conditional bound; Chandee–Soundararajan for the constant). Only inputs with polylogarithmic error terms can matter. In the analytic world these are sums over primes (log-derivatives), as in Li–Liu. The arc count, however, is a count over all lattice points: its generating function is $L(s)$, not $L'/L(s)$.

### 2.1 Transplant H: Hecke Grössencharacters (the angle projection)

Let $\xi(z)=(z/|z|)^4$, a Hecke character of $\mathbb Q(i)$. For $N$ odd with all primes $p\equiv1\ (4)$ and $p^{e_p}\|N$, $\pi_p=\sqrt p\,e^{i\phi_p}$:
$$\lambda_k(N):=\sum_{|z|^2=N}\xi^k(z)=4\prod_{p|N}U_{e_p}(\cos4k\phi_p),\qquad U_e(\cos x)=\frac{\sin((e+1)x)}{\sin x}.$$
This is the $N$-th coefficient of $L(s,\xi^k)$, i.e. of the CM form $\tfrac14\sum_aa^{4k}q^{N(a)}$. (The units contribute $i^{4kj}=1$, and $\sum_{a=0}^{e}e^{4ik(2a-e)\phi}=U_e(\cos4k\phi)$.)
* **Projection.** For any $F\ge\mathbf 1_I$ with $\hat F$ supported in $|k|\le K$, $M_N(I)\le\sum_{|k|\le K}\hat F(k)\lambda_k(N)$. For $F$ of Beurling–Selberg type, $K\asymp\Theta^{-1}=\sqrt R/C$.
* **Diagonal.** The $k=0$ term $\hat F(0)\,r_2(N)\approx(\Theta/2\pi)\,r_2(N)$.
* **Model.** The local densities give again $(\Theta/2\pi)r_2(N)$.
* **Residual.** $\sum_{0<|k|\le K}\hat F(k)\lambda_k(N)$. As a quadratic form it is $\langle Tv,v\rangle$ with $v=(\hat F(k)^{1/2})_k$ and $T=\mathrm{diag}(\lambda_k(N))$, or any Gram rearrangement of it.

**Proposition 2.1 (proved; fake-blindness of H).** Every object in transplant H is a fixed polynomial in $(\cos4k\phi_p,\sin4k\phi_p)_{p|N,|k|\le K}$ with the exponents $e_p$ as parameters. The identity "count $=\sum_k\hat F(k)\lambda_k$" is the Fourier expansion of $\mathbf 1_I$ over the multiset of angles $\{\sum_p(2a_p-e_p)\phi_p\}$, so it holds verbatim for a fake circle (genuine $p,e_p$, arbitrary real $\phi_p$). Hence any bound on the residual, the operator or its spectrum that is derived for all real circles, and is a function of these data alone, either holds for fakes or uses the fact that the $\phi_p$ are Gaussian-prime angles.

*Proof.* The product formula is multiplicativity plus the Chebyshev sum. Nothing else enters. $\square$

GRH for $L(s,\xi^k)$ (and twists) constrains only sums $\sum_{n\le X}$ or $\sum_{p\le X}$. Changing the angles of the $\le\omega(N)\le\log N$ primes of $N$ alters any $\sum_{p\le X}w(p)\xi^k(\pi_p)\log p$ with $|w|\le1$ by at most $2\log^2N$, which is below the GRH error $X^{1/2}$ as soon as $X\ge\log^4N$. So GRH-level information about prime angles is graft-insensitive, the mechanism of record barrier item 2. *Label:* this is a structural argument, not a theorem about all possible uses of GRH. **Outcome under F1: transplant H is fake-blind.**

### 2.2 Transplant Q: Dirichlet characters of $|z|^2$ to an auxiliary modulus (the strip)

Let $S_W(I)=\{z\in\mathbb Z[i]:|\,|z|^2-N|\le W,\ \arg z\in I\}$, with about $\Theta\cdot2\pi W/(2\pi)=\Theta W$ points, and let $q>2W$.

**Lemma 2.2 (exact projection = Parseval; proved).** Let $r_I(n)=\#\{z\in S_W(I):|z|^2=n\}$ and $\Theta_a=\sum_{z\in S_W(I)}e(a|z|^2/q)$. Then
$$M_N(I)=\frac1q\sum_{a\bmod q}e(-aN/q)\Theta_a,\qquad \frac1q\sum_{a\bmod q}|\Theta_a|^2=\sum_nr_I(n)^2=\#S_W(I)+E^W(I),$$
with $E^W(I)=\sum_nr_I(n)(r_I(n)-1)\ge M_N(I)(M_N(I)-1)$. With multiplicative characters, $\frac1{\varphi(q)}\sum_\chi|S_\chi|^2$ gives the same identity on norms coprime to $q$. The "Li–Liu decomposition" here reads: diagonal $=\#S_W(I)$, off-diagonal $=E^W(I)$.

*Proof.* Orthogonality, and $|n-n'|\le2W<q$. $\square$

**Lemma 2.3 (structure of the off-diagonal; proved).** Distinct $z,w$ with $|z|=|w|$ satisfy $w=z\,u\,P/\bar P$ with $u$ a unit and $P$ a primitive Gaussian integer, unique up to sign; $w\in\mathbb Z[i]$ iff $\bar P'\mid z$, where $P'=P/\gcd(P,\bar P)$. Hence
$$E^W(I)=\sum_{u,\pm P}\#\{z\in S_W(I):\bar P'\mid z,\ zuP/\bar P\in S_W(I)\}.$$
For non-trivial $I$, only $u=1$ and $|\arg P|\le\Theta/2$ contribute. So $E^W$ is a sum of lattice-point counts of the sublattices $\bar P'\mathbb Z[i]$ in a thin curved strip. No $L$-function and no character is needed. The characters mod $q$, and the $L$-functions $L(s,\chi\circ\mathrm N\cdot\xi^k)$ that a Mellin–Hecke expansion of $S_\chi$ would produce, re-encode this elementary count.

*Model (exact average; computed).* Averaging over positions, $\mathbb E_\theta E^W=\sum_{n}\sum_{z\ne w,|z|^2=|w|^2=n}(\Theta-\angle(z,w))_+/2\pi$. Exact enumeration (`pairs_window.py`) gives $\mathbb E_\theta E^W=\kappa\,(W/R)\,C^2\log R$ with $\kappa\in[0.25,0.33]$ on five instances ($R\approx7\cdot10^3$–$5\cdot10^4$, $C=4,8$, $W=R,R/10$). This matches the heuristic $\sum_{P}(\Theta-2|\arg P|)_+/|P|^2\asymp\Theta^2\log R$ (**exact** counts; constant **numerical**).

**Proposition 2.4 (rational-direction resonance; proved).** Let $W\ge5R$, $C>0$, $R\ge10$. Put $g=1+2i$, $A=\lfloor\sqrt{N/5}\rfloor$, $z_s=g(s+iA)$ and $z'_s=g(-s+iA)$. Then $|z_s|^2=|z'_s|^2=5(A^2+s^2)\in[N-W,N+W]$ for $1\le s\le\sqrt R$. All pairs with $2s\sqrt5\le C\sqrt R$ lie in one arc of angular width $\Theta$ centred at $\arg(ig)=\pi-\arctan\frac12$, a non-trivial position at angle $\ge18^\circ$ from every axis and diagonal. Hence
$$\sup_{I\ \text{non-trivial}}E^W(I)\ \ge\ 2\big\lfloor\min(C/(2\sqrt5),1)\sqrt R\big\rfloor\ \gg_C\ \sqrt R.$$

*Proof.* $5A^2\in(N-2\sqrt{5N}-5,N]$ and $5s^2\le5R\le W$, so the norms lie in the window. $z_s-z'_s=2sg$, so $|z_s-z'_s|=2s\sqrt5$. Both points have argument within $s\sqrt5/R$ of $\arg(ig)$. $\square$

Verified with exact integer arithmetic for $N=10^8+1,\ldots,10^{16}+61$ (`resonance_and_norm.py`): the count is $1.789\sqrt R=C\sqrt R/\sqrt5$. For $W=R$ the exact sweep gives $\sup E^W\approx5\times$ the mean (238 against 46 at $N=2576450045$, $C=4$; the sector width there also admits other resonant directions), while the central circle has $M(M-1)=12$.

**Corollary 2.5 (window dichotomy).**
* (a) **Proved.** For $W\ge5R$, the residual "$E^W(I)-$model" is not uniformly $O(\mathrm{polylog}R)$ over non-trivial positions. Any operator-norm-type hypothesis for transplant Q, i.e. a bound uniform in the test vector (= arc position), is **false**.
* (b) **Heuristic.** The same lines with $|g|^2=m\approx(R/W)^2$ give $\sup E^W\gtrsim W^{3/2}/R$ for $R^{2/3}\ll W\le5R$.
* (c) **Proved, trivial.** For $W\le R^{2/3}$ the model is $O(R^{-1/3}\log R)$. A bound "$E^W\le$ model $+$ residual" then has residual $\ge M_N(M_N-1)$, so it is the original problem.

So the neighbouring circles enter only monotonically ($E^W\ge M(M-1)$). Large windows are dominated by rational-direction (Jarník-line) resonances, which no hypothesis can remove uniformly. Small windows bring nothing. A non-uniform hypothesis that subtracts the resonances as "major arcs" is still a hypothesis on $\sup_IE^W(I)$ over minor-arc positions, and the cluster's own pairs are in it, so it is a restatement.

### 2.3 Transplant A: Li–Liu verbatim (characters mod $N$ on $N=a^2+b^2$)

For an arc at non-trivial position, the points are $a+bi$ with $a\in J_1$, $b\in J_2$, intervals of length $\ell\asymp C\sqrt R\asymp CN^{1/4}$. Primitive points have $(ab,N)=1$, and $a^2+b^2\in(0,2N)$, so
$$M_N(I)=\frac1{\varphi(N)}\sum_{\chi\bmod N}\chi(-1)A_\chi\overline{B_\chi},\qquad A_\chi=\sum_{a\in J_1}\chi^2(a),\ B_\chi=\sum_{b\in J_2}\chi^2(b).$$
These are character sums of length $N^{1/4}$ to modulus $N$. The critical arc length is exactly the Burgess threshold, a curiosity but a real one.

**Proposition 2.6 (degeneracy; proved).** If $2C\sqrt R\cdot2\sqrt N<N$ (true for $R>16C^2$), then
$$\frac1{\varphi(N)}\sum_\chi|A_\chi|^2=\#\{a\in J_1:(a,N)=1\},$$
exactly. More generally every second moment of $A_\chi\pm B_\chi$, $A_\chi\pm iB_\chi$ equals the diagonal plus the target count itself. Li–Liu's off-diagonal families (which the Beta model approximates) are **empty**. There is nothing to model, and the "residual" is the count.

*Proof.* $a^2\equiv a'^2\pmod N$ means $N\mid(a-a')(a+a')$, with $0<|a-a'|\le\ell$ and $a+a'\le2\sqrt N$. $\square$

To recover Li–Liu's non-trivial structure one must drop the arc localisation from the weights, as Li–Liu do with $\phi(n/N)$. Re-imposing it costs the test-vector norm:

**Lemma 2.7 (uncertainty loss; proved, leading order).** In Li–Liu's Hardy space, the functional of a weight $\Phi(u)=\sqrt u\,\phi((u-1)/\delta)$ of relative width $\delta$ has $\|\ell\|^2=\frac{e^{2}}{2\pi}\delta^{-1}\|\phi'\|_2^2(1+O(\delta)+O(1/\log N))$, by (7.16) and substituting $u=1+\delta x$ in $\int|\Phi-u\Phi'|^2du/u$. For an arc, $\delta\asymp\Theta=CR^{-1/2}$, so $\|v\|^2\asymp R^{1/2}/C$: the same $N^{1/4}$ as the diagonal in Lemma 2.2. Numerically $\delta\int|V|^2du/u\to\|\phi'\|^2=0.40959$ for $\delta=0.1\to0.003$ (`resonance_and_norm.py`).

So in every version the arc localisation puts a factor $\asymp N^{1/4}$ (the number of lattice points of the strip) between the operator norm and the count.

## 3. What GRH and other hypotheses give; filters F1, F2

**3.1 GRH gives holomorphy, not size.** This holds in Li–Liu (§1.2(d)), and it holds in every transplant. In H and Q the relevant $L$-functions are $L(s,\xi^k)$ and $L(s,\chi\circ\mathrm N\cdot\xi^k)$, $\chi\bmod q$, with analytic conductor $Q$ satisfying $\log Q\asymp\log R$. Under GRH the individual sums satisfy $|S_\chi|\ll(\#S_W)^{1/2}Q^{o(1)}$. Summing squares gives $\le\#S_W\cdot Q^{o(1)}$, which is above the diagonal: no information on $E^W$ and, by Prop. 2.0, nothing even before that. The family second moment is the exact identity of Lemma 2.2, so no analytic input is needed for it and none can improve it.

**3.2 Hypotheses that are not restatements: outcome.**

| hypothesis | what it would control here | outcome |
|---|---|---|
| GRH for $L(s,\xi^k)$, $L(s,\chi\circ\mathrm N\,\xi^k)$, Dirichlet $L$ mod $N$ | holomorphy; prime-angle sums with error $X^{1/2}$ | no size information on the residual; graft-insensitive for $X\ge\log^4N$ (§2.1); vacuous by Prop. 2.0 when used through Lindelöf |
| Lindelöf, subconvexity, sup-norm/$L^p$ of toral eigenfunctions | individual sums with $Q^\varepsilon$ or $\lambda^\varepsilon$ loss | vacuous (Prop. 2.0); the $\varepsilon$-free $L^p$ statements are additive-energy statements, which fakes satisfy (F1) |
| pair correlation, GGM-type mean values, ratios/CFKRS moments | family averages at a fixed test function, relative error $O(Q^{-1/2+\varepsilon})$ | for Q with $W\asymp R$ the main term is $\#S_W\asymp C\sqrt R$, so the error is $\asymp CR^{\varepsilon}$: vacuous for growth; the $\varepsilon$-free, position-uniform version is false (Prop. 2.4) |
| Elliott–Halberstam for $\mathbb Z[i]$ | Gaussian primes in sectors/progressions, averaged over moduli | averaged input (record barrier item 2); the arc count is not prime-weighted |
| Li–Liu-type norm bounds (S, R analogues) | residual uniformly over test vectors | restatement (residual at the cluster vector $\ge M(M-1)-$model); false in the uniform form for $W\ge5R$; for small $W$ it is the original problem |

**No hypothesis in the list implies the uniform bound or any $M=o(\log R/\log\log R)$ through this template.** The proofs are Props. 2.0, 2.1, 2.4, 2.6 and Lemma 2.7. The judgements "restatement" and "averaged" are the classification used in the record. The GRH row also rests on the structural argument of §2.1.

**Remark 3.1 (implication proved; hypothesis a restatement).** If for some $W\in[R/10,R]$ one had $E^W(I)\le K(W/R)C^2\log R$ for every position $I$ that carries a cluster, then Lemma 2.2 gives $M\le1+\sqrt K\,C\sqrt{\log R}=o(\log R/\log\log R)$. However:
* the hypothesis contains $M(M-1)$ itself;
* it fails at resonant positions (Prop. 2.4, Cor. 2.5(b));
* a version restricted to "non-resonant" positions is implied by no standard conjecture I know of.

I record it to make the gap explicit: the position-averaged pair count on neighbouring circles is $\asymp\log R$ (exact counts), so this heuristic sits at $M\approx\sqrt{\log R}$.

**3.3 Filter F1.**
* **H:** fake-blind (Prop. 2.1).
* **Q:** not fake-blind, since it uses neighbouring circles and real lattice geometry. But it enters only through the monotone inequality $E^W\ge M(M-1)$, and its uniform form is false (Prop. 2.4).
* **A:** degenerate (Prop. 2.6); its only content is the count.

**3.4 Filter F2.** The template's one inequality is Cauchy–Schwarz $|\langle Tv,v\rangle|\le\|T\|\|v\|^2$. At the cluster vector the quadratic form *equals* the error, so the method is exactly critical in the trivial sense: it neither gains nor loses. Without the norm hypothesis it proves nothing, and so nothing that would extend to the false 4-point analogue of Question B.

**3.5 Literature.** Goldston–Gonek–Montgomery, J. reine angew. Math. 537 (2001) 105–126: under RH, pair correlation is equivalent to an asymptotic for $\int|\zeta'/\zeta(\frac12+a/\log T+it)|^2dt$. This was seen through secondary sources only (search results: arXiv 2211.08571, Gonek's lecture notes arXiv math/0401126); the full text was not read. The conditional Lindelöf constants (Littlewood; Chandee–Soundararajan) are memory-quotes.


## 4. Status labels and files

| item | label |
|---|---|
| §1.1–1.2 account; GRH bound weaker than trivial; S, R are norm bounds | proved (read from the paper; the trivial bound is Lemma 3.1 of the paper) |
| §1.2(c) hypothesis implies all-weights, all-scales HL statement | proved (Cauchy–Schwarz plus Prop. B.4's argument); scale range stated loosely |
| §1.3 kernel cancellation | numerical evidence (double precision, real points) |
| Prop. 2.0, 2.1, Lemma 2.2, 2.3, Prop. 2.4, 2.6, Lemma 2.7 | proved (elementary); 2.4 and 2.7 also checked exactly or numerically |
| model constant $\kappa\approx0.26$–$0.33$; $\sup E^W\approx5\times$ mean at $W=R$ | exact counts / numerical |
| Cor. 2.5(b) $W^{3/2}/R$ resonance scale | heuristic |
| GRH graft-insensitivity (§2.1); "no standard hypothesis" table | structural argument (not a theorem about all uses of GRH) |

Scripts in `spectral_checks/`:
* `ll_kernel.py` (output `ll_kernel_out.txt`);
* `pairs_window.py` (exact sweeps over sector positions);
* `resonance_and_norm.py` (output `resonance_and_norm_out.txt`).
