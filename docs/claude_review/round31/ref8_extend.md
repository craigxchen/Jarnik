# Referee report, block "extend": what VIII.1 buys, and whether it can be pushed

Referee: analyst, round 31, Part VIII of Astra's `answer.md`. In this block VIII.1 is taken as given, but I re-derived
it line by line (§2), because every extension below modifies the proof. Scripts are in `round31/ref8_extend/`. All
are my own code, and none imports Astra's scripts.

Labels: **Proved** (I checked every step), **Computed** (exact finite computation by my code), **Heuristic**,
**Gap**, **False**. "Record" means `docs/claude_review/`.

## 0. Results

1. **VIII.1 is correct (Proved, re-derived), and its constant can be improved (Proved).** The ramified prime divides
   every reduced chord, so the cross-pair average of $\log|b_{ij}|$ is at least $\frac12\log2$ for free. Hence:
   * no residues are needed for $C\le\sqrt2$, so that $M\le1+\omega(N')$ there;
   * one mod-2 indicator, of affine cost 1, suffices for $C<2^{3/4}$;
   * in general $K_C$ can be taken from the condition
     $\sum_{\ell\in S}\log\ell/(\ell+1)>\log C-\frac12\log2$. This is about 3.5 times smaller than Astra's
     $K_C$: for example $K_2=10$ instead of 21, and $K_8=1150$ instead of 3885.
   * Either way $K_C=C^{4+o(1)}$.
   * VIII.1 is an $A_1(Q_0)$ theorem, with $Q_0=\max S$: it uses only pair instances of the residue axiom (R)
     together with the conic fibre counts.
2. **Exact barycentre identity (Proved; it is the averaged pair inequality).** For any two probability measures
   $P,Q$ with disjoint supports on the marks,
   $$2\,\mathbb E_{P\otimes Q}\log|b_{ij}|-2\log C\ \le\ \sum_p e_p\log p\,\big[2t_p^2-2(s_p-\tfrac12)^2-2\,\mathbb E\min(X,Y)(1-\max(X,Y))\big]$$
   (notation in §2). VIII.1 is the case $t\equiv0$. Theorem 1.1's inert part and its part at split $p\nmid N'$ are the
   average of this inequality over random bisections. In Gram language, Theorem 1.1 is the certificate
   $\mathbf 1^TG\mathbf 1\ge0$. A weaker cut-matrix form of VIII.1 (rank counted in cuts rather than primes) is
   Perron–Frobenius simplicity for the same PSD matrix (Prop 2.4).
3. **(1) Combination with Theorem 1.1: the gain in the constant 3 is exactly zero (Proved, as a statement about
   inputs).**
   * There are data (an exponent profile with residue class counts) that satisfy inequality (5.1) of
     `two_thirds/proof.md`, have $\omega(N')\ge M$ with affinely independent rows, and have
     $W=3M\log M+O(M)$ (Prop 3.2).
   * Adding "at least $M-1-K_C$ primes with nontrivial cut" to the cut sum $\sum_p\sum_\tau(2S_{p\tau}-M)^2\log p$
     adds 0 for even $M$. For odd $M$ it adds only the term $\ge W'$ that the record already uses.
   * The two lower bounds are not additive, for the following reason. The arithmetic side of Theorem 1.1 lives at
     primes $q\le M$, while $W$ enters (5.1) only as a budget. The extra primes forced by VIII.1 can sit above $M$,
     with balanced cuts, and carry that budget at no cost.
   * What *does* move the constant is whether small split primes divide $N'$. The graded form
     $W'\ge3M\log M+2M\sum_{p\le M,\,p\nmid N'}\frac{\log p}{p-1}-O_C(M)$ is Proved (Prop 3.1). It is a
     rephrasing of the record's (5.1) and Lemma 4.3, and it has minimal novelty: `growth/auxiliary.md` already
     has the extreme case $4M\log M$.
4. **(2) Strengthening VIII.1.**
   * *Slack.* The barycentre step has slack only at primes whose cross-measure mean is unbalanced, or where both
     points sit at interior levels (Prop 4.1).
   * *Weighted-rank statements.* These reduce to the trivial count, because every prime has weight
     $\ge\log5$ (Prop 4.2).
   * *A robust form (Proved).* On the subspace orthogonal to $\mathbf1$ and to the residue indicators, the weighted
     Gram matrix $\sum_p(\log p/e_p)\,a_pa_p^T$ satisfies
     $c^TGc\ge(G_S-\log C)\,\|c\|_1^2$, where $G_S=\frac12\log2+\sum_S\log\ell/(\ell+1)$ (Cor 4.3).
   * *No $c>3$ for bounded exponents.* No version for sub-clusters or bounded exponents gives $W\ge cM\log M$ with
     $c>3$: Theorem 1.1's extremal data already have exponent cap one. **Gap/negative.**
   * *Computed* on 732,612 genuine windows (split-only norms, all $N\le10^{11}$ with at least 64 representations,
     plus random norms up to $10^{20}$, with $C\le16$):
     - no exponent dependence occurs at $C\le\sqrt2$ (theorem);
     - the smallest $C$ with a dependence is $3.34$;
     - the deficiency $M-1-\mathrm{affrank}$ is at most 4 for $C\le16$;
     - at most 3 inert residue blocks are ever needed, against the theorem's 54;
     - in every deficiency-one window the cross-pair average of $\log|b|$ is exactly $\frac12\log2$, and in the
       twelve smallest-$C$ ones $\mathbb E d=W/2$ exactly. So the barycentre accounting is tight there, and all
       the room is in the arc length.
5. **(3) No integer-only refinement of the barycentre argument separates the fakes (Proved, conditional on the fakes
   being as claimed).**
   * Any refinement that uses pair collisions only (at inert primes, at split $p\nmid N'$, at level-wise data for
     small $p\mid N'$) is a consequence of pair instances of (R). Every $A_1(Q_0)$ model therefore satisfies it.
   * The round-5 residue designs and Astra's VIII.3 models are such models, and they have no prime below $2m$
     dividing $N'$, so level-wise data are vacuous for them.
   * A separating input must use genuine Gaussian prime arguments or prime-level residues, the archimedean order, or
     non-pair data that the fakes are not built to satisfy.
   * The record's Paley profiles satisfy VIII.1 with a factor-2 margin: $\omega=q-1\approx2m$, and the affine rank
     is exactly $m-1$. This was Computed for $q=43,127,383$, and Proved in §5.
6. **(4) Cleanest corollary (Proved).** $M\le1+K_C+\omega_1(N')\le1+K_C+\omega_1(\mathrm{rad}\,N)$, where $\omega_1$
   counts prime factors $\equiv1\pmod4$. Consequently:
   * with Robin's bound (quoted), $M\le1+K_C+1.3841\,\log\mathrm{rad}(N)/\log\log\mathrm{rad}(N)$, with **no**
     threshold $R_0$;
   * for $R=n^j$: $M\le(\frac1j+o(1))\log R/\log\log R$;
   * in general $M\le(\frac2\kappa+o(1))\log R/\log\log R$ with $\kappa=\log N/\log\mathrm{rad}(N)$.
   This beats Theorem 1.1 exactly when $\kappa>3$, so for $R$ a perfect square ($1/2$ against $2/3$). It is
   weaker than Theorem 1.1 for $R\in\mathbb N$ ($\kappa\ge2$ only, giving 1) and for squarefree-type $N$.
7. **(5) Assessment (Heuristic, judgement).**
   * VIII.1 is a correct, clean, new structural theorem and should be recorded prominently. It replaces
     $\omega(N')>M/4$ (`round8/stability.md` Thm 3.3) and $M\le(2+o(1))\pi(Y;4,1)$ (Prop 3.4) by
     $M\le1+K_C+\omega(N')$ and $M\le\pi(Y;4,1)+1+K_C$.
   * It does **not** change the growth rate, or the constant $2/3$, for general $R$. It cannot break the fake
     barrier, because it is an $A_1$ statement.
   * The three most promising extensions are ranked in §8.

## 1. Definitions

The cluster consists of $M\ge2$ distinct lattice points on an arc of length $\le C\sqrt R$ of $|z|^2=N=R^2$. It is
normalised as in `two_thirds/proof.md` Lemma 1.1:
* $N'=\prod_pp^{e_p}$ is odd and split-only, and $W=\log N'$;
* $z_i=\varepsilon_i\prod_p\pi_p^{a_{ip}}\bar\pi_p^{e_p-a_{ip}}$;
* $x_{ip}=a_{ip}/e_p\in[0,1]$ and $w_p=e_p\log p$;
* $g_{ij}=\gcd(z_i,z_j)$, $b_{ij}=(z_i-z_j)/g_{ij}$ (the reduced chord, $c_{ij}$ in the record), and
  $d_{ij}=\sum_pw_p|x_{ip}-x_{jp}|$.

The **pair inequality** is
$$d_{ij}\ \ge\ \tfrac W2+2\log|b_{ij}|-2\log C. \tag{1.1}$$
It holds because $\log|g_{ij}|=(W-d_{ij})/2$ and $|z_i-z_j|\le C\sqrt{R'}=Ce^{W/4}$; I re-derived it.

$S$ is a finite set of inert primes, and $\mathrm{ind}_\ell(i)\in\{0,1\}^{\ell+1}$ is the indicator of the class of
$z_i\bmod\ell$ in the norm-$N'$ fibre, which has $\ell+1$ elements.

$G_S:=\frac12\log2+\sum_{\ell\in S}\frac{\log\ell}{\ell+1}$.

$\omega_1(n)$ is the number of primes $\equiv1\pmod4$ dividing $n$.

## 2. VIII.1 re-derived, the exact identity, and a sharper $K_C$

**Lemma 2.1 (the unit-square inequality; Proved).** For $x,y\in[0,1]$:
$$|x-y|=x+y-2xy-2\min(x,y)(1-\max(x,y)).$$
*Proof.* Take $x\le y$. Then $x+y-2xy-(y-x)=2x(1-y)$. $\square$

**Theorem 2.2 (exact barycentre inequality; Proved).** Let $P,Q$ be probability measures on the marks with disjoint
supports, and let $i\sim P$ and $j\sim Q$ be independent. Put $\mu^P_p=\mathbb E_Px_{ip}$ and similarly $\mu^Q_p$.
Put $s_p=(\mu^P_p+\mu^Q_p)/2$, $t_p=(\mu^P_p-\mu^Q_p)/2$, and $m_p=\mathbb E\min(x_{ip},x_{jp})(1-\max(x_{ip},x_{jp}))$.
Then
$$\mathbb E d_{ij}-\tfrac W2=\sum_pw_p\big[2t_p^2-2(s_p-\tfrac12)^2-2m_p\big],\tag{2.1}$$
and therefore
$$2\,\mathbb E\log|b_{ij}|-2\log C\ \le\ \sum_pw_p\big[2t_p^2-2(s_p-\tfrac12)^2-2m_p\big].\tag{2.2}$$

*Proof.* Take expectations in Lemma 2.1, using independence: $\mathbb E[xy]=\mu^P\mu^Q$. Then
$\mu^P+\mu^Q-2\mu^P\mu^Q=2s-2(s^2-t^2)=\frac12-2(s-\frac12)^2+2t^2$. This gives (2.1). Average (1.1) over $P\otimes Q$;
this is legitimate because $i\ne j$ almost surely, by disjointness. $\square$

**Theorem 2.3 (VIII.1, sharpened; Proved).** Let $S$ satisfy $G_S>\log C$. Then the $M$ rows $(a_i,\mathrm{ind}_\ell(i))_{\ell\in S}$
are affinely independent. Hence $M\le1+\omega(N')+K_C$, with $K_C=\sum_{\ell\in S}\ell$.

*Proof.* An affine dependence $c\ne0$ ($\sum c=0$) splits into $P=c^+/|c^+|$ and $Q=c^-/|c^-|$, with disjoint supports
and equal barycentres. So $t_p=0$ for every $p$, and the residue laws $P_\ell=Q_\ell$ agree. The right side of (2.2) is
then $\le0$.

On the left side:
* $(1+i)\mid b_{ij}$ for every pair (record Lemma 3.3(R)).
* $\Pr(\ell\mid b_{ij})\ge\sum_xP_\ell(x)^2\ge1/(\ell+1)$, by Cauchy–Schwarz. A collision mod $\ell$ gives
  $\ell\mid b_{ij}$ because $\ell\nmid N'\supseteq|g_{ij}|^2$.
* Divisors at distinct primes multiply (record Lemma 3.1).

So the left side is $\ge2G_S-2\log C>0$, a contradiction. An indicator block lies in the hyperplane "sum $=1$", so it
adds at most $\ell$ to the affine dimension. $\square$

Astra's statement is the same with $\frac12\log2$ omitted. Every other step agrees with hers.

*Refinements (Proved).*
* **2-adic block.** Odd Gaussian integers fall into 2 classes mod $2$. Equal laws mod 2 give
  $\Pr(2\mid b_{ij})\ge\frac12$, which adds $\frac14\log2$ for affine cost 1. So $K_C\le1$ for $C<2^{3/4}=1.682$.
  With the $\ell=3$ block as well, $K_C\le4$ for $C<2.213$.
* **Optimal $S$.** Among inert-only sets, the optimum is a prefix of $3,7,11,19,\dots$ (verified by brute force for
  $C\le5$ in `kc.py`). Computed values of the sharpened $K_C$: $K_{1.5}=3$, $K_2=10$, $K_3=40$, $K_4=94$,
  $K_5=184$, $K_8=1150$, $K_{10}=2411$, $K_{20}=31014$. Astra's values are 10, 21, 94, 310, 753, 3885, 8684, 111269.
  Since $\sum_{\ell\le x}\log\ell/(\ell+1)=\frac12\log x+O(1)$, we get $K_C=C^{4+o(1)}$.
* **Scope.** The proof uses only (1.1) with $\log|b_{ij}|\ge\log(\sqrt2D_1(i,j))$ and the fibre counts. In angle
  models this is exactly the pair instance of round 5's axiom (R) (`round5/construct.md` §1.2, step 6 of Lemma 1.2:
  $|\sin|\ge D_1e^{-w/2}/\sqrt2$). **So Theorem 2.3 holds for every arc-admissible model in $A_1(Q_0)$ with
  $Q_0\ge\max S$** (Proved).

**Proposition 2.4 (Gram form; Proved).** Put $G=\sum_{p}\sum_{\tau=1}^{e_p}\log p\,\sigma_{p\tau}\sigma_{p\tau}^T$, where
$\sigma_{p\tau}(i)=\pm1$ according as $a_{ip}<\tau$ or not. Then:
* $G\succeq0$ and $G_{ii}=W$;
* $G_{ij}=W-2d_{ij}\le-4\log|b_{ij}|+4\log C=:-h_{ij}$;
* $\mathbf 1^TG\mathbf 1=\sum_p\log p\sum_\tau(2S_{p\tau}-M)^2=4Q$.

Theorem 1.1's identity (2.3)–(2.5) is the statement $\mathbf 1^TG\mathbf1\ge0$. If $C<\sqrt2$ then every $h_{ij}>0$.
By Perron–Frobenius, the matrix $WI-G$ (entrywise $>0$ off the diagonal) has a simple top eigenvalue with a positive
eigenvector, so $\ker G\cap\mathbf 1^\perp=0$. Hence the **cut** rows are affinely independent: $M\le1+\#\{(p,\tau)\}$.

This is weaker than Theorem 2.3 when exponents exceed 1. The reason is that $c\perp$ all cut vectors forces
$c\cdot a_p=\sum_\tau c\cdot\mathbf 1[a_p\ge\tau]=0$, but not conversely.

So Theorem 1.1 and VIII.1 are two certificates for the same PSD system: the all-ones vector, and the kernel. For
$x\ge0$ the system gives $W\ge x^THx/\|x\|^2$. Its optimum, the Perron root of $H$, exceeds the average row sum used
by Theorem 1.1 only through irregular collision patterns, and that excess is $O(\log M)$ for equidistributed classes
(Heuristic computation: the row-sum variance is $\asymp M\log^2M$).

## 3. Question (1): VIII.1 inserted into Theorem 1.1

**Proposition 3.1 (graded Theorem 1.1; Proved, novelty minimal).** For every normalised cluster with $M\ge2$,
$$W'+(K_0+4\log^+C)M\ \ge\ 3M\log M+2M\sum_{p\le M,\ p\nmid N'}\frac{\log p}{p-1}.$$
*Proof.* In (6.2) of `two_thirds/proof.md`, the split term for $p\nmid N'$ may be taken from the first case of
Lemma 4.3, $(\lambda_pM^2/2-M/2)\log p$, instead of $(c_pM^2/4-M/2)\log p$. After multiplying by $8/M$ the extra term
is $(4\lambda_p-2c_p)M\log p\ge2\lambda_pM\log p$, since $c_p\le\lambda_p$ (Remark 4.2.3). The rest is the proof of
Theorem 7.5 verbatim. $\square$

If no split prime $\le M$ divides $N'$, this gives $W'\ge4M\log M-O_C(M)$, which is `growth/auxiliary.md` l.74. So
**the only obstruction between 3 and 4 is split primes $p\le M$ dividing $N'$ with near-balanced levels.** Those
primes cost only $\theta(M;4,1)\approx M/2$ in $W$.

**Proposition 3.2 (no additivity; Proved).** Let $M\ge M_0$ be even and $C\le1$. There are data with the following
properties:
* exponent levels at all split primes;
* residue class counts at all $q\le M$ meeting the pigeonhole minima of Lemma 3.4;
* $W=3M\log M+O(M)$, satisfying inequality (5.1) of `two_thirds/proof.md`;
* in addition, $\omega(N')\ge M$ and affinely independent exponent rows, so that VIII.1's conclusion holds even with
  $K_C=0$.

Hence no argument whose inputs are (5.1) together with the conclusion of VIII.1 proves $W\ge(3+\eta)M\log M$ for any
$\eta>0$.

*Proof.*
1. At every split $p\le M$ take $e_p=1$ and $n_{p0}=n_{p1}=M/2$. At every $q\le M$ take equidistributed classes. By
   Remark (iv) of the record's review status this attains $F_1(M)\le\frac38M^2\log M+O(M^2)$, so the left side of
   (5.1) is $\le F_1(M)$. These small primes contribute $\theta(M;4,1)\le M\log4$ to $W$.
2. Add $M$ split primes $p>M$ with $e_p=1$ and balanced 0/1 columns forming a basis of $\mathbb R^M$. Such columns
   exist: differences of balanced vectors give every $e_i-e_j$, and a balanced vector has nonzero sum, so the
   balanced vectors span $\mathbb R^M$. Linear independence of the columns makes the rows linearly, hence affinely,
   independent.
3. At these new primes $E(n_{pt},(p-1)p^{a-1})=0$, since $n\le M\le p-1$, and $(S_{p1}-M/2)^2=0$. So the left side of
   (5.1) does not change.
4. Choose the new primes near $M^3$, adjusting one of them, so that $W=8F_1(M)/M+O(\log M)$. Then (5.1) holds, and
   $W=3M\log M+O(M)$. $\square$

*The cut identity.* With $\ge M-1-K_C$ primes carrying nontrivial cuts, the cut sum
$\sum_p\sum_\tau(2S_{p\tau}-M)^2\log p$ gains **nothing** for even $M$: balanced cuts contribute 0, as in step 3 of the
proof. For odd $M$ each cut contributes $\ge\log p$, so the sum is $\ge W'$, i.e. $Q\ge W'/4$. The record already uses
this as Remark 5.3.1, changing $M/4$ to $(M-1)/4$. The prime count of VIII.1 adds nothing beyond it.

*Why the gain is zero.* Theorem 1.1's arithmetic lower bound comes from collisions at $q\le M$, and $W$ enters only
as the archimedean budget $MW/8$. VIII.1 constrains how many primes carry $W$, not which primes or how much. Primes
above $M$ with balanced cuts satisfy VIII.1 and are invisible to (5.1).

*Heuristic (pair class).*
* The full family of per-pair inequalities (1.1), as opposed to their average, is the real content available at
  this level. The truth for $A_1$ lies in $[3,\,23/3]$ (if VIII.3 is right) or $[3,8]$ (round 5).
* Equidistant metrics (Hadamard or Paley rows) need margin $W/(2(M-1))\ge2\max_{ij}\log(\sqrt2D_1(i,j))$. Pigeonhole
  on CRT classes mod $\prod_{q\le Y}q$ forces $\max_{ij}\log D_1\ge(1-o(1))\log M$, even when small split primes
  divide $N'$. Astra's "geometric" residues $s_0\omega_i/\bar\omega_i$ reach $\max D_1\le L_0=M^{1+o(1)}$.
* So pair-only $A_1$ models at about $4M\log M$ are plausible. Whether the per-pair constraints can force more than 3
  when small primes divide $N'$ is the open LP described in §8.

## 4. Question (2): strengthening VIII.1

**Proposition 4.1 (where the slack is; Proved).** By (2.1), with $t\equiv0$, the barycentre step loses
$\sum_pw_p[2(s_p-\frac12)^2+2m_p]$. This is zero if and only if, for every $p$:
* the common cross mean is exactly $\frac12$; and
* almost surely one of the two points lies at an extreme level ($0$ or $e_p$).

For squarefree profiles $m_p=0$, so only an imbalance of the mean helps. Middle levels help only when both points sit
strictly inside $(0,e_p)$.

*Computed (`slack.py`, 1,056 genuine windows with exponent deficiency exactly 1):*
* the cross-pair average of $\log|b_{ij}|$ is exactly $\frac12\log2$ in all of them (every cross chord has
  $b=$ unit$\cdot(1+i)$);
* in the twelve with smallest $C$, $\mathbb Ed-W/2=0$ to machine precision;
* in every case $(2\mathbb E\log|b|-2\log C)-(\mathbb Ed-W/2)\le-0.46$.

So both sides of the accounting are tight there, and the room lies entirely in $\log C$: the first dependence appears
at $C=3.34$, far above the threshold $\sqrt2$.

**Proposition 4.2 (weighted rank is trivial; Proved).** Let $c$ be an affine dependence of the rows restricted to
primes $p\le Y$, with residues at $S$. Then (2.2) gives
$2(G_S-\log C)\le\sum_{p>Y}w_p\,2t_p^2\le W_{>Y}/2$. So such dependences force only $W_{>Y}\ge4(G_S-\log C)$, whereas
the rank count already gives $\mathrm{rank}_{\le Y}\ge M-1-K_C-\omega_{>Y}$ with $\omega_{>Y}\le W_{>Y}/\log Y$. No
weighted improvement is possible this way, because every prime weighs at least $\log5$.

**Corollary 4.3 (robust form; Proved).** Let $U=\{c:\ c\perp\mathbf1,\ c\perp\mathrm{ind}_\ell(\cdot)\ \forall\ell\in S\}$
and $G'=\sum_p(\log p/e_p)\,a_pa_p^T$. Then for every $c\in U$,
$$c^TG'c\ \ge\ (G_S-\log C)\,\|c\|_1^2.$$
*Proof.* Put $P=c^+/|c^+|$, $Q=c^-/|c^-|$, with $|c^\pm|=\|c\|_1/2$. Then $t_p=c\cdot a_p/(e_p\|c\|_1)$ and (2.2) gives
$\sum_pw_p2t_p^2\ge2(G_S-\log C)$, which is the displayed inequality. $\square$

There is no approximate dependence with small coefficients. But since $\mathrm{tr}(G'|_U)\le MW/4$, this lower bound
on the bottom of the spectrum is far below the average and gives no growth information.

**Sub-clusters and bounded exponents (Gap; negative).**
* VIII.1 for a subset $T$ is implied by VIII.1 for all marks, since subsets of affinely independent sets are
  independent.
* Theorem 1.1 applied to sub-clusters is the record's Prop 2.2 recursion (`round8/stability.md`). Uniform pair
  weights are optimal for the averaged method: restricting to a cell of size $n$ raises the diagonal penalty to
  $W/(2(n-1))$.
* Theorem 1.1's extremal data (Prop 3.2) have all exponents equal to 1. So bounded exponents cannot by themselves
  push the constant past 3 in any argument whose inputs are (5.1) and VIII.1.
* With exponents $\le k$, VIII.1 and Theorem 1.1 give $P_+(N')\ge\max(2,6/k)\,M\log M\,(1-o(1))$. For $k\le2$ this is
  Theorem 1.1's bound, not VIII.4's.

*Conjectures tested (Computed, `analyse.py`, 732,612 windows; full run output in §9).*
* (T1) Sharpened Theorem 2.3: 0 violations.
* (T2) "No exponent dependence for $C\le C_1$": true in the data for $C_1=3.34$. Its first instance is $M=4$,
  $\omega'=7$, $N=9862667358769$. The theorem proves $C_1\ge\sqrt2$.
* (T3) Deficiency against $C$: at most 1 on $(3,4]$, at most 2 on $(4,8]$, at most 4 on $(8,16]$. The extreme case is
  the record's 10-point cluster: $N=1176852625$, $C=8.2393$, $\omega'=5$, deficiency 4, made good by $\{3,7\}$ as
  Astra reports. At most 3 inert blocks are ever needed, against the theorem's 18–54.
* The data are consistent with a deficiency of $O(\log C)$, but they are small ($M\le10$) and prove nothing. **Heuristic.**

## 5. Question (3): can a refinement be violated by the fakes?

**Proposition 5.1 (no-go for pair refinements; Proved, conditional on the cited constructions).** Let $\mathcal I$ be
any inequality derived from the arc condition, the pair instances of (R) at odd prime powers $\le Q_0$ (including
level-wise data at split $p\le Q_0$ dividing $N'$), and the conic fibre counts. Then every arc-admissible model of
$A_1(Q_0)$ satisfies $\mathcal I$.

The VIII.3 models are claimed to satisfy (R) at every odd prime power $\le2m$, with all primes of $N'$ near
$q^{1.99}\gg2m$. The round-5 designs have all primes above the modulus bound. In both, the level-wise data are vacuous.
So for $Q_0\le2m$, **no refinement of the barycentre argument through collisions at inert primes, or at split primes
not dividing $N'$, separates these fakes.**

*Proof.* Each such input holds in $A_1(Q_0)$ by definition (`round5/construct.md` Lemma 1.2). Averaging, choosing
measures and pigeonholing preserve validity. $\square$

The correctness of the VIII.3 construction itself is outside this block.

*Paley profiles (Proved and Computed).*
* $a_{ij}=1+\chi(j-i)$ on rows $S\subset\mathcal A$, with $|S|<q$. A dependence would give $K^Tc=0$, so
  $c\in\ker K=\mathrm{span}\,\mathbf1$, which is impossible on a proper subset. Hence the affine rank is $|S|-1$.
* Computed for $q=43,127,383$ (`paley.py`): rank $m-1$, $\omega=q-1\approx2m$, every pair at $\ell^1$ distance
  $q+1$, margin $2\log P$.
* So VIII.1 holds with room to spare. The factor 2 cannot be attacked by a VIII.1-type count, because genuine clusters
  reach $\omega'=M/2$ (the 10-point record).

*Heuristic.* Astra's residues are geometric: $D_1\mid\mathrm{Im}(\omega_i\bar\omega_j)$ and
$D_u\le\sqrt2L_0^{s/2}$, independently of the modulus bound. So her models plausibly satisfy (R) at **all** odd moduli
coprime to $N'\prod\ell_i$, and raising $Q_0$ would not help either. I did not check this.

*What a separating input must use.* It must use one of the following:
* the genuine association $p\mapsto\arg\pi_p$, or prime-level residues (the VIII.4 area test is of this kind);
* the archimedean order beyond the chord bound (triple and curvature inequalities);
* data at primes dividing $N'$ above $2m$.

## 6. Question (4): corollaries for the original conjecture

**Corollary 6.1 (Proved).** For every $C>0$ and every arc of length $\le C\sqrt R$:
$$M\ \le\ 1+K_C+\omega_1(N')\ \le\ 1+K_C+\omega_1(\mathrm{rad}\,N).$$
With Robin's inequality $\omega(n)\le1.3841\log n/\log\log n$ for $n\ge3$ (quoted, not re-checked):
$$M\le1+K_C+1.3841\,\frac{\log\mathrm{rad}(N)}{\log\log\mathrm{rad}(N)}$$
for **all** $R$, with no $R_0$ and with explicit $K_C$ (§2).

Since the product of the first $r$ primes $\equiv1\pmod4$ is $e^{(1+o(1))r\log r}$, we get
$\omega_1(m)\le(1+o(1))\log m/\log\log m$. Write $\kappa=\log N/\log\mathrm{rad}(N)\ge1$. Then:
* $M\le(\frac2\kappa+o_C(1))\log R/\log\log R$ as $\mathrm{rad}(N)\to\infty$; if the radical stays bounded, so does $M$.
* If every prime of $N$ has exponent $\ge k$, then $\kappa\ge k$. This is VIII.5, and I re-derived its proof:
  $\log N\ge k\log((r+1)!)$ with $r\ge M-1-K_C$.
* $R=n^j$ ($j\ge1$): $\mathrm{rad}(N)\mid n$, so $M\le1+K_C+\omega_1(n)\le(\frac1j+o(1))\log R/\log\log R$.
  For $R$ a perfect square this gives $\frac12$ against Theorem 1.1's $\frac23$.

*Check against Theorem 1.1.* The two bounds combine as a minimum: $\min(\frac23,\frac2\kappa)$. VIII.1 wins exactly
when $\kappa>3$, so for every $R=n^j$ with $j\ge2$. It loses for integer $R$ ($\kappa\ge2$ gives only 1), and for
$N$ near squarefree.

As in §3, the two do not interpolate. Data with all $e'_p=3$ at $M$ primes of size about $M$, and
$v_p(N)=\kappa\ge3$, satisfy both with $\log N=\max(3,\kappa)M\log M$. So nothing better than the minimum follows from
these inputs (Proved for the inputs; Heuristic as to realisability).

## 7. Question (5): honest assessment

* **New and worth recording (Heuristic judgement; the mathematics is Proved).** VIII.1 is short and correct, holds
  for every $C$, and is sharp in its accounting: the deficiency-one genuine windows sit exactly at
  $\mathbb E\log|b|=\frac12\log2$ (§4).
  * It upgrades the record's prime-count statements: $\omega(N')>M/4$ for $C\le1$, $\ge(M/2)(1-M^{-1/3})^2$, and
    $M\le(2+o(1))\pi(Y;4,1)$ all become $\omega(N')\ge M-1-K_C$, uniformly in $C$. VIII.4 ($P_+\ge(2-o(1))M\log M$)
    follows.
  * Literature: a web search (snippets only, not read) found only Cilleruelo–Granville's four-point classification
    and the Cilleruelo–Córdoba and Bourgain–Rudnick bounds. I found no prime-count statement of this type there.
    Novelty in the literature is therefore unverified.
* **Not a growth result.** For general $R$ it gives $W\ge(1+o(1))M\log M$ only. It does not interact with the
  collision side of Theorem 1.1 (Prop 3.2). It is an $A_1$ theorem, so every fake family of the record satisfies it.
  The growth rate and $2/3$ are unchanged. The only genuine improvement is in the class $\kappa>3$, which includes
  perfect-square radii.
* It should be recorded with the sharpened $K_C$ (the $\frac12\log2$ term and the 2-adic block), the $A_1$ scope
  remark, and the $\kappa$-corollary (Cor 6.1).

## 8. The three most promising extensions, ranked

1. **The per-pair (cut-cone) problem for $A_1$: does any pair input beat 3?**
   * Theorem 1.1 and VIII.1 are two certificates of one system: $G$ in the cut cone, $G_{ii}=W$, $G_{ij}\le-h_{ij}$
     (Prop 2.4).
   * The constant moves from 3 only through small split primes dividing $N'$ (Prop 3.1).
   * *First step:* solve the exact LP over cut metrics (column generation) for $M\le40$, with adversarial levels at
     $p\le M$ and residue classes. In parallel, build pair-only $A_1$ models with balanced small primes dividing $N'$
     and geometric residues.
   * If $W=(3+o(1))M\log M$ is reached, 2/3 is the limit of all pair inputs, which is a clean barrier. If not, the
     per-pair constraints improve 2/3.
2. **Structure of circuits in genuine clusters; $K_C$ polynomially smaller.**
   * The theorem has $K_C=C^{4+o(1)}$. The data have deficiency at most 4 for $C\le16$, and every cross chord of a
     minimal dependence is $b=$ unit$\cdot(1+i)$.
   * *First step:* prove that in a circuit $(P,Q)$ with $\mathbb E\log|b|<\log C$, most cross pairs have
     $|b_{ij}|\le C$, i.e. $z_i/g_{ij}-z_j/g_{ij}$ is a short vector. Then count such configurations by the number of
     Gaussian integers of norm $\le C^2$, aiming at a deficiency of $O(C^2)$ or $O(\log C)$.
3. **Interpolating the two bounds on non-squarefree norms ($2<\kappa<3$).**
   * *First step:* decide whether Prop 3.2's non-additivity survives when the exponents are forced. Put
     $\kappa=W/\log\mathrm{rad}(N')$, so that the large primes must carry exponent about $\kappa$. Test whether the
     chain lemma at the large primes ($e_p\ge2$, so middle levels have slack $m_p>0$ in (2.1)) together with VIII.1
     gives $W\ge c(\kappa)M\log M$ with $c(\kappa)>\max(3,\kappa)$.

## 9. Files and reproducibility

All files are in `round31/ref8_extend/`. Run with `prlimit --as=3500000000 python3 ...`.

* `kc.py`: Astra's and the sharpened $K_C$, and a brute-force check that the optimal sets are prefixes.
* `gaussx.py`: exact Gaussian arithmetic and exact rational rank.
* `enum_clusters.py`: random split-only norms up to $10^{20}$, and windows with $C\le16$.
* `dfs_norms.py` and `scan.py`: all split-only $N\le10^{11}$ with at least 64 representations and primes $\le200$
  (110,295 norms).
* `analyse.py`: affine deficiencies, residue blocks needed, and the tests T1–T3. Ranks are computed mod
  $2^{61}-1$, with exact rational recomputation whenever the result is not full.
* `slack.py`: barycentre accounting on deficiency-one windows.
* `paley.py`: the Paley-profile checks.

Key output of `analyse.py`:
```text
windows analysed 732612
violations of sharpened VIII.1 (need > |S_C|): 0
windows with C<=sqrt2 and d0>0: 0
def0>=1: min C=3.3406 (M=4, omega'=7, N=9862667358769)
def0>=2: min C=4.0070 (M=6, omega'=6, N=1720010305)
def0>=3: min C=7.5214 (M=8, omega'=5, N=1176852625)
def0>=4: min C=8.2393 (M=10, omega'=5, N=1176852625)
```

Literature: one web search; I saw snippets only, of the Cilleruelo–Granville "Close lattice points on circles" page
(Canad. J. Math.) and of arXiv math/0608109. I did not read either paper.
