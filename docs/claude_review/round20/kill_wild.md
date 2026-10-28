# Adversarial kill report: `round20/wild.md` (unlikely-intersections lens)

Scripts and outputs: `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round20/kill_wild_checks/`
(`gauss.py` exact Gaussian-integer toolkit; `smin_counter.py/.out`, `smin_alln.py/.out`, `smin_explore.py/.out`,
`recheck.py/.out`; all exact integer arithmetic, logs only at the end, run under `prlimit --as=3500000000`).

## 0. Verdict

**Killed.** The proposal already declares itself "reduced to a cap". The re-derivation confirms that cap. The one
statement it leaves open, (S$_{\min}$), is **false**, and false for **every** $n$, so nothing of the lens survives:

* (S$_{\min}$) fails at $n=2$ on the proposal's own example $(1:\gamma:\gamma^2)$, $\gamma=A+i$, $A\to\infty$. The write-up
  says the examples of Prop 6.6 "do not refute it", but that check only let $K\to\infty$; letting $A\to\infty$ refutes it.
* (S$_{\min}$) also fails at $n=2$ on **genuine short-arc circle clusters**: $z=((b-i)^2,\,b^2+1,\,(b+i)^2)$ with $b$ even, which
  satisfy $C\to4$, rank 1 and $R\to\infty$.
* (S$_{\min}$) fails for **every** $n\ge2$ (Prop K.2 below; a proof). Take $u=(1,X_1+i,\dots,X_{n-1}+i,(X_1+i)(X_2+i))$ with
  $X_x=T+x$. Its min-anchor ratio is $\to\frac n2>\frac n3$, its rank deficiency is 1, and there is one multiplicative
  relation. A variant with rank $\le k\approx\sqrt{2n}$ gets ratio $\ge(0.64\text{ to }0.81)n$ (exact evidence). A single
  circuit (an even cycle of products $(T+v+i)(T+v+1+i)$) gets ratio $\to n$ (exact evidence). So no version restricted
  to large $n$, to large rank deficiency, or to relations of full support can be rescued.

Every proved step of §6 (Lemma 6.1, Thm 6.2, Lemma 6.3, Prop 6.4, Prop 6.6, Lemma 6.7, Prop 6.8) re-derives correctly,
and every number in §7 recomputes, apart from one rounding typo. The barrier-map content of the proposal (fixed carriers
vacuous, uniform line comparisons false, ceiling $\omega(N')$) stands, with the corrections in §2 below. In particular the
"ceiling" Cor 6.9 is proved only for subgroup membership ($\varepsilon=0$), not for Rémond cones.

## 1. What the proposal claims

The lens is unlikely intersections in tori: BMZ, Habegger's bounded height theorem, Rémond cones, and anomalous
Zilber–Pink. The write-up proves four things, all unrefereed:
(1) The fixed carriers $X_M$, $L_M$ ($M\ge4$) and $M_{0,k}$ ($k\ge5$) have $V^{oa}=\varnothing$ (Lemma 6.1, Thm 6.2).
(2) The moving carrier $\ell_u$ is never in a coset (Lemma 6.3). A cluster has $h(z)=\log R$,
$h(\ell_z)\le\frac32\log R+\log C$, and chart sums $\ge(M-1)(\frac12\log R-\log C)$ (Prop 6.4).
(3) The uniform-in-the-line comparisons (UI-$\ell$) and (UI-$\ell_\Sigma$) are false (Prop 6.6).
(4) An unlikely-intersection input acts only on rank-deficient clusters, so its ceiling is $M\le\omega(N')+O_C(1)$
(Lemma 6.7, Prop 6.8, Cor 6.9).
The only residue left open is (S$_{\min}$): for some $\delta>0$, every multiplicatively dependent
$u\in\mathbb P^n(\mathbb Q(i))$ whose coordinates are pairwise non-$\mathbb R$-proportional satisfies
$\min_a\sum_{x\ne a}h(u_x/u_a)\le(\frac13-\delta)\,n\,h(\ell_u)+e(n)$ with $e(n)=o(n^2\log n)$.

## 2. Re-derivation of the proved steps

Every item was re-derived by hand. Items marked OK are correct as written.

* **Lemma 6.1: OK.** For a dominant morphism of irreducible varieties $V\to V'$, every component of every nonempty fibre
  has dimension $\ge\dim V-\dim V'$. This holds at all points of the image, not only generically, so the "every point"
  claim is right. A fibre of $\varphi$ is a finite union of cosets of $(\ker\varphi)^0$, and the anomaly inequality
  $\dim Y>\max(0,\dim V-m)$ is exactly the two hypotheses.
* **Thm 6.2: OK.**
  * (a) $X_M$ is, as in overview §4.1, the cone over $\{u_1v_1=\dots=u_Mv_M\}$. Its torus part is the subtorus cut out by
    the characters $e_{u_x}+e_{v_x}-e_{u_1}-e_{v_1}$. That lattice is saturated (each character has its own $u_x$), so
    the subtorus is connected. Its dimension is $M$ in $\mathbb G_m^{2M}/\mathbb G_m$, and $M>\max(0,1)$.
  * (b) $D_{12}-D_{13}+D_{23}=0$. The image is a curve in $\mathbb G_m^2$ and $\dim V=M-2\ge2$.
  * (c) $\dim T_4=2$, the image is $M_{0,4}$, and $1<k-3$. The anomaly test gives $\dim Y\ge k-4>\max(0,k-5)$.
* **Lemma 6.3: OK.** The function $\prod L_x^{m_x}$ with $\sum m_x=0$ has divisor $\sum m_x[\zeta_x]\ne0$, because the
  $\zeta_x$ are distinct.
* **Prop 6.4: OK.**
  * (ii) The finite places contribute 0 and the complex place gives $\log\max|z_x|$.
  * (iii) $|\mathrm{Im}(\bar z_az_b)|\le N\cdot CR^{-1/2}$, and $\gcd\ge1$ is the safe direction.
  * (iv) $h(z_x/z_a)=\log|\beta|=\frac12d_{xa}$, and $d\ge\log R-2\log C$ follows from $\mathrm{Norm}\,g\le|z_x-z_a|^2\le C^2R$.
  * **Correction to the reading, not to the proof.** (iii) is only an upper bound. On genuine short-arc clusters the
    Plücker gcd can be of size $\sqrt R$. For the family $z=((b-i)^2,b^2+1,(b+i)^2)$ ($b$ even, $C\to4$) one has
    $G=2b$ and $h(\ell_z)=\log(2(b^2-1))=\log R+O(1)$, not $\frac32\log R$ (`smin_counter.out` §3).
  * So the §7 evidence "Plücker gcd 1 and ratio equal to the bound on actual clusters" is not representative. On real
    clusters $h(z)/h(\ell_z)$ can tend to $1$, not $\frac23$. The phrase "the moving carrier encodes the arc exactly in
    its height" holds in one direction only.
* **Prop 6.6: OK**, including the algebra of (d). If $\alpha+\beta<\frac{2\alpha+n\beta}3$ then $\alpha<(n-3)\beta$ and
  $\beta>0$, which gives $\frac{n\beta-\alpha}6>\frac\beta2$.
  * (b) and (c) use $h(\alpha)=\log|\alpha|$ for $\alpha\in\mathbb Z[i]$, which is correct.
  * Typo: §6.3 quotes "0.9906 at $(10,20,2000)$". The computed ratio is $0.990755$; $0.990589$ is the rigorous lower bound.
* **Lemma 6.7: OK.** After normalisation $N$ is odd, since $(1+i)$ and $(1-i)$ are associates and gcd 1 removes them.
  The $\pi_j,\bar\pi_j$ are then pairwise non-associate, and $\chi_m(z)\in\mu_4$ gives $z\in\ker\chi_{4m}$.
* **Prop 6.8: OK.** The (UI-$\ell_\Sigma$) branch gives $\frac32\delta\log R\le(\frac43-\delta)\log C+e(M-1)/(M-1)$.
  The step $e(n)/n=o(n\log n)$ is what kills the result against Theorem 2.1, and that is correctly stated.
* **Cor 6.9: overclaimed in scope.** The proposal's summary bullet (4) says "even an ideal unlikely-intersection input
  gives at best $M\le\omega(N')+O_C(1)$". That is proved only for inputs that act through subgroup membership
  ($\varepsilon=0$, via Lemma 6.7). For Rémond cones the write-up's own remark is only **heuristic**.
  * Impact: none on the verdict. The cone versions of (UI-$\ell$) and (UI-$\ell_\Sigma$) are refuted anyway, since cones
    enlarge the class, and (S$_{\min}$) is refuted with $\varepsilon=0$.
  * Correct wording: "every conclusion delivered through subgroup membership".
* The quoted forms of BMZ 1999 (Thm 1: bounded height on $\mathcal C\cap\bigcup_{\mathrm{codim}\ge1}H$ for curves not in a
  translate of a proper subtorus) and Habegger 2009 agree with my recollection of the literature; not re-verified.

## 3. Numerical checks (independent code, `recheck.out`)

| claim in wild.md | recomputed | status |
|---|---|---|
| $F_1$ ratio at $(A,n,K)=(10,20,2000)$ "0.9906" | 0.990755 (bound 0.990589) | typo; claim holds |
| $F_1$ at $(2,30,3000)$; $F_2$ chart ratio 14.93 vs $n/2=15$ | 0.990428; 14.9325 | exact match |
| $F_2$ at $(3,8,50)$: 4.3582, $G_2=2^{21}$ | 4.3582, $\log_2G_2=21$ | match |
| $(1:\gamma:\gamma^2)$: $h(u)=h(\ell_u)=\log(A^2+1)$ | exact for $A=10,\dots,10^{24}$ | match |
| 10-point cluster on $1176852625=5^3 13^2\,17\cdot29\cdot113$: chord $C=8.2386$, rank 5 | $8.2386$, rank 5 | match |
| 7-window on $3159797225$: rank $6=M-1$ (full) | $C=10.5641$, rank 6 | match |

## 4. (S$_{\min}$) is false for every $n$

Notation: $\Sigma_{\min}(u)=\min_a\sum_{x\ne a}h(u_x/u_a)$. (S$_{\min}$) asserts
$\Sigma_{\min}(u)\le(\frac13-\delta)n\,h(\ell_u)+e(n)$, where $e$ depends only on $n$. To refute it, it is enough to find,
for one fixed $n$, a sequence of admissible $u$ with $h(\ell_u)\to\infty$ and
$\Sigma_{\min}(u)\ge(\frac13+\eta)n\,h(\ell_u)-O_n(1)$ for some $\eta>0$.

**Prop K.0 ($n=2$, the proposal's own point; proved).** Let $u=(1:\gamma:\gamma^2)$ with $\gamma=A+i$, $A\ge2$. It lies in
$\mathcal U_2(0)$ (Prop 6.6(a)). The Plücker coordinates are $1$, $2A$ and $A^2+1$, so $h(\ell_u)=L:=\log(A^2+1)$. The
anchor $\gamma$ gives $h(\gamma^{-1})+h(\gamma)=L$, and the other anchors give $\frac32L$. So $\Sigma_{\min}=L=h(\ell_u)$, the
ratio is $1>\frac23$, and $L-(\frac23-2\delta)L-e(2)\to\infty$ as $A\to\infty$.
So the sentence in wild.md §6.3 "(S$_{\min}$) is not refuted by $F_1,F_2$" is **wrong**. It looked only at $K\to\infty$.
$F_1$ with $n=K=2$ and $A\to\infty$ refutes it. The same holds at $n=3$: $(1:\gamma:\gamma^2:\gamma^3)$ has
$\Sigma_{\min}=h(\ell)=2L$ exactly, ratio $1=\frac n3>(\frac13-\delta)n$ (`smin_counter.out` §2).

**Prop K.1 (genuine circle clusters; proved).** Let $b\ge2$ be even and $z=((b-i)^2,\ b^2+1,\ (b+i)^2)$.
* Equal norms and gcd: all three have norm $N=(b^2+1)^2$, so $R=b^2+1$. The Gaussian gcd is $1$, because $b\pm i$ are
  coprime: their gcd divides $2i$ and the norm is odd.
* Arc: the arc length tends to $4\sqrt R$ from above ($4.0066$ at $b=10$, $4.000067$ at $b=100$).
* Rank: $z_1z_3=z_2^2$, so the cluster has rank $1$.
* Plücker coordinates: they are $2b(b^2+1)$, $4b(b^2-1)$ and $2b(b^2+1)$. Their gcd is $2b$, because $\gcd(b^2+1,2(b^2-1))=1$
  for even $b$. So $h(\ell_z)=\log(2(b^2-1))$.
* Anchor sums: the anchor $z_2$ gives $\Sigma=2\cdot\frac12\log(b^2+1)=\log R$, and the other anchors give $\frac32\log R$.
* Conclusion: $\Sigma_{\min}/h(\ell_z)\to1>\frac23$, so (S$_{\min}$) fails at $n=2$ on actual short-arc clusters with
  $C\to4$ and $R\to\infty$ (exact values in `smin_counter.out` §3).
* By Prop 6.8's own logic, any statement that bounds $R$ for rank-deficient $M$-clusters at fixed $M=3$ is false.

**Prop K.2 (every $n\ge3$, deficiency 1; proved).** Let $n\ge3$, $T\ge n^2$, $X_x=T+x$, and
$$u=(1,\ X_1+i,\ \dots,\ X_{n-1}+i,\ (X_1+i)(X_2+i)).$$
(a) *Admissible.*
* Dependence: $u_n u_0=u_1u_2$, so the point is multiplicatively dependent, i.e. it lies in $\mathcal U_n(0)$.
* No real ratio between singles: $(X_a+i)/(X_b+i)\in\mathbb R$ iff $X_a=X_b$.
* No real ratio involving $u_0$ or $u_n$: $\mathrm{Im}\,u_n=2T+3\ne0$, and $\mathrm{Im}(\bar u_xu_n)=T^2+2xT+3x-1\ne0$.

(b) *Height of the line.* The Plücker coordinates are:
* $p_{0x}=1$ and $p_{0n}=2T+3$;
* $p_{xy}=X_x-X_y$;
* $p_{xn}=X_x(2T+3)-((T+1)(T+2)-1)=T^2+2xT+3x-1$.

All are $\le(T+n)^2$. Since $\gcd\ge1$ (the safe direction), $h(\ell_u)\le2\log(T+n)$.

(c) *Every pair ratio is large.* Write $u_x/u_a=\beta/\delta$ in lowest terms and $g=\gcd(u_x,u_a)$. Then
$h(u_x/u_a)\ge\log(\max(|u_x|,|u_a|)/|g|)$. Cases:
* Two singles $x,a\in[1,n-1]$: $g$ divides $X_x-X_a\in\mathbb Z\setminus0$, so $|g|\le|x-a|<n$ and $h\ge\log T-\log n$.
* Either coordinate equal to $u_0=1$: $h=\log|u|\ge\log T$.
* $u_n$ against $u_1$ or $u_2$: the ratio is the other single, so $h\ge\log T$.
* $u_n$ against $u_a$ with $a\ge3$: $|g|\le|\gcd(u_1,u_a)\gcd(u_2,u_a)|<n^2$ and $|u_n|\ge T^2$, so $h\ge2\log T-2\log n$.

So $\Sigma_{\min}(u)\ge n(\log T-2\log n)$.

(d) *Contradiction.* (S$_{\min}$) would give $n(\log T-2\log n)\le(\frac13-\delta)\,2n\log(T+n)+e(n)$, that is
$(\frac13+2\delta)\,n\log T\le O(n\log n)+e(n)$. For fixed $n$ this is false once $T$ is large. $\square$

The limiting ratio $\Sigma_{\min}/h(\ell_u)$ is $\frac n2$. Exact values for $n=2,\dots,40$ and $T=10^{10},10^{30},10^{60}$
are in `smin_alln.out`: ratio$/n=0.494$–$0.500$, $G=1$.

**K.3 (restricted versions also fail; exact evidence, with proofs that follow the pattern of K.2).**
* *Small rank.* Take $k$ generic singles $X_x+i$ plus all $\binom k2$ pairwise products, so the rank is $\le k\approx\sqrt{2n}$.
  Ratio$/n$ is $0.643,\ 0.726,\ 0.775,\ 0.809$ for $n=10,21,36,55$ (`smin_alln.out`).
  Reason: product–product Plücker coordinates cancel at order $T^3$, so $h(\ell)=2\log T+O(\log k)$, while disjoint
  products have ratio height $2\log T$.
* *Single circuit (a relation with full support and rank deficiency exactly 1).* Take an even cycle $C_{2k}$ and
  $u_x=(T+x+i)(T+x+1+i)$, indices mod $2k$. Then $\prod_{x\ \rm even}u_x=\prod_{x\ \rm odd}u_x$.
  Ratio$/n$ is $0.66,\ 0.79,\ 0.84,\ 0.89,\ 0.92,\ 0.94$ for $n=3,5,7,11,19,29$ at $T=10^{30}$ (`smin_circuit.out`), and it
  tends to $1$.

So no rescue of (S$_{\min}$) works: not "for $n\ge n_0$", not "rank $\le n-d$", and not "relation of full support".

**Why it fails (structural).** Multiplicative dependence is one relation. It can be carried by a handful of coordinates
while the others stay generic, or be spread around a cycle of nearly parallel Gaussian integers $T+v+i$.
* For generic tuples close to a common real direction, $h(\ell_u)\approx2\log T$ while every pair ratio has height
  $\gtrsim\log T$. So the generic ratio is $\ge\frac n2$, and dependence costs nothing.
* What makes clusters special (ratio exactly $\frac n3$ at the Dirichlet exponent) is the **equal-norm constraint**, which
  (S$_{\min}$) drops on purpose.
* Restoring that constraint gives back the circle problem itself, i.e. Vojta on $X_M$ / $(L_k)$. This is the proposal's own
  §6.3 reading, now confirmed: the residue was never weaker than the target.

## 4b. Attacks (i)–(iii) on the lens as a whole

* **(i) Fakes.** The lens's best conclusion through Lemma 6.7, "rank-deficient clusters are bounded", is satisfied by
  every fake profile. Those profiles have full rank $M-1$ (overview §3.3), so they lie in no proper subgroup, and the
  proposal concedes this.
* **(i) Algebraic fakes.** K.1 adds genuine short-arc algebraic points (actual clusters, $C\to4$) whose height ratios beat
  the cluster bound of Prop 6.4. So $(\ell_z,\ \text{subgroup membership})$ data cannot separate "few points" from "many
  points" even inside the true circle problem.
* **(ii) Reduction.**
  * The ideal subgroup-membership input caps at $M\le\omega(N')+O_C(1)$. That is the lit18 Prop B ceiling of the QS route,
    i.e. checklist item 5 in effect: the only cluster information used is the multiplicative rank of the exponent rows.
  * The cone version adds nothing proved.
  * (UI-$\ell$), (UI-$\ell_\Sigma$) and (S$_{\min}$) are all false.
* **(iii) Strength.** The only form not refuted is the equal-norm comparison on $X_M$ or $M_{0,k}$, which is Vojta-type /
  $(L_k)$-type (target strength). The fixed-carrier unlikely-intersection versions there are vacuous ($V^{oa}=\varnothing$).

## 5. Checklist placement

| item | placement after this report |
|---|---|
| 1 (fake-satisfied input) | The conclusion is fake-satisfied (full-rank fake profiles are invisible). The inputs are statements about algebraic points, but K.1 shows that the line/subgroup data is also "satisfied" by genuine 3-clusters at ratio 1. |
| 2 (averaged) | not applicable |
| 3 (auxiliary polynomials) | not applicable. On $X_M$ and $M_{0,k}$ the lens is vacuous (Thm 6.2, confirmed). |
| 4 (single form / SLW) | not applicable |
| 5 ($S$ through $\lvert S\rvert$) | The constants are $S$-free, but the input is multiplicative rank, which lands at the lit18 Prop B ceiling $M\le\omega(N')$. |
| 6–8 | not used |
| Dirichlet exponent | (UI-$\ell$) needs $\kappa<\frac23$ exactly at $\alpha=\frac12$, and it is false even at $\kappa<1$. |
| Project constraints | none violated: gcds are kept, no phase is set to zero, and the BMZ constant per line is not treated as uniform |

## 6. What survives

* **Nothing toward the target.** There is no growth improvement and no uniform bound. The single open residue
  (S$_{\min}$) is refuted (Props K.0–K.2, proved; K.3, exact evidence), and so are its natural rescues.
* **What survives as barrier facts** (the re-derivations confirm them; still unrefereed beyond this report):
  * Thm 6.2: $V^{oa}=\varnothing$ for $X_M$, for $L_M$ ($M\ge4$) and for $M_{0,k}$ ($k\ge5$) in the Plücker torus.
  * Lemma 6.3 and Prop 6.4, with the correction that (iii) is only an upper bound: genuine clusters can have Plücker gcd
    $\asymp\sqrt R$ and $h(\ell_z)=\log R+O(1)$.
  * Prop 6.6, Lemma 6.7 and Prop 6.8.
  * Cor 6.9, restricted to subgroup-membership inputs (the cone case is heuristic only).
* **Suggested barrier-map row.** "Unlikely intersections (BMZ / Habegger / Rémond / anomalous ZP)":
  * fixed carriers vacuous;
  * uniform moving-line comparisons false at every useful constant, already on actual 3-clusters (K.1) and for every $n$
    with a single relation (K.2);
  * the ceiling through subgroup membership is $\omega(N')$ (= lit18 Prop B).
* **Status: killed.**
