# Round 29: adversarial check of lane spectral

Scripts and outputs: `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round29/kill_spectral_checks/`
(`kernel_indep.py` + `kernel_indep_out_fixed.txt`, `kappa_indep.py` + `kappa_out.txt`, `resonance_indep.py` + `resonance_out.txt`, `small_checks.py` + `small_checks_out.txt`).
Notes: `.../round29/kill_spectral_notes.md`.

## 0. Verdict

**The lane's obstruction survives, with corrections.** Nothing in the lane proves or suggests a growth improvement. No claimed positive result exists to kill. The negatives are mostly right; four of them need correcting.

1. **Prop. 2.0 is overstated.** It is true for bounded-degree statistics. It is false as a blanket claim that every input with a $Q^\varepsilon$ loss is vacuous: a loss of $Q^\varepsilon$ on a degree-$d$ statistic costs only $Q^{\varepsilon/d}$ in $M$.
2. **Lemma 2.3 is over-parametrised.** The unit $u$ is redundant. The restriction "only $u=1$ for non-trivial $I$" is a normalisation, not a consequence of non-triviality.
3. **The model constant $\kappa\in[0.25,0.33]$ is a finite-size effective value.** The heuristic leading constant is $2/\pi^2\approx0.203$, and the data's slope in $\log R$ matches it. The order $\log R$ stands.
4. **The kernel numerics are thinner than stated.** The $\sigma=0.05$ rows lie outside the Hardy strip. The "cancellation to $10^{-2}$" occurs only at the largest $N$.

There is also a sharpening that the lane missed (§2.4). The "minor-arc lattice counts" in the closing Lead are **not** an $\mathrm{SL}_2(\mathbb Z)\backslash\mathrm{SL}_2(\mathbb R)$ equidistribution problem. Every sublattice $\bar P'\mathbb Z[i]$ is similar to $\mathbb Z[i]$, so $E^W$ is exactly a sum of short-sector counts on the circles of radius $R/|P'|$. The Lead's hypothesis is therefore self-similar: it asks for the original problem at smaller scales, summed. This makes "restatement" a proved structural fact rather than a judgement.

**Novelty:** a repackaged known barrier. The transplant negatives are new packaging of record barriers: fake-blindness, criticality, Jarník lines.

## 1. The account of Li–Liu checked against the PDF

I read the PDF in place (pdftotext to stdout, no copy). Every load-bearing citation in lane §1 is **confirmed verbatim**:

* Lemma 3.1, $G=D-P$ with $0\le P\le D\ll_\phi\log N$. So $|G-c_\phi\mathfrak S|\ll\log N$ trivially, since $\mathfrak S(N)\ll\log\log N$.
* (2.28): $S(N;h)\iff\|T_N\|\le\mathfrak S(N)/(2h\sqrt{1-\eta_N})$.
* (2.30) and Thm 6.1: $R(N)\iff\|T_N\|\le\delta_N$. The proof of Thm 6.1 is the factor-by-factor quadratic, as the lane says.
* Thm 6.2, (2.24): $\|T_N\|\le\|T_N\|_{HS}\ll(\log N)^2$. The p. 12 sentence "they do not measure cancellation between the components" is present.
* §2.3.1: "GRH supplies holomorphy; the further spectral hypotheses will supply smallness." This is the paper's own phrasing of the lane's point (d).
* Prop. A.4: for fixed $N$, $\|T_{N,c}\|\ge z_N^{obs}/c$, and "$R(N,c)$ fails for all sufficiently small $c$".
* App. B: $N=30030$; 905 pairs recovered; $\|C\|/\delta_N\approx0.826864$ with $\delta_N=0.43231$; "they do not certify Hypothesis R".
* §8: "The central advance … is the exact identification $c_\phi\mathfrak S-G=\langle T_Nv_N,v_N\rangle$."
* (7.15)/(7.16): $Q_N=\frac{e^2}{2\pi}\int|V|^2du/u+O_\phi(1/\log N)$ with $\tilde V=(1+s)\tilde\Phi$, i.e. $V=\Phi-u\Phi'$.

**Re-derived:**
* S and R are restatements, by the paper's own equivalences.
* The conclusion is one coordinate, $\langle T v_N,v_N\rangle$, of the hypothesis $\sup_{\|v\|=1}|\langle Tv,v\rangle|$, so the hypothesis is strictly stronger.
* The GRH consequence $|G-c\mathfrak S|\ll(\log N)^2$ is weaker than the trivial $\ll\log N$.

The lane's (c), "all weights, all scales", rests on Prop. B.4 (same kernel for other weights). It is correctly labelled loose for the scale range. I agree with the verdict that the conditional theorem is a formal reformulation with no cancellation-producing step.

**Literature status.** WebSearch did not find "A conditional proof of the Goldbach conjecture" online. The same authors do have arXiv:2606.05224, "Theorem (1+1.9) on the Goldbach Conjecture", which is a different paper. The account therefore rests on the uploaded PDF alone; its publication status is unknown.

## 2. Re-derivation of the transplant propositions

**Prop. 2.0 (scale obstruction).** The inequalities $M\le r_2(N)\le4\tau(N)=R^{o(1)}$ are correct, and so is the comparison with the proved $(\frac23+\varepsilon)\log R/\log\log R$, provided $\log B\ge(1+\delta)\log\log R$ (the lane's "$\gg$" needs this constant).

*Correction.* The Consequence paragraph says every $Q^\varepsilon$-lossy input is vacuous "before any structure is examined". That is true only for statistics of **bounded degree** in $M$. Suppose an input bounds a degree-$d$ statistic, such as the $d$-tuple energy or the $L^{2d}$ norm of the toral eigenfunction restricted to the arc, by $B\,M^{d/2}$ (or $B$ times the trivial size) with $B=Q^{\varepsilon}$. Then $M\le Q^{O(\varepsilon/d)}$. With $d\asymp\log R$ and $\varepsilon$ fixed this is $O(1)$, which is the uniform bound.

So the corrected statement is: an input with a $Q^{\varepsilon}$ loss is vacuous **unless the loss is uniform in a degree $d\to\infty$**. Standard $\varepsilon$-loss results (Bombieri–Pila type, $L^p$ bounds of Bourgain–Rudnick type; memory-quote) carry $C_{\varepsilon,d}$ constants, and the record's criticality (Lemma K1) is where the $d$-uniformity breaks. Prop. 2.0 should be cited only together with this caveat.

**Prop. 2.1 (H, fake-blind).** The formula $\lambda_k=4\prod_pU_{e_p}(\cos4k\phi_p)$ is correct for $N$ composed of primes $\equiv1\ (4)$. In general there is an extra factor $(-1)^{ke_2}$ for $2^{e_2}$, since $\xi(1+i)=-1$, and a factor $1$ for $q^{2f}$ with $q\equiv3\ (4)$. I checked this exactly on 5 values of $N$ (including $2^3\cdot5^2\cdot13$ and $2\cdot9\cdot5\cdot13$) and 4 values of $k$ (`small_checks_out.txt`). Fake-blindness holds: the identity is a trigonometric-polynomial identity in the real angles. The GRH "graft" estimate holds (it is even better: the change is $\le2\sum_{p|N}\log p\le2\log N$). It is correctly labelled structural.

**Lemma 2.2 (Q, Parseval).** Correct, since $|n-n'|\le2W<q$.

**Lemma 2.3.** The conclusion is correct but the parametrisation is over-counted. Every $x\in\mathbb Q(i)$ of norm 1 is $x=P/\bar P$ with $P$ primitive, and $P$ is **unique up to sign**.
* Existence: Hilbert 90 with $\alpha=1+x$ (or $\alpha=i$ for $x=-1$), then remove the rational content.
* Uniqueness: $\alpha/\bar\alpha=\beta/\bar\beta$ forces $\alpha\bar\beta\in\mathbb Q$, so $\beta\in\mathbb Q\alpha$, so $\beta=\pm\alpha$.

So no unit is needed. Summing over $u$ and $\pm P$ counts each $x$ more than once, and "only $u=1$ contributes for non-trivial $I$" is just the choice of the unique representative. Non-triviality of $I$ plays no role. The final formula is right once read with the unique representation and $|\arg P|<\Theta/2$. Also, $\gcd(P,\bar P)\in\{1,1+i\}$, because a split or inert prime dividing both would contradict primitivity.

**Prop. 2.4 (resonance).** The proof is correct: $N-5A^2\le2\sqrt5R+5$ and $5s^2\le5R$, so $W\ge5R$ suffices. The angles are $\arg(ig)\pm s\sqrt5/R$. Recomputed independently (§3).

**Cor. 2.5.**
* (a) Correct.
* (c) Correct.
* (b) The heuristic $m\approx(R/W)^2$, $\sup E^W\gtrsim W^{3/2}/R$ comes from random hits $|mA^2-N|\le W$ of probability $W/(\sqrt mR)$ over $\asymp m$ directions. I re-derived it; it is internally consistent. But my data at $W=R$ do **not** show $\sqrt R$ growth: $\sup/\sqrt R=1.40,1.06,0.74$ for $\sqrt R=100,225,1000$. It must stay labelled heuristic. What the data do show is $\sup\gg$ mean: 740 against 57 at $R=10^6$.

**Prop. 2.6 (A, degenerate).** Correct. I checked it on 3 values of $N$ with $R>16C^2$: zero non-trivial collisions $a^2\equiv a'^2$. Two caveats:
* The identity counts only **primitive** points, since $(ab,N)=1$.
* "Exactly the Burgess threshold" (memory-quote) is loose. Burgess is non-trivial at length $N^{1/4+\varepsilon}$ only for prime or cube-free moduli, and at $N^{3/8+\varepsilon}$ in general. Here $N=R^2$ need not be cube-free, and the length is $CN^{1/4}$ with no $\varepsilon$.

**Lemma 2.7.** Correct to leading order from (7.15)/(7.16). The $O(1/\log N)$ error is uniform in $\delta$, because $(Nu)^{2\alpha_N}=e^2u^{2/\log N}$ on a fixed $u$-range. The lane's numbers match ($\to\|\phi'\|^2=0.40959$). Its relevance to the arc problem is by analogy, since no Hardy space is constructed for transplant A.

### 2.4 The self-similar form of $E^W$ (new here; proved)

Let $z\ne w$, $|z|=|w|$, and $w/z=P/\bar P$ with $P$ primitive. Write $P=gP'$, $g\in\{1,1+i\}$, $\gcd(P',\bar P')=1$. Then $z=\bar P'y$ and $w=\epsilon P'y$ with $y\in\mathbb Z[i]$ and $\epsilon=g/\bar g$ a unit. Hence, exactly,
$$E^W(I)=\sum_{P'\ne1}\#\{y\in\mathbb Z[i]:\ |P'|^2|y|^2\in[N-W,N+W],\ \bar P'y\in I,\ \epsilon P'y\in I\}.$$

Each term counts Gaussian integers $y$ on circles of radius $\approx R/|P'|$, in a norm window of width $2W/|P'|^2$, inside a sector of angular width $\Theta-2|\arg P'|$. That is an arc of length $\le C\sqrt R/|P'|=(C/\sqrt{|P'|})\sqrt{R/|P'|}$.

For $|P'|^2>2W$ the window holds at most one integer, so the term is **a single short-arc count of the original type**, at radius $R/|P'|$ and constant $C/\sqrt{|P'|}$. The central circle's own $M(M-1)$ is the sub-sum with $|P'|^2|y|^2=N$.

Consequences:
* (i) The "minor-arc counts conj(P)Z[i] in the thin strip" of the Lead are short-arc counts for $\mathbb Z[i]$ itself, so no homogeneous-dynamics input over $\mathrm{SL}_2(\mathbb Z)\backslash\mathrm{SL}_2(\mathbb R)$ enters. All the lattices are similar to $\mathbb Z[i]$, because $\mathbb Z[i]$ is a PID.
* (ii) Remark 3.1's hypothesis is a bound on a sum of instances of the problem itself, which is a proved form of "restatement".
* (iii) Large $|P'|$ gives rational-direction resonances (Prop. 2.4 is $P'$ near direction $g$, with many $y$ on a line).

## 3. Independent recomputation

**Li–Liu kernel.** `kernel_indep.py` uses a different implementation:
* mpmath Hurwitz $\zeta$ with the exact derivative instead of a finite difference;
* an explicit character DFT;
* the prime zeta via $\sum_k\mu(k)k(\log\zeta)''(ks)$;
* $P_N$ as primezeta plus a direct correction.

$\mathscr P$, $\mathscr K$ and $\mathscr D$ agree with the lane to about $10^{-5}$. My first $\mathscr D$ omitted prime powers; that bug was found and fixed. The table row $N=2018$, absent from the lane's output file, is confirmed: $E=0.2373$.

**Caveat.** The Hardy strip is $\alpha_N=1/\log N<\Re s<1/8$. For every $N\le2\cdot10^5$, $\alpha_N>0.05$, so **all $\sigma=0.05$ rows are outside the strip**. $\sigma=1/9$ is inside only for $N\ge20014$. The in-strip values are $E=0.234$ (20014), $-0.059$ (60022), $-0.015$ and $-0.0115$ (200006). The components are of size 2–19. So the evidence for R is four points and no visible decay law. It is consistent with R but weaker than the lane's phrase "cancel to $10^{-2}$" suggests.

**Position-averaged pair count** (`kappa_out.txt`; exact enumeration up to $6.3\cdot10^6$ strip points, $R$ from $7\cdot10^3$ to $10^6$, $C=4,8$, $W=R,R/10$):

| | |
|---|---|
| $A/((W/R)C^2)$ | 2.59→3.58 ($C=4$, $W=R$) and 2.79→3.73 ($C=8$, $W=R$) as $\log R$ goes 8.86→13.82 |
| slope | $\approx0.19$–$0.20$ |
| $\kappa_{\rm eff}$ range | 0.24–0.44 (0.44 at $W=R/10$, $R=7\cdot10^3$, outside the lane's range) |

A heuristic that I derived and corrected: density $8/\pi^2$ of primitive $P$ (weighted for $(1+i)\mid P$), the profile $\int(\Theta-2|\varphi|)_+d\varphi=\Theta^2/2$, and $|P|\in[2/\Theta,R]$. It gives
$$\mathbb E_\theta E^W\approx\frac{2}{\pi^2}\frac{W}{R}C^2\big(\log R+2\log\tfrac C2\big)+O\big(\tfrac WRC^2\big).$$
Its slope $2/\pi^2=0.203$ matches the data. So the lane's "$\kappa\in[0.25,0.33]$" is an effective finite-size value, and the leading constant is about $0.20$ (heuristic; the counts are exact).

**Resonance** (`resonance_out.txt`, $C=4$, independent sweep). At $W=5R$, the exact $E^W$ on the predicted arc at $153.43^\circ$ is 484, 980, 3904 for $\sqrt R=100,225,1000$. This is $\ge$ the proved lower bound $2\lfloor C\sqrt R/(2\sqrt5)\rfloor=178,402,1788$, so **Prop. 2.4 is confirmed**. The non-trivial supremum sits at the symmetric images $17.3^\circ$–$26.1^\circ$ of that direction. At $W=R$: $\sup=140,238,740$ against a mean of about 46–57. The lane's figure of 238 at $N=2576450045$ is reproduced.

## 4. Filters

**F1.**
* H: fake-blind, confirmed (§2, Prop. 2.1).
* Q: not fake-blind. $E^W\ge M(M-1)$ uses the fact that the central points are lattice points, which fakes lack. But by §2.4 everything Q supplies is the original problem at smaller radii, so it passes F1 only by containing the target.
* A: uses genuine integers $a,b$, but is degenerate (Prop. 2.6).

**F2.** I agree. The only inequality is Cauchy–Schwarz, and it is an equality at the cluster vector. No mechanism is present that could prove the false 4-point analogue of Question B. My own check: the Parseval identity of Lemma 2.2 holds verbatim for any finite point set with integer norms, so it carries no 3- versus 4-point information.

## 5. Literature (WebSearch)

* **Li–Liu, "A conditional proof…":** not found online. Only arXiv:2606.05224 by the same authors (a "(1+1.9)" paper) appears. The account above rests on the uploaded PDF.
* **Goldston–Gonek–Montgomery, J. reine angew. Math. 537 (2001) 105–126:** confirmed through the reference list of Chirre's paper in Proc. AMS and through the introductions of arXiv 2107.13636 and 2310.15918. The result: under RH, pair correlation is equivalent to the asymptotic for the mean square of $\zeta'/\zeta$ at $\frac12+a/\log T$. I did not read the full text.
* Littlewood, Chandee–Soundararajan, Burgess and Bourgain–Rudnick are memory-quotes.

## 6. Corrected statement and what survives

**Corrected statement.**
* Li–Liu's Hypotheses S and R are, by the paper's (2.28), (2.30) and Thm 6.1, operator-norm bounds strictly stronger than the Goldbach asymptotic. GRH supplies only holomorphy, and its norm bound $(\log N)^2$ is weaker than the trivial $\log N$.
* Transplanted to short arcs, the template gives three versions:
  * H is fake-blind (proved).
  * Q is the Parseval identity $\frac1q\sum|\Theta_a|^2=\#S_W+E^W$ with $E^W\ge M(M-1)$. By the exact self-similar decomposition (§2.4), $E^W$ is a sum of short-arc counts of the same problem on the circles of radius $R/|P'|$. Its supremum over non-trivial positions is $\ge(2/\sqrt5)\min(C,2\sqrt5)\sqrt R\,(1+o(1))$ for $W\ge5R$ (proved). Its position average is $\asymp(W/R)C^2\log R$ (exact counts; heuristic constant $2/\pi^2$).
  * A is exactly diagonal on primitive points for $R>16C^2$ (proved).
* Prop. 2.0 holds for statistics of bounded degree. A $Q^\varepsilon$ loss uniform in degree $d\to\infty$ would not be vacuous; no such input is known (record: Lemma K1).
* No GRH-type, pair-correlation, ratios, subconvexity or sup-norm hypothesis gives the uniform bound or $M=o(\log R/\log\log R)$ through this template. Remark 3.1, $M\ll C\sqrt{\log R}$, is a correct implication from a hypothesis that §2.4 shows to be self-referential.

**What survives:**
* the paper analysis;
* Props. 2.1, 2.4 and 2.6 and Lemmas 2.2 and 2.7;
* Prop. 2.0 in its bounded-degree form;
* the overall obstruction verdict.

**What does not survive:**
* the blanket "every $\varepsilon$-loss is vacuous";
* $\kappa$ as a constant in $[0.25,0.33]$;
* the $\sqrt R$ scaling at $W=R$, which is not supported by the data;
* the $\mathrm{SL}_2$-equidistribution framing of the Lead;
* the strength of the kernel evidence: the $\sigma=0.05$ rows are out of the strip.
