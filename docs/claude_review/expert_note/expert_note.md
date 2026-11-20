# Lattice points on short arcs: three Diophantine questions

*A short note for experts in Diophantine approximation, compiled on 2026-10-28 from the review record
`docs/claude_review/` of the repository `craigxchen/Jarnik` (branch `claude/upbeat-einstein-p0peju`). Each
statement below is the referee-corrected or killer-corrected version, and each carries a status label and a file
pointer. Paths are relative to `docs/claude_review/`.*

*Labels.* **Proved** means checked by an independent adversarial referee. **Proved (killer-checked)** applies to
round 20 only: the result was checked by an adversarial "killer" but not by a separate referee. **Proved by the
killer** marks the killer's own new results (round 20), which no second party has checked. **Conditional**
means proved from a named hypothesis. **Exact** means a finite exact computation. **Heuristic** means a count that
no proof uses. **Unrefereed** means checked by its author only. **Open** means neither proved nor refuted.

## 1. The problem and what is known

Let $N=R^2\in\mathbb N$ and $\mathcal E_R=\{z\in\mathbb Z[i]:|z|^2=N\}$. For $C>0$, let $M(C)$ be the supremum,
over all $R$ and over all arcs $I$ of $|z|=R$ of length $\le C\sqrt R$ **in any position**, of
$\#(I\cap\mathcal E_R)$.

> **Conjecture (Cilleruelo–Granville, critical case).** $M(C)<\infty$ for every $C>0$.

The scale $\sqrt R$ is critical because the sagitta $L^2/8R$ of an arc of length $L$ becomes of order $1$ there.
Arcs that touch an axis point are trivial, so the whole difficulty is uniformity in the position of the arc
(`overview/overview.md` §1.1).

**Known for every $R$.** The literature bound is $\ll_C\log R$ (Cilleruelo–Córdoba 1992; Bourgain–Rudnick).

> **Theorem 1.1** (**proved**, two referees; `two_thirds/proof.md` Thm 7.5, Cor 8.1). Let $M\ge2$ distinct lattice
> points lie on an arc of length $\le C\sqrt R$. Then $3M\log M\le\log R^2+(14+4\log^+C)\,M$. Hence for every
> $C,\varepsilon>0$ there is an explicit $R_0(C,\varepsilon)$ such that $M\le(\frac23+\varepsilon)\log R/\log\log R$
> for all $R\ge R_0$.

> **Theorem 1.2** (**proved**; round 31, Part VIII.1, a theorem of the outside colleague, checked by seven
> referees and a judge, 0 failures in 1.27 million exact window tests; `round31/README_VIII.md` §§1.1, 1.8).
> Divide a cluster of $M\ge2$ distinct points by its Gaussian gcd; let $N'$ be the normalised norm and $C'\le C$
> the chord constant $\max_{i,j}|z_i'-z_j'|/N'^{1/4}$. Let $S$ be a finite set of primes $\ell\equiv3\ (4)$ with
> $\frac12\log2+\sum_{\ell\in S}\frac{\log\ell}{\ell+1}>\log C'$. Then the exponent rows $(a_{ip})_{p\mid N'}$, each
> enlarged by the indicator of the class of $z_i'$ mod $\ell$ for every $\ell\in S$, are affinely independent.
> Hence $M\le1+\omega(N')+K_C$ with $K_C=\sum_{\ell\in S}\ell$ (one can take $K_C=\Theta(C^4/\log C)$), and
> $M\le1+\omega(N')$ when $C'<\sqrt2$.

The proof of Theorem 1.2 is short. An affine dependence splits into two probability measures on disjoint sets of
marks with equal barycentres. For independent marks $i,j$ drawn from them, $|x-y|\le x+y-2xy$ on $[0,1]^2$ gives
$\mathbb E\,d_{ij}\le W/2$ for the weighted exponent distance $d_{ij}=\sum_p|a_{ip}-a_{jp}|\log p$, while collisions
modulo the inert $\ell\in S$ (norm fibres of size $\ell+1$) and the factor $1+i$ give
$\mathbb E\log|b_{ij}|>\log C'$ for the cofactor $b_{ij}=(z_i'-z_j')/\gcd(z_i',z_j')$. This contradicts the pair
inequality $d_{ij}\ge W/2+2\log|b_{ij}|-2\log C'$. Consequences (**proved**, same source):
* $P_+(N')\ge(2-o_C(1))M\log M$ for the largest prime $P_+(N')$ of $N'$, so $M\le M_0(C,\eta)$ on clusters with
  $P_+(N')\le(2-\eta)M\log M$;
* $M\le1+K_C+\omega_1(\mathrm{rad}\,N)$, where $\omega_1$ counts primes $\equiv1\ (4)$; hence
  $M\le(\frac2\kappa+o(1))\log R/\log\log R$ with $\kappa=\log N/\log\mathrm{rad}\,N$. This beats $\frac23$ exactly when
  $\kappa>3$; for squarefull $R$ the coefficient is $\frac12$.

Scope: **for general $R$ nothing improves.** Together with the inputs of Theorem 1.1 the constant $3$ stays exactly
$3$ (`round31/ref8_extend.md` Prop 3.2). Theorem 1.2 replaces the record's $\omega(N')>M/4$ (`round8/stability.md`
Thm 3.3) and turns the conditional ceiling $M\le\omega(N')+O_C(1)$ of §5 items 5 and 8 into a theorem. It holds
verbatim over $\mathbb F_2[t]$ and $\mathbb F_3[t]$ and in every arc-admissible point-residue angle model with
residues up to $\max S$, so every fake of the record satisfies it and it cannot by itself give the uniform bound.
$K_C$ is lossy: on all circles $N\le10^6$ no affine dependence of the exponent rows occurs for $C'<2$, and
$M-1-\omega(N')\le6$ for arc constant $\le9$, against $K_C=1150$ at $C'=8.24$. Novelty relative to the literature was checked only by
snippet-level searches.

> **Theorem 1.3** (**proved**; round 31, Part IX.2, the outside colleague's projection theorem, two referees and a
> judge, 0 failures in 1.68 million exact tests; `round31/README_IX.md` §§1.2, 1.9(a), 3(a)). Let $M\ge3$ distinct
> lattice points lie on $|z|=R$, $N=R^2$. For $0\ne b\in\mathbb Z[i]$ the integers $\mathrm{Re}(\bar bz_i)$ have
> diameter $X_b$ with $2+X_b\gg M^{3/2}$ (absolute constant; $\gg M^2$ if no prime $p\equiv1\ (4)$, $p\le M$, divides
> $N|b|^2$). Hence an arc of length $\le C\sqrt R$ whose midpoint direction $\theta$ satisfies
> $|\sin(\theta-\arg b)|\le\kappa R^{-1/2}$ holds $\ll(|b|(C\kappa+C^2))^{2/3}+1$ points.

Theorem 1.3 is Gallagher's larger sieve applied to the projected points, which lie on the conic $x^2+y^2=N|b|^2$
and so occupy about half the residues modulo each prime. It is **position-dependent**: it gives
$M=o(\log R/\log\log R)$ only when $|b|(1+\kappa)=o((\log R/\log\log R)^{3/2})$ (the classical projection count,
`round4/hardness.md` Lemma 4.2(b), needs $o(\log R/\log\log R)$), and nothing for generic arcs, whose best Gaussian
direction has height about $R^{1/4}$. In row language: for integers $c_i$ with $\sum c_i=2$ and
$T(c)=\prod_pp^{|\sum_ic_ia_{ip}-e_p|}$ one has $(1+\|c\|_1)\sqrt{T(c)}\gg_CM^{3/2}$, so a combination with
$(1+\|c\|_1)^2T(c)=o(\omega(N')^3)$ gives $M=o(\omega(N'))$; two rows with $\sum_p|a_{ip}+a_{jp}-e_p|\le2$ (two marks
that are reflections of each other in a line of small height) force $P_+(N')\gg_CM^{3/2}$. The hypothesis is
exceptional: generic clusters have $T\approx N'^{0.2\text{–}0.3}$, the extremal data of Theorem 1.1 admit no small
certificate, and the record's large clusters are exactly symmetric, where the classical count already gives
$O_C(1)$. Angle models satisfying every point-residue axiom violate the bound in a rational-direction window
(Part IX.6), so it uses the integrality of the projected coordinates. The sieve step is probably folklore; no
source was found.

The proof compares an archimedean upper bound with an arithmetic lower bound for the product of the reduced chords.
After dividing by $\gcd(z_1,\dots,z_M)$ (the *normalised* cluster: $N$ odd, all $p\mid N$ split,
$z_j=\varepsilon_j\prod_p\pi_p^{a_{jp}}\bar\pi_p^{e_p-a_{jp}}$), put $c_{ij}=(z_i-z_j)/\gcd(z_i,z_j)$ (the reduced
chords) and $S_{p\tau}=\#\{i:a_{ip}<\tau\}$ for $1\le\tau\le e_p$. The link is an exact cut identity
(`overview/overview.md` §2.2, Step 2):
$\prod_{i<j}|z_i-z_j|^8=N^{M(M-2)}\prod_{p\mid N}p^{\sum_\tau(2S_{p\tau}-M)^2}\prod_{i<j}|c_{ij}|^8$.
The arithmetic side comes from pigeonhole collisions at primes $q\lesssim M$, together with a chain lemma on
exponent levels. The factor $\log M$ is the mass $\sum_{q\le M}\log q/q$. No bound $o(\log R/\log\log R)$ is
known unconditionally. Conditionally, the only known routes go through Vojta's Main Conjecture and give the
uniform bound directly; no intermediate growth bound is known from a weaker standard conjecture. In the other direction, the sharpest Vojta-type form of the problem, $(L_8)$, implies Rudin's conjecture on squares in arithmetic progressions (§2, "Strength of $(L_8)$"); the uniform bound itself is not known to imply any recognised open problem. They use Vojta's Main Conjecture for **one fixed** variety:
* $D=0$ on a resolution of $\mathrm{Bl}_YX_7$, where $X_7$ is the equal-norm variety and $Y$ its diagonal line
  (`growth/conditional.md` Thm 5.1);
* $D=0$ on $\overline M_{0,8}$ (§2 below);
* $D=E$, $S=\{\infty\}$ on a 4-fold (§3 below).

Explicit examples give $M(C)\ge10$ for $C\ge8.2394$ (one explicit 10-point cluster) and $M_\infty(C)\ge8$ for $C>5.89771$ (a
family) (`growth/families.md`, `growth/referee_families.md`).

**Why $\log R/\log\log R$ is a barrier for local inputs: fake circles** (`overview/overview.md` §§5.1, 5.4).
Divide a cluster by its Gaussian gcd. Then $N=\prod_jp_j^{e_j}$ with $p_j\equiv1\ (4)$, and
$z_x=\varepsilon_x\prod_j\pi_j^{a_{xj}}\bar\pi_j^{e_j-a_{xj}}$. The arc condition becomes a linear system in the
angles $\phi_j=\arg\pi_j$: for all $x$, $\langle2a_x-e,\phi\rangle\in\theta_0+[0,\Delta]+\frac\pi2\mathbb Z$, with
$\Delta\asymp R^{-1/2}$. A **fake circle** (angle model) keeps the genuine primes and the genuine exponent profile
but uses an arbitrary *real* $\phi$. That $\phi$ must solve the system and satisfy every monomial gap
$\mathrm{dist}(\langle v,\phi\rangle,\frac\pi2\mathbb Z)\ge\arcsin e^{-w(v)/2}$, where
$w(v)=\sum|v_j|\log p_j$; these gaps hold on actual circles. Since $W=\log N=cM\log M$ corresponds to
$M\approx(2/c)\log R/\log\log R$, an input class that admits fakes with $W=O(M\log M)$ cannot prove
$M=o(\log R/\log\log R)$. The record builds such fakes for several classes:
* Skew-Hadamard (Paley) profiles give fakes at $W\le(8+o(1))M\log M+8M\log(1/C)$. These satisfy every pair, cut,
  monomial-gap and rational-character inequality (`round8/lw.md` Thm 6.1; **proved**).
* Legendre (Paley) profiles with exponent cap two and randomly assigned primes in $[P,2P]$, $P=M^{2-\xi}$,
  $\xi<1/12$, give fakes **inside the low-height window** at $W\le(23/3+o(1))M\log M$ (along a diagonal sequence).
  They satisfy every monomial gap, every rational character and point residues (axiom (R) of `round5/construct.md`)
  at all odd prime powers $\le2M$, on a positive-measure set of non-arithmetic angles (round 31, Part VIII.3,
  (VIII.10); **proved**, four referee blocks and a judge's recheck of the mod-four support bound,
  `round31/README_VIII.md` §§1.3–1.5, 1.8). They have positive, not majority, Haar measure on the realisation
  component, and they are not checked against Lang–Waldschmidt. All margins are asymptotic ($q>10^{180}$ roughly).
  The same profiles, with primes below $M^{11/6+o(1)}$ and weight $W\le(22/3+o(1))M\log M$, carry prime-residue
  witnesses in $\mathbb Q(i)^*$ of the genuine norms that satisfy the character-tested factored residue axiom
  (P3-char, `round5/round6_sharp/sharp.md` Def. 1.1 with the same-level clause) **at every modulus at once**
  (round 31, Part IX.3–4; **proved**, `round31/README_IX.md` §§1.3–1.5). A separate family at larger prime height
  also satisfies the single-label max-scale cross-gap at constant one (Part IX.5; **proved** modulo one Fourier-decay
  citation not read in full; an elementary twelve-copy variant covers the finite test). A $128H^2$-fold replicated
  version also passes all single-label rays and all words of length $\le H$ with nonzero class sum of the
  all-monomial test, at $W\le(1536H^2/\rho+o(1))M\log M$ (Part X.2; **proved**). For the IX.3 witness design no
  right inverse gives short witnesses: almost every owner twist has log norm
  $\ge(\frac12-o(1))\log\mathrm{rad}\,N$ (Part X.3; **proved**).
* Residue designs extend this to point-level residues at moduli $\le M^A$ with $A$ fixed in advance
  (`round5/construct.md`, `round5/lower.md`). Referee scope limits (removed for zero-sum characters by round 31, Part IX.3): all primes of $N$ must exceed the modulus bound
  $X$, and the twisted extension at $X=M^A$, $A>1$, is open (`round5/README.md`; `overview/overview.md` §5.2 item 2).
* Grafts cap robust averaged analytic inputs (`round10/analytic.md` Thm B; `round20/composite.md`).

Actual clusters satisfy $W\ge3M\log M-O(M)$. Not everything is capped. P3-all$^*$ and the sharp-constant large
sieves (prime-only and non-multiplicative forms; the completely multiplicative, all-ideal, all-height sharp sieve
of round 20 collapses: it forces genuine angles mod $\pi/4$, with $B_0=4$ attained by an irrational graft,
ineffectively in the height; `round31/README.md` §1, sieve) remain uncapped; see `overview/overview.md` §5.5.
Residue caps uniform in the modulus are now proved for one class: residues that factor through prime residues of
the genuine norms and are tested through zero-sum characters (P3-char) are capped at every modulus at once
(round 31, Part IX.3, the bullet above), hence so are point residues at moduli $M^A$ with $A$ chosen after the
cluster, for zero-sum characters at odd prime powers. Within P3-all$^*$ (prime-level residues tested on every
monomial at moduli $3M^A$), one parent sequence at $W=O_H(M\log M)$ passes, at every modulus at once and together
with P3-char, every power of every single-label ray, every pair sum and every signed word of length $\le H$ whose
coefficient sum is nonzero on one of $b\to\infty$ label classes chosen by the construction (round 31, Part X.2;
**proved**, `round31/README_X.md` §1.2). Not capped: the balanced words inside a class of P3-all$^*$, in particular
label differences $e_j-e_k$ chosen from the data, twisted point collisions ($\sum c\ne0$) at moduli $M^A$, the 2-adic clause, the
operator-form union large sieve and the composite cross-gaps, every cross-gap inside the low-height window,
square-root-quality prime sums and box counts, Lang–Waldschmidt, and the Gaussian-integer cocycle of the reduced
chords. At the other end, with one rational lift per prime and all moduli at once, the all-monomial residue test at
the single unit $-1$ on powers of single primes holds only on actual clusters (round 31, Part IX.7; **proved**);
since such lifts exist at every finite modulus, this is the problem itself and says nothing about P3-all$^*$.
Effectively (Part X.1; **proved**): if all primes of $N$ are $\equiv5\pmod8$ and lie in $[P,2P]$, a witness that is
the reduction of $\pi_j\gamma/\bar\gamma$ with $\mathrm{Norm}\,\gamma\le P^h$ and passes the two diagonal tests on
$2e_j$ at a modulus $\ge2\sqrt2P^{1+2h}$ is a genuine prime factor, and with genuine factors the model is an actual
cluster (round 14, Prop. P). Passing data need not have such lifts: on one fixed fake every lift at label $p_j$ has
$\mathrm{Norm}\,\gamma>(X/(\sqrt2p_j))^{1/2}$ at modulus $X$ (Part X.2), so bounded-height rigidity yields no lower
bound for P3-all$^*$. Genuine witnesses pass every residue test on every parent, so residue tests can serve a growth
proof only through rigidity; **the residue route is closed for internal work** (`round31/README_X.md` §3(b)). So the
residue column is settled at both ends and open in the middle (`round31/README_IX.md` §3(c)–(d)). The low-height
window $P_+(N')\le M^{2-\delta}$ ($N'$ the normalised norm, $P_+(N')$ its largest prime factor) was the one regime
with no known fake (`round21/README.md`). It is now **capped for $0<\delta<1/6$** for pair, cut, monomial-gap,
rational-character, point-residue and character-tested prime-residue inputs (round 31, Parts VIII.3 and IX.4;
`round31/README_VIII.md` §3(d), `round31/README_IX.md` §3(b)). It stays **open** for $\delta\ge1/6$ and for the
uncapped inputs just listed. Actual clusters have $P_+(N')\ge(2-o(1))M\log M$ (Theorem 1.2). The record's synthesis is a judge's assessment, **unrefereed**
(`round17/README.md` §7.1; `overview/overview.md` §7.6). It says a new mechanism must be global, multi-place and
tied to the actual arithmetic of $N$. The only statements of that type in the record are Vojta-type, and their
sharpest arithmetic forms are Questions A and B below.

## 2. Question A: $(L_8)$, a statement about integer $2\times8$ matrices

**Data** (`growth/ptolemy.md` §7; `overview/overview.md` §7.1). The data are:
* a matrix $P=(P_0|\cdots|P_7)\in\mathbb Z^{2\times8}$ with primitive columns, $P_0=(1,0)^T$ and
  $P_j=(X_j,Y_j)^T$;
* positive integers $n_T$ for $\emptyset\ne T\subseteq[7]$, such that for all $0\le i<j\le7$
  $$\det(P_i,P_j)=t_{ij}\prod_{T\ni i,j}n_T,\qquad t_{ij}\in\mathbb Z\setminus\{0\}\quad(t_{0j}=Y_j);$$
* the **chain condition**: at each prime $p$, the sets $T$ with $p\mid n_T$ form two chains under inclusion, and
  the unions of the two chains are disjoint.

The *visible* blocks are those with $2\le|T|\le6$. There are $119$ of them, one for each boundary divisor
$D_T$ of $\overline M_{0,8}$. Singletons occur in no minor. The top block $[7]$ divides every minor with
$i,j\ge1$, but it can be inflated by a fresh prime ($z\mapsto nz$) without moving $\Phi$, so imposing balance
only on the 119 visible blocks is equivalent (`round18/README.md` §1, "Balance"). Column $0$ is the **anchor**: it
is normalised to $P_0=(1,0)^T$ and sent to $\infty$ in Kapranov's chart
$\mathbb P^5=\mathbb P(\mathbb Q^7/\mathbb Q\mathbf 1)$. The chain condition does not depend on the anchor, but $Z$
depends on the labelling, so a statement with an arbitrary anchor must symmetrise $Z$ (`round18/README.md` §1,
"Anchor"). Put
$H_\Sigma=\sum_j\log|Y_j|+\sum_{1\le i<j\le7}\log|t_{ij}|$. The datum is **$\eta$-balanced at scale $w$** if
$|\log n_T-w|\le\eta w$ for all 119 visible $T$. Finally, $\Phi\in M_{0,8}(\mathbb Q)$ is the point given by the
eight columns as points of $\mathbb P^1(\mathbb Q)$.

> **$(L_8)$** (**open**). There are $c>0$, $\eta_0>0$ and a proper Zariski-closed $Z\subset\overline M_{0,8}$ such
> that every such datum with $\Phi\notin Z$ and $\eta\le\eta_0$ satisfies $H_\Sigma\ge c\,w-O(1)$.

**Cleanest form** (`round18/theory.md` §1.3, Prop 2.1; `round18/census.md` Prop 1.3; **proved** reductions).
Every $\Phi\in M_{0,8}(\mathbb Q)$ has a primitive representative $x\in\mathbb Z^7$, taken modulo translation, with
point $0$ sent to $\infty$ and
$$|x_i-x_j|=\prod_{T\ \text{visible},\ T\ni i,j}d_T(\Phi).$$
Here $d_T=\prod_pp^{\ell_T(p)}$ is the **contact integer**, the finite part of the height $h_{D_T}(\Phi)$, and
$\ell_T(p)$ is the number of levels $t$ at which $T$ is a residue class of $x\bmod p^t$ (`round18/theory.md`
Lemma 2.2). Let $\mathrm{Bal}(\eta)$ be the set of $\Phi$ with $|\log d_T-w|\le\eta w$ for all 119 $T$. The
strongest natural form, $(L_8^c)^{\rm all}$, which implies $(L_8)$, says that **a rational point of $M_{0,8}$ off a
fixed $Z$ cannot have all 119 boundary contacts of the same size $e^{w}$, $w\to\infty$**; $(L_8)$ itself gives this
only at caterpillar points. Round 18 proved a sandwich, and no arrow is known to reverse:
$$(L_8^c)^{\rm all}\Rightarrow(L_8)\Rightarrow(L_8^f)\Rightarrow(L_8^c)^{\rm cat}\Rightarrow(L_8^{pc})'.$$
* $(L_8^c)^{\rm all}(\eta_1)$: $\mathrm{Bal}(\eta_1)\cap\{w\ge w_1\}\subset Z$. It implies $(L_8)$ with
  $\eta_0=\eta_1/2$, $c=\eta_1/4$. Conversely $(L_8)$ implies $(L_8^c)(\eta_0,Z)$ only for $\Phi$ admissible from
  anchor 0 (any anchor, with symmetrised $Z$); this is not an equivalence (`round18/README.md` §2.1 item 1;
  `overview/overview.md` §7.1 item 2).
* "cat" restricts this to points that are **caterpillar** at every prime. A compatible split system is a
  *caterpillar* if one side of each split can be chosen so that the chosen sides form a chain; equivalently, all
  internal edges of its tree lie on one path. The chain condition is exactly the caterpillar condition on the block
  system at each prime (`round18/theory.md` Lemma 1.1; exact on all 39208 trees with 8 leaves). Caveat (referee
  item H): the block system equals the $p$-adic tree of $\Phi$ only at primes dividing no residue.
* "pc" restricts it to pairwise coprime $d_T$. **As first printed this set is empty**: at $p=2$ exactly five of
  the 119 contact integers of every point are even, and at $p=3$ at least two are divisible by $3$ (round 31, Part
  IX.2 of the Question B thread; **proved**, refereed, `round31/README_IXL.md` §1.5). The corrected restriction
  $(L_8^{pc})'$ asks for $d_T$ pairwise coprime outside $\{2,3\}$, with the $\{2,3\}$-part of
  $\prod_{i<j}(x_i-x_j)$ at most a constant $B_0$.
* $(L_8^f)$ replaces $cw$ by any $f(w)\to\infty$ with $f\le cw$. Every proof of $(L_8^f)$ must prove
  $(L_8^c)^{\rm cat}$, the $H_\Sigma=0$ content of $(L_8)$, so a slowly growing $f$ is no shortcut; but $(L_8^f)$ is
  not known to be equivalent to $(L_8)$ (`round18/referee_theory.md` §2.3, item A).

**Lattice form** (`round18/theory.md` Lemma 2.2, §5(2), as corrected by `round18/referee_theory.md` item G and
by round 31, Part IX.2; **as first printed, vacuous; the corrected form below is open**; a consequence of $(L_8)$,
not known to be equivalent to it).
* In Kapranov's model $\overline M_{0,8}\to\mathbb P^5$: $H\cdot\beta=31$, $E_T\cdot\beta=1$,
  $K=-6H+\sum_{3\le|T|\le6}(|T|-2)E_T$ and $K\cdot\beta=10$. Here $\beta$ is the balanced class, defined by
  $D_T\cdot\beta=1$ for all $T$.
* The congruence lattice $\Lambda(x)$ consists of the $y\in\mathbb Z^7$ with $y_i\equiv y_j\bmod p^t$ whenever
  $i,j$ lie in one residue class of $x\bmod p^t$ of size $\ge3$. It contains $x$ and has index $[\mathbb Z^7:\Lambda(x)]=\prod_{3\le|T|\le6}d_T^{|T|-1}$.
* A balanced point has $|x|\approx e^{31w}$. The Minkowski scale of its own lattice is $e^{49w}$, and the trivial
  (Liouville) bound is $e^{30w}$.

> **$(L_8^{pc})'$, lattice form (corrected).** There are $\eta_1>0$ and a proper closed $Z$, and for every
> $B_0\ge1$ a $w_1=w_1(B_0)$, with the following property. Let $e_T\in[e^{(1-\eta_1)w},e^{(1+\eta_1)w}]$, for
> $T\subseteq[7]$ with $3\le|T|\le6$, be integers that are pairwise coprime outside the primes $2$ and $3$, and let
> $w\ge w_1$. Then the lattice $\Lambda_e=\{y\in\mathbb Z^7:e_T\mid y_i-y_j\ \text{for } i,j\in T\}$ has no primitive
> vector $y$ with distinct entries such that all of the following hold:
> * $[y]\in M_{0,8}\setminus Z$;
> * $|y|\le e^{31w(1+\eta_1)}$;
> * the pair cofactors $(y_i-y_j)/\prod_{T\supseteq\{i,j\}}e_T$ are integers, pairwise coprime outside $2$ and $3$,
>   of size $e^{w(1\pm\eta_1)}$;
> * the $\{2,3\}$-part of $\prod_{i<j}(y_i-y_j)$ is at most $B_0$.

*Correction (round 31, Part IX.2 of the Question B thread; refereed, `round31/README_IXL.md` §1.5).* The lattice
form first printed here asked for all $e_T$ and all pair cofactors to be pairwise coprime, with no prime exempt.
That statement has **empty input**: seven integers have at least nine pairs of equal parity; at most one $e_T$ is
even, and in every case at least five pair cofactors are even (exhaustive check), and at $p=3$ at least two are
divisible by $3$. So it was vacuously true, and $\mathrm{Bal}^{pc}$ as first defined is empty. The corrected form
follows from $(L_8)$ with $\eta_1=\eta_0/2$ (blocks are the $\{2,3\}$-free parts, the $\{2,3\}$-parts go into the
residues, $H_\Sigma\le\log B_0$), and its coprimality pattern is realised by integers, for example
$y=(0,5,10,107,267,69,-229)$ with $e_{\{1,2,3\}}=5$ and all other $e_T=1$. Below, $(L_8^{pc})$ means the corrected
form $(L_8^{pc})'$.

In words, structured lattices have no primitive vector $e^{18w}$ below their Minkowski scale. Nothing between the
trivial bound and Minkowski is known for these lattices. The heuristic count gives $98-108=-10=-K\cdot\beta$ per
unit of $w$: $e^{98w}$ profiles against $e^{-108w}$ for shortness (**heuristic**). $(L_8^{pc})$ alone does **not**
give the uniform theorem; that needs the chain cases, that is, $(L_8)$ itself.

**Theorems around $(L_8)$.**
* **$(L_8)\Rightarrow$ uniform bound** (**proved**; `growth/ptolemy.md` §7 and Steps 0–5 of §6.2, refereed in
  `growth/referee_ptolemy.md`). The proof has four steps.
  1. Subdivide the arc, keep one Gaussian unit class, and divide by its gcd.
  2. A random balanced selection (item 532 of the notes on branch `agent/uniform-gcd-formalization`) shows that at least half of all ordered 8-tuples have
     balance $\eta\le2\sqrt2\cdot255/\sqrt{M-7}$ and $H_\Sigma\le3584\,w/(M-1)$, with $w=W/128$.
  3. At most $c(Z)M^7$ of the 8-tuples have $\Phi\in Z$ (Lemma 6.1).
  4. $(L_8)$ then bounds $M$ independently of $R$.
* **Vojta $\Rightarrow(L_8)$** (**conditional**; Thm D, `growth/ptolemy.md` §6). Assume Vojta's Main Conjecture with
  $D=0$ on $\overline M_{0,8}$ over $\mathbb Q$, with $B=\sum_TD_T$ the boundary and $A=K+B$, which is ample
  (Keel–McKernan). Then $(L_8)$ holds with $Z=Z_{1/43}$ and
  $\eta_0=1/381$. The key estimate is $h_K(\Phi)\ge10w-282\eta w-60H_\Sigma-O(1)$. The sign is
  $K\cdot\beta=2^{k-3}(k-8)+k+2$, which is $-7$ at $k=7$ and $+10$ at $k=8$. So $k=8$ is the first case in which the
  $D=0$ inequality says anything at balanced points.
* **Smaller $k$.** $(L_5)$ and $(L_6)$ are false (**proved**; `growth/referee_ptolemy.md` §3.3, `threshold/L8.md`
  Prop C and Thm B). $(L_7)$ is open and
  expected false. $(L_7)\Rightarrow(L_8)$ (`round18/census.md` Lemma 7.7).
* **Correction B** (`round18/referee_theory.md` §2.7). A refutation of $(L_8)$ refutes Vojta at balanced points, but
  not conversely: a refuter must produce caterpillar points.

**Round-18 facts** (`round18/README.md` §§0–3; all **proved** and refereed unless labelled).
1. **No feasible computation** (`round18/census.md` Thm E). Each prime divides at most $5$ visible blocks, so a balanced datum
   needs at least 24 primes. At $\eta=0.1$ the class is empty for $e^w<59.2$, and all 21 minors
   $\det(P_i,P_j)$ with $i,j\ge1$ are $\ge10^{49.4}$. Block-sampling searches return $H_\Sigma\asymp w$ even at
   $k=5,6$, where bounded-$H$ points are known to be abundant (**exact**). So a failed search says nothing.
2. **No family mechanism.**
   * Density, box and sieve arguments produce only $\kappa\le-\sum\log e_T$, where $\kappa$ is the finite part of
     $h_K$ (Prop 2.3).
   * Bounded-degree covering families of genus-$\le1$ curves whose members *each* carry infinitely many balanced
     small-$H$ points are excluded (`round18/census.md` Prop 7.5). Correction M1: a dense supply with finitely many points on
     each member, for example on the lines of Kapranov's $\mathbb P^5$, is **not** excluded.
   * Class-$\beta_7$ lifts of the $(L_6)$ families are excluded (`round18/lift.md` Thm 4.2). Classes $c\beta_7$
     with $c\ge2$, and point-level lifts, are open.
   * Restricted theorems are known only on special loci with locus-dependent constants, or for a fixed prime set
     (Ridout: $H_\Sigma\ge8(1-\eta)(1-\varepsilon)w-C(S_0,\eta,\varepsilon)$ when all blocks are $S_0$-units;
     theory Thm 2.5). The constant $C$ depends on $S_0$, so this is useless here.
3. **Chains carry little information on profiles pairwise coprime outside $\{2,3\}$** (theory Prop 2.1, corrected
   in round 31, Part IX.2: $\mathrm{Bal}^{pc}$ as first defined is empty, because a single-edge tree is impossible at
   $p=2,3$). On the corrected class each prime $p\ge5$ divides at most one block and the $p$-adic tree is a single
   edge, so the chain condition is vacuous at $p\ge5$; at $2$ and $3$ the structure is carried by residues of total
   size $\le\log B_0$. **Every proof of $(L_8)$ must prove $(L_8^{pc})'$**. The gloss "nothing but
   $K\cdot\beta=10$ survives" is not a theorem (item G).

*Calibration (round 31, Part V).* The analogue of $(L_8)$ over $\mathbb F_2(t)$ with an anisotropic norm fails
asymptotically (plausible, not a theorem: the dictionary of `growth/ptolemy.md` Thm B was not re-derived in
characteristic 2), so a proof of $(L_8)$ cannot be characteristic-free. The same holds over
$\mathbb F_9(t)/\mathbb F_3(t)$ with the form $A^2+B^2$ (round 31, Part VI; same caveat).

**Strength of $(L_8)$** (round 31, Part VIII of the Question B thread; **proved** conditional implications, seven
referee blocks and a judge, `round31/README_VIIIL.md`). $(L_8)$ implies Rudin's conjecture: every nonconstant
$k$-term integer arithmetic progression contains $O(\sqrt k)$ squares. The proof uses $(L_8)$, with its linear
bound $c\,w$, anchor $0$ and unsymmetrised $Z$, only on data whose blocks are pairwise coprime perfect squares, one
block per prime, with arbitrary nonzero cofactors; it does not use $(L_8^f)$ with sublinear $f$,
$(L_8^c)^{\rm cat}$ or the lattice form $(L_8^{pc})$. It passes through three statements.
* (i) If at most two residues are prescribed modulo each exact prime power $p^e\parallel q$, any interval of length
  $q^{1/2+\epsilon}$ contains $O(1)$ integers in the prescribed classes, for one absolute $\epsilon>0$. This is the
  case $d=2$ of Cilleruelo–Granville's Conjecture 10, which they state at every exponent $1-\varepsilon$ (for initial
  segments $[0,b^{1-\varepsilon})$; they prove exponent $(1-\varepsilon_d(k))/d$ and "strongly believe" exponent
  $1/d$), taken at one fixed exponent above $1/2$. The step from
  $(L_8)$ is the colleague's: columns $P_i=(Q/d_i,(y_i-y_0)/d_i)^T$ with $d_i=\gcd(Q,y_i-y_0)$ give an exact
  $(L_8)$ datum, and a random 8-tuple of the integers is balanced on all 128 patterns and lies off $Z$ (a one-sided
  Gram bound, a Walsh expansion, Markov, and Schwartz–Zippel on the grid; every constant is fixed before $q$).
* (ii) Hence any interval of length $\sqrt q$ contains $O(1)$ roots of $x^2\equiv a\pmod q$, for every $a$ (their
  Conjecture 9).
* (iii) Hence Rudin's conjecture (their implication "Conjecture 9 $\Rightarrow$ Conjecture 1"; Theorem 11 in
  arXiv math/0608109v1).

Full $(L_8)$ also implies that an invertible residue class modulo $q\ge n^{1/4-\kappa}$ contains $O(1)$ divisors of
$n$, for one absolute $\kappa>0$: past the exponent $1/4$ of Lenstra and Coppersmith–Howgrave-Graham–Nagaraj (a
cardinality bound, no algorithm). Under Vojta's conjecture on $\overline M_{0,8}$ one may take
$\epsilon=10^{-16}$, $\kappa=10^{-17}$, so **Vojta's conjecture for the single variety $\overline M_{0,8}$ implies
Rudin's conjecture**; we found no earlier derivation of Rudin's conjecture from Vojta's conjecture, Bombieri–Lang,
the uniformity conjecture or abc (literature check of full texts, `round31/lit8L_verify.md`: Bombieri–Lang gives
$\sigma(k)\ll k^{5/7}$, Cilleruelo–Granville §3; the Caporaso–Harris–Mazur uniformity conjecture gives exponent $2/3$,
Shkredov–Solymosi; both are weaker than the unconditional exponent $3/5$ of Bombieri–Zannier; not opened:
Bombieri–Granville–Pintz, Bombieri–Zannier, Caporaso–Harris–Mazur). None of these consequences is
known to fail, and a referee who tried to refute $(L_8)$ through them found nothing: the best integer two-class
families reach exponent $2$ only for at most six points, and the elementary family $\prod_{i\le4}(1+iq)$ gives $16$
divisors at exponents just below $1/4$ (our computation; families of this type give at most
$2^{\lfloor1/\alpha\rfloor}$ at a fixed exponent $\alpha$, and the published examples we read have $6$ divisors, at
$\alpha>1/3$). They measure the strength of $(L_8)$, not its plausibility. The uniform bound itself is not known
to imply Rudin's conjecture, and Cilleruelo–Granville do not derive Rudin from any arc conjecture (their routes are
two congruence conjectures). Literature status of the divisor bound (`round31/lit8L_verify.md`): Hales (arXiv:2410.05030, read in
full) restates Lenstra's and Coppersmith–Howgrave-Graham–Nagaraj's bounds for exponents $>1/4$ and calls $1/4$ "the
threshold of both the Lenstra and Coppersmith et al. methods"; we found no result in either direction for invertible
classes at exponents $\le1/4$. Lenstra's unboundedness remark concerns a non-invertible class. Lenstra 1984 and
Coppersmith–Howgrave-Graham–Nagaraj 2008 themselves were seen only through snippets and Hales's restatement (their
hosts were unreachable); the published version of Cilleruelo–Granville was not read, only arXiv v1, which is the
only arXiv version. Cilleruelo–Granville's Theorem 10 should be cited for $k\ge3$ (it fails as printed for $k=2$).

> **Question A.** Prove $(L_8)$, or at least the corrected lattice form $(L_8^{pc})'$ (coprimality outside the
> primes $2$ and $3$, bounded $\{2,3\}$-part), which every proof must cover; the lattice form as first printed had
> empty input. Any proof
> of $(L_8)$, even restricted to pairwise coprime square blocks, proves Rudin's conjecture, and full $(L_8)$ also
> proves a divisor bound below Lenstra's exponent $1/4$ (see "Strength of $(L_8)$" above); the lattice form
> $(L_8^{pc})'$ is not known to imply either. A proof of
> $(L_8)$ would be an unconditional case of $D=0$ Vojta on the rational 5-fold $\overline M_{0,8}$, for a structured
> but unbounded family of points. Alternatively, find, for every $\eta>0$ and every $w_1$, a Zariski-dense set of
> $\eta$-balanced points with $w\ge w_1$ that are caterpillar up to $o(w)$ of residue. That would refute both $(L_8)$
> and Vojta (Correction B; `round18/referee_theory.md` §2.7, `round18/theory.md` §3.1).

## 3. Question B: $(L_5^\circ)$, five lattice points on a circle

Every theorem in this section is **proved (killer-checked)** or **conditional**, except where labelled otherwise; the
statements $A_5$, $(L_5^\circ)$ and its arc form are **open**. Sources: `round20/gauss.md` §§1–4,
`round20/kill_gauss.md` §§1.2, 6–7, and `overview/overview.md` §§4.3, 7.1A.

**Setting.**
* Let $\zeta_1,\dots,\zeta_k$ be a **primitive** tuple of distinct lattice points of $|\zeta|=R''$, meaning
  $\gcd_{\mathbb Z[i]}=1$.
* Its similarity class $[\zeta]=(\zeta_j/\zeta_1)_j$ is a point of the norm-one torus
  $T_k=(R^1_{\mathbb Q(i)/\mathbb Q}\mathbb G_m)^k/\text{diag}$, which has dimension $k-1$. $T_k$ is the open torus of
  the twisted Losev–Manin space $\overline L{}^{\rm tw}_k$ (the permutohedral fan in $\mathbb Z^k/\mathbb Z\mathbf 1$,
  with Galois acting by $-1$ on cocharacters). $[\zeta]$ is invariant under $\zeta\mapsto\gamma\zeta$,
  $\gamma\in\mathbb Q(i)^*$, which is why the exceptional set $Z$ must live in $T_k$. For $A\subsetneq[k]$ nonempty,
  $E_A$ is the toric divisor of the ray $\mathbf 1_A$.
* **Unit class.** Fix the Gaussian primes $\pi_p$ once. Each point factors as
  $\zeta=\varepsilon\prod_p\pi_p^{a(p)}\bar\pi_p^{e_p-a(p)}$ with a unit $\varepsilon\in\{\pm1,\pm i\}$; points lie in
  one **unit class** when their units agree. Within one class, $\zeta_a/\zeta_b=P/\bar P$ with $P\in\mathbb Z[i]$
  conjugate-primitive (`growth/ptolemy.md` §1, following item 532 on branch `agent/uniform-gcd-formalization`;
  `docs/claude_review/` has no separate formal definition).
* Each threshold layer $(p,t)$ of the exponent profile gives a cut $\{A,A^c\}$ of $[k]$. Let $n_A$ be the product
  of the $p$ over the layers with that cut. Then $h_{E_A+E_{A^c}}([\zeta])=\log n_A$ exactly, and
  $\sum_A\log n_A=\log R''^2$ (`gauss.md` Lemma 1.1).
* A short arc means archimedean proximity of $[\zeta]$ to the identity $e$.
* Put $d_{ab}=\sum_{A\text{ separates }a,b}\log n_A$, define $t_{ab}$ by $|\zeta_a-\zeta_b|=2|t_{ab}|R''e^{-d_{ab}/2}$,
  and put $H_\Sigma=\sum_{a<b}\log|t_{ab}|$. Here $t_{ab}$ is a nonzero integer when $\zeta_a,\zeta_b$ lie in one Gaussian
  unit class. Across unit classes it need not be: $|t_{ab}|\in\sqrt{\text{odd}/2}$ occurs (killer correction, minor 1).
  The record allows either repair: restrict to one unit class (as below), or take $t_{ab}:=|c_{ab}|/2$ real, with
  $c_{ab}=(\zeta_a-\zeta_b)/\gcd(\zeta_a,\zeta_b)$, and read the conclusion as $H_\Sigma\ge cw-O(1)$
  (`round20/kill_gauss.md` §7 item 1; `round20/README.md` §1.1 item 1).

> **$A_k(\delta)$.** There are a proper Zariski-closed $Z\subset T_k$ and an $R_0$ such that every primitive
> $k$-tuple of distinct lattice points of a circle of radius $R''\ge R_0$, with all mutual distances
> $\le R''^{1/2-\delta}$, has $[\zeta]\in Z$.
>
> **$(L_k^\circ)$.** There are $c,\eta_0>0$, $w_0$ and a proper Zariski-closed $Z\subset T_k$ with the following
> property. Every primitive $k$-tuple **within one unit class** with $[\zeta]\notin Z$, $w\ge w_0$ and
> $|\log n_A-w|\le\eta_0w$ for all $2^{k-1}-1$ cuts has $H_\Sigma\ge cw$.

At $k=5$ there is also a sufficient **arc form** (**open**; it implies $(L_5^\circ)$ by the argument of `gauss.md`
Prop 3.1, and the converse is not proved; `overview/overview.md` §7.1A). *Take five primitive lattice points of a circle of radius $R$, in
one unit class, whose Gaussian profile is $\eta_0$-balanced over all 15 cuts and whose class lies outside $Z$. Then,
once $R$ is large, they do not all lie within $R^{7/15+c}$ of each other.* The exponent $7/15$ arises because each
pair is separated by 8 of the 15 cuts. So $|t_{ab}|\ge1$ (valid within one unit class) gives diameter $\ge2R''^{7/15}$ at exact balance, and
this is exactly the Liouville bound. For a balanced $k$-tuple taken from a large cluster on a critical $\sqrt R$-arc, dividing by the tuple gcd moves it
to the "subcritical" exponent $\frac12-\frac1{2^k-2}$ in the reduced radius (`gauss.md` Lemma 1.2).

**Theorems.**
* **Reduction** (`gauss.md` Thm 3.4′, Prop 3.1). Let $k\ge3$ and $0<\delta<1/(2^k-2)$. Then
  $A_k(\delta)\Rightarrow(L_k^\circ)\Rightarrow$ the uniform bound. For $\delta\le2^{-k-1}$, only a convexity
  averaging of the tuple gcd is needed. The converse $(L_5^\circ)\Rightarrow A_5$ is **not** proved.
* **The 4-point version is false** (`gauss.md` Thm 4.1, plus the killer's projection corollary; exact Jacobian rank
  reproduced at 6 rational points).
  * Take Gaussian linear blocks $H_T(s)=q_T(s-r_T)$, with roots $r_T$ on an explicit rational 2-parameter family.
    These give balanced 4-tuples with bounded residues at exponent $3/7=\frac12-\frac1{14}$, dense in $T_4$.
  * Hence $A_4(\delta)$ is false for every $\delta<1/14$, and $(L_4^\circ)$ is false for every $c,\eta_0,Z$.
  * In fact **$A_4(\delta)$ is false for every $\delta<1/6$**, the full non-vacuous range (round 31, Part VI.4 of the
    Question B thread; refereed and reproduced with exact arithmetic, `round31/README_VIL.md` §1). Pell families
    $G_j=A^{r_jn}v_j$ with $A=\begin{pmatrix}5&2\\2&1\end{pmatrix}$ give four points
    $\overline{G_1G_2G_3}\cdot\{1,-ig_1g_2,-ig_1g_3,-ig_2g_3\}$, Zariski dense in $T_4$, on arcs of length
    $O(R^{1/3})$ at rates $(1,1,1)$ and $R^{2/5}$ at rates $(6,7,7)$. The families are Cilleruelo–Granville's; the
    statement that they leave every closed set is new. They are consistent with the three-point theorem (smallest
    limiting constant found $4.60154>4.34737$). The $1/14$ statement remains the best for the balanced
    $(L_4^\circ)$.
  * Likewise $A_3(\delta<1/6)$ and $(L_3^\circ)$ are false (the killer's own projection corollary,
    `kill_gauss.md` §3, `k3_projection.py`; not independently checked).
  * So $k=5$ is the first case where the reduction can work.
* **Vojta predicts it** (**conditional**; `gauss.md` Thm 4.2, re-derived by the killer).
  * Let $Y_k=\mathrm{Bl}_e\overline L{}^{\rm tw}_k$, so $K_{Y_k}=-\pi^*\sum_\rho E_\rho+(k-2)E$.
  * $V_E(k,\varepsilon)$ is Vojta's Main Conjecture for $(Y_k,D=E,S=\{\infty\})$: off a proper closed
    $Z_\varepsilon$,
    $\lambda_{E,\infty}(y)+h_{K_{Y_k}}(y)\le\varepsilon h_A(y)+C_\varepsilon$.
  * Here $A=\pi^*L-\varepsilon'E$ is ample with $L=\sum_\rho a_\rho E_\rho$, $a_\rho>0$, and $c_A=\max_\rho a_\rho$, so
    that $h_A\le c_AW'+O(1)$ with $W'=\log R''^2$ (`gauss.md` Thm 4.2 proof).
  * If $k\ge5$, $\kappa=(k-1)(1+2\delta)/4-1>0$ and $\varepsilon<\kappa/c_A$, then $V_E(k,\varepsilon)\Rightarrow A_k(\delta)$.
  * At $k=5$, $\kappa=2\delta$, so $V_E(5)$ with one $\varepsilon<(1/32)/c_A$ gives the uniform bound (via
    $A_5(1/64)$), and $V_E(5)$ predicts $(L_5^\circ)$.
  * The instance is non-vacuous (`gauss.md` Thm 4.2, Remark (3)). The $D=0$ form would need $k\ge6$.
  * The literature's instance (Vojta on $X$ with a normal-crossings divisor of multiplicity 4 through $e$) gives the
    same $\alpha^{\rm ess}\ge4$ (`papers_check/README.md` §1.3).
* **Covering families are excluded** (`gauss.md` Thm 4.3(ii); $k\ge5$, $\delta>0$). Put $B=\pi^*\sum_\rho E_\rho$, so
  $h_B=W'$. For each $d$, the curves $C$ with $B\cdot C\le d$ that carry infinitely many "A-points" lie in a proper
  closed set. Here an A-point is the class of a primitive tuple with diameter $\le R''^{1/2-\delta}$. (Thm 4.3(i) was
  narrowed by the killer, `kill_gauss.md` §1.2 item 4.) The density mechanism is excluded at $k=5$ only
  **heuristically** (count $e^{-w}$ per scale).
* **Lemma K1: every recorded method that yields an inequality at the actual point is exactly critical**
  (`kill_gauss.md` §6.1; **proved by the killer**, not independently checked; route 5's forced vanishing is
  supercritical but yields no inequality without its open zero estimate). Let
  $F$ be bihomogeneous of bidegree $(a,b)$.
  * Then $F(z)=g^a\bar g^bF(\zeta)$, and the forced divisor scales in the same way.
  * The angular width is unchanged, and so are the pair data and the residues.
  * Consequently, every forced-divisibility, Liouville, pair, cut, character or single-form inequality is the same
    before and after division by the gcd. The move from $\frac12$ to $\frac7{15}$ is a change of units.
  * Auxiliary forms are capped at ratio exactly $1$ (`threshold/construct.md` Cor C3, with Thm B), and single-form
    Lang–Waldschmidt is exactly tight.
  * Only Vojta's $D=E$ term and the heuristic count see the reduction.

**Approximation reading** (`kill_gauss.md` §6.2). Measure height by $h_B=\log R''^2$. An A-point satisfies
$\mathrm{dist}(y,e)\le H_B(y)^{-(1/4+\delta/2)}$, approximating the **rational** point $e$ of the 4-fold $T_5$.
* $\delta=0$ is counting-critical (**heuristic**).
* Liouville excludes balanced points for $\delta>1/30$, and arbitrary points for $\delta>1/10$ (Cilleruelo–Córdoba).
* What is needed is that, for some $\delta<1/30$, **the Liouville exponent at $e$ is not attained by a Zariski-dense
  set of balanced points**.
* At $k=4$ this fails (Thm 4.1). On $\mathbb P^1$, Liouville is attained at every rational point.
* Roth- and Schmidt-type theorems apply but are weaker. McKinnon–Roth (Invent. Math. 200 (2015)) Thm 5.1 gives,
  off a proper closed set, $\alpha^{\rm ess}_x\ge\beta_x(L)$, where
  $\beta_x(L)=\int_0^{\tau}\mathrm{Vol}(\pi^*L-\gamma E)/\mathrm{Vol}(L)\,d\gamma<\tau_x(L)$ (ibid. p. 25). Their Thm 6.3
  gives $\alpha_x\ge\varepsilon_x/2$, strict at rational $x$. At $e$ (on $X'$), $\beta_e(-K)\in[2.31,2.61]$ (inference),
  below Liouville's $10/3$. As McKinnon–Roth note (ibid. p. 24), at a $k$-rational point Liouville $\alpha_x\ge\varepsilon_x$
  is "a stronger inequality than $\alpha_x\ge\frac12\varepsilon_x$" (`papers_check/README.md` §2.3).
* So this is the McKinnon–Roth regime. According to the killer, the only predictions beyond Liouville/Seshadri are
  conjectural: Vojta, and McKinnon's conjecture that best approximations lie on rational curves. The killer also
  knew of no unconditional theorem that beats Liouville at a rational point of such a variety on a Zariski-dense
  set. Round 22 checked this at snippet level (two independent searches plus an offline arXiv mirror): no
  unconditional theorem with $\alpha^{\rm ess}>\tau_x$ at a rational point was found for any variety
  (`round22/README.md` §2, J3). Seven full texts were then read (McKinnon–Roth 2015, 2016; the McKinnon–Satriano
  erratum; Lehmann–McKinnon–Satriano 2024; Huang ANT 2021, Acta Arith. 2019, Bull. SMF 2020). Every unconditional lower
  bound on $\alpha^{\rm ess}$ in them is Liouville with a section, at most $\tau_x$. Every exact value of $\alpha^{\rm ess}$
  they compute equals $\tau$, with $\tau\ge\dim$ (`papers_check/README.md`).

**In the approximation-constant literature's own terms** (round 22; `round22/README.md` §§0, 4; `round22/toric.md`,
`round22/lit.md` with their killer reports). The 4-fold is the twisted Losev–Manin space $X=\overline{L}_5^{\rm tw}$. It
is smooth and **weak Fano**, not Fano: $-K$ is nef and big but not ample, with polytope the $A_4$ root polytope, and
$-K=\phi^*(-K_{X'})$ for the small crepant contraction $\phi$ onto the twisted Gorenstein toric Fano $X'$, with
$X(\mathbb Q)=X'(\mathbb Q)=T_5(\mathbb Q)$ and $(-K)^4=70$ (**proved**). With $H=H_{-K}=R''^2$, an A-point at
$\delta$ has McKinnon–Roth exponent $\gamma=4/(1+2\delta)$, and:
* $\varepsilon_e(-K)=\alpha_e(-K)=2$ (**proved**), attained only on diagonal conics, which carry repeated points;
* $\tau_e(-K)=10/3$ (**proved**, round 23 `width.md` Thm 2.5, refereed). For every $m$, the maximal order at $e$ of a
  section of $m(-K)$ is $\lfloor10m/3\rfloor$, for the split and the twisted form, and it is attained by Vandermonde
  products. The proof slices the lattice points of $a\Delta-b\Delta$ coordinatewise, which gives the same family again,
  and peels the slices. Consequently every finite family of auxiliary Laurent polynomials over $\mathbb Q(i)$, used
  through Liouville with toric forced divisibility, gives for arbitrary profiles exactly the Cilleruelo–Córdoba exponent
  $\delta=1/10$ and never more; for balanced profiles it gives $\delta=1/30$. McKinnon–Roth's base-locus Liouville at $e$
  is therefore exactly Cilleruelo–Córdoba. (McKinnon–Roth, Eur. J. Math. 2 (2016), Thm 3.3. It needs $L$ ample, so it is
  applied to $-K_{X'}$ on $X'$. The values are unchanged.) The open window $\alpha^{\rm ess}\in(15/4,4)$ is untouched.
* **Upper bound: $\frac{10}3\le\alpha^{\rm ess}_e(-K)\le4$** (round 31, Part VII; drawn by two referees and the
  judge from the colleague's polynomial Pell families, then checked by three independent sceptics with exact
  arithmetic, `round31/ref7L_alpha_exact.md`, `ref7L_alpha_density.md`, `ref7L_alpha_record.md`; **proved**).
  Exact statement: for every $n\ge1$ and every proper closed $Z\subset T_5$ there are rational seeds $s$ and a
  constant $C(s)$ such that, for all but finitely many integers $T$, the five marks of the Pell family at $t=T$,
  divided by their Gaussian gcd (which divides a Bezout integer depending on $s$ only), form a primitive
  five-tuple of distinct lattice points of one circle with class outside $Z$, mutual distances
  $\asymp_sR''^{1/2+1/(2(8n+1))}$ and $\mathrm{dist}^{4+1/(2n)}H_{-K}\le C(s)$. Hence **$A_5(\delta)$ is false
  for every $\delta<0$**. The argument is curve by curve (seeds chosen after $Z$, then $T\to\infty$); the
  constants are not uniform in the seeds; the tuples use 5 of the 15 cuts and are in general not in one unit
  class. Nothing follows for $\delta\ge0$, for $(L_5^\circ)$, or for exponent exactly 4 on a dense set with one
  constant. This is on the far side of Question B, which needs $\alpha^{\rm ess}>15/4$; it is a calibration
  (Vojta predicts exactly 4) and a short corollary of the Cilleruelo–Granville Pell mechanism. One sceptic's own
  further observation, **computed, not refereed**: at a fixed integer $T$ and $n\to\infty$ each curve has
  exponent exactly 4 with a seed-dependent constant.
* $A_5(\delta<1/30)$ is exactly $\alpha^{\rm ess}_e(-K)>15/4$ off a proper closed set, and $(L_5^\circ)$ is the same
  over balanced sequences; Vojta predicts $\ge4=\dim X$.
* *Calibration (round 31, Part V; item 12 of §5).* The analogues of $A_5(\delta<1/30)$, $(L_5^\circ)$ and $V_E(5)$
  over $\mathbb F_2(t)$ with the norm of $\mathbb F_4(t)/\mathbb F_2(t)$ are **false** ($15/4$ is attained
  Zariski-densely). The same holds over $\mathbb F_9(t)/\mathbb F_3(t)$ with $A^2+B^2$ (Part VI; every
  5-subtuple of the ternary clusters is balanced, by Weil; the Schwartz–Zippel step was not rechecked in
  characteristic 3). The number-field statements are untouched, but no proof of them can be characteristic-free.
* *Characteristic-0 function field (round 31, Part VII of the Question B thread; `round31/README_VIIL.md`).* The
  analogue of $A_5(\delta)$ over $\mathbb C(t)$ (five distinct primitive $z_i\in\mathbb C[t]$ of degree $D$ with
  one norm $z\bar z$, contact $h=D-\max\deg(z_i-z_j)\ge(\frac12+\delta)D$; is the ratio curve in one proper closed
  subset of $T_5$, uniformly in $D$?) is **open**. *Proved:* (i) it holds for every $\delta>0$ on ratios that are
  quotients of Laurent polynomials with at most $b$ monomials, with one ineffective $\mathbb Q$-closed
  $Z_{b,\varepsilon}$ for all degrees, coefficients, exponents and cancelled common factors (threshold
  $h\ge(\frac14+\varepsilon)W$ in the height $W=2D$; the constant $\frac14$ is optimal as $b\to\infty$); (ii)
  $20h\le9D+18r$, where $2r$ roots of $N$ lie off a circle with real centre (false for lines; no integer analogue);
  (iii) at $2h>D$ every four marks are linearly independent; (iv) a finite-type family dominating $T_5$ has
  $4h+1\le\#\{\text{charged places}\}\le W$ (characteristic 0 only: false over $\mathbb F_4(t)/\mathbb F_2(t)$,
  where $W=30$, $h=8$ dominates and the sparse statement fails with a uniform term budget). *Dense families:*
  polynomial Pell families ($A=\begin{pmatrix}t^2+1&t\\t&1\end{pmatrix}$) give five marks with $D=8n+1$, $h=4n$,
  i.e. $W=4h+2$, one unit short of (iv) and half a degree below the threshold; four marks with $D=6n+1$, $h=4n$;
  three marks with $h=2D/3$. *Data (referee):* rigid real five-mark configurations with $h=(D+1)/2$ exist at
  $D=11,13$ (numerical, 40 digits, not interval-certified), so the strict inequality is not vacuous; none at
  $D=7$, $h=4$ (exact). Published function-field tools (Mason–Stothers, Brownawell–Masser, Corvaja–Zannier,
  Guo–Nguyen–Sun–Wang) pay $\chi_S\asymp W$ and are vacuous here. The question is an open function-field Vojta
  instance of the same generic complexity as Question B; it is **not** a stepping stone.
* *The ten-pair slice $(\ast)$ is not a shortcut.* $(\ast)$: are genuine primitive 5-tuples on the ten-pair-cut
  profile at the Cilleruelo–Córdoba endpoint $\alpha=10/3$ (five lattice points on an arc of length
  $R^{2/5+o(1)}$) Zariski dense in $T_5$? **Open.** Non-density in $(\ast)$ is **equivalent** to $A_5(\delta)$ for
  some $\delta\in(0,1/10)$ (round 31, Part VI.6, **proved**; `round31/README_VIL.md`). It is logically weaker than
  Question B, but no separation is known: it asks to beat the Liouville/Cilleruelo–Córdoba exponent at five points
  at all. The profile needs $H\ge(5\cdot13\cdot17\cdot29\cdot37)^2\approx1.4\cdot10^{12}$; a non-exhaustive
  search to $H\approx4.6\cdot10^{15}$ found at most 8 distinct pair cuts and best width exponent $0.2299$ against
  $0.3$. The heuristic count predicts finitely many such tuples.
* **Unconditional on a thin class** (round 31, Part VI.1 of the Question B thread; **proved**, refereed, the judge
  rechecked the Jacobian step). Fix a number field $K\supset\mathbb Q(i)$, a finite set $S$ of places containing the
  infinite ones, $b\ge1$ and $\varepsilon>0$. There are an **ineffective** proper closed $Z(K,S,b,\varepsilon)\subset
  T_5$, $\rho>0$ and $H_0$ such that every primitive 5-tuple of norm $H\ge H_0$ lies in $Z$ if: its anchored ratios
  are $q_i=\eta_i\prod_e(G_e/\overline G_e)^{b_{ie}}$ with at most one Gaussian block per cut,
  $\prod_e|G_e|^2\le H$, $h(\eta_i)\le\rho\log H$; each $G_e$ is a sum of at most $b$ terms $a\cdot u$ with
  $u\in\mathcal O_{K,S}^\times$ and $h(a)\le\rho\log H$; and $\Delta\le H^{-1/4-\varepsilon}$. So on this class the
  essential exponent is $\ge4$, the value Question B needs.
  * The proof is an Evertse/Corvaja–Zannier-type cancellation estimate from the Subspace Theorem with $K$ and $S$
    **fixed by hypothesis** (not tied to the parent, so the standing caution on fixed-$S$ theorems is respected),
    a formal Laurent model and a logarithmic Jacobian.
  * Scope. It covers fixed-order Pell, trace and golden templates with honest blocks (a bounded number of template
    blocks, pairwise coprime up to $H^{o(1)}$). A large raw Gaussian gcd is covered only if the reduced blocks are
    still bounded $S$-unit sums. The golden 8-point family sits at exactly exponent 4 on a curve
    ($\Delta H^{1/4}=3.70245917364$), consistent with the theorem.
  * It is **false or meaningless** for varying fields or orders ($v_n=(n+\sqrt{n^2-4})/2$ satisfies
    $v_{n+1}+v_{n+1}^{-1}-v_n-v_n^{-1}-1=0$ with no vanishing subsum), for unboundedly many summands and for $S$
    growing with $H$; a degree-uniform version would cover every Gaussian integer. Generic Gaussian blocks, which
    are all but an $X^{o(1)}$ proportion, are untouched.
  * Open, narrow: resonant templates (unbounded products of $a\lambda^r+b\lambda'^r$ in one fixed quadratic order).
* *Closed (round 31, Part VI.3).* General-type arguments on the cofactor cube-class cover
  $w^3=\prod_{a<b}(y_a-y_b)$: it is the triple cover of $\overline M_{0,7}$ branched on $B_2\in|-3K|$, with
  $\kappa=0$; it is not uniruled, and rational carriers of bounded degree lie in a closed set depending on the
  degree. A moving affine pencil of blocks $G_e=a_eB+c_e$ is rigid (tuples with $\Delta\le P^{-6+\nu}$, $\nu<1$,
  $P=H^{1/20}$, lie in a fixed closed set of dimension $\le2$; effective; a codimension-18 toy class).

> **Question B, literature form.** Is the identity $e$ **essentially bounded** on the smooth weak Fano 4-fold
> $X=\overline L{}^{\rm tw}_5$ (so that $X(\mathbb Q)=T_5(\mathbb Q)$ for the anisotropic norm-one torus $T_5$), in the sense
> of Lehmann–McKinnon–Satriano (arXiv 2403.02480, Def 1.3: $\alpha_{\rm ess}(e,-K_X+\varepsilon A)\ge\dim X$ for every
> ample $A$ and rational $\varepsilon>0$)? Since $-K_X$ is big, this is equivalent to $\alpha^{\rm ess}_e(-K)\ge4$
> (inference). The same constant is obtained on $X'$, since $\phi$ is an isomorphism near $e$. Less ambitiously: is
> $\alpha^{\rm ess}_e(-K)>15/4$, at least along sequences whose 15 toric heights are asymptotically equal? A yes gives
> the uniform Cilleruelo–Granville bound (gauss Thm 3.4'). Vojta's Main Conjecture implies it. McKinnon–Satriano's
> erratum (correcting their TAMS 374 (2021) Prop 1.3) has Prop 0.1: "Let X be a smooth projective variety over a number
> field k, and let P ∈ X(k) … If −K_X is big, then Vojta's Main Conjecture implies that X is canonically bounded at
> P." Lehmann–McKinnon–Satriano Prop 1.4 gives the same for essential boundedness on any smooth projective $X$.
> Neither assumes splitness. Both use Vojta on $X$ with a normal-crossings divisor of multiplicity $\dim X$ through
> $P$, a different instance from $V_E(5)$ with the same conclusion. In Lehmann–McKinnon–Satriano Thm 1.5 essential
> boundedness is a hypothesis. No unconditional proof was found for any variety where $\tau_P(-K)<\dim X$. Full texts
> were checked for seven papers; McKinnon–Satriano (TAMS) and Monahan–Satriano remain unread (`papers_check/README.md`).

> **Split calibration** (round 23 `split.md`, refereed). The same window occurs on the *split* Losev–Manin 4-fold at
> $e$. There a primitive positive 5-tuple $u$ has $H_{-K}=(\max u/\min u)\,\mathrm{lcm}(u)$ and
> $\mathrm{dist}\asymp(\max u-\min u)/\min u$. Split Question B asks whether five distinct integers with gcd 1 in
> $[U,U(1+\mathrm{lcm}^{-(1+2\delta)/4})]$, $\delta<1/30$, have their ratio class in a fixed proper closed set once the
> lcm is large. Equivalently, it concerns five factorisations of $N$ within hyperbolic angle $N^{-(1+2\delta)/4}$ on
> $xy=N$. It is **not known**: it implies Erdős Problem #887 in its absolute form ($\exists K\,\forall C$: eventually at
> most $K$ divisors of $n$ in $[X,X(1+Cn^{-1/4})]$, for every $X$), and Ruzsa's #886 for
> $\varepsilon>(15/64)(1+2\delta)$. Both are open (formal-conjectures, October 2026). The proof uses a random 5-subset,
> averaging over threshold layers, and Schwartz–Zippel. The 2026 five-divisor construction for #887
> ($n=(x^2-1)(x^2-4)(x^2-16)(x^2-49)$, $z^2=56x^2+65$) is a rational curve hitting each of the 15 cuts once, with
> $\alpha=15/4$ exactly. It lies in a 2-parameter family sweeping a hypersurface $W$, which every admissible split
> exceptional set must contain. Such curves do not exist on the twisted form: for a curve over $\mathbb Q$ there,
> $E_A\cdot C=E_{A^c}\cdot C$, so a balanced twisted curve at $15/4$ has $B\cdot C\ge30$ and $\mathrm{mult}_e\ge8$. A
> proof of Question B valid for every $\mathbb Q$-form of the 4-fold would therefore solve #887; a proof for $T_5$ alone
> may need the anisotropy. Split Question B is undecided. It would fail if infinitely many such hypersurfaces existed,
> or if absolute #887 failed.

> **Anisotropy** (round 24 `theory.md`, `search.md`, refereed). *Pointwise, it does nothing.* Let $\mathcal L$ be the
> arguments "archimedean Taylor upper bound $\ge$ product-formula lower bound" whose lower bound may use integrality,
> toric forced divisibility at every place of $\mathbb Q(i)$, point parity $v_{\bar\pi}=-v_\pi$ and reality of values,
> with auxiliary objects that are isotypic sections for $\sigma_X=\iota\circ\overline{\phantom{x}}$, polynomials in
> group words of A-points (including $\iota y$), or polynomials in the Hilbert-90 vectors ($x_j=P_j/\bar P_j$). Every
> such inequality is implied by toric Liouville for one section (or its $n$-point form at the generating points), so
> on balanced simple profiles with disjoint primes nothing in $\mathcal L$ beats $15/4$ ($\delta=1/30$), and nothing
> beats $\delta=1/10$ for arbitrary profiles, for every family with $\log H(F)+\sum d_i\log d_i=o(\sum d_iw_i)$. The
> reason is that $\sigma_X$ fixes $e$ and orders at $e$, so Galois averaging is lossless; the isotypic parts of
> $H^0(m(-K))$ have maximal orders $\lfloor10m/3\rfloor$ and $\lfloor10m/3\rfloor-1$ in both forms. Twist-sensitive
> congruence divisibility at inert places contributes $O_F(1)$. *Globally, it is visible only through curves.* By
> the parity lemma a balanced twisted curve at $15/4$ has $B\cdot C\ge30$. Such curves exist over $\mathbb R$: an
> explicit one-parameter real family $M$ (two conjugate branches of contact 4, each ray hit once; Jacobian rank 28,
> so a component). Only finitely many curves of its complex component are defined over $\mathbb Q$ (geometric
> monodromy $S_4$ of the place configuration, Faltings; ineffective; none of height $\le300$). The same finiteness
> holds for the split form, so this is not an anisotropy effect. *Sharpened in round 31 (refereed,
> `round31/ref_local.md`): no nondegenerate ray-generic member of the known component is defined over $\mathbb Q$,
> in either form; each form is excluded at 2 and at 3.* Whether the class has other components (Q24.1) is open.
> Blind numerical solvers all missed $M$, so failed searches in this class are not evidence of emptiness. Consequently
> a proof of Question B for $T_5$ that is not also a proof for the split form (hence of Erdős #887) must use a global
> twist-sensitive input; the only one known is the parity lemma.
>
> The disjoint-primes hypothesis is essential (round 25 `intra.md`, `referee_intra.md`, **proved**). The
> $\binom M5$ gcd-reduced 5-subtuples $y_S$ of one cluster are the projections $\mathrm{pr}_S(Z)$ of a single point
> $Z\in T_M(\mathbb Q)$, $T_M=(R^1_{\mathbb Q(i)/\mathbb Q}\mathbb G_m)^M/\mathrm{diag}$. Every multi-point
> Liouville/Dyson comparison among them, with toric forced divisibility across shared primes, is one-point toric
> Liouville on $T_M$ at $Z$. Conversely, a star partition $T_M\cong T_5^{(M-1)/4}$ identifies route 5's
> forced-divisibility forms with intra-cluster forms. Two results follow.
> * If $Z$ is balanced on all $2^{M-1}-1$ cuts of $[M]$, nothing in this class with a Taylor-type archimedean bound
>   beats $15/4$, uniformly in the auxiliary form.
> * For 5-wise uniform cut designs (e.g. extended BCH, $M\ge8189$) every 5-subtuple is exactly balanced with
>   $\mu<4(1-\eta)$, yet route 5's LP supercriticality forces $F(y)=0$.
>
> So shared primes do beat the one-point cap, but only by producing exact algebraic relations. Using those relations
> requires route 5's open zero estimate at the actual point. This gives nothing towards Question B.

> **Absolute form of the A-branch** (round 23 judge J4; independent referee `round23/referee_J4.md`, **proved**). If
> $A_5(\delta)$ holds for one $\delta\in(0,1/30)$, there is $K=K(\delta,Z)$ such that for every $C>0$ and every
> $R\ge R_1(C)=\max(C^{2/\kappa},(CR_0)^2)$, every arc of length $\le C\sqrt R$ holds at most $K$ lattice points; so
> $\sup_C M_\infty(C)<\infty$ (here $\kappa=1-(15/16+\eta)(1+2\delta)>0$, $\eta=(1/(1+2\delta)-15/16)/2$). The same
> argument gives it from $A_k(\delta)$, $\delta<1/(2^k-2)$, for every $k\ge3$: convexity averaging of the tuple gcd alone
> reaches the sharp range, with no balanced selection and no unit-class restriction. Scope: this is proved from the
> $A_5$ branch only, not from $(L_5^\circ)$ (whose route subdivides the arc and gives $K$ linear in $C$). So a
> construction with $M_\infty(C)\to\infty$ would refute $A_5(\delta)$ for all $\delta<1/30$, and $V_E(5)$, but not
> $(L_5^\circ)$. The absolute form is only formally stronger than $M(C)<\infty$, and it also follows from
> Cilleruelo–Granville's own $R^{1-\varepsilon}$ conjecture ($\varepsilon=1/4$).

> **Question B.** Prove $(L_5^\circ)$, or $A_5(\delta)$ for some $\delta<1/30$. Or refute $(L_5^\circ)$ by a
> Zariski-dense supply of balanced bounded-residue 5-tuples at exponent $\ge7/15$ (for a rational family, one
> dominating the 4-dimensional $T_5$); that would refute $V_E(5)$. Unlike $k=4$, a 2-parameter family plus its
> curve parameter has image of dimension $\le3$ and refutes nothing. By the killer's dimension count, one-variable
> polynomial-block families cannot be dense at $k=5$ (**heuristic**), and bounded-degree covering families are
> excluded (Thm 4.3(ii)) (`kill_gauss.md` §§4, 6.3).
> The two-branch twisted Pell ansatz ($B\cdot C=30$, contact 4, round 24) yields real curves at exactly $15/4$ but
> none over $\mathbb Q$ on the known component (round 31; other components open); it can only enlarge $Z$.

**Round 27 additions to Question B** (`round27/README.md`; all **proved** or **exact** unless labelled).
* *Split exceptional set.* Every admissible split exceptional set contains, besides $W$ and its image under the
  inversion $u_a\mapsto\mathrm{lcm}(u)/u_a$, the 15 quadric cones $(u_2-u_1)(u_3-u_1)=(u_4-u_1)(u_5-u_1)$ (and their
  $S_5$-images), swept by the second 2026 five-divisor construction ($C=7/2$, $\alpha=7/2$), and the 15 cubic threefolds
  obtained from them by inversion. Only $Q\cap\sigma Q$ can carry rational points in the twisted form.
* *A 2-dimensional calibration.* At an explicit rational point $P$ of an explicit degree-2 del Pezzo surface with
  $\mathrm{Pic}=\mathbb ZK$: $\tau_P(-K)=3/2<2$, the $\tau$-curve has a real cubic branch field and $\alpha=2$, and no
  $\mathbb Q$-curve of degree $\le6$ through $P$ has $\alpha<2$. So "$\alpha^{\rm ess}>\tau$ with $\tau<\dim$" is open
  already there. The surface has no balance, so it does not model $(L_5^\circ)$; and the Losev–Manin tower has no
  lower-dimensional balanced analogue, since $\alpha_{\rm bal}(k)=(2^{k-1}-1)/2^{k-3}<k-1$ iff $k\ge5$ (judge check,
  **exact**, unrefereed).
* *Pointer* (literature seen through search results only). $\mathrm{Bl}_eLM_5$ is a contraction of $\overline M_{0,7}$
  (Castravet–Tevelev, Duke 2015), so $\tau_e(-K)$, Theorem B and the section caps are statements about
  $\mathrm{Eff}(\overline M_{0,7})$; Question A lives on $\overline M_{0,8}$.

**Compared with $(L_8)$** (`round20/README.md` §§5, 7). $(L_5^\circ)$ is smaller: a 4-fold and 5 points, against a
5-fold and 8 points. It needs the archimedean divisor $D=E$ rather than $D=0$, and its heuristic margin is much
thinner: $\kappa=2\delta\approx1/15$ and $e^{-w}$ per scale, against $K\cdot\beta=10$ and $e^{-10w}$. Neither
statement is known to imply the other. Fakes violate $(L_5^\circ)$, whereas $(L_8)$ is fake-blind.

## 4. Question C (smaller): positive-gain non-unit pencils on $\overline M_{0,8}$

**Gain** (`round18/theory.md` §2.4; `round19/theory.md` §1). Let $X=\overline M_{0,8}$, let $\partial=\sum_TD_T$
(119 components, snc), and let $\beta$ be the balanced class. For a pencil (a non-constant rational function) $f$,
write $F_p=f^*(p)\sim L$ and $R_p=(F_p)_{\rm red}\cdot\beta$. Define
$$\mathrm{gain}(f)=L\cdot\beta-(R_0+R_1+R_\infty)=\sum_{Y\subset F_0+F_1+F_\infty}\tfrac{m_Y-3}{3}\,Y\cdot\beta .$$
When the base locus lies in codimension-2 boundary strata, this equals $10-K^{\log}_{\tilde X/\mathbb P^1}\cdot\tilde\beta$.

**The abc criterion** (`round18/theory.md` Prop 2.6; **conditional** on abc; implication **proved**).
* A *unit pencil* is $f=\lambda A/C$ with $\lambda\in\mathbb Q^*$, where $A,C$ are coprime products of differences
  $x_i-x_j$ in Kapranov's chart. Put $B=C-\lambda A$.
* Claim: if $\mathrm{gain}(f)>0$, then $(L_8^{pc})$ holds with $Z=\{B=0\}$. As first printed this conclusion is
  vacuous ($\mathrm{Bal}^{pc}=\emptyset$, round 31, Part IX.2); whether the criterion gives the corrected
  $(L_8^{pc})'$, where the gcd identity should hold up to a factor $B_0^{O(1)}$, is **unchecked**.
* Mechanism: on $\mathrm{Bal}^{pc}$ we have $\gcd(A(x),C(x))=\prod d_T^{\mu_T}$, so
  $h(f(\Phi))=(L\cdot\beta)w+O(\eta w)$. Divide $(\lambda A(x),B(x),C(x))$ by its gcd to get a coprime abc triple.
  The log of the radical of that triple is $\le(F_{\rm red}\cdot\beta)w+O(\eta w)$.
* Scope: the refereed abc barrier covers only relations among products of differences (`round8/weakest.md`
  Prop 4.2), so this was a genuine gap in the barrier map.

**Round 19 closed the unit part** (`round19/theory.md`, refereed in `round19/referee_theory.md`; summary in
`round19/README.md` §1).

> **Theorem 2.4** (**proved**; $\beta$-weighted log Riemann–Hurwitz). For every non-constant rational function $f$ on
> $\overline M_{0,8}$,
> $$\sum_{p\in\mathbb P^1}\big(L\cdot\beta-R_p^{\rm int}\big)\ \le\ 2L\cdot\beta+\Lambda_{k(f)}\cdot\beta\ \le\ 2L\cdot\beta+31.$$
> * $R_p^{\rm int}$ is the $\beta$-degree of the reduced non-boundary part of $F_p$.
> * $k(f)$ is the least $|S|$ such that $f$ factors through a forgetful map $\pi_S$.
> * $\Lambda_k=\pi_k^*\psi_k$, with $\Lambda_k\cdot\beta=16,24,28,30,31$ for $k=4,\dots,8$.
>
> Equality holds for the cross-ratio and for $4z(1-z)$.
>
> **Corollary 2.5.** $\mathrm{gain}(f)\le31-n'(f)$, where $n'(f)=\#\{T: f|_{D_T}\text{ constant}\}$.

*Proof idea.* The form $df$ spans a line subsheaf of $\Omega^1_X(\log\partial)$, and its divisor is read off the
fibres. The forgetful filtration
$F_k=\pi_k^*\Omega^1_{\overline M_{0,k}}(\log\partial_k)$ has line-bundle quotients $\Lambda_k$. Every line
subsheaf of $F_k$ that is not in $F_{k-1}$ has $\beta$-degree $\le\Lambda_k\cdot\beta$, because $\beta$ is positive on
every effective divisor (Keel–McKernan, via `threshold/L8.md` Thm A(a)). This is not a curve argument. Restricting to covering curves of
class $d\beta$ only gives gain $\le(2g-2)/d$, a bound that is $\ge10$ and so cannot exclude positive gain.

> **Theorem 3.3, Corollary 3.4** (**proved**). For every non-constant unit $u$ of $M_{0,8}$ and every
> $\lambda\in\mathbb C^*$,
> $$\mathrm{gain}(\lambda u)\le\Lambda_{k(u)}\cdot\beta-N(u)-n_1\le-2^{8-k(u)}\le-1\qquad(\le-16\ \text{if}\ k(u)=4),$$
> where $n_p=\#\{T:D_T\subset\operatorname{supp}F_p\}$ counts the boundary divisors in the fibre over $p$.
> The same bound $\mathrm{gain}\le-1$ holds for every composite $h\circ u$.
>
> The key input is $N(u)=\#\{T:\mathrm{ord}_{D_T}u\ne0\}\ge32$, from multilinear Schwartz–Zippel on $\{0,1\}^7$.

So Prop 2.6 never applies. Exact computation agrees: gain $\le-32$, attained, for all unit pencils of degree
$D\le6$ and all $\lambda\in\mathbb Q^*$ (`round19/search.md`, reproduced by the referee; `round18/theory.md` §2.4).
Morphisms have gain $\le-32$ (Bruno–Mella, **read from the abstract only**).

**What is open** (`round19/theory.md` §4; `round19/referee_theory.md` m2). The remaining case is non-unit pencils in
$$P^\circ=\{f\ \text{regular on}\ X\setminus\textstyle\bigcup_{T\ne T'}(D_T\cap D_T')\}.$$
Two things are needed, and both are open.
* **(i) A gcd extension of Prop 2.6** to $P^\circ$. This means bounding the excess of $\gcd(A(x),C(x))$ over
  $\prod d_T^{\mu_T}$ at balanced pairwise-coprime points. Round 18 asserted, without proof, that for base locus in
  codimension-2 boundary strata the excess gcd vanishes on $\mathrm{Bal}^{pc}$ (`round18/theory.md` §2.4, after
  Prop 2.6). Round 19 §4.1 records this as plausible but **unproved** and gives only an unrefereed sketch of the
  obstruction. (For base locus meeting the interior, outside $P^\circ$, round 18 calls the gcd "a Vojta-type gcd
  problem in its own right".)
* **(ii) A pencil with gain $\ge1$.** Any such pencil must satisfy these proved necessary conditions:
  * $n'(f)\le31-2^{8-k(f)}$ and $k(f)\ge5$;
  * $f$ is not of the form $h\circ u$ with $u$ a unit;
  * all three fibres are non-reduced, with $\beta$-average multiplicities $\mu_p=L\cdot\beta/R_p$ satisfying
    $\sum_{p\in\{0,1,\infty\}}1/\mu_p<1$;
  * $\sum_{p\notin\{0,1,\infty\}}(L\cdot\beta-R_p^{\rm int})\le30-n(f)$, where $n(f)=n_0+n_1+n_\infty$ counts the
    boundary divisors in the three special fibres ($n(f)\le n'(f)$).

  **No example is known**, and Theorem 2.4 alone allows gain $\le31$. Campana-type statements see only gcd
  multiplicities, so they do not decide the question.

Even with both (i) and (ii), the result would be only an abc-conditional $(L_8^{pc})$ (for the corrected form $(L_8^{pc})'$ this is unchecked; see the claim above). Reaching $(L_8)$ also
needs a nested-chain extension. There, for the chain $\{1,2\}\subset\{1,2,3\}$, the gcd exceeds Prop 2.6's value
by a full block (`round18/README.md` §3.3, J2; **unrefereed**). The judge adds that the radical also counts $p$
once, which might compensate, and that this is not analysed; so J2 shows a gcd excess, not that the extension fails.

> **Question C.**
> * (a) Does some $f\in P^\circ$ have $K^{\log}_{\tilde X/\mathbb P^1}\cdot\tilde\beta<10$, that is,
>   $\mathrm{gain}(f)>0$?
> * (b) Does Prop 2.6 extend to $P^\circ$? This would be a gcd inequality for two forms at balanced points.
> * (c) Is Theorem 2.4 known, or useful elsewhere? The record notes that it holds verbatim for any movable class
>   $\gamma$, with $31$ replaced by $\max_k\Lambda_k\cdot\gamma$, the maximum taken over orderings (unused there;
>   `round19/theory.md` §7 item 4, referee m7). Is it known, or sharp, in that generality, and on $\overline M_{0,n}$?

## 5. What any proof must avoid

Each item below is a barrier result in the record (`overview/overview.md` §§5.4, 7.6). Statuses are as given there.
1. **Local inputs.** These are pairs, gcds, cuts, monomial (Liouville) gaps, rational characters, and point-level
   residues at moduli $\le M^A$ with $A$ fixed in advance. Fakes at $W=O(M\log M)$ satisfy all of them (**proved**,
   with the scope limits of `overview/overview.md` §5.2 item 2).
   Caveat: caps uniform in $A$ are proved for zero-sum characters and for factored prime residues tested through
   characters (round 31, Part IX.3), not for twisted collisions or the 2-adic clause; P3-all$^*$ is uncapped only on the balanced words inside
   label classes, in particular data-chosen label differences (round 31, Part X.2), and the residue route is closed
   for internal work (`round31/README_X.md` §3(b)).
2. **Averaged or robust analytic inputs.** This covers large sieves with non-sharp constants, prime sums or box
   counts with power-saving error, and "almost all" statements. Grafted fakes satisfy them (`round10/analytic.md`
   Thm B, **proved**; composite-supported version `round20/composite.md`, killer-checked). Inside the low-height window the
   earlier exceptions (`round21/README.md`) stand only for $\delta\ge1/6$ or for the inputs listed as uncapped in
   `round31/README_IX.md` §3(d). The single-label max-scale cross-gap at constant one is capped outside the window
   (round 31, Part IX.5).
3. **Auxiliary polynomials and single-place Diophantine approximation.**
   * On $X_M$: $\tau^*(M)=2-1/\lceil M/2\rceil<2$, so Liouville at $\infty$, and Ru–Vojta with $S=\{\infty\}$, never
     reach exponent $\frac12$ (`threshold/upper.md`, **proved**).
   * On $\overline M_{0,k}$ for any $k$: $\beta$ lies in the interior of the movable cone (`threshold/L8.md` Thm A,
     **proved** as a formal barrier).
   * On the Gaussian model: the forms that yield an inequality at the actual point are exactly critical (Lemma K1
     with `threshold/construct.md` Cor C3 and Thm B). Grid determinants with free nonvanishing are all subcritical;
     a supercritical design needs an irreducible factor of degree $\ge2$, whose nonvanishing is route 5's open zero
     estimate (`round31/README_II.md` §1.3).
   * Roth- and Schmidt-type bounds (McKinnon–Roth) give $\alpha^{\rm ess}\ge\beta<\tau$ at a rational target: weaker
     than Liouville.
   * Multi-point (Thue/Siegel/Dyson) arguments on the 5-subtuples of one cluster are one-point toric Liouville on
     $T_M$ at the cluster point, which is the determinant method with forced divisibility ("route 5"). In the growth
     regime they force exact algebraic relations and stop at a zero estimate at the actual point. With lacunary degrees
     they reduce to one-point caps (round 25 `intra.md`, **proved**; the comparable-degree dichotomy is **heuristic**).
   * *Zero estimates do not supply route 5's missing lemma* (round 26 `ze.md`, `referee_ze.md`; **proved** at the
     degree level, given the shape of Philippon's 1986 zero lemma, with margins $\ge e^{1.7\cdot10^4}$ that make its
     exact constants irrelevant). Route 5's Liouville step certifies, for one form $P_c$, vanishing at $Z^M$ only for
     integer endomorphisms $M$ with $\Lambda_M(c)<\mathrm{ord}(c)-O(E)$; for the powers $Z^j$ this means
     $|j|<\mathrm{ord}(c)/\Lambda(c)\le1.07\sqrt n/\lambda$, by the spectral cap
     $\mathrm{ord}(c)\le(\sqrt n/\lambda)\Lambda(\mathrm{supp}\,c)$, and only $|j|\le1$–$4$ on regime boxes. A
     Philippon estimate on $\mathbb G_m^d$ needs progressions of length about $d$, so its conclusion is met by the
     trivial subgroup. The dual framing (one exponential polynomial on $\mathbb G_a^d\times\mathbb G_m$ vanishing on
     the slab box) has obstruction subgroups whose counts do not depend on the point, is vacuous at the degree level,
     and is false in Hilbert-sharp form. Moreover it concerns $\mathrm{Bs}(V(A,m))$, while route 5 needs
     $Z\notin\mathrm{Bs}(V_H)$ for the height-cut forced space $V_H$. The global statement any group-theoretic
     argument would have to prove is therefore **ZE$^{\rm all}_H$: $\mathrm{Bs}(V_H)$ contains no nondegenerate point**,
     which is not group-structured. The Hilbert-level claim that no single-form estimate can help is proved only for
     powers-only data and well-conditioned profiles ($\lambda\ge\lambda^*(n,\mu)$, checked at four values of $n$).
     Profile automorphisms are certified too, and once they are counted the argument fails near $m_{\max}$ (the
     claim is then unproved, not refuted), so it should not be cited as general.
     The actual point lies on no proper algebraic subgroup, but that does not affect a count at the trivial subgroup.
   * *A linear-algebra sufficient condition for route 5's target* (round 27 `interp.md`, `kill_interp.md`; **proved**,
     new to the record; the hypothesis is **open**). Let $L_{\rm pair}=\{\pm e_x,\pm(e_x-e_y)\}$. If, for all $e<m$, the
     degree-$e$ least space of the slab box $A$ is the sum of those of its shift-intersections $A\cap(A-l)$,
     $l\in L_{\rm pair}$, then $\mathrm{Bs}\,V(A,m)$ lies in the braid-and-coordinate arrangement (round 11's ZE-braid);
     with a top-degree height condition this gives $M=O_C(\log\log R)$. The hypothesis is at least as hard as
     ZE$^{\rm all}_H$. Evidence on slab boxes is mixed: the criterion beats every product (zonotope) certificate on the
     box M8$'$, but stays below the zonotope order on 20 of 27 random boxes.
     *Round 28 (pre-registered test, `round28/README.md`): unfavourable.* On 105 independently verified slab boxes
     ($d=5$–$9$), the criterion's reach exceeds the zonotope order only on M8$'$ (reach 11, exact over $\mathbb Q$) and
     on one $d=8$ box containing M8$'$'s configuration; the true one-point extension of M8$'$ has reach $9<10$. The
     share does not grow with $d$. Reach, $\rho$ and the zonotope order add under products of boxes, so no bound
     reach $\le$ zonotope order $+\,c$ holds in general, but products violate the character gap. Under the fixed rule
     the interpolation branch of route 5 is closed for internal work, and the corollary $M=O(\log\log R)$ is treated
     as presumptively vacuous.
   * *No tight pencils* (round 27 `stab.md` Thm 3.3, Cor 3.4, `kill_stab.md`; **proved**, killer-checked). Mean-width
     extremal measures are unique up to scalar on their support ($\dim A_t(S)\le1$ for integer $t\ge MW_n(S)$, all $n$).
     Hence no pencil of tight sections exists on $\mathrm{Bl}_eX$ (both forms, every $k$), and "second-order integrality"
     at Liouville-tight points is void. This closes the two-section part of the equality-case loophole of Theorem B;
     the single-form part still rests on the structure half of `threshold/construct.md` Conjecture E.
   * *Integrality at the root arrangement adds nothing at $e$* (round 27 dp2 Thm C and transfer, split Lemma 3.5, the
     dp2 killer's line lemma; **proved**). Bounded pair cofactors make near-tight tuples integral with respect to the
     diagonals, but on $\mathrm{Bl}_e\overline L_k$ that pair has log-canonical degree $-1$ on the lines through $e$ for
     every $k\ge3$, so Lang–Vojta predicts nothing there, and poles along the diagonals amplify nothing beyond sections
     wherever the section cap holds. A genuine fixed-$S$ Subspace argument on $LM_5\setminus\bigcup\Delta_{ab}$ is not
     covered.
4. **Single linear forms in logarithms** (Baker–Matveev, single-form Lang–Waldschmidt). Fakes satisfy them
   (`round8/lw.md` Thms A, B; `round8/stability.md` Prop 6.1; **proved**). A uniform two-logarithm inequality
   $(1+w(v))\log(2D_L(v))\ge c_*W$ for coordinates $v$ of finite order $D_L(v)$ modulo the exponent row lattice
   holds on actual clusters and fails in the binary function-field family (round 31, Part VIII.2; **proved**, low
   novelty, cf. `growth/auxiliary.md` Prop B2). It gives no growth bound and has not been tested against fakes.
5. **$S$-unit and Subspace arguments with thresholds depending on $S$ only through $|S|$.** With an absolute height
   threshold, polynomial or even sub-exponential exceptional-subspace counts are false in every dimension: covering
   needs $\ge\exp(s^{1-\delta/(m-1)-o(1)})$ subspaces, also on exactly equal norms in dimensions $2\le k\le6$
   (round 31, Part VIII.4 of the Question B thread; **proved**, refereed; earlier lit18 Prop A,
   $\exp(s^{1/2-o(1)})$, unrefereed). So round 8's hypothesis QS at a fixed threshold is **false**. With thresholds
   depending on $|S|$ only, lit18 Prop B (**unrefereed**) still caps these methods at $M\le\omega(N')$. A threshold
   relative to $\prod_{p\in P}p$ (QS$^{\rm prod}$) is open and non-standard. On one critical arc it would suffice
   with $f(s)=o(s^2)$ subspaces, because a hyperplane contains at most $\max(2,1+C^2)M^{k-2}$ ordered $k$-tuples of
   the arc (Part VIII.5; **proved**). A polynomial count, with exponent independent of $\delta$, of norm-one
   $P$-units $u$ with $|u-1|^2\le H(u)^{-\delta}$ above the relative threshold would give $M\ll(\log R)^\tau$ for
   every $\tau>0$ (Part VIII.6; **proved** implication); with exponent $0$ that count already implies the uniform
   bound, and no standard conjecture is known to imply it. Uniform unit-equation counts with fixed anchors are closed:
   a rank-$r_0$ multiplicative coset contains the divisor splits of at most $8r_0+4$ third vertices of a fixed edge,
   so the rank needed at a growth failure is $\gg M$; the obstruction holds in angle fakes
   (`round31/ref3_rank.md`; **proved**). Fixed-$S$ theorems such as `round18/theory.md` Thm 2.5 have
   $S$-dependent constants, and $S$ is the unbounded set of primes of $N$. The same holds for the fixed-field
   unit-sum theorem of round 31, Part VI.1 (§3): its degree-uniform or growing-$S$ version is false.
   Since $M\le1+\omega(N')+K_C$ is now proved unconditionally (Theorem 1.2), these methods add nothing unless they
   beat $\omega(N')$.
6. **abc on single relations.** abc through relations among products of chords (pair, triangle, Ptolemy or
   $n$-term) is vacuous on squarefree profiles (`round8/weakest.md` Prop 4.2), which contain every recorded fake
   (Cor 4.3), so it cannot prove growth (**proved**). abc on unit pencils is closed as well (Theorem 3.3).
   abc on character relations (products of points, not of chords) is equivalent to the Lang–Waldschmidt support
   test (round 31, Part VI, `round31/ref6_code.md`): it excludes the fully flipped families with $b\ge4$ and not
   the skew double core. The stronger parity condition that holds over $\mathbb F_2[t]$ (item 12, Part VI) has no
   integer analogue: $z=\gamma^2(2+i)$ with $\gamma=x+iy$, $(x+2y)^2-5y^2=\pm1$, has $|z-\bar z|=2$ and parity
   support $\{5\}$ while $\log N\to\infty$.
   *What is not capped* (round 31, Part IX of the Question B thread; `round31/README_IXL.md` §1.8). abc on the
   integer triangle $t_{12}A_{12}+t_{23}A_{23}=t_{13}A_{13}$ (three marks in one unit class: $z_j-z_i=2i\eta
   H_{ij}t_{ij}$ with $t_{ij}\in\mathbb Z\setminus\{0\}$ and $H_{ij}$ the common Gaussian part; the $A_{ij}$ are
   products of the parent's primes) is vacuous for the same reason. But the exact identity itself, with small
   integer cofactors, is **not capped**: it is the input the fakes lack, and it gives a rate for random parents
   (item 12, Part IX). For a fixed parent it is a squarefree three-term $S$-unit equation with $|S|=s$ and small
   coefficients, whose content is open; known counts are $\exp(O(s))$. Its local consequences at primes not dividing
   $N$, triangle areas and abc remain capped, and the earlier classification of the cocycle as capped
   (`round10/additive.md` §6; `overview/overview.md` §5.4 rows 15–16) holds only for those consequences.
7. **gcd normalisation.** It is invisible to forced-divisibility, auxiliary-form, Liouville, pair, cut, character
   and single-form inequalities (Lemma K1, proved by the killer). Vojta's $D=E$ term and the heuristic count do
   see it.
8. **Genus-0 Arakelov theory and unlikely intersections in tori.** The first reduces to items 1 and 3 (overview
   §5.4 rows 2 and 8) and to the limitation "Lim" of known tools (Mordell–Weil packing gives $e^{O(\mathrm{rank})}$; bounded rank is equivalent
   to the uniform theorem). The second is vacuous, or capped at $M\le\omega(N')+O_C(1)$ (a bound now proved unconditionally, Theorem 1.2); the Rémond-cone case of
   that ceiling is only **heuristic** (`round20/arakelov.md`, `round20/wild.md`; killer-checked;
   `overview/overview.md` §5.4 rows 18–19).
9. **For refutations.**
   * No feasible computation exists (`round18/census.md` Thm E).
   * No density mechanism works: proved at $k=8$, heuristic at $k=5$ for circles.
   * Bounded-degree covering curve families whose members each carry infinitely many balanced small-$H$ points
     are excluded at $k=8$ (genus $\le1$ by `round18/census.md` Prop 7.5; genus $\ge2$ by Faltings), and, for
     circles, those whose members each carry infinitely many A-points are excluded at $k\ge5$, $\delta>0$
     (Thm 4.3(ii)). Families with infinitely many rational points per member are not excluded as such (the lines of
     Kapranov's $\mathbb P^5$ are one; correction M1).
   * What remains are dense supplies with finitely many points per curve, curves of unbounded degree, and sporadic
     points.
10. **Round 29 (four unconventional lanes, `round29/README.md`; killer-checked).** None passes both filters.
   * *Function-field model* (`round29/ff.md`, `kill_ff.md`). Over $\mathbb C(t)$, twisted by conjugation, the arc
     condition is $L\approx D/2$ exact real-linear equations on exponent vectors (already `families.md` Lemma P). The
     model has its own fakes (arbitrary columns in $\mathbb R^L$, linear growth), covering families satisfy
     $r(L+1)+1\le 2D$ (the record's free-curve bound). Bounded clusters are equivalent to bounded rank of the ratio
     lattice; the uniform question was open there until round 31, Part IV, which proves $M(D,0)\le24$ for every $D$
     and $\sup_DM(D,c)<\infty$ for every $c$ over $\mathbb C(t)$, $\mathbb Q(i)(t)$ (item 12). New and exact: the *cube family*, 8-point clusters at $c=0$ for
     every $D\equiv0\ (8)$ over $\mathbb Q(i)(t)$, a non-covering rational family, with an 8-point template in
     $\mathbb Z[i][u]$ (constant $\approx3.2\cdot10^9$). The transfer to $\mathbb Z$ is negative: the $D/2$
     independent equations come from the infinite residue field; a $p$-adic arc gives one congruence. Over the fixed
     finite fields $\mathbb F_2$ (round 31, Part V) and $\mathbb F_3$ with the form $A^2+B^2$ (Part VI) the model is
     nevertheless unbounded (item 12), by additive (linearised-polynomial) mechanisms with no integer
     counterpart.
   * *Operator templates in the style of Li–Liu's conditional Goldbach paper* (`round29/spectral.md`). Their
     Hypotheses S and R are operator-norm bounds equivalent to what the conclusion needs. Characters mod $N$ vanish on
     the circle (every lattice point divides $N$ in $\mathbb Z[i]$); Hecke-angle versions are fake-blind; the
     neighbouring-circle pair count is exactly a sum of short-arc counts of the same problem at smaller radii, so any
     uniform bound on it is a restatement.
   * *Residue-class twin* (`round29/clothes.md`, `kill_clothes.md`). Points of the circle in one class mod
     $q\ge\lambda N^{1/4}$ obey the same $2/3$ bound for every odd $q$ coprime to $N$; the lane's claim that smooth
     $q$ escape it is **false** (killer's Lemma K). No implication between twin and arc problem is known.
11. **Round 30 (the function-field model, three lanes, `round30/README.md`; killer-checked).** The line is
   **closed for internal work**; nothing changes for the integer problem.
   * *Rank.* No rank-5 cluster at $c=0$ was found for $D\le12$ (search evidence, not proof); rank 5 does occur at
     $c=1$ ($D=10,12$, certified). Bounded rank (R$_4$) stays **open**. It has exactly the shape of truncated
     function-field Vojta for the pair $(\mathrm{Bl}_eX,\pi^*\partial+E)$, $r(L+1)\le|S|-1+\varepsilon h$, which is
     itself open, so this is calibration and not evidence. *This concerns only the sharp constant.* A uniform bound
     is a theorem (round 31, Part IV, item 12): $8\le\sup_DM(D,0)\le24$.
   * *The one new tool: the saturated balance lemma* (`round30/search.md` Lemma 1.1 = `prove.md` Lemma 1.2;
     proved twice). Every nonzero integer vector $a$ in the real solution space of the moment equations has
     $\|a\|_1\ge L+1$ and $\min(a^+,a^-)\ge L/2$. It kills the round-29 function-field fakes, and an exhaustive
     search with it gives $M(7,0)\le5$, $M(8,0)=8$, $M(9,0)\le5$ (**judge's computation**; $M(9,0)$ remains
     unrefereed; $M(7,0)=4$ exactly by round 31, Part VII, `round31/ref7L_search.md`). It does not bound the rank by itself: 21 classes of 6-point rank-5 codes survive it at
     $D=8$. *Round 31 (refereed, `round31/ref_ff.md`):* all 21 classes are excluded by non-linear certificates, so
     at $D=8$ the rank is $\le4$ and $M(8,0)=8$ with a unique extremal (cube) code, for squarefree $N$ over
     $\mathbb C(t)$; there are at most 12 constant-offset divisors up to sign, attained only by the cube
     (`round31/ref2_min.md`), and the lattice minimum can exceed $D/2$. Coordinator's remark (unrefereed): the weight half is the function-field form of the monomial-gap
     inequality of §1, which the integer fakes already satisfy; the balance half uses the real place and has no
     integer counterpart.
   * *Bridge to $\mathbb Z$* (`round30/bridge.md`). A primitive $\mathbb Q(i)(t)$-cluster with $D$ even, $c=0$
     specialises to integer clusters with an explicit limiting constant (proved). The best cube-family constant is
     $C_\infty=145877.66$, far above the golden family's $5.89771$, so no lower bound for $M_\infty(C)$ improves.
     All 56 five-subsets of a cube cluster have $\alpha=4$: no A-points, and they avoid the round-27 quadric cones
     and cubic threefolds (exact). Question B is unchanged.
   * *Withdrawn in the round:* "vertical-line loci have rank $\le3$" (certified rank-4 counterexample at $D=10$);
     "balance does not separate real clusters from fakes" (true only for pairwise balance); the "real-algebraic
     obstruction" lead (it was linear algebra, Prop K1).

12. **Round 31 (ten exchanges with an outside colleague, "Astra"; 108 referee blocks and seventeen judge reports:
   `round31/README.md`, `README_II.md`, `README_III.md`, `README_IV.md`, `README_IVL.md`, `README_V.md`,
   `README_VL.md`, `README_VI.md`, `README_VIL.md`, `README_VII.md`, `README_VIIL.md`, `README_VIII.md`, `README_VIIIL.md`, `README_IX.md`, `README_IXL.md`, `README_X.md`, `README_XL.md`).** No referee found a false theorem in the colleague's reports; one case missing from a proof
   (Part V, unit-twisted points) was repaired by the referee. **The growth rate for general $R$, the uniform bound,
   Questions A and B and route 5 keep their status** (judges' scores on the goal: 2/10, 2/10, 1/10, 2/10, 1/10, 2/10,
   1/10, 2/10, 2/10, 1.5/10, 1.5/10, 3/10, 2/10, 2.5/10, 2.5/10, 2/10, 2/10). Parts VIII and IX are about integer circles: Part VIII adds Theorem 1.2 (§1)
   and the first fakes in the low-height window, and on the Question B thread shows that Question A contains
   Rudin's conjecture (§2); Part IX adds Theorem 1.3 (§1), caps character-tested prime residues at every modulus
   and extends the window cap to $\delta<1/6$; on the Question B thread it gives a rate for random parents at the
   smooth critical height and corrects the lattice fallback of Question A (§2).
   The other new theorems are in the function-field model (Parts IV, V and VI below): the uniform statement is true in
   characteristic 0, false over $\mathbb F_2(t)$ (also for squarefree parents) and false over $\mathbb F_3(t)$
   with the Gaussian form $A^2+B^2$. What is now in the record, all
   **proved** and refereed unless labelled:
   * *Question B is the index-one case.* The valuation columns of a five-tuple span a lattice of index $I$.
     If $I\ge5$ the tuple lifts to a degree-$I$ isogeny cover and McKinnon–Roth gives exponent $>10/3$; if $I\ge7$
     it gives $A_5(1/31)$ off one proper closed set, for heights $H_S\le C^{-94/25}N^{47/50}$. This is irrelevant
     to clusters: in a shortened genuine one-unit cluster all but $O(m^3)$ ordered five-tuples have index one
     (absolute constant, via Baker–Wüstholz, confirmed through secondary sources only), and a subcluster of
     $\ge m/145-O(1)$ points has every five-tuple of index one. Hence **$(L_5^\circ)$ restricted to index-one
     tuples suffices for the uniform bound**. It is formally weaker and no easier: every hard model has index one.
     Covers étale above the torus cannot lift index-one points, and the one cover that does (cubic) gains nothing.
   * *First exponent bounds that genuine tuples satisfy and gap-admissible fakes fail, on a thin family* (judge's
     reading of the fake side, unrefereed). On the ten-pair-cut profile, where the gap-only angle model attains
     exactly $10/3$: tuples whose ten block norms lie on a fixed rational plane $U$ have exponent $\ge38/11$
     ($40/11$ for comparable norms, $4$ on arithmetic progressions), with constants depending on $U$, ineffective
     in one stratum (Faltings on a genus-five curve); six norms in $[X,X+cX^{1/6}]$ give exponent $4$. Neither
     hypothesis is supplied by the growth problem, and $40/11<15/4$.
   * *Integrality with all divisibilities is not enough at five marks.* There are arbitrarily large parents with
     a primitive integer $2\times5$ star matrix, index one, all ten divisibilities with bounded cofactors, norms to
     relative error $O(N^{-1/10})$, and width $\asymp N^{-3/10}$ (exponent $10/3$). Exact individual norms together
     with the divisibilities force a genuine cluster. So the gap between an integral model and a circle is exactly
     the norm condition. No integral model with full divisibility is known at growing $m$ (half the log
     divisibility is the trivial CRT ceiling); **open**.
   * *Rigidity.* Under pair separation and for $m\ge16390$, the parent and the labelled separation divisors
     determine all rational star directions up to one element of $\mathrm{PGL}_2(\mathbb Q)$ (robust to deleting a
     tenth of the constraints for $m\gtrsim4\cdot10^4$). Scalar and five-vertex reconstruction prove a round-8
     converse. No count follows.
   * *Random parents (almost-all; not a worst-case statement).* Let $N$ be the product of $r$ distinct split
     primes drawn uniformly, in random order, from $[r^D,2r^D]$, and fix $C$. Outside a set of parents of relative
     size $\exp(-c(D,C)\,r\log r)$, clusters of normalised norm $N$ have at most $M_D$ points on any arc of length
     $C\sqrt R$: $M_D=1$ for $D>4$, $\le2$ for $8/3<D\le4$, $O(1)$ for $D>2$ (elementary pair first-moment count,
     the referee's), and $M_D\le37$ at $D=3/2$, finite for every $D>1$ (quartet exposure). This averages over a
     density-zero family. It says nothing about any individual $N$, about the narrow-bin parents on which fakes
     are built (they lie inside this exceptional set), or about the growth rate.
     *Extension (round 31, Part IX of the Question B thread; refereed, `round31/README_IXL.md` §§1.1–1.3;
     almost-all over a prime pool, not a worst-case statement).* Let $S$ be a uniformly random $s$-subset of any
     set of $K\ge e^{40}s$ primes in $[P,2P]$, $q=\prod_{p\in S}p$, $J=\log(K/s)$ and
     $b-1=\max(2,\lceil100\log(2P)/J\rceil)$. For $s\ge s_0(C)$, outside a fraction $\exp(-sJ/20)$ of parents, every
     arc of length $Cq^{1/4}$ of $x^2+y^2=q$ (pool of split primes) carries at most $4(b-1)$ points, and every
     interval of length $C\sqrt q$ contains at most $b-1$ integers lying in one of two prescribed classes modulo
     each $p\in S$, uniformly in the arc, the interval and the classes. At $P=As\log s$ ($A\ge2e^{40}$) this is
     $O(\log s)$ at height $\log q\sim s\log s$, the endpoint $D=1$ at which $M_D$ above is unbounded and where the
     worst-case bound is of order $s$; at $P=s^{1+\epsilon}$ it is $O_\epsilon(1)$. Exponents $1\le e_p\le E$ chosen
     after the sample are allowed, with constants polynomial in $E$. The proof uses three marks: some triple has
     small integer cofactors, its triangle identity determines a block of at least $s/5$ of the primes from the
     rest, and a union bound over $7^s$ labellings beats the pool entropy. The exceptional fraction is $R^{-o(1)}$
     at the smooth scale, the exceptional set may contain $(\#\text{parents})^{0.95}$ parents, and the dense
     profile extremal for Theorem 1.1 ($K=s$) is not covered. A first-moment heuristic predicts $O(1)$ already at
     $P=As\log s$. *Same norms for fakes and genuine arcs:* on typical $s$-subsets $S$ of a bin of split primes of
     log-width $l^{-3/4}$ in $[l^{1.99},2l^{1.99}]$ ($s=l-1$), the Legendre fakes of Part VIII exist on
     $N=\prod_{p\in S}p^2$ with $\asymp s/2$ marks while every genuine arc of length $C\sqrt R$ carries $O_C(1)$
     points (asymptotic only). So the inputs the fakes satisfy cannot certify $M=o(s)$ even where the truth is
     bounded; the separation is by the truth, not by a test the fakes fail.
   * *Closed routes.* Unit-equation counting with fixed anchors (item 5); grid determinants with free
     nonvanishing (item 3); fixed-cofactor elimination (the fixed-$t$ norm vectors form an irreducible fourfold of
     degree 12 that dominates ratio space); quartet cross-ratio consistency (automatic); template density (ten-pair
     template curves of fixed block degree sweep out dimension $\le2$ of the 4-dimensional ratio space); the
     completely multiplicative sharp sieve (§1); the near-square-divisor count $D_B(N)$ (fakes match it).
   * ***Part IV: a uniform bound in the function-field model, and a new filter*** (`round31/README_IV.md`; seven
     referee blocks; Langer's theorem read in the original, `round31/langer_check.md`).
     * *Inequality (IV.2)* (**proved**; characteristic 0). Let $N,n_i\in\mathbb C[t]$, $L=\deg N>0$, $n_i\mid N$,
       $\alpha_i=\deg n_i$, and let $P_1,\dots,P_m\in\mathbb C[t]^2$ ($m\ge3$) have coprime entries and nonzero minors
       $D_{ij}$. Count one layer $h$ per root of $N$ and depth, with $s_h$ the number of $n_i$ vanishing to that
       depth; put $E_{ij}=\deg\gcd(n_i,n_j)-\deg\gcd(n_i,n_j,D_{ij})$ and suppose
       $\deg D_{ij}\le(\alpha_i+\alpha_j)/2-L/4+\delta$. Then
       $\frac34(m-24)L+\sum_h(s_h-m/2)^2\le m(m+3)\delta+4\sum_{i<j}E_{ij}$. The proof applies Langer's orbifold
       Bogomolov–Miyaoka–Yau inequality to the ruled surface of the column directions, blown up along the contacts,
       with weight $2/3$ on the sections. Two referees recomputed the Chern numbers directly on about 70 surfaces.
     * *Theorem (function-field uniform bound).* For $k=\mathbb C$, $\mathbb R$ or $\mathbb Q(i)$ and the model of
       item 10 (distinct $z_j\in k[t]$ of degree $D$, $z_j\bar z_j=N$, $\deg(z_i-z_j)\le D/2+c$): $M(D,0)\le24$ for
       every $D$ ($\le23$ for odd $D$); $M(D,c)\le47$ for $D>68c$ and $\le24$ for $D>1400c/3$; hence
       $\sup_DM(D,c)\le\max(47,2^{2c+2}(\lfloor68c\rfloor+2))$. With the cube family, $8\le\sup_DM(D,0)\le24$. The
       colleague did not state this corollary; the coordinator and three referees derived it, and the judge
       rechecked the dictionary (stars $P_i=\gcd(z_i,U)$, $E_{ij}=0$, $\delta=\lfloor D/2+c\rfloor-D/2$). The sharp
       statement (R$_4$), $\limsup M=8$, stays **open**. Without the norm condition the maximum at $\delta=E=0$ lies
       in $[12,24]$ (a 12-mark example on an elliptic curve; 24 is the method's limit). Characteristic 0 is
       essential: the statement is false over $\mathbb F_p(t)$, $p\equiv3\ (4)$, with $A^2+B^2$ (Part VI below) and
       over $\mathbb F_2(t)$ with an anisotropic norm (Part V below), where
       (IV.2) fails for a separable configuration on a rational surface over $\overline{\mathbb F}_2$.
     * *Characteristic 2* (**proved**). Over $\mathbb F_q(t)$, $q=2^d$, trace hyperplanes give $m=q-1$ columns at
       weight $m\log(m+1)$ with every divisibility and constant cofactors; (IV.2) fails for $q\ge32$. The
       configuration is a genuine split-form (hyperbola) cluster; no anisotropic norm fits it. Its incidence table
       cannot be realised over $\mathbb Z$ for $d\ge15$ (rigidity). (Part V below removes both limitations: fixed
       field $\mathbb F_2$, anisotropic norm.)
     * *Filter F3* (sharpened to F3+ in Part V below, which removes the second alternative for uniform
       statements). An argument that bounds integer-matrix models without the exact norm must use characteristic 0
       or collisions in residue fields smaller than $m$; every inequality whose characteristic-$p$ analogue the
       trace surface satisfies is excluded (Hodge index, Riemann–Roch, geometry of numbers, Bogomolov on rational
       surfaces; the last via secondary sources).
     * *For the integers: nothing.* One-parameter polynomial templates with bounded $C$ carry at most 24 points;
       the unit (golden/Pell) templates that hold the records sit at $\delta=L/8$, where (IV.2) is vacuous.
       Isolated parents are untouched. The arithmetic analogue is an arithmetic orbifold BMY inequality on
       $\mathbb P^1_{\mathbb Z}$ with $m$ sections, of at least Szpiro strength; no proved substitute was found (any
       constant $\kappa$ in place of 3 would give $m\le2\kappa(\kappa+1)$). Whether integer-matrix models without
       the norm condition exist at $W=O(m\log m)$ is **open** in both directions.
   * ***Part IV, Question B thread: the norm defect*** (`round31/README_IVL.md`). Write $\delta_i=|V_i|^2-n_i$ for an
     integral five-mark model. With each cut marking its two endpoints ($n_i\sim N^{2/5}$) there are models with
     all divisibilities, index one, cofactors $(\log N)^{O(1)}$, $|\delta_i|=N^{1/10+o(1)}$ (every defect exponent
     $>1/4$) and width $N^{-3/10+o(1)}$, avoiding any prescribed proper closed set of ratios. With all
     $n_i\sim\sqrt N$ one pencil gives $o(n_i^{3/5})$ with fixed cofactors. Below $n^{1/4-30\eta}$ (resp.
     $n^{2/5-24\eta}$) at width $N^{-3/10+\eta}$, a model is exact or has $\delta_i=\lambda\ell(V_i)^2$; the latter is
     a genuine cluster for an order $\mathbb Z[\sqrt k]$, $|k|\le N^{3/10}$ (analyst's normal form), neither
     constructed nor excluded. The first family consists of exact clusters for definite forms of discriminant
     $\asymp N^{3/10}$. So a proof of Question B must use local isotropy or the form $X^2+Y^2$ itself (or bounded
     discriminant); archimedean norm accuracy $N^{1/10+\varepsilon}$ and exact values of some form do not suffice.
     Nothing here bears on defect zero. The coordinator's quadratic-form lattice suggestion was wrong (one marked
     line per prime; index $N$).
   * ***Part V: the uniform statement is false over the fixed field $\mathbb F_2$*** (`round31/README_V.md`; five
     referee blocks; the coordinator and the judge reproduced the construction with their own code).
     * *Theorem* (**proved**). Let $\mathbb F_4=\mathbb F_2(\theta)$, $\theta^2+\theta+1=0$, with conjugation
       $\theta\mapsto\theta+1$, so the norm form of $\mathbb F_4(t)/\mathbb F_2(t)$ is $A^2+AB+B^2$ (anisotropic,
       also at infinity; it is not $X^2+Y^2$, which is degenerate in characteristic 2, and the extension is
       unramified everywhere). Let $M_2(D,c)$ be the largest number of distinct $z_j\in\mathbb F_4[t]$ of degree $D$
       with one norm $N\in\mathbb F_2[t]$ and $\deg(z_i-z_j)\le D/2+c$. For $r\ge2$, $m=2^r$, $d=r+1$,
       $V=\mathbb F_2[t]_{<d}$, every integer $L\ge m(r-1)$ and $Y=t^L+\theta$, the $m$ points
       $z_H=\prod_{v\in H\setminus0}(Y+v)$, $H$ running over the hyperplanes of $V$ not containing 1, have the
       common norm $N=\prod_{v\in V\setminus\{0,1\}}(Y+v)$, degree $D=(m-1)L$ and
       $\max\deg(z_H-z_{H'})=D/2+m(r-1)/2-L/2$. Hence $M_2(D,0)\ge2^r$ for every $D\in(2^r-1)\mathbb Z$ with
       $D\ge2^r(2^r-1)(r-1)$, and $\sup_DM_2(D,0)=\infty$. The proof is the identity
       $z_H-z_{H'}=P_J(1)P_J(Y)/Y$ for the subspace polynomial $P_J$ of $J=H\cap H'$; conjugation acts on $Y$ as
       translation by 1. No gcd cost is discarded (the common gcd has degree $\le30$ for $m\le64$ and dividing by it
       tightens the cluster). $N$ is not squarefree here; squarefree parents carrying $m=2^r$ points exist (Part VI.5 below,
       **proved**). The mechanism (subspace
       polynomials agreeing in their top coefficients) is Ben-Sasson–Kopparty–Radhakrishnan (IEEE Trans. Inf.
       Theory 2010); the equal-norm specialisation is the colleague's. The worst-case fixed-field short-arc problem
       was not found in the literature (search snippets only).
     * *Rate.* Along $D_m=m(m-1)(\log_2m-1)$ the cluster has $m\sim\sqrt{2D/\log_2D}$ points; the parent degree is
       $W=2D=\Theta(m^2\log m)$. The record's $2/3$ argument transfers to $\mathbb F_2[t]$ with this norm
       ($W\ge3m\log_2m-O(m)$; referee, not redone by the judge), so the construction sits a factor about $(2/3)m$
       above the growth scale and says nothing about growth. Inside the construction's class (any additive $V\ni1$,
       any $T$ in $Y=T+\theta$, hyperplane sides, a common divisor, unit twists) the parent degree is
       $\Omega(s^2\log s)$ for $s$ points (analyst's Theorem E, **proved**, which supersedes the colleague's V.2;
       the loss sits at small inert primes). This is not a bound for $M_2(D,0)$. Exact values: $M_2(D,0)=4$ for
       $11\le D\le17$; $M_3(D,0)\le4$ for $D\le11$ (form $A^2+B^2$ over $\mathbb F_9(t)/\mathbb F_3(t)$); nevertheless
       $\sup_DM_3(D,0)=\infty$ (Part VI; first 5-point instance found at $D=568{,}386$).
     * *Consequences inside the model* (**proved** unless labelled). (IV.2) fails with $\delta=E=0$ (left side
       62016 at $m=32$). The $\mathbb F_2(t)$ analogues of $A_5(\delta)$ for every $\delta<1/30$, of $(L_5^\circ)$ and
       of $V_E(5)$ are false: the exponent $15/4$ is attained Zariski-densely. The analogue of $(L_8)$ fails
       asymptotically (**plausible, not a theorem**: the star dictionary was not re-derived in characteristic 2). The
       point is separable (all ratios have nonzero derivative), so "Frobenius" is not the explanation, and neither
       is isotriviality (the characteristic-0 model is equally isotrivial and bounded).
     * *Filter F3+.* Call an argument $\mathbb F_2$-valid if every step holds verbatim for $\mathbb F_2[t]$ with this
       norm, all $N$, degree for logarithm, and the inert place at infinity for the archimedean place. Unique
       factorisation, exact norms, residue counts and pigeonhole at small primes, anisotropy, the cut/chord/Plücker
       identities, the chain lemma, Liouville with toric forced divisibility, Schwartz–Zippel and intersection
       theory short of BMY are all $\mathbb F_2$-valid. **No $\mathbb F_2$-valid argument proves the uniform bound,
       nor any bound $M\le f(W)$ with $f(W)=o(\sqrt{W/\log W})$.** In particular F3's alternative "collisions in
       small residue fields" does not suffice for a uniform statement, and no characteristic-free proof of
       $A_5(\delta<1/30)$ or $(L_5^\circ)$ exists. Not affected: growth-rate arguments (fixed-field growth scale is
       **open**); number-field Vojta and the implications $A_5\Rightarrow$ uniform, $(L_8)\Rightarrow$ uniform.
       (As first recorded, F3+ also exempted arguments needing odd characteristic, a reflection-type conjugation
       or squarefree $N$. Part VI removes these: see F3++ below.)
     * *Not evidence against the conjecture over $\mathbb Z$.* The mechanism needs many pairwise cofactor
       differences that are constant, on sets of even size. In characteristic $\ne2$ with a reflection-type
       conjugation a constant difference forces odd size and three points cannot be pairwise of this type
       (analyst's Lemma P, **proved**; it applies to $\mathbb F_3(t)$ and to $\mathbb Z[i]$). Lemma P explains
       only why the characteristic-2 mechanism does not port; it is not a barrier, because over $\mathbb F_p(t)$,
       $p\equiv3\ (4)$, a different mechanism gives unbounded clusters (Part VI below). Do not write "the
       conjecture fails in positive characteristic" ($p\equiv1\ (4)$ and ramified quadratic models are untested).
   * ***Part V, Question B thread: local isotropy*** (`round31/README_VL.md`; three referee blocks). Call a block
     $J_e$ isotropic when $J_e\mid|V_i|^2$ at its marked columns (this propagates along the marked line). Count
     isotropic blocks, not incidences. (i) Without a size bound, full isotropy with fixed cofactors and index one
     allows arbitrarily small width off any closed set ($\Gamma(N)$ shear); the genuine height decouples from $N$.
     (ii) Full isotropy with two-sided width $N^{-3/10\pm o(1)}$ forces subpower multipliers and is the genuine
     problem on the ten-pair-cut profile at the Cilleruelo–Córdoba endpoint $\alpha=10/3$; excluding it off a closed
     set follows from $A_5(\delta)$ for any $\delta<1/10$; by Part VI below it is in fact equivalent to $A_5(\delta)$
     for some $\delta<1/10$ (logically weaker than Question B; no separation known). (iii) In the
     endpoint orientation with $|V|^2/n\to1$ (defect $n^{7/8}$), four isotropic blocks forming a star occur
     densely, and five occur on one curve; six to nine are **open**. (iv) At defect $\le P^{2-\varepsilon}$
     ($P=N^{1/10}$) a nonstar triple of isotropic blocks forces exact norms, and at $\le P^{3-\varepsilon}$ so do
     $\ge8$ blocks. (v) With $k$ isotropic blocks the model yields genuine tuples of exponent $\le(20-k)/3$; a
     dense nine-block family would refute $A_5$. All of this characterises genuine tuples and gives no inequality
     for them. *Closed as routes to Question B:* local isotropy without archimedean size; full isotropy with
     subpower multipliers as a reformulation (it is the Gaussian factorisation); and, on the judge's
     recommendation, the integer-matrix-model programme as a whole (four passes located where models stop being
     fake, and each stopping point is "the tuple is genuine").
   * ***Part VI: unbounded over $\mathbb F_3(t)$ with $A^2+B^2$; squarefree over $\mathbb F_2(t)$***
     (`round31/README_VI.md`; six referee blocks; no referee found an error in a claim of the colleague; the judge
     reproduced the four-mark ternary instance, the squarefree instance at $m=4$ and two referee computations with
     own code).
     * *Theorem (VI.4;* **proved***).* Let $p\equiv3\ (4)$, $\mathbb F_{p^2}=\mathbb F_p(i)$, conjugation
       $i\mapsto-i$ on coefficients: the norm form of the constant extension $\mathbb F_{p^2}(t)/\mathbb F_p(t)$ is
       $A^2+B^2$, anisotropic, with the infinite place inert. Let $M_p(D,c)$ be the largest number of distinct
       $z_j\in\mathbb F_{p^2}[t]$ of degree $D$ with one norm $N\in\mathbb F_p[t]$ (any multiplicities) and
       $\deg(z_i-z_j)\le D/2+c$. For every $m\ge2$ there is $D_1(m)=\exp(O(m\log m))$ with $M_p(D,0)\ge m$ for
       every $D\ge D_1(m)$. Hence $\sup_DM_p(D,0)=\infty$ and
       $\liminf_DM_p(D,0)/(\log D/\log\log D)\ge\frac12\log_2p$. In the model's units ($W=\deg N$ for $\log R$)
       this is $m\asymp\log W/\log\log W$, far below the growth scale.
     * *Construction.* $q=p^r>2^{m-1}$ with $r$ odd, $k=\lceil\log_p(2m+1)\rceil$,
       $V=\{i\sum_{b<k}x_bt^b:x_b\in\mathbb F_q\}$, $n=q^k$, $h=n/q$; $m$ pairwise nonproportional
       $\mathbb F_p$-linear forms $\lambda_j$, signs $\varepsilon_j$, $\chi$ the quadratic character of
       $\mathbb F_q$; $F_j(X)=\prod_{v\in V}(X-v)^{1+\chi(\varepsilon_j\lambda_j(v))}$. The coefficient Frobenius
       fixes $F_j$, so $F_j\in\mathbb F_{p^2}[t,X]$. Conjugation is the $q$-power map on $\mathbb F_{p^2}$, sends
       $v\mapsto-v$, and $\chi(-1)=-1$, so $F_j\overline F_j=P_V(X)^2$. With $A=P_{\ker\lambda_j}$
       ($\mathbb F_q$-linearised), $F_j=A(A^{(q-1)/2}-c^{(q-1)/2})^2=P_V+2c^{q-1}A-2c^{(q-1)/2}A^{(q+1)/2}$, so
       pairwise differences have $X$-degree $\le(n+h)/2$. Signs can be chosen so that the common factor has degree
       $g\ge2^{1-m}n>h$; divide by it and put $X=t^L$, $L=2(k-1)n+1$. No gcd cost is discarded (the full gcd is
       larger and dividing only tightens). First instances for $p=3$: 4 points at $D=6981$ (reproduced by two
       referees and the judge), 5 at $D=568{,}386$, 6 at $D=3{,}582{,}084$ (via a full-gcd formula; certified
       without it at $6{,}880{,}250$).
     * *Caveats, mandatory when citing.* (1) $N$ is far from squarefree; most primes are squared, and for
       $k\ge3$ multiplicities at $\mathbb F_{q^2}$-rational places are about $q^{k-2}$, unbounded in $m$.
       Squarefree or bounded-multiplicity $N$ in odd characteristic is **open**, and this mechanism cannot reach
       it. (2) Only the constant field is fixed: the bulk of $\deg N$ sits on primes whose residue fields contain
       $\mathbb F_{q^2}$ and grow with $m$. (3) Only $p\equiv3\ (4)$ and the constant extension: for
       $p\equiv1\ (4)$ the mechanism fails formally (every half-system $S\subset\mathbb F_q^\times$ has
       $\delta_S\ge(q-1)/4$; referee, **proved**), and ramified quadratic models are untested. (4) The constant
       $\log_2p$ in place of $\frac12\log_2p$ is plausible only.
     * *Why it does not port* (analyst; **proved** at the level of mechanism). The one step without a
       characteristic-0 counterpart is a single parent with many same-degree right composition factors
       ($P_V=A_j^q-c_j^{q-1}A_j$ for every kernel). If $\operatorname{char}k\nmid\deg f$ and $f(A)=f'(A')$ with
       $\deg A=\deg A'$, then $A'=\alpha A+\beta$ (Ritt/Engstrom-type; the $d$-th root of $P$ as a Laurent series
       exists when the characteristic does not divide $d$). Descent, common division, specialisation and character
       exponents all have integer analogues and do not help. The integer counterpart class (polynomial templates
       in $\mathbb Z[i][u]$) is bounded by 24 (Part IV); computed integer analogues have width exponents
       $0.94$–$1.00$ against the $1/2$ needed. This is an obstruction to porting the mechanism, not a theorem about
       integer clusters.
     * *Theorem (VI.5;* **proved***).* For $d\ge3$, $m=2^{d-1}$ and every integer $L\ge4m(d-2)$ there is a
       **squarefree** $N\in\mathbb F_2[t]$ of degree $2(m-1)L$ carrying $m$ distinct points of $\mathbb F_4[t]$ at
       $c=0$; least supplied parent degree $8m(m-1)(\log_2m-1)$. ($V=\{a+(f^2+f)^2\}$, a space of squares with
       trace zero at even-degree primes, and $T$ chosen by CRT with a trace condition; threshold exact for this
       $V$; verified to $m=64$.) So squarefreeness alone does not rescue the uniform statement in characteristic 2.
     * *Parity-code condition (VI.1;* **proved***, binary model only).* On every genuine parent over
       $\mathbb F_2[t]$, every nonzero word of the anchored binary span of the split-prime parity rows has weighted
       degree $>W/4$ (logarithmic derivative; $d/dt$ kills squares), so the rank is $<0.457n$. This excludes from
       the binary model the record's squarefree skew (Paley) double core, all fully flipped Hadamard cores and
       their residue designs. It does not exclude low-rank (Walsh) profiles, the shape of every genuine binary
       cluster; two Walsh cores on $x,x^3$ still allow $W=O(m\log m)$. **No integer analogue** (barrier 6).
       Do not write "the fakes exist because $\mathbb Z$ has no derivation".
     * *Class bounds over $\mathbb F_2[t]$* (**proved**). Complete Hadamard prime profiles force $m$ a power of 2
       and $D\ge m(m-1)(\log_2(m-1)-3)/6$ (so $\Theta(m^2\log m)$ with VI.5); full hyperplane designs with any
       translation $B$ cost $\Omega(s^2)$ (asymptotic, no threshold). Neither bounds $M_2(D,0)$: with $k\le m/2$
       columns the method gives only $k\log k$. Nonconstant $B$ beats Part V by a constant factor (8 points at
       point degree 68 against 104; 16 at 472; 32 at 2607). **The growth scale over $\mathbb F_2[t]$ is open**:
       $3M\log_2M-O(M)\le W(M)\le2M^2\log_2M$.
     * *Filter F3++* (replaces F3+; **extended to F3+++ in Part VII below**, which removes bounded multiplicity
       and finite ramification from the "still available" list). Call an argument function-field valid if every step holds verbatim, with
       degree for logarithm and the inert degree-one place at infinity for the archimedean place, in either
       setting: (i) $\mathbb F_2[t]$ with the norm of $\mathbb F_4(t)/\mathbb F_2(t)$, for all $N$ or for
       squarefree $N$ only; (ii) $\mathbb F_p[t]$, $p\equiv3\ (4)$, with $A^2+B^2$, for all $N$. **No
       function-field-valid argument proves the uniform bound.** One valid in (i) cannot prove $M\le f(W)$ with
       $f=o(\sqrt{W/\log W})$, even for squarefree $N$; one valid in (ii) cannot prove $M\le f(W)$ with
       $f=o(\log W/\log\log W)$ (it could still prove $M\ll W^\varepsilon$). Admissible, hence insufficient:
       unique factorisation, exact norms, anisotropy and inertness at infinity, a reflection-type conjugation and
       a fixed odd characteristic (ii), squarefreeness (i), residue counts and pigeonhole at small primes,
       cut/chord/Plücker identities, the chain lemma, Liouville with toric forced divisibility, Schwartz–Zippel,
       Lemma P. Still available to a proof over $\mathbb Z$ (none known to suffice): characteristic 0 used
       essentially (Newton linearisation, BMY, uniqueness of composition factors); archimedean order and
       curvature; Baker-type input; squarefree or bounded-multiplicity $N$ in odd characteristic (**open**);
       ramification ($\mathbb Q(i)/\mathbb Q$ ramifies at 2, constant extensions nowhere); absence of a
       coefficient Frobenius. Growth-rate arguments are not filtered at the integer growth scale: the record's
       $W\ge3M\log M-O(M)$ transfers to both settings and no function-field family reaches it.
     * *Not evidence against the conjecture over $\mathbb Z$* (judge's ruling). Both mechanisms rest on additive
       polynomials; the characteristic-0 function-field model shares the reflection, anisotropy, exact norms and
       isotriviality with $\mathbb Z[i]$ and is bounded; the failing step is identified and impossible in
       characteristic 0; and the ternary clusters live at $m\asymp\log W/\log\log W$ on square-full parents with
       growing residue fields. Do not write "false in positive characteristic", "false for sums of two squares
       over function fields" (ramified models with inert infinity are unbounded at cap two, Part VII; the model
       $A^2+tB^2$ with ramified infinity is computed at $m=3$ only), or "the growth rate over $\mathbb F_q$ is
       $\log W/\log\log W$" (only a lower bound).
   * ***Part VI, Question B thread: the ten-pair slice*** (`round31/README_VIL.md`; six referee blocks; nothing
     false). Recorded in §3: $A_4(\delta)$ is false for every $\delta<1/6$; the slice $(\ast)$ is equivalent to
     $A_5(\delta)$ for some $\delta<1/10$ and is still undecided; the fixed-field unit-sum theorem (essential
     exponent $\ge4$ on a thin class, ineffective, false for varying orders); the cube-class cover has $\kappa=0$;
     the moving affine pencil is rigid. *Closed:* the $(\ast)$ sub-thread as a route to Question B. *Open, narrow:*
     resonant templates in one fixed quadratic order (numerical test first).
   * ***Part VII: multiplicity two and ramified models*** (`round31/README_VII.md`; six referee blocks; nothing
     false; judge's score on the goal 1.5/10).
     * *Theorem VII.1* (**proved**; the separation step rechecked by the judge with negative controls). Let
       $p\equiv3\ (4)$ and $\mathcal O=\mathbb F_p[t,\eta]$ with $\eta^2=-1$, or with $\eta^2=\Delta(t)$,
       $\deg\Delta$ even and positive with nonsquare leading coefficient (infinity inert; finite ramification;
       conjugation fixes the constants). Critical clusters on parents with **every prime exponent at most two** are
       unbounded: $\liminf M^{(\le2)}(D)\,\log\log D/\log D\ge\frac12\log_2p$, with an $m$-point construction for
       every $D\ge C_{p,\Delta}2^mn^2$. *Device:* labels that are $p$-th powers,
       $V=\{\eta^pC^p\sum x_jt^{pj}\}$, so $dv=0$; $C$ divisible by $\Delta\prod_{j<k}(t^{q^j}-t)$ kills every
       collision; $T=1+C^pt^s$ is $1$ on the collision set; hence $R=\prod_v(T-v)$ is squarefree and the norm
       divides $R^2$. Exact width threshold (referee): $L_{\min}=\max(b+1,\ p\deg C+2,\ \lceil b(n-h)/(g-h)\rceil)$.
       Instances for $p=3$: 4 marks at $D=578{,}396$ (constant model), fully expanded by two referees; ramified
       $\Delta=2t^2+1$ at $D=619{,}801$.
     * *Caveats.* The cap is reached in a thin sense: squared primes carry $0.98$ of $\deg N$ at $m=4$ and a
       fraction tending to 1. $N$ is coprime to the ramified primes, so ramification plays no role in the
       mechanism. With inert infinity the order is never a PID (class number even), unlike $\mathbb Z[i]$.
       Ramified infinity ($A^2+tB^2$) is computed at $m=3$ only. $p\equiv1\ (4)$ is open.
     * *Exactly squarefree $N$ in odd characteristic is open.* The Part VI pattern cannot be specialised to a
       squarefree norm (**proved** for $m\ge16$, $q\ge27$; one family only). Exhaustive tables over
       $\mathbb F_9(t)/\mathbb F_3(t)$: $M_3^{\rm sf}(D,0)=2,4,3,4,3,4,4,4,4,4,4,5,4,5,4,5$ for $D=1,\dots,16$ and
       $4,5,4$ for $D=17,18,19$ (referee); the first five-point cluster for any parent is at $D=12$ and is
       squarefree. The squarefree table equals the all-parent table for $D\le13$, and all-parent clusters are
       unbounded only beyond $D\sim10^5$, so **the tables are no evidence of boundedness**.
     * *Binary growth scale* (**proved**, block form, sharp at $m=8$, false at $m=4$): a single Hadamard profile
       of order $m\ge8$ with all rows and a proper subset of the nonconstant columns never gives a critical
       cluster. The scale $3M\log_2M-O(M)\le W(M)\le2M^2\log_2M$ is unchanged. The colleague's remark that the
       Gram step is available over $\mathbb Z$ is withdrawn (off-diagonal entries are $\le4\log C$, not $\le0$).
     * *Filter F3+++* (replaces F3++). Add to the settings of F3++: (iii) $\mathbb F_p[t]$, $p\equiv3\ (4)$, with
       $A^2+B^2$ or with the order $\mathbb F_p[t,\eta]$, $\eta^2=\Delta$ as above, for $N$ with every prime
       exponent at most 2. An argument valid in (ii) or (iii) cannot prove $M\le a\log W/\log\log W$ for any
       $a<\frac12\log_2p$. Newly admissible, hence insufficient: bounded prime multiplicity (cap two) in odd
       characteristic; finite ramification; a conjugation that is not a coefficient Frobenius. Still available
       to a proof over $\mathbb Z$ (none known to suffice): characteristic 0 used essentially (Newton
       linearisation, BMY, uniqueness of composition factors); archimedean order and curvature; Baker-type and
       transcendence input; exactly squarefree $N$ in odd characteristic (**open**); $p\equiv1\ (4)$ (**open**);
       ramification at infinity (not filtered).
     * *Meaning for $\mathbb Z$* (judge). The closed escape had no client: every squarefree-specific step of the
       record already holds over $\mathbb F_2[t]$ with squarefree $N$, where clusters are unbounded. Not evidence
       against the conjecture, for the same reason as Parts V–VI. Do not write "unbounded for squarefree parents
       in odd characteristic", "false for $A^2+tB^2$", or "the census suggests boundedness". **Judge's
       recommendation: close the function-field programme** and return integer effort to characteristic-0 or
       archimedean input.
   * ***Part VII, Question B thread: five marks over $\mathbb C(t)$*** (`round31/README_VIIL.md`; seven referee
     blocks; nothing false; judge's score on the goal 1.5/10). Recorded in §3: the characteristic-0 function-field
     analogue of $A_5(\delta)$ is **open**; it is a theorem on $b$-sparse presentations (the colleague's Theorem
     VII.2: a moving-coefficient gcd estimate, Dirichlet plus a weight-homogeneous resultant, a special case of
     Guo–Sun–Wang, then one exponent vector and the frozen-support family, which must have **exact** supports);
     the circle bound $20h\le9D+18r$; dense polynomial Pell families at $W=4h+2$ for five marks (reproduced by the
     judge with exact arithmetic; $h=4n$ for **every** choice of distinct seeds, a correction); failure in
     characteristic 2 with a uniform term budget (rank corrected from "$\le1$" to 2, still $<4$). Referee
     additions: $H_5(7)=3$, hence $M(7,0)=4$ **exactly** (Gröbner bases over $\mathbb Q$); rigid five-mark
     configurations with $2h>D$ at $D=11,13$ (numerical); the remark that "arbitrary complex coefficients" are
     height-0 constants, so VII.2 is not Diophantinely stronger than Part VI.1. The referees and the judge also
     draw an integer consequence from the Pell families, $\alpha^{\rm ess}_e(-K)\le4$ ($A_5(\delta)$ false for every
     $\delta<0$), recorded in §3 with its verification by three further checkers. *Closed* (judge's
     recommendation): this thread as internal work; the unrestricted question is one for function-field Vojta
     specialists.
   * ***Part VIII: low-height primes on integer circles*** (`round31/README_VIII.md`; seven referee blocks; no
     error found; judge's score on the goal 3/10, quality as mathematics about 8/10). After the function-field
     programme was closed the colleague was asked to return to the integers and chose the low-height window.
     * *Theorem VIII.1* is Theorem 1.2 of §1, with its consequences $P_+(N')\ge(2-o(1))M\log M$ and the
       coefficient $2/\kappa$. Referee sharpenings: the term $\frac12\log2$ (since $1+i$ divides every cofactor),
       hence $M\le1+\omega(N')$ for $C'<\sqrt2$ where the colleague stated $C<1$; a block mod $2$ gives
       $M\le\omega(N')+2$ for $C'<2^{3/4}$; $K_C=\Theta(C^4/\log C)$. Data (rank referee and judge, exact):
       rank deficiencies $1,\dots,5$ first occur at arc constant $2.83,4.00,5.18,6.64,8.75$ on circles
       $N\ge10^4$; none below chord constant $2$. The 10-point cluster $N=1176852625$ has affine rank $5$, raised
       to $6$ by the mod-3 block and to $9$ by the mod-7 block.
     * *No growth from it for general $R$* (analyst referee, `round31/ref8_extend.md`). Inserting Theorem 1.2 into
       Theorem 1.1 leaves the constant $3$ (Prop 3.2); every inequality derived from pair instances of axiom (R) at
       moduli $\le Q_0$ holds in the point-residue fakes $A_1(Q_0)$ (Prop 5.1), so a separating input must use
       genuine prime arguments or prime-level residues, the archimedean order, or primes of $N'$ above $2M$. The
       only obstruction between $3$ and $4$ in Theorem 1.1 is small split primes dividing $N'$ with balanced
       levels.
     * *Height-order inequality (VIII.6)* (**proved**; low novelty): §5 item 4. It bites only on coordinates of
       finite order modulo the row lattice; in a pinned cluster every coordinate at a prime $\le M^A$ has order
       $\ge\exp(cM)$. Corrections: no window hypothesis and no $\delta$-dependence are needed; the 10-point
       illustration lies outside the lemma's range. The Baker–Wüstholz input was checked only through
       `round31/ref_baker.md`.
     * *Window fakes (VIII.3, (VIII.10))* (**proved**): the bullet in §1. Devices: random assignment of primes in a
       wide interval with martingale control of the Paley margins (bypassing judge J1 of round 21); a square-root
       support bound for $Kv$ mod $4$ and mod odd $\ell$ (the minimum-distance bound for extended
       quadratic-residue codes; corrected to $\sqrt q+1-h$; mod $2$ alone fails), used as an erasure estimate;
       integral moment lifts with base-$B$ digits and a Newton-identity count of words with vanishing initial
       moments; residue witnesses ordered by auxiliary Gaussian prime directions. Wording corrections: at fixed
       $\rho$ the coefficients are $7.96/\rho$ and $23/(3\rho)$, and $23/3$ holds along a diagonal sequence (one
       referee's "false as stated" was overruled by the judge); inside the window only for $\delta<1/12$. It
       answers round 21's question (Q) in that range and contradicts no round-21 theorem; what falls is the
       round-21 heuristic that window fakes admit no positive-probability recipe. Brackets: $[3,23/3]$ for inputs
       up to point residues at $2M$.
     * *Area test and spread fakes (VIII.4)* (**proved**; the test has minimal novelty):
       $M\le1+K_C+O((1+\sqrt P\sigma)(1+\sqrt P\eta))$ when the primes lie in $[P,Pe^\sigma]$ and the chosen
       Gaussian arguments in an interval of length $\eta$. It detects the first window fake; the spread version of
       the fake (prime angles separated by $\ge1/(6q)$ up to conjugation) evades it. Referee addition: the spread
       fakes are $L$-spread in the sense of `round21/prove.md` Def 5.1, which discharges the fake side of round 21
       Cor 5.4; the genuine side (Caveat 3) is **open**.
     * *Judge's ranking of what could still move the integer rate*: whether the per-pair point-residue
       inequalities force $W\ge(3+\eta)M\log M$ (estimated 4%; the only item whose success would improve $\frac23$
       for every $R$), then interpolation between Theorems 1.1 and 1.2 for $2<\kappa<3$ (2%). Do not write "the
       low-height window is uncapped", "Theorem 1.2 improves $\frac23$", "fakes violate the rank bound" or
       "coefficient $\frac12$ for integer radii" (only squarefull radii).
   * ***Part IX: projections, and residue witnesses at all heights*** (`round31/README_IX.md`; eight referee
     blocks; no error found; judge's score on the goal 2.5/10, quality as mathematics about 8/10). The colleague
     was asked to beat the rank bound of Theorem 1.2 by an unbounded factor ($M=o(\omega(N'))$ would improve the
     rate), or to prove the strongest obstruction.
     * *Theorem IX.2* is Theorem 1.3 of §1, with its row certificate. Restricted and position-dependent; see the
       scope paragraph there. Function-field check: proved for constant quadratic extensions only.
     * *All-height factored fakes (IX.3)* (**proved**; the judge rechecked the same-level deletion at parent primes,
       0 failures in 22,992 identities). The anchored exponent matrix $D$ of the Legendre profile has an integral
       right inverse $U$; with $V=U(H-D)/2$ and $F=I+2V$ one has $DF=H$, and the witnesses
       $r_j=i^{u_j}\pi_j\prod_l(\pi_l/\bar\pi_l)^{V_{jl}}$ have norm $p_j$ and make every mark ratio equal to
       $\gamma_i/\gamma_{i'}$ at every odd prime power; at a parent prime the singular factor is deleted. The global
       collision cost of a character is bounded once for all moduli because the colliding powers divide one nonzero
       Gaussian integer. Corrections: $m/\omega(N')\to\frac12$ only along a diagonal; the witness heights depend on
       $U$ and are enormous; "these witnesses fail the all-monomial test" is proved for these and for $S$-unit
       witnesses only, at an astronomical modulus (computed for one $V$: failure already at modulus $3m$ for most
       labels at $q=11,19$).
     * *Window to $\delta<1/6$ (IX.4)* (**proved**): a Halász-type atom bound
       $\sup_b\Pr(\sum w_iX_i=b)\le q^{o(1)}(A\sqrt{qh})^{-1}$ removes the condition that closed the parameter system
       at $1/12$; it now closes at $1/6$, where the square-root support bound and the lift sup norm $q^{3/2}$ bind
       together.
     * *Cross-gap fakes (IX.5)* (**proved** modulo the Bourgain–Dyatlov Fourier-decay theorem, cited, hypotheses
       checked, full text not read): each Legendre column replicated $g$ times, angles drawn from a measure on a
       continued-fraction limit set conditioned on prescribed sums. Outside the window (prime height $m^d$, $d$
       large and unspecified). Referee addition: these fakes satisfy round 21's union large sieve at the genuine
       constant. As written the test includes the owner's own direction, where genuine angles fail it; the record
       forms exclude it.
     * *A fake in a rational-direction window (IX.6)* (**proved**): the integrality of projections is not implied by
       the point-residue axioms. The remark that the Gram identity with real cofactor magnitudes permits the fakes
       is correct and near-tautological; integrality of the cofactors and their cocycle relation stay outside.
     * *Rigidity (IX.7)* (**proved**; judge rechecked line by line): §1. Referee addition: at every finite modulus
       every factored datum is the reduction of some $r_j=\pi_j\gamma/\bar\gamma$ (Hilbert 90 and the Chinese
       remainder theorem), so the rational-lift hypothesis has content only in the limit over all moduli.
     * *Judge's reading*: the thread has capped every local and character-tested input at every height; the residue
       route to the rate reduces to P3-all$^*$ (open since round 6) or to the problem itself. Recommended: one more
       pass on P3-all$^*$, and to stop the residue route if it returns another fake. Do not write "P3-all is
       capped", "every residue-type input is capped or is the problem itself", "point residues at moduli $M^A$ are
       capped uniformly in $A$" (zero-sum characters only), "Part IX proves growth in the window", or "the max-scale
       cross-gap is capped" without "single-label, outside the window".
   * ***Part X: finite cutoff and coordinate rays*** (`round31/README_X.md`; six referee blocks; no error found, one
     sentence overstated; judge's score on the goal 2/10). Offered the cocycle on one fixed parent or the residue
     class P3-all$^*$, the colleague took the second. **P3-all$^*$ is not settled in either direction.** Proved:
     effective finite rigidity for bounded-height rational lifts (§1; polynomial modulus on the dyadic parent of all
     primes $\equiv5\ (8)$; in general the test exponent is $2^{T-1}$, $T=\max v_2(\ell-1)$); a partial
     P3-all$^*$ fake (all single-label rays and all short words with nonzero class sum, every power, unit and
     modulus, together with P3-char), on which the least lift heights diverge with the modulus; no short exact
     witnesses on the IX.3 design (a mod-four bound; exact minima 20, 72, 110 at $q=11,19,23$, judge). The
     conditional weight bound X.4 (if the residue tuples of the labels at inert primes up to $K\log m$ generate a
     group of bounded rank, then $W/(m\log m)\to\infty$) is a reformulation: its hypothesis fails on every genuine
     parent with polynomially bounded primes. Corrections: the sentence that the fake covers "a fraction tending to
     one of the tests" concerns data-independent, class-separable word sets only, and the uncovered words are
     exactly the ones a pigeonhole argument would pick; the table of least cutoffs (439, 3881, 8681, 5783) is for
     lifts with parent-prime denominators; the coordinator's phrase "the twist must be $\pm1$" was wrong (units and
     conjugation are allowed on genuine circles). **Computed (referee, spot-checked by the judge):** on
     $q_P=\prod_{p\equiv1\,(4),\,P\le p\le2P}p$ no two lattice points lie within $8q_P^{1/4}$ of each other for
     $200\le P\le450$ (the closest pair is between $468$ and $8\cdot10^7$ times $q_P^{1/4}$ away), against a proved
     bound of order $s/3$ points. *Judge's reading*: the residue route has been mapped to its end and carries no
     rate; the one live item left on this thread is the cocycle on a fixed parent. Do not write "P3-all$^*$ is
     capped", "Part X gives a finite-cutoff lower bound", "X.4 proves a growth rate", "the omitted difference must
     fail", or "$M_{\max}(q_P,C)=1$ is proved".
   * ***Part VIII, Question B thread: hardness of Question A*** (`round31/README_VIIIL.md`; seven referee blocks;
     nothing false, one citation corrected; judge's score on the goal 2/10). The colleague was asked for the weakest
     hypothesis that gives the conjecture, or for a hardness statement, and took the second. Recorded in §2
     ("Strength of $(L_8)$"): $(L_8)$ restricted to pairwise coprime square blocks implies the two-residue
     statement at length $q^{1/2+\epsilon}$, Cilleruelo–Granville's Conjecture 9 for all $a$, and Rudin's
     conjecture; full $(L_8)$ implies $O(1)$ divisors of $n$ in an invertible class modulo $q\ge n^{1/4-\kappa}$.
     The judge rechecked the two steps the colleague asked to have prioritised (the square-block identity and the
     escape from the fixed exceptional set) and that the data are admissible for the $(L_8)$ that Theorem D derives
     from Vojta's conjecture (blocks need not be contact integers: they differ in about 2% of 2.3 million checked
     instances, within the record's bound (3.2)). **No referee refuted or weakened $(L_8)$.** Recorded in §5 item 5:
     round 8's QS at an absolute height threshold is false in every dimension; the conditional growth statements
     from relative-height counts are correct implications from non-standard, untestable hypotheses, so the sentence
     "no intermediate growth bound is known from a weaker standard conjecture" stands. Referee additions: a
     weighted Prouhet–Tarry–Escott dictionary for two-class families, with exponent-2 families of five and six
     points; a characteristic-2 Hadamard family with unboundedly many points at exponent exactly 2; the coincidence
     of the Liouville exponents of the two-residue problem and the arc problem for every $k$; the correction of
     three sentences of the record that had Cilleruelo–Granville deriving Rudin from arc conjectures. The
     colleague could not derive growth from abc or Bombieri–Lang, nor Rudin from the uniform bound. *Judge's
     reading*: Question A is not a realistic standalone target; this is information about difficulty, not
     progress; Question B is untouched. Do not write "the uniform bound implies Rudin", "Question A is equivalent
     to Rudin", "$(L_8^{pc})$ implies Rudin", or "Vojta's conjecture is known to imply Rudin" (say "we found no
     earlier derivation").
   * ***Part IX of the Question B thread: the critical two-class problem*** (`round31/README_IXL.md`; eight referee
     blocks; nothing false; judge's score 2.5/10). The twin of the arc problem (at most two classes modulo each
     $p^f\parallel q$, intervals of length $C\sqrt q$) is equivalent to Cilleruelo–Granville's Conjecture 9 at the
     uniform scale and at the little-oh scales $\log q/\log\log q$ and $\omega(q)$ (**proved**); so its rate
     question is a named open problem: do intervals of length $\sqrt q$ contain $o(\omega(q))$ roots of
     $x^2\equiv a\ (q)$? No transfer of statements between the twin and circles is known in either direction;
     methods transfer both ways. **Proved:** the random-parent rate for the twin and for circles (Part III bullet
     above), including bounded exponents and the smooth scale; the same-norm comparison of fakes and genuine arcs;
     a squarefree relaxation (an edge table with all pair magnitudes, gcds, valuation patterns at every prime,
     integral cofactors and Sidon signatures) that reaches the larger-sieve constant
     $\frac12\log q/\log\log q$ and has no integer realisation up to length $q^{1-O(1/m)}$: a barrier for
     valuation-based arguments, not for arguments using signed residues of single differences; no affine 3-flat of
     twin signatures at length below $q^{2/3}/4$ (sharp; characteristic 0 is used, which is why the
     Sylvester–Hadamard families of characteristic 2 have no integer counterpart); at most $(1+C^2)(1+2C^2)$ pairs
     per signature difference in the twin; $M,B\ll n^{1-2^{1-2d}}$ ($n^{7/8}$ at $d=2$) for profiles whose bits are
     Boolean functions of $d$ affine forms on $n$ labels, in both problems and in angle models (binary rank
     $O((\log n)^d)$, so no recorded fake is in the class). Not new: $M\ll_C\sqrt{|A+A|}$ on circles (the
     distinct-label lemma of Part II) and the one-large-prime-power criterion (round 8). The printed lattice form
     $(L_8^{pc})$ had empty input (§2). **The growth rate for general $R$, the uniform bound and Questions A and B
     are unchanged**; worst-case parents, even squarefree with comparable primes, are open in both problems.
     *Judge's reading*: at the smooth critical height the truth on typical pool parents is $O(\log s)$, very
     probably $O(1)$, so the order-$s$ barrier is a worst-case phenomenon and not a height phenomenon; the input
     that separates genuine from fake, the exact additive identity with integer cofactors, is identified but has
     been used only through the entropy of the parent. Removing the exceptional parents by counting over exponent
     patterns and cofactor tables is closed as a method for $P\le s^2$. Do not write "the uniform bound holds for
     almost all smooth parents", "Part IX gives the rate $O(\log\log R)$" (without "on random parents of a pool"),
     "the exceptional parents are polynomially rare", "Part IX covers the extremal profile of Theorem 1.1", "the
     relaxation isolates the triangle identity" (without "equivalently, the signed pair residues modulo $q$"),
     "fakes fail the same-norm test", "$(L_8^{pc})$ is open" (printed form), "$(L_8^{pc})$ is Rudin-hard" (either
     form), or "the twin rate implies the circle rate" or its converse.
   * ***Part X of the Question B thread: signed residues on one fixed modulus*** (`round31/README_XL.md`; five
     referee blocks; nothing false; judge's score 2/10). In the critical two-class problem ($q$ squarefree, interval
     $C\sqrt q$, classes at each $p\mid q$ differing by $\delta_p$), a *sign-free table* prescribes for each pair of
     marks an integer $d_{ij}$ with $0<|d_{ij}|\le C\sqrt q$, $\gcd(d_{ij},q)=h_{ij}$ and
     $d_{ij}\equiv\pm\delta_p\pmod p$ at every separating $p$, the sign free per pair and per prime. **Proved:** on
     every triple $q^2g^2$ divides the Heron product $(a+b+c)(-a+b+c)(a-b+c)(a+b-c)$ of the three magnitudes ($g$ the
     product of the primes constant on the triple), so a nondegenerate triple has a side $\ge\sqrt{qg}$ (sharp, also
     at $g=1$). With Menger's line-embedding theorem (1928) this gives
     $B_{\rm int}\le B_{\rm free}\le2\max(4,B_{\rm int})$ whenever $q$ has a prime $p>C^2$ (no loss at $C=1$), and on
     parents with all primes in $[P,2P]$ every sign-free table of $B$ rows becomes a genuine configuration on the same
     modulus and interval after deleting $O_C(1+s/B)$ rows. So the pair-level square congruences already force
     genuineness: among the relaxations studied in the record nothing lies strictly between the valuation data of
     Part IX (which reach the larger-sieve constant $\frac12$) and the genuine problem. No sign-free table exists on
     the Paley cut pattern of Part IX for any modulus and prime assignment: on the full core at length $C\sqrt q$ and,
     with equal numbers of primes per column, at $q^{5/8-\epsilon}$; and at critical length on every set of at least
     $(3/5+\epsilon)m$ rows (two primes per column) or $(2/3+\epsilon)m$ rows (one prime). For genuine configurations
     only the row-subset statement is new. **Open:** a worst-case $o(\omega(q))$ bound; a genuine example of order
     $\omega(q)$; the half-row Paley pattern, which is the row class of the circle fakes (excluded computationally at
     $m=19,23$ for large primes). **The mechanism does not transfer to circles**: separated squared chords modulo a
     split prime are not constant (norm 65: chords of squared length 52 and 4 at the prime 5), circle triangles never
     degenerate, and the only Heron identity there is $\mathrm{area}=2g|t_{12}t_{23}t_{13}|$. Linear
     (column-orthogonality) certificates succeed on a cut pattern exactly when its lattice has no nonconstant short
     vector, a problem of the same type as the twin itself. **The growth rate for general $R$, the uniform bound and
     Questions A and B are unchanged**; on the judge's advice the two-class thread is closed. Do not write "Part X
     proves Conjecture 9 for Paley moduli", "the sign-free relaxation is a new barrier", "there is no relaxation
     between valuation data and the genuine problem" (say "none among those studied"), "Paley profiles are
     impossible in the twin" (without "full core, or at least $3/5$ of the rows"), "the circle fakes are impossible
     in the twin", "column orthogonality is closed as a route" or "gives a rate beyond the sieve constant", or "the
     twin result implies that pair-level circle data force genuineness".
   * *Withdrawn or corrected in the round:* in Part X of the Question B thread, the line lemma is Menger's theorem
     (citation added), and a referee's claim that column orthogonality is closed as a route was overturned by the
     judge (it was shown only for patterns carrying $O(C)$ genuine marks); the claim that the lattice form $(L_8^{pc})$ is open (as first printed
     it has empty input; corrected to $(L_8^{pc})'$, §2); the classification of the cocycle of reduced chords as
     capped (narrowed to its local consequences, triangle areas and abc, §5 item 6); a sentence on "different local obstructions" for the round-24 family
     (both forms are excluded at 2 and at 3); the coordinator's "weight $D/2$ is exactly the lattice minimum" (only
     a lower bound) and "every determinant is critical" (see Lemma K1 as now worded); the headline "17 points" of
     the random-parent theorem (superseded by 2 and 1); "asymptotically exact radii" (relative error
     $O(N^{-1/10})$); in Part V, the colleague's lower bound V.2 omitted unit-twisted points (repaired, loss $\le4$)
     and is empty below $|V|=2^{48}$ (superseded by Theorem E), and her construction is not optimal in its own
     family (8 points at degree 104 against 112); in Part VI, nothing of the colleague's is withdrawn, but her
     five-mark ternary certificate ($D=1{,}231{,}594$) is far from minimal, her rate constant is $\frac12\log_2p$
     as proved, her attribution of the non-transfer to Frobenius descent is corrected (descent ports; the
     composition-factor step does not), "strictly weaker than Question B" is replaced by "logically weaker, no
     separation known", and the coordinator's F3+ exemption list is superseded by F3++. Unreproduced: two eliminations for the strict linear-block system, whose decision is still
     open and peripheral either way.

## 6. Pointers

Repository `craigxchen/Jarnik`, branch `claude/upbeat-einstein-p0peju`, directory `docs/claude_review/`. Each
write-up `X.md` has its referee report `referee_X.md` next to it, or, in round 20, a kill report `kill_X.md`. Its
exact checks are in `X_checks/`. The round READMEs give the referee-corrected statements.

| item | file |
|---|---|
| consolidated account (status labels, corrections) | `overview/overview.md` (§§4, 5.4–5.5, 7.1, 7.1A, 7.6, 8.2) |
| Theorem 1.1, the $(2/3+\varepsilon)$ bound | `two_thirds/proof.md` (Thm 7.5, Cor 8.1, §10) |
| Vojta on $\mathrm{Bl}_YX_7$ | `growth/conditional.md` (Thms 5.1–5.2), `growth/referee_conditional.md` |
| $(L_8)$; Vojta on $\overline M_{0,8}$ implies uniform | `growth/ptolemy.md` (§6 Thm D, §7), `growth/referee_ptolemy.md` |
| $(L_5)$, $(L_6)$ false; Liouville barrier on $\overline M_{0,k}$ | `threshold/L8.md` (Thms A, B, Prop C), `threshold/referee_L8.md`, `growth/referee_ptolemy.md` §3.3 |
| round 18: census, sandwich, lattice form, lifts, Prop 2.6 | `round18/README.md`, `round18/theory.md`, `round18/census.md`, `round18/lift.md` (with `referee_*.md`) |
| round 19: Thm 2.4, Thm 3.3, non-unit residue | `round19/README.md`, `round19/theory.md`, `round19/search.md` (with `referee_*.md`) |
| round 20: $(L_5^\circ)$, $V_E(5)$, Thm 4.1, Lemma K1 | `round20/README.md`, `round20/gauss.md`, `round20/kill_gauss.md` |
| round 21 and round 31 Parts VIII–IX: low-height window (capped for $\delta<1/6$ for inputs up to character-tested prime residues; open otherwise) | `round21/README.md`, `round31/README_VIII.md`, `round31/README_IX.md` |
| Theorem 1.3, orientation (projection) theorem | `round31/README_IX.md` §§1.2, 1.9(a), 3(a); `round31/ref9_proj.md` |
| Theorem 1.2, $M\le1+\omega(N')+K_C$ | `round31/README_VIII.md` §§1.1, 1.8; `round31/ref8_rank.md`, `ref8_extend.md` |
| fake circles and capped classes | `round4/outside.md`, `round5/construct.md`, `round8/lw.md`, `round10/analytic.md` |
| mean-width Thm B, Cor C3 (auxiliary forms capped) | `threshold/construct.md` |
| round 26: zero estimates vs route 5's ZE; spectral LP cap; target ZE$^{\rm all}_H$ | `round26/README.md`, `round26/ze.md` (with `referee_ze.md`) |
| round 27: dP2 calibration, uniqueness of extremal measures, split cones, multi-shift criterion | `round27/README.md` (with `kill_*.md`) |
| round 28: pre-registered test of the multi-shift criterion (unfavourable) | `round28/README.md` (with `verify.md`) |
| round 29: function-field model, Li–Liu template, residue-class twin, 18 screened ideas (no change) | `round29/README.md` (with `kill_*.md`) |
| round 30: function-field model decided as far as internal work goes (rank open, line closed; no change) | `round30/README.md` (with `kill_*.md`) |
| round 31: outside colleague, ten refereed passes (index-one reduction, fixed-plane exponents, integral five-mark models, norm defect and local isotropy, rigidity, random parents; **function-field uniform bound $M(D,0)\le24$ in characteristic 0; unbounded over $\mathbb F_2(t)$ (also squarefree) and over $\mathbb F_3(t)$ with $A^2+B^2$**, **unbounded at prime-exponent cap two, constant and finitely ramified models**; filters F3, F3+, F3++, F3+++; **characteristic-0 function-field $A_5$ open, true on sparse presentations**; $A_4(\delta)$ false for $\delta<1/6$; fixed-field unit-sum theorem; **on integer circles $M\le1+\omega(N')+K_C$ and fakes inside the low-height window**; **$(L_8)$ implies Rudin's conjecture**; random-parent rate at the smooth scale, same-norm separation of fakes, $(L_8^{pc})$ correction; sign-free residue tables collapse onto the genuine two-class problem; finite rigidity and a partial P3-all$^*$ fake, residue route closed; no change to the growth rate for general $R$) | `round31/README.md`, `README_II.md`, `README_III.md`, `README_IV.md`, `README_IVL.md`, `README_V.md`, `README_VL.md`, `README_VI.md`, `README_VIL.md`, `README_VII.md`, `README_VIIL.md`, `README_VIII.md`, `README_VIIIL.md`, `README_IX.md`, `README_IXL.md`, `README_X.md`, `README_XL.md` (with `ref*_*.md`, `langer_check.md`) |
| papers check: seven full texts on approximation constants | `papers_check/README.md` |

## Appendix: notes on the check

All 23 items of the adversarial check were applied. Three were applied in a modified form after re-reading the
sources.
* *Item 6 (unit class).* The check's suggested "equivalently $\zeta_a/\zeta_b=P/\bar P$ with $P\in\mathbb Z[i]$" is
  not an equivalence as stated: $i=(1+i)/(1-i)$ and $-1=i/(-i)$ are themselves of the form $P/\bar P$. The note
  therefore defines the class by equal units for a fixed choice of the $\pi_p$ (as in `growth/ptolemy.md` §1) and
  states the $P/\bar P$ form, with $P$ conjugate-primitive, only as a consequence.
* *Item 8 (labels).* Both remedies were used: the legend gains the label "proved by the killer", and Lemma K1 and
  the $k=3$ corollary carry it.
* *Item 17 ($\hat\gamma$, $D_{PQ}$).* The sentence was dropped in favour of the pointer to `gauss.md` Thm 4.2,
  Remark (3), as the check allowed.
