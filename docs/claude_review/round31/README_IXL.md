# Judge report, round 31: Part IX of the Question B thread (Astra, "The critical two-class problem on smooth parents")

Source: `astra_lemma_IX/answer.md` from line 1931 (`## Part IX.`) to the end (line 2149), with `final9.txt`,
`ix_checks.py` and the question `reply9.txt`.
Referee reports in this directory, **all eight present**: `ref9L_pool.md` (IX.0 squarefree theorem, corollaries,
polynomial cutoff), `ref9L_powers.md` (IX.0 bounded powers, twin and Gaussian), `ref9L_samenorm.md` (same-norm P2
comparison; IX.6 first paragraph), `ref9L_additive.md` (IX.1, IX.3, IX.6 modulus criterion and density),
`ref9L_fallback.md` (IX.2; IX.6 equivalence with Cilleruelo–Granville), `ref9L_relax.md` (IX.4; the m = 43
computation), `ref9L_nonlinear.md` (IX.5; the Prouhet–Tarry–Escott regression; comparison of scales),
`ref9L_meaning.md` (novelty, meaning, next steps).
Record consulted: `expert_note/expert_note.md` §§1, 2, 5; `round31/README_VIIIL.md` (format and the previous rulings on
this thread); `round31/README_VIII.md` (the P2 fakes); `round31/README_III.md` (the earlier random-parent theorem III.9);
`round18/theory.md` §§1.3, 2.1, 5 and `round18/referee_theory.md` (sources of $(L_8^{pc})$);
`astra_free_III/iii_graph_notes.md` lines 128-158 (the parallel thread's unrefereed pool estimate).
Judge scripts (my own code; nothing imported from Astra or any referee), in `README_IXL/`:
* `judge9L_pool.py` (output `judge9L_pool.out`): the union-bound exponents of IX.0 for $J=13,\dots,60$; exact enumeration
  of the conditional law of the hidden block on toy pools; the inequality $\binom{K-s+k}{k}^{-1}\le(s/K)^k$ on all
  $K<120$; the integer triangle identity, the coprime factorisation $q=GA_{12}A_{23}A_{13}$ and the heavy-block count on
  4,846 admissible triples; the Gaussian triangle on 245 random triples.
* `judge9L_powers.py` (output `judge9L_powers.out`): the bounded-exponent Gaussian extraction on 120 random
  configurations with $E\in\{2,3\}$: 4,068 pair factorisations, the per-prime inequalities, the case analysis at the
  hidden primes, and recovery of the hidden prime set from the triangle equation on 12,440 ordered triples.
* `judge9L_pc.py` (output `judge9L_pc.out`): exhaustive local minimisation for the printed lattice form at $p=2$
  ($v_2(e_T)\le3$) and $p=3$ ($v_3(e_T)=1$); the 2- and 3-adic active sets of 20,000 random primitive 7-tuples; the
  referee's example for the corrected form.
* `judge9L_relax.py` (output `judge9L_relax.out`): the Paley relaxation at $m=11,19$, including the pair-level
  signed-residue test.

Claims of my own are labelled **[J]** and are unrefereed.

**Notation and definitions used throughout.**
* **Twin (critical two-class problem).** $q\ge2$; at every exact prime power $p^f\parallel q$ at most two residues
  $r_p,r'_p$ modulo $p^f$ are prescribed (they may coincide, or agree modulo a lower power, or be nonunits). A *mark*
  is an integer $y$ with $y\equiv r_p$ or $r'_p\pmod{p^f}$ for every $p^f\parallel q$. $B_{\max}(q,C)$ is the maximum,
  over all residue choices and all real intervals of length $C\sqrt q$, of the number of marks in the interval.
  $s=\omega(q)$, $W=\log q$.
* **Circle problem.** $M$ is the number of lattice points of $x^2+y^2=N$ on an arc of length $C\sqrt R=CN^{1/4}$;
  $N'$ is the norm after dividing the cluster by its Gaussian gcd. The *growth rate* is the worst-case order of
  $M_{\max}(R)$ as $R\to\infty$; Theorem 1.1 of the record gives $(2/3+\epsilon)\log R/\log\log R$.
* **Pool.** A set $\mathcal P$ of $K$ distinct primes in $[P,2P]$; a *parent* is a uniformly random $s$-subset $S$;
  $J=\log(K/s)$.
* **Signature.** At each $p^f\parallel q$ label the two classes by bits; the signature of a mark is its bit vector in
  $\mathbb F_2^{s}$; $A$ is the set of signatures and $D=|A+A|$ (including 0). On circles the unary signature has bits
  $\mathbf 1_{a_{ip}\ge j}$, $1\le j\le v_p(N')$.
* **Fakes, P2.** The Legendre angle models of `README_VIII.md` §§1.3-1.5: exponent cap two, $\asymp s/2$ marks, all
  inputs of the list $A_{2M}$ (monomial gaps, rational characters, point residues at odd prime powers $\le2M$).
* **$(L_8^{pc})$** is the lattice form printed in `expert_note.md` lines 214-222; $\mathrm{Bal}^{pc}$ is the set of points of
  $M_{0,8}(\mathbb Q)$ whose 119 visible contact integers $d_T$ are pairwise coprime (`round18/theory.md` §1.3, Prop 2.1).
* **Conjecture 9 (Cilleruelo–Granville, arXiv math/0608109v1):** uniformly in $a$ and $b$, every interval of length
  $\sqrt b$ contains $O(1)$ roots of $x^2\equiv a\pmod b$. Its *rate form*: $o(\log b/\log\log b)$, or $o(\omega(b))$, roots.

## 0. Verdict (five lines)

1. Nothing in Part IX is false and nothing of Astra's mathematics is withdrawn. Eight referees returned no FALSE verdict on a claim of hers; the one FALSE (meaning) is on a proof route she did not claim, and I uphold it in the range that matters. I checked myself [J] the three items prioritised: the probability argument of IX.0 (conditional law, union bound, simultaneity), the bounded-exponent Gaussian extraction, and the vacuity of the printed $(L_8^{pc})$.
2. **The random-parent rate is a theorem**, for the twin and for circles, uniformly over all intervals, arcs and residue choices: for a random $s$-subset of any pool of $K\ge e^{40}s$ primes in $[P,2P]$, outside a fraction $e^{-sJ/20}$, $B_{\max}\le b-1$ and $M_{\max}\le4(b-1)$ with $b-1=\max(2,\lceil100\log(2P)/J\rceil)$; this is $O(\log s)$ at the smooth scale $P\asymp s\log s$ and $O(1)$ at $P=s^{1+\epsilon}$; bounded exponents are allowed. It is an average over a density-zero family. The twin statement is new; on circles it adds the endpoint $D=1$ to III.9 (the estimate itself was in the parallel thread's unrefereed notes).
3. **The same-norm comparison is established** (asymptotically, on typical parents of a narrow prime bin): on the same norms $N=\prod p^2$ the P2 fakes have $\asymp s/2$ marks while every genuine critical arc has $O(1)$. The input list $A_{2M}$ therefore cannot certify $M=o(s)$ even where the truth is bounded. Typical dyadic parents and worst-case parents are untouched.
4. **The printed lattice fallback $(L_8^{pc})$ is vacuous** (empty input, at the prime 2 and again at 3), and so is $\mathrm{Bal}^{pc}$ (every point has exactly five even contact integers). The Question A box, fact 3 of §2, the sandwich and the conclusion of the abc criterion must be restated with pairwise coprimality outside $\{2,3\}$ and a bounded $\{2,3\}$-part; that corrected form follows from $(L_8)$ and is open. The critical twin is exactly Cilleruelo–Granville's Conjecture 9 (uniform and little-oh scales).
5. **For the integer problem nothing moves**: the growth rate $(2/3+\epsilon)\log R/\log\log R$ for worst-case $R$, the uniform bound and Questions A and B are unchanged. The deterministic results (no affine 3-flat, $B\ll\sqrt D$, $n^{7/8}$) kill linear/characteristic-2-type profiles over $\mathbb Z$ and no recorded fake; the relaxation IX.4 is a barrier for valuation-based arguments (sharp at the larger-sieve constant $\tfrac12$) but not for signed pair residues. Score on the goal: 2.5/10.

## 1. Per block: referee-corrected statements

### 1.1 IX.0, the squarefree random-parent theorem (referee: pool; probability argument rechecked by the judge)

**Proved (statement as corrected; quoted in full).**
> **Theorem IX.0 (random parents; almost-all over a prime pool, not a worst-case statement).** Fix $C>0$. There is
> $s_0(C)$ such that for all $s\ge s_0(C)$, all $P$ and every set $\mathcal P$ of $K\ge e^{40}s$ primes in $[P,2P]$ the
> following holds, with $J=\log(K/s)$ and $b=\max(3,1+\lceil100\log(2P)/J\rceil)$. Let $S$ be a uniformly random
> $s$-element subset of $\mathcal P$ and $q=\prod_{p\in S}p$.
> 1. *(Twin.)* Outside an event of probability at most $\exp(-sJ/20)$: for every real interval of length at most
>    $C\sqrt q$ and every choice of at most two residues modulo each $p\in S$ (equal and nonunit residues allowed), at most
>    $b-1$ integers of the interval lie in a prescribed class modulo every $p\in S$.
> 2. *(Circle.)* If $\mathcal P$ consists of primes $\equiv1\pmod4$, outside a second event of probability at most
>    $\exp(-sJ/20)$: every arc of length $Cq^{1/4}$ of $x^2+y^2=q$ carries at most $4(b-1)$ lattice points; the same
>    holds for clusters on $x^2+y^2=q|\Gamma|^2$ all of whose points are divisible by one Gaussian integer $\Gamma$.
>
> The residues, the interval and the arc are chosen after $S$. The bounds $b-1$ and $4(b-1)$ are absolute: $C$ enters
> only through $s_0(C)$ (one may take $s_0(C)=\max\{(3.2+9\log C_0)/2.24,\ (\log C)/3.7,\ 1\}$ with $C_0=\max(2,C,C^2)$,
> pool referee). Not covered: norms whose Gaussian-gcd-normalised norm is not a sampled $q$ (for example the cap-two
> norms of P2, for which see §1.2).

**Corollaries (Proved).**
* *Smooth critical scale.* Pool = all primes (all split primes) in $[P,2P]$, $P=As\log s$ with $A\ge e^{40}(1+o(1))$
  ($A\ge2e^{40}(1+o(1))$ for split primes): $J\to\log A$ (resp. $\log(A/2)$), $\log q\sim s\log s$, so
  $\log q/\log\log q\sim s$, and $B_{\max},M_{\max}=O(\log s)=O(\log\log q)$ outside a fraction $e^{-\Omega(s)}$.
  Purely asymptotic: $\log s\gg40$ is needed before "$\log q\sim s\log s$" is accurate.
* *Polynomial scale.* $P=s^{1+\epsilon}$: $J\sim\epsilon\log s$, $b-1\le2+100(1+\epsilon)/\epsilon+o(1)$.
* *Primes below a polynomial cutoff.* For fixed $A>1$ and any pool of $K$ primes below $X\le s^A$ (no lower bound on the
  primes) with $J=\log(K/s)\ge40A$: outside a fraction $\exp(-\Omega_A(sJ))$ the bound is $O_{A,C}(1+\log X/J)$ in either
  problem (heavy block $\ge s/(5A)$ from $W\ge\log s!$; exponent per $sJ/A$ is $-0.106$). A cutoff $X=A's\log s$ needs
  $\log A'\ge80$.
* **[J] Corollary for square roots (Proved, immediate).** For squarefree $q$, the roots of $x^2\equiv a\pmod q$ occupy at
  most two classes modulo each $p\mid q$ (one if $p\mid a$ or $p=2$). Hence, outside the event of item 1, **for every
  $a$**, every interval of length $C\sqrt q$ contains at most $b-1$ roots: almost every parent of the pool satisfies the
  rate form of Conjecture 9 ($O(\log\omega(q))$ at the smooth scale) and, at $P=s^{1+\epsilon}$, its bounded form. With
  exponents (§1.2) this holds for $a$ coprime to $q$ only (for $p\mid a$, $x^2\equiv a\pmod{p^e}$ can have many classes).

**[J] Judge check of the probability argument.** The points that carry the weight:
1. *Deterministic part (holds for every parent).* For $b$ marks, $h_{ij}=\prod\{p:y_i\equiv y_j\}$ divides $y_j-y_i$, so
   $h_{ij}\le C\sqrt q$ and $t_{ij}=(y_j-y_i)/h_{ij}\in\mathbb Z\setminus\{0\}$. A prime separates $k(b-k)\le b^2/4$ pairs, so
   mean $\log|t_{ij}|\le\log C+W/(2(b-1))$; each pair lies in $b-2$ triples, so some triple has
   $\log|t_{12}t_{23}t_{13}|\le3\log C+3W/(2(b-1))$, and each $|t|\le$ that product since $|t|\ge1$. At a prime not
   constant on the triple, "at most two residues" forces exactly one pair to agree, giving $q=GA_{12}A_{23}A_{13}$
   (pairwise coprime) and (IX.0a) $t_{12}A_{12}+t_{23}A_{23}=t_{13}A_{13}$. The primes separating marks 1 and 2 have
   log-mass $\ge W/2-\log C$ and each lies in $A_{13}$ or $A_{23}$, so one block has $\ge(W/4-\tfrac12\log C)/\log(2P)\ge s/5$
   primes once $P>16$ (automatic: $K<P$ forces $P>e^{40}s$) and $s\ge s_0(C)$. Computed [J]: 4,846 admissible triples,
   0 failures of the identity, the factorisation or the block count.
2. *Conditional law.* With the sample in uniformly random order, $(p_1,\dots,p_s)$ is a uniform injective sequence.
   For a labelling $\lambda$ of slots by $\{G,12,23,13\}$ and a forced block $F$ with $k$ slots, both fixed before
   sampling, conditioning on the $s-k$ exposed slots leaves the $F$-slots a uniform injective $k$-tuple from the
   $K-s+k$ unused primes. Given $t$, (IX.0a) fixes the product over $F$ as one number, hence (unique factorisation) at most
   one $k$-set: probability $\le\binom{K-s+k}{k}^{-1}\le(s/K)^k\le e^{-sJ/5}$. Computed [J]: exact enumeration on four toy
   pools gives the maximum conditional probability equal to $\binom{K-s+k}{k}^{-1}$; the inequality with $(s/K)^k$ has
   0 violations for all $K<120$. A sorted order would break step 2 (the hidden slots would be confined to an interval of
   the pool); Astra uses a random order, correctly.
3. *Union bound.* Events are indexed by $(\lambda,F,t)$: $4^s\cdot3\cdot(2Y)^3$ in the twin, $7^s\cdot3\cdot(2Y)^3$ on circles
   ($7=1+3\cdot2$: a $G$-slot needs no orientation because $G$ cancels), $\log Y=3\log C_0+3s\log(2P)/(2(b-1))$. The log is
   $\le O(\log C_0)+s\log4\ (\text{or }\log7)+9s\log(2P)/(2(b-1))-sJ/5$. Per $sJ$ at $J=40$: $-0.120$ (twin), $-0.106$
   (circle), against $-0.05$ needed. Computed [J]: the 7-label case still clears $-0.05$ at $J=18.6$, fails at $J=13$;
   the 4-label case clears it at $J\ge18.6$ and fails at $J=13$. So $J\ge40$ is safe and not sharp.
4. *Simultaneity.* The event $E(\lambda,F,t)$ is a property of the sampled primes alone. Intervals, residues, the triple,
   its labelling and $t$ are existential witnesses and are never enumerated. Hence the exceptional set is one set of
   parents serving all intervals of length $\le C\sqrt q$ and all residue choices at once. The same holds on circles for all
   arcs and unit classes (one unit class retained, factor 4).
5. *Gaussian triangle.* With fixed generators $\pi_p$, a point is $\eta\prod g_{ip}$, $g_{ip}\in\{\pi_p,\bar\pi_p\}$; for two
   points in one unit class $z_i=\eta H_{ij}X$, $z_j=\eta H_{ij}\bar X$, so $z_j-z_i=2i\eta H_{ij}t_{ij}$ with
   $t_{ij}=-\operatorname{Im}X\in\mathbb Z\setminus\{0\}$; the forced Gaussian product has norm equal to the rational product of
   the hidden primes. Computed [J]: 245 random triples, 0 failures (identity, integrality, triangle in $\mathbb Z[i]$).

**Computed (pool referee).** Exact maximal twin counts for $s=6..10$: $B_{\max}=3$ or $4$ (every extremal 4-set is a
2-cube $\{0,q_1,q_2,q_1+q_2\}$, $q_1q_2=q$); circles $M_{\max}\le3$. The hypotheses fail there by about $10^{16}$; these
data test nothing asymptotic. **Heuristic (pool referee, upheld):** a first-moment count predicts $O(1)$ already at
$P=As\log s$, so the factor $\log s$ is very probably an artefact of the entropy count.

**Novelty.** Twin: new (the record had no random-parent statement for the two-class problem). Circles: the endpoint
$P\asymp s\log s$ is new relative to the refereed III.9 (whose $M_D$ blows up as $D\downarrow1$), but the estimate
$O(1+\log P/\log(K/r))$ and the exposure lemma are in `astra_free_III/iii_graph_notes.md` lines 131-135 and 155, proved there
by quartets with a Jensen step and never refereed [J: I read lines 128-158; the statement is there as described]. Astra
cites these notes. New in method: three marks instead of four, and a heavy block that is deterministic on every triple.

**Verdict.** Correct. Gap: none. Withdrawn: none.

### 1.2 IX.0, bounded exponents (referee: powers; Gaussian extraction rechecked by the judge)

**Proved (statements as corrected; quoted in full; constants the referee's).**
> **Theorem T (twin, bounded powers).** Pool of $K$ primes in $[P,2P]$, $s$ drawn uniformly, $J=\log(K/s)$, fixed $C>0$ and
> integer $E\ge1$, $J\ge40E\log(1+3E)$. There is $s_0(C,E)$ such that for $s\ge s_0$, outside a set of samples of
> probability at most $\exp(-sJ/(8E))$: for all exponents $1\le e_p\le E$, all choices of two classes modulo each
> $p^{e_p}$ and all intervals of length $\le C\sqrt q$, $q=\prod p^{e_p}$, all chosen after $S$, the interval contains at most
> $b-1=\max(2,\lceil100E^2\log(2P)/J\rceil)$ admissible integers.
>
> **Theorem G (circles, bounded powers).** Same pool, now of split primes, $J\ge400E^2\log(E(E+1)^3)$. There is $s_0(C,E)$
> such that for $s\ge s_0$, outside probability $\exp(-sJ/(250E^2))$: for every $N=\prod_{p\in S}p^{e_p}$ with $1\le e_p\le E$
> chosen after the sample, every arc of length $CN^{1/4}$ of $x^2+y^2=N$ carries at most $4(b-1)$ lattice points, where
> $b-1=\max(216E^2,\lceil90000E^5\log(2P)/J\rceil)$.

Nothing is uniform in $E$: $J_0(E)\asymp E^2\log E$, exceptional exponent $\asymp sJ/E^2$, bound $O(E^5(1+\log P/J))$.
Corrections (wording): $\mathbb E D_{12}\ge W/3$ needs $b\ge4$ (supplied by $b\ge216E^2$); the label counts $1+3E$ and
$[E(E+1)^3]^s$ are overcounts ($(E+1)^3$ per slot suffices); Astra's "$\le\log C_0+W/(4(b-1))$" can be sharpened to
$\log(C/2)+W/(4(b-1))$.

**[J] Judge check of the Gaussian extraction (the step Astra asked to have prioritised).**
1. *Pair factorisation.* With $z=\eta\prod\pi_p^{a_p}\bar\pi_p^{e_p-a_p}$ and $H_{ij}=\prod\pi_p^{\min}\bar\pi_p^{e_p-\max}$, the
   two quotients $u_i=z_i/(\eta H_{ij})$, $u_j$ have at each $p$ the factors $\bar\pi^{d}$ and $\pi^{d}$ ($d=|a_{ip}-a_{jp}|$) in some
   order, so $u_j=\bar u_i$ exactly, $t_{ij}=-\operatorname{Im}u_i\in\mathbb Z\setminus\{0\}$ and $\log|H_{ij}|=(W-D_{ij})/2$. The chord bound
   gives $D_{ij}\ge W/2-2\log(C/2)+2\log|t_{ij}|$; threshold cuts give mean $D\le W/2+W/(2(b-1))$.
2. *Extraction.* Three rows drawn independently with replacement: $\mathbb E D_{12}\le\sum_pw_p\delta_p$ (because
   $|a_1-a_2|\le e_p\mathbf 1_{a_1\ne a_2}$), $\mathbb E D_{12}\ge(1-1/b)(W/2-O_C(1))\ge W/3$ for $b\ge4$; by Cauchy–Schwarz
   with $\sum_{a=0}^{e_p}(1-\alpha_a)=e_p$, $\Pr(a_1=a_2\ne a_3)=\sum\alpha_a^2(1-\alpha_a)\ge\delta_p^2/E$; weighted
   Cauchy–Schwarz gives expected exclusive weight $\ge W/(9E)$ and count $\ge s/(18E^2)$; repetitions cost $\le3s/b$; for
   $b\ge216E^2$ the mean over distinct ordered triples is $\ge s/(24E^2)\ge2c_Es$, $c_E=1/(100E^2)$, so a fraction $\ge c_E$
   of distinct triples have $\ge c_Es$ exclusive primes, and (nonnegativity of $\log|t|$, Markov on that fraction) one of
   them has $\sum\log|t|\le3(\log C_0+W/(4(b-1)))/c_E$. All directions are right.
3. *Division and hiding.* Dividing by $G_3$ ($\pi$-exponent $\min_ia_{ip}$, $\bar\pi$-exponent $e_p-\max_ia_{ip}$): at
   $a_1=a_2=x\ne a_3=y$ the residual factor sits in the $t_{12}$-coefficient only, $\bar\pi^{y-x}$ if $y>x$ and $\pi^{x-y}$ if
   $y<x$. Three-level slots put factors in two coefficients but are exposed and known from $(e_p,a_{1p},a_{2p},a_{3p})$. The
   triangle equation then fixes the hidden Gaussian product $X=(t_{13}R_{13}-t_{23}R_{23})/(t_{12}U)$, whose norm
   $\prod p^{|x_p-y_p|}$ has every exponent $\ge1$, so it determines the hidden prime set; the conditional probability is
   $\le\binom{K-s+k}{k}^{-1}$ with $k\ge c_Es$. The order of choices is $E,C\to J_0(E)\to b(E,P,J)\to s_0(C,E)$.
4. Computed [J] (`judge9L_powers.out`): 120 configurations ($E=2,3$, 2-4 primes, 4-9 marks, random units): 4,068 pair
   factorisations, 0 failures; the per-prime identities $\Pr(a_1=a_2\ne a_3)=\sum\alpha^2(1-\alpha)$ (direct triple count),
   $\ge\delta^2/e_p$, $\mathbb E|a_1-a_2|\le e_p\delta$: 0 failures; 21,006 ordered triples, of which 12,440 have hidden slots:
   the case analysis at every hidden prime, the triangle identity in $\mathbb Z[i]$, and exact recovery of the hidden set from
   the norm of $X$: 0 failures.

**Verdict.** Correct (Theorems T and G as stated). Gap: none. Withdrawn: none.

### 1.3 Same-norm comparison and the IX.6 Legendre paragraph (referee: samenorm)

**Proved, conditional only on refereed results** (VIII.3 of the parallel thread, refereed in Part VIII; Theorem G, §1.2).
Statement as corrected, quoted in full:
> **Theorem (same-norm comparison).** Fix $\rho\in(0.99,1)$. There is $l_0(\rho)$ such that for every prime $l\ge l_0$,
> $l\equiv3\pmod4$, with $s=l-1$ and $P=l^{1.99}$, some interval $I_l=[P_b,P_be^{\sigma}]\subset[P,2P]$ with
> $\sigma\le l^{-3/4}$ contains $K\ge(1+o(1))\,l^{1.24}/(2\cdot1.99\cdot\log2\cdot\log l)$ primes $\equiv1\pmod4$, and:
> (i) *(fake half, every subset)* for **every** $s$-set $S$ of split primes of $I_l$ and $N=\prod_{p\in S}p^2$, the Legendre
> model of VIII.3 (rows $1+K_i$, the primes of $S$ placed on the columns by a suitable permutation) is an angle model of
> normalised norm exactly $N$ with $(\rho/2+o(1))l$ marks on an arc of constant one, satisfying both monomial gaps for every
> nonzero integral exponent vector (hence all rational characters of the profile), axiom (R) at every odd prime power
> $\le2m$ for all four units, positive Haar measure of good angles, $P_+(N)\le2P$ and $W=(3.98+o(1))l\log l$;
> (ii) *(genuine half, almost every subset)* outside a fraction $\exp(-(0.24/1000+o(1))\,l\log l)$ of the $s$-subsets of the
> split primes of $I_l$, for every $C>0$ every arc of length $C\sqrt R$ ($R^2=N$) of $x^2+y^2=N$ carries at most
> $\lceil C\rceil M_1$ lattice points, $M_1$ absolute.
> Hence there are arbitrarily large norms $N=R^2$, $R$ squarefree, on which (i) and (ii) both hold.

Corrections: the fake half holds for every subset (it must be in the statement); purely asymptotic (no bin holds $s$ primes
for $l\le10^7$, computed; Theorem G's threshold $J\ge1600\log54\approx6383$ with $J\approx0.24\log l$ needs
$\log l\gtrsim2.66\cdot10^4$); the norms do not depend on $C$ (subdivision). [J] Minor: with the constants of §1.2 at $E=2$,
$M_1=4(b-1)\approx4\cdot90000\cdot32\cdot8.3\approx10^8$, not "of order $10^7$". **Not checked:** Astra's remark that the older
quartet exposure would also do the genuine half (III.9 is stated for dyadic pools only); it is not needed.

What the widening needs, checked by the referee and upheld: of all uses of the prime set in VIII.3 only the pair margin
$w(K_i-K_j)-W/2\ge2\log P-o(1)$ depends on the bin width; its deviation under a uniform permutation is
$O(\sigma\sqrt{l\log l})=O(l^{-1/4}\sqrt{\log l})$ (Azuma under the transposition coupling, union over $<l^2$ pairs); VIII.3
tolerates up to $(1.96-\epsilon)\log l$. The normalisation (both exponents 0 and 2 in every column of the selected rows) is
deterministic.

**IX.6 first paragraph (Proved; computed for $q=43..383$).** The P2 parity is the singleton indicator of the owner; so P2 is
not a Boolean function of $d$ affine forms on $2^r=O(s)$ labels (Reed–Muller dimension $\sum_{j\le d}\binom rj<m$); stronger,
any such labelling needs $n\ge\exp(c_dm^{1/d})$; P2 signatures are Sidon, $D=1+\binom m2$. "Evades" means the tests give
bounds of order $m$, not that P2 violates them.

**Verdict.** Correct. Gap: none. Withdrawn: none (one remark not checked).

### 1.4 IX.1, IX.3, IX.6 criterion and density (referee: additive)

**Proved.**
* **IX.1.** For $q>\max(C^2,(4C)^6)$ the twin signature set contains no affine 3-dimensional subspace of
  $\mathbb F_2^{\omega(q)}$; eight marks with signatures forming an affine 3-space force $L\ge q^{2/3}/4$. The exponent $2/3$ is
  sharp ($y_x=\sum x_i\prod_{j\ne i}Q_j$ has span $\approx3q^{2/3}$). The proof is characteristic-zero (Walsh inversion over
  $\mathbb Q$), which is why the characteristic-2 Sylvester–Hadamard families have no counterpart over $\mathbb Z$. New; low.
* **IX.1 sumset estimate** $B(B-1)\le(D-1)(1+\sqrt{1+8(D-1)})$: correct, superseded by IX.3. "Such bounds transfer to the
  log-height rate since $s!\le q$" has content only with an a priori bound on $D$.
* **IX.3 twin.** For every $h\ne0$ at most $(1+\lfloor2L/K\rfloor)(1+\lfloor(L/H)^2/K\rfloor)\le(1+2C^2)(1+C^2)$ unordered pairs
  have signature difference $h$ ($H$ the product of the prime powers where $h=0$, $K=q/H$); hence $B\le1+\sqrt{2c_C(D-1)}$.
  Read correctly this says "$A$ is a $c_C$-fold Sidon set"; it bounds $B$ by $o(s)$ only if $D=o(s^2)$ is known a priori.
  New (twin analogue of the distinct-label lemma).
* **IX.3 circle.** $M\ll_C\sqrt D$ for unary signatures, with subarcs of angular diameter $<2N''^{-1/4}$ ($c<2$ suffices;
  Astra's $c<1/10$ is safe). **Not new in substance**: an immediate corollary of the refereed distinct-label lemma
  (`README_II.md`, `ref2_tri.md` §1.1). Holds in angle models; saturated by the fakes (Sidon).
* **IX.6 twin criterion** $B\ll_C1+W/\log P_*$, $P_*=\max_{p^f\parallel q}p^f$: correct, new as a twin statement, routine;
  on comparable-prime smooth parents it gives the larger-sieve order $s$.
* **IX.6 circle criterion** $M\ll_C1+W/\max_pe_p\log p$: correct, **already in the record** with better constants
  (`round8/stability.md` Thm 3.3(2), Cor 3.2), and dominated by $P_+(N')\ge(2-o(1))M\log M$.
* **IX.6 density** $B_{\max}(q)/\omega(q)\to0$ in natural density: correct and trivial (in fact $B_{\max}$ is bounded in
  probability, Dickman); it concerns typical, non-smooth integers.

**Computed (referee).** Twin, exhaustive for $q\le400$ ($C=1$) and selected $q$ up to $720720$; circles, all $N\le10^6$ with
split primes only and structured $N$: 0 violations of every proved inequality; the circle constant matters (violations at
$c=8$, none at $c\le4$).

**Verdict.** Correct. Gap: none. Withdrawn: none. Two items not new (IX.3 circle, IX.6 circle).

### 1.5 IX.2 and the Cilleruelo–Granville equivalence (referee: fallback; vacuity rechecked by the judge)

**Proved (IX.2, strengthened).** The lattice form $(L_8^{pc})$ printed in `expert_note.md` lines 216-222 (copied from
`round18/theory.md` §5(2); I compared the wording, [J]) has **empty input** at every scale, so it is a theorem with $Z=\emptyset$.
The set $\mathrm{Bal}^{pc}$ is **empty**. Hence the Bal-form $(L_8^{pc})$ of `round18/theory.md` §2.1, the "pc" restriction in
the sandwich, Prop 2.1(4), census Prop 1.3(2) "in particular", fact 3 of §2 and the conclusion of the abc criterion
(theory Prop 2.6) are all vacuously true.

**[J] Judge check.**
1. *Exact wording.* The note requires: $e_T$ ($3\le|T|\le6$) **pairwise coprime integers**; $y\in\Lambda_e$ primitive,
   $y\notin\mathbb Z\mathbf1$, $[y]\in M_{0,8}\setminus Z$ (distinct entries); the 21 pair cofactors
   $c_{ij}=(y_i-y_j)/\prod_{T\supseteq\{i,j\}}e_T$ "pairwise coprime of size $e^{w(1\pm\eta_1)}$". No prime is exempt and
   nothing is said "up to bounded factors". Pairwise coprimality of the $e_T$ makes each $c_{ij}$ an integer.
2. *The prime 2.* At most one $e_T$ is even. If none: seven integers have $\ge\binom43+\binom32=9$ equal-parity pairs, all with
   odd denominators, so $\ge9$ even $c_{ij}$. If $e_T$ is even with $v=v_2(e_T)$: for $|T|\ge4$, four $T$-marks in one class
   mod $2^v$ give $\ge2$ pairs congruent mod $2^{v+1}$, whose denominators have 2-adic valuation exactly $v$; for $|T|=3$ the
   four outside marks give $\ge2$ equal-parity pairs with odd denominators. So $\ge2$ even cofactors: not pairwise coprime.
   Computed [J] (exhaustive, $v_2(e_T)\le3$): the minimum number of even cofactors is 9 (no even block), 7 ($|T|=3$),
   5 ($|T|=4,5$), 6 ($|T|=6$). At $p=3$ ($v_3(e_T)=1$): minimum 5 (no block), 2 ($|T|=3,4,5$), 3 ($|T|=6$). So pairwise
   coprimality fails at 3 independently.
3. *$\mathrm{Bal}^{pc}$.* For the primitive representative $x$ (point 0 at $\infty$), the sets $T$ that are residue classes of
   $x\bmod2^t$ for some $t\ge1$ form a full binary tree on 7 leaves (each class splits into at most two lifts; the root is a
   class only at $t=0$ by primitivity), with 6 internal nodes, hence exactly 5 non-root internal nodes, all visible
   ($2\le|T|\le6$) and all with $2\mid d_T$. Computed [J]: 20,000 random primitive 7-tuples, always exactly 5 active sets at 2
   (at least two of size $\ge3$); at 3 between 2 and 5. So no point has pairwise coprime contact integers.
4. *Is there a non-empty intended reading?* Not in the text: the source's own motivation ("the tree at $p$ is a single edge")
   is impossible at 2 and 3 (a single-edge tree needs $p\ge5$, or $p\ge7$ when $|T|\in\{2,6\}$). The intended statement
   is the corrected one below.

**Proved (corrected fallback; referee's wording, upheld).**
> **$(L_8^{pc})'$, lattice form (corrected).** Let $(c,\eta_0,Z,C_0)$ be as in $(L_8)$ (anchor $0$, unsymmetrised $Z$). For
> every $B_0\ge1$ there is $w_1=w_1(B_0)$ with the following property. Let $e_T\in[e^{(1-\eta_1)w},e^{(1+\eta_1)w}]$
> ($3\le|T|\le6$, $\eta_1=\eta_0/2$) be integers that are pairwise coprime outside the primes $2,3$, and let $w\ge w_1$.
> Then $\Lambda_e$ has no vector $y$ with distinct entries such that $[y]\notin Z$, the pair cofactors
> $c_{ij}=(y_i-y_j)/\prod_{T\supseteq\{i,j\}}e_T$ are integers that are pairwise coprime outside $2,3$, of size
> $e^{w(1\pm\eta_1)}$, and the $\{2,3\}$-part of $\prod_{i<j}(y_i-y_j)$ is at most $B_0$.

$(L_8)\Rightarrow(L_8^{pc})'$ is proved (blocks = $\{2,3\}$-free parts, residues = $\{2,3\}$-parts, $H_\Sigma\le\log B_0$; at
$p\ge5$ the divisible sets are nested or disjoint, the crossing case would force three divisible cofactors;
$w_1=\max(2\log B_0/\eta_0,(\log B_0+C_0)/c)+1$). Its hypotheses are satisfiable (Computed, referee and [J]:
$y=(0,5,10,107,267,69,-229)$, $e_{123}=5$, all other $e_T=1$). **Open:** $(L_8^{pc})'$ itself; whether it implies Rudin;
whether the abc criterion (Prop 2.6) survives on the corrected class (Heuristic: the gcd identity should hold up to
$B_0^{O(1)}$; unchecked).

**Proved (IX.6, equivalence with Conjecture 9; quantifiers as corrected).** With $B_C(q)$ the twin maximum and $N(b)$ the
maximum number of roots of $x^2\equiv a\pmod b$ in an interval of length $\sqrt b$:
(a) $B_C(q)\le(C+2)\,N(4q)$ (CRT integers $a,b$, map $y\mapsto2y-a-b$; no unit assumption); (b) $N(b)\le2B_1(Q_2)$ with
$Q_2\le b$, $\omega(Q_2)\le\omega(b)$ (valuation division of VIII.2). Hence Conjecture 9 $\Leftrightarrow$ $B_C$ bounded for one
$C$ $\Leftrightarrow$ for every $C$; and $N=o(f)\Leftrightarrow B_C=o(f)$ for $f=\log x/\log\log x$ and for $f=\omega$. So the twin's
$o(s)$ question is exactly the $o(\omega)$ rate form of Conjecture 9. Computed (referee): 61,833 admissible $y$ ($q<3000$) and
326,907 square classes ($b\le1500$), 0 failures. **Literature: not checked beyond snippets**; no result improving the trivial
order $\omega(b)$ at length $\sqrt b$ was found (bounded search).

**Verdict.** Correct, and stronger than claimed ($\mathrm{Bal}^{pc}=\emptyset$). This is a correction to the record, not progress.

### 1.6 IX.4, the squarefree relaxation (referee: relax; spot-checked by the judge)

**Proved (statement as corrected).** For primes $m\equiv3\pmod8$, with $P=c_0m(\log m)^2$ ($c_0$ a fixed absolute constant), there
are $m$ pairwise coprime semiprimes $Q_v$ of split primes in $[P,2P]$, in one multiplicative bin of ratio $e^{1/\lceil m\log m\rceil}$
(no prime-tuple hypothesis), $q=\prod Q_v$, $s=2m$, $W=(2+o(1))m\log m$, so $m=(\tfrac12+O(\log\log m/\log m))\log q/\log\log q$;
and an edge table $d_{ij}=h_{ij}t_{ij}$ (Paley cuts $S=K-I$, $t_{ij}=w_j-w_i$, $w_i=ki+\mathbf 1_{i=m-1}$, $k=\lfloor P/(4m)\rfloor$) with:
(a) $0<|d_{ij}|<\sqrt q$; (b) $\gcd(d_{ij},q)=h_{ij}$, two-class cuts at every prime, minimal agreement; (c) integral reduced
cofactors of mean log size $\log m+O(\log\log m)$ forming an exact cocycle; (d) at every prime and depth the valuation
(collision-partition) pattern of a genuine local configuration, simultaneously at finitely many places; (e) Sidon signatures,
$D=1+\binom m2$; (f) an ordered real configuration matching every edge to relative error $O(1/\log m)$. It violates:
(g) the signed residues of single edges modulo every $p\mid q$ (a genuine configuration has $d_{ij}^2\bmod p$ constant over the
pairs separated at $p$; here it takes almost as many values as there are pairs); (h) signed consistency modulo small primes;
(i) the triangle identities over $\mathbb Z$. **No two-class configuration with these cuts exists** in any interval of length
$<q/(16(m+1)P^4)=q^{1-O(1/m)}$, for every residue choice (referee's strengthening of $C\sqrt q$).
**The constant is sharp for the larger sieve** (referee, Proved): when all primes of $q$ exceed $b$, the Vandermonde/pair-valuation
form of the larger sieve gives $b\le(\tfrac12+o(1))\log q/\log\log q$; the relaxation attains it to the leading constant
(asymptotically only: at $m=43$ the sieve still allows about 118 marks).

**Correction of interpretation (upheld, with a [J] precision).** "Isolates the signed triangle identity" is accurate but
misleading. For tables with $|d_{ij}|<q/3$, the triangle identities hold **iff** there are residues $\rho_i$ with
$d_{ij}\equiv\rho_j-\rho_i\pmod q$ [J: the referee's "for any residues $r$" should read "for some residues"; the substantive
direction is "signed pair residues with one consistent $\rho$ $\Rightarrow$ triangle identities"]. Already the residue-free pair
test "$d_{ij}^2\bmod p$ constant over separated pairs" fails maximally. Computed [J] (`judge9L_relax.out`): at $m=11$ and
$m=19$, Gram identity, $\binom m2$ distinct pair sums, $\gcd(d,q)=h$, all $\binom m3$ triangles fail, and at every $p\mid q$ the
squares $d_{ij}^2\bmod p$ take 28-30 (resp. 86-90) distinct values over 30 (resp. 90) separated pairs. So the relaxation is a
barrier for arguments that read the table through magnitudes, gcds, valuation partitions and signatures (the larger sieve in
Gallagher and Vandermonde forms), **not** for arguments using signed residues of single differences (IX.3 uses these).

**Computed (referee, reproduced).** $m=43$: 86 primes, 903 distinct pair sums, 12,341 failed triangles $=\binom{43}3$, i.e.
every triple; the count is generic and carries no information beyond the theorem.

**Verdict.** Correct, interpretation narrowed. Gap: none. Withdrawn: none.

### 1.7 IX.5, the PTE regression and the comparison of scales (referee: nonlinear)

**Proved (statement with the scope remark added).**
> Let $V=\mathbb F_2^r$, $n=2^r$. Prescribe at every exact prime power (twin), or at every exponent threshold of the
> Gaussian-gcd-normalised parent (circle), a bit that is a Boolean function of at most $d$ affine forms on $V$ (forms and
> functions may vary). For any set of labels carrying distinct marks: $B\ll_{C,d}n^{1-2^{1-2d}}$ in the twin and
> $M\ll_{C,d}n^{1-2^{1-2b_d}}$ on circles, $b_d=\lfloor2^{d-1}\sqrt d\rfloor$; for $d=2$ both are $\ll_Cn^{7/8}$. The circle bound
> holds in every angle model with the record's monomial gaps. Partial coverage: in the twin it suffices that the condition
> hold on prime powers of product $q_g\ge q^{1/2+\epsilon}$, with $t=\lceil d(1/2+\epsilon)/\epsilon\rceil$.
> *Scope.* Such a profile has binary signature rank $\le\sum_{j\le d}\binom{\log_2n}{j}$; the bound is $o(s)$ only if
> $n=o(s^{8/7})$ ($d=2$). Every recorded extremal fake (Paley, skew-Hadamard, Legendre P2) and the IX.4 relaxation have
> binary rank of order $m$, so for them the class needs $n\ge\exp(cm^{1/d})$ and the bound is vacuous.

The $p=2$ worry does not arise (unnormalised integer Walsh sums). Exact maxima of $\sum_u|b_{p,j,u}|$: 2, 6, 12 for
$d=2,3,4$ ($b_2$ sharp, $b_4=16$ not). **Heuristic (referee):** no natural reason for worst-case clusters to lie in the class.
**Computed (reproduced):** all PTE counts (3, 15, 6; 1, 1, 3; 14 and 1; 7,315 and 2,520), both inequalities $(2L)^2\ge q$,
$(4L)^3\ge q^2$ globally and per flat. Comparison of scales consistent with `ref8L_attack.md` (the 8- and 10-mark families are
at $L=q^{4/7+o(1)}$, $q^{3/5+o(1)}$, not at critical scale; in characteristic 2 Walsh inversion fails).

**Verdict.** Correct. Gap: none. New to the record; literature (Bonin–Qin type counting) seen in snippets only.

### 1.8 Meaning (referee: meaning)

Upheld, with two judge remarks.
* The random-parent rate is a rate **on average over a pool measure**, not "the uniform bound for almost all smooth parents":
  at the smooth scale the bound is $O(\log s)$; the exceptional fraction is $e^{-sJ/20}=R^{-J/(10\log P)}=R^{-o(1)}$ there; the
  theorem allows $\exp(s(1+J)-sJ/20)=(\#\text{parents})^{0.951}$ exceptional parents at $J=40$ (Computed; [J] rechecked:
  $(41-2)/41=0.951$); it does not cover the dense profile extremal for Theorem 1.1 ($K=s$, $J=0$), and any slot-exposure
  version needs $J>5\log4\approx6.9$.
* The input that the fakes lack is the exact additive identity with integer cofactors (twin: (IX.0a); circle:
  $z_j-z_i=2i\eta H_{ij}t_{ij}$ with the cocycle in $\mathbb Z[i]$). The record's classification "the $\sigma$-cocycle: capped
  ($A_0$)" (`round10/additive.md` §6; overview §5.4 rows 15-16) is correct only for the consequences it analysed (local
  consequences at primes not dividing $N$, triangle areas, abc); the identity itself is not capped.
* **FALSE verdict (on a route, not on Astra): removing exceptional parents by counting over exponent patterns and cofactor
  tables.** Upheld in the range that matters, with a [J] precision. The analyst's encoding count compares
  $(m-1)s\log2+W/2$ with the pool entropy $s(1+J)$ and uses $W/2\ge(s/2)\log P>sJ$. Since $K<P$ gives $J<\log(P/s)$, the
  last inequality holds when $P\le s^2$, which contains the smooth critical scale and every $P=s^{1+\epsilon}$ with
  $\epsilon<1$; for $P\gg s^2$ it is not shown. The conclusion is a statement about that counting method (for a fixed parent
  the union bound is "#clusters ≤ #candidates"), consistent with the record's closure of random-parent density statements
  as a route to the worst case (`README_III.md` §3.1). Label: **Proved for $P\le s^2$ as a statement about the method.**
* Why average and not worst case (Proved, as an analysis of the proof): randomness enters only where the forced hidden block
  costs $e^{-kJ}$; for a fixed parent the equation is a squarefree three-term $S$-unit equation $aX+bY=cZ$ with $|S|=s$ and
  small coefficients, on which abc is vacuous and the known upper counts are $\exp(O(s))$.

### 1.9 Withdrawn or corrected in Part IX (collected)

1. IX.0: the implied constant may be dropped ($b-1$ absolute; $C$ only in $s_0(C)$); circle coverage is norm $q$ and its
   common-Gaussian-factor multiples only; the circle rate at $P\asymp s\log s$ is credited to the parallel notes (unrefereed)
   with the first refereed proof here.
2. Bounded powers: explicit constants; $b\ge4$ needed in one step; nothing uniform in $E$.
3. Same-norm: fake half for every subset; purely asymptotic; $M_1\approx10^8$; the quartet alternative not checked.
4. IX.1: the "transfer to the log-height rate" remark needs an a priori bound on $D$. IX.3 circle and IX.6 circle: not new.
5. IX.2: stronger than stated ($\mathrm{Bal}^{pc}=\emptyset$; at least 5 even and at least 2 multiples of 3 among the
   cofactors); the record, not Astra, is corrected.
6. IX.4: "isolates the triangle identity" read as "equivalently, signed pair residues modulo $q$"; barrier for valuation
   arguments only.
7. IX.5: scope (polylogarithmic binary rank; $n=o(s^{8/7})$) must be stated with the theorem.

## 2. What is new relative to the record and to the literature, and correct

New to the record and correct:
1. **The twin random-parent theorem** (IX.0(1), Theorem T): the first statement of any kind for the two-class problem on a
   natural family, uniform over residues and intervals, with rate $O(\log s)$ where the larger sieve gives order $s$; with the
   [J] corollary that almost every pool parent satisfies the rate form of Conjecture 9 for every $a$ (squarefree case).
2. **The circle random-parent rate at the endpoint $P\asymp s\log s$**, first refereed proof (the statement was in unrefereed
   notes), and **Theorem G**, the bounded-exponent Gaussian extraction by triangles (needed for cap-two norms).
3. **The three-mark proof** with a deterministic heavy block, replacing quartets and the Jensen/matching step.
4. **The same-norm comparison**: the first common norms on which the P2 fakes ($\asymp s/2$ marks) and genuine arcs ($O(1)$)
   coexist; it answers the record's caveat that fakes live only inside the exceptional set (expert note line 940) for bin
   parents.
5. IX.1 (no affine 3-flat at critical length; sharp exponent $2/3$), IX.3 twin (bounded Sidon multiplicity), IX.5
   ($n^{1-2^{1-2d}}$ for prescribed low-complexity profiles, both problems, angle models included).
6. IX.4: the first construction in the record that meets the pair-valuation larger-sieve constant ($\tfrac12$, twin), and the
   unrealisability of Paley cuts in the twin up to length $q^{1-O(1/m)}$ (the circle analogue is open in the record).
7. The correction $(L_8^{pc})$ vacuous, $\mathrm{Bal}^{pc}=\emptyset$, with a proved non-vacuous replacement implied by $(L_8)$.
8. The equivalence of the critical twin with Cilleruelo–Granville's Conjecture 9 at the uniform and both little-oh scales
   (the uniform direction "twin $\Rightarrow$ Conj. 9" is in substance VIII.2's valuation division; the rate equivalences are new).

Not new: the circle pool estimate and the exposure lemma (parallel notes, unrefereed); IX.3 circle (distinct-label lemma);
IX.6 circle criterion (round 8); the counting lemma for sets without affine $t$-spaces (Bonin–Qin type, snippets only); the
twin modulus criterion and density statement are routine. Literature: no full text was read in this round; no result on the
rate form of Conjecture 9 or on random-parent statements of this shape was found (snippets only).

## 3. Rulings

**(a) The random-parent rate is a THEOREM**, in exactly the form of §1.1 (squarefree) and §1.2 (Theorems T and G):
* *Measure on parents:* uniform $s$-subsets of an arbitrary set of $K\ge e^{40}s$ primes in a dyadic interval $[P,2P]$
  (for bounded exponents $J\ge J_0(E)$); also any pool below $X\le s^A$ with $J\ge40A$.
* *Exceptional fraction:* $\exp(-sJ/20)$ (squarefree), $\exp(-sJ/(8E))$ (twin, exponents $\le E$), $\exp(-sJ/(250E^2))$
  (circle, exponents $\le E$); at the smooth scale this is $R^{-o(1)}$ and the exceptional set may contain
  $(\#\text{parents})^{0.95}$ parents.
* *Simultaneity:* one exceptional set serves all intervals of length $\le C\sqrt q$ and all residue choices (twin), and all arcs
  and unit classes (circle); exponents $1\le e_p\le E$ may be chosen after the sample. Checked [J] (§1.1 items 2-4).
* *Bound:* $B_{\max}\le b-1$, $M_{\max}\le4(b-1)$, $b-1=\max(2,\lceil100\log(2P)/J\rceil)$ (squarefree), absolute; $O(\log s)$ at
  $P=As\log s$ with $A\ge e^{40}$ (twin) or $2e^{40}$ (circle), $O_\epsilon(1)$ at $P=s^{1+\epsilon}$.
* *Relation to III.9:* III.9 drew from $[r^D,2r^D]$ at fixed $D>1$ and gave $M_D$ (1 for $D>4$, $\le2$ for $8/3<D\le4$, $O(1)$ for
  $D>2$, 37 at $D=3/2$, finite for $D>1$, unbounded as $D\downarrow1$), circles only. Part IX allows any dyadic pool with
  $K/s\ge e^{40}$, gives the explicit rate $O(1+\log P/J)$, reaches $D=1$ ($P\asymp s\log s$, the height at which the worst case is
  of order $s$), covers the twin, and handles cap-$E$ circles by triangles. It does **not** improve III.9's small constants at
  large $D$, and it does not cover the dense profile of Theorem 1.1. It is an average over a density-zero family; nothing is said
  about any individual $q$ or $N$.

**(b) The same-norm comparison is ESTABLISHED** (a theorem conditional only on refereed results, purely asymptotic), on typical
$s$-subsets of a narrow bin of split primes; the fake half holds on every subset. Meaning for the record's use of fakes as
barriers: until now a barrier read "the inputs $A_{2M}$ cannot prove $M=o(s)$ on the special parents where fakes live, where the
truth is unknown". Now: **on the same norms the truth is $M=O(1)$ while the $A_{2M}$ inputs allow $\asymp s/2$ marks**. So these
inputs carry no information about the genuine count even where it is bounded; they cannot tell such norms apart from any other
subset of the bin. The barrier is strengthened in kind, not extended in range: typical dyadic parents (a dyadic draw lands in
one bin with probability $\approx e^{-0.75s\log s}$, far below the dyadic exceptional fraction) remain open as the Part III judge
posed it, and worst-case parents are untouched (the fake does not see the choice of subset). The separation is by the truth,
not by a test: the fakes fail no test of the list.

**(c) The deterministic results are THEOREMS.** New: IX.1 (twin), IX.3 twin, IX.5 (both problems), the twin modulus criterion,
the density statement (trivial). Not new: IX.3 circle (corollary of the distinct-label lemma, constant improved to $c<2$),
the circle modulus criterion (round 8). What they exclude: signature sets of small doubling; linear (Walsh,
Sylvester–Hadamard) profiles, impossible over $\mathbb Z$ beyond 7 marks at critical length (IX.1); every integer transfer with
linear bits of the characteristic-2 hyperplane mechanisms; prescribed profiles of polylogarithmic binary rank with $n=o(s^{8/7})$
labels (IX.5); parents with one dominant prime power. What they do not exclude: any recorded fake (all Sidon, binary rank of
order $m$), Gold-type full-rank quadratic cores, balanced smooth parents. They are filters, not rates.

**(d) The relaxation IX.4 is a barrier for** every argument that reads the configuration through pair magnitudes, exact gcds,
two-class cuts, valuation (collision-partition) patterns at finitely many places, integral cofactors with a cocycle, signature
structure (Sidon, no affine 3-flat) and approximate real positions; in particular for Gallagher's larger sieve and its
Vandermonde refinement, **sharply** (constant $\tfrac12$). It is **not** a barrier for arguments using signed residues of single
differences modulo the primes of $q$ (already $d_{ij}^2\bmod p$ constant on separated pairs fails), and it is weaker than an angle
model (no exact real or local signed points). It is not a counterexample.

**(e) The printed $(L_8^{pc})$ is VACUOUS** (empty input at 2, and at 3), and $\mathrm{Bal}^{pc}=\emptyset$. Checked [J] against
the note's wording and `round18/theory.md` §5(2). The Question A box must drop "or at least its lattice form $(L_8^{pc})$" and
name the corrected $(L_8^{pc})'$ (pairwise coprime outside 2, 3; $\{2,3\}$-part bounded by $B_0$), state that the printed form had
empty input, and re-pose the judge's Part VIII question "is $(L_8^{pc})$ Rudin-hard?" for $(L_8^{pc})'$ (for the printed form the
answer is trivially no: an unconditional theorem implies Rudin only if Rudin is provable outright). Fact 3 ("every proof of $(L_8)$
must prove $(L_8^{pc})$") holds for the corrected form. The abc criterion's conclusion is vacuous as printed; whether it survives on
the corrected class is unchecked. None of the $(L_8)^{\rm sq}$ results of Part VIII is affected.

**(f) The equivalence with Cilleruelo–Granville's Conjecture 9 is a THEOREM** (§1.5): uniform at one $C$, at every $C$, and at the
little-oh scales $\log q/\log\log q$ and $\omega(q)$. The twin's rate question is therefore a named open problem: do intervals of
length $\sqrt q$ contain $o(\omega(q))$ roots of $x^2\equiv a\pmod q$, uniformly in $a$? No transfer of statements between the twin
and circles is known in either direction (split versus anisotropic torus; the root/representation correspondence does not preserve
short intervals); methods transfer both ways. A rate for the twin is not known to imply Rudin (CG's implication uses the uniform
bound).

**(g) Consequences for the integer problem. Blunt: nothing moves.**
* Growth rate for worst-case parents, $(2/3+\epsilon)\log R/\log\log R$ (Theorem 1.1): **unchanged**. Every rate in Part IX is
  either an average over a pool or for prescribed restricted signature classes that contain no recorded fake.
* Uniform bound: **unchanged**; not known to imply Rudin.
* Questions A and B: **unchanged** in status; the Question A box's fallback target is corrected (it was empty).
* What has changed is diagnostic: at the smooth critical height the truth on typical pool parents is $O(\log s)$ (very probably
  $O(1)$), so the order-$s$ barrier is a worst-case phenomenon, not a height phenomenon; and the input that separates genuine
  from fake, the exact additive identity with integer cofactors, is identified but has been used only through the parent's
  entropy.

## 4. Recommended edits to the expert note

1. **§5 item 12, Part III bullet "Random parents" (lines 934-940; the phrase is on line 940).** Replace "about the narrow-bin parents on which fakes are
   built (they lie inside the exceptional set), or about the growth rate." by "about the narrow-bin parents on which fakes are
   built (they lie inside this exceptional set), or about the growth rate." and append:
   > *Extension (round 31, Part IX of the Question B thread; refereed, `round31/README_IXL.md`; almost-all over a prime pool,
   > not a worst-case statement).* Let $S$ be a uniformly random $s$-subset of any set of $K\ge e^{40}s$ primes in $[P,2P]$,
   > $q=\prod_{p\in S}p$, $J=\log(K/s)$, $b-1=\max(2,\lceil100\log(2P)/J\rceil)$. For $s\ge s_0(C)$, outside a fraction
   > $\exp(-sJ/20)$ of parents, every arc of length $Cq^{1/4}$ of $x^2+y^2=q$ (pool of split primes) carries at most $4(b-1)$
   > points, and every interval of length $C\sqrt q$ contains at most $b-1$ integers lying in one of two prescribed classes
   > modulo each $p\in S$, uniformly in the arc, the interval and the classes. At $P=As\log s$ ($A\ge2e^{40}$) this is
   > $O(\log s)$ at height $\log q\sim s\log s$, the endpoint $D=1$ at which $M_D$ above is unbounded; at $P=s^{1+\epsilon}$ it is
   > $O_\epsilon(1)$. Exponents $1\le e_p\le E$ chosen after the sample are allowed, with constants polynomial in $E$. The
   > exceptional fraction is $R^{-o(1)}$ at the smooth scale and the exceptional set may contain $(\#\text{parents})^{0.95}$
   > parents; the dense profile extremal for Theorem 1.1 ($K=s$) is not covered. On typical $s$-subsets $S$ of a bin of split
   > primes of log-width $l^{-3/4}$ in $[l^{1.99},2l^{1.99}]$ ($s=l-1$), the fakes P2 exist on $N=\prod_{p\in S}p^2$ with
   > $\asymp s/2$ marks while every genuine arc of length $C\sqrt R$ carries $O_C(1)$ points (asymptotic only): the inputs of
   > the fake list cannot certify $M=o(s)$ even where the truth is bounded.
2. **§2, sandwich (line 191).** Replace the right end $(L_8^{pc})$ by $(L_8^{pc})'$, and replace the bullet (line 201)
   "* "pc" restricts it to pairwise coprime $d_T$." by:
   > * "pc" restricts it to pairwise coprime $d_T$. As printed this set is empty: at $p=2$ exactly five of the 119 contact
   >   integers of every point are even, and at $p=3$ at least two are divisible by 3 (round 31, Part IX.2 of the Question B
   >   thread; **proved**, refereed). $(L_8^{pc})'$ restricts instead to $d_T$ pairwise coprime outside $\{2,3\}$, with the
   >   $\{2,3\}$-part of $\prod_{i<j}(x_i-x_j)$ at most an absolute constant.
3. **§2, lattice form (heading line 206, box lines 216-222).** In the heading replace "**open**" by "**as printed, vacuous; corrected form open**".
   Replace the box by the statement $(L_8^{pc})'$ of §1.5 above, and add after it:
   > *Correction (round 31, Part IX.2; refereed, `round31/README_IXL.md`).* The lattice form first printed here (all $e_T$ and
   > all pair cofactors pairwise coprime) has empty input: at most one $e_T$ is even, and then at least five pair cofactors are
   > even; at $p=3$ at least two are divisible by 3. The corrected form follows from $(L_8)$ (blocks are the $\{2,3\}$-free parts,
   > residues the $\{2,3\}$-parts, $H_\Sigma\le\log B_0$), and its coprimality pattern is realised by integers, e.g.
   > $y=(0,5,10,107,267,69,-229)$ with $e_{\{1,2,3\}}=5$.
   In "In words, ..." and the heuristic count, read $(L_8^{pc})'$.
4. **§2, fact 3 (lines 266-269).** Replace by:
   > 3. **Chains carry little information on profiles pairwise coprime outside $\{2,3\}$** (theory Prop 2.1, corrected in round
   >    31, Part IX.2). $\mathrm{Bal}^{pc}$ as first defined is empty. On the corrected class each prime $p\ge5$ divides at most
   >    one block among the $e_T$ and one among the pair cofactors, nested or disjoint, so the chain condition is vacuous at
   >    $p\ge5$, and at 2 and 3 the structure is carried by residues of total size $\le\log B_0$. **Every proof of $(L_8)$ must
   >    prove $(L_8^{pc})'$.** The gloss "nothing but $K\cdot\beta=10$ survives" is not a theorem (item G).
5. **Question A box (line 318).** Replace "Prove $(L_8)$, or at least its lattice form $(L_8^{pc})$, which every proof must cover."
   by "Prove $(L_8)$, or at least the corrected lattice form $(L_8^{pc})'$ (coprimality outside the primes 2 and 3, bounded
   $\{2,3\}$-part), which every proof must cover; the lattice form as first printed had empty input." and replace "the lattice
   form $(L_8^{pc})$ is not known to imply either" by "$(L_8^{pc})'$ is not known to imply either".
6. **§4, abc criterion (line 655).** After "then $(L_8^{pc})$ holds with $Z=\{B=0\}$." add: "As printed the conclusion is vacuous
   ($\mathrm{Bal}^{pc}=\emptyset$, round 31 Part IX.2); whether the criterion gives $(L_8^{pc})'$, where the gcd identity should
   hold up to a factor $B_0^{O(1)}$, is **unchecked**." Treat line 715 ("an abc-conditional $(L_8^{pc})$") the same way.
7. **§5 item 6 (line 823).** Append: "abc on the integer triangle of the two-class twin, $t_{12}A_{12}+t_{23}A_{23}=t_{13}A_{13}$, is
   vacuous for the same reason (squarefree terms). The exact identity itself, with integer cofactors and blocks that are products
   of the parent's primes (twin), resp. $z_j-z_i=2i\eta H_{ij}t_{ij}$ with $t_{ij}\in\mathbb Z$ (circle), is **not capped**: it is the
   input the fakes lack, it gives the random-parent rate (item 12, Part IX of the Question B thread), and for a fixed parent it is a
   squarefree three-term $S$-unit equation with $|S|=s$, whose content is open. Its local consequences at primes not dividing $N$
   (round 9 Prop 5.1(a)), triangle areas and abc remain capped." (Outside the note: qualify `round10/additive.md` §6 and overview
   §5.4 rows 15-16 accordingly.)
8. **§5 item 12, new bullet** (after the bullet on Part VIII of the Question B thread, before "*Withdrawn or corrected in the
   round:*"):
   > * ***Part IX of the Question B thread: the critical two-class problem*** (`round31/README_IXL.md`; eight referee blocks;
   >   nothing false; judge's score 2.5/10). The twin of the arc problem (at most two classes modulo each $p^f\parallel q$, intervals
   >   of length $C\sqrt q$) is equivalent to Cilleruelo–Granville's Conjecture 9 at the uniform scale and at the little-oh scales
   >   $\log q/\log\log q$ and $\omega(q)$ (**proved**). Proved: the random-parent rate for the twin and for circles (Part III bullet
   >   above), including bounded exponents and the smooth scale $P\asymp s\log s$ ($O(\log s)$ marks; the circle estimate was in the
   >   parallel thread's unrefereed notes); the same-norm comparison (P2 fakes with $\asymp s/2$ marks on norms where genuine arcs
   >   carry $O(1)$, typical bin parents, asymptotic only); a squarefree relaxation with all pair magnitudes, gcds, valuation
   >   patterns and Sidon signatures that reaches the pair-valuation larger-sieve constant $\frac12\log q/\log\log q$ and fails the
   >   signed pair residues (a barrier for valuation-based arguments, not for signed residues); no affine 3-flat of twin signatures
   >   at $L<q^{2/3}/4$ (sharp); at most $(1+C^2)(1+2C^2)$ pairs per signature difference in the twin; $n^{1-2^{1-2d}}$
   >   ($n^{7/8}$ at $d=2$) for profiles whose bits are Boolean functions of $d$ affine forms on $n$ labels, in both problems and in
   >   angle models (binary rank $O((\log n)^d)$, so no recorded fake is in the class). The printed lattice form $(L_8^{pc})$ had empty
   >   input (§2). **The growth rate for general $R$, the uniform bound and Questions A and B are unchanged**; worst-case parents, even
   >   squarefree with comparable primes, are open in both problems. Do not write "the uniform bound holds for almost all smooth
   >   parents", "Part IX gives the rate $O(\log\log R)$" (without "on random parents of a pool"), "the exceptional parents are
   >   polynomially rare", "Part IX covers the extremal profile of Theorem 1.1", "the relaxation isolates the triangle identity"
   >   (without "equivalently, the signed pair residues modulo $q$"), "fakes fail the same-norm test", "typical dyadic parents
   >   separate genuine from fake", "$(L_8^{pc})$ is open" (printed form), "the twin rate implies the circle rate" or the converse.
9. **§5 item 12, header and withdrawn list.** Add `README_IXL.md` to the list of judge reports, add 2.5/10 to the scores and eight
   blocks to the count; add to "*Withdrawn or corrected in the round:*" "the claim that the lattice form $(L_8^{pc})$ is open (as
   printed it has empty input; corrected to $(L_8^{pc})'$); the classification of the $\sigma$-cocycle as capped (narrowed to its
   local consequences, triangle areas and abc)".
10. **§6 Pointers (line 1353).** Add `README_IXL.md` and the phrase "random-parent rate at the smooth scale; same-norm separation
    of fakes; $(L_8^{pc})$ correction" to the round-31 row.

Sentences the note must not contain: those listed in edit 8, and "$(L_8^{pc})$ is Rudin-hard" (either form).

## 5. Open sub-questions, ranked by honest probability of progress on the growth RATE for worst-case parents

Probabilities are of moving $o(\log R/\log\log R)$ (or $o(\omega(q))$ in the twin) for worst-case parents in a natural class,
within one more exchange.

1. **Signed pair residues on a fixed squarefree comparable-prime twin** (P ≈ 0.03). The one input that the relaxation lacks and
   that is meaningful for a single parent: at each $p\mid q$ all separated differences satisfy $(y_i-y_j)^2\equiv\delta_p^2$, with a
   sign fixed by one residue vector. IX.3 uses it only for pairs with equal signature difference. *First step:* on the Paley cut
   pattern (unrealisable up to length $q^{1-O(1/m)}$), extract which pair congruences give the contradiction, and test whether the
   same congruences constrain an arbitrary Sidon cut pattern with $B\gg s$. [J] Note on formulation: a relaxation that keeps the
   signed pair residues with one consistent residue vector modulo $q$ and $|d_{ij}|<q/3$ is automatically a genuine configuration,
   so the meaningful relaxation question is the sign-free one: can a table meet $d_{ij}\equiv\pm\delta_p\pmod p$ (sign free per pair)
   with $B\gg s$ and fail the triangle identities?
2. **That sign-free relaxation** (P ≈ 0.02 for the rate; P ≈ 0.3 of a clear answer). A construction would show that pair-level
   residues do not suffice either; an impossibility would reduce the critical twin rate to a question about cut matrices compatible
   with the residue data. *First step:* count, for a Paley cut pattern, the cofactors $t_{ij}$ with $|t_{ij}|\le C\sqrt q/h_{ij}$ in the
   $2^{\#\text{sep}}$ admissible classes modulo $q/h_{ij}$, and search for consistent sign choices at $m=11,19$.
3. **Combination with the parallel thread's projection theorem or with Theorem 1.2's barycentre argument** (P ≈ 0.02). No visible
   mechanism (projected differences of $2i\eta H_{ij}t_{ij}$ do not factor; the barycentre averages lose integrality). *First step:*
   decide whether a large cluster on a fixed parent forces a triangle whose hidden block is constrained by a second, independent
   equation.
4. **Literature: is the rate form of Conjecture 9 known, or the random-parent twin statement?** (P ≈ 0.6 of a clear answer; no rate.)
   *First step:* full texts of Marklof–Welsh (Duke 2023), Shkredov–Shparlinski–Zaharescu (Math. Z. 2024), Cilleruelo–Garaev and
   Chan–Shparlinski on concentration of roots, and the published Cilleruelo–Granville.
5. **Fakes on typical dyadic parents** (P ≈ 0 for the rate; barrier map only; P ≈ 0.15 of a construction). Needs a replacement for
   the pair-margin concentration over windows of width $\gg\sqrt{\log l/l}$ (column-weighted profiles, prime-dependent row
   selection).
6. **Sharpening the random-parent theorem** ($O(1)$ at $P\asymp s\log s$ by a second-moment or longer-cycle exposure, $J_0$ down to
   about 7) (P ≈ 0 for the worst case; P ≈ 0.3 of success). Value low: it remains an average.
7. **Removing the exceptional parents by counting** (P ≈ 0.01): closed for $P\le s^2$ as a method (§1.8).
8. **Bookkeeping** (P ≈ 0.95): the edits of §4.

**Recommendation.** Ask the colleague **one** more question and close the Question B two-class thread if it does not move:
> In the critical two-class problem with $q$ squarefree and all primes in $[P,2P]$ ($P\le s^{A}$), use that at each $p\mid q$ the
> separated differences satisfy $(y_i-y_j)^2\equiv\delta_p^2\pmod p$ with signs given by one residue vector. Either prove $B=o(s)$
> for **every** such $q$ (extending IX.3's sum-and-square device beyond equal signature differences), or exhibit for one fixed dense
> $q$ (for example all primes in $[P,2P]$) configurations with $B\ge cs$, or build an edge table with $|d_{ij}|\le C\sqrt q$, the
> correct gcds and $d_{ij}\equiv\pm\delta_p\pmod p$ on separated pairs (sign free per pair) that reaches $B\gg s$, or prove that no
> such table exists for the Paley cuts. Full proofs or exact computations; say which of these you can do.

If the answer is a fourth random-parent variant or a further restricted class, stop: the thread has then shown what it can show
(the order-$s$ barrier is a worst-case phenomenon, and the separating input is identified but only usable through entropy).

## 6. Scores (0-10, on the goal)

| block | score | reason |
|---|---|---|
| pool (IX.0 squarefree) | 3 | correct, checked [J]; first twin rate and first refereed circle rate at the smooth scale; an average over a density-zero family |
| powers (bounded exponents) | 2 | correct, Gaussian extraction checked [J]; needed for cap-two norms; constants degrade in $E$ |
| samenorm | 2.5 | correct; first common norms separating fakes from genuine arcs; bin parents only, asymptotic |
| additive (IX.1, IX.3, IX.6) | 1.5 | correct; clean characteristic-0 exclusion of Walsh profiles; two items not new; no fake touched |
| fallback (IX.2, CG equivalence) | 1.5 | correct and stronger than claimed; a correction to the record and a named sister problem, not progress |
| relax (IX.4) | 2 | correct; sharp barrier at the larger-sieve constant; interpretation narrowed to valuation-based arguments |
| nonlinear (IX.5) | 1.5 | correct power saving for prescribed low-complexity profiles; outside the class of every recorded fake |
| meaning (referee) | 2 | sound reading; identifies the misfiled cocycle cap; FALSE on the counting route upheld for $P\le s^2$ |
| **Part IX overall** | **2.5** | **nothing false; the first rate $o(\omega)$ at the smooth critical height for a natural measure, and the first same-norm separation of fakes from genuine arcs; the record's $(L_8^{pc})$ fallback was empty; no movement on the worst-case growth rate, the uniform bound or Questions A and B** |

Literature: no full text was read in this round (snippets only, by the fallback, nonlinear and meaning referees). My own checks
are the probability argument of IX.0 (§1.1), the bounded-exponent Gaussian extraction (§1.2), the vacuity of the printed
$(L_8^{pc})$ and of $\mathrm{Bal}^{pc}$ against the note's and round 18's wording (§1.5), the relaxation spot check (§1.6), the
parallel notes' prior statement (§1.1), and the [J] remarks in §§1.3, 1.6, 1.8 and 5.
