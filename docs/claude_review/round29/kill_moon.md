# Round 29: adversarial review of lane moon (kill_moon)

*Reviewer's labels: proved, exact, evidence, heuristic. Scripts are in `kill_moon_checks/`. Working notes are in `kill_moon_notes.md`.*

## 0. Verdict

**The obstruction verdict survives.** moon claims no growth improvement, and nothing in it asserts one. Every
screen in §3 reaches the right *conclusion*. The function-field development (§4) is correct as stated, with
small exceptions. However, the review found the following:

1. **A false data claim (idea 17).** moon says the reduced chords of the five record clusters carry "exactly the
   pigeonhole-forced collisions at ℓ=3 (4,4,4,3,1 pairs)". The actual counts are **4,4,4,4,2**, and the forced
   minima are 4,4,4,3,1. So the 7-point and the 5-point clusters each have **one collision above the forced
   minimum**. moon's own output file shows `3: (4, 8)` and `3: (2, 6)` for these two clusters, so the text
   misreports it. The qualitative point (no hidden identity) still stands as weak evidence.
2. **An invalid proof (idea 1, Erdős–Turán).** "$L(G)\le 2^n$, so the bound is $\le16\sqrt{\log2}>1$, hence vacuous"
   is a non sequitur. An *upper* bound on the ET bound cannot show that the bound is useless. The conclusion is
   rescued below by a lower bound, but only for the record clusters, not "for every N". The F1 kill is unaffected.
3. **A misstated scope (F_3 data).** The scans are not over "all squarefree split N". They cover only those whose
   Gaussian-type prime factors $\pi$ have degree ≤3 (`ff_check`) or ≤4 (`ff_scan`). A complete scan for
   $D\le6$ (106,149 N, against about 18k in moon's scope) gives the same maxima, so the data statement holds for
   $D\le6$ after correction. For $D=7,8$ it is unverified outside the restricted scope.
4. **A wrong heuristic in the lead (open question).** moon says "the generic count predicts 5". That count treats
   the $M-1$ differences of a fibre as independent. Allowing dependent differences ("Sperner blocks", §2.4) gives
   **6-point clusters over $\mathbb R(t)$ at $D=8$**. This is proved: an explicit algebraic construction plus
   Lemma 4.5. There is also a 6-point cluster at $D=12$ (numerical evidence, residual $10^{-41}$). The record's
   integer-box template search found only 4 at $D=8$. The open question must be restated (§5).
5. **Novelty is overstated.** Lemma 4.5 is `growth/families.md` Lemma P (moon says so). The record already has a
   function-field analogue (round 7, `L7_prev_attempt1.md` Corollary 2: "$(L_7)$ is false over $\mathbb C(t)$").
   Observation 4.4 is the Plotkin-extremal case of a one-line Plotkin bound $M_{\mathbb R}(0)\le 2D$ that moon
   does not state.
6. **A sharper idea-14 kill (new).** Every lattice point $z$ of $x^2+y^2=N$ divides $N$ in $\mathbb Z[i]$
   ($N=z\bar z$). So every character modulo $N$ (or modulo any multiple of a prime of $N$ that $z$ contains)
   vanishes identically on the circle. The Li–Liu device of "fixing the modulus at the integer $N$ itself" has
   no non-trivial analogue here.

## 1. Screens of §3 re-derived

* **Obs 1.1 (product set): correct**, by unique factorisation. The 2-part and the inert part multiply all points by
  one fixed Gaussian integer.
* **Idea 1 (ET): conclusion correct on the data, but the stated proof is invalid; the F1 kill stands.**
  Let $G=R^{-n}F_N(RX)$, with $|a_0|=|a_n|=1$. By Parseval, $\sup_{|X|=1}|G|\ge(|a_0|^2+|a_n|^2)^{1/2}=\sqrt2$. So
  every ET-type bound $C\sqrt{\log(\|G\|/\sqrt{|a_0a_n|})/n}$ is at least $C\sqrt{\log\sqrt2/n}$. A cluster gives an
  excess of at most $M/n$, so ET can detect it only if $n<M^2/(C^2\log\sqrt2)$. That is $n<M^2/88.7$ for
  $C=16$, and $n<M^2/2.27$ for Ganelius's $C=\sqrt{2\pi/k}\approx2.56$ (literature found, §4).
  * Exact, `nf_clusters.py`: $n=r(N)=192,1152,2048,288,768$ against $M^2/2.27\le28.2$. So ET is useless on all
    five clusters.
  * In general, $r(N)\ge cM^2$ is not proved by moon or here.
  * ET holds for arbitrary polynomials, so fakes satisfy it, and the F1 kill is valid regardless.
* **Idea 2: correct.** $|\mathrm{disc}F_N|=\prod|z-z'|^2\ge1$ already follows from $|z-z'|\ge|g|\ge1$ (pair data).
* **Idea 3: formula off by a factor 4.** $s_m=4R^{4m}(-1)^{me_2}\prod_p\frac{\sin(4m(e_p+1)\phi_p)}{\sin4m\phi_p}$. This
  was verified numerically for $m=1,2,3$ on all five $N$ (`nf_clusters.py`). The kill is unaffected.
* **Idea 4:** it is a screen, not a theorem, and moon says so. Agreed.
* **Idea 5: correct.** Re-derivation: with $u=z_x-z_w$ and $v=z_y-z_w$, $|z_w+u+v|^2=N+2\mathrm{Re}(u\bar v)$, which
  equals $N+n_{xw}+n_{yw}-n_{xy}$. Andrews' theorem ($\mathrm{area}\gg n^3$) is confirmed in the literature.
  The segment area is $\approx C^3\sqrt R/12$, which gives $M\ll CR^{1/6}$. That is weaker than known results;
  it is harmless.
* **Ideas 6–9, 11–13, 15, 18: no objection** beyond the remark that they are verbal screens. Idea 9 is right,
  and it is even sharper in the FF setting (§2.3).
* **Idea 14: moon's account of the paper is accurate at abstract level, but the screen is imprecise.**
  * What the PDF says (read in place, pp. 1–3): "Jiamin Li & Jianya Liu, *A conditional proof of the Goldbach
    conjecture*". It assumes **GRH(N)**, that is, GRH only for primitive Dirichlet characters *whose conductor
    divides N*, together with one of two spectral hypotheses on the zeros of $\mathfrak X_N(q)$, a function on
    the quaternions built from the eigenvalues of a compact self-adjoint Hilbert–Schmidt operator $T_N$. The
    abstract says that odd characters mod $N$ isolate the correlation and that an exact quadratic-form identity
    expresses the weighted Goldbach error.
  * moon's phrase "GRH concerns all primes, not the ω primes of N" misses that the paper's hypothesis is tied to
    $N$ through the modulus.
  * The correct kill is item 6 of §0. Characters mod $N$ vanish on every point because $z\mid N$, which is proved
    and trivial. Characters mod divisors $\pi_p^k$ see only the profile and the residues mod the primes of $N$
    (idea 10).
  * The arc count is a finite exact sum with no prime sum in it, so GRH has nothing to act on.
  * Web check: a search did not find the conditional paper online. It found the same authors' "Theorem (1+1.9) on
    the Goldbach Conjecture" (arXiv 2606.05224), which the PDF cites as [9]. So the uploaded paper's status is
    unverified.
* **Idea 10:** not capped, no mechanism. Agreed. It remains the only screened structure that is not angle data.

## 2. Function-field section re-derived

**2.1 Lemma 4.1: correct.** In (iv), $h$ must be monic so that monic points stay monic. The claim "$M_k(c)$ is
constant on $c<0$" holds: monotonicity in $c$, plus composition with $\kappa\ge-1/c$, gives $M(c)\le M(\kappa c)\le
M(-1)\le M(c)$. Calling (iv) "the FF form of Lemma K1" is an analogy (**heuristic**), not a theorem. The exact FF
analogue of the NF scaling $z\mapsto gz$ is multiplication by $g\in B$. That sends $c\mapsto c+\deg g/2$, whereas
composition has no NF counterpart.

**2.2 Props 4.2 and 4.3: correct.** I re-derived both.
* Prop 4.2: with $\sum_\pi e_\pi\deg\pi=D$, $\binom M2(s-D)+M^2D/4\le\binom M2c+MD/4$.
* Prop 4.3, residue counts:
  * odd $d$: inert, residues in a norm-one coset of size $q^d+1$;
  * even $d$: $\ell=\lambda\bar\lambda$ with $\deg\lambda=d/2$ and $B/\lambda=\mathbb F_{q^d}$; a collision mod
    $\lambda$ implies a collision mod $\bar\lambda$.
* Prop 4.3, constants: the bound $dI_q(d)\ge q^d-2q^{d/2}$ holds, since
  $\sum_{e\mid d,e<d}eI_q(e)\le\frac{q}{q-1}q^{d/2}$, and $K_q$ is generous.
* (H) is a real restriction, and it is labelled as one.
* Exact check (`ff3_scan.py`, independent code): 23,238 clusters at $c=0$, 0 violations of Prop 4.2. This
  reproduces moon's count, which shows moon used the same restricted scope.
* Prop 4.3 is not tested by the data, because $M\ge q^2=9$ never occurs.

**2.3 Obs 4.4: correct, and stronger than stated.** For squarefree $N$, a pair needs $\deg g_{xy}\le s$, so the
Hamming distance is $\ge\lceil D/2\rceil$. Plotkin then gives **$M_{\mathbb R}(0)\le2D$ for even $D$ and $\le D+1$
for odd $D$** (proved, for squarefree $N$). So "local" FF inputs give $M=O(\log R)$, and Hadamard profiles show this
is sharp for local data. moon's idea 9 (Plotkin) is thus exactly the FF local bound.

**2.4 Lemma 4.5 and Prop 4.6: correct.** Lemma 4.5 is `families.md` Lemma P. For Prop 4.6, the DBTK level sizes
6, 4, 3, 2, 1 were checked by hand. DBTK itself is confirmed in the literature: Nieuw Arch. Wisk. 23 (1952).
**The "generic count predicts 5" heuristic is wrong. A new construction (Sperner blocks):**
* Split the $D=4b$ roots into four disjoint blocks $B_j$ with equal vectors
  $T=(\sum_{k\in B_j}\mathrm{Im}\,a_k^r)_{r\le m}$, where $m=2b-1$.
* The six sign patterns that are $+$ on two blocks and $-$ on two blocks satisfy $\sum\varepsilon_kw_k=0$.
  Pairwise they differ on at least two blocks (distance $\ge D/2$). By Lemma 4.5 they form a 6-point ball cluster
  at $c=0$.
* Count: each block has $2b$ real unknowns and $m=2b-1$ equations, so the condition is not over-determined.
* $D=8$ (`rt_sperner_blocks.py`, **proved**): the blocks are $\{i,1+2i\}$ with $T=(3,4,-3)$. Three further pairs
  have $x_1\in\{-1,2,-2\}$ and $y_1$ a real root in $(0,3)$ of an explicit polynomial, isolated exactly by sympy.
  The 8 roots are distinct and lie in $H$. A 60-digit check confirms that coefficients 5–7 agree and coefficient 4
  differs. So $M_{\mathbb R}(0)\ge6$ at $D=8$; the record's integer box gave 4.
* $D=12$ (`rt_blocks12.py`, **evidence**): four triples found by mpmath. Coefficients 7–11 agree to
  $5\cdot10^{-41}$.
* The correct heuristic is rank-based. The differences of a fibre span a space of rank $\rho$ with
  $\rho m\le2D-2$, so $\rho\le4$ for $D\ge7$. Antichain and distance constraints inside a rank-4 structure
  suggest a bounded answer of about 6. Larger $\rho$-block designs (5+ blocks, Hamming-code patterns) fail the
  count or the distance. This is **heuristic**; whether 7 is possible at large $D$ is open.

**2.5 F_3 data, re-run.** `ff3_complete.py` covers **all** squarefree split $N$ with $D\le6$: 106,149 N, with prime
pairs of every degree, found by python-flint factorisation over $\mathbb F_9$. The maxima at $c=-1,0,1$ are
$(1,2,2),(2,4,4),(2,3,8),(2,4,8),(2,3,6),(3,4,8)$ for $D=1..6$. This agrees with moon. moon's $D=7,8$ values
hold only within its restricted scope.

## 3. Recomputation of numerical claims (summary)

| claim | status |
|---|---|
| ℓ=3 collisions = forced min (4,4,4,3,1) | **false**: actual 4,4,4,4,2 |
| every prime ℓ>4k divides ≤1 chord | confirmed, exact |
| odd ℓ≤4k extras {5:1},{7:2,31:1} | confirmed |
| C values 5.8977, 5.6739, 5.7542, 5.2699, 1.9962 | confirmed |
| ET vacuous for every N | proof invalid; true on the 5 clusters |
| $s_m$ formula | off by factor 4 |
| 23,238 F_3 clusters, 0 violations | confirmed, but scope is π-degree ≤3, not "all N" |
| F_3 maxima 3 (c=−1), 4 (c=0) | confirmed for all N, D≤6 |
| generic count 5 over R(t) | refuted at D=8 (6 points, proved) |

## 4. Literature checked by web search

* Erdős–Turán constant 16; Ganelius $\sqrt{2\pi/k}$; Amoroso–Mignotte lower bound $\sqrt2$. These come from
  secondary sources (arXiv 1307.5734 and others).
* Andrews 1963, Trans. AMS 106 (area $\gg n^3$).
* de Bruijn–Tengbergen–Kruyswijk 1952.
* No function-field Cilleruelo–Granville analogue was found. Related work found: Sedunova (Bombieri–Pila over
  $\mathbb F_q(t)$) and Yiasemides (lattice points in annuli over $\mathbb F_q[t]$). So Prop 4.3 appears new, but
  it is a direct transfer of Theorem 1.1's mechanism.
* The Li–Liu conditional paper was not found online (see §1).

## 5. Novelty and corrected statement

**Novelty: repackaged known barrier, with a small new calibration.**
* Repackaged: the screens restate expert-note §5. Lemma 4.5 is Lemma P.
* New, but small and direct:
  * the F_q analogue Prop 4.3 under (H);
  * the R(t) Sperner case Prop 4.6;
  * from this review: the Plotkin bound $M_{\mathbb R}(0)\le2D$, the 6-point clusters at $D=8$ (and $D=12$), and
    the character-vanishing kill of the Goldbach device.

**Corrected statement.**
* No growth improvement; $M\le(2/3+\varepsilon)\log R/\log\log R$ stands.
* Over $\mathbb F_q(t)$ under (H): $M\le(\frac12+o(1))\log_qR/\log_q\log_qR$.
* Over $\mathbb R(t)$ at $c=0$:
  * $M\le2D$ (Plotkin, sharp for local data, by Obs 4.4);
  * $M=6$ at $D=4$;
  * $M\ge6$ at $D=8$ (proved) and at $D=12$ (evidence).
* Idea 17 is weak evidence only: the 8-point clusters sit at the forced ℓ=3 minimum, and the 7- and 5-point
  clusters sit one above it.
* The open question becomes: *is $\sup_D M_{\mathbb R}(0)<\infty$, and is it 6?* "Generic count 5" must be
  dropped. A rank heuristic suggests about 6. Any proof must allow the 6-point Sperner-block clusters, which exist
  at $D=4,8$ and, numerically, at $D=12$.
