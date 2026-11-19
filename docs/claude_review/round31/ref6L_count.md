# Referee report, round 31, Part VI, block "count": VI.2 and the extension test of VI.4

Source: `astra_lemma_VI/answer.md`, Part VI (line 1319 on), §VI.2 (l.1402-1412) and the "extension test"
paragraphs of §VI.4 (l.1479-1487). Scripts: `round31/ref6L_count/` (all mine, exact Gaussian/integer
arithmetic, python3 + sympy for factoring only). Record compared: `growth/cotangent.md` (§1 D1-D4, §3 Theorem B,
block balance (B3), the CG reformulation, heuristic and data), `growth/cotangent_checks/search_cg.py` and
`cg_K8_L120_A1e5.txt`, `round31/README_VL.md`, `round31/ref5L_full.md`.

## 0. Verdict in six lines

1. Every claim in my block is **correct**. No false statement, no gap that affects a conclusion.
2. **(a) volume and (b) heuristic**: correct. They are not new. The heuristic $P^{20-4s}$ at $s=6$ is
   $H^{-1/5}$, which is exactly the record's heuristic $N^{1-6/5}=N^{-1/5}$ (`cotangent.md` l.260).
3. **(c) inverse bound $X^{3/10+o(1)}$**: correct, and $P^6=H^{3/10}$ is right. In substance it is **not new**.
   Astra's word $C$ for the extreme pair is the record's endpoint-chart vector $X+iL$ (`cotangent.md` D1,
   $q_j=(X_j+iL)/(X_j-iL)$). Theorem B gives $L=H^{o(1)}$, the width gives $A\approx LH^{3/10}$, and the
   record's divisibility $X_j-A\mid A^2+L^2$ is a shorter proof of the reconstruction than Astra's
   direction-spacing step. The record never states the count. The explicit statement is new and elementary.
   It does **not** give a closed exceptional set (confirmed). It beats the easy pair-count bound
   $X^{7/10+o(1)}$ by $X^{2/5}$. It is $X^{1/2}$ above the heuristic $X^{-1/5}$. Non-density or finiteness
   would need the count pushed below $X^{0}$.
4. **(c) enumeration (real data)**: I wrote an independent enumerator over $C=A+iL$ with a Gaussian-lcm
   rebuild of each primitive tuple. It reproduces the record's run exactly: 407 ratio classes, minimum
   $K_5=4.9475$, largest $H=683{,}617{,}125{,}425$. Extending $A$ to $10^6$ at $L=6$ finds **no new cluster
   with $K_5\le8$**. This covers $H$ up to about $10^{19}$ at that $L$. **No cluster found has all ten pair
   cuts present.** I prove that this is forced at these heights: a ten-pair profile needs $H\ge1.4\cdot10^{12}$.
   The best width exponent found is $e=0.2246$ at $H=6.8\cdot10^{11}$, against the target $0.3$. The data say
   nothing about (*) itself, except that nothing small exists. With $K_5\le20$ up to $H=4.6\cdot10^{15}$, there is again no ten-pair cluster, and the best exponent is $e=0.2299$.
5. **(d) extension test**: correct. The bound $|q_i-q_j|\ge e^{-d_{ij}/2}$, the identity
   $\sum d_{ij}=6W-2W_{\rm sing}$ and the rigidity $w_{ij}=W/10+o(W)$ are **in the record**:
   (B3) in `cotangent.md` with explicit constants, and Part V's converse (`README_VL.md` §1.1). Astra says so.
   The deletion statements are new, elementary and correct. These are the four-mark profile (singleton $W/10$,
   $2|2$ cuts $W/5$), primitive triples of norm $H^{9/10+o(1)}$ that are three-point sharp, and the
   incompatibility with the Pell four-point profiles (checked numerically).
6. **Significance for the goal: about 0.5/10.** These are bookkeeping and calibration. Nothing here moves
   Question B, (*), the growth rate or $M(C)$.

## 1. Standing cautions, checked explicitly

* *Uniform exceptional sets from fixed-$S$ Subspace*: not used in my block. VI.2 and VI.4 are elementary.
* *Gaussian gcd / denominator costs*: handled correctly. The inverse bound enumerates the singleton and unit
  corrections $A$ (cost $P^{o(1)}$). The deletion statements subtract the common Gaussian gcd prime by prime,
  and I checked this. The four Pell examples were measured after dividing the gcd.
* *Drifting common phase replaced by zero*: no. The phases $\xi$ are carried into $A$, with $A/\bar A=\xi$, and
  the argument of $U$ is taken modulo $\pi$ relative to the known factors.
* *Template-dependent threshold treated as uniform*: the $o(1)$ in the count is uniform given uniform bounds
  on the profile errors. The divisor bound $\tau(n)\le\exp(O(\log n/\log\log n))$ and the direction spacing
  are uniform. Legitimate.
* *Ineffectivity*: none. VI.2 is fully effective. (VI.1's $K,S,b$ are not in my block.)

## 2. (a) Volume and (b) the heuristic

**Rank-4 claim: CORRECT.** $\arg z_k=\sum_e\theta_e-2\sum_{e\ni k}\theta_e$ (+ const), so the anchored map is
$\theta\mapsto-2(\mathbf 1_{k\in e}-\mathbf 1_{1\in e})_{k,e}$. Its rank is 4, because the $K_5$ incidence matrix
has rank 5 (the graph is non-bipartite). Its Smith form is $\mathrm{diag}(2,2,2,2)$ (`volume_check.py`). As a
homomorphism of tori $T^{10}\to T^4$ it is surjective, so the push-forward of Haar measure is Haar. Fixed phases
translate.

**Haar fraction $5(\delta/2\pi)^4$, $\delta<\pi$: CORRECT.** This is the standard $n x^{n-1}$ for $n$ iid marks
in an arc of fraction $x\le1/2$. Monte Carlo with $4\cdot10^5$ samples at $\delta=0.3,1,2,3$ agrees, both for
iid marks and through the actual block map. Chord version: $\Delta=2\sin(\delta/2)$. Correct.

**Heuristic $P^{20-4s}$: CORRECT as a heuristic, NOT NEW.** There are $(3\pi)^{10}P^{20}$ assignments with
$|G_e|\in[P,2P]$, so $H\in[P^{20},4^{10}P^{20}]$ (constants only). Width $P^{-s}$ gives $P^{20-4s}$. At $s=5$
($\Delta=H^{-1/4}$, generic exponent 4) the count is 1. At $s=6$ it is $P^{-4}=H^{-1/5}$ per dyadic window. This
is summable, so the heuristic predicts **finitely many** (*) tuples, not merely non-density. Including $H^{o(1)}$
slack changes nothing. The record has the same count: "$N$ choices of blocks, four angle conditions of
precision $N^{-3/10}$, expected $N^{-1/5}$" (`cotangent.md` l.260-264). Astra's caveat (no discrepancy estimate at
this scale; an arithmetic family would invalidate it) is right. The four-point analogue (count $N^0$ at
$R^{1/3}$, where the Pell families exist) shows that the heuristic is only a heuristic.

## 3. (c) The inverse bound $X^{3/10+o(1)}$

### 3.1 Line-by-line check: CORRECT

* **Word for pair 12.** The blocks separating 1 and 2 are $\{1,j\},\{2,j\}$ for $j=3,4,5$, six in all. From
  $z_k=u\prod_e(\bar G_e$ if $k\in e$, else $G_e)$ we get
  $z_2/z_1=\xi\,\prod_{e\ni1,\not\ni2}G_e/\bar G_e\prod_{e\ni2,\not\ni1}\bar G_e/G_e$. The singleton and unit
  phase $\xi$ is a reduced Gaussian ratio of norm $P^{o(1)}$ (and $i=(1+i)/\overline{(1+i)}$), so $\xi=A/\bar A$
  with $|A|=P^{o(1)}$. Then $C=A\prod G_e\prod\bar G_e$ has $C/\bar C=z_2/z_1$ exactly, and
  $|C|=P^{6+o(1)}$ on the profile.
* **Thin strip.** $C/\bar C-1=2i\,\Im C/\bar C$, so $|\Im C|=\tfrac12|C|\,|z_2/z_1-1|\le\tfrac12P^{6+o(1)}\Delta=P^{o(1)}$.
  Also $\Im C\ne0$. So $C$ lies in $\{|x|\le P^{6+o(1)},\,1\le|y|\le P^{o(1)}\}$, which holds $P^{6+o(1)}$
  points, not the $P^{12}$ of the whole disc. **$P^6=H^{6/20}=H^{3/10}$: right.**
* **Divisor step.** $N(C)\le P^{12+o(1)}$. The number of ordered factorisations into $A$ and six labelled
  factors is at most $\tau_7(N(C))\cdot4^7=H^{o(1)}$. Correct.
* **Other pairs.** The separating sets of 12 and 13 share $\{14,15,23\}$ (3 edges). Those of 12 and 34 share
  $\{13,14,23,24\}$ (4 edges). Correct. After pair 12 the unknown edges are $12,34,35,45$. The 13-word contains
  $12,34,35$, and the 14-word contains $45$. So "every unknown edge occurs" holds, and two further words suffice.
* **Argument of $U$.** $|z_3/z_1-1|\le P^{-6+o(1)}$ gives $\arg C_{13}\equiv0$ modulo $\pi$ to within
  $P^{-6+o(1)}$. Subtracting the known exact arguments places $\arg U$ in a known interval of length
  $P^{-6+o(1)}$ modulo $\pi$. Correct.
* **Primitivity.** At a split $p$, $z_j/z_i$ has $\pi$-valuation $a_j-a_i$. The layers separating $i,j$ are the
  thresholds strictly between $a_i$ and $a_j$, and all of them have the larger one high. So they all point the
  same Gaussian way (monotonicity), and no $\pi\bar\pi=p$ occurs in any sub-product. Inert primes and $1+i$ are
  never layers. Hence every sub-product of the cut blocks of a reduced word is a primitive vector. I checked this
  exactly: **31,744 sub-products from 400 random primitive 5-tuples, 0 non-primitive**
  (`primitivity_check.py`).
* **Spacing.** Two non-parallel primitive $u,u'$ with $|u|,|u'|\le P^{3+o(1)}$ have
  $|\sin\angle|=|\det|/(|u||u'|)\ge P^{-6-o(1)}$. So an interval of length $P^{-6+o(1)}$ modulo $\pi$ contains
  $P^{o(1)}$ directions, each giving $\pm u$. Correct.
* **Total.** $P^{6+o(1)}\cdot P^{o(1)}=X^{3/10+o(1)}$ ordered primitive ratio classes, i.e. primitive tuples up
  to units. The uniformity claim is right. Shared rational primes between layers cause no loss, by the
  monotonicity argument.

Minor point: "ratio classes" and primitive tuples are in bijection up to units and the anchoring. The primitive
tuple is the least scaling, as in V.1. Harmless.

### 3.2 Novelty: the record's endpoint chart in other words

The record's chart (`cotangent.md` D1) writes $q_j=(X_j+iL)/(X_j-iL)$, with $L=L_0$ the lcm of the residues.
So $C_{0j}=X_j+iL$, up to the bounded factor $L/s_{0j}$: Astra's "word with tiny imaginary part" is this
vector. Theorem B with $K=H^{o(1)}$ gives $L_0\le\prod s_{ij}\le K^{10}=H^{o(1)}$. The width
$\Delta=2L/\sqrt{A^2+L^2}$ and the pair lower bound give $A=L\,H^{3/10\pm o(1)}$. The record's integrality
$Q_{0j}=(AX_j+L^2)/(X_j-A)\in\mathbb Z$ gives $X_j-A\mid A^2+L^2$. So:

> *Corrected/simplified proof.* Choose $L\le H^{o(1)}$ and $A\in[LX^{3/10-o(1)},LX^{3/10+o(1)}]$, and then the
> three $X_j$ among the divisors of $A^2+L^2$. That is $X^{3/10+o(1)}$ choices in all. The profile hypothesis is
> not even needed: the width $\le H^{-3/10+o(1)}$ forces it by the extension test.

The record's `search_cg.py` is exactly this enumeration, and the CG equivalence ("bounded $L$, 4-cliques with
$N\le c(A/L)^{10/3}$") is stated there. The record does not state the counting bound as such (grep: no "3/10"
count). So: **correct; the statement is new and elementary; the method is the record's.** Astra's U-spacing step
is a valid but longer substitute for the divisor step.

### 3.3 Comparison with trivial bounds and with what density needs

* Profile alone, without width: $\asymp P^{20}=X^{1+o(1)}$ assignments.
* **Pair count** (my easy bound, any profile). Fix $w=z_2-z_1$ with $|w|\le X^{1/2}\Delta=X^{1/5+o(1)}$. Equal
  norms put $z_1$ on a line with $\ll\sqrt X\,\gcd(w)/|w|$ lattice points in the disc. Summing gives
  $\ll\sqrt X\,X^{1/5}\log X=X^{7/10+o(1)}$, and the other three points then cost $r(H)^3=X^{o(1)}$. **So the
  easy bound is $X^{7/10+o(1)}$, and the one-chord bound gains $X^{2/5}$.**
* Heuristic: $X^{-1/5}$. The gap between proved and heuristic is $X^{1/2}$.
* **No closed set: confirmed.** A count of $X^{c}$ per window, $c\ge0$, says nothing about Zariski closure. Even
  $X^{o(1)}$ per window allows a dense sequence. To exclude (*) by counting one would need a bound $X^{-\eta}$,
  which would give finiteness, i.e. a gain of $X^{3/10+\eta}$. That is far stronger than non-density, and stronger
  than the open CG five-point question (bounded $K$). Astra says this correctly.

### 3.4 The enumeration, as real data on (*)

`enum_chart.py Kmax Amax Ls` enumerates the word $C=A+iL$ of the extreme pair. For each $L$ it sieves the
factorisations of $A^2+L^2$. The candidates are $X=A+d$ with $d\mid A^2+L^2$ and $X\le2L\sqrt T$, where
$T=(K\sqrt{A^2+L^2}/2L)^{10/3}$. This bound is complete: $n_{0X}\ge(X^2+L^2)/(2L^2)$ and $n_e\le H\le T$. Cliques
come from pairwise integrality, and pruning uses the lcm of the edge norms. **Every found tuple is rebuilt
independently** as $Z=\mathrm{lcm}_G(\bar\beta_j)$, $z_j=Z(X_j+iL)/(X_j-iL)$. The script asserts
$N(Z)=\mathrm{lcm}_e n_e$ (record D2, verified on all tuples), equal norms, Gaussian gcd 1, and width
$=2L/\sqrt{A^2+L^2}$. `cutprofile.py` factors $H$ and computes all fifteen cut weights.

Runs (A ≥ L, i.e. arc ≤ π/2, unless noted):

| run | range | unique ratio classes | min $K_5$ | largest $H$ |
|---|---|---|---|---|
| run2 | $K\le8$, **all** $L\le60$, $A\le2\cdot10^4$ | 285 (180 with $A\ge L$) | 4.9475 | $1.6\cdot10^{10}$ |
| run1 | $K\le8$, $L\in6\mathbb Z$, $L\le120$, $A\le10^5$ | 407 (244) | 4.9475 | 683,617,125,425 |
| run3 | $K\le8$, $L\in\{6,12\}$, $A\le10^6$ | 48 and 45, the same as for $A\le2\cdot10^4$ | | |
| run4 | $K\le20$, $L\in\{6,12,18,24,30\}$, $A\le10^5$ | 2218 (1928) | 4.9475 | $4.55\cdot10^{15}$ |

* run1 reproduces the record's run (`cg_K8_L120_A1e5.txt`) exactly: the same minimum 4.9475 at $N=65$ and the
  same largest $N$. The record's count of 1476 includes repeats over multiples of $L$.
* run2 independently confirms $6\mid L_0$. Over all $L\le60$, only multiples of 6 produce clusters, and every
  found $L_0\in6\mathbb Z$.
* run3: **no cluster with $K_5\le8$ and $L\in\{6,12\}$ has $2\cdot10^4<A\le10^6$**. This range contains every such
  cluster with $H$ up to about $(8\cdot10^6/12)^{10/3}\approx10^{19}$.

Best width exponent $e=-\log\Delta/\log H$ per decade (run1; target $0.3$; $0.3-e=\log K_5/\log H$):

| $H\sim$ | $10^3$ | $10^5$ | $10^7$ | $10^9$ | $10^{10}$ | $10^{11}$ |
|---|---|---|---|---|---|---|
| best $e$ | 0.097 | 0.170 | 0.195 | 0.212 | 0.214 | 0.225 |
| its $K_5$ | 5.76 | 5.99 | 5.42 | 6.62 | 7.56 | 7.80 |
| nonzero pair cuts | 4 | 5 | 5 | 6 | 7 | 7 |

The exponent creeps toward 0.3 only because $K_5$ stays bounded (about 5-8) while $\log H$ grows. That is the
trivial behaviour of a constant. No found cluster has all ten pair cuts present (0 of 244 in run1). The best
"near-profile" cluster at $10^7$ has 5 pair cuts at 0.16-0.21 and singleton weight 0.10.

**Lemma (mine, elementary).** A primitive 5-tuple has at most two distinct pair cuts per prime. The high sets of
one prime form a nested chain, so only sizes 2 and 3 give pair cuts. A ten-pair profile therefore needs at least
5 split primes, each with $e_p\ge2$. So
$H\ge(5\cdot13\cdot17\cdot29\cdot37)^2\approx1.4\cdot10^{12}$, and approximate balance needs far more. **Hence all
data in the record and in run1/run2 lie below the first height at which a (*)-type tuple can exist.** The record's
"data on CG" and mine bear on (*) only through run3 (no cluster at all with $K_5\le8$ up to about $10^{19}$ at
$L=6$) and run4.

### 3.5 run4 (larger $K$)

run4 uses $K_5\le20$, $L\in\{6,\dots,30\}$ and $A\le10^5$, and reaches $H$ up to $4.6\cdot10^{15}$, beyond the
ten-pair threshold. It finds 1928 classes with $A\ge L$ and **0 with all ten pair cuts present**. The largest
number of distinct nonzero pair cuts in the per-decade best clusters is 8 (at $H=8.5\cdot10^{12}$, $K_5=10.0$).
Singleton weight stays at 0.1-0.28 of $W$, far from $o(W)$. Best exponents per decade:

| $H\sim$ | $10^{11}$ | $10^{12}$ | $10^{13}$ | $10^{14}$ | $10^{15}$ |
|---|---|---|---|---|---|
| best $e$ | 0.205 | 0.2225 | 0.216 | 0.218 | **0.2299** |
| min $K_5$ | 11.2 | 10.0 | 13.8 | 15.2 | 11.9 |

The best exponent found anywhere is $e=0.2299$ ($H=2{,}086{,}319{,}768{,}148{,}425$, $L=6$,
$X=(39858,40708,68808,630438)$, $K_5=11.87$, 6 pair cuts, singleton weight $0.20W$). The minimal $K_5$ per decade
does not decrease with $H$. The counts per decade fall roughly like the heuristic $K^4H^{-1/5}$, but the data are
not a test of (*). Caveat: completeness holds only for the listed $L$. Theorem B allows $L_0$ up to
$K^{10}/128$ ($8\cdot10^6$ at $K=8$), so none of these runs (nor the record's) is exhaustive in $L$.
To push the enumeration to the (*) regime ($H\gtrsim10^{20}$ with near-balanced weights) one would need
$A\approx LH^{3/10}\gtrsim10^{7}$ and $K$ of order $H^{o(1)}$. The cost is about $A\cdot\tau(A^2+L^2)^2$ per
$L$, which is feasible in C but not in this session's Python. I recommend this as the one concrete computational
follow-up, with the expectation (heuristic) of finding nothing.

## 4. (d) The extension test

* **$|q_i-q_j|\ge e^{-d_{ij}/2}$: CORRECT.** $z_j/z_i=\xi C/\bar C$ with $C=z_j/\gcd(z_i,z_j)$ and
  $|C|^2=e^{d_{ij}}$, and $\xi C-\bar C$ is a nonzero Gaussian integer. I checked this exactly on 273 random
  primitive 5-tuples, including $|C|^2=H/N(\gcd)=e^{d_{ij}}$: 0 violations (`extension_check.py`).
* **$\sum_{i<j}d_{ij}=6W-2W_{\rm sing}$: CORRECT.** A $2|3$ cut separates 6 pairs and a $1|4$ cut separates 4, and
  every layer of a primitive tuple is a nontrivial cut. Checked exactly on all 273 tuples.
* **Rigidity.** $d_{ij}\ge(3/5-o(1))W$ and the sum give $W_{\rm sing}=o(W)$ and $d_{ij}=(3/5+o(1))W$. With
  $d_{ij}=w_i+w_j-2w_{ij}+o(W)$ and $\sum_{j\ne i}d_{ij}=w_i+2W_{\rm pair}$, we get $w_i=2W/5$ and
  $w_{ij}=W/10$. CORRECT. **Not new**: this is `cotangent.md` (B3), $|\log b_T-W/10|\le60\log K$ via
  $(I+\mathrm{Pet})^{-1}$, and Part V's converse (`README_VL.md` §1.1, `ref5L_full.md` §1.7). Astra cites the
  record.
* **Deleting one mark: CORRECT, new (elementary).** Of the cuts, $\{5\}$ becomes constant, $\{a,5\}$ becomes the
  singleton $\{a\}$, and the pairs $\{a,b\},\{c,d\}$ with $a,b,c,d\ne5$ merge into one $2|2$ cut. So the four-mark
  norm is $H^{1-o(1)}$, the four singletons have $W/10+o(W)$ and the three $2|2$ cuts have $W/5+o(W)$ (total $W$).
  The primitive norm drops exactly by the constant layers, which is the gcd, as Astra says.
* **Deleting two marks: CORRECT, new (elementary).** $\{4,5\},\{4\},\{5\}$ become constant, so the norm is
  $H'=H^{9/10+o(1)}$. Each triple singleton gets $3W/10$. Normalised widths are unchanged by removing a common
  factor, so $\Delta'=H^{-3/10}=H'^{-1/3}$. This is the sharp three-point exponent (area $\ge1/2$). So all ten
  triples are simultaneously three-point sharp. Correct.
* **Pell incompatibility: CORRECT (checked).** A four-mark sub-tuple of (*) has two-sided width
  $H'^{-3/10\pm o(1)}$, so an $R^{1/3}$-sharp four-tuple cannot be one. The mixed-rate $(6,7,7)$ Pell example has
  only three $2|2$ cuts. With my cut code (`pell_profile.py`) at $n=1,2$ the weights are
  $(0.28\text{-}0.30,0.34\text{-}0.35,0.34\text{-}0.35)$ with singleton $\le0.04$ (bounded), and the width exponent
  is $0.304, 0.299$. The required profile is $2|2$ cuts at 0.2 and singletons at 0.1. Wrong profile, as Astra says.
* **Comparison with the record.** The record has the rigidity with constants (B3) and the sub-tuple projections
  (expert note §on $\mathrm{pr}_S$). The record has no statement of the four-mark and triple profiles of a (*)
  tuple, nor of the simultaneous three-point sharpness. These are new, but each is a few lines. They give a
  necessary condition, not an obstruction: three-point-sharp triples are abundant.

## 5. Significance (blunt)

* VI.2 contains no new idea. The volume count is the record's heuristic, and the "proved inverse bound" is the
  record's endpoint chart plus a divisor bound, now stated as a count. The count is $X^{1/2}$ away from anything
  that would decide (*). It cannot yield a closed set.
* The extension test is the record's rigidity plus elementary sub-tuple bookkeeping.
* The data are honest negatives. No five-point cluster with $K_5\le8$ exists beyond $H=6.8\cdot10^{11}$ in the
  ranges searched (to about $10^{19}$ at $L=6$), and nothing searchable at these heights can have the ten-pair
  profile. This is consistent with the heuristic $K^4H^{-1/5}$ but proves nothing.
* Score: **0.5/10** on the goal. Labels: all PROVED claims in my block confirmed. The HEURISTIC labels are
  appropriate.

## 6. Files

`ref6L_count/volume_check.py`, `glib.py`, `extension_check.py`, `primitivity_check.py`, `enum_chart.py`,
`cutprofile.py`, `analyse.py`, `pell_profile.py`; outputs `run1_K8_A1e5_L6-120.txt`, `run2_K8_A2e4_Lall60.txt`,
`run3_K8_A1e6_L6-36.txt` (partial: $L=6,12$ completed; killed for time), `run4_K20_A1e5_L6-30.txt`.
