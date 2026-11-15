# Round 30 -- adversarial check of lane `search`

## 0. Verdict

**Status of lane `search`: survives in part, with corrections. One proved claim is false, and the lane's "lead" mostly collapses.**

1. **Confirmed (independent code).**
   * Lemma 1.1, Prop. 2.1 and Thm 3.1 hold. Thm 3.1 has one misstated degeneracy line, but the conclusion holds.
   * The horizontal-line bound $M\le2$ holds.
   * Balanced clique numbers: 4, 7, 7, 14, 6 at $c=0$ and 14 at $D=7$, $c=1$, so $M(8,0)\le14$ and $M(9,0)\le6$.
   * Circle-locus numbers match.
   * All rank-$\le3$ and rank-$\le4$ code maxima match.
   * Every Krawczyk certificate was reproduced with my own certifier and gauge ($D=9$; $D=10$; five $D=12$ half-cube solutions; the three $c=1$ rank-5 clusters). In every case the maximum cluster over **all** $2^D$ points of $N$ equals the cluster through $Q$.
2. **False claim.** "Vertical lines give rank $\le3$ at $c=0$" (lane §3.2, repeated in the summary). It is a parameter count, not a proof, and it is false. I give **Krawczyk-certified rank-4, 5-point clusters on the vertical line at $D=10$, $c=0$** (§4.1).
3. **Overstatements.**
   * "8-point rank-$\le4$ needs $4\mid D$ for $D\le16$ by computation": the lane had no data at $D=11,13,15$. I supply it (max 5), so the claim now stands for $D\le16$.
   * "Along $D\equiv2\ (4)$, (R$_4$) gives $M\le6$": shown only for $D=10,14$.
   * "Circle locus $M\le6$" does not exclude rank 5.
4. **Lead instance 1 is not real-algebraic.** The "141 six-point rank-4 codes at $D=10$, expected dimension 2, never realised" are **one** code up to symmetry, and its span contains a weight-2 integer vector. The lane's own Lemma 1.1 kills it. The same holds at $D=14$ (2 classes, both killed). New **Proposition K1 (proved)**: at $D=10$ and $D=14$, $c=0$, every cluster with $\ge6$ points has rank $\ge5$. Hence $M(10,0)=5$ unless rank $\ge5$ occurs.
5. **Lead instance 2 is weakened.** Of the "16 of 46 $D=8$ five-point types never realised", my solver realises 35 of 46 with the prescribed orientation (lane: 30), leaving 11. Two of them, types 9 and 12, which the lane's log lists as unrealised, I certified.
6. **Stop rule.**
   * No rank-5 cluster at $c=0$ for $D\le12$ was found. That is correct, and it is evidence only.
   * Proved rank bounds on positive-dimensional loci: the cube locus (rank exactly 4) and the horizontal line ($M\le2$, hence rank $\le1$). The vertical-line bound is withdrawn.
   * So the second condition is formally not met, but only through weak bounds. The real-algebraic "lead" offered as the reason to continue is, in its cleanest instance, linear algebra.
   * **I recommend closing the function-field line for internal work.**


## 1. Proofs re-derived

All re-derived by hand, with all quantifiers.

* **(S1), (S3), (S5): correct.** The (S1) sign is immaterial. (S5) is triangular with leading term $mx^{m-1}$.
* **Lemma 1.1 (balance on $V\cap\mathbb Z^D$): correct.**
  * The equations are real-linear in $a$.
  * $p_m(H)-\overline{p_m(H)}=2i\sum_k a_k\,\mathrm{Im}\rho_k^m$, so $e_j(H)$ is real for $j\le L$.
  * $\mathrm{Im}H\not\equiv0$, because $H$ uses at most one of $\rho_k,\bar\rho_k$ for each $k$.
  * $H/\bar H=\prod((t-\rho_k)/(t-\bar\rho_k))^{a_k}=e^{i\Theta}$, and $\Theta$ is continuous from $-2\pi(a^+-a^-)$ to $0$.
  * It meets the $|a^+-a^-|-1$ interior multiples of $2\pi$ at distinct real $t$, each a zero of $\mathrm{Im}H$.

  This gives $\min(a^+,a^-)\ge L/2$ and $\|a\|_1\ge L+1$. It is an elementary extension of (S2), but it is the lane's most useful tool (see §2.3).
* **Prop. 2.1 (half-cube needs $4\mid D$): correct.**
  * The six nontrivial characters other than $\chi$ vanish.
  * Fourier inversion gives $m_u=(D+\lambda\chi(u))/8$.
  * $D\pm\lambda\equiv0\ (8)$ gives $2D\equiv0\ (8)$.

  It is proved only for the half-cube combinatorial type.
* **Thm 3.1 (cube locus): correct, with one wrong line.** The moment equations on the locus ($\mathrm{Im}\,r^2=2\varepsilon$, $\mathrm{Im}\,r^3=3/y-y^3$) and the case analysis by $\mathrm{supp}\,s$ check out. Odd support is excluded by parity of $\sum t\varepsilon$.
  * The third degeneracy case is misstated. $y_j=-y_i$ and $\varepsilon_j=\varepsilon_i$ give $x_j=-x_i$, so $r_j=-r_i$, not $-r_j=\bar r_i$. It is still a degeneracy, since $-r_i$ is a root of $Q$.
  * Scope: the theorem bounds the cluster **through $Q$**. I checked the stronger statement numerically (`cubeloc.py`, 60 random nondegenerate locus points, all 256 points of $N$, exact max clique). The largest cluster on $N$ is 8 at every sample. This is evidence, not proof, for clusters not containing $Q$.
* **§3.2, horizontal line: correct.**
  * The rows are $x^0,\dots,x^{L-1}$, so $a^+=a^-$.
  * $a^+\le L-1$ would make two disjoint multisets of size $\le L-1$ share their first $L-1$ power sums, which is impossible.
  * So distances are $\ge2L$, which is even. Plotkin, for even $d$ with $2d>n$, gives $A(n,d)\le2\lfloor d/(2d-n)\rfloor$ (checked against the literature, §5), hence $M\le2$ for $D\ge7$.
* **§3.2, circles: the lemma is correct.** $\sin m\theta=\sin\theta\,U_{m-1}(\cos\theta)$, and the standard sign-change argument applies. But "$M\le6$" bounds the size, not the rank: a 6-point cluster may have rank 5. So the circle result **does not exclude rank 5** for $D\le10$.
* **§3.2, vertical lines "$r\le3$ at $c=0$": FALSE as stated and not proved.**
  * "$r\lceil L/2\rceil\le D-1$" is a parameter count, the same heuristic kill_ff called "no excess". It is not a theorem.
  * At $D=10$ the count itself allows $r=4$ ($4\cdot2=8\le9$).
  * Counterexample, Krawczyk-certified (§4): rank-4 clusters on the vertical line at $D=10$, $c=0$.


## 2. Independent recomputation of combinatorial bounds

All code is mine and is in `kill_search_checks/`. None of the lane's scripts was imported.

**2.1 Balanced clique numbers** (`balclique.py`, networkx exact clique, split by an $S_D$-representative vertex).

| $D$, $c$ | 5,0 | 6,0 | 7,0 | 8,0 | 9,0 | 7,1 |
|---|---|---|---|---|---|---|
| clique number | 4 | 7 | 7 | 14 | 6 | 14 |

These match the lane. $M(8,0)\le14$ and $M(9,0)\le6$ are **confirmed** (proved, given Lemma 1.1).

**Circle locus** (`circ.py`): the values 6, 4, 6, 4, 6 for $D=6,\dots,10$ match.

**2.2 Rank-$\le r$ codes** (`rankcodes.py`). This is a different formulation from the lane's:
* base row, basis rows $e_1,\dots,e_r$, and extra rows given by functionals $\lambda\in\Lambda_r$, with $|\Lambda_4|=575$;
* depth-first search over the extra rows, with a HiGHS MILP for the oriented column multiplicities (distance and pairwise balance).

Maximal sizes at $c=0$:

| $D$ | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| rank $\le3$ | 4 | 4 | 6 | 4 | 4 | 4 | 6 | 4 | 4 | 4 | 6 |
| rank $\le4$ | 6 | 5 | 8 | 5 | 6 | **5** | 8 | **5** | 6 | **5** | 8 |

All of the lane's entries match. Without balance, $D=10,14$ also give 6.

**A gap in the lane's claim.** The lane's table had no rank-4 entry at $D=11,13,15$, so its claim "8-point rank-$\le4$ clusters need $4\mid D$, $D\le16$ by computation" was **unsupported** there. My values (5) close the gap: the claim is now true for $D\le16$ (exhaustive, floating MILP). For general $D$ only the half-cube case (Prop. 2.1) is proved. That every 8-point rank-4 code is half-cube is not proved.

**2.3 The D=10 "141 codes" are one code, and Lemma 1.1 kills it (proved; this changes the lane's lead).** `six10.py` enumerates **all** 6-point rank-4 codes at $D=10$:
* all multiplicity vectors, exact integer distance filter, then a balance-feasibility MILP;
* reduced modulo XOR-translation (change of base point), coordinate permutations and row permutations.

The result is exactly **one class**. The lane sampled up to 6 MILP solutions per labelled $S$ and counted labellings and orientations separately; its 141 "codes" are this one class. Its rows are
$R_1=\{5,..,9\}$, $R_2=\{2,3,4,8,9\}$, $R_3=\{0,1,4,6,7\}$, $R_4=\{0,1,2,3,5\}$, $R_5=[10]$.

Then $\tfrac12\big((R_1-R_2)-(R_3-R_4)\big)=e_5-e_4$ is an integer vector of weight $2<L+1=5$ in the real span. Lemma 1.1 makes the code unrealisable for **every** orientation and every choice of roots.

At $D=14$ the same enumeration (`sixD.py`) gives 2 classes. Their spans contain vectors of weight 2 and 4, both $<7$, so both are dead (`spanweight.py`). Rank-$\le3$ codes have $\le4$ points at $D=10$ and $D=14$.

> **Proposition K1 (proved: exact enumeration plus Lemma 1.1).** At $D=10$ and $D=14$, $c=0$, every cluster with at least 6 points has rank $\ge5$. With the certified 5-point cluster: $M(10,0)=5$, unless there is a cluster of rank $\ge5$ (6 points, 5 independent differences).

This upgrades the lane's "evidence $M(10,0)=5$ if rank $\le4$" to a theorem. It also removes the lane's main example of a "real-algebraic obstruction beyond dimension counting". The 6-point codes are never realised for a linear reason the lane had already proved, but did not apply to these codes.

**2.4 Types at D=7.** `ktypes.py` finds 2 five-point balanced types up to $S_7\times S_5\times$flip. One is killed by the span vector $e_4-e_1$ (weight 2). The survivor equals $\mathcal T_7$ (canonical forms compared). Corollary 2.2 is confirmed for $D=7$. The $D=8$ (46) and $D=9$ (6) type counts were not recomputed: both enumerations were stopped, for time and the process limit. The $D=8$ realisation test used the lane's type list (data only).


## 3. Independent certification of the claimed clusters

`kcert.py` is my own certifier, and it uses a different gauge from the lane's.
* It fixes the $2D-rL$ coordinates selected by QR column pivoting to 8-digit rationals.
* It runs Newton at 60 digits.
* It runs a Krawczyk test with `mpmath.iv` on a box of radius $10^{-25}$.
* It then recounts at the certified point: all $2^D$ base-form sets, and also the exact **maximum cluster over all $2^D$ points of $N$**, built from the kernel vectors in $\{0,\pm1\}^D$.

| object | Krawczyk | cluster through $Q$ | rank | max cluster on all of $N$ | min nonmember residual | root separation |
|---|---|---|---|---|---|---|
| $D=9,c=0$ (`cert_D9_M4_c0_t0`) | PASS | 4 | 3 | 4 | 0.214 | 0.123 |
| $D=10,c=0$ (`rtq_D10_c0_M5_1`) | PASS | 5 | 4 | 5 | 0.195 | 0.088 |
| $D=12,c=0$ half-cube, 5 solutions (`hc12_lam4` rows 0–3, `hc12_lam-4`) | PASS ×5 | 8 | 4 | 8 | $\ge0.015$ | $\ge0.0074$ |
| $D=10,c=1$ (`rtq_D10_c1_M6_1`, `_10`) | PASS ×2 | 6 | 5 | 6 | 0.075, 0.10 | 0.048, 0.096 |
| $D=12,c=1$ (`rtq_D12_c1_M6_1`) | PASS | 6 | 5 | 6 | 0.021 | 0.0087 |

All claimed certified objects survive. For each, $N$ is squarefree (separation $\gg$ box radius), and no root is real ($\min|\mathrm{Im}|\ge0.016$). For the $D=12$ cluster I recomputed the column multiset: it is antipodal, all non-antipodal distances are 6, and the 8 column types have multiplicities $2,2,2,2,1,1,1,1$. That is the half-cube type with $(2,1)$, as claimed.

The "not a composition" argument (a composed cluster has $\le M(12/\deg h,0)<8$ points) is valid. The 8-point cluster is certified over $\mathbb R$ only, so it is not shown to be a $\mathbb Q(i)$ template; the lane says so itself.


## 4. Counterexample hunt for the claimed rank bounds

**4.1 Vertical lines: the lane's "rank $\le3$" is refuted (certified).**
* On $q_k=iy_k$ ($y_k$ real, $|y_k|$ distinct) the even moments vanish identically. A set $F$ is zero-sum iff $\sum_Fy_k=\sum_Fy_k^3=0$ (here $D=10$, $L=4$).
* `vert.py` takes 4 random sets with a square-norm gauge. Its 17 nondegenerate solutions are all 5-point clusters of rank 4.
* `vcert.py` certifies two of them in $\mathbb R^{10}$ (system: 8 moments, $\sum y^2=10$, one coordinate fixed to a rational).
  * Krawczyk PASS at radius $10^{-25}$.
  * Solution 2: $\min|y|=0.055$, $\min$ gap of $|y|$ $=0.0157$, so the 20 roots $\pm iy_k$ are distinct and non-real.
  * Full-model recount with all moments $m\le4$: the cluster through $Q$ is 5, the rank is 4, and the maximum over all points of $N$ is 5.
  * Certified point, solution 2: $y\approx(1.36636,\,0.60445,\,-1.35061,\,1.09705,\,-0.26318,\,-1.52475,\,1.13694,\,0.05528,\,-0.77066,\,-0.67549)$ (20 digits in `vcert.py` output).
* So "Vertical lines give rank $\le3$" (lane §3.2 and the coordinator summary) is **false**.
* The parameter-count inequality $r\lceil L/2\rceil\le D-1$ is unproved in general; at $D\equiv0\ (4)$ it is an open heuristic.

**4.2 Cube locus.** I found no counterexample. Theorem 3.1 holds (§1), and 60 random locus points show max cluster 8 on all of $N$.

**4.3 Lane's lead: "types of non-negative expected dimension are never realised".**
* *$D=10$ six-point codes:* explained by Lemma 1.1 (§2.3). This is **not** a real-algebraic phenomenon.
* *$D=8$ five-point types:* I ran my own solver (`realise.py`). It uses free orientation, random coordinate slices to make the system square, and LM with 300 starts (plus 150 with free orientation), on the lane's 46 types (data file only). Result:
  * **35/46** types realised with the prescribed orientation, and 40/46 with some orientation. The lane got 30/46.
  * The lane's log (`ext8_a.txt`) lists types 9 and 12 as "not realised in 30 starts". I **Krawczyk-certified** both (`real_one.py`, `kcert.py`):
    * type 9: sep 0.025, nonmember residual $\ge0.054$;
    * type 12: sep 0.15, nonmember residual $\ge0.28$;
    * in both, the cluster through $Q$ is exactly 5, the rank is 4, the max on $N$ is 5, and the orientation matches the type.
  * So the lane's "16 of 46 never realised" is at least partly a **solver artefact**. At most 11 types (6 with any orientation) remain unrealised under my solver, and that is evidence only.
* *$c=1$, $D=14$:* the lane tried 15 random types × 20 starts. That is too thin to show an obstruction. I did not rerun it.

**4.4 Rank 5 at $c=0$.** No rank-5 claim was made, so there was nothing to destroy. The D=8 "near-hits" were sent to degenerate limits by the lane itself. I did not repeat that hunt.


## 5. Literature and novelty

* The lane cites no literature. It uses Plotkin's bound and the Chebyshev-system sign-change lemma as memory-quotes.
* I checked Plotkin by WebSearch:
  * $A(n,d)\le2\lfloor d/(2d-n)\rfloor$ for **even** $d$ with $2d>n\ge d$ (de Launey–Gordon, dmgordon.org/papers/plotkin.pdf; MIT lecture notes give $|C|\le2d/(2d-n)$ for $d>n/2$).
  * The lane applies it with $d=2L$, which is even, so the use is legitimate.
* The sign-change lemma has a one-line proof, which the lane gives.

**Novelty.**
* Lemma 1.1 is a routine extension of (S2) (argument principle / Hermite–Biehler type).
* The genuinely new items are internal to the record:
  * the certified non-composed $D=12$ 8-point cluster;
  * the $D=10$, $D=12$ rank-5 clusters at $c=1$;
  * the exhaustive rank-$\le4$ code tables;
  * Proposition K1 (mine, from the lane's lemma).
* None of this is a new barrier type. The overall line is a **repackaging of the known rank/exceptional-set barrier** in the function-field model: everything still hinges on the open (R$_4$).


## 6. Corrected statement

At $c=0$, no rank-5 cluster is known for $D\le12$. That is a search result (evidence), not a proof. The following are **proved** (exhaustive computation where stated):

* Lemma 1.1.
* Prop. 2.1.
* Thm 3.1: the cluster through $Q$ on the cube locus has exactly 8 points.
* $M\le2$ on horizontal lines for $D\ge7$.
* $M(8,0)\le14$ and $M(9,0)\le6$.
* Rank-$\le4$ code maxima for $6\le D\le16$. In particular, 8-point rank-$\le4$ clusters need $4\mid D$ for $D\le16$; for all $D$ this is proved only for the half-cube type.
* **Prop. K1:** at $D=10,14$, $c=0$, clusters of rank $\le4$ have $\le5$ points.

**Certified (Krawczyk):**
* $M(9,0)\ge4$ and $M(10,0)\ge5$;
* a non-composed 8-point rank-4 half-cube$(2,1)$ cluster at $D=12$, $c=0$, over $\mathbb R$ (not shown to be a $\mathbb Q(i)$ template);
* rank-5 clusters at $c=1$ for $D=10,12$;
* rank-4 clusters on the vertical-line locus at $D=10$, $c=0$. This refutes the lane's vertical-line rank bound.

**Evidence only:** $M(7,0)=4$, $M(9,0)=4$, the absence of rank 5 at $D=8$, and the residual unrealised $D=8$ types.

**Withdrawn:** "vertical lines give rank $\le3$"; "141 never-realised codes of dimension 2" (one code, killed linearly); "16/46 never realised" (at most 11/46 under a stronger solver; two certified counterexamples).

**Files** (`kill_search_checks/`):
* `kcert.py`, `vcert.py`, `real_one.py` — certificates;
* `balclique.py`, `circ.py`, `rankcodes.py`, `six10.py`, `sixD.py`, `spanweight.py`, `ktypes.py` — exact or exhaustive computations;
* `vert.py`, `realise.py`, `cubeloc.py` — numerics;
* logs: `*.txt`.

