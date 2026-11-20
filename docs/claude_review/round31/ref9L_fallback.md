# Referee report ref9L_fallback: Astra Part IX.2 (lattice fallback) and IX.6 "relation to a named hypothesis"

Scripts and outputs are in
`/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref9L_fallback/`
(`local.py`, `trees.py`, `examples.py`, `cg9.py`). They are my own code, with exact integer arithmetic. I did not import
or run Astra's `ix_checks.py`.

## 0. Results

1. **IX.2 is correct, and the defect is larger than Astra says** (Proved, Computed). The lattice form $(L_8^{pc})$
   printed in `expert_note.md` §2 (lines 181-188) has no admissible input at any scale. The reason is local, at the
   primes 2 and 3. So it is a theorem, with $Z=\emptyset$, every $\eta_1$ and $w_1=0$. It is not open. The same holds
   one level up: the set $\mathrm{Bal}^{pc}$, the points of $M_{0,8}(\mathbb Q)$ whose 119 contact integers are
   pairwise coprime, is **empty**. At $p=2$ exactly 5 visible contact integers are even, for every point. At $p=3$ at
   least 2 are divisible by 3. Every statement in the record of the form "restricted to $\mathrm{Bal}^{pc}$" or
   "conclusion $(L_8^{pc})$" is therefore vacuous. This includes theory Prop 2.1(4) and the reduction "every proof
   must prove $(L_8^{pc})$", census Prop 1.3(2) "in particular", expert-note fact 3, and the conclusion of the abc
   criterion Prop 2.6.
2. **The intended statement is recoverable** (Proved). The intended fallback is "coprime outside 2 and 3, with a
   bounded $\{2,3\}$-part". This version follows from $(L_8)$, and its coprimality hypotheses are satisfiable by
   integers. §3 gives exact replacement wording.
3. **IX.6, "relation to a named hypothesis", is correct** (Proved, Computed), once the quantifiers are written out.
   - Cilleruelo–Granville's Conjecture 9 is equivalent to the uniform critical two-class bound at one length
     $C\sqrt q$, and then at every $C$.
   - The little-oh versions at scale $\log q/\log\log q$ are equivalent.
   - So are the little-oh versions at scale $\omega(q)$.
4. **Literature** (snippets only). I found no result giving $o(\log q/\log\log q)$, or any improvement on the trivial
   order $\omega(q)$, for the number of roots of $x^2\equiv a\pmod q$ in an interval of length $\sqrt q$.

## 1. The literal statement and its sources

**Verbatim requirements** (`expert_note.md` lines 181-188, copied from `round18/theory.md` §5(2), lines 576-582):

- Integers: seven finite marks $y=(y_1,\dots,y_7)\in\mathbb Z^7$. The eighth point, $0$, is sent to $\infty$ in
  Kapranov's chart.
- The vector $y$ is primitive, $y\notin\mathbb Z\mathbf 1$, and $[y]\in M_{0,8}\setminus Z$ (so the $y_i$ are
  distinct).
- Height: $|y|\le e^{31w(1+\eta_1)}$.
- Blocks: integers $e_T\in[e^{(1-\eta_1)w},e^{(1+\eta_1)w}]$ for $T\subseteq[7]$ with $3\le|T|\le6$ (98 blocks).
  They are **pairwise coprime** among themselves, and $y\in\Lambda_e$, that is, $e_T\mid y_i-y_j$ for $i,j\in T$.
- Cofactors: $c_{ij}=(y_i-y_j)/\prod_{T\supseteq\{i,j\}}e_T$ for the 21 pairs. They are "pairwise coprime of size
  $e^{w(1\pm\eta_1)}$". Coprimality is required among the $c_{ij}$ themselves. Nothing is said between $c$'s and
  $e$'s; referee item G (`round18/referee_theory.md` §2.8) explicitly allows $c_{ij}$ to share primes with $e_T$
  when $T\supseteq\{i,j\}$.
- No prime is exempt, and nothing says "up to bounded factors".

The Bal-form source (`round18/theory.md` §2.1, line 178) is
$(L_8^{pc})$: $\mathrm{Bal}^{pc}(\eta_1)\cap\{w\ge w_1\}\subset Z$. Here $\mathrm{Bal}^{pc}$ is the set of points
"whose 119 finite boundary heights are pairwise coprime". Prop 2.1(4) motivates it: "every prime divides exactly one
block, the tree at $p$ is a single edge". Census Prop 1.3(2) also uses it: "in particular for every $\Phi$ whose 119
contact integers are pairwise coprime".

## 2. Vacuity (Proved; Computed)

**Lemma 1 (contact splits at 2 and 3).** Let $x\in\mathbb Z^7$ have distinct entries, primitive modulo translation.
Call $T$ ($2\le|T|\le6$) *active at $p$* if $T$ is exactly a residue class of $x\bmod p^t$ for some $t\ge1$; then
$p\mid d_T$ (Lemma 2.2 / Lemma 1.1 of round 18).
- At $p=2$, exactly 5 sets are active, and at least 2 of them have $|T|\ge3$.
- At $p=3$, at least 2 sets are active.

*Proof.* The sets of size $\ge2$ that occur as classes at some level form a laminar family with root $[7]$; the root
is a class only at level 0, by primitivity. When a class first fails to be a single class, it splits into at most $p$
parts. These parts are its children, and singletons are leaves. This gives a rooted tree on 7 leaves in which every
internal node has between 2 and $p$ children.
- For $p=2$ the tree is a full binary tree. It has 6 internal nodes, so 5 non-root active sets. At most 3 of these
  are pairs, since 7 leaves give at most 3 cherries, so at least 2 have $|T|\ge3$.
- For $p=3$, $6=\sum(\text{children}-1)\le2\,(\#\text{internal})$, so there are at least 3 internal nodes and at
  least 2 active sets. $\square$

(This is the record's own bound "each prime divides at most 5 visible blocks", census Thm E, which is attained at
$p=2$ for every point.)

**Corollary.** $\mathrm{Bal}^{pc}=\emptyset$. The Bal-form $(L_8^{pc})$, the restriction "pc" in
$(L_8^c)^{pc}$, and every "in particular" for pairwise coprime contact integers are vacuous.

**Lattice form, Astra's argument checked.** By pairwise coprimality, $\prod_{T\supseteq\{i,j\}}e_T\mid y_i-y_j$, so
the $c_{ij}$ are integers. At most one $e_T$ is even.
- *All $e_T$ odd.* Seven integers have at least $\binom43+\binom32=9$ pairs of equal parity, so at least 9 $c_{ij}$
  are even.
- *$e_T$ even, $|T|\ge4$.* Four $T$-marks lie in one class mod $2^v$, $v=v_2(e_T)$, which has only two lifts mod
  $2^{v+1}$. That gives at least 2 pairs with $2\mid c_{ij}$.
- *$e_T$ even, $|T|=3$.* The four marks outside $T$ give at least 2 equal-parity pairs whose block products are odd.

In every case at least two cofactors are even, which contradicts pairwise coprimality. The argument uses neither
height, balance, primitivity nor $Z$. **Correct.**

**Computed** (`local.py`): an exhaustive local minimisation over the $p$-adic configurations, with at most one
$e_T$ divisible by $p$ and arbitrary valuation $a\le3$.
- The minimum number of even $c_{ij}$ is **5**.
- The minimum number of $c_{ij}$ divisible by 3 is **2**.

This confirms Astra's sharper count. **Computed** (`trees.py`, 20,000 random primitive 7-tuples):
- At $p=2$ the number of active sets is always 5, with $\ge2$ of size $\ge3$.
- At $p=3$ it is $\ge2$, with distribution 2..5.
- At $p=5$ it can be 1.

The canonical-contact version is vacuous too, since $c_{ij}=\pm d_{\{i,j\}}$ and $e_T=d_T$ there.

**The census does not contradict vacuity.** `round18/census.md` §6.2 prescribes pairwise coprime prime *blocks*
$n_T$ of $(L_8)$-**data**. Those data have residues $t_{ij}$, with $H_\Sigma\approx464w$, and the residues absorb
the 2- and 3-adic trees. They are not points of $\mathrm{Bal}^{pc}$, and they are not inputs of the lattice form.

**Computed** (`examples.py`): census-type vectors of $\Lambda_e$, with 98 distinct primes $e_T\in[1000,3000]$ and
random CRT lattice vectors. They have 9-11 even cofactors and 5-7 cofactors divisible by 3. Likewise
`README_VIIIL.md`'s $(L_8)^{\rm sq}$ concerns blocks with arbitrary cofactors, which is a non-vacuous class, and is
unaffected.

**What the source intended.** The source meant that "each prime divides one block, so the chain condition and the
residues are vacuous". A single-edge tree at $p$ requires both vertex degrees to be $\le p+1$:
- $|T|\in\{3,4,5\}$ needs $p\ge5$;
- $|T|\in\{2,6\}$ needs $p\ge7$.

So the intention is realisable only away from 2 and 3. Coprimality must be imposed outside $\{2,3\}$, and the
$\{2,3\}$-part must be paid for (bounded, or carried by residues).

## 3. Corrected fallback (Proved) and replacement wording

**Replacement for `expert_note.md` lines 181-188:**

> **$(L_8^{pc})$, lattice form (corrected, round 31 IX.2).** Let $(c,\eta_0,Z,C_0)$ be as in $(L_8)$ (anchor $0$,
> unsymmetrised $Z$). For every $B_0\ge1$ there is $w_1=w_1(B_0)$ with the following property. Let
> $e_T\in[e^{(1-\eta_1)w},e^{(1+\eta_1)w}]$ ($3\le|T|\le6$, $\eta_1=\eta_0/2$) be integers that are pairwise
> coprime outside the primes $2,3$, and let $w\ge w_1$. Then $\Lambda_e$ has no vector $y$ with distinct entries such
> that all of the following hold:
> * $[y]\notin Z$;
> * the pair cofactors $c_{ij}=(y_i-y_j)/\prod_{T\supseteq\{i,j\}}e_T$ are integers that are pairwise coprime outside
>   $2,3$, of size $e^{w(1\pm\eta_1)}$;
> * the $\{2,3\}$-part of $\prod_{i<j}(y_i-y_j)$ is at most $B_0$.
>
> The literal version with full pairwise coprimality has empty input: at $p=2$ at least five $c_{ij}$ are even, and
> at $p=3$ at least two are divisible by $3$. Likewise $\mathrm{Bal}^{pc}=\emptyset$, since every point has exactly
> five even contact integers.

*Proof that $(L_8)$ implies it.*
- **Datum.** Take $P_0=(1,0)^T$ and $P_j=(y_j,1)^T$. Then $t_{0j}=Y_j=1$.
- **Blocks.** $n_T$ is the $\{2,3\}$-free part of $e_T$ for $|T|\ge3$; $n_{ij}$ is the $\{2,3\}$-free part of
  $|c_{ij}|$; and $n_{[7]}=1$, as do singletons.
- **Residues.** Then $t_{ij}=\pm(\{2,3\}\text{-part of }y_i-y_j)$, so $H_\Sigma=\log(\{2,3\}\text{-part of}\prod(y_i-y_j))\le\log B_0$.
- **Chain condition.**
  - At 2 and 3 no block is divisible, so the condition is empty.
  - At $p\ge5$ at most one $n_T$ with $|T|\ge3$ and at most one $n_{ij}$ are divisible. These are nested or
    disjoint, which is allowed (two chains with disjoint unions).
  - The crossing case $i\in T$, $j\notin T$ is impossible. It gives $y_k\equiv y_i\equiv y_j\bmod p$ for all
    $k\in T$, so $p\mid c_{kj}$ for the $\ge2$ indices $k\in T\setminus\{i\}$, and then $p$ divides at least three
    cofactors.
- **Balance.** Every visible $n_T$ satisfies $|\log n_T-w|\le\eta_1w+\log B_0\le\eta_0w$ once
  $w\ge2\log B_0/\eta_0$.
- **Conclusion.** $(L_8)$ gives $\log B_0\ge H_\Sigma\ge cw-C_0$. So
  $w_1(B_0)=\max(2\log B_0/\eta_0,(\log B_0+C_0)/c)+1$ works.

The height bound and primitivity are not needed. $\square$

**Not locally empty** (Computed, `examples.py`). Take $y=(0,5,10,107,267,69,-229)$, $e_{\{1,2,3\}}=5$ and all other
$e_T=1$. The 21 cofactors are pairwise coprime outside $\{2,3\}$. So the coprimality pattern is satisfiable by
integers; balance at large $w$ is the content of the conjecture. At $p=5$ every input is forced to have a class of
size 3 or 4 carried by one $e_T$, since 7 marks in 5 classes otherwise give two divisible cofactors. This is a
constraint, not an obstruction. The heuristic count $98-108=-10$ is unchanged up to $O(1)$.

**Other places needing rewording** (Gap in the record, not in Astra).
- **Sandwich, line 156.** The last arrow $\Rightarrow(L_8^{pc})$ is trivially true. Replace it with the corrected
  statement; $(L_8)\Rightarrow$ corrected is proved above.
- **Line 166.** "pc restricts to pairwise coprime $d_T$" should read "pairwise coprime outside $\{2,3\}$, bounded
  $\{2,3\}$-part".
- **Line 172, "open".** For the literal statement, "open" should read "vacuous (proved)". The corrected statement is
  open.
- **Fact 3, lines 231-233, and theory Prop 2.1(4).** On a corrected $\mathrm{Bal}^{pc}$ the chain condition is
  vacuous at $p\ge5$, and at 2 and 3 the structure is absorbed by residues of total size $\le\log B_0$. The
  reduction "every proof of $(L_8)$ must prove [corrected]" then holds.
- **Question A, lines 283-286.** Replace "its lattice form $(L_8^{pc})$" with "the corrected lattice form (coprime
  outside 2, 3)".
- **Abc criterion, line 620, and theory Prop 2.6; overview lines 13, 81, 1035-1039, 1058, 1302-1312.** As stated,
  the conclusion is vacuously true, so the "uncapped conditional lane" has an empty target. Whether Prop 2.6 survives
  on the corrected class is **Heuristic**: its gcd identity $\gcd(A,C)=\prod d_T^{\mu_T}$ should hold up to a factor
  dividing $B_0^{O(1)}$, which costs $O(\log B_0)$. I did not check this.
- **Census Prop 1.3(2) "in particular".** This is vacuous; delete it.
- **Unaffected.** Theory Thm 2.5, Prop 2.3 and the $(L_8)^{\rm sq}$ results of `README_VIIIL.md` do not use
  $\mathrm{Bal}^{pc}$.

**"Is $(L_8^{pc})$ Rudin-hard?"** (judge question 4, `README_VIIIL.md` line 467). For the literal statement the
answer is trivial: it is a theorem. Any implication from it to Rudin would be an unconditional proof of Rudin.
Astra's remark is correct. The question must be re-posed for the corrected form. That form only covers
$H_\Sigma\le\log B_0=O(1)$, while the $(L_8)^{\rm sq}$ route to Rudin needs the linear bound $cw$ against arbitrary
cofactors. So no implication is known (**Open**). Items in `README_VIIIL.md` lines 321 and 399 should name the
corrected statement.

## 4. IX.6: equivalence with Cilleruelo–Granville Conjecture 9 (Proved; Computed)

**Definitions.**
- $N(b)$ is the maximum, over $a$ and over half-open intervals of length $\sqrt b$, of the number of roots of
  $x^2\equiv a\pmod b$ in the interval. Conjecture 9 (verbatim in `lit8L_verify.md` C1) says $\sup_b N(b)<\infty$.
- $B_C(q)$ is the maximum, over choices of at most two residues per $p^f\parallel q$ and over intervals of length
  $C\sqrt q$, of the number of admissible integers.

**(a) $B_C(q)\le(C+2)\,N(4q)$.**
- By CRT choose $a,b$ with $a\equiv r_1$ and $b\equiv r_2\pmod{p^f}$, repeating a residue when only one is
  prescribed. Every admissible $y$ has $q\mid(y-a)(y-b)$.
- So $(2y-a-b)^2-(a-b)^2=4(y-a)(y-b)\equiv0\pmod{4q}$.
- The map $y\mapsto2y-a-b$ is injective and sends the interval to one of length $C\sqrt{4q}$. Cover that by
  $\lceil C\rceil+1$ half-open intervals of length $\sqrt{4q}$.

Nothing about units or coalescence is used. **Computed:** all 61,833 admissible $y$ for $q<3000$ (random residue
choices, including coalescent ones) satisfy the congruence.

**(b) $N(b)\le2\,B_1(Q_2)$, where $Q_2\le b$ and $\omega(Q_2)\le\omega(b)$.** This is the valuation division of
VIII.2 (refereed in `ref8L_rudin.md` §5).
- $u$ divides every root.
- Modulo each $p^f\parallel Q_2$ the roots divided by $u$ lie in $\le2$ classes.
- The interval shrinks to length $\le\sqrt{Q_1}\le\sqrt{2Q_2}$, which two intervals of length $\sqrt{Q_2}$ cover.

$Q_2=1$ is trivial. **Computed** (`cg9.py`): all 326,907 square classes $(b,a)$ with $b\le1500$ were checked for the
divisibility, the $\le2$ classes and the length bound; there were 0 failures.

**Uniform equivalence.**
- Conjecture 9 $\Rightarrow$ $B_C$ bounded for every $C$, by (a).
- $B_c$ bounded for one $c>0$ $\Rightarrow$ Conjecture 9, by (b) together with
  $B_1\le(\lceil1/c\rceil+1)B_c$.

So Conjecture 9 $\Leftrightarrow$ TR(1/2) at one $C$ $\Leftrightarrow$ TR(1/2) at all $C$. Here TR(1/2) is in
`README_VIIIL.md`'s notation, at length $\sqrt q$.

**Rate equivalence.** Let $f(x)=\log x/\log\log x$, which is increasing for $x\ge16$.
- (a) gives $B_C(q)\le(C+2)N(4q)$ and $f(4q)\sim f(q)$.
- For (b), given $\delta>0$ pick $Q_\delta\ge16$ with $B_1(Q)\le\delta f(Q)$ for $Q\ge Q_\delta$, and put
  $K_\delta=\max_{Q<Q_\delta}B_1(Q)<\infty$. Then $N(b)\le2\max(K_\delta,\delta f(b))$, so
  $\limsup N/f\le2\delta$.

Hence $N=o(f)\Leftrightarrow B_C=o(f)$ for one, equivalently every, $C$. Astra's monotonicity remark is exactly this.

The same argument works at scale $\omega$. It needs $B_1(Q)\ll1+\omega(Q)$ for bounded $\omega$. That follows from
Astra's $P_*$ criterion $B\ll_C1+\log q/\log P_*$, which I re-derived: the averaging over the $b\ge B/2$ marks that
agree at $P_*$ gives $b(\log P_*-2\log C)\le W-2\log C$, and $\log q/\log P_*\le\omega(q)$. So
$N=o(\omega)\Leftrightarrow B_C=o(\omega)$, which is the project's $o(s)$ question for the twin.

No separation from $(L_8)$, and no circle consequence, is claimed, and none follows.

## 5. Consequence for the record and the literature

The split twin's uniform bound **is** Cilleruelo–Granville's Conjecture 9, a published conjecture that implies Rudin's
conjecture (their Thm 11). Its rate version is the natural growth question there. The same reduction gives the trivial
bound $N(b)\ll1+\omega(b)\ll\log b/\log\log b$.

**Literature search** (WebSearch result snippets only; no paper was downloaded or read in full). Three searches were
run:
- on roots of $x^2\equiv a\bmod q$ in intervals of length $\sqrt q$ with a $\log q/\log\log q$ factor;
- on follow-ups to Cilleruelo–Granville (arXiv math/0608109);
- on Cilleruelo–Garaev, Chan–Shparlinski and Cilleruelo–Garaev–Ostafe–Shparlinski concentration results.

They returned:
- equidistribution results for $x^2\equiv p\pmod q$ on average over $q$ (arXiv 2009.03460; Math. Z. 2024);
- the root count $\rho(q)\le2^{\omega(q)+2}$ (arXiv 1401.2355);
- modular-hyperbola concentration bounds, mainly for prime moduli and boxes below $p^{1/4}$ (arXiv 1007.1526,
  1103.2879, 1203.0017).

None states an $o(\log q/\log\log q)$ or $o(\omega(q))$ bound at length $\sqrt q$ for general $q$, and none states
anything beyond Cilleruelo–Granville Theorem 10. Theorem 10 gives the Liouville-type span bound
$x_k-x_1>b^{1/2-1/(2\ell)}$; it gives $O(1/\delta)$ roots at length $b^{1/2-\delta}$ and nothing at length
$\sqrt b$. **I found no result on the rate question.** This is a bounded search over snippets, not a proof of absence.

## 6. Verdict table

| claim | verdict |
|---|---|
| IX.2: literal lattice $(L_8^{pc})$ has empty input (2 and 3) | CORRECT (Proved, Computed), and stronger: $\mathrm{Bal}^{pc}=\emptyset$ |
| IX.2: $\ge5$ even, $\ge2$ multiples of 3 | CORRECT (Computed minimum 5 / 2; tree proof) |
| IX.2: canonical-contact version empty | CORRECT |
| IX.2 remark: literal $\Rightarrow$ Rudin would prove Rudin | CORRECT; the judge's question must be re-posed for the corrected form |
| IX.6: Conjecture 9 $\Leftrightarrow$ uniform critical twin | CORRECT with quantifiers (§4) |
| IX.6: little-oh equivalence at $\log q/\log\log q$ | CORRECT (monotonicity for $Q\ge16$, finite $K_\delta$); also at scale $\omega$ |
