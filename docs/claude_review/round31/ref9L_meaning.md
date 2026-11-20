# ref9L_meaning: what Part IX of the two-class thread changes, and what to ask next

Analyst block for round 31, Part IX of `astra_lemma_IX/answer.md` (lines 1931-2149, "The critical two-class problem on
smooth parents"), read in full, with `final9.txt`, `reply9.txt`, and the parallel thread's Part IX
(`astra_free_IX/answer.md` from line 1955, `final9.txt`; **unverified**, refereed separately). Record read: expert note
§§1, 2, 5 (items 1-12 in full), `overview/overview.md` §5.4, `round31/README_III.md` §§1, 3.1-3.2, `README_VIII.md`
(P2 statement), `README_VIIIL.md` in full, `round9/L1L2.md` §5 and `round9/referee_L1L2.md` §4, `round10/additive.md`
§§5-6. I also read, as factual input for the analysis and not as my own verification, the block reports already
present in this directory: `ref9L_pool.md`, `ref9L_powers.md`, `ref9L_samenorm.md`, `ref9L_relax.md`,
`ref9L_additive.md`, `ref9L_fallback.md` (`ref9L_nonlinear.md` was a skeleton when I read it). I did not re-referee
the proofs. My own code (no import of Astra's script) is in `ref9L_meaning/`: `budget.py` (entropy budgets,
output `budget.out`) and `pc2.py` (the parity obstruction of IX.2, output `pc2.out`).

Labels: **Proved** (checked by me), **Computed** (reproduced by my code), **Heuristic**, **Gap**, **False**,
**Record** (taken from a refereed block report or the record, not rechecked by me).

---

## 0. Results (answers to the seven questions)

1. **Random-parent rate (Record: proved in `ref9L_pool.md`, `ref9L_powers.md`).** What may be recorded is a rate
   *on average over a prime pool*: for a uniformly random $s$-subset of a pool of $K\ge e^{40}s$ primes in $[P,2P]$,
   outside a fraction $\exp(-sJ/20)$, $J=\log(K/s)$, every critical two-class configuration and every critical arc of
   the normalised norm $q$ has $O(1+\log P/J)$ marks. At $P=As\log s$ ($A\ge 2e^{40}$ for split pools) this is
   $O(\log s)=O(\log\log R)$ where the worst-case bound is of order $s\asymp\log R/\log\log R$. It is **not** "the
   uniform bound for almost all smooth parents": at the smooth critical scale the bound is $O(\log s)$, not $O(1)$;
   $O(1)$ needs $P\ge s^{1+\epsilon}$, which the record already had on circles (Part III, III.9, $D>1$). The new items
   are the endpoint $D=1$ (where III.9's $M_D$ blows up), the twin, and the simpler three-mark proof.
   The exceptional fraction is only $R^{-O(J/\log P)}=R^{-o(1)}$ at the smooth scale, and the exceptional *number*
   may be as large as $(\#\text{parents})^{0.95}$ (Computed). The theorem does **not** cover the profile extremal for
   Theorem 1.1 (all split primes up to $\approx6M\log M$): that profile is a single dense parent with $K=s$, $J=0$,
   and the method needs $K/s\ge e^{c}$ with $c\ge5\log4\approx6.9$ even after optimisation (§1.4). It does not
   reduce the growth problem to anything easier: for comparable-prime squarefree parents the worst-case rate is
   *exactly* the statement that the exceptional parents carry no large cluster, and nothing makes that set empty.
2. **The input fakes lack.** Both the same-norm comparison and the IX.4 relaxation isolate the same input: the
   *exact* additive identity among genuine differences with integer cofactors (twin: $t_{12}A_{12}+t_{23}A_{23}=t_{13}A_{13}$;
   circle: reduced chords are Gaussian integers satisfying the cocycle). The record **had** the cocycle as an input
   (round 9 `L1L2.md` Prop 5.1(a); round 10 `additive.md` §6, classified "capped ($A_0$)"; overview §5.4 rows 15-16;
   expert note §5 item 6 for abc on triangle relations), but it only analysed its consequences at primes $q\nmid N$
   (transitivity, i.e. residue partitions) and abc on it (vacuous). The round-9 referee had flagged that "the cocycle
   is also an exact identity in $\mathbb Z[i]$ ... which (a) does not address". Part IX is the first result that
   uses that exact content, and the first that shows it separates genuine from fake on a common norm. The classification
   "cocycle: capped" must be narrowed (§2.2). In the twin, ref9L_relax proves the triangle identities are equivalent
   (for $|d_{ij}|<q/3$) to the signed pair residues modulo $q$; so the missing input is "local signed data at the
   parent's own primes, coming from one global integer".
   **Why average but not worst case (Proved, §2.3):** the union bound gains $k\log(K/s)$ with $k\ge s/5$ from the
   fact that the triangle *determines* a block of the parent; for a fixed parent the same equation is a three-term
   squarefree $S$-unit equation $aX+bY=cZ$ with $S$ the parent's primes and $|abc|\le\exp(O(W/b))$. There is no
   probability to spend; abc is vacuous (squarefree terms); the known upper counts (Evertse, Beukers-Schlickewei) are
   $\exp(O(s))$, no better than the $4^s$ trivial count; and the known lower constructions (Erdős-Stewart-Tijdeman,
   Konyagin-Soundararajan, snippets only) show that for suitable dense $S$ the unrestricted equation has
   $\exp(s^{2-\sqrt2-\epsilon})$ solutions. The triangle therefore gives a rate exactly as far as the parent carries
   entropy ($J\gg1$), and nothing when $J=0$.
3. **Deterministic results (Record: proved in `ref9L_additive.md`).** $B\ll_C\sqrt{|A+A|}$ (twin, new) and
   $M\ll_C\sqrt{|A+A|}$ (circle, not new in substance) are saturated by every fake (Sidon signatures), so they kill no
   fake. IX.1 (no affine 3-flat of signatures at $L<q^{2/3}/4$) and IX.5 ($n^{1-2^{1-2d}}$, $n^{7/8}$ at $d=2$) kill the
   *linear / low-complexity* families: Sylvester-Hadamard (Walsh) signature profiles in the twin, and every integer
   transfer of the characteristic-2 hyperplane mechanisms (Part V, Part VI binary families, the attack referee's char-2
   Hadamard family) that keeps linear bits: such a transfer contains affine 3-flats and is impossible over $\mathbb Z$
   beyond 7 marks at critical length. They do not reach Paley/Legendre profiles (Sidon, high degree), nor Gold-type
   quadratic bits of full rank ($x^3$ cores). The modulus-alone criterion is new for the twin, old for circles
   (round 8 `stability.md` Cor 3.2), and does not reach balanced smooth parents. New relative to the record: the twin
   statements and IX.5; literature not checked.
4. **IX.2.** Under the outcome "vacuous" (supported by `ref9L_fallback.md` and by my check, `pc2.out`: every
   admissible choice forces $\ge5$ even pair cofactors), the Question A box must drop "or at least its lattice form
   $(L_8^{pc})$"; fact 3 of §2 ("every proof of $(L_8)$ must prove $(L_8^{pc})$") becomes true but empty; the judge's
   open question "is $(L_8^{pc})$ Rudin-hard?" is answered **no** in the trivial sense (a theorem with empty input
   implies Rudin only if Rudin is provable outright) and must be re-asked for the corrected fallback ("pairwise coprime
   outside $\{2,3\}$, with bounded $\{2,3\}$-part"). Under the outcome "not vacuous" (if another referee finds an
   intended reading with nonempty input), the box and the question stand as printed, with that reading made explicit.
5. **The twin.** The record gains a clean sister problem with a name: the **rate** form of Cilleruelo-Granville's
   Conjecture 9 (roots of $x^2\equiv a\bmod q$ in an interval of length $\sqrt q$: $o(\log q/\log\log q)$, or
   $o(\omega(q))$). The equivalence (IX.6, Record: `ref9L_fallback.md`) holds at the uniform and both little-oh
   scales. No transfer in either direction is known, and I see none (the two problems live on the split and the
   anisotropic form of one torus; there is no $\mathbb Q$-homomorphism, and the classical root/representation
   correspondence does not preserve short intervals). Every *method* of Part IX ports both ways; no *statement* does.
6. **Next steps** (§6): the honest probabilities of moving the worst-case rate are all $\le5\%$. Counting over
   exponent patterns cannot remove the exceptional parents: for a fixed parent the union bound is the trivial
   inequality "#clusters ≤ #candidates", and the pattern-plus-cofactor entropy, $(m-1)s\log2+W/2$, already exceeds the
   entropy of the whole pool, $s(1+J)$, for every $m\ge2$ (Proved, §6.1).
7. **Wording** (§7). The thread should get **one** more targeted question (worst-case twin for squarefree
   comparable-prime moduli, or an explicit obstruction there), and then close unless that question moves.

---

## 1. The random-parent rate

### 1.1 What may be recorded (Record: proved by `ref9L_pool.md` and `ref9L_powers.md`; quantifiers as I read them)

> **Theorem (random parents, round 31 Part IX; almost-all over a prime pool, not a worst-case statement).** Fix
> $C>0$. There is $s_0(C)$ such that for all $s\ge s_0(C)$, all $P$, and every set $\mathcal P$ of $K\ge e^{40}s$
> primes in $[P,2P]$, the following holds with $J=\log(K/s)$ and $b-1=\max(2,\lceil100\log(2P)/J\rceil)$. Let $S$ be
> a uniformly random $s$-element subset of $\mathcal P$ and $q=\prod_{p\in S}p$. Outside an event of probability at
> most $\exp(-sJ/20)$: (i) for every choice of at most two residues modulo each $p\in S$ and every interval of length
> $C\sqrt q$, at most $b-1$ integers lie in the prescribed classes; (ii) if $\mathcal P$ consists of primes
> $\equiv1\ (4)$, every arc of length $Cq^{1/4}$ of $x^2+y^2=q$ carries $O(b)$ lattice points, and the same holds for
> clusters on $x^2+y^2=q|G|^2$ all of whose points are divisible by one Gaussian integer $G$. The residues, the
> interval and the arc are chosen after $S$. Variants: primes below $X\le s^A$ (bound $O_{A,C}(1+\log X/J)$, $J\ge
> J_0(A)$); exponents $1\le e_p\le E$ chosen after $S$ (twin: $J\ge40E\log(1+3E)$, failure $\exp(-sJ/(8E))$, bound
> $O(E^2(1+\log P/J))$; circle: $J\gtrsim E^2\log E$, failure $\exp(-\Omega(sJ/E^2))$, bound $O(E^5(1+\log P/J))$;
> `ref9L_powers.md`). Nothing is uniform in $E$.

Corollaries. With the pool all split primes in $[P,2P]$, $P=As\log s$, $A\ge2e^{40}$ fixed: $J=\log(A/2)+o(1)$ is a
constant, $\log q\sim s\log s$, and almost every parent has $B_{\max},M_{\max}\ll_C\log s$. With $P=s^{1+\epsilon}$:
$J\sim\epsilon\log s$ and $B_{\max},M_{\max}=O_{C,\epsilon}(1)$.

What it is: a rate $o(\omega(q))$, in fact $O(\log\log q)$, for a measure on parents of the right height.
What it is not: a statement about any individual $q$ or $N$; a statement about clusters whose normalised norm is not a
member of the family (an extra non-constant split prime destroys coverage, as in III.9); a uniform bound at the smooth
scale.

### 1.2 Comparison with the record's earlier random-parent theorem (III.9)

| | III.9 (round 31 Part III, free thread) | Part IX.0 |
|---|---|---|
| parents | $r$ split primes from $[r^D,2r^D]$, $D>1$ fixed | $s$ primes from any pool of $K\ge e^{40}s$ primes in $[P,2P]$ |
| height | $\log N\sim Dr\log r$ | $\log q\sim s\log P$, including $P=As\log s$ |
| bound | $M_D$: 1 ($D>4$), 2, $O(1)$ ($D>2$), 37 ($D=3/2$), finite for $D>1$, $\to\infty$ as $D\downarrow1$ | $O(1+\log P/J)$: $O(1)$ at $D>1$, $O(\log s)$ at $D=1$ |
| exceptional fraction | $\exp(-cr\log r)=N^{-c/D}$, $c=c(D,M,C)$ | $\exp(-sJ/20)$ |
| problems | circles | twin (all residues) and circles; bounded exponents both |
| proof | pair count ($D>2$); quartet exposure | one triangle with a deterministic heavy block |

So for circles Part IX adds the endpoint $D=1$ and a cleaner proof; the circle statement and the exposure lemma were
already in the parallel thread's unrefereed notes (`astra_free_III/iii_graph_notes.md`, per `ref9L_pool.md`). For the
twin it is new.

### 1.3 How sparse are the exceptional parents (Computed, `budget.out`)

The number of parents is $\binom Ks$, with $\log\binom Ks=s(1+J)-O(\log s)$ for $K\gg s$. The theorem allows
$\exp(s(1+J)-sJ/20)$ exceptional parents, i.e. $(\#\text{parents})^{0.951}$ at $J=40$ and at $J=60$, for every $s$
tested ($10^3$ to $10^9$). In terms of $R=\sqrt q$: $\exp(-sJ/20)=R^{-J/(10\log P)}$, which is $R^{-0.06}$ to $R^{-0.08}$ in
the computed range and $R^{-o(1)}$ asymptotically at $P=As\log s$ (polynomially small in $R$ only at $P=s^{1+\epsilon}$).
So "almost all" means "all but a fraction that is small on the scale of the pool, but subpolynomial in $R$, and
numerous". Nothing in the proof makes the exceptional set structured (it is defined by the existence of a solution to a
triangle equation), and nothing suggests it is empty.

The union-bound exponent per slot (with $b-1=100\log(2P)/J$) is $\log(\text{labels})+9J/200-J/5$; it is $\le-J/20$
exactly when $J\ge\log(\text{labels})/0.105$, i.e. $J\ge18.6$ for 7 labels (circle) and $J\ge13.2$ for 4 labels (twin)
(Computed for $J=18,18.5,20,40$ in `budget.out`). Astra's $J\ge40$ is safe. Note also that $J\ge40$ forces
$A\approx2e^{40}$, so "$\log q\sim s\log s$" holds only once $\log s\gg40$: the smooth-scale corollary is purely
asymptotic.

### 1.4 Does it cover the profile extremal for Theorem 1.1?

The extremal profile is squarefree, $s\approx3M$ split primes, all split primes up to $X\approx6M\log M$, so that
$\log N\approx3M\log M$. Part IX covers the *height scale* of that profile ($\log q\sim s\log s$) but **not the profile**:

* That profile is one parent with $K=s$ (the pool is exhausted), so $J=0$ and the theorem is empty.
* Even as a member of a larger pool (e.g. all split primes below $e^{40}s\log s$) it is one individual subset, about
  which an almost-all statement says nothing.
* **Lower limit of the method (Proved).** Any version of the argument that exposes slot labels (at least 4 per slot in
  the twin, 7 on circles) and gains only on a heavy block of $k\ge s/5$ slots needs
  $s\log4<kJ\le sJ$ at least, and with the guaranteed $k\ge s/5$, $J>5\log4\approx6.9$. Even the best conceivable
  heavy block $k\le s$ needs $J>\log4$. So dense pools ($K\le e^{6.9}s$ with the $s/5$ guarantee) are out of reach of
  the exposure method, whatever the cofactor bookkeeping.

### 1.5 Does it turn the growth problem into a problem about exceptional parents?

Only tautologically. For the class "squarefree, $s$ primes in $[P,2P]$", the worst-case rate is equivalent to "no
exceptional parent carries a cluster of size $\gg\log s$". That set may have $(\#\text{parents})^{0.95}$ elements, it
contains every dense configuration (J small), and the general problem also contains unbounded exponents, mixed prime
sizes, and the Theorem 1.1 extremal profile. No reduction in difficulty results. The useful content is diagnostic: on
typical parents of the smooth critical height the truth is $O(\log s)$ (and on typical bin-parents $O(1)$), so the
barrier at order $s$ is a worst-case phenomenon, not a height phenomenon.

---

## 2. The input that fakes lack

### 2.1 What the two results identify

* **Same-norm comparison** (Record: proved, `ref9L_samenorm.md`). On arbitrarily large norms $N=\prod_{p\in S}p^2$
  ($S$ a typical $s$-subset of a bin of multiplicative width $l^{-3/4}$ of split primes in $[P,2P]$, $P=l^{1.99}$), the
  P2 fakes have $\asymp s$ marks satisfying every $A_{2M}$ input, and every genuine critical arc has $O_C(1)$ points.
  The fake half holds for **every** $s$-subset of the bin; the genuine half for all but $\exp(-\Omega(s\log s))$. The
  genuine half is purely asymptotic ($\log l\gtrsim2.7\cdot10^4$ with the constants as written). These are typical
  bin-parents, not typical dyadic parents (a dyadic sample lands in one bin with probability $\approx e^{-0.75s\log s}$),
  so the Part III judge's question about typical dyadic parents remains open as asked.
* **Relaxation IX.4** (Record: proved, `ref9L_relax.md`). A table of small signed integer edges with correct gcds,
  integral cofactors forming a cocycle, all local valuation partitions, Sidon signatures and an approximate real
  configuration reaches $m\sim\frac12\log q/\log\log q$, the larger-sieve constant to leading order, and has no integer
  realisation. All $\binom{43}{3}=12{,}341$ triangle identities fail at $m=43$ (the count is "all triples" and carries no
  information by itself). For $|d_{ij}|<q/3$ the triangle identities are **equivalent** to the signed pair residues
  $d_{ij}\equiv\rho_j-\rho_i\pmod q$; the relaxation fails these maximally ($d_{ij}^2\bmod p$ takes 426-461 values over
  462 separated pairs; a genuine configuration has one).

Together: inputs that read a configuration through magnitudes, gcds, valuation partitions, cut signatures, monomial
gaps, rational characters and point residues at moduli $\le2M$ cannot see the difference between $O(1)$ and $\asymp s$;
the exact additive identity among genuine differences, with integer cofactors, can, on average.

On circles the corresponding input is: the reduced chords $b_{ij}$ are Gaussian integers and
$z_j-z_i=2i\eta H_{ij}t_{ij}$ with $t_{ij}\in\mathbb Z$, with the cocycle exactly. Equivalently, the local data at the
parent's own primes are reductions of one Gaussian integer per mark. The parallel thread's Part IX (unverified)
approaches the same point from the fake side: its IX.3 builds fakes whose residues factor through norm-correct prime
residues at every odd prime power, including same-level pairs at parent primes (so per-prime witnesses are imitable),
and its IX.7 shows that *if* the witnesses come from one rational Gaussian lift per label, one-unit tests force
genuine angles. Its closing question ("can the local witnesses on a worst parent be forced to share a rational lift?")
is the circle form of the same missing input. Both threads converge on: **local data must come from one global
integer**, and no local test certifies that.

### 2.2 Had the record identified this input class?

Partly, and it misfiled it.

* `round9/L1L2.md` Prop 5.1(a): the cocycle $g_{xw}\beta_{xw}=g_{xy}\beta_{xy}+g_{yw}\beta_{yw}$; "its only consequence at
  an odd prime $q$ not dividing $N$ is transitivity", i.e. point-residue partitions. The referee
  (`referee_L1L2.md` §4) added: "The cocycle is also an exact identity in $\mathbb Z[i]$ with archimedean content, which
  (a) does not address."
* `round10/additive.md` §6 (refereed): "sums, midpoints, Sidon, chord lengths, $\sigma_{xy}$ cocycle, triangle areas,
  4-point relations: capped ($A_0$)". Overview §5.4 row 16 lists "the $\sigma$-cocycle" among tools that "reduce to pair
  data on one circle".
* Expert note §5 item 6 / overview row 15: abc on triangle relations among chords is vacuous on squarefree profiles.

Part IX shows that the exact identity with integer cofactors at the parent's primes is **not** captured by those
consequences: fakes satisfy all of them and fail the identity (they have no integer cofactors). The "capped" entry is
correct only for the listed consequences (local consequences at $q\nmid N$, triangle areas, abc). Proposed correction
(§7.3): the cocycle *as an exact identity in $\mathbb Z[i]$ (resp. $\mathbb Z$) with integer cofactors at the primes of
$N$* is **not capped**; it gives a rate on random parents and its fixed-parent content is open. Row 16's "reduce to
pair data on one circle" remains true as a statement about known worst-case arguments (Part IX gives nothing on one
fixed circle), but should no longer be read as a cap on the identity.

### 2.3 Why the triangle gives a rate on average and not in the worst case (Proved, as an analysis of the proof)

*Where randomness enters.* Steps (a) averaging to find a triple with $\log|t_{12}t_{23}t_{13}|\le3\log C+3W/(2(b-1))$, (b) the
deterministic heavy block of $\ge s/5$ primes, (c) the identity (IX.0a), are all deterministic and hold for every
parent. Randomness enters only in (d): conditional on the identities of the exposed slots, the hidden block's product is
forced by (IX.0a), and a uniformly random completion hits that product with probability $\binom{K-s+k}k^{-1}\le
e^{-kJ}$. The union over $4^s$ (or $7^s$) slot labels and $(2Y)^3$ coefficient triples, $\log Y=O(W/b)$, is beaten
because $kJ\ge sJ/5$ grows faster than $s\log7+O(W/b)$ once $J\ge40$ and $b\gg\log P/J$.

*What fails for a fixed parent.* The probability space is a point. For fixed $q$, the candidates (label assignment,
coefficient triple) either solve (IX.0a) or not, and the union bound returns "#solutions ≤ #candidates". The question
becomes deterministic: does the three-term equation
$$aX+bY=cZ,\qquad X,Y,Z \text{ products of disjoint subsets of } S,\ |abc|\le Y^3,\ \omega(\text{one of }X,Y,Z)\ge s/5,$$
with the size constraints $|aX|,|bY|,|cZ|\le C\sqrt q/G$, have solutions? This is a squarefree three-term $S$-unit
equation with small coefficients and $|S|=s$ growing.

* abc: the radical of $aX\cdot bY\cdot cZ$ is at least the product of the primes of $XYZ$, which is at least each
  term up to $Y^3$; abc is vacuous (record §5 item 6, now for the twin as well).
* Upper counts: Evertse and Beukers-Schlickewei give $\exp(O(s))$ solutions per coefficient triple (rank of the
  $S$-unit group $\asymp s$); this is no better than the trivial $4^s$, and to exclude a cluster one needs "no
  solution", not a count (record §5 items 5, 14: "S-unit counts exponential in $s$; exceptional sets depend on $S$").
* Lower constructions (literature, **snippets only**, not read in full): Erdős-Stewart-Tijdeman (Compositio 1988) give
  sets $S$ of $s$ primes for which $x+y=z$ has $\exp((4-\epsilon)(s/\log s)^{1/2})$ coprime $S$-unit solutions;
  Konyagin-Soundararajan improve this to $\exp(s^{2-\sqrt2-\epsilon})$; for the first $s$ primes a lower bound
  $\exp(s^{1/8-\epsilon})$ is attributed to Lagarias-Soundararajan. These solutions are not squarefree with balanced
  sizes and a block of $s/5$ primes, so they are not twin clusters; but they show that no argument that ignores the
  squarefree-balanced structure can give "no solution" uniformly in $S$, and that dense initial segments are where
  $S$-unit equations are richest.
* Information view: on average, (IX.0a) carries $\ge(s/5)J$ nats of information about the parent, while its
  description costs $O(s)+O(W/b)$ nats. For a fixed parent the information is already spent (the adversary chose $q$),
  so the equation carries none. That is the precise sense in which the triangle gives a rate on average but not in the
  worst case. It also explains why the dense profile ($J=0$) is the hardest case.

---

## 3. Deterministic results

| result | status (Record: `ref9L_additive.md`, `ref9L_relax.md`) | new? | what it kills |
|---|---|---|---|
| IX.1: no affine 3-flat of twin signatures at $L<q^{2/3}/4$ | proved; exponent 2/3 sharp (cube family) | new (twin) | linear (Walsh, Sylvester-Hadamard) signature profiles with $\ge8$ marks; every integer transfer of the char-2 hyperplane mechanism with linear bits (Part V, Part VI binary families, attack referee's char-2 Hadamard family) |
| IX.3 twin: $\le(1+C^2)(1+2C^2)$ pairs per signature difference; $B\ll_C\sqrt D$ | proved | new (twin analogue of the distinct-label lemma) | nothing among fakes: P2 and Paley/skew-Hadamard fakes are Sidon, $D=1+\binom M2$ |
| IX.3 circle: $M\ll_C\sqrt D$ | proved | not new in substance (README_II distinct-label lemma) | same: saturated by fakes |
| IX.5: bits Boolean in $\le d$ affine forms on $n$ labels: $B\ll n^{1-2^{1-2d}}$, $M\ll n^{1-2^{1-2b_d}}$; $n^{7/8}$ at $d=2$; partial coverage $q_g\ge q^{1/2+\epsilon}$ | twin parts per `ref9L_additive`/`ref9L_nonlinear` (the latter not finished when read) | new to the record; literature not checked | low-complexity (degree-$\le d$ Reed-Muller-type) profiles with $n=O(s)$; **not** Legendre/P2 (parity singletons need degree $>d$, IX.6), not full-rank quadratic (Gold $x^3$) cores, not Paley |
| IX.6 modulus-alone: $B\ll1+W/\log P_*$ | proved | new for twin; old for circles (round 8 `stability.md` Cor 3.2, Thm 3.3) | parents with one dominant prime power; not balanced smooth parents |
| IX.6 density: $B_{\max}/\omega(q)\to0$ in natural density | proved, weak | new, trivial | typical integers, which are not smooth; irrelevant to the extremal regime |

Summary: the deterministic results give the first integer theorems that kill the characteristic-2 (linear/Walsh)
mechanisms *by an argument valid in characteristic 0 at critical length*; this is a concrete instance of the filter
F3+++ list "characteristic 0 used essentially" (here: Walsh inversion over $\mathbb Q$ and the divisibility
$q/Q_u\mid F_u$, which in characteristic 2 fails because the sign characters become 1, IX.6). They leave every
recorded integer-circle fake untouched. Their value for the rate is as filters: a proposed worst-case input that does
not beat Sidon-type signature conditions is useless.

---

## 4. IX.2 and the Question A box

**My check (Computed, `pc2.py`).** With seven finite marks and pairwise coprime $e_T$ ($3\le|T|\le6$), at most one $e_T$
is even. Enumerating residues modulo $2^{v+1}$ ($v=v_2(e_T)\in\{1,2\}$, all sizes $|T|=3,\dots,6$) and the case "no
even $e_T$", the minimum number of even pair cofactors $c_{ij}$ is 9 (no even block), 7 ($|T|=3$), 5 ($|T|=4,5$), 6
($|T|=6$). So at least 5 cofactors are even and pairwise coprimality fails. The argument depends only on lifts modulo
$2^{v+1}$, so it holds for every $v$. This agrees with Astra ("at least five") and `ref9L_fallback.md` (which also shows
$\mathrm{Bal}^{pc}=\emptyset$ for contact integers: exactly 5 even at $p=2$, at least 2 divisible by 3).

**Consequences if the fallback is vacuous (the expected outcome):**

1. Question A box: delete "or at least its lattice form $(L_8^{pc})$, which every proof must cover". Replace by the
   corrected fallback with an explicit small-prime exemption (wording in §7.2), and say plainly that the printed form
   had empty input.
2. §2 fact 3 ("every proof of $(L_8)$ must prove $(L_8^{pc})$"; theory Prop 2.1(4); census Prop 1.3(2); the abc
   criterion Prop 2.6's conclusion) are true but vacuous; restate them for the class "pairwise coprime outside
   $\{2,3\}$, with $\{2,3\}$-part bounded". The heuristic count $98-108=-10$ is about balanced profiles generally and
   survives, but must be attached to the corrected class.
3. The sandwich $(L_8^c)^{\rm all}\Rightarrow(L_8)\Rightarrow(L_8^f)\Rightarrow(L_8^c)^{\rm cat}\Rightarrow(L_8^{pc})$
   stays formally true with a vacuous right end; replace the right end by the corrected form.
4. Judge's open question 4 of `README_VIIIL.md` §5 ("Is $(L_8^{pc})$ Rudin-hard?", P≈0.25): answered **no** for the
   printed form (an unconditional theorem with empty input implies Rudin only if Rudin is provable outright). Re-ask
   it for the corrected form. The obstacle the judge named stays: VIII.2's square-block data have small cofactors
   sharing primes, so they are far from the corrected class; and the corrected class now has a mandatory
   $\{2,3\}$-part that the VIII.2 data must be made to respect.
5. The record's "must not contain" list gains: "$(L_8^{pc})$ is open" (for the printed form).

**Consequences if the fallback is not vacuous** (if another referee finds an intended reading with nonempty input):
the box and the judge's question stand as printed, but the intended reading (which primes are exempt, whether $e_T$
for $|T|=2$ are blocks) must be written into the box, since the literal reading is empty.

---

## 5. The twin as a sister problem

* **Gain (Record: proved, `ref9L_fallback.md`).** The uniform critical two-class bound is equivalent to
  Cilleruelo-Granville's Conjecture 9, and the little-oh versions at scales $\log q/\log\log q$ and $\omega(q)$ are
  equivalent too. So the record gains a sister problem that is (i) named in the literature, (ii) free of fakes (every
  admissible configuration is genuine), and (iii) has its own rate question: *do intervals of length $\sqrt q$
  contain $o(\log q/\log\log q)$ roots of $x^2\equiv a\pmod q$?* The fallback referee found no literature result
  improving the trivial order (snippets only). My own search (snippets only: Marklof-Welsh, Duke 2023, fine-scale
  limit laws for roots of quadratic congruences with squarefree modulus; Shkredov-Shparlinski-Zaharescu, Math. Z. 2024,
  short intervals for prime modulus) found distributional results, not worst-case counts. **Gap:** whether
  "almost all squarefree $q$ of a smooth pool have $O(\log\omega(q))$ square roots of every $a$ in every interval of
  length $\sqrt q$" (the twin part of IX.0) is in the literature was not determined.
* **Transfer.** None known in either direction, and I see none. The twin is a statement on the split torus, the circle
  problem on the anisotropic torus $x^2+y^2$ (the Liouville exponents coincide for every $k$, `README_VIIIL.md` §1.7, but
  that is a coincidence of local data). The classical bijection between primitive representations of $N$ and roots of
  $\nu^2\equiv-1\pmod N$ does not send short arcs to short intervals. What transfers is method: the triangle
  extraction, the signature tests, and the modulus criterion all run in both settings (IX.0, IX.3, IX.5, IX.6); the
  twin is therefore a good **laboratory** (a method that fails there fails on circles), but a result there is not a
  result on circles.
* Calibration: a twin *rate* is not known to imply anything about Rudin's conjecture. Cilleruelo-Granville's
  implication (Conjecture 9 implies Rudin) uses the uniform bound; no rate version of it is recorded, and I did not
  derive one.

---

## 6. Next steps toward the growth rate (ranked)

Probabilities are of moving the worst-case rate $o(\log R/\log\log R)$ (or $o(\omega(q))$ in the twin), unconditionally,
in a natural class, within one more exchange.

### 6.1 (a) Removing exceptional parents by many triangles / counting over exponent patterns — P ≈ 1%

Explicit count (Proved). Fix a parent $q$ of $s$ primes, $W=\log q$, and a cluster of $m$ marks.
* Sign patterns: $2^{(m-1)s}$ (modulo global flips). Even up to permuting primes the column multiset alone has
  $\log\binom{s+2^{m-1}}{s}$ nats, e.g. $9.6\cdot10^2$ at $s=10^3,m=10$ (Computed); but for a fixed parent the
  assignment of columns to the actual primes is what matters, and that costs $(m-1)s\log2$.
* Cofactors: the triangle identities make the table determined by a star; a good star has
  $\sum_i\log|t_{1i}|\le W/2+O_C(m)$, so $\exp(W/2+O(m))$ tables.
* Constraints: $\binom{m-1}2$ divisibilities $h_{ij}\mid h_{1j}t_{1j}-h_{1i}t_{1i}$, each of "probability" $1/h_{ij}\approx
  q^{-1/2}$. Heuristic expected count $2^{(m-1)s}e^{W/2}q^{-(m-1)(m-2)/4}$, negligible for $m\ge4$ (**Heuristic**; this
  is exactly the independence heuristic that the IX.4 relaxation shows no local input can certify).
* **Fixed parent, union bound over cofactors only:** the events are deterministic, so the bound is
  "#clusters ≤ $2^{(m-1)s}e^{W/2+O(m)}$", which excludes nothing.
* **Counting parents:** encoding a parent by (pattern, star table, constant part) costs $(m-1)s\log2+W/2$ nats, and
  $W/2\ge(s/2)\log P>sJ\ge\log\binom Ks-s$ because $K\le P$. So the encoding is *longer than the parent itself* for every
  $m\ge2$: an entropy count over patterns and cofactor tables cannot bound the number of parents admitting a cluster
  below the trivial $\binom Ks$. Astra's single triangle is cheaper precisely because it pays only $3\log Y=O(W/b)$ and
  $s\log4$, and it gains at most $sJ$. Additional triangles add cost at least $s\log2$ per mark and their total gain is
  capped by the parent's entropy $s(1+J)$; so at best they improve the exceptional exponent by a constant factor, and
  at $J=0$ they give nothing.
* Conclusion: (a) is not a route to the worst case. It could sharpen the measure-theoretic statement (e.g. $J$ down to
  about 7), which has low value. *First step if pursued anyway:* write the exposure with heavy blocks from $m$ marks and
  compute the optimal $J_0$.

### 6.2 (e) Worst-case twin for squarefree comparable-prime moduli through signed pair residues — P ≈ 4%

The relaxation and `ref9L_relax.md` show that the missing input in the twin is equivalent to the signed pair residues
mod $q$, a pair-level datum. IX.3 already exploits it for pairs with equal signature difference (sum and squared
difference congruences). *Question:* can the sum/square device be extended from "equal signature difference" to
"signature differences in a small set" or to triples, giving $B=o(s)$ for all squarefree $q$ with comparable primes?
*First step:* on the Paley cut pattern (where the cuts are unrealisable even up to length $q^{1-O(1/m)}$, IX.4 strengthened),
isolate which pair-residue congruences produce the contradiction, and test whether the same congruences constrain any
Sidon pattern. This is the only item that uses a property of genuine integers that the relaxation demonstrably lacks
**and** that is meaningful for a fixed parent.

### 6.3 (c) The parallel thread's projection theorem (unverified) — P ≈ 3%

If correct, $2+X_b\gg M^{3/2}$ for every Gaussian integer projection, and the exponent-row certificate gives $M=o(s)$
when some $c$ with $\sum c_i=2$ has $(1+\|c\|_1)^2T(c)=o(s^3)$. This is a deterministic, parent-specific test, the circle
counterpart of IX.5's restricted classes. It does not combine with the triangle in any visible way (projected
differences $\mathrm{Re}(\bar b\,2i\eta H_{ij}t_{ij})$ do not factor). *First step:* decide, once it is refereed, whether
typical-parent clusters of size $\gg\log s$ would automatically carry a short certificate; if yes, the projection gives
an alternative proof of the random-parent rate, not a worst-case one.

### 6.4 (b) Triangle identity plus Theorem 1.2's barycentre method — P ≈ 2%

The barycentre argument averages the pair inequality over two disjoint measures with equal barycentre; the triangle
is an exact identity on three marks. A combination would need an averaged *exact* identity, e.g. a weighted sum of
triangle equations over a measure on triples, which loses integrality. *First step:* check whether the affine
dependence produced in the proof of Theorem 1.2 can be chosen supported on a triangle-rich subset; I do not see why.

### 6.5 (d) Derandomisation: fix the parent, vary the cluster — P ≈ 2%

The cluster is already varied (the good triple is chosen by averaging). Varying the arc on a fixed circle averages
$\binom M3$ over arcs and cannot bound a single arc. A genuine derandomisation would need a family of auxiliary parents
built from one cluster (e.g. multiplying by Gaussian integers, or deleting primes constant on a sub-cluster) on which a
random-parent statement applies; the derived parents are not uniformly distributed in any pool. *First step:* none
worth spending effort on unless a construction of such derived parents with genuine entropy is found.

### 6.6 Bookkeeping and literature — P ≈ 0.9 (of a clear answer; no rate)

Record §7 below; confirm in full text whether the twin random-parent statement or any rate for roots of quadratic
congruences in intervals of length $\sqrt q$ is in the literature (Marklof-Welsh; Shkredov-Shparlinski-Zaharescu;
Cilleruelo-Garaev-Ostafe-Shparlinski on concentration of points), and read Erdős-Stewart-Tijdeman and
Konyagin-Soundararajan for the exact exponents quoted in §2.3.

---

## 7. Proposed wording

### 7.1 Expert note §1 (after the fake-circle paragraph, before "Actual clusters satisfy")

> *Random parents* (round 31, Parts III and IX; almost-all over a prime pool, **not** a worst-case statement). Draw $s$
> primes uniformly without replacement from a pool of $K\ge e^{40}s$ primes in $[P,2P]$ and let $q$ be their product,
> $J=\log(K/s)$. Outside a fraction $\exp(-sJ/20)$ of parents, every arc of length $Cq^{1/4}$ of $x^2+y^2=q$ (pool of
> split primes) and every interval of length $C\sqrt q$ with at most two prescribed classes modulo each $p\mid q$ carries
> $O(1+\log P/J)$ marks, uniformly in the arc, the interval and the classes. At $P=As\log s$ ($A\ge2e^{40}$) this is
> $O(\log s)$ at height $\log q\sim s\log s$, where the worst-case bound is of order $s$; at $P=s^{1+\epsilon}$ it is
> $O(1)$. Bounded exponents chosen after the sample are allowed, with constants degrading in $E$. The exceptional
> fraction is $R^{-o(1)}$ at the smooth scale and the exceptional set may contain $(\#\text{parents})^{0.95}$ parents.
> Nothing is said about any individual $N$, about the dense profile extremal for Theorem 1.1 ($K=s$), or about the
> growth rate for general $R$. On typical parents of narrow bins the P2 fakes exist with $\asymp s$ marks while
> genuine arcs carry $O(1)$ points (Part IX, same-norm comparison): the inputs of the fake list cannot certify $M=o(s)$
> even where the truth is $O(1)$.

### 7.2 Expert note §2

* Replace the Question A box's first sentence by: "**Question A.** Prove $(L_8)$, or at least its lattice form with
  the blocks $e_T$ and the pair cofactors pairwise coprime **outside the primes 2 and 3**, with $\{2,3\}$-parts bounded
  by an absolute constant, which every proof of $(L_8)$ must cover." (Exact replacement text for the lattice box: as in
  `ref9L_fallback.md` §3.)
* After the lattice box add: "*Correction (round 31, Part IX.2; refereed).* As printed above (all $e_T$ and all pair
  cofactors pairwise coprime), the lattice form has empty input: among seven marks at least five pair cofactors are
  even. The same holds for pairwise coprime contact integers ($\mathrm{Bal}^{pc}=\emptyset$: exactly five visible
  contact integers are even at every point). Statements restricted to $\mathrm{Bal}^{pc}$ are vacuous and are restated
  with the $\{2,3\}$ exemption."
* Fact 3 of "Round-18 facts": insert "outside $\{2,3\}$" after "pairwise-coprime".

### 7.3 Expert note §5

* New item (or extension of item 1): "**The exact additive identity is the input the fakes lack** (round 31, Part IX;
  refereed). Genuine differences satisfy $t_{12}A_{12}+t_{23}A_{23}=t_{13}A_{13}$ with integer cofactors and the blocks
  products of the parent's primes (twin), resp. $z_j-z_i=2i\eta H_{ij}t_{ij}$ with $t_{ij}\in\mathbb Z$ and the cocycle in
  $\mathbb Z[i]$ (circle). A relaxation with every pair magnitude, gcd, valuation partition, cut signature and Sidon
  signatures, but without these identities, reaches the larger-sieve constant $\frac12\log q/\log\log q$ in the twin
  (IX.4); in the twin the identities are equivalent to the signed pair residues modulo $q$. The identity gives a rate on
  random parents (§1) because it determines a block of $\ge s/5$ parent primes; for a fixed parent it is a squarefree
  three-term $S$-unit equation with $|S|=s$, where abc is vacuous and known counts are exponential in $s$. Its
  consequences at primes not dividing $N$ (transitivity, round 9 Prop 5.1(a)) and its triangle-area and abc
  consequences remain capped; **the identity itself is not capped**, and its fixed-parent content is open."
* Item 6: append "abc on the integer triangle of the two-class twin is vacuous for the same reason (squarefree terms)."
* Overview §5.4 row 16 / round 10 §6: qualify "$\sigma$-cocycle ... capped ($A_0$)" to "its local consequences and
  triangle areas are capped; the exact identity with integer cofactors is not (Part IX)".

### 7.4 Expert note §5 item 12 (new bullet)

> * ***Part IX: the critical two-class problem*** (`round31/README_IXL.md`, to be written; refereed blocks ref9L_*).
>   The twin of the arc problem (at most two classes modulo each $p^f\Vert q$, intervals of length $C\sqrt q$) is
>   equivalent to Cilleruelo–Granville's Conjecture 9 at the uniform and little-oh scales. Proved: the random-parent
>   rate (§1) for the twin and for circles, including the endpoint $P=As\log s$ ($O(\log s)$ marks) and bounded
>   exponents; the same-norm comparison (P2 fakes with $\asymp s$ marks on norms where genuine arcs have $O(1)$);
>   a squarefree relaxation reaching the larger-sieve constant without the triangle identities; no affine 3-flat of
>   signatures in the twin at $L<q^{2/3}/4$, $B\ll\sqrt{|A+A|}$, and $n^{1-2^{1-2d}}$ bounds for signatures of
>   complexity $d$ (these kill linear/Walsh and characteristic-2-type profiles over $\mathbb Z$, and no recorded fake);
>   $B\ll1+\log q/\log P_*$. The printed lattice form $(L_8^{pc})$ has empty input (§2). **The growth rate for general
>   $R$, the uniform bound and Questions A and B are unchanged.** Worst-case parents, even squarefree with comparable
>   primes, are open in both problems.

Add `README_IXL.md` and the judge's score to the list.

### 7.5 Sentences the note must not contain

* "The uniform bound holds for almost all circles" / "for almost all smooth parents" (without the pool measure, and
  at the smooth scale the bound is $O(\log s)$, not $O(1)$).
* "Part IX improves the growth rate", "the rate is $O(\log\log R)$" (without "on random parents of a pool").
* "The exceptional parents are polynomially rare" (at $P=As\log s$ they are $R^{-o(1)}$-rare and may number
  $(\#\text{parents})^{0.95}$).
* "Part IX covers the extremal profile of Theorem 1.1".
* "The triangle identity is capped" and, conversely, "the cocycle is capped (A_0)" without the qualification of §7.3.
* "The relaxation isolates the triangle identity" without "equivalently, the signed pair residues modulo $q$".
* "The twin rate implies the circle rate" or the converse; "the uniform circle bound implies Rudin".
* "$(L_8^{pc})$ is open" (printed form); "$(L_8^{pc})$ is Rudin-hard" (either form).
* "Fakes fail the same-norm test" (the fakes satisfy every $A_{2M}$ input; the separation is by the truth, not by a
  test).
* "Typical dyadic parents separate genuine from fake" (shown only for typical bin-parents).

### 7.6 The question to put to the colleague

> Part IX is refereed and nothing is false (see the block reports). Its content for the rate is an average over
> parents, and the analysis shows why: the triangle identity determines a block of $\ge s/5$ primes, which is
> information only if the parent carries entropy; for a fixed parent it is a squarefree three-term $S$-unit equation,
> and counting over exponent patterns and cofactor tables is longer than the parent itself. The one input that is
> both missing from your relaxation and meaningful for a single parent is the **signed pair residue modulo $q$**
> (for $|d_{ij}|<q/3$ the triangle identities are equivalent to $d_{ij}\equiv\rho_j-\rho_i\pmod q$), which your IX.3 already
> uses for pairs with equal signature difference. **Question.** In the critical two-class problem with $q$
> squarefree and all primes in $[P,2P]$ (any $P\le s^{A}$), prove $B=o(s)$ for **every** such $q$, or exhibit, for one
> fixed dense $q$ (e.g. all primes in $[P,2P]$), a family of configurations with $B\ge cs$ at length $C\sqrt q$, or a
> relaxation that keeps all signed pair residues modulo $q$ and still reaches $B\gg s$. Say which of these three you
> can do; full proofs or exact computations; test against the IX.4 Paley cuts.

### 7.7 Should the thread continue?

Frankly: Part IX is the most informative pass of this thread on the rate. It gives the first rate $o(\omega)$ at the
smooth critical height for a natural measure, and the cleanest identification so far of what the fakes lack. But
every rate in it is either measure-theoretic or for restricted signature classes that contain no recorded fake, and
the remaining step (a fixed parent) is a uniform-in-$S$ $S$-unit problem that the record (§5 items 5, 14) and the
literature show to be hard. I recommend **one** more exchange with the question in §7.6, and closing the thread if it
returns neither a worst-case theorem for the comparable-prime squarefree twin nor an obstruction there. Score on the
goal: 3/10 (sound, new, diagnostic; worst-case rate unchanged).

---

## Sources (web; snippets only, nothing read in full)

* Erdős, Stewart, Tijdeman, Compositio Math. 66 (1988): [uwaterloo PDF](https://uwaterloo.ca/pure-mathematics/sites/default/files/uploads/documents/cm_1988_66_1_37_0_0.pdf), [Rényi archive](https://www.renyi.hu/~p_erdos/1988-31.pdf)
* Lagarias, Soundararajan, smooth solutions to $a+b=c$: [arXiv 0911.4147](https://arxiv.org/pdf/0911.4147), [arXiv 1102.4911](https://arxiv.org/pdf/1102.4911)
* Many S-unit solutions on average: [arXiv 1108.3819](https://arxiv.org/pdf/1108.3819)
* Marklof, Welsh, fine-scale distribution of roots of quadratic congruences: [Duke Math. J. 172 (2023)](https://projecteuclid.org/journals/duke-mathematical-journal/volume-172/issue-12/Fine-scale-distribution-of-roots-of-quadratic-congruences/10.1215/00127094-2022-0081.full)
* Shkredov, Shparlinski, Zaharescu, Math. Z. (2024): [Springer](https://link.springer.com/article/10.1007/s00209-024-03436-5)
