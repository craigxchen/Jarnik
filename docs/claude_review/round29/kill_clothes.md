# Round 29: adversarial review of lane "clothes"

Reviewed: `round29/clothes.md`, plus its notes and scripts. Independent scripts are in `round29/kill_clothes_checks/`:
* `kc.py`: Gaussian point enumeration via Gaussian factorisation, arc maxima, twin maxima, residue classes;
* `runtwin.py`: the twin tables;
* `smoothcls.py`: smooth-modulus classes;
* `check44.py`: the refutation of Prop 4.4.

Labels are as in the lane. **exact** means a finite exact computation. **evidence** means data that proves nothing.

## 0. Verdict

**Status: obstruction, survives in part. One headline item is false.**

* **Killed: Proposition 4.4 ("smooth moduli defeat the engine").** It is false, and with it the stated hypothesis
  of Theorem 4.3(ii) and the lane's *Lead*.
  * Inside one class mod $q$, with $\ell^t\Vert q$, the residues mod $\ell^{t+b}$ lie in exactly $\ell^b$
    classes.
  * So the depth-one pigeonhole that Theorem 7.5 of `two_thirds/proof.md` actually uses survives at every
    $\ell\mid q$, with no loss of constant.
  * Consequence: **the twin obeys Theorem 7.5 with $C=2/\lambda$ for every odd $q$ coprime to $N$.** There is no
    "smooth-conductor regime where the pigeonhole gives nothing".
* **Survives (re-derived, correct):**
  * the dictionary Theorem 4.1 (i)–(iv);
  * the form of discriminant $-4q^2$ (with a converse the lane did not prove, §1);
  * Prop 4.2 (twin Liouville gap);
  * Theorem 4.3(i) (twin cut inequality);
  * Props 5.1 and 5.2, Cor 5.3, and the negative decision on peeling;
  * forms (a) and (e).
* **Minor errors:**
  * an arithmetic slip in the Cilleruelo–Córdoba exponent;
  * a misdescription of the modular-concentration range in (b);
  * an unnecessary factor $2$;
  * the evidence table is thinner than stated (§3–4).
* **Novelty:** repackaged known barrier. `round4/hardness.md` already names "the Lenstra problem at exactly
  $\alpha=1/4$" as the non-archimedean twin. Peeling is the record's Thm 3.10, as the lane itself says. What is
  new is only the explicit single-form / ring-class dictionary, and (after this correction) the fact that the
  twin carries the full $\frac23$ theorem. Neither changes the growth rate.

## 1. What was re-derived and holds

* **Form (a).**
  * $m=|z-z_1|^2/2=N-\mathrm{Re}(z\bar z_1)$, and $v^2=m(2N-m)$.
  * $\gcd(m,2N-m)=\gcd(m,2N)=g$, so $m=ga^2$ and $2N-m=gb^2$.
  * Lattice condition $\bar z_1\mid -m+iv$: with $\iota=x_1/y_1$ ($\iota^2\equiv-1$) this is
    $v\equiv m/\iota=-\iota m\pmod N$. This matches the lane's $\rho=-x_1/y_1$. Correct.
  * F1 failure confirmed: the only angle information is the choice of root.
* **Theorem 4.1.**
  * (ii)⇒(iii) uses $\gcd(q,2N-m)=1$ to get $q^2\mid m$, and $\gcd(g,q)=1$; then $q^2\mid a^2\Rightarrow q\mid a$.
    Correct.
  * $L=\mathbb Zz_1+q\mathbb Z[i]$ has index $q$, because $z_1$ has additive order $q$ mod $q$ (a prime dividing
    $q$ and both coordinates would divide $N$).
  * The Gram discriminant is $-4q^2$. The form is primitive, since it represents $N$ and $q^2$ and
    $\gcd(N,q)=1$; the lane did not state this.
  * The class factor is $\#\{c:c^2\equiv1\ (q)\}=2^{\omega(q)}$ for odd $q$, so the lane's $2^{\omega(q)+1}$ is
    harmless but unnecessary.
  * *Converse (missing from the lane, needed for "equivalent to $r_f$ for single forms").* Let $f$ be any
    primitive form of discriminant $-4q^2$ and let $\gcd(n,q)=1$.
    * Realise $f$ as $|\cdot|^2$ on a proper $\mathcal O_q$-lattice $L_f\subset\mathbb Q(i)$ of covolume $q$.
    * Since $h(\mathbb Z[i])=1$, $L_f\mathbb Z[i]=\gamma\mathbb Z[i]$, and $[L_f\mathbb Z[i]:L_f]=q$ forces
      $|\gamma|=1$.
    * So $L=\gamma^{-1}L_f\subset\mathbb Z[i]$ has index $q$ and $L\supseteq q\mathbb Z[i]$.
    * If $z\in L$ has $|z|^2=n$ coprime to $q$, then $\mathbb Zz+q\mathbb Z[i]\subseteq L$ has the same index,
      so equality holds.

    Hence $r_f(n)\le2^{\omega(q)}T(n,q)$. With this the equivalence claimed in §4 is proved for $\gcd(n,q)=1$.
* **Prop 4.2.** Correct. The case $v=0$, $\varepsilon\ne1$ is vacuous, since $q\ge3$ cannot divide
  $\varepsilon-1$.
* **Theorem 4.3(i).**
  * $q\mid c_{ij}$ because $\gcd(q,g_{ij})=1$.
  * Identity (2.3) is algebraic, and $|z_i-z_j|\le2R'$.
  * $K\log q\ge K\log\lambda+\frac{M(M-1)}4\log R'$ leaves exactly $\frac M4\log R'+K\log(2/\lambda)$. Correct.
  * **exact**: 1396 classes (sizes up to 8) on 16 norms at $\lambda\in\{1,\frac12,\frac14\}$ with smooth $q$,
    normalised by the class gcd and with $Q$ computed from the $\pi_p$-valuations. 0 failures (`check44.py`).
* **Props 5.1, 5.2, Cor 5.3.**
  * 5.1: the algebra $p^{e(r+1)}\ge(C/\sqrt2)^{2(2r+1)}R'$ ⇔ the stated form with $k=2r$ was re-derived. The
    Cilleruelo–Córdoba exponent $r/(2r+1)=\frac12-\frac1{4r+2}$ is correct here.
  * 5.2: Plotkin, with $\langle h_x,h_y\rangle_w=W_1-2d_w\le-(e\log p-4\log C)$. Correct.
  * These are the record's Thm 3.1/Cor 3.2 of `round8/stability.md`, $M<2\log N/(\log P_+-4\log^+C)$,
    confirmed in the file. The lane's "rediscovered" label is accurate.
* **Theorem 5.4.** Correct as a statement about inputs: the inequalities 5.1–5.3 are pair data, and every
  recorded fake has $\log P_+=O(\log M)$ (`round8/stability.md` Cor 3.6).
* **Form (e).** The Dirichlet strip count $h|v|+L/Q+1$, 2 points per line, and $O(C^{3/2}N^{1/8})$ are correct.
  The bound is weaker than Cilleruelo–Córdoba.

## 2. Kill: Proposition 4.4 is false, and Theorem 4.3(ii) holds for all odd $q$ coprime to $N$

**Lemma K (confinement; proved).** Let $\ell$ be an odd prime with $\ell\nmid N$, $\ell^t\Vert q$ with $t\ge1$, and
$b\ge0$. For any $z_1$ with $|z_1|^2=N$,
$$\#\{z\bmod\ell^{t+b}:\ z\equiv z_1\ (\ell^t),\ |z|^2\equiv N\ (\ell^{t+b})\}=\ell^b.$$

*Proof.* By the lifting step of `two_thirds/proof.md` Lemma 3.2, each solution of $|z|^2\equiv N$ mod $\ell^a$
has exactly $\ell$ lifts mod $\ell^{a+1}$, because the gradient $(2x,2y)\not\equiv0$. Start from the single class
$z_1$ mod $\ell^t$ and apply this $b$ times. $\square$
(**exact**: 4528 base classes, $\ell\le13$, $t,b\le2$; 0 failures.)

**Corollary (proved).** Let $z_1,\dots,z_M$ be one twin class mod $q$ (odd, $\gcd(q,N)=1$), and $\ell^t\Vert q$.
* For $i\ne j$, $v_\ell(c_{ij}/q)=v_\ell(c_{ij})-t=\sum_{b\ge1}\mathbf 1[z_i\equiv z_j\ (\ell^{t+b})]$. This uses
  $\gcd(\ell,g_{ij})=1$, as in Lemma 3.3 (I)/(S∤).
* By Lemma K and (3.5),
  $$\sum_{i<j}v_\ell(c_{ij}/q)\ \ge\ \sum_{b\ge1}E(M,\ell^b)\ \ge\ E(M,\ell)\ \ge\ \frac{M^2}{2\ell}-\frac M2 .$$
* Section 6 of `two_thirds/proof.md` (the only part Theorem 7.5 uses) needs, per prime $\le M$, only:
  * inert $\ell$: $\frac{M^2}{2(\ell+1)}-\frac M2$. This holds, since $\frac1{2\ell}>\frac1{2(\ell+1)}$.
  * split $p\nmid N$: $\frac{c_pM^2}4-\frac M2$ (Lemma 4.3). This holds, since $c_p\le\lambda_p=\frac1{p-1}\le\frac2p$.
  * ramified: $(1+i)\mid c_{ij}/q$, since $q$ is odd.
  * $p\mid N$: unchanged, since $p\nmid q$.

  So $D(M)\le\sum\log|c_{ij}/q|+\frac Q2$, and Theorem 4.3(i) gives (6.2) with $C=2/\lambda$.
* **Hence Theorem 7.5 and Cor 8.1 hold for the twin, $3M\log M\le W+(14+4\log^+(2/\lambda))M$, for every odd $q$
  coprime to $N$.** The hypothesis "no prime factor $\le M$" is unnecessary.

**Where Prop 4.4 went wrong.** It counts only the collision levels $a\le t$. Those are indeed absorbed into $q$.
It misses that confinement to one class mod $\ell^t$ cuts the number of residues available at depth $t+b$ from
$(\ell\mp1)\ell^{t+b-1}$ to $\ell^b$. Equivalently, apply the engine to $\delta_j=(z_j-z_1)/q$ rather than to $z_j$.
- **exact** (`check44.py`): in all 1396 classes and at every $\ell\mid q$, the residues mod $\ell^{t+1}$ occupy
  $\le\ell$ classes, and $\sum_{i<j}v_\ell(c_{ij}/q)\ge E(M,\ell)$. 0 failures.

**Consequences.**
* The lane's *Lead* question (b), "can smooth-conductor ring class groups carry large fibres, where the Theorem 1.1
  engine is void?", has a false premise. The engine is not void there.
* The claim that the twin is "harder by known methods on smooth $q$" is withdrawn. On every odd $q$ the twin has
  exactly the arc's proved bound.
* What remains open is twin boundedness, or growth below $\frac23\log R/\log\log R$, which is open for the same
  reasons as for arcs.

## 3. Smaller errors

1. **§2(b): Cilleruelo–Córdoba exponent.** At $k=4$ the exponent is $\frac12-\frac1{4\cdot2+2}=\frac25$, not
   $\frac12-\frac18$. The conclusion "$\le4$ points on arcs $R^{3/8}$" still holds, since $\frac38<\frac25$.
2. **§2(b): the concentration range is misdescribed.**
   * The identity "cluster = box solutions mod $p$" holds for **every** prime $p>6LR$, not only for
     $p\asymp N^{3/4}$.
   * Taking $p>H^4$ puts any box in the range $H<p^{1/4}$. So $H^{o(1)}$ is available for arcs of every length.
     It is simply no stronger than the divisor bound $R^{o(1)}$.
   * The correct reason (b) gives nothing is that $H^{o(1)}$ is divisor-bound strength, not that it is confined to
     $R^{3/8}$.
   * Also, the circle analogue is not uncharted. Zumalacárregui, arXiv 1012.3569, "Concentration of points on
     modular quadratic forms", extends Cilleruelo–Garaev to $Q(x,y)\equiv\lambda$ (literature, abstract-level
     snippet; exact thresholds not verified).
3. **Theorem 4.1:** the factor $2^{\omega(q)+1}$ can be $2^{\omega(q)}$.
4. **§3 table:**
   * Blomer–Granville (Duke 135 (2006) 261–302) is stated for fundamental discriminants (snippet). $-4q^2$ is
     non-fundamental, and the extension is due to Toma (snippet).
   * Cilleruelo–Garaev's GAFA venue was not confirmed by search; arXiv 1007.1526 is.
5. **Evidence description (§4).**
   * Of the 24 distinct norms in the two lists, 12 have $r(N)>1100$, and the full twin search was skipped there.
   * The twin maxima come from small norms. The value 8 at $\lambda=\frac12$ is $N=1105$, $q=3$: the trivial
     pigeonhole $32$ points / $4$ residue classes mod $3$.
   * "Smooth-modulus classes $\le4$" holds at $\lambda=1$ only. At $\lambda=\frac12$ the smooth modulus $q=3$ gives 8.

## 4. Numerical claims recomputed (independent code)

| claim | recomputed | status |
|---|---|---|
| twin max at $\lambda=1$ is 5 | 5 ($N=2576450045$, $q=269$) | agrees |
| arc max at $C=2$ is 3 | 3 | agrees |
| twin max at $\lambda=\frac12$ is 8 | 8 ($N=1105$, $q=3$; $N=32045$, $q=7$) | agrees; trivial small-$q$ pigeonhole |
| arc max at $C=4$ is 5 | 5 | agrees |
| smooth-$q$ classes $\le4$ | $\le4$ at $\lambda=1$; 8 at $\lambda=\frac12$ | holds only at $\lambda=1$ |
| dictionary identities | re-proved by hand (§1) | proof suffices |

Outputs: `kill_clothes_checks/twin_1.out`, `twin_0.5.out`.

## 5. Filters applied independently

* **F1.**
  * (a), (b), (e), (f): the filter fails, as the lane says. All four are pair data plus the angle/root choice,
    satisfied by fakes.
  * Twin: after §2, the twin's proved inequalities are the arc's (5.1) with $C=2/\lambda$, and they are
    angle-blind. Group fakes (arbitrary classes in a group of order $\asymp q$) satisfy them.
    * The existence of group fakes at fixed $\lambda$ is heuristic only, as the lane says.
    * Note the extra freedom: $q$ ranges over $\asymp N^{1/4}$ values. The expected-count heuristic gives
      twin($M+1$) $\approx$ arc($M$), which is consistent with the data offset of about 2.
    * So the twin is not a harder problem in any proved sense, and no reduction in either direction exists.
* **F2.**
  * Lemma K1 transfers, as the lane says.
  * The Lenstra/CHN constant is $2+2\pi\alpha/(\alpha-\frac14)^{3/2}+4\alpha/(\alpha-\frac14)$, stated only for
    $\alpha>\frac14$ (Hales, arXiv 2410.05030, snippet). It diverges at $\frac14$.
  * No source found treats $\alpha=\frac14$ exactly. The twin is exactly critical. Confirmed.

## 6. Literature checks (WebSearch, snippets only)

| item | outcome |
|---|---|
| Cilleruelo–Córdoba (PAMS 115 (1992)) | statement: arc $\sqrt2R^{1/2-1/(4\lfloor k/2\rfloor+2)}$ has $\le k$ points. Confirmed via arXiv 1603.09646 |
| CHN, Math. Comp. 77 (2008) | $O((\alpha-\frac14)^{-3/2})$ for $\alpha>\frac14$. Confirmed; Lenstra's original bound was $O((\alpha-\frac14)^{-2})$ |
| Cilleruelo–Garaev, arXiv 1007.1526 | $M<p^{1/4}\Rightarrow M^{o(1)}$ for $xy\equiv\lambda$. Confirmed; venue not confirmed |
| Blomer–Granville, Duke 135 (2006) | averaged $r_f(n)^\beta$, uniform in $D=o(x)$, fundamental $D$. Confirmed |
| Ellenberg–Venkatesh | Lemma 2.3 concerns small primes and $\ell$-torsion. The precise "$<\sqrt{|D|}/2$ ⇒ non-principal" form was not confirmed (memory; elementary in any case, and Prop 4.2 is self-contained) |

The uploaded PDF concerns the spectral lane. This lane makes no claim about it, so there was nothing to check.

## 7. Novelty and corrected statement

**Novelty: repackaged known barrier.**
* The twin as Lenstra at $\alpha=\frac14$ is in `round4/hardness.md` §6–7.
* Peeling is `round8/stability.md` Thm 3.1.
* New but modest: the explicit dictionary class ↔ $r_f(N)$ for one form of discriminant $-4q^2$, with the converse
  of §1; the twin Liouville gap; and the transfer of the full $\frac23$ theorem to every odd modulus (corrected here).

**Corrected statement.**
* *(Twin dictionary.)* For odd $q$ with $\gcd(q,N)=1$, a class $\{z:|z|^2=N,\ z\equiv z_0\ (q)\}$ is a coset of
  $\ker(a\mapsto\prod u_p^{a_p})$ in the exponent box. Moreover
  $T(N,q)\le\max_f r_f(N)\le2^{\omega(q)}T(N,q)$, over primitive forms $f$ of discriminant $-4q^2$.
* *(Twin bound.)* For every such $q\ge\lambda N^{1/4}$, a twin class of size $M$ satisfies
  $3M\log M\le\log N+(14+4\log^+(2/\lambda))M$. So $M\le(\frac23+\varepsilon)\log R/\log\log R$, with no
  condition on the prime factors of $q$.
* *(Open.)* Boundedness of $T(N,\lambda N^{1/4})$, equivalently of $r_f(n)$ for $n\asymp D^2$, is open. It is
  equally blocked by F1 (angle-blind inequalities) and F2 (critical exponent). No implication to or from the arc
  problem is known.
* *(Peeling.)* Peeling gives exactly $M<2W/(\log P_+-4\log^+C)+O(1)$, the record's Thm 3.10. Every fake is
  $M^{O(1)}$-smooth, so the dichotomy route is closed.
* *(Withdrawn.)* Prop 4.4 and the smooth-conductor lead.
