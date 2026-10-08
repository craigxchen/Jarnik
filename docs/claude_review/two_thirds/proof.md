# A $(2/3+\varepsilon)\log R/\log\log R$ bound for lattice points on short arcs

> **Review status (2026-10-08).** Two independent adversarial referees re-derived every step and
> could not refute the result; neither found a mathematical error. Their minor points: (i) the
> strict inequality in Lemma 7.1 step 3 should be non-strict at `y=1` (the final bound is unaffected);
> (ii) the floating-point rounding claim for `sigma, a_3, a_1` is asserted, but a referee recomputed
> all three with correctly rounded 34-digit logarithms over every prime up to `10^7` and confirmed
> the stated bounds; (iii) the tail bound for `a_1` uses an extra factor `1/(P-1)` that `constants.py`
> applies but the text does not state; (iv) the remark that 2/3 is the limit of these ingredients
> needs the easy upper bound `F*(M) <= (3/8)M^2 log M + O(M^2)` (take `e=1`, `n_0=n_1=M/2`);
> (v) "only `l < M` contribute" should read `l <= M-2`; (vi) the `F*(48)` table entry rounds to
> 1686.6. **This improves only the leading constant; the growth rate is still `log R/loglog R`.**
> Checkers: `check_chain_lemma.py`, `check_identity.py`, `constants.py`, `finite_inequality.py`
> (run from this directory; `gauss.py` is a shared helper).


**Status.** This note gives a complete, self-contained proof of the target theorem
(Corollary 8.1 below):

> For every $C>0$ and $\varepsilon>0$ there is an explicit $R_0(C,\varepsilon)$ such that for all
> $R\ge R_0$ with $R^2\in\mathbb N$, every arc of the circle $x^2+y^2=R^2$ of length at most
> $C\sqrt R$ contains at most $(2/3+\varepsilon)\log R/\log\log R$ lattice points.

It follows from an effective inequality (Theorem 7.5): every $M\ge 1$ distinct lattice points on such
an arc satisfy

$$3M\log M\;\le\;W+K_C\,M,\qquad W=\log(R^2),\qquad K_C=14+4\log^+C ,$$

which in turn comes from an exact finite inequality (Theorems 5.1 and 5.2). That inequality extends
inequality (6) of `inert_prime_cofactor_bound.md` to split primes, at every prime-power depth.
The route is the one in `two_thirds_sketch.md`. It uses:

- the gcd normalisation;
- the pair chord identity;
- collision gains at inert primes, at split primes not dividing $N$, and at split primes dividing $N$
  for pairs on the same exponent level;
- the cut identity;
- a chain lemma (Lemma 4.1), which couples the cut slack at a split prime with its level collisions;
- one-sided Mertens estimates in the classes mod 4 (proved here from scratch, with explicit constants).

Every step was re-derived and checked. Section 9 lists the specific audit points. **The sketch is
correct.** No step needed repair, beyond making all constants explicit.

This does **not** prove the uniform (bounded) conjecture. It improves the leading constant of the
growth rate from $1$ to $2/3$, and the growth rate itself is still $\log R/\log\log R$ (Section 10).

What is proved, what is quoted, and what is numerical:

* **Proved here, in full:** all algebra and arithmetic (Sections 1–6), the chain lemma (Section 4),
  the prime-sum bounds (Lemmas 7.1–7.3), the effective inequality (Theorem 7.5), and the asymptotic
  corollary (Section 8).
* **Quoted classical facts, with short proofs or standard references:**
  - unique factorisation in $\mathbb Z[i]$ and the classification of Gaussian primes;
  - Leibniz's series $\sum_{k\ge0}(-1)^k/(2k+1)=\pi/4$;
  - the elementary alternating-series bound.
* **Numerical (rigorous finite computation with explicit tail bounds):** only the three constants
  $\sigma,a_3,a_1$ of Lemma 7.4, which enter $K_C$. All other numerics in Section 11 are
  *evidence and sanity checks only*. No proof depends on them.

---

## 0. Notation and conventions

* $N\ge1$ is an integer and $R=\sqrt N$. A *configuration* is a set of $M\ge1$ distinct Gaussian
  integers $z_1,\dots,z_M$ with $|z_j|^2=N$, all lying on one arc $\Gamma$ of the circle $|z|=R$,
  of length $L\le C\sqrt R$. Then $|z_i-z_j|\le L$, because a chord is never longer than the arc
  joining its endpoints inside $\Gamma$. This chord bound is the **only** geometric input, so
  nothing depends on the position of the arc.
* $W:=\log N=2\log R$, $\;K:=\binom M2$, and $\log^+C:=\max(0,\log C)$.
* $p$ always denotes a rational prime $\equiv1\pmod 4$ ("split"), $\ell$ a rational prime
  $\equiv 3\pmod 4$ ("inert"), and $q$ an arbitrary rational prime. For each split $p$ we fix one
  Gaussian prime $\pi=\pi_p$ with $\pi\bar\pi=p$. Then $\pi$ and $\bar\pi$ are non-associate primes
  of $\mathbb Z[i]$. The prime $1+i$ is the unique prime above $2$, and inert $\ell$ stay prime in
  $\mathbb Z[i]$.
* "$w\mid z$" means divisibility in $\mathbb Z[i]$, and $\gcd$ is a Gaussian gcd, defined up to a
  unit. Every modulus, norm and divisibility statement below is invariant under changing associates.
* $\chi=\chi_4$ is the non-trivial character mod 4. $\theta(y)=\sum_{q\le y}\log q$ and
  $\psi(y)=\sum_{n\le y}\Lambda(n)$.
* For integers $n\ge0$ and $Q\ge1$,
  $$E(n,Q):=\min\Big\{\sum_{r=1}^{Q}\binom{m_r}{2}:\ m_r\in\mathbb Z_{\ge0},\ \sum_r m_r=n\Big\}$$
  is the least number of unordered pairs in a common class when $n$ objects are put into at most
  $Q$ classes (Lemma 3.4).

---

## 1. Gcd normalisation

**Lemma 1.1.** Take a configuration (radius $R$, arc length $L\le C\sqrt R$), let
$g=\gcd(z_1,\dots,z_M)$, and put $z_j'=z_j/g$.

1. The $z'_j$ are distinct Gaussian integers with $|z'_j|^2=N':=N/|g|^2\in\mathbb Z$, so $N'\le N$.
2. They lie on an arc of length $L'=L/|g|$ of the circle of radius $R'=R/|g|$, and
   $L'\le C\sqrt{R'}$.
3. $\gcd(z'_1,\dots,z'_M)=1$. Moreover:
   * (a) $N'$ is odd;
   * (b) no prime $\ell\equiv3\ (4)$ divides $N'$;
   * (c) if $N'=\prod_p p^{e_p}$, then every $z'_j=\varepsilon_j\prod_{p\mid N'}\pi_p^{a_{jp}}\bar\pi_p^{\,e_p-a_{jp}}$ with $\varepsilon_j\in\{\pm1,\pm i\}$ and $0\le a_{jp}\le e_p$;
   * (d) for each $p\mid N'$, $\min_j a_{jp}=0$ and $\max_j a_{jp}=e_p$.

*Proof.*

1. Multiplication by $1/g$ is a similarity of $\mathbb C$ with ratio $1/|g|$. It is injective and maps
   the circle and arc onto a circle of radius $R/|g|$ and an arc of length $L/|g|$. Since $|g|\ge1$,
   $$L'=\frac{L}{|g|}\le \frac{C\sqrt R}{|g|}=\frac{C\sqrt{R'}}{\sqrt{|g|}}\le C\sqrt{R'}.$$
   This proves 1 and 2.
2. *Gcd.* $\gcd(z_j/g)=1$ by the definition of $g$.
3. *(a)* If $2\mid N'=z'_j\bar z'_j$, then $1+i$ divides $z'_j$ or $\bar z'_j$. Since
   $\overline{1+i}=-i(1+i)$, in either case $(1+i)\mid z'_j$. This holds for every $j$, which
   contradicts $\gcd=1$.
4. *(b)* An inert $\ell$ is a Gaussian prime. From $\ell\mid z'_j\bar z'_j$ we get
   $\ell\mid z'_j$ or $\ell\mid \bar z'_j$, and since $\bar\ell=\ell$, in both cases
   $\ell\mid z'_j$. Again this holds for every $j$, a contradiction.
5. *(c)* By (a) and (b), every Gaussian prime dividing $z'_j$ is some $\pi_p$ or $\bar\pi_p$ with
   $p\mid N'$. Put $a=v_{\pi_p}(z'_j)$ and $b=v_{\bar\pi_p}(z'_j)$. Then
   $v_{\pi_p}(z'_j\bar z'_j)=a+b$, while $v_{\pi_p}(N')=e_p$. Unique factorisation gives (c).
6. *(d)* If $\min_j a_{jp}\ge1$, then $\pi_p$ divides every $z'_j$. If $\max_j a_{jp}\le e_p-1$,
   then $\bar\pi_p$ divides every $z'_j$. $\square$

Because $N'\le N$ and $L'\le C\sqrt{R'}$, **every bound below that is proved for normalised
configurations and is monotone in $W$ transfers to the original configuration with $W'=\log N'$
replaced by $W=\log N$.** The number $M$ is unchanged.

**Standing assumption for Sections 2–6.** The configuration is normalised:
$\gcd(z_1,\dots,z_M)=1$, with $N=\prod_p p^{e_p}$ and $z_j=\varepsilon_j\prod_p\pi_p^{a_{jp}}\bar\pi_p^{\,e_p-a_{jp}}$
as in Lemma 1.1. The units $\varepsilon_j$ are arbitrary. No common unit class is imposed.

---

## 2. Pair factorisation and the exact cut identity

For $i\ne j$ let $g_{ij}=\gcd(z_i,z_j)$ and
$$c_{ij}:=\frac{z_i-z_j}{g_{ij}}\in\mathbb Z[i]\setminus\{0\}.$$
This is integral because $g_{ij}$ divides both $z_i$ and $z_j$, and non-zero because
$z_i\ne z_j$. (In the sketch's notation, $c_{ij}=u_i-u_j$ with $u_i=z_i/g_{ij}$.)

**Lemma 2.1 (pair gcd).**
$$g_{ij}\sim\prod_p \pi_p^{\min(a_{ip},a_{jp})}\bar\pi_p^{\,\min(e_p-a_{ip},\,e_p-a_{jp})},$$
hence
$$\log|g_{ij}|^2=W-d_{ij},\qquad d_{ij}:=\sum_{p\mid N}|a_{ip}-a_{jp}|\log p .$$

*Proof.* The first formula is unique factorisation. The exponent of $p$ in $|g_{ij}|^2$ is
$\min(a_i,a_j)+e-\max(a_i,a_j)=e-|a_i-a_j|$. $\square$

**Lemma 2.2 (cut identity).** For $p\mid N$ and $1\le\tau\le e_p$ put
$S_{p\tau}:=\#\{i:a_{ip}<\tau\}$. Then
$$\sum_{i<j}|a_{ip}-a_{jp}|=\sum_{\tau=1}^{e_p}S_{p\tau}(M-S_{p\tau})=\frac{e_pM^2}{4}-\sum_{\tau=1}^{e_p}\Big(S_{p\tau}-\frac M2\Big)^2 .$$

*Proof.* For integers $0\le a,b\le e$, $|a-b|$ is the number of $\tau\in\{1,\dots,e\}$ such that
exactly one of $a<\tau$ and $b<\tau$ holds. (These are the $\tau$ with
$\min(a,b)<\tau\le\max(a,b)$, all of which lie in $[1,e]$.) The pairs separated at $\tau$ number
$S_\tau(M-S_\tau)$. Finally, $S(M-S)=M^2/4-(S-M/2)^2$. $\square$

**Proposition 2.3 (exact identity).** Put
$$Q:=\sum_{p\mid N}\log p\sum_{\tau=1}^{e_p}\Big(S_{p\tau}-\frac M2\Big)^2\;\ge 0 .$$
Then
$$\sum_{i<j}\log|c_{ij}|+\tfrac12 Q=\sum_{i<j}\log|z_i-z_j|-\frac{M(M-2)}{4}\log R. \tag{2.3}$$
Equivalently, in integers,
$$\prod_{i<j}|z_i-z_j|^8=N^{M(M-2)}\prod_{p\mid N}p^{\sum_\tau(2S_{p\tau}-M)^2}\prod_{i<j}|c_{ij}|^8. \tag{2.4}$$

*Proof.*

1. By Lemma 2.1, $\log|z_i-z_j|=\log|g_{ij}|+\log|c_{ij}|=\log R-\tfrac12d_{ij}+\log|c_{ij}|$.
2. By Lemma 2.2,
   $$\sum_{i<j}d_{ij}=\frac{M^2}{4}\sum_p e_p\log p-Q=\frac{M^2}{2}\log R-Q.$$
3. Summing step 1 over pairs and substituting step 2:
   $$\sum_{i<j}\log|z_i-z_j|=K\log R-\frac{M^2}{4}\log R+\frac Q2+\sum_{i<j}\log|c_{ij}|.$$
4. Since $K-M^2/4=M(M-2)/4$, this is (2.3).
5. Multiplying (2.3) by $8$ and exponentiating gives (2.4). $\square$

**Corollary 2.4.** For a configuration on an arc of length $\le C\sqrt R$,
$$\sum_{i<j}\log|c_{ij}|+\tfrac12Q\;\le\;\frac M4\log R+K\log C. \tag{2.5}$$

*Proof.* We have $\log|z_i-z_j|\le\log C+\frac12\log R$. The right side of (2.3) is therefore at
most $K\log C+\log R\,\big(\frac{M(M-1)}{4}-\frac{M(M-2)}{4}\big)=K\log C+\frac M4\log R$. $\square$

(2.5) is the sketch's "summed pair inequality", divided by 2: there
$\gamma_{ij}=2\log|c_{ij}|$. Dropping $Q\ge0$ and using $S(M-S)\le\lfloor M^2/4\rfloor$ recovers (2)
of `inert_prime_cofactor_bound.md`, including its $b_M$.

---

## 3. Arithmetic lower bounds for the reduced chords

**Lemma 3.1 (distinct prime divisors multiply).** For $c\in\mathbb Z[i]\setminus\{0\}$ and an odd
rational prime $q$, let $v_q(c)=\max\{a:q^a\mid c\}=\sum_{a\ge1}\mathbf 1[q^a\mid c]$. Then
$$\log|c|\;\ge\;\tfrac12\log2\cdot\mathbf 1[(1+i)\mid c]\;+\sum_{q\text{ odd}}v_q(c)\log q .$$

*Proof.* The elements $(1+i)^{s}$, with $s=\mathbf 1[(1+i)\mid c]$, and $q^{v_q(c)}$ for odd $q$ all
divide $c$. Their norms $2^s$ and $q^{2v_q(c)}$ are pairwise coprime, so no Gaussian prime divides
two of them. By unique factorisation their product $P$ divides $c$. Hence
$|c|^2=|P|^2|c/P|^2\ge|P|^2=2^s\prod_q q^{2v_q(c)}$. $\square$

In particular, collisions at distinct primes multiply. If $c_{ij}$ is divisible by distinct odd
primes $q_1,\dots,q_r$, then $|c_{ij}|^2\ge\prod q_k^2$. This is the sketch's
$\gamma_{ij}\ge\sum 2\log q$.

**Lemma 3.2 (conic counts).** Let $q$ be an odd prime, $a\ge1$ and $n\in\mathbb Z$ with $q\nmid n$.
Then
$$\nu(q^a,n):=\#\{(x,y)\in(\mathbb Z/q^a)^2:\ x^2+y^2\equiv n\}=(q-\chi(q))\,q^{a-1}.$$
That is $(q+1)q^{a-1}$ for $q\equiv3$, and $(q-1)q^{a-1}$ for $q\equiv1\pmod 4$.

*Proof.* We treat $a=1$ first, then lift.

* *Level $a=1$, $q\equiv1$.* Pick $\iota\in\mathbb F_q$ with $\iota^2=-1$. The change of variables
  $u=x+\iota y$, $v=x-\iota y$ is invertible ($2\iota\neq0$), and $x^2+y^2=uv$. The equation
  $uv=n\neq0$ has exactly $q-1$ solutions.
* *Level $a=1$, $q\equiv3$.* Here $\mathbb F_q[i]=\mathbb F_{q^2}$ and
  $x^2+y^2=N_{\mathbb F_{q^2}/\mathbb F_q}(x+iy)=(x+iy)^{q+1}$. This norm map
  $\mathbb F_{q^2}^*\to\mathbb F_q^*$ is a surjective homomorphism (its kernel, the $(q+1)$-th roots
  of unity, has $q+1$ elements). So each fibre over $n\ne0$ has $(q^2-1)/(q-1)=q+1$ elements.
* *Lifting from $a$ to $a+1$.* Let $x^2+y^2\equiv n\pmod{q^a}$. Then
  $(x+q^as)^2+(y+q^at)^2\equiv x^2+y^2+2q^a(xs+yt)\pmod{q^{a+1}}$, because $2a\ge a+1$. So the
  lifts are the solutions $(s,t)\in\mathbb F_q^2$ of the linear equation
  $2(xs+yt)\equiv-(x^2+y^2-n)/q^a\pmod q$. Since $x^2+y^2\equiv n\not\equiv0$, the vector
  $(2x,2y)\not\equiv 0$, so there are exactly $q$ solutions. Distinct classes mod $q^a$ have
  distinct lifts. $\square$

A Gaussian integer $z=x+iy$ satisfies $q^a\mid z$ iff $q^a\mid x$ and $q^a\mid y$. So "$z$ mod $q^a$"
is the pair $(x,y)$ mod $q^a$, and $|z|^2=x^2+y^2$.

**Lemma 3.3 (collision criteria; normalised configuration).** Let $i\ne j$ and $a\ge1$.

* **(R) Ramified prime.** $(1+i)\mid c_{ij}$ for every pair.
* **(I) Inert $\ell$.** $\ell^a\mid c_{ij}\iff z_i\equiv z_j\pmod{\ell^a}$. All $z_j$ mod
  $\ell^a$ lie in a set of $(\ell+1)\ell^{a-1}$ classes.
* **(S$\nmid$) Split $p\nmid N$.** $p^a\mid c_{ij}\iff z_i\equiv z_j\pmod{p^a}$. All $z_j$ mod $p^a$
  lie in a set of $(p-1)p^{a-1}$ classes.
* **(S$\mid$) Split $p\mid N$, $e=e_p$, $\pi=\pi_p$.** Put
  $r_i:=z_i/(\pi^{a_{ip}}\bar\pi^{\,e-a_{ip}})\in\mathbb Z[i]$, so $|r_i|^2=N/p^e$ and
  $p\nmid N/p^e$. If $a_{ip}=a_{jp}$, then $p^a\mid c_{ij}\iff r_i\equiv r_j\pmod{p^a}$. For each
  level $t$, the $r_i$ with $a_{ip}=t$ lie mod $p^a$ in a set of $(p-1)p^{a-1}$ classes.
  If instead $a_{ip}\ne a_{jp}$, then neither $\pi$ nor $\bar\pi$ divides $c_{ij}$. This last fact
  is not used in the proof; it shows that different levels never collide at $p$.

*Proof.*

* *(R).* $N$ is odd, so no $z_j$ is divisible by $1+i$. Since $\mathbb Z[i]/(1+i)\cong\mathbb F_2$,
  every $z_j\equiv1\pmod{1+i}$, hence $(1+i)\mid z_i-z_j=g_{ij}c_{ij}$. Also
  $|g_{ij}|^2\mid|z_i|^2=N$ is odd, so $(1+i)\nmid g_{ij}$, and therefore $(1+i)\mid c_{ij}$.
* *(I) and (S$\nmid$).* In both cases $q\nmid N$ (for $q=\ell$ by Lemma 1.1(b)). Since
  $|g_{ij}|^2\mid N$, no prime above $q$ divides $g_{ij}$, so $\gcd(q^a,g_{ij})=1$. Then
  $q^a\mid g_{ij}c_{ij}=z_i-z_j$ is equivalent to $q^a\mid c_{ij}$ (Euclid's lemma in the UFD
  $\mathbb Z[i]$). The class counts follow from Lemma 3.2 with $n=N$, since $|z_j|^2=N$ and $q\nmid N$.
* *(S$\mid$), same level.* Let $a_{ip}=a_{jp}=t$ and $D=\pi^t\bar\pi^{e-t}$. Then $z_i=Dr_i$ and
  $z_j=Dr_j$, so $g_{ij}\sim D\,h_{ij}$ with $h_{ij}=\gcd(r_i,r_j)$, and
  $c_{ij}\sim(r_i-r_j)/h_{ij}$. Now $|h_{ij}|^2\mid N/p^e$ is prime to $p$, so $\gcd(p^a,h_{ij})=1$.
  As before, $p^a\mid c_{ij}\iff p^a\mid r_i-r_j$. The class count is Lemma 3.2 with $n=N/p^e$.
* *(S$\mid$), different levels.* Say $t=a_{ip}<a_{jp}=t'$. The $p$-part of $g_{ij}$ is
  $\pi^{t}\bar\pi^{e-t'}$. So $z_i/g_{ij}$ has $p$-part $\bar\pi^{t'-t}$, which is divisible by
  $\bar\pi$ but not $\pi$, and $z_j/g_{ij}$ has $p$-part $\pi^{t'-t}$. Their difference $c_{ij}$
  is therefore $\not\equiv0$ modulo $\pi$ and modulo $\bar\pi$. $\square$

*Units.* The units $\varepsilon_j$ sit inside $z_j$ and $r_j$. The residue classes above are those
of the actual $z_j$ and $r_j$, and the class counts cover **all** residues of the given norm. So
units affect neither the criteria nor the counts, and no "common unit class" normalisation is
needed.

**Lemma 3.4 (pigeonhole).** Let $n\ge0$ and $Q\ge1$, and write $b=\lfloor n/Q\rfloor$,
$r=n-Qb$. Then
$$E(n,Q)=Q\binom b2+rb\;\ge\;\frac{n^2}{2Q}-\frac n2 .\tag{3.5}$$

*Proof.* If two occupancies differ by $\ge2$, moving one object from the larger class to the smaller
strictly lowers $\sum\binom{m_r}2$. So the minimum is attained by $r$ occupancies $b+1$ and $Q-r$
occupancies $b$, which gives the formula. Also, by Cauchy–Schwarz,
$\sum\binom{m_r}{2}=\frac12(\sum m_r^2-n)\ge\frac12(n^2/Q-n)$. $\square$

**Corollary 3.6.** Write $n_{pt}:=\#\{i:a_{ip}=t\}$ for $p\mid N$ and $0\le t\le e_p$. Then
$$\sum_{i<j}\log|c_{ij}|\ \ge\ \frac{K\log2}{2}+\sum_{\ell}\log\ell\sum_{a\ge1}E\big(M,(\ell+1)\ell^{a-1}\big)
+\sum_{p\nmid N}\log p\sum_{a\ge1}E\big(M,(p-1)p^{a-1}\big)
+\sum_{p\mid N}\log p\sum_{a\ge1}\sum_{t=0}^{e_p}E\big(n_{pt},(p-1)p^{a-1}\big). \tag{3.6}$$

*Proof.* Apply Lemma 3.1 to each $c_{ij}$ and use $v_q=\sum_a\mathbf 1[q^a\mid\cdot]$:
$$\sum_{i<j}\log|c_{ij}|\ge\frac{K\log2}{2}+\sum_{q\text{ odd}}\log q\sum_{a\ge1}\#\{i<j:q^a\mid c_{ij}\},$$
using (R) for the first term. By Lemma 3.3, the pairs with $q^a\mid c_{ij}$ include all pairs in a
common residue class among at most the stated number of classes. For $p\mid N$ this uses only pairs
on a common level, and pairs from different levels are distinct pairs. Lemma 3.4 bounds each count
from below. $\square$

All sums in (3.6) are finite, because $E(n,Q)=0$ for $Q\ge n$.

---

## 4. The chain lemma

**Lemma 4.1 (chain lemma).** Let $\lambda>0$ and let $c=c(\lambda):=\frac{\sqrt{1+4\lambda}-1}{2}$ be
the positive root of $c^2+c=\lambda$. Let $e\ge0$ be an integer, $h>0$, and let
$x_0,x_1,\dots,x_{e+1}\in\mathbb R$ with $x_0=-h$ and $x_{e+1}=h$. Then
$$\Phi:=\sum_{\tau=1}^{e}x_\tau^2+\lambda\sum_{t=0}^{e}(x_{t+1}-x_t)^2\ \ge\ 2c\,h^2 .$$

In particular, take a distribution $(n_0,\dots,n_e)$ of $M$ points into levels, with
$S_\tau=n_0+\dots+n_{\tau-1}$. Put $x_\tau=S_\tau-M/2$, so that $x_0=-M/2$, $x_{e+1}=M/2$ and
$n_t=x_{t+1}-x_t$. Then
$$\Phi=\sum_{\tau=1}^{e}\Big(S_\tau-\frac M2\Big)^2+\lambda\sum_{t=0}^{e}n_t^2\ \ge\ c(\lambda)\frac{M^2}{2}. \tag{4.1}$$

*Proof.* The proof has four steps: a one-sided chain bound, a monotone recursion, the median split,
and the final estimate.

*Step 1: one-sided chains.* For $n\ge0$ and reals $y_0,\dots,y_n$ put
$$F_n(y_0,\dots,y_n):=\sum_{k=1}^{n}y_k^2+\lambda\sum_{k=0}^{n-1}(y_k-y_{k+1})^2+\lambda y_n^2 .$$
We claim that $F_n\ge c\,y_0^2$.

Define $\gamma_0=\lambda$ and $\gamma_n=T(\gamma_{n-1})$, where $T(\gamma):=\frac{\lambda(1+\gamma)}{1+\gamma+\lambda}$.
We show by induction that $F_n(y_0,\dots)\ge\gamma_n y_0^2$ for **all real** $y_0$.

- For $n=0$, $F_0=\lambda y_0^2$.
- For $n\ge1$, $F_n(y_0,\dots,y_n)=\lambda(y_0-y_1)^2+y_1^2+F_{n-1}(y_1,\dots,y_n)$. By induction,
  applied with starting value $y_1$ of either sign, this is
  $$\ge\lambda(y_0-y_1)^2+(1+\gamma_{n-1})y_1^2\ \ge\ \frac{\lambda(1+\gamma_{n-1})}{\lambda+1+\gamma_{n-1}}\,y_0^2 .$$
  Here we used $\min_x[\alpha(y-x)^2+\beta x^2]=\frac{\alpha\beta}{\alpha+\beta}y^2$ for
  $\alpha,\beta>0$.

*Step 2: the recursion stays above $c$.* $T$ is increasing on $[0,\infty)$, since
$T'(\gamma)=\lambda^2/(1+\gamma+\lambda)^2>0$. Its fixed point is $c$:
$T(c)=c\iff c(1+c+\lambda)=\lambda(1+c)\iff c^2+c=\lambda$. Since $c^2+c=\lambda$ and $c>0$, we have
$c<\lambda=\gamma_0$. By induction, $\gamma_{n-1}\ge c$ implies $\gamma_n=T(\gamma_{n-1})\ge T(c)=c$.
Hence $F_n\ge c\,y_0^2$ for every $n\ge0$.

*Step 3: split at a median.* Let $m:=\max\{t\in\{0,\dots,e\}:x_t\le0\}$. This exists because
$x_0=-h<0$, and since $x_{e+1}=h>0$ we have $x_{m+1}>0$. Because $x_m\le0<x_{m+1}$,
$$(x_{m+1}-x_m)^2=x_m^2+x_{m+1}^2-2x_mx_{m+1}\ \ge\ x_m^2+x_{m+1}^2 .$$
Split the sums of $\Phi$ at the step $t=m$, and replace the middle step term
$\lambda(x_{m+1}-x_m)^2$ by $\lambda x_m^2+\lambda x_{m+1}^2$. This gives
$$\Phi\ \ge\ F_m(y_0,\dots,y_m)+F_{e-m}(w_0,\dots,w_{e-m}),$$
with $y_k=-x_k$ ($0\le k\le m$, so $y_0=h$) and $w_k=x_{e+1-k}$ ($0\le k\le e-m$, so $w_0=h$).

To see this, note that the left block contains $x_1^2,\dots,x_m^2$, the steps $t=0,\dots,m-1$, and
the extra $\lambda x_m^2$. The right block contains $x_{m+1}^2,\dots,x_e^2$, the steps
$t=m+1,\dots,e$, and the extra $\lambda x_{m+1}^2$. The degenerate cases are included:

- if $m=0$, the left block is $F_0=\lambda h^2$;
- if $m=e$, the right block is $F_0=\lambda h^2$;
- if some $x_\tau=0$ (possible for integer $S_\tau$ with $M$ even), it is assigned to the left.

*Step 4: conclude.* By Steps 1–2, $\Phi\ge c h^2+c h^2=2ch^2$. $\square$

**Remarks 4.2.**

1. **Generality.** Nothing about integrality, the parity of $M$, the value of $e$, or empty levels
   (including levels $0$ and $e$) is used. The bound holds for all real profiles.
2. **Sharpness.** $\inf_e\min_{x\in\mathbb R^e}\Phi=c\,M^2/2$ exactly.
   For $e=2k+1$, take $x_{k+1}=0$ and mirror-symmetric halves, each optimal for $F_k$; this gives
   $\Phi=2\gamma_kh^2$. Moreover $\gamma_k\downarrow c$, since $0<T'(c)<1$. Check (C) of
   `check_chain_lemma.py` confirms the convergence numerically.
   So the chain constant cannot be improved. The integer minima exceed $cM^2/2$ only by lower-order
   amounts (checks (A) and (B)).
3. **Size of $c$.** $c\le\lambda$, and $c=\lambda-c^2\ge\lambda-\lambda^2$. For
   $\lambda=1/(p-1)$:
   $$c_p\log p\ \ge\ \frac{(p-2)\log p}{(p-1)^2}=\frac{\log p}{p}-\frac{\log p}{p(p-1)^2}.\tag{4.2}$$
4. **What the lemma fixes.** For $e=1$ (the squarefree case),
   $\Phi=\lambda M^2/2+(1+2\lambda)(n_0-M/2)^2\ge\lambda M^2/2$. So the chain lemma loses only the
   factor $c/\lambda=1-O(\lambda)$ against the squarefree model of
   `inert_prime_residual_growth.md` §3.1, uniformly in the exponent pattern. This removes the
   "dilution by exponent levels" of §5 of that note: spreading the points over many levels makes
   the cut slack $\sum_\tau(S_\tau-M/2)^2$ large.

**Lemma 4.3 (per split prime).** Let $p\equiv1\ (4)$, $\lambda_p=1/(p-1)$ and $c_p=c(\lambda_p)$.
In a normalised configuration define the $p$-contribution
$$\Gamma_p:=\log p\sum_{i<j}v_p(c_{ij})\;+\;\mathbf 1[p\mid N]\,\frac{\log p}{2}\sum_{\tau=1}^{e_p}\Big(S_{p\tau}-\frac M2\Big)^2 .$$
Then
$$\Gamma_p\ \ge\ \Big(\frac{c_pM^2}{4}-\frac M2\Big)\log p .$$

*Proof.* There are two cases.

* *If $p\nmid N$:* by (S$\nmid$) and (3.5),
  $$\Gamma_p\ge E(M,p-1)\log p\ge\Big(\frac{\lambda_pM^2}{2}-\frac M2\Big)\log p\ge\Big(\frac{c_pM^2}{4}-\frac M2\Big)\log p .$$
* *If $p\mid N$:* by (S$\mid$), (3.5) on each level, and (4.1),
  $$\Gamma_p\ \ge\ \log p\Big[\sum_tE(n_{pt},p-1)+\tfrac12\sum_\tau x_\tau^2\Big]
  \ \ge\ \log p\Big[\frac{\lambda_p}{2}\sum_tn_{pt}^2-\frac M2+\frac12\sum_\tau x_\tau^2\Big]
  =\frac{\log p}{2}\big(\Phi_p-M\big)\ \ge\ \Big(\frac{c_pM^2}{4}-\frac M2\Big)\log p. \qquad\square$$

---

## 5. The exact finite inequality

**Theorem 5.1 (exact finite inequality, configuration form).** Let $z_1,\dots,z_M$ ($M\ge2$) be
distinct lattice points on an arc of length $\le C\sqrt R$ of $x^2+y^2=R^2=N$. Normalise as in
Lemma 1.1, with new norm $N'=R'^2\le N$, levels $a_{jp}$, level sizes $n_{pt}$ and cut counts
$S_{p\tau}$. Then
$$\boxed{\begin{aligned}
&\frac{K\log 2}{2}+\sum_{\ell\equiv3\,(4)}\log\ell\sum_{a\ge1}E\big(M,(\ell+1)\ell^{a-1}\big)
+\sum_{p\equiv1\,(4),\ p\nmid N'}\log p\sum_{a\ge1}E\big(M,(p-1)p^{a-1}\big)\\
&\quad+\sum_{p\mid N'}\log p\Big[\sum_{a\ge1}\sum_{t=0}^{e_p}E\big(n_{pt},(p-1)p^{a-1}\big)+\frac12\sum_{\tau=1}^{e_p}\Big(S_{p\tau}-\frac M2\Big)^2\Big]
\ \le\ \frac M4\log R'+K\log C\ \le\ \frac M4\log R+K\log C .
\end{aligned}}\tag{5.1}$$

*Proof.* Combine (3.6) with (2.5) for the normalised configuration (Lemma 1.1, $L'\le C\sqrt{R'}$),
then use $R'\le R$. $\square$

To remove the dependence on the configuration, define for each split $p$ and each $M\ge2$
$$\Psi_p(M):=\min_{e\ge1}\ \min_{\substack{n_0,\dots,n_e\ge 0,\ \sum n_t=M\\ n_0\ge1,\ n_e\ge1}}\Big[\sum_{a\ge1}\sum_{t=0}^{e}E\big(n_t,(p-1)p^{a-1}\big)+\frac12\sum_{\tau=1}^{e}\Big(S_\tau-\frac M2\Big)^2\Big],$$
$$G_p(M):=\min\Big\{\sum_{a\ge1}E\big(M,(p-1)p^{a-1}\big),\ \Psi_p(M)\Big\},$$
$$F^*(M):=\frac{K\log2}{2}+\sum_{\ell}\log\ell\sum_{a\ge1}E\big(M,(\ell+1)\ell^{a-1}\big)+\sum_{p}G_p(M)\log p .$$

**Theorem 5.2 (configuration-free finite inequality).** For every configuration of $M\ge2$ points on
an arc of length $\le C\sqrt R$,
$$F^*(M)\ \le\ \frac M4\log R+\binom M2\log C. \tag{5.2}$$
$F^*(M)$ is a finite, exactly computable quantity:

- only $\ell<M$ and $p\le M$ contribute;
- the minimum in $\Psi_p$ exists and is attained with $e\le M-1$ and no empty level.

*Proof.*

1. *Each term is at least its minimum.* In (5.1), the bracket for $p\mid N'$ is an admissible
   value in the definition of $\Psi_p(M)$ (Lemma 1.1(d) gives $n_0,n_e\ge1$). For $p\nmid N'$, the
   $p$-term is the first entry of $G_p$. Hence the left side of (5.1) is $\ge F^*(M)$.
2. *Finiteness.* The bracket takes values in $\frac18\mathbb Z_{\ge0}$, so the minimum exists.
3. *No empty levels needed.* If $e\ge2$ and some middle level is empty ($n_t=0$, $0<t<e$), delete
   it. Then $e$ drops by one, $n_0,n_e$ are unchanged, the $E$-terms are unchanged, and exactly one
   repeated value $S_{t+1}=S_t$ disappears from the cut sum. So the bracket does not increase, and
   the minimum is attained with all $n_t\ge1$, hence with $e\le M-1$.
4. *Only small primes count.* $E(n,Q)=0$ for $Q\ge n$, and $G_p(M)=0$ once $p-1\ge M$. $\square$

**Remarks 5.3.**

1. **This extends (6) of the inert note.** Drop all split-prime terms and the term $K\log2/2$. For
   $M$ odd, also use $(S-M/2)^2\ge\frac14$, which gives $\frac12Q\ge\frac18\sum_pe_p\log p=\frac14\log R$.
   Then (5.1) becomes exactly $F(M)\le b_M\log R+K\log C$, inequality (6) of
   `inert_prime_cofactor_bound.md`. In general, for $M$ odd one may replace $(S_\tau-M/2)^2$ by
   $(S_\tau-M/2)^2-\frac14$ inside $\Psi_p$ and $M/4$ by $(M-1)/4$ on the right.
2. **Numerics** (evidence only; `finite_inequality.py`).

   | $M$ | $F^*(M)$ | $3\log M-8F^*(M)/M^2$ | $\log R_{\min}=4F^*(M)/M$ ($C=1$) |
   |---|---|---|---|
   | 16 | 108.8 | 4.92 | 27.2 |
   | 48 | 1686.7 | 5.76 | 140.6 |
   | 96 | 9022.3 | 5.86 | 375.9 |
   | 200 | 49482.5 | 6.00 | 989.7 |

   At $M=200$ the bound reads $3M\log M\le W+6.0\,M$. That is, for $C=1$, at least 200 points
   need $\log R\ge 989.7$. At $M=200$, for 21 split primes $p$ the $p\mid N$ branch $\Psi_p$ is the
   smaller one: the worst case has many small split primes dividing $N$ with balanced levels.

---

## 6. Depth-one consequence

Keep only the following terms in (5.1):

- the ramified term;
- inert $\ell\le M$ at depth $a=1$, via (3.5);
- split $p\le M$, via Lemma 4.3.

Every discarded term is $\ge0$. The result is
$$D(M):=\frac{K\log2}{2}+\sum_{\ell\le M}\Big(\frac{M^2}{2(\ell+1)}-\frac M2\Big)\log\ell+\sum_{p\le M}\Big(\frac{c_pM^2}{4}-\frac M2\Big)\log p\ \le\ \frac M4\log R'+K\log C .\tag{6.2}$$
Individual summands of $D(M)$ may be negative; (6.2) remains valid, since each summand is a lower
bound for a non-negative quantity.

Multiplying (6.2) by $8/M$, with $W'=\log N'=2\log R'$, gives
$$M\big(4A(M)+2B(M)\big)-4\theta_{\rm odd}(M)+2(M-1)\log 2\ \le\ W'+4(M-1)\log C ,\tag{6.3}$$
where
$$A(x):=\sum_{\ell\le x}\frac{\log\ell}{\ell+1},\qquad B(x):=\sum_{p\le x}c_p\log p,\qquad \theta_{\rm odd}(x):=\theta(x)-\log 2 .$$

*Choice $x=M$.* Primes $q>M$ are dropped (their terms are $\ge0$). Inside $q\le M$ every prime
keeps its $-\frac M2\log q$ correction, which costs $-\frac M2\theta_{\rm odd}(M)=O(M^2)$, the same
order as the other error terms. Any $x=\kappa M$ with fixed $\kappa>0$ gives the same leading term.

---

## 7. Explicit prime sums and the effective inequality

Define
$$O(x):=\sum_{3\le q\le x}\frac{\log q}{q},\qquad X(x):=\sum_{q\le x}\chi(q)\frac{\log q}{q},\qquad S_1=\tfrac12(O+X),\quad S_3=\tfrac12(O-X).$$
Here $S_1(x)=\sum_{p\le x}\frac{\log p}{p}$ and $S_3(x)=\sum_{\ell\le x}\frac{\log\ell}{\ell}$.
These are the Mertens sums in the two odd classes mod 4. We only need one-sided explicit bounds.
In particular we need no asymptotic for $X$, which in fact converges to a negative constant.

**Lemma 7.1 (Chebyshev-type bounds).** For real $y\ge1$:
$$\theta(y)<y\log4\qquad\text{and}\qquad \psi(y)<y\log 4+\sqrt y\log y\le y\big(\log4+2/e\big).$$

*Proof.*

1. *Integer bound.* We show $\prod_{q\le n}q\le4^{n-1}$ for integers $n\ge1$, by induction.
   - $n=1,2$: clear.
   - Even $n>2$: the product equals the one for $n-1$.
   - $n=2m+1$: every prime in $(m+1,2m+1]$ divides $\binom{2m+1}{m}$, and
     $2\binom{2m+1}{m}=\binom{2m+1}{m}+\binom{2m+1}{m+1}\le2^{2m+1}$. So
     $\prod_{q\le 2m+1}q\le4^{m}\cdot4^{m}=4^{n-1}$.
2. *$\theta$.* For real $y\ge1$, $\theta(y)=\theta(\lfloor y\rfloor)<\lfloor y\rfloor\log4\le y\log 4$.
3. *$\psi$.* $\psi(y)-\theta(y)=\sum_{q\le\sqrt y}(\lfloor\log y/\log q\rfloor-1)\log q\le\sum_{q\le\sqrt y}(\log y-\log q)<\sqrt y\log y$.
   Finally, $\log y/\sqrt y\le 2/e$. $\square$

**Lemma 7.2 (Mertens, lower bound).** For integers $n\ge2$:
$$O(n)\ \ge\ \log n-\beta+\frac1n,\qquad \beta:=1+\sigma+\tfrac12\log2,\qquad \sigma:=\sum_q\frac{\log q}{q(q-1)} .$$

*Proof.* By Legendre's formula,
$$\log n!=\sum_{q\le n}\log q\sum_{k\ge1}\lfloor n/q^k\rfloor\le n\sum_{q\le n}\frac{\log q}{q-1}=n\sum_{q\le n}\Big(\frac{\log q}{q}+\frac{\log q}{q(q-1)}\Big).$$
Also $\log n!\ge\int_1^n\log t\,dt=n\log n-n+1$. Hence $\sum_{q\le n}\frac{\log q}{q}\ge\log n-1+\frac1n-\sigma$.
Subtract the $q=2$ term $\frac12\log2$. $\square$

**Lemma 7.3 ($\chi_4$-twisted sum, upper bound).** For real $y\ge1$:
$$X(y)\le\frac4\pi\cdot\frac{\psi(y)}{y}<\xi:=\frac4\pi\Big(\log4+\frac2e\Big)<2.70189 .$$

*Proof.*

1. *Factorisation.* Since $\log n=\sum_{d\mid n}\Lambda(d)$ and $\chi$ is completely multiplicative,
   $$T(y):=\sum_{n\le y}\frac{\chi(n)\log n}{n}=\sum_{d\le y}\frac{\chi(d)\Lambda(d)}{d}\,\mathcal L(y/d),\qquad \mathcal L(t):=\sum_{m\le t}\frac{\chi(m)}{m}.$$
2. *The partial sums $\mathcal L(t)$.* $\mathcal L(t)$ is a partial sum of Leibniz's series
   $1-\frac13+\frac15-\dots=\frac\pi4$. By the alternating-series bound,
   $\big|\mathcal L(t)-\frac\pi4\big|\le\frac1{m_0}<\frac1t$ for $t\ge1$, where $m_0$ is the least
   odd integer $>t$.
3. *Substitute.* With $Y(y):=\sum_{d\le y}\chi(d)\Lambda(d)/d$,
   $$\frac\pi4\,Y(y)=T(y)-\sum_{d\le y}\frac{\chi(d)\Lambda(d)}{d}\Big(\mathcal L(y/d)-\frac\pi4\Big)\le T(y)+\sum_{d\le y}\frac{\Lambda(d)}{d}\cdot\frac dy\le T(y)+\frac{\psi(y)}{y}.$$
4. *Sign of $T$.* $T(y)=\sum_{k\ge1,\,2k+1\le y}(-1)^k\frac{\log(2k+1)}{2k+1}$, since $\chi(2k+1)=(-1)^k$.
   The terms are strictly decreasing in absolute value for $2k+1\ge3>e$ and start with a negative
   sign. So every partial sum lies in $[-\frac{\log3}{3},0]$, and $T(y)\le0$. Hence
   $Y(y)\le\frac4\pi\psi(y)/y$.
5. *Prime powers.* $X(y)=Y(y)-\sum_{q\ \rm odd}\log q\sum_{k\ge2,\ q^k\le y}\chi(q)^kq^{-k}$. For each
   $q$ the inner sum is $\ge0$. For $\chi(q)=1$ every term is positive. For $\chi(q)=-1$ it is a
   partial sum of an alternating series with decreasing terms that starts with $+q^{-2}$. Hence
   $X(y)\le Y(y)$.
6. *Conclude* with Lemma 7.1. $\square$

**Lemma 7.4 (numerical constants).** Let
$$\sigma=\sum_q\frac{\log q}{q(q-1)},\qquad a_3:=\sum_\ell\frac{\log\ell}{\ell(\ell+1)},\qquad a_1:=\sum_p\frac{\log p}{p(p-1)^2}.$$
Then
$$\sigma\in[0.7553665,\,0.7553683],\qquad a_3\le0.1743774,\qquad a_1\le0.0225136 .$$

*Proof (finite computation, `constants.py`).* We take partial sums over primes $\le P=10^7$ and add
explicit tail bounds. The function $f(t)=\log t/(t(t-1))$ is decreasing for $t\ge2$, so
$$\sum_{n>P}f(n)\le\int_P^\infty f\le\frac{P}{P-1}\cdot\frac{\log P+1}{P}<1.8\cdot10^{-6}.$$
Double-precision rounding over fewer than $10^6$ terms is $<10^{-9}$, and this is added as a margin. $\square$

**Theorem 7.5 (effective inequality).** Let $N\ge1$, $R=\sqrt N$, $W=\log N$ and $C>0$. Every set of
$M\ge2$ distinct lattice points on an arc of length $\le C\sqrt R$ of $x^2+y^2=N$ satisfies
$$3M\log M\ \le\ W'+(K_0+4\log^+C)\,M\ \le\ W+(14+4\log^+C)\,M,$$
where $W'=\log N'$ is the normalised norm and
$$K_0:=3\beta+\xi+2(2a_3+a_1)+6\log2<13.9092 .$$
Equivalently, $M\log M\le(W+K_CM)/3$ with $K_C=14+4\log^+C$. (For $M=1$ the inequality is
trivial.)

*Proof.*

1. *Rewrite $A$ and $B$.* $A(x)=S_3(x)-\sum_{\ell\le x}\frac{\log\ell}{\ell(\ell+1)}\ge S_3(x)-a_3$.
   By (4.2), $B(x)\ge S_1(x)-a_1$. Therefore
   $$2A+B\ \ge\ 2S_3+S_1-c_0=\tfrac32O-\tfrac12X-c_0,\qquad c_0:=2a_3+a_1 .$$
2. *Apply the prime-sum lemmas at $x=M$.* Insert this into (6.3), together with
   - $O(M)\ge\log M-\beta+1/M$ (Lemma 7.2),
   - $X(M)\le\xi$ (Lemma 7.3),
   - $\theta_{\rm odd}(M)<M\log4-\log2$ (Lemma 7.1).

   This gives
   $$3M\log M-3\beta M+3-\xi M-2c_0M-4M\log4+4\log2+2(M-1)\log2\ <\ W'+4(M-1)\log C ,$$
   that is,
   $$3M\log M<W'+M(3\beta+\xi+2c_0+6\log2)-3-2\log2+4(M-1)\log C .$$
3. *Simplify.* $4(M-1)\log C\le4M\log^+C$, and $-3-2\log 2<0$. This proves the first inequality.
   The second follows from $W'\le W$ (Lemma 1.1) and $K_0<14$ (Lemma 7.4:
   $3\beta\le6.30583$, $\xi\le2.70189$, $2c_0\le0.74254$, $6\log2=4.15888$). $\square$

*How sharp is $K_0$ (evidence only).* The proof's lossy steps are Lemma 7.1, which is crude
($\theta(y)/y<1$ in reality), and Lemma 7.3, where in reality $X(n)\le0$ for all
$n\le2\cdot10^6$. Evaluating $D(M)$ directly, the constant actually needed for $C=1$ is
$\max_{2\le M\le2\cdot10^6}\big(3\log M-8D(M)/M^2\big)=7.77$, and the full $F^*(M)$ needs about
$6.0$ at $M=200$.

---

## 8. The asymptotic statement

**Corollary 8.1.** Let $C>0$ and $0<\varepsilon\le\frac13$, and put $K=K_C=14+4\log^+C$. Let
$\ell_0=\ell_0(C,\varepsilon)\ge1$ be such that for all $\ell\ge\ell_0$:
$$\text{(i)}\ \ \frac{\log\ell+\log\frac32+\frac K3}{\ell}\le\frac{3\varepsilon}{2+3\varepsilon},\qquad
\text{(ii)}\ \ \frac23\cdot\frac{e^{\ell}}{\ell}\ge e^{K/3}.$$
Both conditions are monotone in $\ell\ge1$, because $(\log\ell+a)/\ell$ decreases for
$a\ge1,\ \ell\ge1$, and $e^\ell/\ell$ increases. So it suffices to check them at $\ell_0$. Put
$R_0:=\exp(\exp\ell_0)$.

Then for every $R\ge R_0$ with $R^2\in\mathbb N$, every arc of $x^2+y^2=R^2$ of length at most
$C\sqrt R$ contains at most
$$M\le\Big(\frac23+\varepsilon\Big)\frac{\log R}{\log\log R}$$
lattice points. For $\varepsilon>\frac13$, use $\varepsilon=\frac13$.

*Proof.* Write $\Lambda=\log R$ and $\ell=\log\Lambda\ge\ell_0$, and put
$M_0:=(\frac23+\varepsilon)\Lambda/\ell$. Suppose $M>M_0$.

1. By (ii), $M_0\ge e^{K/3}$, so $M\ge2$ and Theorem 7.5 applies:
   $h(M):=M(3\log M-K)\le2\Lambda$.
2. Since $h'(M)=3\log M+3-K>0$ for $M>e^{K/3-1}$, $h$ is strictly increasing on $[M_0,\infty)$.
   So $h(M)>h(M_0)$.
3. Using $\log M_0=\ell-\log\ell+\log(\frac23+\varepsilon)$ and $\log(\frac23+\varepsilon)\ge-\log\frac32$,
   $$h(M_0)=(2+3\varepsilon)\Lambda\Big(1-\frac{\log\ell-\log(\frac23+\varepsilon)+K/3}{\ell}\Big)\ge(2+3\varepsilon)\Lambda\Big(1-\frac{3\varepsilon}{2+3\varepsilon}\Big)=2\Lambda,$$
   by (i).
4. Hence $2\Lambda\ge h(M)>h(M_0)\ge2\Lambda$, a contradiction. $\square$

For example, with $C\le1$ (so $K=14$), one can take
$\ell_0=24.9,\ 71.7,\ 795.2$ for $\varepsilon=\frac13,\ 0.1,\ 0.01$ (`constants.py`).
These thresholds are astronomically large but explicit. Qualitatively, Theorem 7.5 gives
$$M\le\frac23\,\frac{\log R}{\log\log R}\Big(1+O_C\Big(\frac{\log\log\log R}{\log\log R}\Big)\Big).$$

Uniformity: $R_0$ depends only on $(C,\varepsilon)$. The proof uses only $|z_i-z_j|\le C\sqrt R$,
so the bound is uniform in the position of the arc and needs no condition on the factorisation of
$R^2$.

---

## 9. Audit checklist (the points the task asked to verify)

1. **Residue-class counts $(\ell+1,\ p-1,\ p-1)$.** These are Lemma 3.2 at $a=1$, together with the
   reduction in Lemma 3.3.
   - Inert $\ell$: the norm fibre in $\mathbb F_{\ell^2}^*$ has $\ell+1$ elements.
   - Split $p\nmid N$: $\{(X,Y):XY=N\}$, i.e. $p-1$ classes.
   - Split $p\mid N$, same level: the residues $r_i$ mod $p$ of norm $N/p^{e}\not\equiv0$, i.e.
     $p-1$ classes.

   Depth $a$ multiplies each count by $p^{a-1}$ or $\ell^{a-1}$ (Hensel). For split $p\nmid N$,
   collision mod $\pi$ already forces collision mod $p$, because $\bar z\equiv N/z$ mod $\pi$.
   So no half-gains at $\pi$ alone are lost. Different levels never collide at $p$ (Lemma 3.3).
   *Verified exactly* in `check_identity.py`, which brute-forces the conic counts for all
   $q^a\le50$.
2. **Collisions at distinct primes multiply into $|u_i-u_j|^2$.** This is Lemma 3.1, by coprimality
   and unique factorisation. The ramified prime contributes a further factor $2$ to
   $|c_{ij}|^2$ for every pair (Lemma 3.3(R)).
3. **Chain lemma.** Proved for all real profiles (Lemma 4.1).
   - Degenerate medians ($x_\tau=0$, $m=0$, $m=e$) are handled explicitly.
   - $M$ odd needs nothing special.
   - Levels $0$ and $e$ are the boundary steps $n_0=x_1-x_0$ and $n_e=x_{e+1}-x_e$.
   - The constant $c(\lambda)$ is sharp (Remark 4.2).

   *Verified exactly* in `check_chain_lemma.py`:
   - 2,034,690 compositions for $M\le12$, $e\le7$ and seven values of $\lambda$;
   - a DP up to $M=40$, $e=24$;
   - the exact continuous optima;
   - the recursion;
   - the integer bracket bound $\ge c_pM^2/4-M/2$;
   - the restricted $\Psi_p$ including prime powers.

   The smallest ratio $\min\Phi/(cM^2/2)$ found was $1.000999$.
4. **Units.** These are irrelevant (Section 3, end). No unit-class normalisation is used, so no
   factor 4 is lost.
5. **Choice $x\sim M$ and error terms.** Section 6: every error term is explicit (Lemmas 7.1–7.4).
   The only $O(M)$ losses in Theorem 7.5 come from:
   - the Mertens constant;
   - the $\chi_4$ bound;
   - the $-\frac M2\log q$ pigeonhole corrections, totalling $-4\theta(M)$ after scaling;
   - the $C$-term.
6. **Arc position.** Only the chord bound is used (Section 0).
7. **Normalisation.** Dividing by the global gcd keeps the points distinct, does not increase the
   radius, and keeps $L'\le C\sqrt{R'}$ (Lemma 1.1). The final inequality is monotone in $W$.
   No per-point conjugate-primitivity and no single unit class is assumed.
8. **Exact identity.** (2.4) was *verified exactly* (big-integer equality) on 6,176 configurations,
   both arcs and random subsets, of eight norms. These norms include inert and ramified factors
   (exercising the normalisation) and repeated split primes up to $5^4$, $13^2$ and $17^2$. The
   same run checked every collision criterion of Lemma 3.3 and inequality (5.1) with the actual
   $C=L/\sqrt R$.

---

## 10. What this does not do

* **Not uniform.** The method gives a total lower bound of
  $(\frac14+\frac18)M^2\log M=\frac38M^2\log M$ (inert $+$ split) against the diagonal cut slack
  $\frac M4\log R$. Configurations with $M=o(\log R/\log\log R)$ are not excluded.
* **2/3 is the limit of these ingredients.** The obstacles are:
  - the chain constant is sharp;
  - the pigeonhole bounds are sharp up to $O(M)$ per prime;
  - deeper prime powers add only $O(M^2)$ in total, since
    $\sum_q\sum_{a\ge2}\frac{\log q}{q^{a-1}(q\pm1)}<\infty$.

  So with every small split prime dividing $N$ and levels following the optimal chain profile, the
  arithmetic lower bound is $\frac38M^2\log M+O(M^2)$. Whether actual configurations realise this
  worst case was **not** examined.
* No Lean formalisation is claimed.

---

## 11. Files and reproducibility

All files are in `scratchpad/two_thirds/` and use pure Python 3 with integer or `Fraction`
arithmetic, without sympy.

* `gauss.py`: Gaussian-integer helpers (Euclidean gcd, split-prime factor, valuations,
  brute-force lattice points).
* `check_chain_lemma.py`: the exact checker for Lemma 4.1, checks (A)–(F) as described in its
  docstring. Output ends with `ALL CHECKS PASSED: True` (about 2.5 min).
* `check_identity.py`: exact checks of Lemmas 1.1, 2.1 and 3.3, identity (2.4) and inequality (5.1)
  on actual configurations. Output ends with `ALL IDENTITY/COLLISION CHECKS PASSED` (about 25 s).
* `constants.py`: the constants of Lemma 7.4 with tail bounds, and $K_0<13.9092$; a sanity check
  that $3\log M-8D(M)/M^2\le K_0$ for all $2\le M\le2\cdot10^6$ (maximum 7.772); numerical
  illustration of Lemmas 7.1–7.3; and the thresholds $\ell_0$ of Corollary 8.1.
* `finite_inequality.py`: evaluates the exact configuration-free $F^*(M)$ of Theorem 5.2.

Excerpt of the `check_chain_lemma.py` output:

```text
compositions enumerated: 2034690;  min ratio minPhi/(cM^2/2) = 1.000999
(A)+(C)+(E) all exact checks passed: True
DP == brute force on M<=10, e<=5, p in {5,13}: True
   5   40   1..16         167.0000  165.6854  1.0079       165.6854  1.0000
 101   31   1..24           5.0100    4.7579  1.0530         4.7579  1.0000
 lam=1/4    c=0.207107  Phi_cont/(cM^2/2) for e=1,2,3,5,9,17,...: 1.207107 1.034663 1.005922 1.000174 1.000000 ...
ALL CHECKS PASSED: True
```

## References (classical inputs)

* Unique factorisation in $\mathbb Z[i]$ and the splitting of rational primes.
* Hardy & Wright, *An Introduction to the Theory of Numbers*, §22: $\theta(n)<n\log4$; the proof is
  reproduced in Lemma 7.1.
* Leibniz's series and the alternating-series remainder bound.
* For context only, not used: K. S. Williams, *Mertens' theorem for arithmetic progressions*,
  J. Number Theory 6 (1974) 353–359, gives the two-sided asymptotics
  $S_{1},S_3=\frac12\log x+O(1)$.
