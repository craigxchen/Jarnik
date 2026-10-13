# Referee report on `round9/L4.md` (L4 = (I1), the inverse step of P3-all\*)

Referee: adversarial, independent. Scripts: `round9/referee_L4_checks/` (my own code; nothing imported from
`L4_checks/`). Each script runs with `python3 <script>` from that folder; outputs are stored next to it as
`*_out.txt`. Working notes: `referee_L4_checks/notes.md`. The repository was not modified.

## 0. Verdict

**The proved theorems survive; two of the headline consequences do not.**

*Survives (re-derived step by step, with exact checks where possible):*
* the reduction (T') to regime (1.4) (`p_j <= M^(h_0)`, `X = 3 M^(A(c_0))`), from `round8/stability.md` Cor. 3.6(1);
* Lemma 2.1, Theorem A and Corollary A.1 (unconditional, via Linnik). Corollary A.1 can in fact be strengthened, see §3;
* Theorem B (via Heath-Brown's quadratic large sieve), with two harmless slips;
* Lemma 4.2 (the duality step was checked exactly on 3000 instances, 819 of them with relations), Theorem C, and Theorem C' (GRH);
* Lemma 5.1 (exhaustive for `16 <= N <= 1200`), Proposition D as a divisibility statement, and Corollary D.1.
  My own implementation reproduces `check_D.py`: no failure at norm `<= 2`, 202-446 failures at norm 3, none on
  the moduli of cases 1-3.

*Does not survive as stated.*
* **(M1) Regime mismatch.** The summary says Proposition D is "uncertified in every sense of Theorems A-C". It
  says D has "slack (3/2) log M in the construction regime". It concludes "Consequently L4 as stated is false up to
  norm 2, and any inverse theorem must use monomials of norm >= 3 and only the moduli `q = +-1 mod 2^B`".
  * The two properties of D hold in **disjoint** regimes.
    * Theorems A-C need `p_j Norm zeta <= X^(1/L)` (or `<= D <= X`), so they need `p_j << X`: that is regime (1.4).
    * The slack of Corollary D.1 needs `p_j >= M^(16h+3) > X`. There Theorems A-C are vacuous.
  * Corollary D.2's clause "or in the setting of Corollary D.1" is therefore wrong.
  * In regime (1.4), which by the write-up's own (T') is the only regime growth needs, Construction D fails the
    per-unit height condition. In exact instances it fails for 93-134 of the 100-144 norm-2 monomials, and for 8-10
    of the 10-12 norm-1 monomials, while the genuine labels on the same profile have zero violations (§6.3).
  * So nothing in the write-up shows that an inverse theorem in regime (1.4) must use norm `>= 3` or deep 2-adic
    moduli. Those "must" statements are proved only for inverse theorems that also have to hold when `p_j > X`.
* **(M2) (I1\*) is not "the precise remaining inverse step".** Its agreement notion is per column, with a unit
  allowed to vary with `(j, q)` (Def. 3.1). Its threshold is `1/2 + delta` of the `Q_7` log-mass.
  * Against the **genuine** labels, a uniformly random parity-respecting column already agrees, in expectation, on
    0.540-0.555 of the `Q_7` log-mass. This is exact per modulus, with 40 columns, `X = 2*10^5` (§7.2).
  * So for `delta < 0.04` the conclusion of (I1\*) is met by structureless labels. Calibrating it would need a
    random baseline above 1/2, not the parity baseline 1/2.
  * Per-column agreement does not control the collision sets of multi-column monomials. So "(I1\*) + Theorem B,
    then only (I3) remains" has a missing transfer step.
  * The formal merge "(I1) + (I2) ⟺ (I1\*) given Theorem B" is correct, but only for this weak notion of (I1).

*Medium.*
* The bold claim of §1 (also in §0 and the task summary) omits the norm bound: "no column is certified ... by
  Gaussian integers with `p_j` not dividing `Norm zeta_j`". Without the bound it is false: an explicit
  `zeta` over a prime `p' = p mod 4 prod_(q <= X, q = 7 (8)) q` is sub-genuine and parity-exact on all of
  `Q_7 ∩ [1, X]` (§7.3).
* Proposition D in the construction regime is certified in the generalised sense `m_v | Im(Z_v^16)`, by Gaussian
  integers of controlled height `M^h << p_j`. That is the plain-language conclusion of round 7's L4 ("induced by
  Gaussian integers `zeta_j` of controlled height"). Whether D refutes L4 there depends on a definition of
  "certified" that the write-up fixes without flagging it.
* The §4 remark on twisted coherence: Theorem C is vacuous for `p3Dlower.md` Def. 4.1, which requires
  `Norm zeta_j = p'_j > X`, i.e. `h > A`. Only Theorem C' (GRH) applies there.

*Minor.* Listed in §§2-7. They concern constants, a GRH exponent, the mass `2^(1-B)` (not `2^(-B)`),
"independent modulo squares" (should be "modulo squares and 2"), unlabelled heuristics in §5 and §6.3, and the
definition of "genuine height".

**Actual circles.** Nothing in the write-up bounds actual circles beyond `two_thirds`, and it says so. No growth
claim is made. Every theorem is about label systems (designs) in P3-all\*. This scope statement is correct.

## 1. What was checked and how

* Read in full: `L4.md`, `L4_notes.md`, `round7/README.md` §§0-3, `round8/stability.md` §§1, 3, `sharp.md` Def. 1.1,
  Prop. 1.2 and Prop. 6.5, and `p3Dlower.md` §§3-4.
* Every numbered statement was re-derived by hand (notes in `notes.md`).
* Exact computations, all with my own library `glib.py`. It does Gaussian arithmetic, works in the norm-one torus
  `T_(q^a)`, finds a generator by `w/conj w`, and takes discrete logs by Pohlig-Hellman, each asserted by
  exponentiation. `t_glib.py` tests the logs and the parity fact `log(w/conj w) = [(Norm w/q) = -1] mod 2`.
  * `consD.py` implements Construction D literally, from the text of §5.
  * The scripts `c1`-`c7` are listed in §8.
* Not checked: the literature citations. Linnik's theorem with `L = 5` (Xylouris) is quoted from memory by the
  author; so is the exact form of Heath-Brown's quadratic large sieve (Acta Arith. 72 (1995), Thm 1). My own
  recollection matches the form (3.1): `sum*_(m<=M) |sum*_(n<=N) a_n (n/m)|^2 <<_eps (MN)^eps (M+N) sum* |a_n|^2`,
  sums over odd squarefree positive integers. Only the size of `A` depends on these citations.

## 2. Section 1 (setting, F1, F2, regime reduction T')

**(F1) Correct.** Fix `ell' = u ell + (m_Q/4) w` with `u in (Z/m_Q)^*` and `q in Q_1`.
* For every `v`, `<v,ell'> in (m_Q/4)Z` iff `<v,ell> in (m_Q/4)Z`.
* At a lower level `a`, `m_Q/4 = (m_(q^a)/4) q^(A-a)`, so the collision sets agree at every level.
* `u` is odd because `m_Q` is even, and `m_Q/4` is even because `8 | q - chi(q)`. So the parities agree.

Remark: (F1) is a *system-level* operation, with one `u` for all columns. Definition 3.1 and (2.1) instead allow
the unit to depend on the column. That weaker notion of agreement is the one used in Theorems A-C. It makes those
theorems stronger, but it matters for (I1\*) (§7.2).

**(F2) Correct.** Take `q = 7 mod 8` with `q` not dividing `n`. Then `(2/q) = 1`, and Jacobi reciprocity with
`(q-1)/2` odd gives `(m/q) = (-1)^((m-1)/2) (q/m)`. Every `q = 3 mod 4` is coprime to `N`.

**(F3), the reduction (T'). Verified.**
* Suppose `W <= c_0 M log M` with `M -> infinity`.
  * Then `log R <= (c_0/2) M log M`, and `log M <= loglog R` for large `R` (stability (T3)).
  * Hence `M >= (2/c_0) log R/loglog R`, i.e. `H(eps)` holds with `eps = 2/c_0 <= 1`. This needs `c_0 >= 2`,
    which is automatic since `c_0 >= 3 - o(1)` by `two_thirds`.
  * Stability Cor. 3.6(1) then gives `P_+(N) <= M^(4/eps - 3 + o(1)) = M^(2c_0 - 3 + o(1))`.
* Conversely, if `M = o(log R/loglog R)` failed, there would be `eps > 0` and actual clusters with
  `M >= eps log R/loglog R`. These have `M -> infinity` and `W <= (2/eps + o(1)) M log M`.
  * Normalised clusters are members of `F^all_X` for every `X` (`sharp.md` Prop. 1.2).
  * Normalisation does not increase `W` or `C`.
  * The primes of `N'` are among those of `N`.
* So (T'), quantified as `forall c_0 exists A`, implies growth. This is a correct and useful observation: the
  quantifier order `forall c_0 exists A` is all that growth needs, and it puts the problem in the regime
  `X >> p_j^C`. There the unconditional least-prime theorems bite, unlike the regime `p_j > X` of rounds 5-7.

*Consequence the write-up does not draw:* in regime (1.4) the per-unit height condition is already a strong
constraint at norm 1.
* For `v = e_j` it says that the moduli `q <= X` where column `j`'s label lies in `mu_4` have product
  `<= nu_s^(-1) sqrt(p_j)`, where `X >> p_j^C`.
* A random parity-respecting column has expected log-demand about `4 sum_(q <= X) log q/q ~ 4 log X >> log p_j`.
* A column avoids every hit at the primes `q in (p_j^2, X]` with probability only about
  `prod (1 - 4/q) ~ (2 log p_j/log X)^4`. So among `r ~ M` independent columns, some column fails at norm 1 with
  probability `-> 1`.
* This is a routine remark, used nowhere below as a proof. In regime (1.4), independent random labels fail already
  at norm 1. That changes which constructions are relevant there (§6.3).

## 3. Theorem A / Cor A.1 (Linnik exact parity lock)

**Lemma 2.1: correct.**
* `n = 2^e m s^2` with `m > 1` odd squarefree. Then `k = lcm(8, rad n) <= 8n`.
* The Jacobi symbol `(./m)` is non-principal mod `m`, so some `b` has `(b/m) = -sigma_m`.
* By CRT choose `a = 7 (8)`, `a = b (m)`, and `a = 1` at the odd primes of `s` that do not divide `m`.
* Every prime `q = a mod k` is `7 mod 8`, does not divide `n`, and has `(n/q) = (m/q) = sigma_m (q/m) = -1`.
  Linnik then gives such a `q <= C_L k^L`.

**Theorem A: correct.**
* At `q = 7 mod 8` with `q` not dividing `n`, (2.1) together with the forced parity gives `(n/q) = 1`. The parity of
  `u log(zeta/conj zeta)` is that of `log(zeta/conj zeta)`, which is `[(Norm zeta/q) = -1]` (re-tested exactly in
  `t_glib.py`).
* `q` divides `n` iff `q` divides `Norm zeta`, since `q != p` and `q` is odd.
* Lemma 2.1 then forces `n in {t^2, 2t^2}`, so `v_p(Norm zeta)` is odd. In particular `p | Norm zeta` and
  `|zeta|^2 >= p`.
* *Precision.* The invariant is "`v_(p_j)(Norm zeta)` odd", not "`p_j | Norm zeta`". Take
  `zeta = p_j pi_(sigma(j))`. It has `p_j | Norm zeta`, but its labels are those of the sub-genuine `pi_(sigma(j))`,
  because `p_j/conj p_j = 1`. Theorem B and Theorem C extend verbatim to "`v_(p_j)(Norm zeta)` even", because
  `p_j` still divides the odd squarefree part `m`. So only the wording of §6.1 needs fixing: "genuine height" should
  mean `v_(p_j)(Norm zeta)` odd (equivalently, the primitive part of `zeta` is divisible by `pi_j` or `conj pi_j` to
  an odd power), or simply `|zeta_j|^2 >= p_j`.

**Corollary A.1: correct, and it can be strengthened.**
* With `A >= L(h_0 + h) + 2` the hypothesis `X >= C_L (8 p_j Norm zeta_j)^L` holds for large `M`.
* The restriction `Norm zeta_j <= M^h` is unnecessary for the "no slack" conclusion. If `Norm zeta_j < p_j` then
  `Norm zeta_j <= M^(h_0)`, and Theorem A applies with `h = h_0` once `A >= 2 L h_0 + 2`. Otherwise
  `|zeta_j|^2 >= p_j` trivially.
* So in regime (1.4) with `A >= 2 L h_0 + 2`, **every** exactly coherent certificate (2.1) has `|zeta_j|^2 >= p_j`.
  This is the unconditional half of L6 for exact coherence, as claimed.

*Scope note.* L6 was posed for the construction regime `p_j > X`. There "genuine size" means `|zeta|^2 > X`. In
regime (1.4) the conclusion `|zeta_j|^2 >= p_j` is much weaker in absolute terms (`p_j << X`). It still gives zero
slack in the certificate, which is all Corollary A.1 claims.

## 4. Theorem B (Hamming lock via quadratic large sieve)

**Proof re-derived; correct.**
* *Step 1.* `Ag_j(zeta)` is contained in `{q <= X, q = 7 (8), q not dividing n, (n/q) = 1}`. So
  `theta(Ag) <= (theta_7' + S_n)/2`. This gives `S_n >= (1+delta) theta_7 - theta_7' >= delta theta_7`. That is
  slightly better than the `delta theta_7 - log n` written in the proof.
* *Step 2.* `S_n = sigma_m (S_m - E)` with `|E| <= log s`. With the sharper Step 1,
  `|S_m| >= delta theta_7 - log D >= delta X/5 - log D >= delta X/6` when `delta X >= 30 log D`.
  * *Slip:* with the proof's own accounting (`- 2 log D`) the hypothesis would have to be `delta X >= 60 log D`.
    This is harmless.
* *Step 3.* The bound is `|E_bad| (delta X/6)^2 <= C (DX)^(eps/2) (D + X) sum a_n^2` with
  `sum a_n^2 <= X log X`, giving `|E_bad| <= 72 C X^eps log X/delta^2`. Correct for `D <= X`. The `n` are primes
  `= 7 mod 8` and the `m` are odd squarefree, as (3.1) requires.
* *Step 4.* An odd `m <= D` has at most `log_3 D` prime factors, and the `p_j` are distinct. Correct.

**Corollary B.1: correct.** With `eps = 1/(4A)`, `D = M^(h_0+h) <= X` and `delta = M^(-1/10)`, the exceptional set
is `O(M^(0.45) log^2 M) = o(M)`.

**Remarks.**
* *GRH remark (slip).* Under GRH, `|S_m| << X^(1/2) log^2(mX)`. Excluding `|S_m| >= delta X/6` therefore needs
  `X >> (log D)^4/delta^2` (up to `loglog` factors), not `X >= (log D)^(2+eps)`. The `(log)^2` threshold is the one
  for exact coherence (least prime in a Chebotarev class), not for a bias `delta`.
* *Twisted-coherence remark (out of scope).* The remark says twisted coherent labels "agree with `zeta_j` exactly
  on the set where the parities match, which by Theorem B has relative mass `1/2 + o(1)`".
  * Theorem B gives only an upper bound `(1+delta)/2`.
  * Def. 4.1 of `p3Dlower.md` needs `p'_j > X`, so `p_j Norm zeta_j > X >= D`, and Theorem B does not apply.
  * In regime (1.4) the statement can be repaired by using a sub-genuine `zeta_j` of norm `<= M^h`.
* *What Theorem B measures.* It measures agreement **in parity**. For sub-genuine `zeta` it bounds per-column
  agreement (Def. 3.1) by `(1+delta)/2`, which is correct. It says nothing about genuine-height `zeta`, for which
  parity agreement is automatic (see §7.2).

## 5. Theorems C, C' (piecewise covering)

**Lemma 4.2: correct.**
* *Duality.* The image of `Psi: G -> {+-1}^P` is the annihilator of the characters `prod_(c in S) x_c` that are
  trivial on the image. This is standard for finite abelian groups.
* *Relations.* `prod_S psi_c = (./prod_S m_c)` is principal on `(Z/K)^*` iff `prod_S m_c` is a square. Since
  `p || m_c` for every `c`, this forces `|S|` even.
* *Signs.* `(ab-1)/2 - (a-1)/2 - (b-1)/2 = (a-1)(b-1)/2` is even for odd `a, b`, and an odd square is `1 mod 8`.
  So `prod_S t_c = +1`.
* *Exact check (`c5_duality.py`).* Realisability of a sign vector is an `F_2` linear system in the quadratic
  characters at the odd primes `l | K`, and every character vector is realised by CRT. I solved that system for 3000
  random instances: `p` in `{5, ..., 73}`, `P <= 8` pieces, norms `a^2 + b^2 < 1800`.
  * 819 of the instances have a nontrivial relation.
  * Every instance is solvable.
  * In every instance an actual prime `q = 7 mod 8` with all `(n_c/q) = -1` was found below `2*10^6`; the largest
    least witness was 51,871. This last part is evidence.

**Theorem C: correct.** At the prime of Lemma 4.2, every piece has `(Norm zeta_c/q) = -(p_j/q)`, so by (F1) it has
the wrong parity, while `ell_j(q)` has the forced one. Hence `X < C_L (8K)^L`. (4.1) follows for large `M`.

**Theorem C': correct (GRH).**
* The counting identity is right: `prod_c (1 + t_c psi_c(q))/2`.
* Principal `psi_S` have `prod_S t_c = 1`.
* A non-principal `psi_S` twisted by a character mod 8 stays non-principal, because the moduli are coprime.
* The GRH bound `<< sqrt X log(XK)` gives `2^P >= c sqrt X/(log X log(XK))`, and hence (4.2) in regime (1.4).

**Remark on twisted coherence: partly vacuous.** The remark says "Theorem C gives `T(X) + 1 > (A-1)/(L(h_0+h))`
unconditionally, which is unbounded as `A` grows".
* `p3Dlower.md` Def. 4.1 requires `p'_j > X`, so the pieces `zeta_j`, `zeta_j rho_t` have norm `> X = 3M^A`, i.e.
  `h > A`. Then `(A-1)/(L(h_0 + h)) < 1/L`, and Theorem C says nothing.
* Only Theorem C' (GRH), whose `log K = O(log^2 X)` survives `h ~ A`, gives `T(X) >= (1/2 - o(1)) log_2 X`.
* The unconditional statement holds for a regime-(1.4) variant of twisted coherence with base `zeta_j` of norm
  `<= M^h`, `h` fixed. In regime (1.4) such a `zeta_j` is a monomial in genuine primes (§6.1, second bullet).

**Scope of §§2-4 (medium).** §1 states in bold, and the task summary repeats, that "no column is certified (exactly,
approximately in Hamming distance, or piecewise with few pieces) by Gaussian integers `zeta_j` with `p_j` not
dividing `Norm zeta_j`".
* As written this is **false**. The theorems all need `p_j Norm zeta <= X^(1/L)/C` (A, C) or `<= D <= X` (B).
* `c7_nonorm.py` gives explicit examples. A prime `p' = p mod 4 prod_(q <= X, q = 7 (8)) q` exists, with 24, 57 and
  111 digits for `(p, X) = (13, 200), (101, 500), (1009, 1000)`.
  * It satisfies `(p'/q) = (p/q)` at every prime `q = 7 mod 8` below `X`.
  * So the Gaussian prime `zeta` over `p'` is sub-genuine and **exactly** parity-coherent with column `p` on all of
    `Q_7 ∩ [1, X]`.
  * In regime (1.4) such a `zeta` has `|zeta|^2 > X >> p_j`, so its certificate has no slack anyway.
* The correct dichotomy in regime (1.4) is therefore:

> a certificate (exact, Hamming above `(1+delta)/2` for all but `o(M)` columns, or with fewer than
> `(A-1)/(L(h_0+h))` pieces) by Gaussian integers of norm `<= M^h` has `v_(p_j)(Norm zeta)` odd; certificates of
> norm `> M^h` have `|zeta|^2 > p_j` once `h >= h_0`. In both cases `|zeta_j|^2 >= p_j`.

Under this wording the conclusion "certified implies genuine size" stands. The phrase "certification by sub-genuine
Gaussian integers is impossible" needs the norm bound.

## 6. Proposition D and corollaries

### 6.1 What is correct

**Lemma 5.1: correct.**
* For `b = 1` the targets are the `N/4` orbits of odd residues. For `b = 0` they are `N/4 - 1` orbits.
* There are at most `N/8 + 1` source orbits.
* The case analysis for `s = +-1` is right.
* `c1_lemma51.py` builds the map exactly as in the proof and verifies the property exhaustively for every
  `N = 0 mod 4`, `16 <= N <= 1200`, `b in {0, 1}`.
* *Cosmetic.* `N = 12` (`q = 47`) also admits such a map for both parities, so `47` need not be in `S_0`.

**Proposition D: correct.** Each case was re-derived.
* *Case 1.* `iota = log i` is odd when `m_q = 4 mod 8`.
* *Case 2.* The parity holds by the choice of `p'_j` on `S_0`.
* *Case 3.* The 2-part is a direct factor `Z/8`. `mu_4 = {0,2,4,6} x {0}`, so a collision forces `u^8 = 1`.
* *Case 4.*
  * Different parities never collide at `Q_1`.
  * Equal parities use the same `iota_b`, and Lemma 5.1 gives `x_j = s x_l`. So `pi_1(<v,L>)` lies in the subgroup
    of order 16 and `pi_q(<v,L>) = 0`.
  * `+-e_j` and `+-2e_j` never collide.
* *Integrality.* `u^16 = 1` in `T_(q^a)` iff `q^a | Im(W^16)`. `Im(W^16) != 0` because the `p'_j` are distinct and
  split.
* *Moduli dividing N.* The free values at `q | N` do not matter: `sharp.md` (e-all) takes `m_(v,s)` over
  `q not dividing N` only.

**Exact reproduction (`c2_propD.py`).** This uses the write-up's parameters: `p_j in (10^4, 10^6)`, `p'_j` the next
prime `= p_j mod 4*17*31*47`, `r = 12`, `X = 400, 1000, 2000, 3000`.
* The case counts at `X = 3000` match (219, 3, 105, 102).
* Every label has the forced parity.
* Norm 1 and norm 2: 0 failures of `m_v | Im(Z_v^16)`.
* Norm 3: 202, 288, 371, 446 failures out of 1156.
* Norm 3 restricted to the moduli of cases 1-3: 0 failures.

**Corollary D.1: correct.**
* The arithmetic: `8hk log M <= k(tau/2 - (3/2) log M)` iff `tau >= (16h + 3) log M`.
* The comparison with independent labels is `sharp.md` Prop. 6.5. That proof already works on the window of
  `q = +-3 mod 8`, where `ell_j - ell_l` is uniform on its parity class, and disjoint pairs are independent.

### 6.2 Corollary D.2: scope error (major, part of M1)

D.2 claims that the system is uncertified "in regime (1.4) with `p'_j | N` ..., **or in the setting of
Corollary D.1**", by Theorems B and C.
* In the setting of D.1, `p_j >= M^(16h+3) > X`.
* Theorem B needs `p_j Norm zeta <= D <= X`. So no `zeta` qualifies, and the statement is vacuous.
* Theorem C needs `X >= C_L (8 prod_c p_j Norm zeta_c)^L`, which is false.
* So in the regime where D has slack, the write-up has no unconditional proof that D is uncertified.
* Under GRH, `p3Dlower.md` Thm 3.2 (which assumes `p_j > X`) excludes exact coherence with slack, and a GRH Hamming
  analogue of Theorem B is routine. Unconditionally the question is open, because sub-genuine `zeta` of norm
  `<= X^(O(1))` with `p_j Norm zeta` a pseudosquare on `Q_7 ∩ [1, X]` are not excluded.

### 6.3 Construction D in regime (1.4): the height condition fails (major, part of M1)

In regime (1.4) the auxiliary Gaussian integers must be units at every `q <= X` that does not divide `N`. So
`Norm zeta_j` either exceeds `X` (much larger than `p_j`) or is a product of genuine primes (§6.1 of the write-up).
D.2's choice `zeta_j = pi_(sigma(j))` is of the second kind.
* The certificate `|Im(Z_v^16)| <= e^(8 w(sigma v))` then exceeds the per-unit budget `e^(w(v)/2)` by a factor 16
  in the exponent on average.
* `c3_regime14.py` shows that this is not just slack in the bound: the actual demand exceeds the budget.

The test: profile primes `p_j = 1 mod 4` (not 17) of size `<= 400`, and `X = 2*10^4` or `5*10^4`, so
`X >= p_max^(1.65)`. `sigma` cycles each Legendre-pattern class on `S_0`. Exact per-unit moduli `m_(v,s)` were
compared with `nu_s^(-1) e^(w(v)/2)` as integers. "max ratio" is `max log m_(v,s)/budget`.

| instance | system | norm-1 monomials violating | norm-2 monomials violating | max ratio |
|---|---|---|---|---|
| `X = 20000`, `r = 12`, `p in [73, 397]` | genuine `G` | 0/12 | 0/144 | 1.00 |
| | D(`sigma = id`) | 2/12 | 124/144 | 6.55 |
| | D(`sigma`), 10 of 12 columns moved | 10/12 | 134/144 | 8.15 |
| `X = 50000`, `r = 10`, `p in [29, 277]` | genuine `G` | 0/10 | 0/100 | 1.00 |
| | D(`sigma = id`) | 4/10 | 93/100 | 8.03 |
| | D(`sigma`), 7 of 10 moved | 8/10 | 93/100 | 10.13 |

These are exact computations on explicit instances. They prove that Construction D, with the parameters of D.2,
does not in general satisfy the norm-`<= 2` height condition in regime (1.4).
* The mechanism is visible in the proof of Proposition D.
  * At `b_q = 3` the column `e_j` collides wherever `u^8 = 1` and `nu_j(q) = 0`. This brings in the prime factors of
    `Re(zeta_j^4)` (size `p^2`), and all of them lie below `X` in regime (1.4).
  * At `b_q >= 4`, pairs collide wherever `L_j = +-L_l mod N_q/4`, which is a 16-torsion condition.
* So in regime (1.4), Construction D is only `beta`-good with `beta = 8 h_0`. Members have `beta <= h_0/2`
  (since `tau_j <= h_0 log M`).
* Proposition D therefore says nothing about what an inverse theorem in regime (1.4) must use.

### 6.4 "Uncertified" in the construction regime is a matter of definition (medium)

In the construction regime D is built from `zeta_j` of norm `<= M^h << p_j`.
* It is system-level coherent with them, in the sense of (F1) with `u = 1`, at every `q in Q_3 ∪ S_0`, and at the
  `b_q = 3` moduli where the parities match.
* Its norm-`<= 2` demand is bounded by the Gaussian certificate `Im(Z_v^16)`.

Round 7 phrases the conclusion of L4 as "induced by Gaussian integers `zeta_j` of controlled height". D satisfies
that conclusion, through a fixed non-linear encoding at each `q`.
* The write-up counts D as uncertified because it uses the narrower notion "coherent up to (F1) on most of `Q_7`".
* That is a legitimate choice, but it should be stated as a choice. Under the broader notion, which is the one a
  route to L6/L7 can use (demand bounded by a fixed power of a Gaussian certificate), D at norm `<= 2` is a
  certified system and not a counterexample.

## 7. Section 6 claims (where L4 sits)

### 7.1 §6.1, second bullet: correct

* A Gaussian integer of norm `<= M^h < X` that is a unit at every `q <= X` not dividing `N` has all its prime factors
  above 2 or above the `p_l`. A prime factor above an inert `l <= X`, or above a split `l <= X` not dividing `N`,
  would make it a non-unit.
* So its labels are `B G` modulo `mu_4`. Theorem A applied to each row forces `B = I mod 2`.
* Then `||B^T v||_1 >= #{j : v_j odd}`, so there is no slack at `e_j`.

### 7.2 (I1\*) is too weak to be "the precise remaining inverse step" (major, M2)

(I1\*) asks for agreement with a genuine-height coherent system on `> 1/2 + delta` of the `Q_7` log-mass, for a
positive proportion of columns. Two problems.

**(a) Calibration.** Agreement in the sense of Def. 3.1 at `q` (some unit `u`, modulo `mu_4`) holds iff
`gcd(ell_j(q), n) = gcd(L(q), n)` with `n = m_Q/4`. Units of `Z/m_Q` surject onto units of `Z/n`, and in `Z/n`
associates are exactly the elements with equal gcd.
* For a genuine-height `zeta`, parity agreement is automatic. So the natural baseline is not 1/2.
* Take `x` uniform in the parity class of `Z/m_Q`, and `d = gcd(G_j(q), n)`. Then
  `P(agree) = phi(n/d)/(n/2)` exactly.
* `c4b_baseline.py` sums this over all 4511 primes `q = 7 mod 8` below `2*10^5`, log-weighted, against the genuine
  labels of 40 random columns `p_j in [10^3, 10^6]`.
* **The expected agreement of a uniformly random parity-respecting column with the genuine column is between 0.540
  and 0.555 (median 0.548).** The value per modulus is exact; the aggregate is a floating-point sum.
* `c3_regime14.py` agrees with this: sampled random labels give a median of 0.545-0.554. Construction D(`sigma`)
  gives 0.556-0.564, and D(`id`) gives 0.70.
* So (I1\*) with `delta < 0.04` is satisfied, in expectation, by structureless labels against `G` itself. For larger
  `delta` the threshold is arbitrary.
* The `1/2` in (I1\*) is inherited from Theorem B, where it is the right baseline **for sub-genuine `zeta`**. There,
  parity agreement is about 1/2 and per-column agreement is at most parity agreement. It is the wrong baseline for
  the genuine-height side, which is the side (I1\*) asserts.
* This is evidence about the calibration, not a theorem that (I1\*) is vacuous: an adversarial system can avoid the
  gcd class of `G_j` at every `q`. But it shows that (I1\*) as stated cannot carry the information an inverse
  theorem has to deliver.

**(b) Transfer.** In Def. 3.1 the unit `u` may depend on the column. So agreement on a set `Q'` of moduli does not
make the collision sets of multi-column monomials on `Q'` equal to those of the certificate. Only a common `u`
(F1) does that.
* (I3)/L7 is a statement about collisions of `D_s(Z_v)` for multi-column `v`. So "(I1\*) plus Theorem B ... what
  then remains is (I3) for genuine-height certificates" skips a step.
* Either (I1\*) must be strengthened to **system-level** agreement, with a common unit per modulus, on a set of
  moduli large enough to carry L7-type tails. Or a separate transfer lemma is needed.
* With system-level agreement the random baseline is negligible. Theorem B still applies, since system-level
  agreement implies per-column agreement. But the threshold that Theorem B allows is still `(1+delta)/2`. So the
  strengthened (I1\*) needs system-level agreement on more than half of `Q_7`, which is a much stronger demand than
  the write-up suggests.

**(c) "Equivalent".** §0(b) and §6.4(b) say that the inverse step "is equivalent to (I1\*)". Only one direction is
argued: (I1) in the weak per-column Hamming sense, plus Theorem B, gives (I1\*). Nothing argues that (I1\*), plus
L7, gives growth, and (b) shows that it does not as stated.

### 7.3 §6.2: the "must" statements hold only for regime-independent strategies (major, M1)

* *"Norm `<= 2` is free ... So (I1\*) restricted to `K = 2` is false."* Not established.
  * The slack is proved only for `p_j > X`. In regime (1.4), where (I1\*) lives, Construction D fails the height
    condition (§6.3).
  * Refuting (I1\*) would also require D to be far from **every** genuine-height coherent system. Theorems A-C
    concern sub-genuine `zeta` only, and the write-up never addresses this.
  * In the instances of §6.3, D(`sigma`) agrees with the genuine `G` on 0.50-0.71 of `Q_7`, median about 0.56. So
    for small `delta` it satisfies the conclusion of (I1\*).
* *"Shallow 2-adic classes are free for every norm."* (6.1) is a correct certificate bound,
  `2^(B+1) log|Z_v|`. It is "free" only if `2^(B+1) log|Z_v| <= w(v)/2`, i.e. (column by column) in the
  construction regime with `|zeta_j|^2 <= p_j^(2^(-B-1))`.
  * In regime (1.4) an auxiliary `zeta` of norm `< X` is a monomial `B G` in the genuine primes. If `B = I mod 2` it
    has the right parity and needs no repair.
  * Otherwise the repaired certificate is `2^(B+1) w(B^T v)/2`. For a permutation matrix `B` this exceeds the budget
    `w(v)/2` by the factor `2^(B+1)` on average over the columns.
  * Auxiliary `zeta` of norm `> X` exceed the budget by a factor `> A/h_0` already before repair.
  * So nothing shows that obstructions in regime (1.4) live only on the deep moduli.
* *Mass (slip).* `{q : b_q > B} = {q = +-1 mod 2^(B+1)}` has relative density `2/2^B = 2^(1-B)`, not `2^(-B)`.

### 7.4 §6.3: heuristic content (minor)

* The section is titled "Why norm 3 is different, precisely". Its content is heuristic and should be labelled so:
  "any `E` linear in `L` turns (6.2) into the twisted shape", "any `E` non-linear in `L` loses the multiplicative
  certificate", "no intermediate mechanism is known".
* Identity (6.2) itself is correct.
* *"Multiplicatively independent modulo squares"* should read "modulo squares and 2". On `Q_7`, `(2/q) = 1`, so for
  example `n_1 = 2 n_2` gives `delta_1 = delta_2`, and not all `2^s` sign vectors occur. The proof via Lemma 4.2 uses
  only the odd squarefree parts `m_j`, and needs those to be independent.

### 7.5 §5, remark on additive labels (minor)

* The divisor count is heuristic, and is marked so only in parentheses.
* *Gap in the pigeonhole.* For additive labels a pair collides at `q in Q_1` iff
  `nu_j(q) = nu_l(q)` and `(q - chi)/8 | a_j - a_l`. So `a_j = a_l mod L'` gives collisions only at those `q | ...`
  where the two Legendre symbols agree. That is about half of them for a generic pair, not "all of them".
* The conclusion survives at half strength. The remark should be labelled heuristic. The proved statement about
  additive labels remains `p3Dlower.md` Prop. 3.1.

### 7.6 §7 evidence: reproduced independently

* *Lemma 2.1.* For 3000 sub-genuine `n = p Norm zeta`, `p in (10^5, 10^7)`, the least `q = 7 mod 8` with
  `(n/q) = -1` is at most 239 (write-up: 191 on its own sample).
* *Lemma 4.2.* See §5.
* *Theorem B.* For 30 columns and 652 norms, the maximal `S_n/theta_7(10^6)` is 0.0247 (write-up: 0.0241).

All of this is evidence, and correctly labelled as such in the write-up.

## 8. Independent computations

All in `round9/referee_L4_checks/`, written from scratch.

| script | output | what it does | status |
|---|---|---|---|
| `glib.py`, `t_glib.py` | `t_glib_out.txt` | Gaussian arithmetic, the torus `T_(q^a)`, Pohlig-Hellman logs (asserted); tests the parity fact `log(w/conj w) = [(Norm w/q) = -1] mod 2` | exact, OK |
| `c1_lemma51.py` | `c1_lemma51_out.txt` | Lemma 5.1 map built and verified for all `N = 0 (4)`, `16 <= N <= 1200` | exact, ALL OK |
| `consD.py`, `c2_propD.py` | `c2_propD_out.txt` | Construction D from the text, construction regime; parity; `m_v \| Im(Z_v^16)` at norms 1-3 | exact: norm `<= 2` 0 failures; norm 3: 202-446 of 1156 fail; cases 1-3 never fail |
| `c3_regime14.py` | `c3_regime14_out.txt` | regime-(1.4)-type instances: per-unit height condition for genuine, D(id) and D(`sigma`); per-column agreement with `G` on `Q_7` | exact per instance (§6.3, §7.2) |
| `c4b_baseline.py` | `c4b_baseline_out.txt` | expected per-column agreement of random parity-respecting labels with genuine labels on `Q_7 ∩ [1, 2*10^5]` | exact per modulus; 0.540-0.555 |
| `c5_duality.py` | `c5_duality_out.txt` | Lemma 4.2 sign system as an `F_2` system (3000 instances, 819 with relations); witnesses; Lemma 2.1 least non-residue | exact (duality); evidence (witness sizes) |
| `c6_bias.py` | `c6_bias_out.txt` | Theorem B bias `S_n/theta_7(10^6)` | evidence; max 0.0247 |
| `c7_nonorm.py` | `c7_nonorm_out.txt` | sub-genuine `zeta` parity-exact on all of `Q_7 ∩ [1, X]` (norm bound necessary) | exact |

Not done: I did not search for a construction that is norm-`<= 2`-good in regime (1.4). A column analogue of the
`D_A` design of `sharp.md` (classes `j mod p'(q)` inside parity classes at small `q`, distinct `+-`-orbits at
`q > 16 r`) is a natural candidate. Whether its `O(log M)` constant fits inside the per-unit budget `w(v)/2` of
near-extremal profiles is open (§9, item 3).

## 9. Corrected statement

**L4 is undecided.** Nothing here bounds actual circles beyond `M <= (2/3 + eps) log R/loglog R`. All statements
are about label systems (designs) in P3-all\*; actual clusters satisfy them with `zeta_j = pi_j`.

1. **Reduction (proved, actual circles).** By `stability.md` Cor. 3.6(1), growth follows from (T'):
   for every `c_0` there is `A(c_0)` such that, for large `M`, no member of `F^all_(3M^A)(M, C)` with all
   `p_j <= M^(2c_0 - 3 + o(1))` has `W <= c_0 M log M`.
2. **Exact coherence (proved, unconditional given Linnik; designs).** In regime (1.4) with `A >= 2 L h_0 + 2`,
   every column whose labels agree, at all `q = 7 mod 8` up to `X`, with `u_q log(zeta/conj zeta)` modulo `mu_4`
   (where `u_q` is a unit that may depend on the column) has `|zeta_j|^2 >= p_j`.
   * If `Norm zeta <= M^h` with `A >= L(h_0 + h) + 2`, then `v_(p_j)(Norm zeta)` is odd.
   * Coherent certificates therefore have zero slack. This is the unconditional half of L6 for exact coherence.
3. **Hamming coherence (proved given Heath-Brown's quadratic large sieve; designs).** Fix `delta = M^(-1/10)` and
   `A >= h_0 + h`. For every parity-respecting system with primes (1.4), all but `O(M^(0.45) log^2 M)` columns agree
   **in parity** on at most `(1+delta)/2` of the `Q_7` log-mass with every `zeta` such that `Norm zeta <= M^h` and
   `v_(p_j)(Norm zeta)` is even. The norm bound is necessary (§5, `c7`).
4. **Piecewise coherence (proved; designs).** Covering a column by `P` pieces of norm `<= M^h` with
   `v_(p_j)(Norm)` even needs `P > (A-1)/(L(h_0+h))` unconditionally, and `P >= (1/2 - o(1)) log_2 X` under GRH.
   For `p3Dlower.md` Def. 4.1 (norms `> X`) only the GRH statement applies.
5. **Proposition D (proved, any regime; designs).** Construction D is parity-respecting, and `m_v | Im(Z_v^16)`
   for every `v` with `||v||_1 <= 2`.
   * *In the construction regime* (`p_j >= M^(16h+3)`, `p'_j in (X, M^h]`) this gives the norm-`<= 2` height
     condition with slack `(3/2) log M` per unit, where independent labels fail.
     * There D is coherent, in the system-level (F1) sense, with the `zeta_j` of height `M^h` on `Q_3 ∪ S_0` and on
       part of the `b_q = 3` moduli.
     * It is not shown to be uncertified unconditionally: Theorems A-C do not apply when `p_j > X`.
     * It is certified in the generalised sense of a 16th-power Gaussian certificate.
   * *In regime (1.4)* D is uncertified in the sense of items 2-4. But it does **not** satisfy the norm-`<= 2`
     height condition in general: the exact instances of §6.3 have 93-134 of the 100-144 norm-2 monomials in
     violation, against 0 for the genuine labels.
   * Consequently:
     * L4 restricted to norm `<= 2` is refuted only for the narrow notion "certified = (F1)-coherent on most of
       `Q_7`", and only in the construction regime, where the uncertification is itself unproved without GRH.
     * In regime (1.4), the only regime growth needs, it is **open** whether norm-`<= 2` data, or the shallow
       moduli `b_q <= B`, can drive an inverse theorem.
     * The statements "any inverse theorem must use norm `>= 3`" and "must use only the moduli
       `q = +-1 mod 2^B`, `B -> infinity`" hold only for inverse theorems that must also cover `p_j > X`.
6. **(I1\*) (reformulation, not a reduction).** Formally, (I1) in the per-column Hamming sense plus items 2-3 is
   equivalent to (I1\*). As stated, (I1\*) is not a usable inverse step, for two reasons.
   * Random parity-respecting columns already agree with the genuine column on 0.540-0.555 of the `Q_7` log-mass in
     expectation (exact baseline, `c4b`), so the threshold `1/2 + delta` is mis-calibrated.
   * Per-column agreement does not transfer the collision structure of multi-column monomials, which (I3)/L7
     needs.
   * A usable version needs system-level agreement, with one unit per modulus for all columns. Given item 3 it must
     hold on more than half of `Q_7`. Formulating that version is left open.
7. **Unchanged.** The bracket `3 M log M - O(M) <= W_min^all <= (66 A + o(1)) M log^2 M/loglog M` is unchanged.
   No growth improvement and no new cap for P3-all\* is claimed or follows.

**Suggested open items** (replacing §8 of the write-up):
1. In regime (1.4), find a parity-respecting system that is uncertified (items 2-4) and meets the per-unit height
   condition for all `||v||_1 <= 2`. Or show that none exists. This is the regime-(1.4) version of Proposition D.
2. Formulate (I1\*) with system-level agreement and a threshold calibrated against random labels. Then prove or
   refute it at norm 3.
3. Test the column analogue of the `sharp.md` `D_A` design at norm `<= 2` in regime (1.4) (§8, last paragraph).

## 10. Issue list (severity)

| # | severity | issue | where |
|---|---|---|---|
| 1 | major | Prop. D's slack (construction regime) and its uncertification (regime (1.4)) are proved in disjoint regimes; in regime (1.4) D violates the height condition (exact instances); "L4 false up to norm 2" and "inverse theorem must use norm `>= 3` / deep moduli" are unproved for the regime growth needs | §0, §5 Cor D.2, §6.2, §6.4 |
| 2 | major | (I1\*) is mis-calibrated (random baseline 0.540-0.555 > 1/2 against genuine labels) and per-column agreement does not transfer multi-column collisions; "what remains is (I3)" skips a step | §6.1, §6.4(b) |
| 3 | medium | Cor. D.2 "or in the setting of Cor. D.1": Theorems B and C are vacuous or false there (`p_j > X`) | §5 |
| 4 | medium | bold claims "no column is certified by sub-genuine `zeta`" omit the necessary norm bound (explicit counterexample) | §1, §0, summary |
| 5 | medium | "uncertified" for D in the construction regime depends on an unflagged narrow definition; D is certified by `Im(Z_v^16)` with height `M^h << p_j` | §5, §6.4(a) |
| 6 | medium | twisted-coherence remark: Theorem C is vacuous for Def. 4.1 (`h > A`); only C' (GRH) applies; Thm B remark likewise out of scope | §3, §4 remarks |
| 7 | minor | "genuine height" should be `v_(p_j)(Norm zeta)` odd (or `|zeta|^2 >= p_j`), not `p_j | Norm zeta` | §1, §6.1 |
| 8 | minor | Theorem B step 2 constant (30 vs 60 `log D`; harmless); GRH remark needs `X >> (log D)^4/delta^2` | §3 |
| 9 | minor | density of `Q^(B)` is `2^(1-B)` | §6.2 |
| 10 | minor | §6.3 heuristic content labelled "precisely"; independence must be modulo squares **and 2** | §6.3 |
| 11 | minor | additive-label remark: heuristic count, and collisions only where Legendre symbols agree | §5 |
| 12 | minor | `47` need not be in `S_0` (Lemma 5.1 holds at `N = 12`) | §5 |
| 13 | positive | Cor. A.1 holds for all `zeta` once `A >= 2 L h_0 + 2` | §2 |
