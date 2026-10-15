# Referee report on `round14/c1construct.md` (C1: a fake circle for P3-all\* in regime (1.4)?)

Referee: adversarial, independent computation. Scripts and outputs are in `referee_c1construct_checks/`. None of
them imports the author's code.

## 0. Verdict

**The negative headline survives: no C1 fake, no lower bound, growth status unchanged.** The individual theorems
(B1-B3, P0, P, P1, P2, L, G, C, R, S, S1) are correct as stated, with minor corrections to B1's constant and B3's
"norm `K + 1`" label.

**One MAJOR error invalidates several of the claimed consequences.** Condition (i) of Prop. P (no `p_l` divides
the pinning numbers `D_s(pi_j)`, `s in S_j`) is presented as "rare (Cor. P1)". It is not:
* The evidence for rarity is vacuous by construction. Primes drawn from one decade `[10^k, 10^(k+1)]` cannot
  divide numbers `< 10^k`.
* With the smallest split primes, (i) fails for 80-91% of columns (exact, `ref_cond_i.py`).
* A builder who selects primes can make it fail in all but one column.
* Even with **genuine labels at every modulus** (P3-G) and all four units used, 50-77% of the columns of the
  first `r` split primes are not pinned at norm 1 (`ref_pin_full.py`). That includes 2 of the 4 columns of the
  round-10 flexible-top instance `{5, 13, 17, 29}`.

Hence the following claims do not hold as stated:
* "P3-G collapses at norm 1 in regime (1.4)".
* "The round-10 flexible tops cannot seed a fake on top of a genuine core".
* "Rig_pin + Prop. P would make P3-all\* exact in regime (1.4)".
* "The label-angle step is settled at norm 1".

All four hold only for prime sets satisfying (i), a condition the fake controls.

**A second, smaller overclaim.** "Every known coordination mechanism (incl. Gaussian coherence with any repair)
is excluded by N2-N6" is unproved. Root-repaired coherence on the deep classes is uncertified, not excluded.

## 1. Class definition and regime check

* **Class: correct.** §1 reproduces `sharp.md` Def. 1.1 (P3-all\* `= F^all_X`) exactly:
  * factored residues at every odd `q <= X`, `q != p_j`, with `A_q = max{a : q^a <= X}`;
  * (G), (e-char), and the slab conditions (1.2) for **every** `v != 0` and every unit;
  * `m_(v,s)` over `q` not dividing `N`;
  * parity `nu_j(q)` forced by Prop. 2.1/Cor. 2.2;
  * per-unit budget `B(v,s) = w(v)/2 + [s odd](log 2)/2`, with `w(v) = sum |v_j| log p_j`.

  The proofs of P, G and C use (1.2) itself, with angles, which the class contains. They do not use any
  stronger condition. P3-G is used in round 7's definition (`r_j(q) = pi_j mod q^(A_q)` with `pi_j` fixed), and
  this matches Prop. P(ii)'s single `(eps_j, t_j)`. Actual clusters are members for every `X` (`sharp.md`
  Prop. 1.2).
* **Regime: correct.**
  * Prop. P and Cor. P1 need `X >= a_j + b_j`. Since `a_j + b_j <= sqrt(2 p_j) <= sqrt 2 M^(h_0/2)`, this holds in
    regime (1.4) once `A > h_0/2`.
  * Lemma L and Thm C need `A >= L(h_0 + h) + 2`. This is chosen after `c_0`, which is legitimate.
  * Thm S is `X`-uniform. Thm R needs `A >= 4`.
  * B1 needs `A > h_0/2`. B2 needs `A > h_0`. B3 needs `A > 1`.
  * No claim uses moduli `q | N` in (1.2).
* **Mod-`p` directions.** The only mod-`q` inferences are "collision at `q^a` ⇒ `q^a | D_s`" (or
  `| Im(gamma^(2^k))`). These are used as upper bounds on `m`. The reverse direction is used only for genuine or
  exactly coherent labels, where it is an equivalence. Both directions are used correctly.

## 2. Break theorems B1-B3

**Thm B1 (random fibre fails at norm 1): survives, with a corrected constant (minor).** The mechanism is right:
at `q = +-3 mod 8`, `m_q = 4 mod 8`, `iota = m/4` is odd, the 4-torsion `{0, m/4, m/2, 3m/4}` has two elements of
each parity, so a label uniform on its parity class is torsion with probability exactly `4/m_q`. A torsion label
is a level-1 collision of `e_j` with some unit. A collision at `q > sqrt(2 p_j)` exceeds `B(e_j, s) <= (1/2) log(2 p_j)`.
Independence over `q` (sharp.md Lemma 6.3(2)) gives the product. Two gaps:
* (B1-a) The product must run over `q` **not dividing `N`**, i.e. it skips `q = p_l = 5 mod 8`. The displayed
  inequality `prod <= (1 + O(1/log M)) (log sqrt(2p_j)/log X)^2` is then false when the prime set contains small
  primes. Exact (`ref_B1.py`): `M = 30`, `A = 3`, `X = 81000`, prime set containing the 90 smallest split primes.
  The product is inflated by factors 2.32, 1.84, 1.59, 1.28 for `p_j = 101, 1009, 5009, 26993`.
  * Correct form: the skipped `q` are at most `rho M` primes `= 5 mod 8` above `Y = sqrt(2p_j)`, so the inflation is
    at most `(1 + o(1)) max(1, log(M log M)/log Y)`.
  * Hence `P(survive) <= (1 + o(1)) log Y max(log Y, log M)/(log X)^2 <= (1 + o(1)) (h_0/2) max(h_0/2, 1)/A^2`.
  * This still tends to 0 as `A -> infinity`, so the qualitative statement (failing fraction `>= 1 - O(1/A^2)`)
    survives. The constant `((h_0/2)/A)^2` and the table column "`(log sqrt(2p)/log X)^2`" do not survive for
    such prime sets.
* (B1-b) `1 + O(1/log M)` needs `p_j -> infinity`. For bounded `p_j` the Mertens error is `O(1)` relative, but the
  bound is then `O(1/log M)`, which is harmless.

**Thm B2 (torsion-avoiding fibre fails at norm 2): survives.** Re-derived:
* Given `ell_k` (not torsion), the four collision values `ell_k + T_4` contain exactly two of parity `nu_j(q)`.
* Neither of these two is torsion, so the conditional probability is `2/(m/2 - 2) >= 4/m_q`.
* `Y = sqrt(2 p_j p_k) <= sqrt 2 M^(h_0) < X` needs `A > h_0`.
* The same `q | N` correction as (B1-a) applies.

**Thm B3 (L4K Thm 2 systems fail at norm 4, deterministically): proof correct. The label "norm `K + 1`" is
wrong for `K = 2` (minor).**
* `<v, a> = 2(K+1)(sigma-sum difference) + (nu-sum difference)`, and `|nu-difference| <= 2 < 2(K+1)`. Pairs with
  equal `sigma`-sum are disjoint. Convexity over the three `nu`-classes gives `Z_r`. All of this re-derived.
* Independent exact check (`ref_B3.py`):
  * `Z_r` reproduced: 1.667 at `r = 10`, 44 at `r = 16`, 223.3 at `r = 24`. `|F| = r(r-1)(r-2)(r-3)/8`.
  * 400 random `(sigma, nu)` per `r`: the number of vanishing `v in F` was always `>= Z_r`. Minima 5, 13, 81, 322 for
    `r = 8, 10, 16, 24`. At `r = 8`, `Z_8 = -2`, so the bound is vacuous there.
  * An instance of my own: `K = 3`, `r = 10`, ten random split primes in `[1000, 6000]`, `sigma_q` random and
    `nu` the actual Legendre parities at every `q in (59, 60000]`. **630/630** monomials of `F` exceed their unit-0
    budget (min demand 1434.7, max budget 16.95).
* `Z_r/|F| = (2/(9r))(1 + O(1/r))` is correct asymptotically, but the `O(1/r)` constant is large: the ratio to
  `2/(9r)` is 0.12, 0.58, 0.76 at `r = 10, 16, 24`.
* For `K = 2`, the vanishing monomials have norm 4 = `K + 2`, not `K + 1`. Norm-3 monomials keep the `2 Lambda` term
  and do not vanish identically. So the summary line "fail at norm `K + 1`" and the §2.6 entry "`K + 1 = 3, 4`"
  should read "norm 4 (= `K + 1` for `K = 3`, `K + 2` for `K = 2`)". Only `h = 1` and `beta_j = j` are proved, not
  general `B_h` data.
* B3 refutes `L4K.md` Remark 2. It is not a statement about fakes, since those systems were never proposed as fakes.

## 3. Norm-1 pinning (Lemma P0, Prop. P, Cor. P1-P2)

**Lemma P0: correct.**
* Re-derived by hand for both cases. For `theta < pi/4`, `s = 2` gives `x in [-theta, theta]` and `s = 3` gives
  `x >= theta`. For `theta > pi/4`, `s = 2, 1, 3` leave `{theta, pi/2 - theta}`.
* Independent arc-union check (`ref_P0.py`, 50-digit `mpmath`) for all 4783 split primes `< 10^5`: the arcs for
  `s in S(theta)` leave no free set of positive length.
* The identity `|sin(theta - s pi/4)| = nu_s |D_s(pi)|/sqrt p` was re-derived.

**Prop. P as a theorem: correct, with its hypotheses (i)-(ii).** Steps checked:
* `m_(e_j, u) >= |D_s|`: every exactly-dividing prime power is odd, prime to `N` by (i), and `<= a_j + b_j <= X`,
  so it lies within `A_q(X)`.
* `nu_u = nu_s`, since `u = +-s mod 2`.
* `|sin(phi - u pi/4)| = |sin(x - s pi/4)|` for `x = eps(phi - t pi/2)`.
* Unwinding to `xi_j in {i^t pi, i^(t+1) conj pi, i^t conj pi, i^(t-1) pi}`.
* `R e^(i <2a_x - e, phi>) = +- prod xi^a conj(xi)^(e - a)`. The `mod pi` ambiguity gives a common sign.
* Distinctness by unique factorisation, and arc length `C sqrt R`.
* (ii) uses one `(eps_j, t_j)` for all pinning moduli of column `j`. That is right for P3-G as defined in round 7
  (`r_j(q) = pi_j mod q^(A_q)`, `pi_j` fixed), so **Cor. P1 is correct as stated, including its hypothesis (i)**.
  As the author says, it is the norm-1 case of `p3G.md` Thm 3.1. The new content is the observation that
  `X >= a_j + b_j` makes the norm-1 relations tight.

**MAJOR: hypothesis (i) is not "rare". It fails for most columns of the prime sets that matter, and the
claims built on its rarity do not hold.**
* The only evidence offered for rarity (`pin_check.py` (c)) is **vacuous by construction**. For primes in
  `[10^k, 10^(k+1)]` with `k >= 2`, every pinning number is `<= sqrt(2 * 10^(k+1)) < 10^k`, so no `p_l` of the set
  can divide it. The "0.00 failing columns" there is automatic.
* Exact (`ref_cond_i.py`), prime set = the first `r` split primes (the smallest budgets, the sets used in rounds
  9-10 and the natural profiles of near-extremal clusters): **(i) fails for 8/10, 16/20, 42/50, 87/100, 180/200,
  448/500, 913/1000 columns.**
* A fake chooses its own primes. About 70% of all split primes `p <= 10^6` have `5` or `13` dividing
  `prod_(s in S_j) D_s(pi_j)`. A set containing `{5, 13}` and `r - 2` such primes fails (i) in **49/50, 199/200,
  999/1000** columns. Regime (1.4) allows such sets: all `p_j <= M^(h_0)`, nothing bounds them from below.
* (i) is only a sufficient condition. Using all four units, the genuine collisions also pin through `{1, 3}`. So I
  tested the strongest form available: genuine labels at **every** odd `q <= X`, i.e. P3-G (`ref_pin_full.py`,
  arc-union in 50 digits).
  * With the first `r` split primes, **5/10, 13/20, 34/50, 74/100, 153/200, 388/500 columns are not pinned at
    norm 1**. The free set is an interval of length about 0.3 around the genuine angle.
  * For the round-10 flexible-top instance `{5, 13, 17, 29}`, columns 13 and 17 are unpinned at norm 1. For 13:
    `D_3(3 + 2i) = 5`, and the free set is `[0.281, 0.588] ∪ [2.554, 2.862]`, total 0.307.
  * The random sets in `[10^4, 10^5]` give 0/200, as the author found.

Consequences:
* **"P3-G collapses at norm 1 in regime (1.4)"** (summary, §0 item 2) is false as a statement about the regime. It
  holds only for prime sets satisfying (i), and P3-G with the smallest primes is **not** collapsed at norm 1.
  Collapse there would need the long closures of `p3G.md` (open for selected sets, round 7).
* **§2.5 and the §2.6 row "genuine core + flexible top: no fake at all, pinned at norm 1, Prop. P"** are false for
  the instance cited. Prop. P does not apply to columns 13, 17, 29 of `{5, 13, 17, 29}`: (i) fails, since 5 divides
  `D_3(13)`, `D_3(17)` and `D_2(29)`. Even with full genuine labels, 13 and 17 are not pinned at norm 1.
  * The claim "a flexible top is useful to C1 only together with a non-genuine core" is therefore unproved.
  * It holds when every column satisfies (i). For the flexible-top instance it would need norm-2 pinning (the
    author's open point 2), which is not proved.
* **§5 item 3, "(b) follows (up to the columns failing (i), which are rare, Cor. P1) from Rig_pin"**: the
  qualifier is false. With `U` = the columns failing (i), Cor. P2 gives at most `prod_(j in U)(e_j + 1)` actual
  clusters of a sub-circle. That is vacuous when `|U| ~ 0.9 r`, and a fake builder controls `U` by choosing primes.
  * So **Rig_pin + Prop. P does not make P3-all\* exact in regime (1.4)**.
  * It does so only on members whose prime sets satisfy (i) at every column (or a weaker all-unit pinning
    condition). That restriction is under the fake's control and fails for the smallest primes.
* **"The label-angle step is settled at norm 1 in regime (1.4)"** (§5 item 3) is overstated. It holds for prime
  sets satisfying (i) and per-unit genuine labels, which the author's own caveat about system-level coherence
  already restricts further.
* **N1 (§3.5)** is correct as a necessary condition. The informal reading ("a C1 fake must be non-genuine below
  `M^(h_0/2)`") applies only to columns satisfying (i). A fake whose primes make (i) fail everywhere is not
  constrained by N1.

**Cor. P2: correct as stated.** The `U`-columns contribute a common constant per `tau`-class. The circle is
`N_(U^c)` with constant `C N_U^(-1/4)`. It is only informative when `|U|` is small.

## 4. Lemma L, Prop. G, Thm C

**Lemma L: correct.** It is `round9/L4.md` Lemma 2.1 with `8 -> 2^(E+1)`, and the proof transfers verbatim:
* `q = -1 mod 2^(E+1)` with `E >= 2` gives `q = 7 mod 8`, `(2/q) = 1` and `(q-1)/2` odd.
* The reciprocity step and the CRT choice are as in Lemma 2.1. `q` does not divide `n` because `q = b mod m` and
  `q = 1 mod l`.

On the gloss "any odd units `u_q`": in `L4.md` (2.1) these are multiplicative odd scalars plus `(m_Q/4) Z`
shifts. Both preserve parity on classes with `k(q) >= 3`, so the gloss is right.

Novelty is modest. Theorem A of `L4.md` already gives this unconditionally on `7 mod 8`. Lemma L only extends it
to each deeper class separately. "Unconditional form of round 7's GRH statement" is fair, for **exact**
(unit-scaled) coherence only.

**Prop. G: correct.** Steps checked:
* (a) The parity of `log(xi/conj xi)` is `[(Norm xi/q) = -1]`. Theorem A applies on `7 mod 8`. The `3 mod 8`
  class excludes `2t^2`, and those `q` never divide `N`. Primitive factorisation then gives
  `xi = u varpi^(2w+1) sigma^2`.
* (b) `D_1`, `D_3` are odd. A prime of `D_s(xi)` cannot divide `Norm xi` (by primitivity). Collision at `q^a`
  holds iff `q^a | D_s(xi)`, by `sharp.md` Prop. 1.2. With the budget, `|D_s| <= sqrt(2p) g`, and with
  `D_1^2 + D_3^2 = 2|xi|^2` this gives `|xi|^2 <= 2 p g^2`.

**Thm C: correct**, given its hypotheses: exact coherence on **all** odd `q <= X`, `Norm xi_j <= M^h`, no `p_l`
dividing `D_1 D_3(xi_j)`, and (i) for every `j`.
* The same MAJOR caveat as §3 applies: the cofactor conditions are "generic" only for prime sets without small
  primes.
* The phrase "every exactly coherent member in regime (1.4) is an actual cluster" must carry "whose prime set
  satisfies the cofactor conditions". A builder choosing the smallest primes violates them in most columns.

## 5. Thm R (pairs, Paley profile)

**Correct.** Re-derived:
* Steps 1-4 of `p3Dconstruct.md` Thm 2 give `|T*| <= 4(C + 2Hbar)`.
* If `|T*| <= (M-1)/2`, Siegel's lemma (`n = M-1`, `m = |T*|`, coefficients `<= bC`) gives `beta != 0` with
  `||beta||_inf <= (M-1) b C`.
* Then `v = sum beta_a s_a` has `Gv = 0` and `Z_v = 1`. By Prop. 3.3(4), `v` collides at every `q^(A_q)`,
  `q = 3 mod 8`, with a single unit.
* Summing the four per-unit budgets gives `2w(v) + log 2 <= 2(M-1) b C W + log 2`. Contradiction with (3.1).
* The regime check (`A >= 4`, `C <= M`, and `C > M` trivial) is right.

Inherited minor point (from `p3Dconstruct.md` Def. 4.1): "`q = 3 mod 8` never divide `Norm gamma_t`" is false if
some `gamma_t` is an inert prime `q` itself. Those finitely many `q` must be dropped from `theta(X; 8, 3)`, which
costs `<= sum_t log|gamma_t|` and is harmless.

Scope: Paley profile and `q = 3 mod 8` coherence only.

## 6. Thm S, Cor. S1 (shallow classes)

**Correct.** Re-derived:
* `k(l_j) >= E_0 + 1`, so the `zeta_j` are units at shallow `q`.
* Lemma 3.2 of `p3Dconstruct.md` gives validity at every level.
* Raising a collision to the `2^(k(q))`-th power kills `f` and `i^s`.
* `q^a | Im(gamma_v^(2^k))`, and coprimality across `q` gives the bound per class. Summing gives `(2^(E_0) - 2) sum |v_j| log l_j`.
* In Cor. S1, `prod (1 + 2x/(1-x)) - 1 <= exp(4 sum x) - 1` with `sum x <= 1/40`. The exponent `2^(E_0+1) - 2` is right.

It is a construction-side localisation, as the author says. It covers monomials outside `L` only, and **all**
`p_j >= M^(2^(E_0+1) - 2 + eps)`.

## 7. §3.5 (N1-N6) and the barrier logic (§5)

* **MAJOR (overclaim): "every mechanism of cross-modulus coordination known from rounds 5-14 (... Gaussian
  coherence with any repair) is excluded by N2-N6".** Root-repaired coherence on the deep classes is not excluded
  by anything proved:
  * N2 and N3 concern **exact** coherence, and root repair is not exact.
  * N6 concerns only the Paley profile and only `q = 3 mod 8`.
  * The only evidence against deep-class root repair is §2.4 ("certified cost `2^k log|gamma_v|`", i.e.
    uncertified) and `deep_demand.py`. That is one random choice of `zeta`'s, level 1 only, evidence.
  * The text in §3.5 ("Gaussian coherence ... is excluded by parity (N2, N3)") contradicts §2.4/Thm S, where root
    repair defeats parity.
  * Correct statement: root-repaired coherence is **uncertified** (not excluded) on the deep classes and is a
    candidate DeepBal mechanism. DeepBal as defined (non-exactly-coherent balanced labels) does include it.
* **N4 is not a property of a fixed label system.** "Not independent over moduli" is a property of a random law.
  B1 and B2 show that samples of the stated laws fail w.h.p. They do not show that no member of the support is
  good. N4 should be restated as a statement about random construction methods. N5 is fine as stated (deterministic, B3).
* **N1**: see §3. It constrains only columns satisfying (i).
* **Barrier implication (§5 item 1): correct.**
  * `W = 2 log R`. `round8/stability.md` Cor. 3.6 under `H(eps)` gives primes `<= M^(4/eps - 3 + o(1))` and
    `r <= (2/eps + o(1)) M`. `sharp.md` Prop. 1.2 makes actual clusters members for every `X`.
  * So a theorem "for every `h_0` there are `A` and `F > (h_0+3)/2` with `W >= F M log M` on members with
    `p_j <= M^(h_0 + o(1))`" gives (T').
* **Dichotomy (§5 item 3): not valid as stated.**
  * (MAJOR) "(b) follows from Rig_pin + Prop. P up to the rare columns failing (i)" is false. See §3: those columns
    are the majority for the smallest primes, and the builder controls them.
  * (minor) "C1 or (b)" is not a strict dichotomy. The negation of C1 ("for some `c_0`, every `A`, all large `M`")
    yields exactness only for infinitely many `M` and some `A`, not "no non-actual members in regime (1.4)".
  * The author's correct and important point stands: exactness of P3-all\* in regime (1.4) would settle the
    barrier question but would give no growth without L7/(I3) for genuine data.
* **No claim of growth or of a cap is made, and none is warranted.** The bracket
  `3 M log M - O(M) <= W_min^all <= O(M log^2 M/loglog M)` and `M <= (2/3 + eps) log R/loglog R` are unchanged.

## 8. Independent computations (referee scripts, `referee_c1construct_checks/`)

| script | output | result | label |
|---|---|---|---|
| `ref_P0.py` | `ref_P0_out.txt` | Lemma P0: for all 4783 split `p < 10^5`, the arcs for `s in S(theta)` leave no free set of positive length (50-digit arc union) | evidence (proof re-derived by hand) |
| `ref_cond_i.py` | `ref_cond_i_out.txt` | (i) fails in 8/10, 16/20, 42/50, 87/100, 180/200, 448/500, 913/1000 columns of the first `r` split primes; selected sets with 5, 13: 49/50, 199/200, 999/1000 | exact (integer divisibility) |
| `ref_pin_full.py` | `ref_pin_full_out.txt` | genuine labels at all odd `q <= X`, all four units: unpinned at norm 1 in 5/10, 13/20, 34/50, 74/100, 153/200, 388/500 columns (first `r` split primes); `{5,13,17,29}`: 13 and 17 unpinned (free length 0.307, 0.368); random `[10^4, 10^5]`: 0/200 | interval endpoints in 50 digits (positive lengths far from 0; evidence-grade exactness, each positive length is a concrete witness interval) |
| `ref_norm2.py` | `ref_norm2_out.txt` | on `{5,13,17,29}`, adding the norm-2 relations `e_j +- e_k` with the other angles fixed at genuine values pins 13 and 17 | evidence only (others fixed, not joint) |
| `ref_B1.py` | `ref_B1_out.txt` | B1 product vs `(log Y/log X)^2`; skipping `q | N` (first 90 split primes in the set) inflates it by 1.28-2.32 | exact products (floating) |
| `ref_B3.py` | `ref_B3_out.txt` | `Z_r` reproduced, `>= Z_r` vanishing monomials in 1600 random `(sigma, nu)`; own instance 630/630 over budget | exact integer counts |

Witness for the main finding (exact, integer data): `p_j = 13`, `pi_j = 3 + 2i`, prime set containing 5.
* `D_3(pi_j) = 5 = p_1`, so (i) fails.
* Even with genuine labels at every odd `q <= X`, the set of `phi_13` allowed by (1.2) at norm 1 is the
  complement of the four arcs with radii `arcsin(1/sqrt 13)` (`s = 0`), `arcsin(1/sqrt 26)` (`s = 1`),
  `arcsin(3/sqrt 13)` (`s = 2`), `arcsin(1/sqrt 26)` (`s = 3`, the 5 removed).
* That complement is `[0.2810, 0.5880] ∪ [2.5536, 2.8606] (mod pi)`, a positive-length set. Column 13 is not
  pinned at norm 1.

## 9. Corrected statement

> **Status: restricted theorem; no C1, no lower bound, growth unchanged.** In regime (1.4) (`p_j <= M^(h_0)`,
> `X = 3M^A`, `A` after `c_0`), unrefereed:
> * (B1) For the random fibre of `sharp.md` Thm 6.4, each non-twin column survives its norm-1 per-unit budget with
>   probability `<= (1 + o(1)) log Y max(log Y, log M)/(log X)^2`, `Y = sqrt(2 p_j)`. The product skips `q | N`.
>   This tends to 0 as `A -> infinity`.
> * (B2) The torsion-avoiding fibre fails at norm 2 (`A > h_0`).
> * (B3) For `K in {2, 3}`, `h = 1`, `beta_j = j`, the `L4K.md` Thm 2 systems have `>= Z_r = r^3/36 - O(r^2)` norm-4
>   monomials vanishing identically at every `q > q_T`, for all permutations. Norm 4 is `K + 1` for `K = 3` and
>   `K + 2` for `K = 2`.
> * (P, P1, P2) Norm-1 pinning holds **for columns satisfying (i)** with per-unit genuine collisions at their
>   pinning prime powers. P3-G members whose prime set satisfies (i) at every column are actual clusters.
>   * (i) fails for most columns when the prime set contains small split primes, e.g. the first `r` split primes,
>     or any set a builder selects.
>   * There even genuine labels leave most columns unpinned at norm 1.
>   * So these results do **not** show that P3-G collapses in regime (1.4).
>   * They do not exclude flexible tops over genuine cores. On `{5, 13, 17, 29}` columns 13 and 17 are unpinned at
>     norm 1; norm-2 pinning there is evidence only.
>   * They do not show that Rig_pin makes P3-all\* exact.
> * (L, G, C) Exact coherence with the forced parity on a deep class forces `p_j Norm zeta in {t^2, 2t^2}` (Linnik).
>   Exactly coherent members with `Norm xi_j <= M^h` and the cofactor conditions are actual clusters.
> * (R) On the Paley profile, G-coherent designs at `X >= 3M^4` have `C + 2 Hbar > (M - 1)/8`.
> * (S, S1) Shallow classes `k(q) <= E_0` can be served `X`-uniformly by root-repaired coherence for the monomials
>   outside `L`, at cost `(2^(E_0) - 2) sum |v_j| log l_j`.
>
> The residual construction problems are DeepBal and PairCompat. Root-repaired coherence on the deep classes is
> not excluded (only uncertified, with one-instance evidence against it), so "every known mechanism is
> excluded" is withdrawn.
>
> The implication "lower bound in regime (1.4) with (T') quantifiers ⇒ growth" is correct. A rigidity answer would
> make the class exact only on prime sets where norm-1 pinning applies. A pinning statement covering prime sets
> with small primes (norm-2 or higher closures, as in `p3G.md`) is needed before "exact in regime (1.4)" can be
> claimed. Even then growth would still need L7/(I3) for genuine data.
