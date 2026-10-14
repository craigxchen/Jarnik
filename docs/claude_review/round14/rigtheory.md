# Round 14 — The theory of Rig(K)

Status: complete draft, unrefereed. Scripts and outputs: `rigtheory_checks/`. Working notes: `rigtheory_notes.md`.

## 0. Summary of results

**No growth improvement and no cap is claimed.** `M <= (2/3 + eps) log R/loglog R` stands; the P3-all\* bracket
`3 M log M - O(M) <= W_min^all <= O(M log^2 M/loglog M)` is unchanged. Rig(K) is a statement about label systems;
actual circles and every fake with genuine labels satisfy all its conclusions (§1). **Rig(K) is not decided**: growth
versus stalling of the *joint* core remains open. What is new:

**(1) Formulation and proved structure (§2-§3).**
* Rig(K) made precise (Def. 2.1-2.4): one-modulus and joint cores; three coherence notions (genuine up to conjugation,
  coherent, strict, §4.1). Proved reductions (Thm 2.5): only the exact forced set `F_K^ex` matters, cores are monotone
  in `K`, and `q_core(K) <= 4A_K + 5 = M^(floor(K/2) + o(1))` (Bose-Chowla; `M^(2 floor(K/2) - 1 + o(1))` with greedy `B_h`).
* **Slack identity** (Lemma 3.1): the genuine slack of `(v, s)` is `log U(D_s(Z_v)) + log(1/|sin(arg Z_v - s pi/4)|)`;
  extra collisions are affordable only on monomials with a large uncharged cofactor or a near-axis `Z_v` (Cor. 3.2).
* **Anticode rigidity** (Theorem 3.4, proved, purely lattice-theoretic): if an `l_1`-anticode `D` of diameter
  `<= K - 1` is mapped onto `Z/n_q` by the genuine label, every label outside the genuine class has a non-genuine
  collision of norm `<= K` (of norm `<= K - 1` inside `D - D` if `Lambda_c` is generated there). Hence the **F1/F2
  dichotomy** (Cor. 3.5): deformations either reuse only genuine collisions (F2, budget-free) or pay for an extra one
  (F1); F2 is impossible on the covered range, jointly as well; a joint deformation must transport budget (Cor. 3.5(3)).
  The covered range lies inside the pigeonhole-forced range of norm `K - 1` (up to its endpoint) and grows like `(2r)^(h')/h'!`,
  `h' = floor((K-1)/2)` (Prop. 3.6): the structural reason for round 10's "rigid range at `K` ≈ forced range at `K-1`".
* What is **not** proved, and why: F1 cannot be excluded by local certificates (50-90 affordable extras at every covered
  modulus, §3.4); rigidity there is a counting fact about lattices (open point 1). Per-modulus packing (Prop. 3.7) gives
  a joint, all-norm rigidity only on a tiny, non-growing range (§5.8).

**(2) Consequences (§4).** Theorem 4.1 (proved, conditional): if RigS holds (strict, coherent rigidity on the window
`Q_eta = (M^(1-eta), M]`), then members of P3-all\* in the regime have the **genuine** unit-1 pair statistics, so
H^gen (L7 for genuine pair integers) plus `K_u`-uniform margins give `W >= ((8(1-eta) - eps)/K_u) M log^2 M/loglog M`,
i.e. growth for the class and for such actual clusters. The size needed is only `q_core >= M`, reached by the anticode
range at bounded `K` (`K >= 2/gamma + 1` if `r >= M^gamma`). For actual circles alone RigS is vacuous: it is a
class-level transfer of H from designed to genuine labels, and the fake-violated input must be H^gen + L1. Stalling would
refute RigS at every bounded `K` but would **not** by itself cap P3-all\* (§4.3: all norms, pair demands and slab
completion are still needed).

**(3) Data explained and predicted (§5; exact, `r = 4`, `{5,13,17,29}` unless stated).**
* Forced sets certified for `K = 2..7`; our code reproduces every round-10 number it was compared with.
* One-modulus core `n_core(K) = 2, 10, 20, 39, 66, 83` (`q_core = 7, 41, 79, 157, 263, 331`) for `K = 2..7`. Round 10's
  observations: the norm-4 core/top split and the norm-5 behaviour (107-151 rigid, 167/199 flexible with 96/40 F1
  survivors, 211 flexible) are exactly the anticode threshold `a(4,4) = |S_2| = 41`: rigid for `n <= 39`, flexible from
  `n = 41` (163) on, with two isolated exceptions.
* Predictions recorded before computing: (A) anticode law, first flexible `n ≈ a(4, K-1)`: **exact at `K = 6`**
  (rigid through `n = 66`, first flexible `n = 67`), **too high at `K = 7`** (first flexible `n = 84` vs 129);
  (B) round-10 law `n ≈ n*(K-1)`: wrong at both.
* In log scale the core is a stable power of the forced range (`log n_core/log n*(K) ≈ 0.8` for `K = 4..7`, and for
  `r = 5, 6`); with `r` it grows: `K = 3`: 10, 15, 18 (`r = 4, 5, 6`); `K = 4`: 20, 32; `K = 5`: 39, 73 (`r = 4, 5`).
  **Evidence for growth of the one-modulus core** in `K` and `r`; not a proof, tiny `r`.
* **Joint cores, corrected for per-column conjugation** (`pi_j <-> conj pi_j` gives other genuine systems; round 10's
  joint searches mixed them): `Core^eps_J(3) = F_3^ex ∩ [3, 37]` **exactly** (rigidity by ILP infeasibility, flexibility
  by exact all-level witnesses outside every conjugation class at every forced `q_0 >= 41`), i.e. the joint core equals
  the one-modulus core minus one modulus. At
  `K = 4`, conditional on genuine labels at all levels of 3 and 7, every good system is genuine up to conjugation at
  every forced `q_0 in [11, 61] ∪ {73}` (exact ILP infeasibility), and exact witnesses at 67 and 71 give
  `q_core^(J,eps)(4) = 61` (`n = 15`, against 20 for the one-modulus core); moreover every such system is **coherent**
  (one conjugation pattern) on all of `F ∩ [11, 61]`, whereas at `K = 3` exact good systems mix patterns inside the
  core. The joint core grows with `K` (`n = 9, 15`) and becomes coherent at `K = 4`. At `K = 5`, with the moduli above
  163 also fixed genuine, every good system is genuine up to conjugation on the whole one-modulus core `[3, 157]`.
* The all-norm reference (every modulus charged, no free moduli) is more rigid: core `n = 27` vs 20 (`K = 4`), 48 vs 39
  (`K = 5`).

**Answer to the round's question.** *Growing* is supported, not proved: the one-modulus core grows with `K` and `r`, the
joint core (up to conjugation) tracks it where it could be computed, the all-norm version is more rigid still, and no
flexible family at all `K` is known. The proved part is the lattice half (no F2 on the covered range, which grows
polynomially); the unproved half is F1 exclusion (a counting lemma), plus covering, coherence and strictness. Stalling is
not supported by any computation once conjugation is accounted for.

## 1. Setting, regime (1.4), and the exact class

Labels: *proved* (all quantifiers, proof here), *exact* (finite exact computation, script named), *evidence*
(computation that no proof uses), *heuristic*, *conditional* (proved from a named open hypothesis). Nothing here
is refereed. **Nothing here bounds actual circles beyond `M <= (2/3 + eps) log R/loglog R`.**

**Conventions** (those of `round10/L4K.md` §1 and `round7/README.md`).
* Columns `j = 1..r` carry distinct primes `p_j = 1 mod 4`, Gaussian primes `pi_j` above them, `N = prod p_j^(e_j)`,
  `w(v) = sum_j |v_j| log p_j = log Norm(Z_v)`, where `Z_v = prod pi_j^(v_j^+) conj(pi_j)^(v_j^-)`.
* For an odd prime `q` not dividing `N`: `T_(q^a)` is the norm-one torus of `(Z[i]/q^a)^*`, cyclic of order
  `m_(q^a) = (q - chi_(-4)(q)) q^(a-1)`; `m = m_q`, `n = n_q = m/4`. A generator is fixed with `log i = m/4`.
* A **label system** is a family `ell(q^a) in (Z/m_(q^a))^r`, compatible under reduction, for `q` in a set of moduli,
  `q^a <= X`. It is **parity-respecting** if `ell_j(q) = nu_j(q) mod 2`, `nu_j(q) = [(p_j/q) = -1]`.
  By `sharp.md` Cor. 2.2 these are exactly the prime-level residue data of P3-char / P3-all\*.
* `v` **collides at `q^a` with unit `s`** if `<v, ell(q^a)> = s m_(q^a)/4 mod m_(q^a)`. At level 1, `v` collides with
  some unit iff `<v, c> = 0 mod n`, `c = ell(q) mod n`; the collision lattice is `Lambda_c = {v : <v,c> = 0 mod n}`.
* **Demand** of `(v, s)` on a set `S` of moduli: `log m^S_(v,s) = sum_(q in S) a_(q,s)(v) log q`, `a_(q,s)(v)` the
  largest colliding level. **Budget** `B(v,s) = w(v)/2 + [s odd] (log 2)/2`, i.e. `m_(v,s) <= nu_s^(-1) e^(w(v)/2)`.
* **Genuine labels** `L_j(q^a) = log(pi_j/conj pi_j)`: `v` collides at `q^a` with unit `s` iff `q^a | D_s(Z_v)`,
  `(D_0, D_1, D_2, D_3) = (Y, X - Y, X, X + Y)` for `Z_v = X + iY` (re-checked exactly here: `check_genuine.py`,
  57,400 checks, 0 mismatches).
* `V_K = {v : 1 <= ||v||_1 <= K}`, `S_h = {v : ||v||_1 <= h}`.
* **Regime (1.4)** (`round9/README.md` §0, §1 item 1): near-extremal `W <= c_0 M log M`, primes `p_j <= M^(h_0)`,
  `h_0 = 2c_0 - 3 + o(1)`, `r <= rho M`, moduli up to `X = 3 M^A` with `A` chosen **after** `c_0`.

**The class.** P3-all\* `= F^all_X(M, C)` (`sharp.md` Def. 1.1): genuine primes and profile, parity-respecting
prime-level residues, and the slab condition (1.2) for every monomial and unit. Every member meets the per-unit
budget for **every** `v` at every modulus `<= X` (`sharp.md` Prop. 6.1). Actual clusters are members with the
genuine labels (`sharp.md` Prop. 1.2). Rig(K) concerns only the budget part of the class (labels, no angles).

**Why the genuine labels are the reference point (fake test, `referee_L4K.md` §7, proved there).** Every fake circle
has genuine primes; given the genuine labels of its primes it satisfies every budget at every norm and every
level and is coherent with `pi`. So **every statement of the form "good systems are genuine on a range of moduli"
is satisfied by every fake whose primes are in the regime.** Rig(K) is a structural lemma; any barrier-beating power
of a route through it must come from a later step (§4).

## 2. Precise formulation of Rig(K)

Fix primes `p_1, ..., p_r` as in §1, a norm bound `K >= 1` and a finite set `S` of odd primes not dividing `N`.

**Definition 2.1 (good systems).** A parity-respecting label system on `S` (all levels `q^a <= X`, `X` arbitrary)
is **`K`-good on `S`** if `log m^S_(v,s) <= B(v,s)` for every `v in V_K` and every unit `s`. The task's
Rig(K) takes `S = S_K := {q <= Q_K}`, `Q_K = 8K(2r+1)^K + 1` (`round9/README.md` Prop. H).

**Definition 2.2 (forced set).** `q` is *free for norm `K`* if some parity-respecting label at `q` has no collision
of any `v in V_K` with any unit at any level, *forced* otherwise; `F_K^ex` is the set of forced `q` (`round10/L4K.md`
§2). Every `q` with `n_q > A_K = K(2(K+1)^2 B + 2K + 1)` is free (`L4K.md` Thm 1(b)).

**Definition 2.3 (coherence classes at one modulus `q`).** Let `L = L(q)` be the genuine labels, `c = L mod n`.
A label `ell'` at `q` (level 1, `c' = ell' mod n`) is
* in the **genuine class** if `c' = u c mod n` for a unit `u` of `Z/n` (round 10's notion; equivalent to
  `Lambda_(c') = Lambda_c` when `c` is onto, Lemma P);
* **strict** if it has the genuine per-unit level-1 signature on `V_K`: `<v, ell'> = s n` iff `<v, L> = s n`, for all
  `v in V_K` and all `s`. (Strict implies genuine class whenever `Lambda_c` is generated by `Lambda_c ∩ V_K`, Lemma P.)
  This is the computable version; §4.1 uses the structural one (`ell' = u ell`, `u = 1 mod 4`, all `v in Z^r`), and
  adds per-column conjugation.

Strictness is what the pair statistics of `p3lower.md` Thm C see (unit-1 collisions, §4); the genuine class only
fixes collisions up to the unit.

**Definition 2.4 (rigid cores).** Let `F = F_K^ex`.
* **Joint core** `Core_J(K)`: the set of `q in F` such that every `K`-good system on `F` is in the genuine class at `q`.
* **One-modulus core** `Core_1(K)`: the set of `q_0 in F` such that every label at `q_0` that keeps the system
  `K`-good on `F`, when every `q in F \ {q_0}` carries the genuine labels at all levels, is in the genuine class.
  The **relaxed** version charges only level 1 at `q_0` (this is what `round10` and this note compute).
* `q_core^J(K)` (resp. `q_core^1(K)`): the largest `Q` with `F ∩ [3, Q] ⊆ Core_J(K)` (resp. `Core_1(K)`).
* In joint questions the genuine class must be enlarged by per-column conjugation (`pi_j <-> conj pi_j`), which maps
  good systems to good systems; §4.1 defines `Core^eps_J(K)` (genuine up to conjugation, per modulus) and the coherent
  and strict versions. All joint data in §5.6 use `Core^eps_J`.

**Rig(K) (the question, made precise).** Regime (1.4), `r <= rho M`. *How do `q_core^J(K)` and `q_core^1(K)`
grow with `K` and with `M`?* "Growing" means `q_core^J(K) >= M^(g(K))` with `g(K) -> infinity`; "stalling" means
`q_core^J(K) <= M^(O(1))` uniformly in `K` (or bounded).

**Theorem 2.5 (reductions; proved).** Fix `p_1..p_r`, `K`.
1. *(Only the forced set matters.)* If `S ⊇ F_K^ex`, the restrictions to `F_K^ex` of the `K`-good systems on `S` are
   exactly the `K`-good systems on `F_K^ex`. Hence `Core_J(K)` computed on `S_K` equals `Core_J(K)` computed on
   `F_K^ex` whenever `Q_K >= max F_K^ex`. This holds when `4A_K + 5 < Q_K`, e.g. for `r >= r_0(K)` (since
   `A_K = O_K(r^(2 floor(K/2) - 1))` while `Q_K >= 8K(2r)^K`) and in every instance computed here; `Q_K` plays no
   other role.
2. *(Monotonicity.)* `Core_J(K) ⊆ Core_J(K+1)`, the relaxed and exact `Core_1(K) ⊆ Core_1(K+1)`, and `F_K^ex ⊆ F_(K+1)^ex`.
3. *(Upper bound.)* `Core_J(K) ⊆ Core_1(K) ⊆ F_K^ex ⊆ {q : n_q <= A_K}`; so `q_core^J(K) <= q_core^1(K) <=
   4A_K + 5 = O_K(r^(floor(K/2)))` (`A_K = O_K(r^(floor(K/2)))` via Bose-Chowla, or `O_K(r^(2 floor(K/2) - 1))`
   greedily). In regime (1.4): `q_core(K) <= M^(floor(K/2) + o(1))`.

*Proof.* (1) Restriction lowers every demand, so restrictions are good. Conversely, a good system on `F_K^ex`
extended by a collision-free parity-respecting label at every `q in S \ F_K^ex` (which exists by the definition of
free; take every level to be its lift, any lift has no collision because a collision at `q^a` implies one at `q`)
has the same demands. For `max F_K^ex <= 4A_K + 5`: every `q > 4A_K + 5` has `n_q >= (q-1)/4 > A_K`, so it is free by
`L4K.md` Thm 1(b).
(2) `V_K ⊆ V_(K+1)`, so a `(K+1)`-good system on `F_(K+1)^ex` restricted to `F_K^ex` is `K`-good there (by (1)
applied to `S = F_(K+1)^ex ⊇ F_K^ex`); if every `K`-good system is genuine at `q`, so is every `(K+1)`-good one.
For `Core_1`: a label surviving at `K+1` has every norm-`<= K` collision `(v, s)` within
`B(v,s) - (genuine demand on F_(K+1)^ex \ {q_0}) <= B(v,s) - (genuine demand on F_K^ex \ {q_0})`, so it survives at
`K`. A forced `q` for norm `K` is forced for `K+1` since `V_K ⊆ V_(K+1)`.
(3) `Core_J ⊆ Core_1`: the one-modulus systems are `K`-good systems on `F` (exact version); for a free `q` see (1).
The bound on `max F_K^ex` is (1). QED

*Remark 2.6 (free moduli).* At `q` free for norm `K` the core question is degenerate: a collision-free label keeps
every system good, and it lies outside the genuine class as soon as the genuine labels have a norm-`<= K`
collision at `q`. So every rigidity statement about Rig(K) lives on `F_K^ex`, inside `n_q <= A_K`.

So Rig(K) as posed can never give more than `M^(floor(K/2)+o(1))`; growth means `q_core(K)` reaching a positive
fraction of that exponent for every `K`. The all-norm class P3-all\* has no free moduli; its version of the
question is Rig(∞) on `[3, X]` (§4.3).

## 3. Theorems on the rigid core

### 3.1 Where the slack of the genuine system sits

**Lemma 3.1 (slack identity; proved).** Let `S` be a set of odd primes not dividing `N`, `v != 0`, `s in Z/4`, and let
`m^(S,gen)_(v,s)` be the genuine demand on `S` (all levels). Write `D = D_s(Z_v) != 0` and `U_S(D) = |D| / (S-part of D)`
(the part of `D` coprime to `S`). Then

```text
B(v,s) - log m^(S,gen)_(v,s)  =  log U_S(D)  +  log( 1 / |sin(arg Z_v - s pi/4)| )      (both terms >= 0).      (3.1)
```

*Proof.* Genuine collisions: `q^a` collides with unit `s` iff `q^a | D` (§1), so `m^(S,gen)_(v,s)` is the `S`-part of `D`.
`|D| = nu_s^(-1) |Z_v| |sin(arg Z_v - s pi/4)|` (`sharp.md` Prop. 1.2, proof of (e-all)) and `|Z_v| = e^(w(v)/2)`, so
`log(nu_s^(-1) e^(w/2)) = log|D| + log(1/|sin|)`. `D != 0` for `v != 0` by unique factorisation (`L4K.md` Thm 2(i)). QED

**Corollary 3.2 (what a one-modulus deformation can afford; proved).** In the one-modulus setting of Def. 2.4 at `q_0`
(genuine on `F \ {q_0}`), a level-1 collision of `(v, s)` at `q_0` with `v` **not** a genuine collision at `q_0` fits the
budget iff `U_(F)(D_s(Z_v)) / |sin(arg Z_v - s pi/4)| >= q_0`. Two sources of such slack exist:
* **large cofactor**: `D_s(Z_v)` has prime factors outside `F` (above `q*(K)`, or free moduli, or `2`, or the `p_l`);
* **archimedean**: `Z_v` lies within angle `~ U/q_0` of the axis `i^(s/2) R`.

In the all-norm class with every prime `<= |D_s(Z_v)|` charged (P3-all\* with `X >= e^(w(v)/2)` and the genuine system
as reference), only the `2N`-part survives in `U`, and the condition becomes `|sin(arg Z_v - s pi/4)| <= U_(2N)(D_s)/q_0`:
**an extra collision at `q_0` is affordable iff the genuine Gaussian integer `Z_v` is anomalously close to an axis.**
(The `q_0`-adic event `q_0 | D_s` and the archimedean event `|D_s| <~ e^(w/2)/q_0` have the same heuristic probability
`~ 1/q_0`: this is the product formula, and it is why local certificates cannot exclude extra collisions, §3.4.)

*Exact example* (`analyze.py 6 5,13,17,29 269`): at the first `K = 6` flexible modulus `q_0 = 269`, a non-genuine
survivor collides on `(0,2,1,-3)`, `(1,-3,2,0)`, `(1,1,1,2)`, `(2,-3,0,-1)` with `D_s = -2309, 1291, 1303, -1621`
(primes above `q*(6) = 1039`, uncharged) and on `(2,1,-1,1)` with `D_s = 32` (a power of 2). Every extra collision
it uses is of the large-cofactor kind.

### 3.2 Anticode rigidity

An **anticode of diameter `d`** is a finite `D ⊂ Z^r` with `||x - y||_1 <= d` for all `x, y in D`. Examples:
`S_h` (diameter `2h`) and the double ball `S_h ∪ (e_j + S_h)` (diameter `2h + 1`). Let `a(r, d)` be the largest size of
an anticode of diameter `d` (for `r = 4`: the balls/double balls give `2, 9, 16, 41, 66, 129` for `d = 1..6`).

**Lemma 3.3 (short generation; proved).** Let `c' : Z^r -> Z/n` be a homomorphism and `D ⊂ Z^r` with `0 in D` and
`c'(D) = Z/n`. Then `Lambda_(c') = ker c'` is generated by the vectors `rho(x) + e_j - rho(x + c'(e_j))`
(`x in Z/n`, `j = 1..r`), where `rho : Z/n -> D` is any section of `c'|_D` with `rho(0) = 0`. All of them lie in
`Lambda_(c') ∩ (D + S_1 - D)`.

*Proof.* Let `H` be the group they generate; `H ⊆ Lambda_(c')` since `c'` kills each generator. For `z in Z^r` put
`f(z) = z - rho(c'(z))`. Then `f(0) = 0`, and for every `z` and `j`,
`f(z + e_j) - f(z) = rho(c'(z)) + e_j - rho(c'(z) + c'(e_j)) in H`. By induction along a lattice path from 0
(steps `+-e_j`; a step `-e_j` is handled by the same identity at `z - e_j`), `f(z) in H` for all `z`. For
`lambda in Lambda_(c')`, `c'(lambda) = 0`, so `lambda = f(lambda) in H`. QED

**Theorem 3.4 (anticode rigidity; proved).** Let `q` be an odd prime not dividing `N`, `n = n_q`, `c = L(q) mod n` the
genuine level-1 label, `K >= 2`, and `D` an anticode of diameter `<= K - 1` with `c(D) = Z/n`. Let `c' in (Z/n)^r`
be **any** label outside the genuine class.
1. `Lambda_(c')` contains a vector of `((D - D) ∪ (D + S_1 - D)) \ Lambda_c`; in particular a **non-genuine** collision
   `v` with `1 <= ||v||_1 <= K`.
2. If moreover `Lambda_c` is generated by `Lambda_c ∩ (D - D)`, then `Lambda_(c')` contains a vector of
   `(D - D) \ Lambda_c`: a non-genuine collision with `1 <= ||v||_1 <= K - 1`.

No budget, parity or level hypothesis is used: this is a statement about index-`n` lattices.

*Proof.* Translating `D` changes neither the hypotheses nor `D - D`, `D + S_1 - D` (`c(D + t) = c(D) + c(t)`), so
assume `0 in D`. Suppose `Lambda_(c')` contains no vector of `(D - D) \ Lambda_c`. For `d, d' in D` with
`c'(d) = c'(d')`, `d - d' in Lambda_(c') ∩ (D - D) ⊆ Lambda_c`, so `c(d) = c(d')`: the map `g : c'(D) -> c(D)`,
`c'(d) -> c(d)`, is well defined and onto `c(D) = Z/n`. Hence `|c'(D)| = n`, `c'(D) = Z/n`, `g` is a bijection, and
`c'(d) = c'(d') <=> c(d) = c(d')` on `D`, i.e. `Lambda_(c') ∩ (D - D) = Lambda_c ∩ (D - D)`.
(2) The generators of `Lambda_c` lie in `Lambda_c ∩ (D - D) ⊆ Lambda_(c')`, so `Lambda_c ⊆ Lambda_(c')`; then `c'`
factors through `Z^r/Lambda_c ≅ Z/n` (the iso given by `c`, which is onto), so `c' = u c` for some `u in Z/n`, and `u`
is a unit because `c'` is onto. So `c'` is in the genuine class: contradiction.
(1) If also `Lambda_(c') ∩ (D + S_1 - D) ⊆ Lambda_c`, Lemma 3.3 (applicable since `c'(D) = Z/n`) gives
`Lambda_(c') ⊆ Lambda_c`. Both have index `n` (`c, c'` onto), so they are equal, and `c' = u c` with `u` a unit as
before: contradiction. So `Lambda_(c')` meets `((D-D) ∪ (D+S_1-D)) \ Lambda_c`. QED

**Corollary 3.5 (the two flexibility mechanisms; proved).** Call a one-modulus survivor outside the genuine class
**F2** if all its norm-`<= K` collisions are genuine collisions (`Lambda_(c') ∩ V_K ⊆ Lambda_c`), and **F1** otherwise
(it uses an *extra* collision, which must then be affordable in the sense of Cor. 3.2).
1. Under the hypothesis of Theorem 3.4, **F2 survivors do not exist** at `q`, at any `K' >= K`, in the one-modulus and
   in the joint setting alike.
2. Under hypothesis 3.4(2), every non-genuine label (survivor or not, in any system) collides at `q` on an extra vector
   of norm `<= K - 1` lying in `D - D`. If no vector of `(D - D) \ Lambda_c` is affordable at `q` (Cor. 3.2), then
   `q in Core_1(K)` (relaxed, hence exact).
3. *(Joint transport inequality.)* Let `ell'` be `K`-good on `F`, non-genuine at such a `q`, with an extra collision
   `(v, s')` at `q`. Then `sum_(q' in F \ {q}) (a^gen_(q')(v,s') - a'_(q')(v,s'))^+ log q' >= log q - sigma(v,s')`, where
   `sigma = log U_F + log(1/|sin|)` is the genuine slack (3.1) on `F`. So a joint deformation is non-genuine at a
   covered modulus only by spending slack or by *dropping* genuine `(v,s')`-collisions at other moduli, which are
   then non-genuine (or unit-shifted) themselves.

*Proof.* (1) is Theorem 3.4(1) (`D + S_1 - D ⊆ V_K`), (2) is 3.4(2) plus Cor. 3.2. (3) The budget gives
`sum_(q') a'_(q') log q' <= B(v,s') = sigma(v,s') + sum_(q') a^gen_(q') log q'` (Lemma 3.1 with `S = F`); use
`a'_q >= 1`, `a^gen_q = 0` and rearrange. QED

*Exact consistency check.* In every scan of §6 (`r = 4`, `K = 2..7`), F2 survivors occur only at moduli where no
ball/double-ball anticode of diameter `<= K-1` covers `Z/n` (`cert_K*_first.txt`, `deform_K*_first.txt`): the first F2
modulus is 19 (`K = 2`), 47 (`K = 3`), 127 (`K = 4`), 197 (`K = 5`), 347 (`K = 6`), 421 (`K = 7`), all outside the
covered range (at 421 and 431, where `K = 7` has F2 survivors with `n = 105, 108 < |S_3| = 129`, the genuine label maps
`S_3` onto only 77 resp. 67 classes: no covering, as the theorem requires).

### 3.3 Range of the criterion, and why it grows with `K`

**Proposition 3.6 (proved).** Let `h' = floor((K-1)/2)`.
1. Theorem 3.4 can apply at `q` only if `n_q <= a(r, K-1)`, and `a(r, K-1) >= |S_(h')| >= (2r)^(h')/h'! (1 - o(1))`
   (`r -> infinity`, `K` fixed); for even `K` the double ball gives about twice this.
2. If `n_q < a(r, K-1)` then `q` is forced for norm `K - 1` (every `c'` collides on `D - D ⊆ V_(K-1)` by pigeonhole).
   So **the anticode range lies inside the pigeonhole-forced range of norm `K - 1`** (apart from the endpoint
   `n_q = |D|`, where the genuine map is a bijection on `D`): this is the structural reason for
   round 10's empirical pattern "one-modulus rigid range at norm `K` ≈ forced range at norm `K - 1`".
3. In regime (1.4) with `r >= M^(gamma)`, the range `n_q <= |S_(h')|` reaches `M^(gamma h' (1 - o(1)))`, which tends to
   infinity with `K`. If the covering hypothesis held on most moduli of this range, F2-rigidity (Cor. 3.5(1)) would hold
   on a range growing polynomially in `M`, with exponent `~ gamma (K-1)/2`.

*Proof.* (1) `|c(D)| <= |D|`. The ball count is `|S_h| = sum_i 2^i C(r,i) C(h,i) >= 2^h C(r,h)`. (2) `|D| > n` values
`c'(d)` collide; the difference of two colliding points lies in `D - D ⊆ V_(K-1)`. (3) is (1) with `r >= M^gamma`. QED

**The covering hypothesis is not automatic.** It asks that the genuine Gaussian integers `Z_v`, `v in D`, meet every
class of `(Z[i]/q)^*/((Z/q)^* mu_4)` (a cyclic group of order `n`). For a fixed adversarial prime set this can fail
(exact, `anticode2.py 5 5,13,17,29`: `S_2` does not cover at `q = 73, 83, 97, 103` although `|S_2| = 41 > n`; generation
fails at `89, 113, 127, 131`). On average over `q` it is an equidistribution statement for products of `<= h'`
of the given Gaussian primes modulo `q`; for the primes of a near-extremal cluster it is an `L7`/`H`-type input (open),
and like every statement about the genuine primes it is shared by the fakes.

### 3.4 What the theory does not exclude: the F1 mechanism

The sharp form (Cor. 3.5(2)) would make a covered modulus rigid as soon as no extra vector of norm `<= K-1` in `D - D`
is affordable. That certificate **never** applies in the computed range beyond `n = 1`: at `K = 5`, `r = 4` every
covered modulus has 50-90 affordable extra vectors of norm `<= 4` in `S_2 - S_2 = S_4` (`anticode2.py`, exact). The same
is expected in the regime: by Cor. 3.2 the affordable extras at `q_0` are the `v` whose `D_s(Z_v)` has a large uncharged
cofactor or an anomalously small archimedean part, an event of probability `>~ 1/q_0` per vector, against
`|S_(2h')| ~ n^2` candidates. So **rigidity at a covered modulus is not a local certificate but a counting fact**: no
index-`n` lattice other than `Lambda_c` has all its short vectors inside the (sparse) affordable set. This is where
every proof attempt stops; §5 quantifies it heuristically.

### 3.5 A per-modulus packing constraint (proved)

**Proposition 3.7.** (1) In every `K`-good system, every collision `v in V_K` at `q` has `w(v) >= 2 log q - log 2`
(`w(v) >= 2 log q` if the unit is even). (2) The genuine map `v -> L(q)·v mod n` is injective on
`{v : Norm(Z_v) < q/sqrt 2}`; hence `#{v in Z^r : Norm(Z_v) < q/sqrt 2} <= n_q` for every odd `q` not dividing `N`.

*Proof.* (1) One collision costs `log q <= B(v,s)`. (2) If `v != v'` have the same image, `u = v - v'` is a genuine
collision, so `q | D_s(Z_u) != 0`, while `|D_s(Z_u)| <= sqrt2 |Z_u| <= sqrt2 |Z_v||Z_(v')| < q`. QED

So every good label at `q` is a code of minimum `w`-distance `~ 2 log q` and index `~ q/4`, and the genuine one is
such a code for arithmetic reasons. With Theorem 3.4(2), (1) gives the **joint** statement "every `K`-good system is
genuine at `q`" whenever every vector of `(D - D) \ Lambda_c` has `w < 2 log q - log 2`. This is stronger than the
one-modulus certificate of Cor. 3.5(2), which already fails at every computed modulus with `n >= 2`, so this certificate
route gives nothing. The direct exhaustive version (every label at `q` all of whose collisions are individually
affordable is genuine up to conjugation, §5.8) does give a joint, all-norm core, `[3, 31]` for `{5,13,17,29}`, but it
shrinks when columns are added (`{3, 7, 11, 23}` for `r = 5, 6`).

## 4. Conditional consequences: growth / cap dichotomy

### 4.1 Three rigidity notions (a correction to the round-10 bookkeeping)

The genuine labels depend on the choice of `pi_j` versus `conj(pi_j)` in each column: replacing `pi_j` by its conjugate
negates column `j` (`L_j -> -L_j`), and the result is again the genuine system of a valid choice of Gaussian primes
(with the profile column reflected, `a_(xj) -> e_j - a_(xj)`). Write `eps ∘ L` for the labels of the choice
`eps in {+-1}^r`. At one modulus `q`, a level-1 label `c'` is
* **`eps`-genuine** if `c' = u (eps ∘ L) mod n` for a unit `u`; **genuine up to conjugation** if this holds for some `eps`;
* **strict for `eps`** if `ell' = u (eps ∘ L) mod m_q` with `u in (Z/m_q)^*`, `u = 1 mod 4` (then every per-unit level-1
  collision of every `v in Z^r` is the `eps`-genuine one).

A system is **coherent** on a set `S` of moduli if one `eps` works at every `q in S`. Accordingly:
* `Core^eps_J(K)`: moduli where every `K`-good system is genuine up to conjugation (the `eps` may depend on `q`);
* `Core^coh_J(K)` (resp. strict): sets `S` on which every `K`-good system is coherent (resp. strict and coherent).

`Core^coh_J ⊆ Core^eps_J`, and round 10's fixed-`eps` class sits in between for one-modulus tests (where the other moduli
fix `eps = 1`; there all three agree). **In joint tests the distinction is essential:** the round-10 `K = 3` joint
searches ("outside the genuine class at every forced modulus except `{3, 7}`") used the fixed class. Our exact `K = 3`
witnesses (§5.3) are genuine up to conjugation at 19-41, with `eps` changing from modulus to modulus, and outside
every conjugation class only from 41-43 on; and the `eps`-ILP proves `Core^eps_J(3) ⊇ {3, ..., 37}`.

### 4.2 What growth of the core gives (conditional theorem)

Let `Q_eta = {q prime : M^(1-eta) < q <= M, q not dividing 2N}` (the window of `p3lower.md` Thm C).

> **(RigS(eta, c_0)).** There are `K = K(c_0)` and `A_1(c_0)` such that for `A >= A_1` and `M >= M_0`: every member of
> `F^all_(3M^A)(M, C)` with `p_j <= M^(h_0)`, `r <= rho M`, `W <= c_0 M log M` is **strict and coherent** on `Q_eta`
> (one `eps` for all `q in Q_eta`).
>
> **(H^gen(eta, delta)).** For every such member, with its profile `(a_x)`, windings `(k_x)` and the Gaussian primes
> `pi^eps` of its `eps`, the integers `c^gen_(xy) = #{q in Q_eta : q | D_(k_x - k_y)(Z^eps_(a_x - a_y))}` satisfy the
> factorial-moment inequalities of Hypothesis `H(eta, delta)` (`p3lower.md` §4).
>
> **(L1(K_u, lambda_0)).** `K_u`-uniform margins and `lambda_Q >= lambda_0`.

**Theorem 4.1 (transfer; conditional on RigS, H^gen, L1).** Fix `c_0, eta, delta, lambda_0, K_u` and `eps' in (0,1)`.
Under RigS, H^gen and L1, every member as in RigS has `W >= ((8(1-eta) - eps')/K_u) M log^2 M/loglog M` for
`M >= M_0`; hence for large `M` there is no such member with `W <= c_0 M log M`. Since actual clusters are members
(`sharp.md` Prop. 1.2), this holds for near-extremal actual clusters with `K_u`-uniform margins, and gives
`M <= (K_u/(4(1-eta)) + o(1)) log R logloglog R/(loglog R)^2` for them (`p3lower.md` Cor. C.4).

*Proof.* Pair `(x, y)` collides with unit 1 at `q` iff `kappa_x = kappa_y`, `kappa_x = <a_x, ell(q)> - iota k_x`
(`sharp.md` Cor. 2.2(2-3) with `c = e_x - e_y`), i.e. `<a_x - a_y, ell(q)> = iota (k_x - k_y) mod m_q`. Under strict coherence
`ell(q) = u (eps ∘ L(q))` with `u = 1 mod 4`, so `u^(-1) iota = iota mod m_q` and the condition becomes
`<a_x - a_y, eps ∘ L(q)> = iota (k_x - k_y)`, which by the genuine collision identity (§1, for the primes `pi^eps`) is
`q | D_(k_x - k_y)(Z^eps_(a_x - a_y))`. So the member's `c_(xy)` equal `c^gen_(xy)`, the member satisfies `H(eta, delta)`, and
`p3lower.md` Thm C (a member of `A_X` with `X >= M`) gives the bound. QED

**What this does and does not say.**
1. *For actual circles alone RigS is vacuous* (their labels are genuine), so Theorem 4.1 adds nothing to Thm C applied to
   actual clusters. Its content is at the **class level**: it lets a growth argument for P3-all\* in regime (1.4) use
   `H` for **genuine** pair integers (the open L7 = (I3) of route 3) instead of `H` for **designed** labels, which is
   false in general (`sharp.md`'s `D_A` violates `H`, `p3lower.md` Lemma C.2). In barrier-map terms: under RigS,
   P3-all\* in the regime has the pair statistics of P3-G, and every cap of P3-all\* there must be a P3-G fake (genuine
   labels, fake angles/windings). P3-G fakes are known only for `log p >= 2.08 X`.
2. *The fake test* (§1): every fake with genuine labels satisfies RigS's conclusion. The fake-violated input, if any,
   is H^gen + L1 (statements about genuine pair integers with the member's profile and windings): that is where the
   barrier must be beaten. RigS only routes P3-all\* to it.
3. *Size needed.* RigS asks for a strict coherent core containing `Q_eta`, i.e. `q_core >= M` (exponent 1, not `L(h_0+h)`;
   neither Linnik nor Heath-Brown is needed because the conclusion is the genuine class itself, `referee_L4K.md` M1).
   By Prop. 3.6 the anticode range at norm `K` reaches `n_q ~ (2r)^(h')/h'!`, `h' = floor((K-1)/2)`; with `r >= M^gamma`
   it contains `Q_eta` once `gamma h' >= 1`, i.e. `K >= 2/gamma + 1` (`K = 3` already if `r >= M/8`, since
   `|S_1| = 2r + 1`). So **the polynomial size that route 3 needs is reached at bounded `K`, provided rigidity holds on the
   anticode range jointly, coherently and strictly**: the open part is not the size but the joint/coherent/strict
   upgrade and the F1 mechanism (§3.4).
4. *Units.* Theorem 4.1 needs strictness (unit-level agreement), not just the level-1 class: a per-modulus shift
   `ell' = u ell + n w` keeps the class but moves pair collisions to units `i^(<a_x - a_y, w>)`, where they cost nothing
   (`|sin((delta_x - delta_y - arg eta)/2)| >= 1/sqrt 2` for `eta != 1`). In the one-modulus scans (§5.2) the class
   survivors inside the rigid range are strict at every modulus **except** `q_0 = 41` (`K = 3`), `79` (`K = 4`), `251`
   (`K = 6`): there exactly half of them need a unit shift `w != 0` (`nonstrict.py`, exact: 8 of 16, 16 of 32, 36 of 72).
   So the per-unit budgets pin the units at most, not all, rigid moduli; strictness is a further condition.

### 4.3 What stalling would give

> **(Stall(eta)).** For every `K` there are `K`-good systems (all levels) that are not coherent on a subset of
> `Q_eta` carrying a fixed positive proportion of its `theta`-mass.

1. Stall refutes RigS for every bounded `K` (definitions). Then no argument whose label input is budgets of norm
   `<= K` can reach Theorem 4.1; a growth proof through Thm C would have to establish `H` for designed labels, which
   fails for `D_A`.
2. **Stall does not by itself cap P3-all\*.** A cap (C1 in regime (1.4)) needs, simultaneously and with constants
   independent of `A`: (i) budgets of **all** norms at all moduli `<= X = 3M^A` (bounded-norm systems typically
   violate the next norm at a positive proportion of moduli, `referee_L4K.md` §4 Remark 2); (ii) pair demands
   `O(log M)` on a near-extremal profile with `p_j <= M^(h_0)`; (iii) the slab completion, e.g. `m_v <= M^(beta ||v||_1)`
   for all `v` outside the character lattice (`sharp.md` Prop. 6.6). The flexible tops are seeds for (i) at small norms
   only. The all-norm, all-modulus version of the question is Rig(∞) on `[3, X]`.
3. **Rig(∞) on the regime.** P3-all\* has no free moduli: at `q` with `n_q > A_K` the norm-`> K` budgets bite. The
   all-norm analogue of the one-modulus test (reference: genuine everywhere, `X >= e^(w/2)`) has affordable extra
   collisions only on archimedean near-axis monomials and on monomials whose `D_s` has a large `2N`-part (Cor. 3.2),
   a much smaller set than in Rig(K) (where the uncharged large-cofactor monomials dominate, §3.1). So Rig(∞) is
   expected to be **more** rigid than every Rig(K): the large-cofactor slack exists in Rig(K) only because the moduli
   above `F_K^ex` are left free. Confirmed at `K = 4, 5` (§5.7).

## 5. Explanation of round-10 data and predictions for norms 6, 7

All data in this section: `r = 4`, `ps = {5, 13, 17, 29}` (the round-10 set) unless stated; one-modulus scans with
the reference system genuine (all levels) on `F_K^ex \ {q_0}`, collision-free above, level-1 relaxation at `q_0`
(the round-10 convention); **exact** counts (`deform.py`, C enumerator, weighted orbit counting under
`{u = 1 mod 4}`). Our code reproduces every round-10 number it was compared with: `K = 4` at 47/83/97/101/103/107
(8/0, 36/24, 16/0, 20/0, 48/0, 216/144) and `K = 5` at 107/113/127/151/167/199/211 (18/0, 24/0, 32/0, 36/0, 120/96,
120/40, 884/832).

### 5.1 The forced sets (exact, certified)

`forced_cert.py`, output `forced_cert_5_13_17_29.txt`: every prime up to `4A_K + 5` (Theorem 1(b) of `L4K.md` with a
verified greedy `B_h` set) was decided by exhaustive DFS (C), every free vector re-verified.

| `K` | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|
| `#F_K^ex` | 7 | 17 | 40 | 65 | 115 | 173 |
| `q*(K)` (`n*`) | 41 (10) | 103 (26) | 199 (50) | 349 (87) | 1039 (260); body to 647 (162) | 2089 (522); body to 1097 (274) |
| anticode `a(4, K)` (ball/double ball) | 9 | 16 | 41 | 66 | 129 | 192 |

### 5.2 The one-modulus core for `K = 2..7` (exact)

`q_core^1(K)` is the largest `Q` such that every forced `q_0 <= Q` admits only genuine-class survivors.

| `K` | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|
| `q_core^1(K)` (`n`) | 7 (2) | 41 (10) | 79 (20) | 157 (39) | 263 (66) | 331 (83) |
| first flexible `q_0` (`n`) | 11 (3) | 43 (11) | 83 (21) | 163 (41) | 269 (67) | 337 (84) |
| isolated rigid moduli above it | -- | -- | 89-103 | 181, 191 | 271-283 | 347, 353, 359, 373, 379, 389, 401, 409, 439 |
| `a(4, K-1)` | 2 | 9 | 16 | 41 | 66 | 129 |
| `log n_core / log n*(K)` | 0.30 | 0.71 | 0.77 | 0.82 | 0.82 | 0.79 |

(Outputs `deform_K{2..6}_first.txt`, `deform_K7_lo.txt`, `deform_K7_hi.txt`; `K = 6` scanned at every forced `q_0`
except 1039, `K = 7` at every forced `q_0 <= 691`.) In the rigid range every survivor is **strict** (genuine per-unit signature) except at `q_0 = 41`
(`K = 3`), `79` (`K = 4`), `251` (`K = 6`), where half the class survivors carry a unit shift; above it the flexible
survivors are overwhelmingly F1 (extra collisions, large-cofactor type, §3.1), with F2 survivors appearing only at
uncovered moduli (§3.2 check).

### 5.3 The round-10 observations, explained

1. **Norm-4 rigid core and flexible top.** At `K = 4` the covered range of Theorem 3.4 (double ball `S_1 ∪ (e_j+S_1)`,
   16 points) is `q <= 53` (`cert_K4_first.txt`): there no F2 deformation exists (proved), and F1 deformations are
   excluded by the counting effect of §3.4. The flexible top `[107, 199]` (`n = 27..50`) lies far above `a(4,3) = 16`,
   where both mechanisms are available; the first F2 survivor is at 127 (`n = 32`).
2. **Why norm 5 rigidifies 107-151 but not 167/199.** At `K = 5` the anticode is `S_2` (`a(4,4) = 41`). The moduli 107,
   113, 127, 151 have `n = 27, 28, 32, 38 <= 41`; 167 and 199 have `n = 42, 50 > 41`; 163 (`n = 41`) is the first flexible
   one. The full `K = 5` scan (§5.2) shows the transition sits exactly at the anticode size: rigid for every forced `q_0`
   with `n <= 39`, flexible at 163, 167, 173, 179, 193, 197, 199 and every forced `q_0 > 199` (181, 191 are isolated
   rigid points). The 96/40 survivors at 167/199 are all F1 (exact), i.e. they use extra collisions on monomials whose
   `D_s` has a prime factor above `q*(5) = 349` (uncharged, Cor. 3.2).
3. **211 flexible.** `n = 53 > 41`, free for norm 4 (so its forced collisions are all of norm 5, where the
   large-cofactor slack is plentiful). 832 F1 survivors (exact).
4. **The round-10 pattern "rigid range at `K` ≈ forced range at `K - 1`"** is the inclusion of Prop. 3.6(2) (anticode
   range ⊆ pigeonhole-forced range of norm `K - 1`) together with the coincidence `a(4, K-1) ≈ n*(K-1)` for `K <= 4`
   (2/2, 9/10, 16/26). It separates for larger `K`: `a(4,4) = 41` vs `n*(4) = 50`, `a(4,5) = 66` vs `n*(5) = 87`.
5. **Joint versus one-modulus at `K = 3`.** Round 10's joint searches destroyed the one-modulus core at `K = 3`. With
   per-column conjugation taken into account (§4.1) this is no longer so: the exact `eps`-ILP (levels at 3 and 7 exact)
   proves `Core^eps_J(3) ⊇ {3, ..., 37}` (`n <= 9`), and exact all-level witnesses outside every conjugation class exist at
   41, 43, 47, 53, 71, 103 (`jointG.py`). So `q_core^(J,eps)(3) = 37` against `q_core^1(3) = 41`: **at `K = 3` the joint
   core (up to conjugation) is the one-modulus core minus one modulus.** What round 10 saw was conjugation mixing.

### 5.4 Predictions for norms 6 and 7, and their outcome

Two predictions were recorded (in `rigtheory_notes.md`) **before** the `K = 6, 7` scans were run:
* (A) *anticode law*: first flexible `n ≈ a(4, K-1)`: `n ≈ 66` (`q ≈ 265`) at `K = 6`, `n ≈ 129` (`q ≈ 517`) at `K = 7`;
* (B) *round-10 law*: first flexible `n ≈ n*(K-1)`: `n ≈ 87` (`q = 349`) at `K = 6`, `n ≈ 162` (`q = 647`) at `K = 7`.

Outcome (exact): `K = 6`: rigid for every forced `q_0 <= 263` (`n <= 66`), first flexible `q_0 = 269` (`n = 67`):
**(A) exact, (B) wrong.** `K = 7`: rigid for every forced `q_0 <= 331` (`n <= 83`), first flexible 337 (`n = 84`),
mostly flexible from 419 (`n = 105`): **(A) and (B) both too high** (by factors 1.5 and 1.9).

*Reading (evidence).* The anticode law is a ceiling for the F2-free range where covering holds, not a theorem about F1;
at `K = 7` (`r = 4`) the norm-6/7 monomials have `|D_s|` up to `29^(3.5) ~ 1.3·10^5`, far above `q*(7)`, so the
large-cofactor slack is plentiful and F1 deformations appear below the anticode size. In logarithmic scale the core is
a stable fraction of the forced range, `log n_core/log n*(K) ≈ 0.8` for `K = 4..7`, so the one-modulus core grows with
`K`; whether it grows **with `M`** is a question about `r`, §5.5.

### 5.5 Dependence on `r` (exact; evidence about `M`)

| prime set | `K = 3` core (`n`) / first flex | `K = 4` core / first flex | `K = 5` |
|---|---|---|---|
| `{5,13,17,29}` (`r = 4`) | 41 (10) / 43 | 79 (20) / 83 | 157 (39) / 163 |
| `{5,13,17,29,37}` (`r = 5`) | 61 (15) / 67 | 127 (32) / 131 | (running) |
| `{5,13,17,29,37,41}` (`r = 6`) | 71 (18) / 73 | -- | -- |

(`deform_K3_r5.txt`, `deform_K4_r5.txt`, `deform_K3_r6.txt`, `deform_K5_r5.txt`.) The `r = 5`, `K = 5` scan: rigid for
every forced `q_0 <= 293` (`n <= 73`), first flexible 307 (`n = 77`), 311-317 rigid, flexible from 331 on (forced range
`q* = 719`, `n* = 180`). At `K = 3` the core grows like `~3r` (forced range `~5r-8r`, anticode `2r+1`); at `K = 4` from
`n = 20` to 32 and at `K = 5` from 39 to 73 when `r` goes from 4 to 5. The ratio `log n_core/log n*(K)` is 0.71, 0.74,
0.79 at `K = 3` (`r = 4, 5, 6`), 0.77, 0.81 at `K = 4` (`r = 4, 5`) and 0.82, 0.83 at `K = 5`. Note that at `r = 5` the
core **exceeds** the anticode size (`a(5,2) = 11`, `a(5,3) = 20`, `a(5,4) = 61` against 15, 32, 73): the anticode law of
§5.4 is specific to `r = 4`. *Evidence only*: `r <= 6`, smallest primes.

### 5.6 Joint cores (exact ILP over minimal signatures; witnesses verified in integers)

Method (`joint.py`, `joint2.py`, `jointG.py`): at every forced modulus, all distinct signatures of labels whose
collisions are individually affordable, reduced to inclusion-minimal ones (a subset signature never costs more);
binary choice of one signature per modulus; budget rows for all `(v, s)`; target: a label at `q_0` outside every
conjugation class. Uncharged levels make it a relaxation, so infeasibility is exact rigidity; feasibility is reported
only with a verified all-level witness.

* **`K = 3`** (levels exact at 3 (`<= 3^4`) and 7 (`<= 7^2`), unconditional): every forced `q_0 <= 31` admits **no**
  individually affordable label outside the conjugation classes (packing, §5.8); `q_0 = 37`: ILP infeasible; every forced
  `q_0 >= 41` has an exact witness (41, 43, 47, 53, 59, 61, 67, 71, 83, 103; genuine at all levels of 3 and 7, collision-free
  level-2 lifts elsewhere, all 256 budgets verified in integers). **So `Core^eps_J(3) = F_3^ex ∩ [3, 37]` exactly**, against
  the one-modulus core `[3, 41]`.
* **`K = 4`, conditional on the genuine labels at all levels of 3 and 7** (`G = {3, 7}`, the other 38 moduli free at
  level 1, about `10^5` binary variables): ILP infeasible at every forced `q_0 in {11, ..., 61} ∪ {73}`, i.e. every
  `K = 4`-good system that agrees with the genuine labels at 3 and 7 is genuine up to conjugation there (exact). The ILP
  is feasible (relaxation) at 67, 71, 79, 83, 89, 97, 101, 103; with `G = {3, 7}` the level-2 lift at 11 blocks a witness
  (`n_121 = 30` forces level-2 collisions), so witnesses are built with `G = {3, 7, 11}`: **exact all-level witnesses
  outside every conjugation class at `q_0 = 67, 71, 79, 83, 89, 97, 101, 103, 107, 127, 151, 199`** (`jointG_K4_eps_G3711.txt`, re-verified by
  `verify_witness.py`: 1280 budgets, 0 violations; they deviate from every conjugation class at 67/71 and at most moduli
  above). These witnesses are unconditional (a witness is a good system, whatever `G` was used to find it). The fixed-class
  test with the same `G` gives the same infeasible set (`jointG_K4_G37.txt`). So **given 3 and 7, `q_core^(J,eps)(4) = 61`
  (`n = 15`)**, against the one-modulus core 79 (`n = 20`) and round 10's `{3, ..., 47}` with all other moduli fixed;
  without the condition, `q_core^(J,eps)(4) <= 61` (the witness at 67).
* **`K = 5`, conditional on genuine labels at all levels of 3, 7 and of every forced modulus above 163** (the
  "free core, fixed top" test of round 10 §4.4, now with the core `F ∩ [11, 163]` (29 moduli) entirely free; signatures
  per modulus up to `~3·10^4`): ILP infeasible at **every** forced `q_0 in [11, 157]`; relaxation feasible at 163 (which
  is one-modulus flexible). So with the top fixed, the joint core (up to conjugation) at `K = 5` is the whole
  one-modulus core `[3, 157]` (exact, conditional; `jointG_K5_eps_top.txt`). The unconditional `K = 5` joint problem
  (`~10^5` signatures per modulus near `q = 200-350`) is beyond this method.
* **`K = 4` unconditional** (levels exact at 3 (`<= 81`) and 7 (`<= 49`)): the ILP at `q_0 = 37` did not terminate within
  the time available (HiGHS, `~10^5` binaries; the LP relaxation is feasible at 37, 41, 47). Open.
* **Coherence** (`coh.py`: columns = (signature, set of conjugation classes containing the label); for every `eps` some
  `q in S` must miss `eps`). `K = 3`, `G = {3, 7}`, `S = F ∩ [11, 37]`: feasible, and the exact `K = 3` witness for
  `q_0 = 41` is genuine up to conjugation at 19, 31, 37 with three **different** patterns `eps` (`(1,-1,1,-1)`,
  `(1,1,-1,1)`, `(1,1,1,-1)`): at `K = 3` the joint core is **not** coherent. `K = 4`, `G = {3, 7}`, `S = F ∩ [11, 61]`
  (11 moduli, 112,654 columns): **ILP infeasible**, so every `K = 4`-good system agreeing with the genuine labels at all
  levels of 3 and 7 is coherent (one `eps`) on all of `F ∩ [11, 61]` (exact, given 3, 7; `coh_K4_G37.txt`). Coherence
  appears between norm 3 and norm 4.
* **Joint core versus one-modulus core** (`r = 4`, conjugation allowed): `n = 9` vs 10 at `K = 3` (exact), `n = 15` vs
  20 at `K = 4` (conditional on 3, 7): the joint core grows with `K` and stays a large fraction of the one-modulus core.

### 5.7 The all-norm reference is more rigid (exact)

`deform_full.py`: the reference system is genuine at **every** odd prime `q != q_0` not dividing `N` (all levels), so the
slack is only the archimedean term and the `2N`-part of `D_s` (Lemma 3.1); this is the norm-`<= K` slice of P3-all\* with
`X >= max |D_s|`. One-modulus core: `K = 4`: every forced `q_0 <= 107` rigid (`n <= 27`, against `<= 79`, `n <= 20` in
Rig(4)); `K = 5`: every forced `q_0 <= 193` rigid (`n <= 48`, against `n <= 39`). At the first flexible modulus the
mechanism is F2 at `K = 5` (`q_0 = 197`: survivors with 3 collisions, all genuine) and F1 via an uncharged `N`-part at
`K = 4` (`q_0 = 109`: extra collision `v = (0,3,0,1)`, `D_1 = 75 = 3·5^2`; the `5^2` is level data at `p_1 = 5 | N`, which
(1.2) does not charge). Outputs `deform_full_K4_first.txt`, `deform_full_K5_first.txt`, `analyze_full.py`.

### 5.7a A counting heuristic for the one-modulus threshold (heuristic, calibrated on the data)

Model the non-genuine labels at `q_0` as `~ n^(r-1)` random lattices (labels up to units), each with `~ |V_K/±|/n`
short collisions, each affordable with probability `alpha = #allowed(v,s)/(4 |V_K/±|)` (the density printed by
`deform.py`). Rigidity is predicted while `(r - 1) log n < (|V_K/±|/n) log(1/alpha)`. With the measured `alpha` at the
first flexible modulus (`0.133, 0.133, 0.165, 0.20, 0.27` for `K = 3..7`) this predicts first flexible
`n ≈ 15, 31, 52, 79, 105`, against the observed `11, 21, 41, 67, 84`: an overestimate by a factor `1.25-1.5`, but it
captures the trend, including the `K = 7` slowdown (driven by the growth of `alpha`: more norm-6/7 monomials have
uncharged large cofactors), and `n = 105` is where `K = 7` flexibility becomes dense. The overestimate is the
structured low-collision lattices that exist near the top of every forced range (Theorem 1(b) type), which the model
ignores. *Asymptotics of the model (heuristic):* `|V_K| ~ (2r)^K/K!` gives a predicted threshold `~ r^(K-1)/log r`,
beyond the forced range `r^(K/2)` for `K >= 3`, as long as `alpha` stays bounded away from 1, i.e. as long as
`log p_j/log r` stays bounded (Dickman: the large-cofactor density at norm `K` is about `1 - rho(log p/log r)`). So in the
model the one-modulus core is limited only by the structured lattices near the top of the forced range, which is what
the stable exponent `log n_core/log n* ≈ 0.8` suggests.

### 5.8 The packing core does not grow (exact)

`packing.py`: `q` is packing-rigid if every label at `q` all of whose collisions are individually affordable
(`q^2 <= Norm(Z_v)(2 if s odd)`, all `v` with `Norm(Z_v) < 2q^2`) is genuine up to conjugation. This uses only
Prop. 3.7(1), so packing-rigid moduli are in every joint core, for every `K` and in P3-all\* itself. Result: `[3, 31]`
for `{5,13,17,29}`; `{3, 7, 11, 23}` (19 not rigid) for `{5,...,37}` and for `{5,...,41}` (scanned to 151). Adding columns
shrinks it: single-collision budgets alone give no growing core.

## 6. Computations (scripts, exact certificates)

All scripts are in `rigtheory_checks/` (written for this round; nothing is imported from the repository; the logic of
`rlib14.py` re-derives round 10's definitions). Heavy loops are in C (`*.c`, compiled to `*.so`, called via `ctypes`).
Each runs as `python3 <script> <args>` from that directory.

| script | what it computes | output | label |
|---|---|---|---|
| `rlib14.py` | tori, genuine labels (normalised `log i = m/4`), parity, `Z_v`, `D_s`, `V_K` | -- | -- |
| `check_genuine.py` | genuine collision iff `q | D_s(Z_v)`; genuine parity | stdout (57,400 checks, 0 mismatches) | exact |
| `freec.c`, `forced.py`, `forced_cert.py` | exact forced sets, certified complete up to `4A_K + 5` (greedy `B_h`) | `forced_cert_5_13_17_29.txt` | exact |
| `deformc.c`, `deform.py` | one-modulus deformations: survivors, genuine class, strict, F1/F2 | `deform_K{2..6}_first.txt`, `deform_K7_{lo,hi}.txt`, `deform_K{3,4,5}_r5.txt`, `deform_K3_r6.txt` | exact (level-1 relaxation at `q_0`) |
| `analyze.py` | anatomy of survivors (collisions, `D_s`, factorisations, slack) | stdout | exact |
| `anticode_cert.py`, `anticode2.py` | covering, generation (Smith form), affordable extras for Theorem 3.4 | `cert_K*_first.txt`, `cert2_K*_first.txt` | exact |
| `nonstrict.py` | class survivors needing a unit shift | stdout | exact |
| `sigc.c`, `joint.py` | all level-1 signatures per modulus, minimal antichains, joint ILP (HiGHS via scipy) | `joint_K4_first.txt` | exact relaxation |
| `sigc2.c`, `joint2.py` | joint ILP with exact multi-level labels at chosen primes and per-column-conjugation classes | `joint2_K4_eps1.txt` | exact relaxation |
| `jointG.py`, `eps_class.py` | witness construction: ILP, actual vectors, collision-free level-2 lifts, integer verification of every budget | `witness*_K*.pkl`, `jointG_K*.txt` | exact (per witness) |
| `verify_witness.py` | independent re-check of a witness pickle (parity, lifts, no level-2 collisions, all budgets for `±v`, class at `q_0`) | stdout | exact |
| `joint2_lp.py` | LP relaxation of the joint ILP (feasible: no shortcut certificate) | stdout | exact |
| `deform_full.py`, `analyze_full.py` | one-modulus deformations with the full genuine reference (all moduli charged) | `deform_full_K{4,5}_first.txt` | exact |
| `sigc3.c`, `coh.py` | coherence ILP: columns (signature, set of conjugation classes), one `eps` must fail somewhere on `S` | `coh_K4_G37.txt` | exact relaxation |
| `packing.py` | packing core (single-collision budgets only; valid for all norms) | `packing_*.txt` | exact |

**Exactness conventions.**
* *Infeasibility* (rigidity) claims come from relaxations: the level-1 relaxation at uncharged levels only lowers
  demand, so "no relaxed system" implies "no system". The ILP uses floating logarithms with a `1e-9` tolerance; no budget
  can be met with equality by a product of primes coprime to `2N` (equality would force `|sin| = 1` or
  `prod q = sqrt(Norm Z_v)`, both impossible), so the tolerance cannot create or destroy feasibility unless two sums of
  logs agree to `1e-9`, which we did not check exhaustively (*caveat*).
* *Flexibility* claims from the joint ILP are made only with an explicit witness verified in exact integer arithmetic
  at all levels (`jointG.py`): parity-respecting vectors at every forced modulus, level-2 lifts with no level-2
  collision (so no deeper ones), genuine all-level data on the fixed set `G`, and `m_(v,s)^2 <= Norm(Z_v)(2 if s odd)`
  for all 640 (`K = 3`: 256) pairs `(v, s)`.
* One-modulus flexibility at `q_0` is reported at the level-1 relaxation (round-10 convention). For the claims used in
  the text the relaxation is immaterial: a non-genuine level-1 survivor lifts to a collision-free level-2 label whenever
  a random `t in (Z/q_0)^r` avoids the `<= |V_K|` affine conditions, which `lift2` checks (*evidence* that the one-modulus
  numbers are exact at all levels; not run for every survivor).

## 7. Open points

1. **F1 exclusion (the missing theorem).** Theorem 3.4 removes F2 on the covered range; nothing proved removes F1.
   The needed statement is a counting lemma: *an index-`n` lattice other than `Lambda_c` cannot have all its short
   vectors inside the affordable set*, where the affordable set has density `~ 1/q` among short vectors (Cor. 3.2) and
   is defined by the genuine arithmetic (large uncharged cofactors, `2N`-parts, near-axis angles). A proof must use the
   structure of the affordable set (it is not random: §3.1 anatomy), since local certificates provably fail (§3.4).
2. **Covering.** Theorem 3.4 needs `c(D) = Z/n` for the genuine labels; this fails at some moduli for a fixed prime set
   and is, on average, an equidistribution statement for products of `<= floor((K-1)/2)` of the given Gaussian primes
   modulo `q` (an L7/H-type input, open, shared by the fakes).
3. **Joint, coherent, strict.** RigS (§4.2) needs all three upgrades on `Q_eta`. Data: joint up to conjugation holds on
   the one-modulus core at `K = 3` minus one modulus (exact) and, given genuine labels at 3 and 7, on `[11, 61]` at
   `K = 4` (75% of the one-modulus core in `n`; witnesses at 67, 71, 79). Coherence (one `eps` on the whole core) fails
   at `K = 3` and holds at `K = 4` on `[11, 61]` (given 3, 7; exact). Joint strictness was not tested. Strictness (no unit shifts) fails at isolated moduli (§4.2 item 4). The unconditional `K = 4` joint ILP
   (levels at 3 and 7 exact) did not finish, and the unconditional `K >= 5` joint test is beyond the method (`~10^5`
   signatures per modulus at `q ~ 200`); with the top above 163 fixed, `K = 5` is jointly rigid on all of `[3, 157]`.
4. **Growth in `M`.** The one-modulus core is a stable power of the forced range, `log n_core/log n*(K) ≈ 0.8`, for
   `r = 4..6`, `K = 3..7` (§5.2, §5.5): **evidence for growth**, at tiny `r`, with the smallest primes (tightest budgets).
   Untested: `r >= 8`, larger primes (`log p/log r` bounded away from 1: `L4K.md` found loose sets far less rigid), and
   the joint core at `K >= 5`.
5. **Rig(∞).** With the full genuine reference (all moduli charged, no free moduli) the one-modulus core grows from
   `n = 20` to `27` (`K = 4`) and from `39` to `48` (`K = 5`) (§5.7): the all-norm class is more rigid than Rig(K), as
   predicted; its first flexible moduli use F2 or `N`-part slack (the `p_l`-parts of `D_s` are uncharged in (1.2)).
6. **The packing core** (per modulus, single-collision budgets only, valid for all norms and in P3-all\* directly) is
   `[3, 31]` for `r = 4` but only `{3, 7, 11, 23}` for `r = 5, 6` (§5.8): it does not grow, so every growing rigidity
   must be a budget-sharing (multi-collision) effect.
7. **What stalling would need.** No flexible family at every `K` is constructed here; the computed flexible tops
   (F1 via uncharged large cofactors) are artefacts of Rig(K)'s free moduli above `F_K^ex` and shrink under the full
   reference (item 5). A C1 seed must work at all norms simultaneously (§4.3).
