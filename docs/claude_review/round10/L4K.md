# L4(K): the scale-dependent inverse step -- exact tests and an inverse/non-inverse verdict

(round 10; unrefereed draft; scripts in `L4K_checks/`, each runs as `python3 <script> [args]` from that directory;
outputs stored next to them)

## 0. Summary and verdict

Labels: *proved* (all quantifiers, here), *exact* (finite computation, script named), *evidence* (computation no
proof uses), *heuristic*. Nothing here is refereed. **No growth improvement and no new cap is claimed**; the
bracket for P3-all\* and the bound `M <= (2/3 + eps) log R/loglog R` are unchanged. Everything concerns label
systems (designs); actual circles give the genuine system, which satisfies every statement here.

**Verdict on L4(K) as posed in `round9/README.md` §6: false, for every `K`** (Theorem 2, proved, regime (1.4) with
`r = Theta(M)`; exact finite instances for `r = 4..8`, `K = 2, 3`). The reason is structural:

1. **Theorem 1 (proved).** Norm-`<= K` budgets act only at moduli with `n_q = m_q/4 <~ r^(floor(K/2))`: below
   `|S_(floor(K/2))| ~ (2r)^(floor(K/2))/floor(K/2)!` every label vector collides (pigeonhole), and above
   `A_K = K(2(K+1)^2 B + 2K + 1)` (`B = r` for `K <= 3`, `B < (2r)^(floor(K/2))` in general) an explicit
   parity-respecting vector has **no** collision at any level. This sharpens README Prop. H from `M^(K + o(1))` to
   `M^(floor(K/2) + o(1))`. Exact boundaries: `q* = 31-89` for `K = 2`, `83-199` for `K = 3` (`r = 4..8`,
   primes `<= 200`), against `Q_2 = 1297-4625` and `Q_3 = 17497-117913`; 97-99.9% of the `Q_7` log-mass below
   `Q_K` is free.
2. **Theorem 2 (proved).** Genuine labels below the free threshold and the vectors of Theorem 1 with independent
   random column permutations above give a parity-respecting system that is good for every monomial of norm
   `<= K` at every level and every `X` (per-unit demand `<=` genuine demand), yet is not system-level coherent with
   any Gaussian integers of norm `<= M^h` on any `delta`-fraction of the `Q_7` log-mass below `Q_K` (or below any
   `Y <= M^(O(1))`). Exact instances: coherent share `<= 4.4-13%` against all `zeta` of norm `<= 2000` or `8000`.
3. **The surviving question** is L4^forced(K) (§4.1): the inverse step restricted to the forced range
   `n_q <~ A_K`, with `K >= 2L(h_0 + h)` (twice the README estimate) if it is to feed Thm A. Lemma P (proved)
   reduces coherence at a forced modulus to equality of collision lattices. **The prescribed first test
   (`K = 2, 3`, `r <= 8`) cannot decide it**: its forced range is `q <= 199`, where coherence is weak or vacuous.
   It does show (exact) that budgets do not force genuine collision patterns: for `{5, 13, 17, 29}`, `K = 2`, there
   are `5.07 * 10^15` good level-1 systems on the forced primes, of which `2^37` (`2.7 * 10^-5`) have the genuine
   signature everywhere; sampled good systems (all levels) are almost never coherent (0-1 of 100 at norm
   `<= 2000`) on all forced primes, and at `K = 3` jointly non-genuine good systems exist on every forced modulus
   except `{3, 7}`.
4. **Norm 4 (exact, evidence for the large-`K` picture; §4.4).** With the smallest primes the forced range
   (`q <= 199`, 40 moduli) splits into a **rigid core** — exhaustively, every good system on `{3, ..., 47}` is in
   the genuine class there, with the rest of the forced range genuine *or* replaced by a non-genuine top; one
   modulus at a time, every `q_0 <= 97` except 83 is rigid — and a **flexible top** `[107, 199]` where good
   systems outside the genuine class at all 18 moduli exist, not coherent (norm `<= 2000`) on top `Q_7` moduli
   carrying 50% of the forced `Q_7` mass. So L4^forced(4) fails at this test in its "most of the mass" form; a
   rigid-core statement is consistent with all data, and rigidity strengthens with `K` and with smaller primes.

**Status for the program.** Route 3 as formulated (L4(K) + Thm A) is refuted; L4^forced(K) is open and is the
corrected missing lemma. The systems of Theorem 2 are not fakes (they violate norm-`> K` budgets), so the barrier
map is unchanged: P3-all\* stays open, uncapped, with no growth.

## 1. Setting and precise definitions

Conventions are those of `round9/L4.md` §1 and `round7/p3Dlower.md` §0.

* Columns `j = 1..r` carry distinct primes `p_j = 1 mod 4`; `N = prod p_j^(e_j)`; `w(v) = sum_j |v_j| log p_j`.
* For an odd prime `q` not dividing `N`: `m = m_q = q - chi_(-4)(q)`, `n = n_q = m_q/4`. A label vector at level
  `q^a` lives in `(Z/m_(q^a))^r`, `m_(q^a) = m_q q^(a-1)`; the top level is `Q = q^(A_q) <= X`.
* **Parity-respecting** (P): `ell_j(q) = nu_j(q) mod 2`, `nu_j(q) = [(p_j/q) = -1]`.
* `v` **collides at `q^a` with unit `s`** if `<v, ell> = s iota mod m_(q^a)`, `iota = log i` (`iota = m/4` or
  `3m/4`; changing it swaps `s = 1, 3`, which have the same budget). At level 1, `v` collides with *some* unit iff
  `<v, c> = 0 mod n`, where `c = ell mod n` (`L4_notes.md` N1). A collision at `q^a` implies one at `q^(a-1)`.
* `m_(v,s) = prod_q q^(a_(q,s)(v))` (largest colliding level). **Per-unit height condition / budget**:
  `log m_(v,s) <= B(v,s) := w(v)/2 + [s odd] (log 2)/2`, i.e. `m_(v,s) <= nu_s^(-1) e^(w(v)/2)`. Genuine labels
  `L_j = log(pi_j/conj pi_j)` meet it because `v` collides at `q^a` with unit `s` iff `q^a | D_s(Z_v)`
  (`D_0 = Y, D_1 = X - Y, D_2 = X, D_3 = X + Y` for `Z_v = X + iY`), and `|D_s(Z_v)| <= nu_s^(-1)|Z_v|`.
  "`beta`-good with the genuine `beta`" in README Prop. H means exactly this; we call such a system **good up to
  norm `K`** (for all `1 <= ||v||_1 <= K` and all `s`). A system good on a set `S` of moduli means the budgets
  are imposed on the part of `m_(v,s)` coming from `q in S`.
* **System-level coherence** (`round9/README.md` §6): at `q`, the system is coherent with Gaussian integers
  `zeta = (zeta_j)` (units mod `q`) if `ell(q) = u_q L^zeta(q) + (m_Q/4) w` for one `u_q in (Z/m_Q)^*` and some
  `w in Z^r`, `L^zeta_j = log(zeta_j/conj zeta_j)`. At level 1 this says `c(q) = u_q L^zeta(q) mod n_q`, one unit for
  all columns. `Q_7 = {q prime, q = 7 mod 8}`; `theta_7(S) = sum_(q in S ∩ Q_7) log q`.
* **L4(K)** (`round9/README.md` §6, verbatim content): regime (1.4) (`p_j <= M^(h_0)`, `r <= rho M`), a
  parity-respecting system good up to norm `K` at the moduli `q <= Q_K := 8K(2r+1)^K + 1`; is it system-level
  coherent with `zeta_j` of norm `<= H = M^h` on all but a small part of `theta_7([1, Q_K])`?

Two sets of monomials are used: `V_K = {v in Z^r : 1 <= ||v||_1 <= K}` and the ball
`S_h = {v : ||v||_1 <= h}`, `|S_h| = sum_(i>=0) 2^i C(r,i) C(h,i) >= (2r)^h/h! (1 - o(1))` (`r -> infinity`).

**Validity for actual circles.** Everything below concerns label systems (designs) in P3-all\*. Actual clusters
give the genuine system, which is good at every norm and coherent with `zeta_j = pi_j`. Nothing in this note
bounds actual circles beyond `two_thirds`.

## 2. Where bounded norms act: forced and free moduli (Theorem 1, exact thresholds)

Call `q` **free for norm `K`** (given the parity pattern `nu(q)`) if some parity-respecting label vector at `q`
has no collision of any `v in V_K` with any unit at any level `q^a`; otherwise **forced**. Since a collision at
`q^a` implies one at `q`, `q` is free iff some parity-respecting `c in (Z/n_q)^r` (with `c_j = nu_j mod 2` when
`n_q` is even; no condition when `n_q` is odd, `round7/p3Dlower.md` Thm 1.1(1)) has `<v, c> != 0 mod n_q` for all
`v in V_K`.

**Theorem 1 (the forced/free boundary is `n_q ~ r^(floor(K/2))`; proved, all quantifiers).** Let `K >= 1`,
`h = floor(K/2)`, `r >= 1`, and let the `p_j` and `q` be arbitrary.

(a) *(forced)* If `n_q < |S_h|`, every label vector at `q` (parity-respecting or not) has a level-1 collision of
some `v` with `1 <= ||v||_1 <= 2h <= K`.

(b) *(free, explicit)* Let `beta_1 < ... < beta_r` be positive integers forming a `B_h` set (all sums of `h`
elements, repetition allowed, are distinct as multisets; for `h <= 1` take `beta_j = j`) with `beta_r <= B`. Put
`Lambda = K(K+1)B + K` and `A_K = K(2(K+1)^2 B + 2K + 1)`. If `n_q > A_K`, then for **every** permutation `sigma`
of `{1..r}` the vector

```text
ell^sigma_j(q) := a_j mod m_Q,      a_j = 2(Lambda + (K+1) beta_sigma(j)) + nu_j(q)                         (2.1)
```

is parity-respecting and has no collision of any `v in V_K` with any unit at any level `q^a`.

(c) Hence `q` is forced whenever `n_q < |S_h|` (`>= (2r)^h/h! (1 - o(1))`) and free whenever `n_q > A_K`. For
`h <= 1` (`K = 2, 3`) one has `B = r`, `A_2 = 36r + 10`, `A_3 = 96r + 21`. For `h >= 2` a `B_h` set with
`B < (2r)^h` exists by Bose-Chowla (cited from memory: for a prime power `P`, there are `P` integers in
`[1, P^h - 1]` with distinct `h`-fold sums modulo `P^h - 1`; take `P` prime in `[r, 2r)`), so
`A_K <= 2K(K+1)^2 (2r)^h + O(K^2)`. Without the citation, the greedy `B_h` set gives
`B <= h (r + h)^(2h-1) + 1` (each new element must avoid at most that many values), which suffices for every
qualitative statement below with `h` replaced by `2h - 1`.

*Proof.* (a) The `|S_h| > n_q` values `<v, c> mod n_q`, `v in S_h`, cannot be distinct; if `<v, c> = <v', c>`
with `v != v'`, then `u = v - v'` has `1 <= ||u||_1 <= 2h` and `<u, c> = 0 mod n_q`, a collision.
(b) Parity: `m_Q` is even, so `ell_j = a_j = nu_j mod 2`. Let `v in V_K`, `k = ||v||_1`, `sigma_v = sum_j v_j`.
Then `<v, a> = 2 Lambda sigma_v + 2(K+1) <v, beta^sigma> + <v, nu>`. If `sigma_v != 0`, then
`|2 Lambda sigma_v| >= 2 Lambda > 2(K+1) K B + K >= |2(K+1)<v, beta^sigma> + <v, nu>|`, so `<v, a> != 0`.
If `sigma_v = 0`, `<v, a> = 2(K+1)<v, beta^sigma> + <v, nu>` with `|<v, nu>| <= K < 2(K+1)`; so `<v, a> = 0`
forces `<v, beta^sigma> = 0`. Then `v^+` and `v^-` are multisets of the same size `k/2 <= h` with equal
`beta`-sums; padding both with the same element up to size `h` contradicts `B_h`, unless `v^+ = v^-`, i.e.
`v = 0`. So `<v, a>` is a non-zero integer with `|<v, a>| <= K max_j a_j <= A_K < n_q`, hence `n_q` does not divide
it, i.e. `<v, ell> = <v, a> mod m_Q` is not in `(m_q/4) Z/m_q` at level 1, so no unit collides at any level.
QED

*Exact check* (`L4K_checks/bh_check.py`, output `bh_check_out.txt`): for `K = 1..5`, `r = 2..6` (`r <= 5` for
`K >= 4`, Mian-Chowla `B_2` sets), every parity pattern and (for `r <= 4`) every permutation, `<v, a> != 0` for all
`v in V_K` and `max |<v, a>| = A_K` exactly. `ALL OK`.

**Correction to README Prop. H.** Prop. H says bounded norms force nothing above `Q_K = M^(K + o(1))`. Theorem 1
improves this to: **norms `<= K` force nothing above `n_q ~ A_K = O_K(r^(floor(K/2)))`, i.e. above
`M^(floor(K/2) + o(1))`**, and they do force collisions below `(2r)^h/h!`. The README §6 sentence "pigeonhole
forces about `#{||v||_1 <= K}/q` collisions at each `q <= Q_K`" describes the *random* expectation, not a forced
quantity; forcing stops at `n_q ~ |S_h|`.

**Exact thresholds for the first test** (`L4K_checks/free.py`, `free_out.txt`, `free_full3.py`,
`free_full3_out.txt`, table `thresholds_table_out.txt`). Exhaustive DFS (column 1 reduced modulo units, every
other column exhaustive, every `v in V_K` and every unit), at every prime `q <= Q_K`, `q` not dividing `N`
(`K = 3`: exhaustive DFS on `q <= 3000`, first-found search on `(3000, Q_3]`, completed for all 15 sets: no forced prime above 3000). Prime sets: the first `r`, the last `r`, and a random `r` of the split primes `<= 200`. Every
returned free vector was re-verified independently.

| K | r | `q*` = largest forced prime (first / last / random set) | `n*` | pigeonhole: forced if `n <` | Thm 1(b): free if `n >` | `Q_K` | `theta_7` share of `(q*, Q_K]` |
|---|---|---|---|---|---|---|---|
| 2 | 4 | 41 / 31 / 41 | 10 / 8 / 10 | 9 | 154 | 1297 | 0.973 |
| 2 | 5 | 41 / 47 / 47 | 10-12 | 11 | 190 | 1937 | 0.975-0.983 |
| 2 | 6 | 47 / 71 / 73 | 12-18 | 13 | 226 | 2705 | 0.975-0.982 |
| 2 | 7 | 73 / 79 / 79 | 18-20 | 15 | 262 | 3601 | 0.976-0.981 |
| 2 | 8 | 79 / 89 / 71 | 18-22 | 17 | 298 | 4625 | 0.981-0.985 |
| 3 | 4 | 103 / 137 / 83 | 21-34 | 9 | 405 | 17497 | 0.993-0.995 |
| 3 | 5 | 151 / 137 / 109 | 27-38 | 11 | 501 | 31945 | 0.996-0.997 |
| 3 | 6 | 151 / 137 / 151 | 34-38 | 13 | 597 | 52729 | 0.997-0.998 |
| 3 | 7 | 167 / 167 / 151 | 38-42 | 15 | 693 | 81001 | 0.998 |
| 3 | 8 | 173 / 199 / 167 | 42-50 | 17 | 789 | 117913 | 0.998-0.999 |

(*exact*.) The true boundary `n*` is about `2.5r` (`K = 2`) and `5r-8r` (`K = 3`): linear in `r` for both, as
Theorem 1 predicts for `h = 1`, and far below `Q_K` (which is cubic in `r` for `K = 3`). The forced primes are not
an initial segment (e.g. `r = 4`, first set: 37 is free, 41 is forced, because `q = 41` is in `Q_1` and the parity
pattern `nu(41)` removes the free vectors).

## 3. L4(K) as stated is false for every K (Theorem 2, exact finite instances)

**Theorem 2 (bounded-norm goodness does not force coherence on most moduli; proved).** Fix `K >= 2`, `h >= 0`
and `delta in (0, 1)`. Let `p_1, ..., p_r` be any distinct primes `= 1 mod 4`, `H >= 2`, and `Y` with
`q_T := 4 A_K + 5 < Y` (`A_K` from Theorem 1(b)). Define a label system by:
* at every odd prime `q <= q_T` not dividing `N`: the genuine labels `log(pi_j/conj pi_j)` at every level;
* at every prime `q > q_T`: the vector (2.1) with a permutation `sigma_q`, the `sigma_q` independent and uniform.

Then:

(i) *(deterministic)* For every `X`, the system is parity-respecting and good up to norm `K` at all moduli
`<= X`: for every `v in V_K` and unit `s`, `m_(v,s)` divides `D_s(Z_v)` (it is the `q_T`-smooth part of the
`X`-bounded part of `D_s(Z_v)`), so `log m_(v,s) <= B(v,s)`, with per-unit demand at most the genuine demand.

(ii) *(non-coherence)* Let `E` be the event that for some `zeta = (zeta_j)` with `1 <= Norm zeta_j <= H` and some
set `S` of primes `q in Q_7 ∩ (q_T, Y]` with `theta(S) >= delta theta_7([1, Y])`, the system is coherent with
`zeta` at every `q in S`. Then

```text
P(E) <= (13 H)^r 2^(pi_7(Y)) ((Y+1)/(4 r!))^(k_0),      k_0 = ceil(delta theta_7([1, Y]) / log Y).                (3.1)
```

(iii) *(regime (1.4))* Let `r >= c M` (near-extremal clusters have `r = Theta(M)` columns, `round8/stability.md`),
`H = M^h`, and `Y` with `Y >= Q_K` and `log Y = O(log M)` (e.g. `Y = Q_K = 8K(2r+1)^K + 1`, or `Y = X = 3M^A`).
Then for `M >= M_0(K, h, delta, c, A)` the right side of (3.1) is `< 1`, and `theta_7([1, q_T]) <= delta
theta_7([1, Y])`. Hence there is a parity-respecting system, good up to norm `K` at every modulus, which is not
system-level coherent with any `zeta` of norm `<= M^h` on any set of `Q_7`-moduli carrying a fraction `2 delta`
of `theta_7([1, Y])`. In particular **L4(K) is false for every `K >= 2`** (and trivially for `K = 1`), for every
choice of primes in regime (1.4).

*Proof.* (i) Parity: genuine labels have the forced parity (`sharp.md` Cor. 2.2); (2.1) does by Theorem 1(b).
At `q > q_T` we have `n_q >= (q - 1)/4 > A_K`, so by Theorem 1(b) no `v in V_K` collides there at any level, for
any `sigma_q`. At `q <= q_T` the labels are genuine, so `v` collides at `q^a` with unit `s` iff `Z_v = i^s conj Z_v` in
`(Z[i]/q^a)^*`, iff `q^a` divides `Z_v - i^s conj Z_v`, which equals `2i Y`, `(1 - i)(X - Y)`, `2X`, `(1 + i)(X + Y)`
for `s = 0, 1, 2, 3` (`Z_v = X + iY`); as `q` is odd this is `q^a | D_s(Z_v)` (exact check:
`L4K_checks/genuine_check.py`). Hence `m_(v,s)` divides `D_s(Z_v)`. And `D_s(Z_v) != 0` for `v != 0`: otherwise the
ideals `(Z_v)` and `(conj Z_v)` coincide, i.e. `prod (pi_j)^(v_j^+) (conj pi_j)^(v_j^-) = prod (conj pi_j)^(v_j^+)
(pi_j)^(v_j^-)`, which by unique factorisation over distinct split primes forces `v^+ = v^-`, i.e. `v = 0`. So
`m_(v,s) <= |D_s(Z_v)| <= nu_s^(-1) |Z_v| = nu_s^(-1) e^(w(v)/2)`.

(ii) Fix `zeta` and `q in Q_7 ∩ (q_T, Y]`. If some `zeta_j` is not a unit mod `q`, coherence at `q` fails by
definition. Otherwise coherence at the top level implies coherence at level 1, which means `c^sigma(q) = u L^zeta(q) mod n_q` for some unit `u` (the `(m/4)w` term
vanishes mod `n_q`). For a fixed `u` put `y = u L^zeta(q)`. Then `c^sigma_j = y_j` says
`2(K+1) beta_sigma(j) = y_j - 2 Lambda - nu_j mod n_q`. The `r` integers `2(K+1) beta_i` are distinct and lie in
`[0, n_q)` (as `2(K+1) B < A_K < n_q`), so `sigma(j)` is determined by `y_j`: at most one `sigma` per `u`. With
`phi(n_q) <= n_q = (q+1)/4 <= (Y+1)/4` units, `P(coherent at q) <= n_q/r! <= (Y+1)/(4 r!)` (`< 1` once
`r! > Y`), independently over `q`, so a fixed `S` is coherent throughout with probability `<= ((Y+1)/(4r!))^(|S|)`. A set `S` with
`theta(S) >= delta theta_7([1, Y])` has `|S| >= k_0`. Summing over the `<= 2^(pi_7(Y))` sets `S` (each needs
coherence on at least `k_0` of its moduli) and the `<= (pi (sqrt H + 1)^2)^r <= (13 H)^r` tuples `zeta` gives (3.1).

(iii) By the prime number theorem in progressions, `theta_7([1, Y]) ~ Y/4` and `pi_7(Y) ~ Y/(4 log Y)`, so
`k_0 >= (delta/4)(1 - o(1)) Y/log Y`. Since `log Y = O(log M) = O(log r)` and `log(r!) >= r log r - r`, we have
`log(4r!/(Y+1)) >= (1/2) r log r` for large `M`. Hence

```text
log P(E) <= r log(13 M^h) + (Y/(4 log Y)) (1 + o(1)) [log 2 - (delta/2) r log r].
```

The bracket is `<= -(delta/4) r log r` for large `M`, and `Y/log Y >= 16 r/log(16 r)` (as `Y >= Q_K >= 16r` and
`t/log t` increases), so the second term is `<= -(delta/4)(1 - o(1)) r^2 log r/log(16r)`, which beats
`r (h log M + log 13)` once `r >= c M` and `M` is large. So `P(E) < 1`. Finally `q_T = O(K^3 (2r)^(floor(K/2)))`
while `Y >= Q_K >= 8K(2r)^K`, so `theta_7([1, q_T])/theta_7([1, Y]) = O_K(r^(-K/2)) -> 0`; a coherent set
therefore has mass `< delta` in the free part plus `<= delta` in the genuine part. QED

*Remarks.*
1. The construction is Prop. H moved down from `Q_K = M^(K + o(1))` to `q_T = M^(floor(K/2) + o(1))`, plus a
   permutation device that makes non-coherence provable. The README's own Prop. H already gave counterexamples
   that are "genuine below `Q_K`"; Theorem 1 shows the genuine part may be confined to `n_q <= A_K`, which is a
   vanishing fraction of the moduli below `Q_K`, and that is what L4(K) asked about.
2. The system meets only the norm-`<= K` conditions. Monomials of norm `> K` collide at the free moduli about as
   often as for random labels, so the system is **not** in P3-all\*; it is a counterexample to L4(K), not a fake.
3. (iii) uses `r -> infinity`. For bounded `r`, `r!` is bounded and the union bound fails; the finite instances
   below are therefore checked directly.

**Finite instances (exact; `L4K_checks/instance.py`, outputs `instance_out.txt`, `instance_H8000_out.txt`).**
Genuine labels at `q <= q*` (the exact boundary of §2) and a random collision-free parity-respecting vector at
every `Q_7` prime in `(q*, Y]` (random greedy, falling back to the exhaustive witness times a random unit). Every
free vector was re-verified against all `v in V_K` and all units. For two pairs of columns and **all** pairs
`(zeta_j, zeta_k)` of Gaussian integers with norm `<= H` (first-quadrant representatives: units quotiented,
conjugates kept; `H = 2000`, i.e. 1573 per column, and `H = 8000`, 6284 per column), the script computes exactly the
`Q_7` log-mass of `(q*, Y]` on which `(c_j, c_k) = u (L^zeta_j, L^zeta_k) mod n_q` for some unit `u` (moduli where a
`zeta` is not a unit are counted as coherent, which only helps coherence). Pair coherence is necessary for system
coherence, so these are upper bounds.

| K | r | set | `Y` | `H` | genuine share `theta_7([1,q*])/theta_7([1,Y])` | max pair-coherent share of the free part | upper bound for any coherent share |
|---|---|---|---|---|---|---|---|
| 2 | 4 | first / last | 1297 = `Q_2` | 2000 | 0.027 | 0.076-0.083 | 0.103-0.110 |
| 2 | 4 | first / last | 1297 | 8000 | 0.027 | 0.086-0.102 | 0.113-0.130 |
| 2 | 6 | first / last | 2705 = `Q_2` | 2000 | 0.018-0.025 | 0.041-0.048 | 0.063-0.069 |
| 2 | 8 | first / last | 4625 = `Q_2` | 2000 | 0.019 | 0.025-0.027 | 0.044-0.045 |
| 2 | 8 | first / last | 4625 | 8000 | 0.019 | 0.031-0.032 | 0.049-0.051 |
| 3 | 4 | first / last | 4700 (`Q_3 = 17497`) | 2000 | 0.022-0.027 | 0.025-0.031 | 0.047-0.056 |
| 3 | 6 | first / last | 4700 (`Q_3 = 52729`) | 2000 | 0.027-0.031 | 0.026-0.032 | 0.053-0.063 |
| 3 | 8 | first / last | 4700 (`Q_3 = 117913`) | 2000 | 0.035-0.045 | 0.022-0.026 | 0.057-0.068 |

So in every test case the system is good up to norm `K` and **no** system-level coherence with Gaussian integers
of norm `<= 8000` holds on more than 13% of the `Q_7` log-mass below `Y`. For `K = 3`, `Y` was capped at 4700
for cost; above `Y` the moduli are free too (`free_full3_out.txt`), so the coherent share below `Q_3` is smaller
still in the same construction. (*exact* for the stated instances; the instances themselves are random draws.)

## 4. The forced region: the exact first test and classification

### 4.1 The corrected question

By Theorems 1-2 the only place where bounded-norm goodness can force anything is the forced range
`F_K = {q : n_q <= A_K}` (in practice the exact forced set of §2). The question that survives is:

> **L4^forced(K).** Regime (1.4), `r = Theta(M)`, `K >= 2L(h_0 + h)`. Suppose a parity-respecting system is good up
> to norm `K` (all levels) at the moduli of `F_K`. Must it be system-level coherent with some `zeta` of norm
> `<= M^h` on all but a small part of `theta_7(F_K)`?

Coherence on `F_K` is a non-trivial condition only when the column patterns outnumber the candidates: per column
there are `prod_(q in F_K) n_q = e^((1+o(1)) theta(F_K))` patterns modulo `n_q` against `<= 13 H` Gaussian
integers and `prod phi(n_q)` shared units, so (heuristic count) coherence is restrictive iff, roughly,
`(r - 1) theta(F_K) >> r log H`. In regime (1.4) with `K >= 2L(h_0+h)` this holds by a huge margin
(`theta(F_K) ~ M^(K/2)`). **At the scale of the prescribed first test (`K = 2, 3`, `r <= 8`, `F_K ⊂ [3, 200]`) it
fails for moderate `H`**, and on the `Q_7` part of `F_K` (three to six primes) it is vacuous; see the data below.

**Lemma P (the label class is determined by the collision lattice; proved).** Let `c, c' in (Z/n)^r` with
`gcd(c_1, ..., c_r, n) = 1`, and let `Lambda_c = {v in Z^r : <v, c> = 0 mod n}`. If `Lambda_c` is generated by
`Lambda_c ∩ V_K` and every `v in Lambda_c ∩ V_K` satisfies `<v, c'> = 0 mod n`, then `c' = u c` for some `u in Z/n`,
and `u` is a unit iff `gcd(c', n) = 1`. *Proof.* `Lambda_c ⊂ Lambda_(c')`, so `v -> <v, c'>` factors through
`Z^r/Lambda_c`, which `c` identifies with `Z/n`; a homomorphism `Z/n -> Z/n` is multiplication by some `u`. QED

So at a forced modulus where the genuine collision lattice is generated by its short vectors (expected when
`#V_K/n_q` is large, i.e. well inside the forced range) and `L(q)`, `c(q)` are both surjective onto `Z/n_q`,
**coherence with `pi` at `q` is the same as having the genuine all-unit collision set there** (one direction is
trivial; the other is Lemma P with `c = L(q)`, `c' = c(q)`). An inverse theorem must therefore show that the budgets force the genuine
collision sets (or those of one other `zeta`) at most forced `q`. Pigeonhole forces only their *number*.

### 4.2 Exact first test, `K = 2`, `r = 4`, tightest prime set (exhaustive)

`L4K_checks/forced_sig.py` (output `forced_sig_2_4_first_out.txt`). Primes `{5, 13, 17, 29}` (the first four split
primes, the smallest budgets), forced set `F = {3, 7, 11, 19, 23, 31, 41}` (exact, §2). For each `q in F` every
parity-respecting `ell in (Z/m_q)^4` was enumerated and grouped by its exact per-unit collision signature
`{(v, s)}`; a system on `F` is good at level 1 iff the signature costs (`log q` on each `(v, s)`) sum to at most
`B(v, s)` for all 20 monomials (up to sign) and 4 units. Exhaustive DFS over signature tuples:

| `q` | 3 | 7 | 11 | 19 | 23 | 31 | 41 |
|---|---|---|---|---|---|---|---|
| parity vectors | 16 | 256 | 1296 | 10^4 | 20736 | 65536 | 160000 |
| signatures | 16 | 128 | 336 | 568 | 640 | 616 | 990 |
| vectors with the genuine signature | 1 | 2 | 4 | 32 | 128 | 1024 | 4096 |
| good systems with a non-genuine signature at `q` | 94% | 98% | 98% | 94% | 88% | 75% | 50% |

* **Good level-1 systems on `F`: 5,066,549,580,791,808** (36,864 signature tuples).
* With the genuine signature at every `q`: `2^37 = 137,438,953,472`, a fraction `2.7 * 10^-5`.

(*exact*.) So even with the tightest budgets the norm-`<= 2` conditions do not force the genuine collision
pattern at any forced modulus. Moreover the genuine signature does not pin the class above `q = 11`: of the
vectors with the genuine signature, only 16/32, 64/128, 256/1024, 512/4096 are in the genuine class
(`c = u L mod n`, `u` a unit) at `q = 19, 23, 31, 41` (all of them at `q <= 11`; `sigclass_check.py`, `sigclass_2_4_first.txt`). So there the genuine collision lattice is not generated by its norm-`<= 2` vectors and
Lemma P does not apply.

### 4.3 Coherence of good systems in the forced range (exact per system; sampling labelled)

`L4K_checks/forced_coh.py`. Good systems are generated either by exact weighted sampling from all good level-1
systems (§4.2 set; signature tuple with probability proportional to its multiplicity, then a uniform vector in each
class) or, for other sets, by a randomized DFS over `q in F` in increasing order. In both cases lifts through **all
levels** are then chosen by an exact search (at each level `t in (Z/q)^r`, exhaustively when `q^r <= 3*10^5`), and
the system is accepted only if every `(v, s)` stays within budget with all levels counted; a level with no
collision ends the search (deeper levels then have none). For each generated system, system-level coherence at
every `q` of a set `S` with some `zeta`-tuple of norm `<= H` is decided **exactly** (DFS over the shared units
`u_q`, with per-column candidate sets of `zeta`).

| K | r | set | generation | `H` | good systems (all levels) | coherent on all of `F` | coherent on `F ∩ Q_7` | genuine signature everywhere |
|---|---|---|---|---|---|---|---|---|
| 2 | 4 | first | exact weighted | 2000 | 100 | 1 | 100 (`{7,23,31}`) | 0 |
| 2 | 4 | first | exact weighted | 20000 | 100 | 68 | 100 | 0 |
| 2 | 4 | last | DFS | 2000 | 100 | 0 | 100 (`{7,23,31}`) | 0 |
| 2 | 4 | random | DFS | 2000 | 100 | 0 | 100 (`{7,23,31}`) | 0 |
| 2 | 6 | first | DFS | 2000 | 50 | 0 | 7 (`{7,23,31,47}`) | 0 |
| 2 | 6 | last | DFS | 2000 | 50 | 0 | 0 (`{7,23,31,47,71}`) | 0 |
| 3 | 4 | first | DFS | 2000 | 30 | 0 | 0 (`{7,23,31,47,71,103}`) | 0 |
| 3 | 4 | last | DFS | 2000 | 30 | 0 | 0 (`{7,23,31,47,79}`) | 0 |

(Outputs `fc_*.txt`; per-system coherence is exact, the proportions are *evidence* about the population of good
systems; the first row is an exactly weighted sample of all good level-1 systems.) Every sampled level-1 good
system admitted a good lift through all levels. Reading: at `K = 2` good systems are almost never coherent on the
whole forced set at `H = 2000`, but coherence restricted to `F ∩ Q_7` is vacuous for `r = 4` (three moduli) and
weak for `r = 6`; raising `H` to 20000 already makes 68% of the `r = 4` systems coherent on all of `F`. This is the
scale problem of §4.1: **the prescribed first test is answered (non-genuine, generically non-coherent good systems
abound), but it is too small to bear on L4^forced at the scale route 3 needs.**

### 4.4 Norm 4: a rigid core and a flexible top (exact; evidence for the large-`K` picture)

At `K = 4`, `r = 4` the forced range is `q <= 199` (first set) / `271` (last set) (`free_K4_out.txt`), 40
moduli, and coherence is restrictive (`prod_(q in F ∩ Q_7) n_q ~ 10^16`). Three exact tests on the tightest set
`{5, 13, 17, 29}`:

1. **One-modulus deformations** (`deform.py`; outputs `deform_*.txt`). Keep the genuine labels (all levels) at every
   forced modulus but one, `q_0`, and the collision-free labels above; enumerate **all** parity-respecting labels
   at `q_0` (up to `10^8` vectors) and keep those for which the system stays good (exact per-unit test against the
   remaining budget; for a sample of the survivors outside the genuine class an exact all-level lift is also
   found). Result: at every `q_0 <= 103` except `q_0 = 83` the **only** surviving labels are in the genuine
   class `c = u L mod n` (`u` a unit); at `q_0 = 83` and at every `q_0 in [107, 199]` non-genuine survivors exist
   (24 at 83; 144-500,160 at 107-199). For comparison: `K = 3` (same primes) is rigid exactly at `q_0 <= 41` and
   flexible from 43 (forced range `<= 103`); `K = 2` is rigid only at `q_0 in {3, 7}` (forced range `<= 41`); the
   loose set (last four primes) at `K = 2` is rigid only at `{3, 5, 7}`, and at `K = 4` it is rigid at every
   `q_0 <= 31` except 23 and flexible from 37 on (forced range `<= 271`; `deform_4_4_last.txt`, complete up to
   `q_0 = 151` when this note was written, flexible at every checked `q_0 >= 37`). *Empirical pattern (evidence):* **the
   one-modulus rigid range at norm `K` is about the forced range at norm `K - 1`** (`41 = q*(2)`, `103 = q*(3)`).
2. **Joint rigidity of the core** (`core_sig.py`, output `core_4_4_first_rest.txt`). With the forced moduli above 47
   genuine, an exhaustive enumeration of **all** good level-1 systems on the core `{3, ..., 47}` (10 moduli,
   relaxation: core levels `>= 2` not charged) finds exactly one signature tuple, the genuine one: all 1,966,080
   good systems are in the genuine class at every core modulus. Without the outer demand (`core_4_4_first.txt`)
   the same core has `1.0 * 10^9` good systems of which only 12.5% are genuine-class everywhere, so the rigidity
   comes from the budgets already spent by the genuine labels further up.
3. **Flexible top** (`forced_prefix.py`, output `prefix_4_4_first_103.txt`). With genuine labels at all forced
   `q <= 103`, a search over the 18 forced moduli in `[107, 199]` finds good systems (all levels) that are
   outside the genuine class at **all 18** of them (3 of the 5 systems found; `prefix_4_4_first_103.pkl` stores the
   first 4, same seed) and are not system-level coherent with any
   `zeta` of norm `<= 2000` on the top `Q_7` moduli `{127, 151, 167, 191, 199}`. These carry
   `25.5/51.1 = 50%` of `theta_7(F)`.

4. **Joint searches over all forced moduli** (`joint_search.py`, outputs `joint_3_4_first.txt`,
   `joint_4_4_first.txt`). At `K = 3` (same primes) 12 of 12 randomized searches found good systems (all levels) on
   all 17 forced moduli that are outside the genuine class at every forced modulus except `{3, 7}` (or `{3, 7, 11}`)
   — including the whole one-modulus-rigid range `11..41` — and at all 9 moduli above 41. So **at `K = 3` the
   one-modulus rigidity is not joint rigidity**: it disappears as soon as the outer labels may also move.
   At `K = 4` the same randomized search from scratch (12 runs, node cap 1500) found no good system on the 40 forced moduli within the time
   available (earlier runs of the same generator stalled at `q = 53`; *evidence*, a search failure proves nothing). Two exact joint tests
   were possible instead:
   * the exhaustive core enumeration of item 2 with the **non-genuine top** of item 3 (stored system index 2 = system 3, outside the
     genuine class at all 18 top moduli) in place of the genuine top (`core_4_4_first_rest_top2.txt`): again
     exactly one signature tuple, all 1,966,080 good core systems in the genuine class at every core modulus;
   * one-modulus deformations at every core/mid modulus `q_0 <= 103` with that non-genuine top
     (`deform_4_4_first_core_given_top.txt`): rigid at every `q_0 <= 97` except 83 (36 survivors), flexible at
     101 (40) and 103 (96).

   So at `K = 4` the rigidity of `{3, ..., 47}` (and, one modulus at a time, of `q_0 <= 97`) survives a
   non-genuine top; at `K = 3` it does not survive joint deformation. *Evidence only* (one prime set, `r = 4`,
   tiny primes; the joint test of the whole range `<= 103` with everything free was not completed exactly).

**Reading (evidence; heuristic beyond the data).**
* At norm 4, with the smallest primes, the forced range splits into a *rigid core* (here `{3, ..., 47}` jointly,
  and `q_0 <= 97` one modulus at a time), where the budgets admit only the genuine class, and a *flexible top*
  (`[107, 199]`), where non-coherent deformations exist. The top carries `50%` of the forced `Q_7` mass, so
  **L4^forced(4) in its "most of the mass" form fails at this test**.
* At norm 3 the analogous one-modulus core (`q_0 <= 41`) is not jointly rigid; at norm 2 only `{3, 7}` are.
  Rigidity therefore strengthens with `K` (it also appears only with small primes: §5 item 4).
* The statement consistent with every computation is a **rigid-core** form: goodness up to norm `K` forces the
  genuine class on an initial segment of the forced range whose length grows with `K`. For route 3 it would have
  to reach `M^(L(h_0+h))` (Thm A) with `K` fixed; the data say nothing about growth in `M`.

## 5. Consequences for route 3; relation to the fakes

1. **Route 3 as formulated fails (proved).** Route 3 was "L4(K) for `K ~ L(h_0 + h)`, then Thm A". Theorem 2
   refutes L4(K) for *every fixed* `K`, in particular for the `K = K(c_0)` route 3 needs, in regime (1.4), with
   `r = Theta(M)` as near-extremal clusters have. The refutation is uniform in the choice of primes and in
   `A`: the system is good up to norm `K` at every `X`.
2. **Where an inverse step can live (proved).** Bounded-norm goodness constrains nothing at moduli with
   `n_q > A_K = O_K(r^(floor(K/2)))` (Theorem 1(b)); it forces collisions only below `|S_h|` (Theorem 1(a)). An
   inverse step that feeds Thm A must certify coherence on `Q_7` up to `C_L (8 p_j Norm zeta_j)^L ~ M^(L(h_0+h))`,
   so it needs `floor(K/2) >= L(h_0 + h)(1 + o(1))`, i.e. **`K >= 2L(h_0 + h)`** (twice the README's estimate),
   and it must be stated on the forced range only. The corrected statement is L4^forced(K) of §4.
3. **What the inverse step would have to show (Lemma P).** At a forced modulus the label class is determined by
   the collision lattice: coherence with `zeta` at `q` is the identity `Lambda_(c(q)) = Lambda_(L^zeta(q))`, where
   `Lambda_(L^zeta(q)) = {v : q | Im((Z^zeta_v)^4)}`. The budgets bound collisions from above and pigeonhole
   forces about `|S_h|^2/(2 n_q)` of them at each forced `q`; an inverse theorem must show that the forced ones
   are, at most forced `q`, the genuine-type ones. The exact data of §4 show that at `K = 2, 3` they are not
   forced to be, and that at `K = 4` (smallest primes) they are on a core of the forced range but not on its top.
4. **Load accounting (heuristic, labelled).** In the forced range `q <= q_T ~ (2r)^(K/2)`, genuine labels put on
   `(v, s)` the `q_T`-smooth part of `D_s(Z_v)`, of mean logarithm `~ log q_T ~ (K/2) log r` (Mertens), against a
   budget `w(v)/2 ~ (K/2) log p`. So the forced range is **tight exactly when `log p_j/log r -> 1`** (primes as
   small as `r` allows) and has a slack factor `log p/log r` otherwise. Independent random labels in the forced
   range fail the union bound at every slack factor `O(1)` (each `(v, s)` exceeds its budget with probability
   bounded below; there are `(2r)^K` of them), so a counterexample to L4^forced(K) needs structure, and an inverse
   theorem must exploit the tightness or the arithmetic of the forced collisions. The tightness is per `(v, s)`,
   not aggregate: on the test sets the level-1 forced mass is only 29-41% of `sum_(v,s) B(v,s)` and the genuine
   demand 43-60% (`tightness_out.txt`), yet the `K = 4` small-prime system is locally rigid (§4.4).
5. **Relation to the fakes.** The systems of Theorem 2 are counterexamples to a lemma, not fakes: they violate the
   budgets of monomials of norm `> K` at the free moduli (Remark 2). The input they falsify, and that actual
   circles satisfy, is goodness of monomials of norm `> K` at moduli above `M^(floor(K/2) + o(1))`; P3-all\*
   contains it. The fakes of rounds 5-7 all have `p_j > X` (round 9 §1 item 1) and are outside regime (1.4); no
   fake of size `O(M log M)` is known inside regime (1.4). Nothing here changes the barrier map: P3-all\* stays
   open at fixed `A` and `A`-uniformly; no cap and no growth are proved.
6. **What is not affected.** Thm A, B, C (round 9) are untouched; L6 stays proved in regime (1.4). Prop. H stands
   and is sharpened by Theorem 1. Construction D is irrelevant here.

## 6. Open points

1. **L4^forced(K)** (§4.1) for `K >= 2L(h_0 + h)` in regime (1.4): open in the form "coherent on most of the
   forced `Q_7` mass", and the `K = 4` test (§4.4) already contradicts that form on one instance (flexible top with
   50% of the mass). The form still consistent with all data is the **rigid-core** statement: *goodness up to norm
   `K` forces the genuine class (Lemma P: the genuine collision lattice) on an initial segment `q <= q_core(K)` of
   the forced range.* For route 3 it would need `q_core(K) >= C_L (8 M^(h_0+h))^L` for some fixed `K`, i.e. a core
   growing polynomially in `M`; the data (`r = 4`, `K <= 4`) cannot say whether it does. Thm A would then be
   applied on the core only.
2. **Mechanism of rigidity (heuristic, labelled).** Rigidity occurs where the genuine labels have nearly used
   the budgets: with small primes `D_s(Z_v)` is mostly smooth over the forced range, the slack of most `(v, s)` is
   `< log q`, and a deformation must place its forced collisions on the few `(v, s)` with slack, which by
   Lemma P pins the class once the short kernel vectors generate the kernel. This is the load accounting of §5
   item 4: the forced range is tight when `log p_j/log r -> 1`, which is also where near-extremal clusters (small
   `W`) would put their primes. A proof would combine (i) a saturation statement for genuine-type budgets and
   (ii) a short-generation statement for kernels of index `n_q <~ |S_h|`. Neither is proved.
3. **Joint rigidity at `K = 4`.** Exhaustive only on `{3, ..., 47}` with the rest fixed (genuine, or one non-genuine
   top). A full joint test of `[3, 103]` with everything free needs an ILP over signatures (the randomized search
   failed to produce any good system, which is evidence for, not proof of, joint rigidity).
4. **Exact forced/free constants.** Theorem 1 brackets the boundary within `r^(floor(K/2))` times constants
   `1/h!` and `2K(K+1)^2 2^h`; the data give `n* ~ 2.5r` (`K = 2`) and `5r-8r` (`K = 3`). The parity pattern at
   `q in Q_1` moves the boundary (e.g. 37 free, 41 forced for `{5, 13, 17, 29}`).
5. **Levels.** Every generated level-1-good system admitted an all-level good lift; a general lemma was not proved
   and is not needed for Theorems 1-2.
6. **Referee items.** Theorem 2(iii) uses `r = Theta(M)` for near-extremal clusters (`round8/stability.md`, as
   summarised in `docs/claude_review/README.md`); Bose-Chowla is cited from memory but only affects constants for
   `K >= 4` (the greedy `B_h` bound suffices qualitatively); `K = 2, 3` use no citation. The `r = 8`, `K = 2`
   coherence sampling (§4.3) was abandoned for run time; `r = 8` enters only through §2 and §3.

## 7. Files (`L4K_checks/`)

| script | output | content | label |
|---|---|---|---|
| `lib.py` | -- | labels, parity, discrete logs in `T_q`, `D_s(Z_v)`, budgets | -- |
| `genuine_check.py` | `genuine_check_out.txt` | collision at `q`, `q^2` with unit `s` iff `q^a | D_s(Z_v)` (209,712 checks) | exact |
| `bh_check.py` | `bh_check_out.txt` | Theorem 1(b) construction, `K = 1..5` | exact |
| `free.py`, `free_full3.py`, `free_K4.py`, `thresholds_table.py` | `free_out.txt`, `free_full3_out.txt`, `free_K4_out.txt`, `thresholds_table_out.txt` | forced/free status of every prime `q <= Q_K` | exact |
| `instance.py` | `instance_out.txt`, `instance_H8000_out.txt` | Theorem 2 finite instances, pair-coherence maxima | exact (per instance) |
| `forced_sig.py` | `forced_sig_2_4_first_out.txt` | all good level-1 systems on the forced primes, `K = 2`, `{5,13,17,29}` | exact |
| `forced_coh.py` | `fc_*.txt` | generated good systems (all levels) and exact coherence tests | exact per system; proportions evidence |
| `deform.py` | `deform_*.txt` | one-modulus deformations of the genuine system | exact (lift column: sampled) |
| `core_sig.py` | `core_*.txt` | all good systems on a core with the rest fixed | exact (level-1 relaxation on the core) |
| `forced_prefix.py`, `joint_search.py` | `prefix_4_4_first_103.txt`, `joint_*.txt` | non-genuine top; joint searches | exact per system; search failures are evidence |
| `tightness.py` | `tightness_out.txt` | aggregate forced mass vs capacity | exact |
| `sigclass_check.py` | `sigclass_2_4_first.txt` | genuine-signature vectors vs genuine class | exact |
