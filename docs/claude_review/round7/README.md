# Round 7: factored residue data (P3-all, P3-gen) and (L_7) (refereed)

Every write-up of this round was attacked by an adversarial referee whose report sits next to it
(`referee_*.md`). **Every statement in §1 is the referee-corrected version.** Anything that no
referee has checked is marked *unrefereed*. That includes the referees' own new sketches, and the
remarks this README adds (marked *README remark*).

The question is the growth rate only. Is there a proof that `M = o(log R/loglog R)` for actual
circles, or a barrier showing that a class of inputs cannot give one? Constants are quoted only
where they measure a barrier.

**Conventions.**
* `N = R^2`, `W = log N = 2 log R`. An arc of length `<= C sqrt R` has angular width
  `Delta = C e^(-W/4)`.
* `w(v) = sum_j |v_j| log p_j`. `d_xy = w(a_x - a_y)`. The pair margin is `mu_xy = d_xy - W/2`.
  `m_xy` is the unit-1 pair collision modulus, `D_xy = 2 log m_xy` the pair demand,
  `D_x = sum_y D_xy` the point total, and `c_C = 2 log C - log 2`.
* *Fake circles* (angle models) have genuine primes, a genuine exponent profile and real angles
  (`round4/outside.md`, `round5/construct.md`).
* **Labels.**
  * *proved*: proved with all quantifiers and refereed;
  * *GRH*: proved under GRH;
  * *model*: proved, but in a probability model, not for actual circles;
  * *exact*: a finite exact computation, independently reproduced;
  * *evidence*: a finite computation that no proof uses;
  * *heuristic*: no proof;
  * *unrefereed*: nobody but the author has checked it.
* Paths of earlier rounds are relative to `docs/claude_review/` in the repository. Round-7 files are
  in this folder.

**Dictionary.** For near-extremal clusters `log M = (1 + o(1)) loglog R`, so `W = c M log M`
corresponds to `M = (2/c + o(1)) log R/loglog R`. A growth improvement means
`W/(M log M) -> infinity` along near-extremal actual clusters. A class is *capped* when it contains
fakes with `W <= c M log M` for some fixed `c`. Then no argument using only its inputs proves
`M = o(log R/loglog R)`.

---------------------------------------------------------------------------------------------

## 0. Answer

**No growth improvement for actual circles was proved this round, either unconditionally or
conditionally on a standard conjecture.** `M <= (2/3 + eps) log R/loglog R` (`two_thirds/proof.md`)
is still the best unconditional bound. The only conditional route known to work is still Vojta's
conjecture for one fixed variety (`growth/conditional.md`), which gives the uniform bound.

**P3-all\*, the open local class singled out in rounds 6 and 8, is still open** (the other open local
classes are subclasses or neighbours of it, §2.1 rows 5-9). P3-all\* is `F^all_X(M, C)`: genuine primes and profile,
prime residues with the norm constraint (hence the forced Legendre parity), and the slab condition
(1.2) of `sharp.md` imposed on every monomial `A_v`. The proved bracket is unchanged:

```text
3 M log M - O(M)  <=  W_min^all  <=  (66 A + o(1)) M log^2 M / loglog M        (X = 3 M^A).
```

The upper end is `sharp.md` Thm 6.4. Every member of the class satisfies the per-unit height
condition `m_(v,s) <= nu_s^(-1) e^(w(v)/2)`. The union-bound completion sum of `sharp.md` §6.1 is
**not** a membership condition. Genuine residue data violate it (sums 4052 and 15190 on `p3G.md`'s
verified prime sets), yet the genuine angles satisfy every slab condition (`p3G.md` Lemma 2.1). It is
a construction tool only (`referee_p3Dlower.md` issue 6).

**What the round did establish** (details in §1, class map in §2, missing lemmas in §3):

1. *The pair channel cannot give growth by itself.* With real weights the pair channel is a cut LP,
   and `2 lambda_max(D) <= W_cut(D) <= pi lambda_max(D)` (proved). The P3-all members of
   `sharp.md` Thm 6.4 have `lambda_max(D) <= 36 M log M + 111 M`. So no certificate that sees only `D`
   through the pair inequalities and the real-weight cut structure proves growth in P3-all (proved).
   What is left at pair level is discreteness: prime weights force `r = O(M)` columns (Question (Q),
   open).
2. *Parity, heights and PSD give constant factors only*, in the precise forms proved by `p3lower.md`
   (B.3, B.5, Cor. A.1). Two of its glosses were false as phrased: "PSD = Perron" (an `O(M)` effect)
   and "parity gains at most a factor 2 per modulus".
3. *Conditional growth on a restricted class* (`p3lower.md` Thm C). Assume `K`-uniform margins and one
   heavy pair (hypothesis `H_1`). Then `M <= (K/(4(1-eta)) + o(1)) log R logloglog R/(loglog R)^2`.
   The referee notes that this is Prop. B.1 applied to an assumed heavy pair: the growth sits in the
   hypothesis. Neither hypothesis is known for any actual cluster.
4. *Construction methods for P3-all\** (`p3Dconstruct.md`, `p3Dlower.md`). Several design families
   provably cannot certify an `O(M log M)` fake:
   * independent random labels on a set of primes of reciprocal mass `(log M)^(-C)`;
   * additive labels;
   * exact Gaussian-coherent labels (under GRH);
   * twisted-coherent labels (certificate caps only);
   * G-coherent designs on the Paley profile, unless their column and pair heights are
     `>= M^(1/2 - o(1))`.

   None of this bounds the class from below. The design suggested in the task (coherent auxiliary
   columns, designed twin/repair columns) is covered by neither theorem and is open.
5. *Genuine residues (route G, `p3G.md`).* The mathematics is correct, but **the headline "P3-G has
   no fakes at `O(M log M)`" did not survive refereeing.** The collapse "P3-G = actual clusters" is
   proved conditional on a complete pinning closure. It is verified exactly on six *random* prime sets,
   and only for members using exactly those primes. Greedily *selected* prime sets defeat every closure
   the paper computes, at the same parameters (referee evidence). Whether P3-G is capped, collapses,
   or gives growth is open.
6. *Discreteness data (`p3data.md`).* For genuine demands, `W_cut/2 lambda_max(D) = 1.00`-`1.05` at
   `12 <= M <= 24` (certified per instance; only evidence about `M -> infinity`). The prime-floor problem is undecided, and the solvers
   are unreliable at `M >= 20`. The claim that adapted profiles are infeasible for their own residues
   is withdrawn: it was an artefact of fixing the windings `k = 0`, which P3-gen leaves free.
7. *(L_7)* is still undecided. It is a global route to the uniform theorem, not a local input class,
   and it does not affect the growth question.

**Correction to earlier READMEs.** `round8/README.md` §3 ("Genuine residues are not a shortcut ...
P3-G there coincides with the actual problem") and the round-8 summary of `p3G.md` rely on the
`p3G.md` headline. They must be read with the scope in item 5: proved on six explicit random prime
sets, and in a non-adaptive model; open for selected prime sets. See §4.

---------------------------------------------------------------------------------------------

## 1. Referee-corrected statements

The classes are nested: `P3-G = P3-gen ⊂ P3-all ⊂ P3-char ⊂ A_X`. Actual normalised clusters
belong to all of them for every `X >= 3` (`sharp.md` Prop. 1.2, re-read in `p3lower.md` Prop. 1.2
and `p3G.md` Prop. 1.3).
* `A_X` is `lower.md` Def. 1.1: point-residue data, with each character's actual collision modulus.
* P3-char and P3-all are `sharp.md` Def. 1.1, with (e-char) and (e-char)+(e-all) respectively.
* P3-gen (`p3lower.md` Def. 1.1) and P3-G (`p3G.md` Def. 1.1) are the same class. It is P3-all with
  genuine residues `r_j(q) = pi_j mod q^(A_q)` for Gaussian primes `pi_j` of norm `p_j`. The angles
  `phi`, and hence `delta`, `theta_0` and the windings `k`, stay free.

### 1.1 `p3lower.md`: the pair channel is a cut LP (refereed: survives, eleven minor issues)

Report: `referee_p3lower.md`. Scripts: `p3lower_checks/`, `referee_p3lower_checks/`.

1. **Theorem A (proved; every member of `A_X`, hence actual clusters, every P3 class, every `X`).**
   * Margins are nonnegative combinations of cuts of total weight `W`.
   * `W >= W_cut(D - c_C(J - I)) >= 2 lambda_max(D) - 2(M-1) c_C^+`.
   * For `B >= 0`: `2 lambda_max(B) <= W_cut(B) <= pi lambda_max(B)`. The lower bound is a Perron or
     hypermetric dual certificate. The upper bound is Goemans-Williamson rounding / Nesterov's `pi/2`
     theorem; the cut cone and hypermetric inequalities are Deza-Laurent. Cite these.
   * *Barrier (A.2).* A certificate that sees only `D`, through the pair inequalities and the
     real-weight cut structure, cannot exceed `pi(lambda_max(D) + (M-1) c_C^-)`.
   * *PSD (corrected).* PSD arguments give `2 lambda_max(D + |c_C|(J - I))` when `C <= sqrt 2`. In
     general they give a doubly-nonnegative value within `2(M-1)|c_C|` of `2 lambda_max(D)`. This is
     not always the Perron bound: the referee's Horn-matrix example at `M = 5` gives `12/5` against 2.
     The effect is `O(M)` and does not change the growth question.
2. **Cor. A.5 (proved barrier).** The `sharp.md` Thm 6.4 members of P3-all (`X = 3M^A`) have
   `lambda_max(D) <= max_x D_x <= 36 M log M + 111 M`. So no `D`-only pair-cut certificate proves
   `W_min^all >= f(M) M log M` with `f -> infinity`.
3. **Elementary bounds (proved, class `A_X`).**
   * `K`-uniform margins: `W >= (2(M-1)/K)(D_max - c_C)` (B.1).
   * Balanced profiles: `W >= 2 max_x D_x - 2(M-1) c_C` (B.2).
   * *Parity (B.3, P3-char and below).* Unit-1 collisions respect the forced Legendre parity. The
     parity-refined pigeonhole raises the bound on `sum D_xy` additively by at most
     `(1/2 + o(1)) M^2 log M`: a factor `<= 2 + O(1/log M)` in aggregate. For a single modulus with
     `M <= m_Q < 2M` the factor is unbounded. With all `p_j > X >= M` it gives `(4 - o(1)) M log M`,
     and `(8 - o(1))` under full bias.
   * *Heights (B.5).* The single-modulus two-term height pigeonhole (`v = +-e_j +- e_l`) gives only
     `r log min(r, X) - O(r)`, which distinct primes already give. Longer monomials, CRT across moduli
     and the parity/size link are **not** covered.
4. **Theorem C (conditional, class `A_X` with hypotheses).** Assume:
   * `X >= M`;
   * `K`-uniform margins;
   * `lambda_Q >= lambda_0`;
   * `H_1`: one pair has unit-1 collisions at `>= (2 - eps/4) log M/loglog M` primes in
     `(M^(1-eta), M]`. This is implied by the factorial-moment Hypothesis `H(eta, delta)`.

   Then `W >= ((8(1-eta) - eps)/K) M log^2 M/loglog M`, so
   `M <= (K/(4(1-eta)) + o(1)) log R logloglog R/(loglog R)^2`. The proof is B.1 applied to the
   assumed heavy pair.
   * The round-6 design `D_A` violates `H` (Lemma C.2).
   * Independent uniform classes satisfy `H` with probability `-> 1` (Lemma C.3). That model ignores
     the forced parity, so it is an `A_X` model, not a P3 model.
   * No actual cluster is known to satisfy `K`-uniformity or `H_1`.
   * The description of `H` through `Im(z_x conj z_y)`-type integers captures units `+-1`, not
     unit 1. Unit 1 also needs `Re = +|B|^2 mod q`.
5. **Discreteness (proved).**
   * Lemma D.0: `W >= (1 + o(1)) r log r`, and `W >= r log X` if all `p_j > X`. So `W = O(M log M)`
     forces `r = O(M)`.
   * Prop. D.1: equal weights need `r >= lambda_max(k)` sign columns, and clique-partition patterns
     attain this.
   * **Question (Q), open.** Can prime weights with `O(M)` cuts absorb the heavy pairs?
6. **Not proved** (glosses withdrawn by the referee):
   * `lambda_max(D) = O(M log X)` for random or actual demands (evidence only);
   * absorption of heavy pairs by actual P3 or P3-gen members: there `D` depends on the profile
     (`kappa = A ell`), so "profiles adapted to `D`" is a statement about the relaxation;
   * "`H` gives nothing without uniform margins";
   * "all height pigeonholes are subsumed by the count of primes";
   * "Theorem L already penalises actual Hadamard-core clusters". This needs the round-8
     qualification: `round8/lw.md` Thm 5.1 shows that the skew double core is immune to every single
     LW form.

### 1.2 `p3data.md`: genuine demands and the discreteness question (refereed: survives as a data paper; two major issues)

Report: `referee_p3data.md`. Scripts: `p3data_checks/`, `referee_p3data_checks/`.

No growth result, no new barrier for actual circles, and (Q) is undecided.

**The problems** (`p3data.md` §3):
* `W_cut` is the cut LP with real weights.
* `W_sc` is the **prime floor**: cut weights in `{0} u [w_0, inf)`, `w_0 = log p_min(X)`.
* `W_prime` uses actual distinct primes `p_j = 1 mod 4`, `p_j > X`, `e = 1`.

Prop. 3.1 (proved) makes `W_sc` a valid relaxation for every member of `A_X` **whose primes all
exceed `X`**. (Q) is the question whether `W_prime(D) = O(M log M)` for the demand matrices `D` that
genuine residues produce.

**Proved, with all quantifiers.**
1. *For every member of `A_X`, hence actual clusters:* Lemma 1.1(1), and Prop. 3.2 (separation):
   `W >= (4/M)[sum_(x<y) D_xy + 2 S_M] - 2(M-1) c_C`, with `S_M = (1/4 + (log 2)/2 + o(1)) M^2`.
   This is an `O(M)` gain over `two_thirds`.
2. *For (Q) and its relaxations:*
   * Prop. 3.1;
   * Lemma 6.1: margins and character margins are additive under concatenation;
   * Prop. 6.2: separation costs `O(M)` at real weights and `O(M log M)` at the prime floor. The
     short-interval prime count it needs is available unconditionally (Huxley-type);
   * Lemma 6.3 with identity (6.1): this is the residual of the Perron dual certificate of
     `p3lower.md` Thm A(2);
   * Cor. 6.4.
3. *Prop. 4.1, a lower bound for equal-weight profiles only:* `W >= min_(tau >= w_0) tau r_P(tau)`,
   from a Plotkin count with ceilings. It does **not** show that equal-weight granularity costs only a
   bounded factor. The strictly equal-weight upper bounds found rise from 1.59 to 1.79 to 1.96 times
   `2 lambda_max` at `M = 12, 20, 24` (medians).
4. *Remark 4.2, at `C = sqrt 2` only:* two `D = 0` pairs whose margins cancel for every weight vector
   make the sign matrix unrealisable by distinct primes.
5. *Example 6.2A, for the relaxation `W_sc` only:* a single collision layer costs
   `W_sc <= 2 lambda_max(D) + 2 M w_0`. Its profiles always contain a Remark 4.2 tie, so they are
   **not** realisable by distinct primes as built.
   * The referee sketches a realisable variant, using a tie-breaking Sylvester block and short-interval
     primes, with `W_prime <= 2 lambda_max + O(M log M)` (*unrefereed sketch*).
   * "The whole difficulty is merging `~ X/log X` layers" is interpretation, not a theorem.
6. *Theorem 5.1 and Cor. 5.2, in the model `U_X`:* squarefree profile, arbitrary real weights,
   windings `k` fixed, labels independent and uniform on `(Z/m_q*)^r`. Then every profile with
   `W <= (2 - eps) M log^2 M/loglog M` satisfies all unit-1 pair inequalities with probability
   `<= 4 M^(-eps/16)`.
   * This holds only for windings fixed independently of the residues.
   * If the windings may be chosen after the residues, the referee sketches that the star argument
     survives with `G(t) = P(min over the four units of D > t)`, at a constant about 4 times smaller
     (*unrefereed sketch*).
   * *README remark:* `U_X`, like Lemma C.3 of `p3lower.md`, draws labels without the forced Legendre
     parity. So it is an `A_X`-type model, not a P3 model.
   * This extends `sharp.md` Prop. 6.5, which should be cited.

**Exact per instance, but only evidence about `M -> infinity`.** These are genuine unit-1 demands of
the fully flipped Paley profiles with actual primes.
* The demands match the independent residue model (`M <= 84`). `D_max/H` runs from 2.2 to 3.2 at
  `X = 2M` (2.9 to 3.7 at `X = 4M`), with `H = log^2 M/loglog M`.
* `W_cut/2 lambda_max = 1.00`-`1.05` for `12 <= M <= 24` (certified, gap `<= 4.3e-10`).
* `LB_eq/2 lambda_max` has median 1.12-1.22 for `M <= 32`. This is a lower bound, so it shows that
  the Plotkin obstruction is a flat factor and nothing more.
* Prime-floor upper bounds, as medians of `W_sc/2 lambda_max`, are 1.31 and 1.21 at `M = 8` (all 12
  instances solved to optimality by the referee's MILP) and 1.37-1.82 at `M = 12`-24. Local search
  misses known Plotkin-tight codes at `M = 20, 24`, so these upper bounds say nothing about the
  optimum there.
* 64 actual-prime profiles were verified exactly. Every character minimum (`|c_x| <= 12`) is attained
  at a pair.
* **Withdrawn:** "the fitted profiles are infeasible for their own genuine residues". With the
  windings free, as P3-gen allows, 23 of the 24 fitted `M = 8` profiles are pair-feasible for their
  own residues for some conjugation choice. The separation profiles are feasible in a median of 23% of
  the 400 choices, and in up to 90%.

**Open.**
* Do adapted near-tight `+-1` designs with `O(M)` prime columns exist for genuine `D`?
* A growth proof through pair data needs an input on how actual residues relate to the exponent
  profile, now including the windings.

### 1.3 `p3G.md`: route G, genuine residues pin the angles (refereed: mathematics correct, headline does not survive)

Report: `referee_p3G.md`. Scripts: `p3G_checks/`, `referee_p3G_checks/`.

Class: `P3-G_X` (genuine residues, every monomial strengthened). It contains every actual cluster for
every `X`, with `phi = theta`. If `X' >= X` then `P3-G_(X') ⊂ P3-G_X`.

**Proved.**
1. *Mechanism.*
   * `m_(v,s) = gcd(D_s(A_v), L'_X)`, where `D_0, ..., D_3 = Y_v, X_v - Y_v, X_v, X_v + Y_v` for
     `A_v = X_v + i Y_v`.
   * A **tight** `(v, s)`, meaning `m_(v,s) = |D_s(A_v)|`, is exactly the one-sided constraint
     `<v, phi - theta> in sigma[0, pi - 2d] + pi Z`, with exact integer data `sigma` and `sin^2 d`.
     The genuine angles `theta` are a planted solution of every monomial condition.
   * Circuits of tight relations with arc sum `< pi` pin `<v_i, phi - theta>` to `pi Z`.
   * The certified closure `Pi*(pi, X)` depends only on the rational prime set and `X`.
   * *Referee Lemma R.1 (proved in the report; not re-refereed).* Relative to `{0}`, arc sums of zero
     combinations lie in `(pi/2) Z`. So level-0 certificates have arc sum exactly `pi/2`.
2. *Collapse (Thm 3.1).* Suppose `2 e_j in Pi*(pi, X)` for every column `j`. Then every member of
   `P3-G_(X')`, `X' >= X`, whose prime set is **exactly** `{p_j}` is an actual cluster after
   `pi_j -> i^(t_j) pi_j`. This holds for any `M`, `C`, exponents and rows.
3. *Dichotomy (Thm 3.2, Prop. 3.3, Cor. 3.4; squarefree members, `X >= 2M`).* Every character that
   vanishes on the unpinned columns is a genuine near-relation. Since `rank V >= (1/2 - o(1)) M`, a
   squarefree member has at least `(1/2 - o(1)) M - k` unpinned columns, or `k + 1` independent
   genuine near-relations.
4. *Upper bound (Prop. 6.1).* `P3-G_X` has non-actual members with `W <= (66 + o(1)) M X` for every
   `X >= 2M` (the paper's constant 68.7 is valid but not sharp).

**Exact (independently re-verified).** The six random prime sets of Thm 4.1 (primes of size `M^3` and
`M^4`, `M = 44`-644, `X = 2M`) have complete closures. So `P3-G` restricted to members with exactly
these prime sets consists of actual clusters.
* Members on **proper subsets** are not covered. 163-subsets of the verified `M = 164`, `b = 8`, `K = 4`
  set leave 35-54 columns free.
* The round-5/6 constructions are not covered either: at `M = 44` they need 1419 primes, and the
  verified set has 86.
* The genuine union-bound sums on these sets are 4052 and 15190. This is how we know that the union
  bound is not necessary for membership.

**Model (`p3G.md` §5).** Assume independent tightness, **non-adaptive** primes, the heuristic count of
`O(M/log M)` near-relations, and squarefree members. Then columns with `log p <= (1/2 - eps) H_X` are
pinned with high probability, and every squarefree **non-actual** member has
`W >= (1/4 - eps - o(1)) M H_X`, where `H_X = log X log M/loglog M`. Actual clusters are untouched by
this bound.

**Heuristic.** Fakes at about `(66 + eps) M H_X`.

**Withdrawn or open.**
* *Withdrawn:* "Route G does not give fakes at `O(M log M)`", "prime sets of polynomial size are
  complete", and "a fake needs `~ M/2` columns of size `c H_X`" as established statements.
* *Withdrawn:* the §4.2 gloss that Thm 4.1 covers Paley profiles with `b(M-1) < r` and transplanted
  round-5/6 constructions.
* *Unproved:* "P3-G-char does not collapse". What is shown is only that the §2 mechanism is
  unavailable there.
* *Referee evidence, exact for each stated sub-closure.* Greedily selected prime sets have an **empty**
  computed closure where random sets collapse:
  * at `M = 644` and `M = 1284` with `K = 4`;
  * at `M = 164` down to `K = 2.75` (`lambda/H = 0.88`, against the measured 1.29 for random sets).

  "Empty" covers the paper's `k <= 2` closure, all closed walks of singles and pairs, and, after at
  most 3 column deletions, all certificates with `<= 3` terms among weight-`<= 2` relations.
  Weight-3 relations re-pin these sets for `K <= 3.5`. Pruning 4-15 more columns leaves 148-159 of
  163 columns free in a weight-`<= 3`, `k <= 3` closure at `K = 3.0`-3.5, where random sets are fully
  pinned. The full closure (all weights) on selected sets has **not** been computed.
* The measured thresholds `lambda* ~ (1.23-1.34) H` are properties of random sets. The premise of
  Heuristic 5.3(3), "an independent set in the tight graph", is the wrong condition: builders keep
  pair degree 26-82 with no certificates.
* **Status:** it is open whether P3-G has non-actual members at `W = O(M log M)`. Route G neither
  decides P3-all\* nor says anything about actual circles beyond `two_thirds`.

### 1.4 `p3Dlower.md`: can designed residues be forced to cost `>> M log M`? (refereed: survives as a modest partial barrier)

Report: `referee_p3Dlower.md`. Scripts: `p3Dlower_checks/`, `referee_p3Dlower_checks/`.

Class: P3-all\* `= F^all_X(M, C)`, with the slab conditions (1.2) for every `v`, hence the per-unit
height condition. The union-bound completion sum is a construction tool, not a membership condition
(see §0). With that definition actual clusters are members. No lower bound `W >= f(M) M log M` with
`f -> infinity` is proved for P3-all\* or for any natural subclass. The bracket of §0 is unchanged.

**Proved unconditionally.**
1. *Parity localisation (Thm 1.1).*
   * At `q = +-3 mod 8`, unit shifts make any label vector parity-respecting. Collision sets and the
     all-unit moduli `m_v` are unchanged, and the units are permuted by `<v, u(q)>`. The parity does
     fix the parity of the collision unit, `s = <v, nu(q)> mod 2`.
   * At `q = +-1 mod 8`, a collision needs `<v, nu(q)>` even. If `gcd(v, m_q/2) = 1`, a uniform
     parity-respecting label collides with probability `8/m_q` or 0, against `4/m_q` without the
     parity.
   * *Withdrawn:* "at `q = +-1 mod 8` the parity acts only as a veto". The parity forces torsion-type
     collisions, for example `(m_q/8) e_j` when `(p_j/q) = 1`. At `q = 7` this is the norm-1 monomial
     `e_j` (exact on 553 actual primes). They are harmless for heights.
2. *Twin decomposition (Thm 2.1, Cor. 2.2).* This restates `sharp.md` Lemma 6.3(2) and §6.6, which
   should be cited. On twin profiles, in the basis `(ell_J, U, Lambda)`:
   * point and pair data depend only on `kappa = F(ell_J, U) + Lambda`;
   * every parity-compatible `kappa` occurs with every parity-respecting `(ell_J, U)`.

   So pairs act on monomials only through the `M` twin labels `Lambda`. **This localises (M1) to
   (M1'); it does not refute it** (§3, item L6).
3. *Additive labels (Prop. 3.1).*
   * On a set `Q` of moduli, additive labels force `W >= (r/(8k))(log Pi_Q - log 2)` whenever
     `binom(ceil(r/2), k) > 2kH + 1`.
   * At `q = +-1 mod 8` they are parity-admissible only for exceptional prime sets. Under GRH this
     never happens once `X >> log^2 p_max`.
   * The admissible version uses `Q ⊂ {q = +-3 mod 8}` and gives `W >= (1/8 - o(1)) M X/(h+2)`.
4. *Twisted coherent labels (Def. 4.1 - Cor. 4.4).*
   * They are parity-respecting residue data, with the explicit cap (4.2), quadratic in `k` and
     valid at every level. The referee checked this at every prime-power level up to `X = 200`:
     20,398 collisions.
   * The cap certifies heights only for `k <= k_1 = O(tau/(T(X) log M))`.
5. *Lemma 5.1.* At one modulus, fewer than `m_q/8` monomials with `gcd(v, m_q/2) = 1` can be avoided.
   This is nearly vacuous for the question.

**Under GRH (Thm 3.2, Cor. 3.3).**
* Assume coherence with `log(zeta/conj zeta)` (up to `mu_4`) at all `q = +-1 mod 8`, `q <= X`, with the
  forced parity. Then `Norm zeta >= p` once `X >= c_1 (log 2n)^2 (loglog 3n)^2`. This is the standard
  GRH/Chebotarev pseudosquare argument.
* Given coherence at the top level and non-degeneracy, the crude certificate `m <= |D|` then has zero
  slack and cannot prove the union-bound completion.
* This does **not** show that such designs, or genuine labels, fail to be members at `O(M log M)`.
  Genuine labels are certified with zero slack, so "certified designs cannot reach `O(M log M)`" read
  literally would prove growth for actual circles. It is withdrawn.

**Heuristic.**
* "(M2): pigeonhole or entropy cannot force a height violation" is a heuristic. It is the remark of
  `sharp.md` §6.7 that counting forces only averages.
* Every analysed design needs `tau >~ log^2 M/loglog M`. This is:
  * proved for the random fibre (constant between 2 and `2A`) and for additive designs;
  * proved at the certificate level under GRH for exact coherence;
  * heuristic for twisted coherence and for `D_A` with coherent free labels;
  * for genuine labels, resting only on random prime sets or a model (`p3G.md`, §1.3).

**Evidence.** The draft's §6 table is withdrawn: `designs.py` mixed the units `i` and `-i` across
moduli. The corrected values are:
* genuine labels: about `1.6 k log M`, at their own cap and far below random;
* coherent labels: about `0.8 k log M`;
* twisted labels: 2-3.4 times coherent at `k = 1`, near random at `k = 2, 3`.

No rate is read off.

**Roadmap (corrected).** (I1) inverse step, (I2) certified implies genuine size, (I3) genuine size
implies Poisson tails. (I2) is proved only for exact coherence under GRH. It is **false** for twisted
coherence, which is certified with `Norm zeta ~ X << p_j`. So the roadmap needs a separate branch for
approximately certified designs (§3, items L5-L7).

### 1.5 `p3Dconstruct.md`: structured fibres for P3-all\* (refereed: survives as a partial barrier for construction methods)

Report: `referee_p3Dconstruct.md`. Scripts: `p3Dconstruct_checks/`, `referee_p3Dconstruct_checks/`.

No construction was found, and nothing is claimed about actual circles. The best rate is still
`sharp.md` Thm 6.4. All results below are barriers to **methods of constructing** fakes. None of them
is a lower bound for the class.

1. **Theorem 1, the Poisson barrier (proved, after the repair (Ind') = unconditional independence).**
   * *Setting.* Fix the primes first. Let `Q* ⊂ [M^(1-eta), X]`, `X <= M^(A_0)`. Suppose the labels of
     a set `J` of columns are independent and uniform in their parity classes on `Q*`, and that
     `>= eps M` columns of `J` have `lambda_j >= (log M)^(-C)`.
   * *Conclusion.* For `M >= M_0`, the design completes with probability `< exp(-M^(eta/4))` unless
     `W >= (eps(1-eta)^2/(2(C+1)) - o(1)) M log^2 M/loglog M`.
   * *Scope.* This is a probability statement, not "every completion". It does not cover
     deterministic designs, column-dependent laws, or weights chosen after the labels.
   * *Counterexample to the summary.* Labels uniform on the parity class **minus the 4-torsion** have
     no single-column collisions. So "independent randomness on any prime set is fatal" is false as
     written.
   * *Referee sketch (unrefereed).* Two-column monomials `e_j - e_k` repair this. They replace
     `lambda_j` by the reciprocal mass of `Q*`, and uniform labels by i.i.d.-within-parity labels, at
     constant `eps/12`.
   * The original conditional (Ind) fails in the write-up's own applications (`sharp.md`'s random
     fibre, the twin columns). The proof uses only unconditional independence.
2. **Root-of-unity parity repair (Lemma 3.2, Prop. 3.3; proved, exact check on 89 prime powers
   `<= 400`).**
   * Repair by a generator of the 2-Sylow subgroup gives `q^a | Im(gamma_v^(2^k(q)))`. So
     size-certified coherent designs pay a factor `2^(K+1)` on the classes `k(q) <= K`. On
     `q != +-1 mod 16` the parity costs only a constant. Prop. 3.3(2) and Def. 4.1 need **split**
     Gaussian primes.
   * *Exact.* Matching the parities through `Norm zeta_j in {1, 2} p_j Q*^2` forces `|zeta_j|^2 >= p_j`,
     so certification by size is lost.
   * *GRH.* Any other matching of a whole 2-adic class on `[M^(1-eta), X]` is impossible.
     Unconditionally it is excluded only polynomially.
   * *Withdrawn:* "matching costs `2^K >~ X/(log X log p_j)` (exact)". It is a pseudosquare density
     heuristic.
   * *Open:* whether bounded-cost repairs exist on `k(q) >= 4` (Heuristic H1).
3. **Theorem 2, the Hadamard height barrier (proved).**
   * Take G-coherent designs on the `sharp.md` §4 Paley profile under (4.1), with `C` the maximum
     column height and `Hbar` the average pair height. Then
     `W >= (M-1)(theta(X; 8, 3) - 2 log 2)/(8 ceil(8(C + 2 Hbar)))`.
   * So at `W = O(M log M)`, `C + 2 Hbar >= M^(1/2 - o(1))`. The referee sharpens this to
     `Hbar >= M^(2/3 - o(1))` when `C = O(log M)`.
   * Pairs and columns cannot both be certified by size. Step 3 of the proof is the Donoho-Stark /
     Elad-Bruckstein `l1` support uncertainty.
4. **Identity (1.4)** restates `sharp.md` Lemma 6.3(2) and Cor. 2.2(5): the twin labels are solved
   coordinates, `t_x = eps_x(kappa_x - P_x)`.
5. **The design suggested in the task is not covered** (coherent auxiliary columns, designed twin or
   repair columns, structured point classes). Its twin labels are deterministic. They are either not
   coherent (e.g. `D_A` point classes) or coherent only with height `>= M^(1/2 - o(1))`, so neither
   Theorem 1 nor Theorem 2 applies. Its fate is **open**. It is blocked only heuristically, at the
   same rate (H2, pair-of-columns Poisson tails).
6. *Scope corrections.*
   * "Any construction at `O(M log M)` must control the smooth parts of `binom(M,2)` integers of size
     `exp(M^(1/2-o(1)))`" holds only for G-coherent Paley designs, and only for the *average* pair
     height.
   * Consequence 4's surjectivity hypothesis holds on `M/2` rows for every Paley order tested
     (`q0 <= 107`, exact SNF), but the consequence is redundant.
   * "Prop. 6.5 is the case `C = 0`" holds in rate only (factor 4 lost).
7. *Evidence.* Single columns at `M <= 2048`: coherent max about `7.3 log M`, random about `11 log M`.
   Coherent pair demands look like random ones, with `lambda_max(D) ~ 0.94`-`0.96 M log M`. Maxima vary
   by about `2 log M` across seeds. At these sizes the data cannot separate `log M` from
   `log M/loglog M`.

### 1.6 `L7.md`: (L_7), second attempt (refereed: survives as a partial barrier with smaller scope)

Report: `referee_L7.md`. Scripts: `L7_checks/r8/`, `referee_L7_checks/`. The first attempt is kept as
`L7_prev_attempt1.md` (cited as r7a).

**Background.** (L_7) is the arithmetic lemma of `threshold/L8.md` §1 and `growth/ptolemy.md` §7 on
`Mbar_{0,7}`. It is a global (Diophantine) route to the **uniform** theorem, not a local input class.
Round 4 Thm 2(a) (refereed) gives a refutation criterion: (L_7) is false if good curves whose classes
lie within `eps` of the balanced ray have Zariski-dense union, for every `eps > 0`.

**(L_7) is still not decided.**

1. **Theorem A (proved).** Seven planes in `P^3` over `Q` satisfying (A1) and (A2) give a
   `Q`-birational map `G(2,4) ⇢ Mbar_{0,7}`, `R -> (R ∩ H_a)_a`. The literature was not checked; it is
   probably a special case of the Gale / Gelfand-MacPherson correspondence.
2. **Theorem B (proved; the class is exact Data; the referee added irreducibility of the parameter
   space).**
   * The rulings of the web quadrics through six cube points of the star arrangement form an
     irreducible 3-dimensional `Q`-family `V_A`, with exactly one member through a general point.
   * Its general member is a good curve of class `gamma_A`: pair 6, triple 12, `(B+8K) = +2`.
   * So rich good `Q`-curves are Zariski dense.
3. **Proved (Lemma 4.2, Prop. 4.4, Cor. 4.5).**
   * Let `C'` be a good curve avoiding `Sing(F_V)`. Suppose `C'` is a section of a `V`-ruled surface
     that moves, or that meets every fibre in the same component as the base section. Then
     `lambda_V(C') >= 0`.
   * So r7a's (G1) fails for good curves `C` avoiding `Sing(F_V)` with `lambda_V(C) < 0`. r7a's
     Remark 1 is false, and its Theorem 3 is only a conditional statement. It is not "incorrect".
4. **Prop. 5.4 (proved; exact computation, independently confirmed).** `lambda_A(beta) <= -23` and
   `lambda_A(mu) <= -10`.
   * *Evidence:* equality at -23 (three members give identical values on all 1260 classes); other
     event sets give -31 and -32.
5. **Consequence (proved).** No approximately balanced good curve **avoiding `Sing(F_A)`** is a moving
   or same-component section of a `V_A`-ruled surface. The same holds for its `S_7`-translates and
   for the forgetful families, where avoidance is automatic.
   * For `V_A`, avoidance is **not** automatic. `Sing(F_A)` contains the point `x_m`, which lies on
     `D_{O b1 b2}` alone (referee Data).
   * Whether passages through `Sing(F_A)` can make up the `23L` deficit is open. The heuristic answer
     is no.
6. **Not proved.**
   * Prop. 5.2 as stated: it needs the leaves and the images of the ruling lines to avoid
     `Sing(F_{V_k})`. A 2-dimensional analogue gives `lambda = 1`, not 2.
   * "`lambda(gamma_V) = 2` in general". It is true for `V_A` and for forgetful fibres.
   * Anything about round 4 Thm 2(b) beyond birational `R = 4` charts.
   * "`V_A` cannot finish ANY addition or section construction".
   * The §9.1 conditional statement (rational connectedness plus Colliot-Thélène implies (L_7) false)
     skips the smooth model, the field-of-moduli twist, and openness of goodness.
7. **Meaning.** Nothing changes for the uniform theorem, for growth, or for actual circles.

### 1.7 Round-8 inputs used by this round

`round8/lw.md` and `round8/weakest.md` were refereed in round 8 (`round8/README.md` §1, with
`referee_lw.md` and `referee_weakest.md`). They are therefore quoted below as refereed, not as drafts.
* *Single-form Lang-Waldschmidt.* No growth, in either the max or the product form, for any
  `eps >= 0`. The Paley skew double core of `round5/exact.md` is immune to every single LW form
  (`lw.md` Thm 5.1). This qualifies `p3lower.md`'s remark on `flipped.md` Theorem L.
* *Simultaneous Lang-Waldschmidt.* The `SLW^log` hypothesis of `weakest.md` is false for genuine
  algebraic numbers. The repaired, coefficient-weighted form is new and unsupported.

The point-level residue version of the double core still rests on the unrefereed `exact.md` Thm 5.1.
No one has given the double core P3-char or P3-all data. That is the "decisive first test" of
`round8/README.md` §3, and **this round did not carry it out**.

---------------------------------------------------------------------------------------------

## 2. Status of every local input class

A *local input class* is a set of conditions that every normalised actual cluster satisfies and
that are read off the cluster's archimedean and `q`-adic data. Status is one of three:
* **capped**: the class contains fakes at `W = O(M log M)`;
* **open**;
* **gives growth**: a proved implication to `M = o(log R/loglog R)`, possibly from extra hypotheses,
  which are then named.

Larger input sets give smaller classes, so a cap on a class says nothing about its subclasses.

### 2.1 Classes that contain every actual cluster

| # | class (inputs) | status | witness or reason | source |
|---|---|---|---|---|
| 1 | `A_0`: exponent profile, pair chords, cut/Plotkin, all monomial (Liouville) gaps, rational characters, GAP/Fourier, energy, incidences, Haar-generic angles | **capped** | skew double core at `(8 + o(1)) M log M + 8 M log(1/C)` | `round8/lw.md` Thm 6.1 (proved, uses `exact.md` Thm 4.1(iii), re-derived by `referee_lw.md`) |
| 2 | `A_X = A_1(X)`: `A_0` plus point-level residues at moduli `<= X`, each character with its actual collision modulus; includes every round-7 inequality valid on `A_X` (Thm A, B.2, separation Prop. 3.2) | **capped** for `X = M^K` | `construct.md` Paley cores at `(27 + o(1)) M log M` (`X = 2M`), constant `9(2 + max(K,1))` at `X = M^K`; `lower.md` `(16 + 32A)` at `M^A` | round 5 (proved); the round-7 inequalities hold for all members, hence for the witnesses |
| 3 | `A_X` plus round-8 structure (sub-cluster gcds, smoothness, prime counts, Plotkin balance, collision and level stability) | **capped** | `construct.md` models via transfer (`stability.md` Prop. 7.1) | round 8 (proved; renormalisation after row deletion still to be written out) |
| 4 | **P3-char** `= F^char_X`, `X <= 3 M^A`: residues factor through prime residues with the genuine norm, so the forced Legendre parity holds; imposed on characters | **capped** | `sharp.md` parity-class designs at `(b max(A,2) + o(1)) M log M`, `b >= 33` | round 6 (proved) |
| 4a | P3-char with `K`-uniform margins | **capped** | the same designs are `K`-uniform, `K <= (b+4)/(b-4) + o(1)` (`p3lower.md` §6) | proved |
| 5 | **P3-all\*** `= F^all_X`, `X = 3 M^A`: P3-char plus the slab condition (1.2) for **every** monomial `v in Z^r` (so the per-unit height condition holds). The union-bound completion sum is **not** part of the class | **open** | bracket `3 M log M - O(M) <= W_min^all <= (66A + o(1)) M log^2 M/loglog M`; with all `p_j > X >= M`, B.3 raises the lower end to `(4 - o(1)) M log M` | `sharp.md` Thm 6.4; `p3lower.md` B.3(4); §1.1-§1.5 |
| 5a | P3-all\* with `K`-uniform margins | **open** | same bracket: the Thm 6.4 members are `K`-uniform with `D_max = O(log M)`, so B.1 gives only `O(M log M)`. Their cost comes from monomials | `p3lower.md` §6 |
| 5b | `A_X` (or any subclass) with `K`-uniform margins, `lambda_Q >= lambda_0`, `X >= M`, and `H_1` (one pair with `>= (2 - eps/4) log M/loglog M` unit-1 collisions at primes in `(M^(1-eta), M]`) | **gives growth** (conditional): `W >= ((8(1-eta) - eps)/K) M log^2 M/loglog M`, `M = O(log R logloglog R/(loglog R)^2)` | Prop. B.1 applied to the assumed heavy pair. No actual cluster is known to satisfy `K`-uniformity or `H_1`; the designed members of rows 4a and 5a violate `H` | `p3lower.md` Thm C (proved, conditional) |
| 6 | **P3-G = P3-gen**: P3-all\* with genuine residues `r_j(q) = pi_j mod q^(A_q)`; angles and windings free | **open** | non-actual members at `W <= (66 + o(1)) M X` for every `X >= 2M` (proved); heuristic fakes at `(66 + eps) M H_X`; lower end `two_thirds` only. On the six verified random prime sets the class **is** the actual problem (collapse); selected sets evade every computed closure | `p3G.md` Prop. 6.1, Thm 3.1, Thm 4.1; `referee_p3G.md` §4 |
| 7 | P3-G-char: genuine residues, characters only | **open** | not shown to collapse; at pair level it is (Q) plus `H` | `p3G.md` §7, `referee_p3G.md` issue D |
| 8 | P3 classes at super-polynomial moduli (`X` up to `N^A`) | **open** | no fake is constructed there; one a-priori bound is useless at `X >= N^(1/4)` | `sharp.md` §7 |
| 9 | additive (archimedean) data, e.g. `Log((z_x + z_y)/(2 z_u)) = O(Delta)` | **open** | not modelled by angle models; natural uses reduce to single-place auxiliary polynomials (barred) | `round8/README.md` §2.5 item 5 |

Rows 5-9 (and 5a) are the only valid local classes not known to be capped. A growth theorem for any of
them is a growth theorem for actual circles. Row 6 is not a shortcut: where the closure is complete it
coincides with the actual problem.

### 2.2 Methods, relaxations and models inside the open classes

These rows are not classes of circles. They say which **arguments** can or cannot prove growth in
row 5 (lower-bound methods), which **designs** can or cannot cap it (construction methods), and what
happens in probability models.

| # | method / relaxation / model | status w.r.t. growth | source |
|---|---|---|---|
| M1 | pair-cut certificates: pair inequalities plus cut/Plotkin/PSD structure, real weights, seeing only `D` | **capped** in P3-all: no such certificate exceeds `pi(lambda_max(D) + (M-1) c_C^-)`, and the Thm 6.4 members have `lambda_max(D) <= 36 M log M + 111 M` | `p3lower.md` Thm A, Cor. A.2, A.5 (proved) |
| M2 | parity pigeonhole (B.3), single-modulus two-term height pigeonhole (B.5), PSD | **constant factors only**: B.3 adds `<= (1/2 + o(1)) M^2 log M` to `sum D_xy`; B.5 gives `r log min(r, X) - O(r)`; PSD is within `O(M)` of Perron | `p3lower.md` (proved) |
| M3 | pair channel with discrete prime weights (prime floor `W_sc`, actual primes `W_prime`; Question (Q)) | **open**. Relevant only for demand matrices that residues cannot make small, i.e. genuine or `H`-type `D`: designed residues in P3-char already make `D_max = O(log M)` | `p3lower.md` §5, `p3data.md` §3-§6 |
| M4 | separation (`delta`-placement) | **`O(M)` only**: `W >= (4/M)[sum D + 2 S_M] - 2(M-1) c_C`; costs `O(M log M)` at the prime floor | `p3data.md` Props. 3.2, 6.2 (proved) |
| M5 | monomial covering: the slabs `{|sin(<v,phi> - s pi/4)| < nu_s m_(v,s) e^(-w(v)/2)}`, `v` outside `L`, covering the torus of free angles | **open**. No necessary condition beyond the heights is known; counting forces only averages (heuristic, `sharp.md` §6.7) | `p3lower.md` §7 item 4; `p3Dlower.md` §5 |
| M6 | genuine-residue pinning (tight relations, closure `Pi*`) | proves collapse on verified prime sets (exact). It bounds only **non-actual** members, so it can show that P3-G is not capped, never growth. Its bound `W >= (1/4 - eps - o(1)) M H_X` for squarefree non-actual members holds only in the non-adaptive model `T`, with the heuristic near-relation count. **Open** for selected primes | `p3G.md` §2-§5; `referee_p3G.md` §4 |
| D1 | uniform random fibre; independent labels uniform in parity classes on primes of reciprocal mass `(log M)^(-C)` | cannot cap row 5: needs `W >~ M log^2 M/loglog M` (with probability `1 - exp(-M^(eta/4))`) | `sharp.md` Prop. 6.5 (per-unit repair); `p3Dconstruct.md` Thm 1 with (Ind') (proved) |
| D2 | additive labels | cannot cap: `W >= (1/8 - o(1)) M X/(h+2)` | `p3Dlower.md` Prop. 3.1 (proved) |
| D3 | exact Gaussian-coherent labels with the forced parity | crude certificate has zero slack (GRH); membership at `O(M log M)` **not** excluded | `p3Dlower.md` Thm 3.2, Cor. 3.3 |
| D4 | twisted-coherent labels | certificate caps quadratic in `k`; certifies only `k <= O(tau/(T(X) log M))`; rate `log^2 M/loglog M` is heuristic | `p3Dlower.md` §4 |
| D5 | G-coherent designs on the Paley profile | need `C + 2 Hbar >= M^(1/2 - o(1))` at `O(M log M)` | `p3Dconstruct.md` Thm 2 (proved) |
| D6 | coherent auxiliary columns plus designed twin/repair columns plus structured point classes (the task's suggestion) | **open**; heuristic H2 says same rate | `referee_p3Dconstruct.md` M1 |
| D7 | skew double core with P3-char or P3-all data | **open, untested** (the round-8 "decisive first test") | `round8/README.md` §3 |
| P1 | model `U_X`: labels uniform and independent of the profile, windings fixed | **gives growth in the model**: every profile with `W <= (2 - eps) M log^2 M/loglog M` fails a pair inequality with probability `>= 1 - 4 M^(-eps/16)`. Windings chosen afterwards: referee sketch, constant about 4 times smaller (unrefereed). *README remark:* no forced parity, so an `A_X`-type model | `p3data.md` Thm 5.1 (model) |
| P2 | independent uniform classes (no parity) | satisfy `H(eta, delta)` with probability `-> 1`, hence growth via Thm C **for `K`-uniform profiles only** | `p3lower.md` Lemma C.3 (model) |

**Reading.**
* Lower-bound methods M1, M2 and M4 are capped.
* M3 and M5 are open. Each needs an arithmetic input about **actual** residues: `H`, or a statement
  that the residues cannot follow the profile. M6 can decide whether row 6 is capped, but by itself it
  cannot give growth.
* Construction methods D1-D5 cannot cap P3-all\*. D6 and D7 are untested or open.
* In the models, growth appears: for every profile in P1, and for `K`-uniform profiles in P2. That
  is exactly why the missing input is a statement that actual residues behave like the model relative
  to their own profile.

### 2.3 Non-local inputs (pointer)

These rows are unchanged from `round8/README.md` §2.2-§2.4 and are not repeated here:
* capped: single-form Lang-Waldschmidt (and the fixed-`eps` max form with `i` listed), Baker/Matveev,
  single-place auxiliary polynomials;
* refuted: `SLW`, `SLW^log`;
* no growth: abc, Bombieri-Lang, the Subspace Theorem with known `S`-dependence;
* gives the uniform bound conditionally: Vojta's conjecture for one fixed variety.

(L_7) (§1.6) is a global route to the uniform theorem and is still undecided.

---------------------------------------------------------------------------------------------

## 3. The precise missing lemmas

Throughout:
* "near-extremal" means `W <= c_0 M log M` for a fixed `c_0`;
* `X = 3 M^A` for a fixed `A >= 1`;
* by `stability.md` Cor. 3.6 (round 8, proved), every lemma below is needed only for profiles whose
  primes are all `<= M^(O(1))`.

Each item says what it would give.

**The target.**
* **(T)** For some fixed `A >= 1`, every member of `F^all_(3M^A)(M, C)` has `W >= f(M) M log M` with
  `f(M) -> infinity`. This gives `M <= (2/f(M) + o(1)) log R/loglog R` for actual circles. Every route
  below is a route to (T), or to its negation.

**Route 1: pairs, with uniform margins** (`p3lower.md` §4, §7). L1 and L2 together give growth
through Theorem C.
* **L1 (uniformity, U).** A valid input forcing near-extremal actual clusters to have
  `max_(x,y) mu_xy <= K mubar`, with `mubar` the mean margin and `K` fixed.
  * The alternative form: force `lambda_max(D) >> M log M`. This needs irregular point totals,
    `max_x D_x >> M log M`, by Cor. A.3.
  * *Status:* open. Theorem A shows that the pair inputs themselves allow adapted margins at real
    weights. Pigeonhole and parity force only averages.
* **L2 (`H_1` for actual clusters).** There are fixed `eta, eps > 0` such that every near-extremal
  actual cluster has a pair `(x, y)` that collides with unit 1 at `>= (2 - eps/4) log M/loglog M`
  primes `q in (M^(1-eta), M]`.
  * Concretely: the reduced chord `(zhat_x - zhat_y)/g_xy` of the normalised points, a nonzero
    Gaussian integer of modulus `<= 2 e^(d_xy/2)`, is divisible by that many window primes. This is
    exactly the unit-1 condition (`p3lower.md` §1.1; checked on 215,670 instances in
    `referee_p3lower.md` §3.3). The `Im(z_x conj z_y)` form captures units `+-1` only.
  * *Status:* open. `D_A` violates `H`, so L2 has to use the genuineness of the residues.
  * It is a joint smooth-part statement about `binom(M,2)` dense Gaussian integers.

**Route 2: pairs, with discrete weights** (`p3lower.md` §5, `p3data.md`).
* **L3 (prime-floor lower bound, (Q) answered "no").** For every demand matrix `D` that an actual
  cluster can produce, every profile with distinct primes `p_j > X` satisfying `mu >= D - c_C` has
  `W >= f(M) M log M`.
  * It cannot hold for all `D`. A single collision layer is cheap: at `W_sc` by Example 6.2A
    (proved), and at `W_prime` by the referee's realisable variant (unrefereed sketch).
  * So L3 needs an `H`-type genericity input as well.
  * *Status:* open. The evidence at `M <= 24` is inconclusive: the solvers miss known Plotkin-tight
    codes at `M >= 20`.

**Route 3: monomials** (`p3Dlower.md` §5.3, `p3lower.md` §7 item 4).
* L4 + L6 + L7, together with L5 on twin profiles, would give growth for the near-uniform-margin
  members of P3-all\*.
* The general case also needs L1 or L3.
* L8 is an independent alternative.
* **L4 (inverse step, I1).** Suppose parity-respecting factored residue data on a profile with
  per-column weight `tau` have per-unit demand `log m_(v,s) <= beta ||v||_1 log M` for every sparse `v`
  outside `L` with `||v||_1` up to the completion range `~ 2.08 X/tau`. Then the labels are close to
  certified: induced by Gaussian integers `zeta_j` of controlled height.
  * *Status:* open. It is of discrepancy type, a statement about concentration over
    `e^(Theta(k log M))` monomials at once.
  * It cannot come from single-modulus counting (Lemma 5.1).
  * A counterexample would be a "multiplicative code with prescribed Legendre parities", and so a
    first step towards capping P3-all\*.
* **L5 (twin channel, M1').** On twin profiles, the twin labels `Lambda_y(q) = kappa_y(q) - F_y(q)`
  cannot both:
  * avoid, at all but `O(||w||_1 + 1)` moduli, the values `-+<w, (ell_J, U)(q)> + (m_q/4) Z` for
    every sparse free vector `w`; and
  * make `kappa = F + Lambda` a point code with maximal pair agreement `O(log M)`.

  *Status:* open. Heuristically the twin channel alone gives the random rate.
* **L6 (I2, certified implies genuine size).**
  * Proved for exact Gaussian coherence under GRH (Thm 3.2), and for additive certification
    unconditionally (Prop. 3.1).
  * **Missing:** a branch for approximately certified designs. Twisted coherence is certified with
    `Norm zeta ~ X << p_j`, so (I2) as stated is false there.
  * Also missing: an unconditional version for exact coherence. Unconditional Chebotarev bounds are
    useless here because `p_j > X`.
* **L7 (I3, genuine size implies Poisson tails).** Collision counts of dense coherent integers (the
  pair integers `B_xy`, the twin values) in `(M^(1-eta), M]` have Poisson-type tails. This is
  Hypothesis `H` of `p3lower.md`, or the smooth-part statistics of `p3G.md`.
  * *Status:* open even for genuine primes.
* **L8 (slab covering).** Any necessary condition on the data of P3-all\* beyond the per-unit heights:
  a statement that the slabs
  `{phi : |sin(<v, phi> - s pi/4)| < nu_s m_(v,s) e^(-w(v)/2)}`, `v` outside `L`, cover the free
  torus for data of a stated kind.
  * *README remark:* such a statement cannot be a threshold on the **sum** of the slab measures. On
    `p3G.md`'s verified prime sets the genuine sums are 4052 and 15190, yet the genuine angles
    `theta` lie outside every slab (planted point, `p3G.md` Lemma 2.1, valid for every set of genuine
    primes). So a covering lemma must use how the slabs overlap, e.g. the tight-relation circuits of
    `p3G.md` §2, or the measure of their union, not the sum of their measures.
  * *Status:* open. The union bound is sufficient for non-covering estimates in constructions, not
    necessary for membership.

**Route 4: genuine residues (P3-G). This route decides row 6 only, not growth.** Where the closure
is complete, P3-G **is** the actual problem (Thm 3.1). So these lemmas would show that P3-G is no
easier than actual circles, or else that it is capped. Neither gives growth.
* **L9 (selected prime sets).** For every set of `r ~ M` primes of size `M^(O(1))`, chosen
  adversarially, the full closure `Pi*(pi, X)` at `X = M^(O(1))` leaves only `o(M)` columns unpinned.
  * *Status:* open. Selected sets evade every closure computed so far (`referee_p3G.md` §4). The
    weight-`>= 4` closure has not been computed.
* **L9' (near-relation count).** An unconditional form of `p3G.md` Remark 3.5: a set of `r ~ M`
  genuine primes of size `M^(O(1))` carries `O(M/log M)` independent genuine near-relations.
  * L9 + L9' and the dichotomy (Thm 3.2, Cor. 3.4) would exclude squarefree non-actual members on
    such prime sets.
  * *Status:* heuristic only.
* **L10 (smooth values at primes).** A lower bound, of order `tau_X(M^K)`, for the density of
  `X`-admissible values `D_s(pi_i pi_j^(+-))` at random Gaussian primes. These are smooth values of
  bilinear forms at primes. L10 would make `p3G.md` Cor. 5.4 unconditional for **random** prime sets.
  *Status:* open.
* **L11 (the window).** Decide whether P3-G has non-actual members with
  `H_X <~ log p <~ 2 H_X`. Fakes are proved only from `log p >= 2.08 X` on (Prop. 6.1).
  *Status:* open.

**Construction side: capping P3-all\*.**
* C1 would make P3-all\* capped (row 5). Rows 6-9 of §2.1 would still be open: genuine residues,
  super-polynomial moduli and additive data.
* C2 is the first finite test.
* C3 is a sub-question for coherent designs.
* **C1 (low-demand factored data with completion).** On the `sharp.md` §4 Paley profile with
  `log P = O(log M)`, construct factored residue data at the prime powers `<= 3M` such that:
  * every pair demand is `log m_xy = O(log M)`;
  * `m_v <= M^(beta ||v||_1)` for every `v` outside `L`, with `beta` fixed.

  `sharp.md` Prop. 6.6 then gives `W = O(M log M)`.
  * Theorems 1-2 of `p3Dconstruct.md` and §3-§4 of `p3Dlower.md` exclude random, additive,
    exact-coherent (at the certificate level, GRH) and G-coherent-with-small-heights data.
  * The task's design D6 is the open candidate.
* **C2 (skew double core).** For each `M`, minimise `D_max` over P3-char, respectively P3-all, residue
  data on the skew double core. Its margins are uniform, so `W >= (2(M-1)/K)(D_max - c_C)` with
  `K = 1 + o(1)`.
  * If `D_max >> log M` is forced, the double core is the model case for (T).
  * Otherwise, the double-core analogue of `sharp.md` Prop. 6.6 has to be written and checked.

  This is a finite exact computation per `M` (evidence). It is the round-8 "decisive first test",
  still not done.
* **C3 (parity repair on `q = +-1 mod 16`).** Is there a deterministic 2-Sylow repair whose certified
  cost is linear in `||v||_1` and independent of `k(q)`? A negative answer turns Heuristic H1 of
  `p3Dconstruct.md` into a theorem for size-certified coherent designs.

**(L_7)** (uniform theorem, not growth; `L7.md` §9, `referee_L7.md` §9).
* **Q1.** Is `lambda_V(beta) < 0` for every covering family with one member through a general point?
  Known: forgetful `-9`, Kapranov `-26`, `V_A <= -23`. Lemma 5.3 makes the test a finite exact
  computation for each explicit family.
* **Q1'.** Can good curves through `Sing(F_A)` (e.g. through `x_m`) make up the `23L` deficit?
  Heuristically no, but not proved.
* **Q2.** Do flexible rich families exist (a positive-dimensional `Q(x)`-rational set of members
  through `x`)?
* **Q3.** Do class-`beta` rational curves exist, and is their 8-dimensional family unirational over
  `Q`? Numerical continuation reaches 32 of 35 triple events.
* Also open: Prop. 5.2 (product charts), which needs leaves avoiding `Sing(F_{V_k})`, and the
  rational-connectedness route of `L7.md` §9.1, which needs a smooth proper model, the
  field-of-moduli twist and openness of goodness.

---------------------------------------------------------------------------------------------

## 4. Corrections to earlier summaries

The repository is read-only for this round. These corrections are recorded here, to be applied when
round 7 is merged.

1. **`round8/README.md` §3, "Genuine residues are not a shortcut".** The sentence "the certified
   closure pins every fake angle to the genuine one, so P3-G there coincides with the actual problem"
   holds only for members whose prime set is **exactly** one of the six verified random sets
   (exact), and in the non-adaptive model `T`. For selected prime sets it is open (§1.3).
   * The conclusion drawn there, "so (T) must be proved for P3-all", should be restated.
     `P3-G ⊂ P3-all`, so a lower bound for P3-G is a *weaker* target than (T), and it still gives
     growth for actual circles.
   * Whether it is a genuinely simpler problem than actual circles is exactly the open collapse
     question (§3, L9-L11). Where the collapse holds, it is not simpler.
2. **`round8/README.md` §2.2, "Certificate-level cap (draft)".** This is now refereed (§1.1, M1).
   * The PSD part is within `O(M)` of the Perron bound, not equal to it.
   * In P3 classes `D` depends on the profile, so "adapted profiles" is a statement about the
     relaxation.
3. **`round8/README.md` §2.3(a), last row (`p3lower.md` Thm C).** This is now refereed. The hypothesis
   can be weakened to `H_1` (one pair, one `k`). The theorem is Prop. B.1 applied to that assumed
   pair; the growth sits in the hypothesis.
4. **`round8/README.md` §3, "decisive first test".** Not carried out in round 7 (§3, C2).
5. **Top-level `README.md` §3.** "Round 7 ... (in progress)" should read "refereed". Its growth-status
   paragraph is unchanged: P3-all is still the open local class. Add genuine residues (P3-G,
   P3-G-char) and super-polynomial moduli to the list of open classes (§2.1, rows 6-8).
6. **Repository `round7/README.md`** (one-paragraph stub). Superseded by this file. Its sentence "no
   pair-cut certificate can prove growth in P3-all" is correct (Cor. A.5, refereed).
7. **`round5/README.md`, `exact.md`.** "Referee running in round 7": no referee report for `exact.md`
   exists among the round-7 files, so `exact.md` is still a draft. Results that use it are marked as
   such in `round8/README.md`.

---------------------------------------------------------------------------------------------

## 5. Files and caveats

| write-up | referee report | verdict | scripts |
| --- | --- | --- | --- |
| `p3lower.md` (also in the repository) | `referee_p3lower.md` | survives; 11 minor issues (PSD gloss, parity gloss, factor-2 slip, scope of Thm C and of "adapted profiles", citations) | `p3lower_checks/`, `referee_p3lower_checks/` |
| `p3data.md` | `referee_p3data.md` | survives as a data paper. Major issues: windings; lower bound read as cost. Medium: Example 6.2A is `W_sc` only | `p3data_checks/` (incl. `r2/`), `referee_p3data_checks/` |
| `p3G.md` | `referee_p3G.md` | mathematics correct; **headline not established** (selected primes; subsets; model) | `p3G_checks/`, `referee_p3G_checks/` |
| `p3Dlower.md` | `referee_p3Dlower.md` | survives as a modest partial barrier; evidence table withdrawn; four headline overstatements; class definition fixed | `p3Dlower_checks/`, `referee_p3Dlower_checks/` |
| `p3Dconstruct.md` | `referee_p3Dconstruct.md` | survives as a barrier for construction methods; the answer to the task's suggestion is not proved; (Ind) repaired | `p3Dconstruct_checks/`, `referee_p3Dconstruct_checks/` |
| `L7.md` (r7a: `L7_prev_attempt1.md`) | `referee_L7.md` | survives as a partial barrier with smaller scope; Prop. 5.2 unproved; `Sing(F_A)` loophole | `L7_checks/r8/` (r7a: `L7_checks/r7/`), `referee_L7_checks/` |

Every script runs from its own directory with `python3 <script>`, and its output is stored next to it
as `*_out.txt` or `*.log`. `patch_qrun.py` and `patch_qrun2.py` in this folder are working patches for
the `p3data.md` runs. Working notes for this README are in `readme_work/notes.md`.

**Caveats.**
* *Unrefereed material in this file.* The following have been checked by nobody but their author:
  * the referees' own sketches: the two-column repair of Theorem 1; the realisable variant of
    Example 6.2A; the winding-robust form of Thm 5.1;
  * referee Lemma R.1, proved in `referee_p3G.md` but not re-refereed;
  * every *README remark*: the parity-free nature of `U_X`; the impossibility of a pure-measure
    covering lemma (L8); the formulation of the missing lemmas in §3.
* *Evidence is not used in any proof.* This covers the tight densities, the thresholds, the solver
  upper bounds, the `lambda_max` and `D_max` statistics, the selected-prime closures, and the (L_7)
  continuation runs.
* *Literature.* The following were cited by the authors or referees and not re-checked here:
  Goemans-Williamson / Nesterov, Deza-Laurent, Donoho-Stark / Elad-Bruckstein, Hildebrand,
  Hildebrand-Tenenbaum, Huxley, Janson-Luczak-Turova-Vallier, Rosser-Schoenfeld, GRH/Chebotarev
  (pseudosquares), Gale / Gelfand-MacPherson, Keel-McKernan and Colliot-Thélène. The novelty of
  Theorem A of `L7.md` was not checked.
* *Scope of every barrier.* Each barrier is relative to the class or method it names (§2). None of
  them covers P3-all\*, P3-G, super-polynomial moduli or additive data. Nothing in round 7 changes the
  bound for actual circles.

---------------------------------------------------------------------------------------------

**Repository note.** The raw trial log `p3data_checks/genuine_trials.jsonl` (156 MB) is not committed;
`p3data_checks/genuine_stats.py` regenerates it. Every table built from it is reproduced in
`p3data.md` and `referee_p3data.md`.
