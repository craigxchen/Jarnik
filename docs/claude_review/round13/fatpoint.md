# Round 13 — fatpoint: base loci of the local system under Loc(b0), and the barrier test

Status: COMPLETE (unrefereed). Write-up of the round-13 task "fatpoint" (route 5, Loc(b0), base loci, barrier test).

## 0. Summary and verdict

**No growth improvement and no uniform bound is proved.** `M <= (2/3 + eps) log R/loglog R` stands. Everything below is
at `T = 1`, `w = sqrt(8 d ln 4r)`, `j = floor(w/2)`, `m' = ceil(4w)`, regime `d >= 316 ln 4r`, unless a model is named.

**Question.** Under locality `Loc(b_0)` the tangent cone `TC(A,m')` is the real base locus of an explicit local linear
system `V_{𝓕_{b_0}}(m') = ∩_F (I_F)^perp_{m'}`. Can that base locus contain nondegenerate, collision-poor points (task 1)?
Do fake pattern boxes have such points (task 2)? What conditional theorem results (task 3)?

| # | result | status |
| --- | --- | --- |
| 1 | **Structure of the local base locus.** `Bs_R(V_𝓕(m)) ⊂ TC(A,m)` for every family set, without `Loc` (so every nondegenerate local base point is a genuine tangent-cone point and feeds fake Thm 2.1); `=` under `Loc_m`. One family: its base locus is `V_F TC(B_F, m)`, a piece of a rational subspace. Joins of such pieces with total tangent degree `< m` lie in the base locus. **Span principle:** a nondegenerate base point needs a `Q`-irreducible base component whose linear span is everything; an algebraic one has degree `>= d` over `Q`. (Lemma 2.1, Props 2.2-2.4) | proved |
| 2 | **Product certificates** (generalised zonotopes): `prod <z_i,u> in J(A)_m^perp` iff the polarisations of `J(A)_m` vanish at `(z_1..z_m)`; then `TC` has no nondegenerate point and `kappa >= |prod <z_i,u>|/(|c|_1 d^m)`. Local criterion: width counts `#{i : <v,z_i> != 0} <= tau_A(v)` (exact for `b = 1`), and `#{i touching F} < alpha_F` sufficient, `<= rho(B_F)` necessary. Under `Loc(b_0)`, local certificates are certificates. (Props 2.5-2.6, Cor 2.7) | proved |
| 3 | **No collision certificates:** for every sign profile with the pair inequalities and `W > 8 ln C'`, no product of `m >= 4w` collision forms satisfies even the `b_0 = 1` counts; certificates must use heavy characters (`|z(P)| >= 2` on average). (Thm 3.1, Prop 3.2) | proved |
| 4 | **Uncertainty budget:** for profiles with a Hadamard core, a certificate in the width system needs `m^2 <= (2d + min(tau_0, m)) T_core/n`. Explains all small exact models (`sdc7`, `pal7`, `M8`, `M8'`); not binding in the regime (`n = 2^{12}`: within a factor 1.15). (Thm 3.3) | proved; exact numbers |
| 5 | **`b_0 >= 2`:** every parity (binary) family, hence every Sylvester binary column family, imposes only **counts**: `<= j|s| - 1` generators per circuit `s` suffices, `<= j|s|` is necessary. Complete XOR blocks reduce self-similarly to the class-sum box. (Prop 4.1, Cor 4.2, Lemma 4.3) | proved |
| 6 | **Complete Fano block `Y_3(2)`, degree `7 = 3.5j`:** of the 91,272 products of 7 heavy restricted generators that satisfy every grid and every line condition, **none** with an axis touched `> j` times lies in `J(Y_3(2))^perp` (0/88,200), while those with zonotope cost `<= 2j` do (1536/3072); mixed light/heavy samples agree (0/4705 above zonotope cost). In degree `<= 3j` costs are counts. So heavy generators pay the **zonotope cost** once the non-parity generators enter. (§5.3) | exact (negatives rigorous), plus sampled rows |
| 7 | **Dichotomy (heuristic):** if complete blocks of growing dimension impose the zonotope cost, the averaging of Thm 3.1 returns and coset-difference certificates are impossible for gap-admissible data; with `Loc(b_0)` for fixed `b_0` only very heavy generators can escape. The certificate route needs `Loc(b_0)` **and** a failure of the zonotope-cost pattern in large blocks, or generators heavier than every active block. (§5.3) | heuristic, one exact data point |
| 8 | **Conditional theorem A:** `Loc(b_0)` + a local certificate + quantitative locality (height of `c`) + tight characters `=> W/4 <= E/2 + psi + m' ln(Cd) + ...`, contradiction, `n = O(loglog R)` through a character form `MIN^char` of the final hypothesis (weaker than `MIN^pair` as a hypothesis, but with non-pair arithmetic inputs). (Thm 6.1, Cor 6.2) | proved implication; hypotheses open |
| 9 | **Conditional theorem B (no locality needed), Sylvester-type profiles:** the box contains direct sums of `~ m'/m_0` coset blocks `Syl^box(N_0, j_0)` sharing the slab budget; with random coset decompositions, `MIN^pair` holds at every large `n` as soon as **one fixed block** satisfies `ZE-TC^braid` at an LP-supercritical order `m_0 > 8 j_0` (constants by Lojasiewicz, ineffective but fixed). Then Sylvester-profile clusters have `n < n_0`. Needs `N_0 >= 32`: beyond exact computation. Does **not** apply to fake (Paley-type) or generic profiles, whose small blocks have full traces. (Lemma 6.3, Lemma 6.5, Thm 6.4, Cor 6.6) | proved implication; `BLOCK` open |
| 10 | **Barrier test (task 2).** No exact fake model reaches the regime ratio (`rho/(2L) <= 2` in all computable models). At their top orders the `n = 8` fake models are **local**: `Loc_8(1)` holds exactly for `sdc7` (`L = 2`), `Loc_7(3)` for `pal7`, `Loc_{11}(1)` for `M8'`. For `sdc7` the width system's real base locus therefore equals `TC`, and it contains the certified nondegenerate, well-separated point of fake §3.2: a pure fat-point system with collision-heavy integer fat points **can** have nondegenerate real base points (here by a dimension count, `H_A(8) = 1`). Product certificates exist in the exact models exactly where the box contains a zonotope of length `m`, and all of them are zonotopes; `M8` shows that certificates are not necessary for the absence of nondegenerate points. The Diophantine condition (2.2) is vacuous at fixed `n`. (§5) | exact / evidence as labelled |

**Answer to task (1).** Not decided. Proved: every nondegenerate point of the local base locus is a global phenomenon
(span principle), every certificate excluding such points must be built from heavy characters, never from collision
forms (Thm 3.1), and for `b_0 = 2, 3` (in fact for all binary families) the local conditions on such certificates are
plain counts. Whether `Bs_R(L_{m'})` itself has nondegenerate points (`b_0 = 1`) remains open; under the counterfactual
`Loc(1)` the question is a combinatorial design problem whose averaged constraints are satisfiable and whose
uncertainty budget is tight only near `n = 2^{12}`.

**Answer to task (2).** Not decidable at regime ratios with exact computation. At the computable ratio 2 the answer is
**yes** for the skew double core `sdc7`: its tangent cone at the top order is exactly the base locus of the width
system, and that base locus has certified nondegenerate, well-separated real points (so fake Thm 2.1 applies). This is
an underdetermined system (`H_A = 1`); the overdetermined exact models (`M8`, `M8'`, `pal7` at `m = 7`) have no
nondegenerate tangent-cone points. The decisive object for the regime is the complete-block image `Y_t(j)` (an
`n`-independent lattice set), and its first exact case shows the zonotope cost, which is the behaviour that would make
the certificate mechanism fail on gap-admissible data, fakes and actual clusters alike.

**Answer to task (3).** Two conditional theorems (8, 9). Theorem A needs `Loc(b_0)` plus three further hypotheses that
locality does not supply (local certificate, height, tight characters); item 7 makes their joint truth doubtful.
Theorem B bypasses locality for Sylvester-type profiles and reduces them to one fixed finite qualitative statement,
but it cannot reach generic or fake profiles.

**Fake-circle test.** Rigidity (through Step 0) remains the only input fakes violate. Every bound proved here is about
the box and holds verbatim for fake pattern boxes; in particular Theorem 3.1 and Proposition 3.2 apply to fakes, and
Theorem 6.4 does not apply to the fake profiles (full traces). Nothing here distinguishes an actual cluster from a fake
on the same box.

## 1. Setting, notation, the hypothesis Loc(b0)

**Labels.** *Theorem / Proposition / Lemma*: proved here with all quantifiers, or quoted with an exact source.
*Exact*: integer or exact-rational computation (python-flint), or a mod-`p` computation used only in the direction in
which it is rigorous. *Evidence*: everything else. Scripts are in `fatpoint_checks/` (`x.py` writes `x_out.txt`),
working notes in `fatpoint_notes.md`. Sources: round 12 (`README.md`, `fake.md`, `zetc.md`, `t1rigour.md`,
`nearstrata.md`) and the round-13 siblings `loctheory.md`, `loccompute.md`, `refJ.md` (all unrefereed).

**1.1 The box.** `T = 1`, `n` points, base point `0`, `d = n - 1`, `k' = (k_1..k_d)`, `k_0 = -sum k'`,

```text
A = { k in Z^n_0 : k' in {-1,0,1}^d, |<k,beta_p>| <= w (all columns p) },   w = sqrt(8 d ln 4r),  j = floor(w/2),
m' = ceil(4w)   (so 8j <= m' < 8j + 9),   regime: d >= 316 ln 4r (then m' <= m_max).
```

(t1rigour §1: the base coordinate `k_0 = -sum k'` is constrained only through the columns.) For a sign column
`beta_p = 2*1_{S_p} - 1` one has `<k,beta_p> = 2k(S_p)` on `Z^n_0`, so the slab is `|k(S_p)| <= j`.
In the small exact models (§5) the half-width is an ad hoc `L` instead of `j`, and the order `m` is stated per model.

**1.2 Ideals, inverse systems, tangent cone.** `S = Q[x_1..x_d]` acts on forms in `u = (u_1..u_d)` by
`f -> f(partial_u)`; for `v in Q^d` write `l_v(x) = <v,x>` and `D_v = l_v(partial_u)`. For `f in S_m` and a form
`Phi` of degree `m`, `<f,Phi> := f(partial)Phi` is a perfect pairing, and `<l_u^m, Phi> = m! Phi(u)`. `I(A)` is the
vanishing ideal, `J(A)` its ideal of leading forms, `H_A`, `rho(A)` as usual. By compute R2 / refJ Lemma 2.0,

```text
span{ Phi_c : c in V(A,m) } = J(A)_m^perp ,     TC(A,m) := { u : Phi_c(u) = 0 for all c in V(A,m) } = { u : l_u^m in J(A)_m } ,
nu_A(u) = max{ e : l_u^e notin J(A) } ,   TC(A,m) = { nu_A < m }   (zetc Lemma 1.1).
```

`TC_C` denotes the complex points. `kappa`, `V_H`, `eps_0`, `W`, `abar`, `eps_xy`, `Gamma`, `P(L)` are as in round-12
README §1. A real point `u` is **nondegenerate** if `<z,u> != 0` for every `0 != z in Z^d` (it lies on no rational
hyperplane), and **collision-poor** if `Gamma(u) = O(log n)` (equivalently, polynomial pair gaps on average).

**1.3 Families.** A *family* is a tuple `F = (v_1..v_b)` of linearly independent integer vectors (directions in
`k'`-space); `pi_F(k') = (<v_q,k'>)_q`, `B_F = pi_F(A) ⊂ Z^b`, `y_q = l_{v_q}`. The **local ideal** of `F` is the ideal
`I_F ⊂ S` generated by `pi_F^* J(B_F)`; `I_F ⊂ J(A)` because pulled-back vanishing polynomials vanish on `A`. By span
invariance (loccompute §1) `I_F` depends only on `span_Q F`. The grid hull `Gamma_F = prod_q I_q` (`I_q` the values
of `<v_q,.>` on `A`), `s_F = sum (|I_q| - 1)`, and `alpha_F = s_F - rho(Gamma_F \ B_F)` are as in loctheory §1.3, §2: by
loctheory Thm 2.2, `J(B_F)_e = J(Gamma_F)_e` exactly for `e < alpha_F`. The **local system** of `F` in degree `m` is

```text
V_F(m) := (I_F)_m^perp = { Phi in Sym^m : g(D_{v_1}, ..., D_{v_b}) Phi = 0 for every g in J(B_F) } ,
```

and for a set `𝓕` of families, `J_𝓕 := sum_{F in 𝓕} I_F` and `V_𝓕(m) := (J_𝓕)_m^perp = ∩_{F in 𝓕} V_F(m)`.
`𝓕_{b}` is the set of all families with `b' <= b` directions. For `b = 1` the local ideals are the width powers
`l_v^{tau_A(v)+1}` (`tau_A(v)` = number of values of `<v,.>` on `A` minus one; refJ Prop 2.1), and `V_{𝓕_1}(m)` is the
fat-point system `L_m` of refJ §2.2 (multiplicity `>= m - tau_A(v)` at `[v]`).

**1.4 The hypothesis.** README §7.1 restricts to families with `rho(pi_F A) < ceil(4w)`; loctheory Prop 2.4 and
loccompute §3.2 show that this restriction loses generators (the correct cost is `alpha_F`). I therefore state locality
**without** a cost restriction; by loctheory Cor 2.3 a family with `alpha_F > i` contributes in degree `i` only width
powers of its own directions, which are `b = 1` families, so nothing is lost or added by restricting to
`alpha_F <= i`.

> **`Loc_D(b_0; A)`.** For every `i <= D`: `J(A)_i = (J_{𝓕_{b_0}})_i = sum_{F in 𝓕_{b_0}} (I_F)_i`.
>
> **`Loc(b_0)`** (regime form). There is `n_L` such that `Loc_{m'}(b_0; A)` holds for every admissible profile with
> `n >= max(n_L, n_R(r))`, at `T = 1`, `m' = ceil(4w)`.

Only the degree `m'` part enters the tangent cone: `TC(A,m')` depends on `J(A)_{m'}` alone (§2.1). Lower degrees
matter for transversal orders and for the quantitative step (§6). `b_0` is independent of `n`. Note that loctheory
Cor 5.3 forces `b_0 >= 8` for Sylvester profiles (circuits of size `<= 8` are essential among column families), and
loccompute Prop 4.4 shows `Loc_7(5)` fails over `Q` in the exact model Paley-12 (`L = 2`). So `Loc(b_0)` is a
hypothesis, possibly false for every `b_0` (complete XOR blocks, loctheory §5.3); everything below that uses it says
so.

## 2. The local linear system and its base locus

Throughout §2: `A ⊂ Z^d` is any finite set with `0 in A` (in particular any `T = 1` box), `m >= 1`, and all statements
hold over `C` and, by taking real points, over `R`. No hypothesis on `n` or the profile.

### 2.1 Base locus = power membership; the unconditional inclusion

**Lemma 2.1.** For every homogeneous ideal `K ⊂ S`,
`Bs_C(K_m^perp) = { u in C^d : l_u^m in K_m }`. If `K ⊂ J(A)` then `Bs_C(K_m^perp) ⊂ TC_C(A,m)`; if moreover
`K_m = J(A)_m` the two sets are equal. In particular, for every set `𝓕` of families,

```text
Bs_R( V_𝓕(m) )  ⊂  TC(A,m)        (unconditional),          Bs_R( V_{𝓕_{b_0}}(m') ) = TC(A,m')   under Loc_{m'}(b_0; A).
```

*Proof.* `Phi(u) = <l_u^m, Phi>/m!`. So `u` is a common zero of `K_m^perp` iff `l_u^m` pairs to zero with `K_m^perp`, iff
`l_u^m in (K_m^perp)^perp = K_m` (finite dimension, perfect pairing). If `K ⊂ J(A)` then `K_m^perp ⊃ J(A)_m^perp`, a
larger system, with a smaller base locus, and `TC_C(A,m) = Bs_C(J(A)_m^perp)` by §1.2. QED

**Consequence for the barrier test.** A nondegenerate real point of `Bs(V_𝓕(m))`, for *any* `𝓕`, is a nondegenerate
real point of `TC(A,m)` without any locality hypothesis, and fake Thm 2.1 applies to it directly. Locality is needed
only in the other direction (to know that `TC` is not larger than the local base locus).

### 2.2 One family: the base locus is a piece of its span

**Proposition 2.2.** Let `F = (v_1..v_b)` be a family, `B = B_F`, `V = (v_1 | ... | v_b)` (a `d x b` matrix).
* (a) `Phi in V_F(m)` iff for every `u_0 in C^d` every homogeneous component (in `t`) of `t -> Phi(u_0 + V t)` lies in
  the inverse system `J(B)^perp` of `B` in `C[t_1..t_b]`.
* (b) `Bs_C(V_F(m)) = { V a : a in TC_C(B, m) }`, where `TC_C(B,m) = { a in C^b : (a.y)^m in J(B)_m }`.
* (c) Hence `Bs_R(V_F(m)) ⊂ span_R(v_1..v_b)`, a rational subspace; for `b < d` it consists of degenerate points.

*Proof.* (a) `(g(D_{v_1},...,D_{v_b}) Phi)(u_0 + V t) = g(partial_t)[Phi(u_0 + Vt)]` by the chain rule. `J(B)` is
homogeneous and `g(partial_t)` maps `t`-homogeneous components to `t`-homogeneous components, so `g(partial_t)P = 0`
for all `g in J(B)` iff each component of `P` is annihilated. Finally `Phi in (I_F)_m^perp` iff `g(D_V)Phi = 0` for every
`g in J(B)`: if `g(D_V)Phi` (a form of degree `m - deg g`) is nonzero, it pairs nontrivially with some `h`, and then
`(h g)(partial)Phi != 0` with `h g in (I_F)_m`; the converse is the same identity read backwards.

(b) Extend `y_q = l_{v_q}` (`q <= b`) by linear forms `y_{b+1}..y_d` with rational coefficients to a basis of `S_1`.
Then `S = Q[y_F] ⊗ Q[y_rest]` and `I_F = J(B) ⊗ Q[y_rest]` is bigraded. Write `l_u = a.y_F + beta.y_rest`
(unique `a, beta`). Then `l_u^m = sum_e binom(m,e) (a.y_F)^e (beta.y_rest)^(m-e)` is the bigraded decomposition, so
`l_u^m in I_F` iff every term lies in `J(B)_e ⊗ Q[y_rest]_{m-e}`. If `beta != 0`, the term `e = 0` is a nonzero element of
`Q[y_rest]_m`, while `J(B)_0 = 0` (`B` is nonempty); so `l_u^m notin I_F`. If `beta = 0`, the only term is `e = m`, and it lies
in `I_F` iff `(a.y)^m in J(B)_m`. Finally `beta = 0` iff `l_u in span(y_q : q <= b)` iff `u = V a`. Lemma 2.1 converts this
into the statement about `Bs`. (c) is immediate. QED

### 2.3 Joins and the span principle

For `a in C^b` put `nu_B(a) = max{ e : (a.y)^e notin J(B) }` (`nu_B(0) = 0`). The set `{e : (a.y)^e in J(B)}` is upward
closed, since `J(B)` is an ideal.

**Proposition 2.3 (joins).** Let `F_1..F_s` be families with matrices `V_i`, `a_i in C^{b_i}`, and `u = sum_i V_i a_i`.
If `sum_i nu_{B_{F_i}}(a_i) < m`, then `u in Bs_C(∩_i V_{F_i}(m)) ⊂ TC_C(A,m)`.

*Proof.* `l_u = sum_i l_{V_i a_i}`, and `l_{V_i a_i} = pi_{F_i}^*(a_i.y)`. In the multinomial expansion of `l_u^m`, every term
`prod_i l_{V_i a_i}^{e_i}` has some `e_i > nu_i` (else `m = sum e_i <= sum nu_i`). For that `i`,
`(a_i.y)^{e_i} in J(B_{F_i})`, so the term lies in `I_{F_i}`. Hence `l_u^m in sum_i I_{F_i}`; apply Lemma 2.1. QED

This contains the cheap-span statement (refJ Prop 2.2: `b_i = 1`, `nu = tau`), the projection-cheap spans
(`s = 1`), and their joins. All points it produces lie in the rational subspace `sum_i span F_i`.

**Proposition 2.4 (span principle).** Let `𝓕` be any set of families and `W = sum_{F in 𝓕} span_Q F`. Then
`Bs_C(V_𝓕(m)) ⊂ W ⊗ C`. Consequently:
* (a) a nondegenerate real point of `Bs(V_𝓕(m))` exists only if `W = Q^d`;
* (b) if `u in Bs_R(V_𝓕(m))` is nondegenerate, the `Q`-Zariski closure `Z_u` of `u` is a `Q`-irreducible subvariety of
  `Bs_C(V_𝓕(m))` whose complex linear span is `C^d` (its geometric components are Galois conjugate and together span
  `C^d`); in particular, if `u` is algebraic, then `[Q(u) : Q] >= d`.

*Proof.* Let `z in Q^d` with `<z, v> = 0` for every `v in W`, and `ell(u) = <z,u>`. For `v in W`, `D_v ell = <z,v> = 0`,
so `g(D_F) ell^m = 0` for every homogeneous `g` of positive degree and every `F in 𝓕`; `J(B_F)` has no constants, so
`ell^m in V_F(m)` for every `F`, hence `ell^m in V_𝓕(m)`, and `Bs ⊂ {ell = 0}`. Intersecting over all such `z` gives
`Bs ⊂ W ⊗ C`. (a) follows. (b) `V_𝓕(m)` is spanned by rational forms (the `J(B_F)` are rational), so `Bs_C` is defined over
`Q` and contains `Z_u`. The linear span of `Z_u` is stable under `Gal(Qbar/Q)`, hence defined over `Q`; it contains `u`, which
lies on no rational hyperplane, so it is `C^d`. If `u` is algebraic, `Z_u` is the set of Galois conjugates of `u`; they
span `C^d`, so there are at least `d` of them. QED

So every classical mechanism that builds base components from a *subset* of the fat points or families (linear cycles,
secant varieties and rational normal curves through those points, Cremona-special loci; Brambilla-Dumitrescu-
Postinghel) produces nondegenerate points only when the subset spans: a nondegenerate base point is necessarily a
**global** phenomenon of the system, involving families whose spans add up to all of `Q^d`.

### 2.4 Product certificates (generalised zonotopes)

The only device known to bound a tangent cone *from above* with explicit constants is the zonotope (fake Lemma 5.1(1),
nearstrata Thm 2.4), and zonotopes are excluded in the regime by the gap inequalities (fake Prop 5.3). The natural
generalisation is a product of rational linear forms lying in the tangent space.

**Proposition 2.5 (polarisation criterion).** Let `z_1..z_m in Z^d \ 0`, `P(u) = prod_i <z_i,u>`, and `K` homogeneous.
Then `P in K_m^perp` iff `Pol(f)(z_1, ..., z_m) = 0` for all `f in K_m`, where `Pol(f)` is the symmetric `m`-linear form
with `Pol(f)(x,...,x) = f(x)`. If `P in K_m^perp`, then `Bs_C(K_m^perp) ⊂ ∪_i {<z_i,u> = 0}`. In particular:
* if `P in J(A)_m^perp` (a *product certificate* for `(A,m)`), then `TC(A,m)` contains no nondegenerate point, and for
  every `c in V(A,m)` with `Phi_c = P`: `kappa_c(u) = |prod_i <z_i,u>|/(|c|_1 d^m)`;
* if `P in V_𝓕(m)` (a *local certificate*), then `Bs(V_𝓕(m))` contains no nondegenerate point; under
  `Loc_m(b_0)` with `𝓕 ⊃ 𝓕_{b_0}`, `P` is a product certificate.

*Proof.* Both sides are linear in `f`, and powers `l_w^m` span `S_m`; for `f = l_w^m`,
`<l_w^m, P> = D_w^m P = m! prod_i <z_i, w>` and `Pol(l_w^m)(z_1..z_m) = prod_i <w, z_i>`. So `<f,P> = m! Pol(f)(z)` for all
`f in S_m`. The rest is Lemma 2.1 and §1.2. QED

A zonotope `k_0 + sum_i [0,1] z_i ⊂ A` gives `P in J(A)_m^perp` (fake Lemma 5.1(1)); Proposition 2.5 allows any
`c in V(A,m)` with `Phi_c = P`. For `b = 1` the criterion is a pure count.

**Proposition 2.6 (local criterion).** Let `F = (v_1..v_b)` be a family, `z_1..z_m in Z^d \ 0`, `P = prod <z_i,.>`, and
`a_i = pi_F(z_i) = (<v_q, z_i>)_q in Z^b` (the restriction of `<z_i,.>` to `span F`). Put `N_F = #{i : a_i != 0}` and
`N_{F,q} = #{i : <v_q, z_i> != 0}`.
* (a) (*sufficient*) If `N_{F,q} <= |I_q| - 1` for every `q` and `N_F < alpha_F`, then `P in V_F(m)`.
* (b) (*necessary*) If `P in V_F(m)`, then `prod_{i : a_i != 0} (a_i.t) in J(B_F)^perp`; in particular
  `N_{F,q} <= |I_q| - 1` for every `q` and `N_F <= rho(B_F)`.
* (c) For `b = 1`: `P in V_{(v)}(m)` iff `#{i : <v, z_i> != 0} <= tau_A(v)`. Hence `P in L_m = V_{𝓕_1}(m)` iff the
  **count condition** `#{i : <v,z_i> != 0} <= tau_A(v)` holds for every `v in Z^d` with `tau_A(v) < m`; and the count
  condition is **necessary** for every product certificate (because `J_tau ⊂ J(A)`, refJ Prop 2.1).

*Proof.* On the flat `u_0 + V t`, `<z_i, u_0 + Vt> = c_i + a_i.t` with `c_i = <z_i,u_0>`. Let `N = {i : a_i != 0}`. The
`t`-homogeneous component of degree `e` of `P(u_0 + Vt)` is `prod_{i notin N} c_i` times a linear combination of the
products `prod_{i in I}(a_i.t)`, `I ⊂ N`, `|I| = e`. Each such product has `t_q`-degree `<= N_{F,q}`.
(a) `J(Gamma_F) = (y_q^{|I_q|})` (grid lemma, loccompute Lemma 4.1), so `J(Gamma_F)^perp_e` is the space of forms of
`t_q`-degree `<= |I_q| - 1` for all `q`; and `J(B_F)_e = J(Gamma_F)_e` for `e < alpha_F` (loctheory Thm 2.2(c)). All
components have `e <= N_F < alpha_F`, so they lie in `J(B_F)^perp`, and Proposition 2.2(a) applies.
(b) Choose `u_0` with `c_i != 0` for all `i` (possible: finitely many nonzero linear forms). The top component
(`e = |N|`) is `prod_{i notin N} c_i prod_{i in N}(a_i.t) != 0` and must lie in `J(B_F)^perp ⊂ J(Gamma_F)^perp`; and
`J(B_F)^perp_e = 0` for `e > rho(B_F)`. (c) For `b = 1`, `J(B) = (y^{tau+1})`, `alpha = infinity`, and (a), (b) coincide. QED

**Corollary 2.7 (what locality would buy).** Assume `Loc_{m'}(b_0; A)`. If `z_1..z_{m'}` satisfy the count condition and,
for every family `F` of `2 <= b <= b_0` directions with `alpha_F <= m'`, the bound `N_F < alpha_F` (or, more generally,
`P in V_F(m')`), then `P` is a product certificate and `TC(A,m')` has **no nondegenerate real point**. Without
`Loc`, the same conclusion holds for `Bs(V_{𝓕_{b_0}}(m'))`.

*Proof.* Proposition 2.6(a) for each `F`, `V_{𝓕_{b_0}} = ∩ V_F`, Lemma 2.1 and Proposition 2.5. QED

So under locality the qualitative barrier question becomes a **combinatorial design question** about the box: are there
`m'` integer vectors whose "touch counts" respect, in every direction, the width `tau_A`, and, in every small family,
the first non-grid degree `alpha_F`? §3-§4 study it in the regime; §5 tests it in exact models.

## 3. b0 = 1: the width (fat-point) system

Setting of §3: `T = 1`, any profile; for the quantitative statements, **sign columns** (`e_p = 1`, so
`d_xy = sum_p log p [x,y separated by S_p]`, `W = sum_p log p`), order `m = m' = ceil(4w)` unless stated.

### 3.1 What the base locus is known to contain

`Bs_R(L_m) ⊂ TC(A,m)` (refJ Prop 2.1; Lemma 2.1). It contains:
* every cheap span (`sum tau_A(v_i) < m`; refJ Prop 2.2), i.e. the BDP linear cycles;
* rational subspaces that are not cheap spans (the Fano spans of `M8`, refJ Prop 2.3);
* all joins of these (Proposition 2.3 with `b_i = 1`).

All of these lie in proper rational subspaces. By Proposition 2.4, the base locus of the subsystem imposed by any set
of fewer than `d` fat points lies in their span, so a nondegenerate point of `Bs(L_m)` needs a `Q`-irreducible base
component whose linear span is everything, forced by the system as a whole.

**Such components do occur (exact).** In the skew-double-core model `sdc7` (`T = 1`, `L = 2`, `m = 8 = rho`),
`locsys.py` certifies `J(A)_8 = (J_tau)_8` exactly (`dim = 356` in `R_8`, sandwich): **`Loc_8(1)` holds**, so
`TC(A,8) = Bs_R(L_8)` is the zero set of the single tangent octic. That octic has a certified regular collision-free real
zero `u*` with all gaps `>= 0.2496` (fake §3.2). Near `u*` the real zero set is a 6-dimensional analytic piece of an
octic that is irreducible over `Q` (R11 Prop 6.2), hence lies in no rational hyperplane; the points of that piece on
rational hyperplanes form a countable union of proper analytic subsets, so nondegenerate real zeros exist arbitrarily
close to `u*`. So the
real base locus of a pure width fat-point system **can** contain nondegenerate, well-separated points. The mechanism
is a dimension count: `H_A(8) = 1`, the system has one form, its base locus is a hypersurface. In the regime the system is
massively overdetermined (`dim V_H >= |A|/2`), and no dimension count produces base points.

Beyond dimension counts I could not decide the question by classical tools: Bezout-type forcing
(`sum_p mult_p(C)(m - tau_p) > m deg C`) is not excluded numerically, because two-class fat points (columns,
`m - tau = m - 2j ~ 3m/4`) are efficient; a curve forced in this way would have to pass through about `(4/3)(d-1)` column
points, i.e. the profile's columns would have to lie on a curve of degree `~d` in very special position. This is not a
theorem and I record it only as the reason why the question was attacked through certificates (§2.4) instead.

### 3.2 No certificate by collision forms

**Theorem 3.1 (braid products never satisfy the width counts).** Let the profile have sign columns only, let
`d_xy >= W/2 - 2 ln C'` for all pairs (true for every actual cluster, nearstrata Lemma 3.1, and for every gap-admissible
fake), and let `m >= 4w`. If `W > 8 ln C'` then there is no multiset of `m` pairs `{x_i, y_i}` (`x_i != y_i`; repetitions
allowed) such that `P = prod_i (u_{x_i} - u_{y_i})` satisfies the count condition of Proposition 2.6(c). In particular:
* no product of collision forms lies in `L_m`, in `V_𝓕(m)` for any `𝓕 ⊃ 𝓕_1`, or in `J(A)_m^perp`;
* `TC(A,m)` and `Bs_R(L_m)` cannot be certified to lie in the braid arrangement by a product certificate.

*Proof.* The direction of column `p` in `k'`-space is `1_{P_p}` (relative pattern), and `<z, 1_{P_p}> = z(P_p)`, which for
`z = e_x - e_y` is nonzero iff `S_p` separates `x` and `y`. On `A`, `k(P_p)` takes at most `2j + 1` values, so
`tau_A(1_{P_p}) <= 2j`, and the count condition gives `#{i : S_p separates x_i, y_i} <= 2j` for every `p`. Multiply by
`log p` and sum over `p`: `sum_i d_{x_i y_i} <= 2j W`. The left side is `>= m (W/2 - 2 ln C') >= 4w(W/2 - 2 ln C') >=
8j (W/2 - 2 ln C')`, using `w >= 2j`. Hence `4jW - 16 j ln C' <= 2jW`, i.e. `W <= 8 ln C'`, a contradiction. The second
and third bullet follow from Proposition 2.6(c) (necessity of the count condition) and `J_tau ⊂ J(A)`. QED

This is the count version of fake Prop 5.3. For pairs the count and the absolute sum coincide (`|z(P)| <= 1`), so the
zonotope obstruction transfers verbatim. For general generators it does not:

**Proposition 3.2 (certificates need heavy generators).** Same setting. Let `z_1..z_m` (`m >= 4w`) satisfy the count
condition at the column directions, and put `X_z = sum_p log p |z(P_p)|`, `Y_z = sum_p log p [z(P_p) != 0]`. Then
`sum_i Y_{z_i} <= 2jW`. If moreover each `z_i` satisfies the character gap inequality
`X_{z_i} >= W/2 - g_i` (fake Thm 5.2 / exact.md Prop 1.1: true for gap-admissible data with
`g_i = 2 log(C|z_i|_1/(2 sqrt 2))`, and for an actual cluster whenever `<z_i, delta> != 0`), then

```text
sum_i X_{z_i}  >=  2 sum_i Y_{z_i}  -  sum_i 2 g_i ,
```

i.e. on average (weighted by `log p`) a generator has `|z(P_p)| >= 2` on the columns it touches.

*Proof.* `sum_i Y_{z_i} = sum_p log p #{i : z_i(P_p) != 0} <= 2j W`, and `sum_i X_{z_i} >= m W/2 - sum g_i >= 4jW - sum g_i`. QED

**Reading.** Any product certificate in the regime is built from *heavy* characters (examples: coset differences
`1_{x+H} - 1_{y+H}` with `|H| >= 2` in Sylvester profiles, `|z(P)| = |H|` on touched columns). Its hyperplanes
`<z_i,u> = 0` are not collision hyperplanes, so:
* a certificate gives `TC ⊂ ∪_i {<z_i,u> = 0}` (no nondegenerate points), but **not** a collision structure;
* to turn it into a lower bound for `kappa(hat u)` one needs lower bounds for `|<z_i, hat u>|`, which the pair data do
  not give. For gap-admissible fakes the character inequalities give them; for actual clusters they hold unless
  `<z_i, delta> = 0` exactly, an arithmetic coincidence that must then be treated separately (§6).

### 3.3 A purely combinatorial budget: uncertainty

The averaging of Theorem 3.1 uses the arithmetic (pair distances). For Hadamard-type profiles there is a budget that
uses only the box.

**Theorem 3.3 (support-uncertainty budget).** Let the profile contain a *Hadamard core*: `n - 1` sign columns
`beta_1..beta_{n-1}` which together with the all-ones vector form an `n x n` Hadamard matrix (Sylvester, Paley, and the
fake profiles `sdc_q`, `pal_q` at `q = 7`; see `uncert_budget_out.txt` for which models have a full core). Let
`T_core = sum_{core} tau_A(beta)` and `tau_0 = tau_A(-1)` (the base coordinate). If `prod_{i=1}^m <z_i,u> in L_m` with
`z_i != 0`, then

```text
m^2  <=  (2d + min(tau_0, m)) T_core / n .
```

In regime boxes (`T_core <= 2j(n-1)`) this reads `m'^2 <= (2d + m') 2j (n-1)/n`.

*Proof.* Let `tilde z_i in Z^n_0` be the balanced extensions and `sigma_i = |supp tilde z_i|`. Count condition at
`v = e_x` (`tau_A(e_x) = 2`, `x != base`) and at the base direction: `sum_i sigma_i <= 2d + min(tau_0, m)`. Count condition
at the core columns: `sum_i #{core beta : <tilde z_i, beta> != 0} <= T_core`. With `Hd` the Hadamard matrix, `Hd tilde z_i`
has its all-ones entry `0` (balance) and its other entries `<tilde z_i, beta>`. The uncertainty principle for Hadamard
matrices, `|supp x| |supp Hd x| >= n` (from `|Hd x|_inf <= |x|_1 <= sqrt(|supp x|) |x|_2` and
`|Hd x|_2^2 = n |x|_2^2 <= |supp Hd x| |Hd x|_inf^2`), gives `#{core beta touched by z_i} >= n/sigma_i`. Hence
`sum_i 1/sigma_i <= T_core/n`, and `m^2 <= (sum sigma_i)(sum 1/sigma_i)` (Cauchy-Schwarz) gives the claim. QED

Values (`uncert_budget_out.txt`, exact): certificates in `L_m` are impossible for `m > 7` (`sdc7`), `> 8` (`pal7`,
`M8 = Syl(8,2)`), `> 10` (`M8' = Syl(8,3)`), `> 9` (`pal11`, Paley-12; `L = 2`). In the regime the budget is not binding:
for Sylvester `n = 2^12..2^20` the right side exceeds `m'^2` (by a factor `1.15` at `2^12`, growing like
`sqrt(n/ln n)`), because the base coordinate of the slab box is constrained only through the columns (t1rigour §1:
`k' in {-1,0,1}^d`, `k_n = -sum k'`). So Theorem 3.3 explains the small models (§5.2) but leaves the regime open.

### 3.4 The b0 = 1 design question in Sylvester profiles

For `b_0 = 1` (that is, under the counterfactual `Loc(1)`: loctheory Cor 5.3 shows that lines contribute at degree
`3j < m'`, so `Loc(1)` is false for Sylvester profiles) the question "does `Bs_R(L_{m'})` have nondegenerate points?" is
answered negatively by any family `z_1..z_{m'}` satisfying the count condition. Theorem 3.1 excludes pairs. The
averaged constraints for coset differences are satisfiable (column load `m'/2^{s+1}` against `2j`, subgroup-indicator
directions `1_V` with `tau_A(1_V) <= 2 floor(2j(1 - 2^{-t}))` for `codim V = t`, by `k(V) = 2^{-t} sum_{a in V^perp \ 0} hat k(a)`),
but the count condition quantifies over **all** integer `v` with `tau_A(v) < m'`, and a proof requires lower bounds on
`tau_A(v)` for every such `v`. I have no proof, and Theorem 3.3 shows that the budget is tight within a constant at
`n = 2^{12}`: there (`m' = 2256`, `j = 281`, `T_core <= 2j(n-1)`) the two budgets give
`9058 <= m'^2 n/T_core <= sum_i sigma_i <= 2d + m' = 10446`, i.e. the balanced generators have average support between
`4.01` and `4.63` (`uncert_budget.py`). §4 shows what the families of `b_0 >= 2` add, and §5.3 tests the decisive block condition.

## 4. b0 >= 2: circuits are count conditions; complete blocks are the crux

### 4.1 Parity images (any profile; in Sylvester profiles all binary column families)

**Proposition 4.1.** Let `F` be a family whose image is a parity section `B_F = P_C` (loctheory §3.3) with all
`|I_q| = 2j+1`, and let `z_1..z_m`, `a_i = pi_F(z_i)`. For a circuit `s` of `C` put `N_s = #{i : a_i|_s != 0}`.
* (a) If `N_{F,q} <= 2j` for all `q in F` and `N_s <= j|s| - 1` for every circuit `s` of `C`, then `P in V_F(m)`.
* (b) If `P in V_F(m)`, then `N_s <= j|s|` for every circuit `s`.

*Proof.* (a) By loctheory Thm 3.4(iii) (pulled back along `pi_F`, Cor 3.5), `I_F = (y_q^{2j+1}) + sum_s I_{F_s}`, so
`V_F(m) = ∩_q V_{(v_q)}(m) ∩ ∩_s V_{F_s}(m)`. The first factors are the counts `N_{F,q} <= 2j` (Prop 2.6(c)). For a circuit
`s`, `B_{F_s} = P_s` (Thm 3.4(ii)) and `alpha(P_s) = j|s|` (Thm 3.3), so Prop 2.6(a) applies to `F_s`. (b) Prop 2.6(b) with
`rho(P_s) = j|s|` (Thm 3.3), applied to `F_s`, using `V_F ⊂ V_{F_s}` (because `I_{F_s} ⊂ I_F`: `pi_{F_s}` factors through
`pi_F`). QED

**Corollary 4.2 (Sylvester, binary families).** In the setting of loctheory Cor 5.3 (Sylvester, `n` large, `j >= 9`),
let `𝓕_bin` be the binary column families of `F_2`-rank `t` with `2^t ln 2n < j`. A product `P = prod <z_i,.>` lies in
`V_{𝓕_bin}(m')` if (i) every column is touched by `<= 2j` generators and (ii) every `F_2`-circuit `s` of columns with
`|s| <= 8` is touched (some `a in s` with `z_i(S_a) != 0`) by `<= j|s| - 1` generators; it lies in `V_{𝓕_bin}(m')` only if
(ii) holds with `<= j|s|`.

*Proof.* Proposition 4.1 for each family, with loctheory Cor 5.3 (circuits of size `>= 9` have `alpha >= 9j > m'`). QED

So **all binary families together impose only count conditions**: the cost of a generator in a circuit is one unit
of `j|s|`, independent of how heavy the generator is. On average this is half as demanding as the column condition
(a generator whose Fourier support has size `sigma` touches a fraction `~ |s| sigma/n` of the `|s|`-circuits and a fraction
`sigma/n` of the columns; the budgets are `j|s|` and `2j`). The heavy generators forced by Proposition 3.2 are therefore
compatible with every binary family. What remains are the non-binary families, i.e. those containing complete XOR blocks.

### 4.2 Complete blocks: self-similarity

Let `U ≤ F_2^k`, `dim U = t`, `F_U = U \ 0` (all `2^t - 1` columns of `U`), and assume `2^t ln 2n < j` so that
`B_{F_U} = Y_t(j)` (loctheory Lemma 5.1). Put `V = U^perp`; its cosets `C_zeta` (`zeta in U^* = F_2^t`) are the classes of
`F_U`. For a generator `z` let `c(z) = (z(C_zeta))_zeta in Z^{2^t}_0` be its class-sum vector; then
`pi_{F_U}(z) = M c(z)` with `M_{a,zeta} = [a.zeta = 1]`.

**Lemma 4.3 (restriction to a block).** Let `N_U = {i : c(z_i) != 0}` and `g_i(t) = <M c(z_i), t>`.
(a) If `prod_{i in I} g_i in J(Y_t(j))^perp` for every `I ⊂ N_U`, then `P in V_{F_U}(m)`. (b) If `P in V_{F_U}(m)`, then
`prod_{i in N_U} g_i in J(Y_t(j))^perp`; if moreover the forms `<z_i,.>` (`i <= m`) are linearly independent, then every
sub-product lies in `J(Y_t(j))^perp`. So membership in `V_{F_U}` depends only on the **class-sum vectors** of the
generators, i.e. on a product of `|N_U|` linear forms on the class space of the block. For a coset difference `z = 1_{x+H} - 1_{y+H}`:
* if `H ⊂ V`: `c(z) = |H| (e_{[x]} - e_{[y]})` is a multiple of a class pair ("light in `U`");
* if `H ∩ V = H_0`: `c(z) = |H_0| (1_{[x]+H'} - 1_{[y]+H'})`, a coset difference in `U^*` with `|H'| = |H : H_0|` ("heavy in `U`").
The minimal multiple of `M c(z)` in the lattice `Lambda_U` has entries `+-|H'|` on the touched axes.

*Proof.* (a), (b) are the proof of Proposition 2.6 written for `F_U`: the homogeneous components of `P(u_0 + Vt)` are
combinations of the sub-products with coefficients `prod_{i notin I} <z_i,u_0>`; the top one is the full product, and
with independent forms the coefficients are independent monomials in `u_0`, so each sub-product can be isolated. `pi_{F_U}(z)_a = z(S_a) = sum_{zeta : a.zeta = 1} z(C_zeta) = (M c(z))_a`. The class sums of
`1_{x+H}` are `|H ∩ V|` on the cosets of `V` meeting `x + H` and `0` elsewhere; those cosets form `[x] + H'` in `U^* = F_2^k/V`.
`M` is injective on `Z^{2^t}_0` (loctheory Prop 5.5, `|det M| = I_t`), so `M c' in Lambda_U` forces `c' = c(z)/g` integral,
and `c(z)/|H_0|` is primitive. QED

So the block condition is **the same problem one level down**: products of linear forms on the class space
`Z^{2^t}_0`, tested against the box `C_t(j) = {c : |(M c)_a| <= j}`, which is a Sylvester `T = 1` box on `2^t` points
**without** the cube constraint. Light-in-`U` generators enter as class pairs (braid forms of the block), heavy ones as
coset differences of the block.

**The dichotomy.** For a zonotope `sum_{i in I} [0,1] g_i` with `g_i = M c(z_i)/|H_0|` to fit in `Y_t(j)` one needs, on every
axis `a`, `sum_i |H'_i| [a touched] <= 2j`. If membership of the products in `J(Y_t(j))^perp` in degrees `>= alpha(Y_t(j))`
required this **weighted** count (the zonotope cost), then heavy generators would pay `|H'|` per axis in every block in
which they are heavy, and the averaging argument of Theorem 3.1 / Proposition 3.2 (Prop 5.3 of fake.md) would come back:
certificates would be impossible for gap-admissible data, fakes and actual clusters alike. If instead membership needs
only the **plain** count (as for widths, Prop 2.6(c), and for parity circuits, Prop 4.1), heavy generators are as cheap as
light ones and the design question of §3.4 stays open. **Which of the two holds is a finite computation per `(t, j)`.**
§5.3 does it for `t = 3`, `j = 2`.

## 5. The barrier test on fake pattern boxes and on the complete Fano block

### 5.1 What can and cannot be computed

The task asks for the base locus of the local system of fake pattern boxes "at regime-like ratios `m'/(2j) ~ 4`".
**No exact fake model reaches that ratio.** The order must satisfy `m <= rho(A)` (else `TC = C^d`), and in every
computable model `rho/(2L) <= 2`:

| model (`T = 1`) | `n` | `L` | `rho` | `rho/(2L)` | source |
| --- | --- | --- | --- | --- | --- |
| `sdc7`, `pal7`, `Syl(8)` | 8 | 2 | 8, 7, 8 | 2, 1.75, 2 | `t_hilb.py` (reproduces fake §3) |
| `Syl(8)` = `M8'` | 8 | 3 | 11 | 1.83 | same |
| `sdc11`, `pal11`, Paley-12 | 12 | 1 | 3 | 1.5 | fake `sizes_out.txt` |
| `sdc11`, `pal11`, Paley-12 | 12 | 2 | 8, `>= 8`, 8 | 2 | fake §3.3, zetc §5.5 |
| Sylvester `n = 16` | 16 | 2, 3 | `<= 15`, `<= 22` | `<= 3.75` | zetc §5 (no LP-supercritical Sylvester model below `n = 32`) |

The regime ratio is reached only asymptotically (`m' ~ 8j`, `rho ~ 0.64 d`). Two substitutes are used below:
1. the exact small models at their own orders (LP-subcritical), where the questions of §2-§4 can be answered exactly
   (§5.2);
2. the complete-block images `Y_t(j)`, which by loctheory Lemma 5.1 are the **exact, `n`-independent** images of `A`
   under complete XOR blocks for every large Sylvester `n`, and which are the place where §4.2 locates the decisive
   question (§5.3). Here the relevant ratio is internal to the block (`e/j` between `alpha = 3j` and `rho(Y_t(j))`), and
   it is regime-faithful at any `j`.

Symmetry (zetc Lemma 4.3: symmetry-fixed subspaces see only degenerate points) is used only to *organise* the search;
no search below is restricted to a fixed subspace: every multiset of generators in the stated pool is tested.

**Certification.** A product `P` is in `J(A)_m^perp` iff `(p_a/multinom(m;a))_a` lies in `T_m = {mu(c) : c in V(A,m)}`
(§1.2). `T_m` mod `p` is computed from the mod-`p` kernel of the evaluation matrix; since that kernel contains the
reduction of the rational kernel, `T_m(mod p) ⊃ T_m(Q) mod p`, so **a mod-`p` non-membership is a rigorous
non-membership over `Q`** (all products have `p`-integral coefficients). Memberships are re-checked exactly over `Q`
(`fmpq` rank of the augmented system). An `int64` overflow in a first version of `tangent_space_mod` was found by the
validation script `zcheck.py` (explicit zonotopes must be members) and fixed; all numbers below are from the fixed code.

### 5.2 Exact models: local systems and product certificates

**(a) Locality of the fake models at their own orders** (`locsys.py`; widths over all directions in `{-2..2}^d`, column
families of size `<= 2, 3`; `dim J(A)_m` in `R = S/(x^3)` from exact ranks; "holds" = exact sandwich, "deficit" = mod `p`,
evidence):

| model | `m` | `H_A(m)` | `dim J(A)_m` | `b_0 = 1` | `b_0 = 2` (columns) | `b_0 = 3` (columns) |
| --- | --- | --- | --- | --- | --- | --- |
| `Syl(8,2) = M8` | 8 | 42 | 315 | deficit 21 | deficit 21 | **holds** (reproduces judge J6(c)) |
| `sdc7`, `L = 2` | 8 (`= rho`) | 1 | 356 | **holds** | holds | holds |
| `sdc7`, `L = 2` | 7 | 43 | 350 | deficit 56 | deficit 56 | **holds** |
| `pal7`, `L = 2` | 7 (`= rho`) | 95 | 298 | deficit 56 | deficit 50 | **holds** |
| `pal7`, `L = 2` | 6 | 241 | 116 | deficit 34 | deficit 33 | deficit 12 (mixed families not tried) |
| `Syl(8,3) = M8'` | 11 (`= rho`) | 7 | 70 | **holds** | holds | holds (agrees with loccompute) |
| `Syl(8,3) = M8'` | 10 | 63 | 98 | **holds** | holds | holds |

So in the two `n = 8` fake models the tangent cone at the top order **is** the base locus of an explicit local system
(columns triples; for `sdc7` at `m = 8` even the width system alone).

**(b) Product certificates** (`prodsearch.py`; pool: balanced generators of support `<= 4`, entries `±1`; every
multiset passing the width counts over all test directions is tested; `zonofit.py` classifies certificates):

| model | `m` | multisets passing the counts | certificates | of which zonotopes | `TC` has nondegenerate points? |
| --- | --- | --- | --- | --- | --- |
| `sdc7` | 8 (`= rho`) | 0 | 0 | — | **yes** (certified, fake §3.2; = base locus of the width system by (a)) |
| `sdc7` | 7 | 24 | 24 | 24 | no (certificates) |
| `pal7` | 7 (`= rho`) | 57 | 57 | 57 | no (certificates) |
| `Syl(8,2) = M8` | 8 (`= rho`) | 0 | 0 | — | no (R11 J4: all real `TC` points degenerate), **without** a certificate |
| `Syl(8,2)` | 7 | 696 | 696 | 696 | no |
| `Syl(8,3) = M8'` | 11, 10 | 0, 0 | 0 | — | no at 11 (R11 J4) |
| `pal11` (pairs only) | 8 | 0 | 0 | — | open (fake §3: floating evidence only) |

(Certificates: mod-`p` membership for all, exact `fmpq` membership re-checked for the first five of each model; every one
is a genuine zonotope by an exact vertex check, hence an exact certificate by fake Lemma 5.1(1). `pal11`: pool restricted
to the 66 pairs because the support-4 pool did not finish in 50 minutes; its 8-pair products all violate a width count.)

**Reading.**
* In the computable range, product certificates exist exactly where the box contains a zonotope of length `m`
  (`sdc7`, `pal7`, `Syl(8,2)` at `m = 7`: the braid zonotopes found by referee_fake), and every certificate found is a
  zonotope. Nothing beyond zonotopes appears.
* Where the box is zonotope-free at order `m`, no certificate exists in the pool. Theorem 3.3 excludes all certificates
  for `m > 7` (`sdc7`), `m > 8` (`M8`), `m > 10` (`M8'`); the search adds `M8` at `m = 8` and `M8'` at `m = 10`.
* Certificates are sufficient, not necessary: `M8` at `m = 8` has no nondegenerate real tangent-cone point (R11 J4,
  through a saturation argument on irreducible factors) but no product certificate.
* `sdc7` at `m = 8` is the barrier-test outcome in miniature: the local (here: width) system's real base locus contains
  certified nondegenerate, well-separated points, so fake Thm 2.1 yields gap-admissible `sdc7` fakes whose directions
  have `kappa <= 2 m theta e^{-W/4}` (fake Cor 2.3 even puts fakes on `Z_P`). This happens in an **underdetermined** system (`H_A(8) = 1`) at
  ratio `m/(2L) = 2`.

**The Diophantine condition (2.2).** (2.2) quantifies over the relevant `lambda` with `||lambda||_1 <= r/rho_1 + 1`, a finite
set for each fixed box; every nondegenerate point therefore satisfies (2.2) with some finite `a = a(u*, n)`. The content of
(2.2) is that `a` stays bounded along a sequence of boxes with `n -> infinity`, which no computable model can test.
For the certified `sdc7` point the question is empty at `n = 8`, and I make no claim about it.

### 5.3 The complete Fano block `Y_3(2)`: heavy generators pay the zonotope cost

`Y_3(2)` (`|Y| = 2605`, 7 axes = the 7 columns of a 3-dimensional `U`) is the image of `A` under a complete Fano block for
every large Sylvester `n` at `j = 2` (loctheory Prop 5.4). Exact Hilbert function (mod-`p` ranks, which are rigorous
upper bounds for `rho`, and agree with loctheory): `H = [1,7,28,84,210,455,812,756,252]`, `rho = 8`, against the grid
`[1,7,28,84,210,455,875,...]`: `alpha = 6 = 3j` (lines), and the non-parity generators start at `7 = 3.5j` (loctheory
Prop 5.4). Restricted generators (§4.2): **light** = class pairs, entries `±1` on 4 axes; **heavy** = coset differences
with `|H'| = 2`, entries `±2` on 2 axes, primitive in `Lambda_U`, zonotope cost 2 per touched axis. Scripts `yblock.py`
(stratified sampling), `yblock_exh.py` (exhaustive), `yblock_mix.py` (mixed).

| degree `e` | generators | condition | in `J(Y_3(2))^perp_e` / tested | status |
| --- | --- | --- | --- | --- |
| 5 (`< alpha`) | heavy, any axis count `<= 4` | — | 3000 / 3000 | sampled; agrees with Prop 2.6(a) |
| 6 (`= 3j`) | heavy, lines `<= 5`, axis count 3 or 4 (zonotope cost 6, 8 `> 2j`) | parity-feasible | 819 / 819 | sampled: **count cost suffices** |
| 7 (`= 3.5j`) | heavy, lines `<= 5`, max axis count 2 (zonotope cost `= 2j`) | parity-feasible | 1536 / 3072 | **exhaustive** |
| 7 | heavy, lines `<= 5`, max axis count 3 or 4 (zonotope cost `> 2j`) | parity-feasible | **0 / 88200** | **exhaustive; rigorous** |
| 7 | light, max axis count 4 (zonotope cost `= 2j`) | — | 24 / 394 | sampled |
| 8 (`= rho`) | heavy, any | — | 0 / 3000 | sampled |
| 7 | mixed: `h` heavy + `7 - h` light (`h = 0..7`), lines `<= 3j` | zonotope cost `<= 2j` | 25 / 324 | sampled (`yblock_mix.py`) |
| 7 | mixed, same | zonotope cost `> 2j` | **0 / 4705** | sampled; each non-membership rigorous |

(The exhaustive row: all 91272 multisets of 7 heavy restricted generators with every axis touched `<= 2j = 4` times and
every line touched `<= 3j - 1 = 5` times, i.e. every grid and every line (parity) condition of Proposition 4.1 satisfied.
Non-memberships are certified mod `p` (§5.1), hence over `Q`; memberships are mod-`p` evidence.)

**Control (`t = 2`, a single line, `yblock_t2_j2_out.txt`, `yblock_t2_j3_out.txt`).** For the line image `P_s` Proposition 4.1
proves count costs below degree `jc = 3j`. Sampled heavy products (axis-parallel, entries `±2`, zonotope cost 2 per axis)
with axis counts up to `2j` are all in `J(P_s)^perp` below `3j` (`j = 2`, `e = 5`: 2000/2000; `j = 3`, `e = 8`: 2000/2000),
and at `e = 3j` only the centred monomial `t^{(j,j,j)}` survives, for even `j` (`j = 2`: 252/252 of the count-2 samples,
`0/1748` otherwise; `j = 3`: `0/2000`, the orbit being killed by the sign `lambda_j = (-1)^j` of loctheory Lemma 3.2).
The method thus reproduces the proved behaviour of parity families.

**Reading.**
* Up to degree `3j` the block behaves like a parity family: costs are counts (Proposition 4.1).
* At degree `3.5j`, where the order-4 relation of the complete block enters, **no** product of heavy generators with an
  axis touched more than `j` times survives, although every grid and line condition holds; heavy products survive only
  when their zonotope cost `sum_i |H'_i| [a touched] <= 2j`. Light generators (zonotope cost = count) still reach count
  `2j`. So in this block the cost of a generator is its **zonotope cost**, not its count, from degree `3.5j` on.
* This is one block at one `j`; whether the same holds for `Y_t(j)`, all `t >= 3`, in the degree range
  `[alpha'_t, rho(Y_t(j))]` with `alpha'_t ~ 2jb/2^{t-1}` (loctheory §5.3) is open. `Y_3(3)` (`~5 * 10^4` points, degrees up to
  `~13` in 7 variables) is the next exact case and needs the `GL(3,2)` isotypic decomposition.

**Consequence (heuristic, conditional on the pattern persisting).** Consider (`W_t`): *for every complete `t`-block `U` with
`2^t ln 2n < j` and every product of `N >= alpha'_t` restricted generators lying in `J(Y_t(j))^perp`, the zonotope cost on
each axis is `<= 2j`.* Let all generators be coset differences `z_i = 1_{x_i+H_i} - 1_{y_i+H_i}`, `|H_i| = 2^{s_i}`. A column `a`
touched by `z_i` has `a in H_i^perp`, and for a uniformly random `t`-block `U` containing `a`, `z_i` is heavy in `U` with
`|H'| = |H_i| = |z_i(S_a)|` as soon as `H_i ∩ U^perp = 0`, which has probability `>= 1 - 2^{s_i - t + 1}` (`U^perp` is a random
subspace of codimension `t - 1` of `a^perp ⊃ H_i`). For `t >= s_i + 1` almost all generators touch `U`, so `N_U ~ m' > alpha'_t`
and the block is active. Averaging (`W_t`) over the blocks containing `a` gives
`sum_i |z_i(S_a)| (1 - 2^{s_max - t + 1}) <= 2j` for every column, hence `sum_i X_{z_i} <= (2 + o(1)) j W` when
`2^{s_max} = o(2^t)`, against `sum_i X_{z_i} >= m'(W/2 - g) >= (4 - o(1)) j W` from the gap inequalities (Proposition 3.2):
**coset-difference certificates would be impossible for gap-admissible data** if (`W_t`) holds for some `t` with
`2^{s_max} = o(2^t)` and `2^t ln 2n < j`. With `Loc(b_0)` for a *fixed* `b_0`, however, only blocks with `2^t - 1 <= b_0` are
local families, and generators with `2^{s} >> b_0` touch them rarely (Fourier support of size `n/2^{s+1}`, and
`U ∩ H^perp = 0` generically), so `Loc(b_0)` alone allows certificates made of very heavy generators. The logical
position is therefore:

```text
Loc(b_0)  and  (W_t) for some t with 2^{s_max} = o(2^t)   ==>   no coset-difference certificate with those s_i ,
```

because under `Loc(b_0)` a local certificate is a certificate (Corollary 2.7), hence lies in every `V_{F_U}`. So the
certificate route of Corollary 6.2 needs **either** `Loc(b_0)` together with a failure of the zonotope-cost behaviour in
large complete blocks, **or** generators heavier than every active block. The `Y_3(2)` computation is the first data
point, and it shows the zonotope cost.

## 6. Conditional theorems and what remains

Setting of §6.1-6.2: an actual cluster of `n` points on `x^2 + y^2 = R^2` (arc length `C sqrt R`, `C' = max(C,1)`),
`d = n - 1`, sign columns, `T = 1` slab box `A`, `m' = ceil(4w)`, regime `d >= 316 ln 4r`, `W > W_1` (so rigidity and
Step 0 hold: `kappa_c(hat u) <= eps_0 = C d e^{-W/4}/(m'+1)` for every integer `c in V(A,m')` with `|c|_inf <= H = |A|^2`;
t1rigour Thm 4.2). `Theta <= C e^{-W/4}` is the angular width (nearstrata Lemma 3.1), `u_x = delta_x - delta_0`,
`hat u = u/|u|_inf`.

### 6.1 Certificates imply the contradiction (Theorem 6.1)

For `z in Z^d` let `tilde z in Z^n_0` be its balanced extension (`tilde z_0 = -sum z`), so `<z,u> = <tilde z, delta>`, and
`X_z = sum_p log p |tilde z(S_p)|`.

**Theorem 6.1.** Suppose there are `z_1..z_{m'} in Z^d`, an integer `c in V(A,m')` with `|c|_inf <= H` and a rational
`lambda != 0` such that
* (C) `Phi_c(u) = lambda prod_i <z_i,u>` (a product certificate with explicit `c`), with `psi := ln(|c|_1/|lambda|)`;
* (Ar) `|<tilde z_i, delta>| >= e^{-X_{z_i}/2}` for every `i` (the character gap inequality; it holds for every
  gap-admissible fake by definition, and for pairs it is nearstrata Lemma 3.1; for general characters of an actual
  cluster it is the same Gaussian-integer argument applied to `prod z_x^{tilde z_x^+} - prod z_x^{tilde z_x^-}`, valid
  only when `<tilde z_i,delta> != 0`, and I take it here as a **hypothesis**).

Put `E = sum_i (X_{z_i} - W/2)`. Then

```text
W/4  <=  E/2 + psi + m' ln(C d) + ln(C d/(m'+1)) .
```

So a cluster with `W/4 > E/2 + psi + m' ln(Cd) + ln(Cd)` and such a certificate does not exist.

*Proof.* `c in V_H`, so Step 0 gives `|lambda| prod_i |<z_i,hat u>| = |Phi_c(hat u)| <= eps_0 |c|_1 d^{m'}`. By (Ar),
`|<z_i, hat u>| = |<tilde z_i,delta>|/|u|_inf >= e^{-X_{z_i}/2}/Theta >= C^{-1} e^{-(X_{z_i} - W/2)/2}`. Multiply over `i`,
take logarithms and insert `eps_0`. QED

**Corollary 6.2 (the conditional growth statement).** Assume:
1. `Loc(b_0)` (§1.4);
2. (*local certificate*) for every regime profile with `n >= n_1` there are `z_1..z_{m'}` with
   `P = prod <z_i,.> in V_{𝓕_{b_0}}(m')` (sufficient: the counts of Proposition 2.6(a) for every family of `<= b_0`
   directions, e.g. Proposition 4.1 for parity families);
3. (*quantitative locality*) the `c` that `Loc(b_0)` then provides (`P in J(A)_{m'}^perp = span{Phi_c}`, Lemma 2.1) can be
   chosen integral with `|c|_inf <= H` and `psi <= (1 - eta) W/8`;
4. (*tight characters*) `E <= (1 - eta) W/4`, and (Ar) holds for every `z_i` and every cluster (or: a finite family of
   certificates is given whose hyperplane arrangements intersect inside the braid arrangement, and (Ar) holds for the
   certificate whose factors do not vanish at `delta`).

Then every regime cluster with `n >= n_1` has `eta W/4 <= m' ln(Cd) + ln(Cd) = O(sqrt(n log n) log n)` (for `C <= 1`;
for `C > 1` add `m' ln C`), which is false because `W = log N = 2 log R >= (3 - o(1)) n ln n` by the proved bound
`n <= (2/3 + eps) log R/loglog R`. So regime clusters have `n < n_1`, and the regime condition `d >= 316 ln 4r` gives
`n = O(loglog R)` on all norms, exactly as in nearstrata Cor 4.2 (with `MIN^char` in place of `MIN^pair`: that proof uses
only Step 0 and a lower bound for `kappa(hat u)` valid at the actual `hat u`).

*Proof.* Items 1-2 and Corollary 2.7 give `P in J(A)_{m'}^perp`; item 3 gives (C) with the stated `psi`; Theorem 6.1 with
items 3-4 and `m' ln(Cd) = O(sqrt(n log n) log n)`. QED

**Remarks.**
* Items 3 and 4 are **not** consequences of locality: `Loc` is a statement about spans and gives no height for `c`, and
  item 4 is an arithmetic property of the profile. In Sylvester-type (Plotkin-tight) profiles coset differences have
  `X_z = (W/2)(1 + O(1/n))` when the column weights are equal, so `E = O(m' abar)`; for general profiles heavy characters
  can have `X_z` up to `|z|_1 W/4`, and item 4 then fails.
* By Theorem 3.1 the generators cannot all be pairs, so any proof along these lines uses **non-pair characters**, i.e. a
  hypothesis of type `MIN^char` (lower bound for `kappa` on the region cut out by character inequalities) rather than
  `MIN^pair`. `MIN^char` is implied by `MIN^pair` only if the character region is inside the pair region, which it is
  (every pair is a character), so `MIN^char` is **weaker** than `MIN^pair` as a hypothesis and also suffices; but its
  inputs (Ar) are not pair data and can fail exactly for actual clusters with `<tilde z_i,delta> = 0`.
* For gap-admissible fakes (Ar) holds by definition, and the bound on `E` in item 4 is a property of the fake profile.
  If items 2-4 held for a fake pattern box, its fakes would violate Step 0 and `TC_H` would have no nondegenerate point:
  the fake test would be passed for that box. This is consistent with fake Thm 2.1 (no contradiction either way).

### 6.2 An unconditional reduction for Sylvester-type profiles: block products

This subsection does not use `Loc`. It shows that for profiles containing many disjoint copies of a small Sylvester
structure, `MIN^pair` follows from **one fixed finite statement** about a single small box.

**Lemma 6.3 (block sub-boxes).** Any sign profile. Let `X_1..X_B ⊂ [n]` be disjoint, `j_1..j_B >= 0` with
`sum_b j_b <= j`, and `A_b = { k in {-1,0,1}^n : supp k ⊂ X_b, sum k = 0, |k(P_p ∩ X_b)| <= j_b for all p }`. Then
* (a) `A ⊃ A_1 ⊕ ... ⊕ A_B := { sum_b k_b : k_b in A_b }`, and `nu_A(u) >= sum_b nu_{A_b}(u)` for every `u`;
* (b) for `c_b in V(A_b, m_b)` the convolution `c = c_1 * ... * c_B` (`c(sum_b k_b) = prod_b c_b(k_b)`) lies in `V(A, sum m_b)`,
  `|c|_1 = prod |c_b|_1`, `|c|_inf = prod |c_b|_inf`, and `Phi_c(u) = (m!/prod_b m_b!) prod_b Phi_{c_b}(u)`;
* (c) for the Sylvester profile `n = 2^k`, if the `X_b` are cosets of a subgroup `U` of order `N_0` not containing the base
  point, then `A_b ≅ Syl^box(N_0, j_b)` (defined below; after translating `X_b` to `U`).

*Proof.* (a) For `k = sum k_b` (disjoint supports): `k in {-1,0,1}^n`, `sum k = 0`, and
`|k(P_p)| <= sum_b |k_b(P_p ∩ X_b)| <= sum_b j_b <= j`, so `k in A`. `nu` is monotone under inclusion (zetc §2.2(1):
`A' ⊂ A` gives `J(A) ⊂ J(A')`) and additive on products (the direct sum is a product in disjoint coordinates; zetc §2.2(1),
or the bigraded argument of Proposition 2.2(b)). (b) is fake Lemma 5.1(3) (the map `(k_b) -> sum k_b` is injective, so
`|c|_1` and `|c|_inf` multiply). (c) For `x = x_b + y` (`y in U`), `[a.x = 1] = [a.x_b = 1] + [a.y = 1] mod 2`; so
`P_a ∩ X_b` is `∅`, `X_b`, or `x_b +` (an affine hyperplane of `U` given by `a|_U != 0`), and for balanced `k_b` the
constraint `|k_b(P_a ∩ X_b)| <= j_b` is the Sylvester column constraint of `a|_U` (the complement has the opposite sum). QED

Write `Syl^box(N_0, j_0) = { k in {-1,0,1}^{F_2^{t_0}} : sum k = 0, |hat k(a)| <= 2 j_0 (a != 0) }` (all `N_0` coordinates
cube-constrained; for a coset `X_b` not containing the base point of the big box, `A_b ≅ Syl^box(N_0, j_b)` by (c)).

**Hypothesis `BLOCK(t_0, j_0, m_0; K, c_0, H_0)`.** `N_0 = 2^{t_0}`, `m_0 = (8 + delta_0) j_0` with `delta_0 > 0`, and for every
`v in R^{N_0}` with `max_x |v_x - v_{x_0}| = 1` (for some base `x_0`) and `min_{x != y} |v_x - v_y| >= sigma > 0`, there is an
integer `c in V(Syl^box(N_0, j_0), m_0)` with `|c|_inf <= H_0` and `|Phi_c(v)| >= c_0 sigma^K |c|_1`.

**Lemma 6.5 (`BLOCK` is qualitative).** For fixed `(t_0, j_0, m_0)` with `m_0 <= rho(Syl^box(N_0,j_0))`, there are
`K, c_0, H_0` with `BLOCK(t_0, j_0, m_0; K, c_0, H_0)` **iff** `ZE-TC^braid` holds for the block: every real point of
`TC(Syl^box(N_0, j_0), m_0)` (in the coordinates relative to `x_0`, all `x_0`) has two equal coordinates among the `N_0` positions.

*Proof.* "Only if": a collision-free real tangent-cone point has `sigma > 0` and `Phi_c(v) = 0` for all `c`. "If": fix a
finite integer basis `c^(1)..c^(D)` of `V(Syl^box, m_0)`, `H_0 = max |c^(i)|_inf`, and `F = sum_i Phi_{c^(i)}^2`. On the compact
set `K_1 = {max_x |v_x - v_{x_0}| = 1}` Lojasiewicz's inequality for the polynomial `F` gives
`F(v) >= C dist(v, F^{-1}(0))^{K'}` (`v in K_1`). The real zero set `F^{-1}(0)` is the cone `TC`, contained in the braid
arrangement `Br` by hypothesis, so `F(v) >= C dist(v, Br)^{K'}`, and `dist(v, Br) >= sigma/2` (sup norm) if all gaps are
`>= sigma`. Then some `i` has
`|Phi_{c^(i)}(v)| >= (C/D)^{1/2} (sigma/2)^{K'/2}`, and `|c^(i)|_1 <= max_i |c^(i)|_1`. QED

So the constants in `BLOCK` are ineffective but harmless: they are attached to one fixed box. `ZE-TC^braid` for the block
is false if `m_0 > rho_A(coordinate columns of U)` (`<= 2 j_0 t_0`), because the `t_0` coordinate columns span a
collision-free rational stratum (README R11 item 7); so `t_0 >= 5`, `N_0 >= 32`. Blocks with `N_0 <= 16` are
LP-subcritical (zetc §5), so `BLOCK` is beyond exact computation.

**Theorem 6.4 (block products give `MIN^pair` for Sylvester profiles).** Assume `BLOCK(t_0, j_0, m_0; K, c_0, H_0)`
(equivalently, by Lemma 6.5, `ZE-TC^braid` for the single box `Syl^box(2^{t_0}, j_0)` at order `m_0 = (8 + delta_0) j_0`). Let
`B = floor(m'/m_0)`, `r = m' - B m_0 < m_0`, and `Q_1 = (m_0 + K) B binom(N_0, 2)/2 + 2.1 r`. Then there is `n_0` such that
for every Sylvester profile with `n = 2^k >= n_0` (sign columns, any multiplicities), every `C`, every pair data `eps_xy >= 0`
with `sum eps_xy <= N_2 abar`, and every `u` in the gap region `G = {|u|_inf = 1, |u_x - u_y| >= e^{-eps_xy/2}}`:

```text
ln(1/kappa(u))  <=  m' ln d + B ln(2^K/c_0) + r ln 2 + Q_1 abar (1 + O(N_0/n)) .
```

Consequently `MIN^pair` holds for `(A, m')` whenever `W (1/4 - Q_1/(2d)) > m' ln d + ln(Cd) + B ln(2^K/c_0) + r ln 2 + 2 Q_1 ln C'`;
since `B = O(j/m_0)`, `Q_1 = O(j N_0^2 (1 + K/m_0))` and `j = O(sqrt(n log n))`, this is the range
`W >= W_2 := (4 + o(1)) m' ln d + O_{BLOCK}(j (1 + ln C'))`, which contains the near-extremal range `W >= 3 n ln n - 14 n`.

*Proof.* *Room.* `B j_0 + r <= j` and `B N_0 + 2r <= n - N_0` for large `n`: `B j_0 <= m' j_0/m_0 <= (8j + 9)/(8 + delta_0)`, and
`(8j+9)/(8+delta_0) + m_0 <= j` once `j delta_0 >= 9 + m_0(8 + delta_0)`.
*Decomposition.* Choose a subgroup `U` of order `N_0`, `B` distinct cosets `X_1..X_B` of `U` not containing the base point,
and `r` disjoint pairs `{x_l, y_l}` of non-base points outside them. Apply Lemma 6.3 with the blocks `A_b ≅ Syl^box(N_0, j_0)` and with the pair zonotope
`Z = {sum_{l in I} (e_{x_l} - e_{y_l})}` (slab budget `r`; its alternating indicator `c_Z` has `|c_Z|_1 = 2^r`,
`Phi_{c_Z} = ± r! prod_l (u_{x_l} - u_{y_l})`, fake Lemma 5.1(1)).
*Blocks.* For `u in G` and a block `X_b` with base `x_0`, put `s_b = max_{x in X_b} |u_x - u_{x_0}| in (0, 2]`,
`g_b = min_{x != y in X_b} |u_x - u_y|`, and `v = (u|_{X_b} - u_{x_0})/s_b`. Then `v` is normalised with gaps `>= g_b/s_b >= g_b/2`,
and `BLOCK` gives `c_b` with `|Phi_{c_b}(u)| = s_b^{m_0} |Phi_{c_b}(v)| >= c_0 2^{-K} g_b^{m_0 + K} |c_b|_1` (`s_b >= g_b`;
`Phi_{c_b}(u)` depends on `u|_{X_b}` modulo constants since `sum k_b = 0`).
*Product.* `c = c_1 * ... * c_B * c_Z in V(A, m')` has `|c|_inf <= H_0^B <= H` for large `n`, so `c in V_H`, and
`kappa(u) >= |Phi_c(u)|/(|c|_1 d^{m'}) >= d^{-m'} prod_b (c_0 2^{-K} g_b^{m_0+K}) 2^{-r} prod_l |u_{x_l} - u_{y_l}|` (the
multinomial is `>= 1`).
*Gaps.* On `G`, `ln(1/g_b) <= (1/2) sum_{x<y in X_b} eps_xy` and `ln(1/|u_{x_l} - u_{y_l}|) <= eps_{x_l y_l}/2`.
*Averaging.* Take `U` uniform among subgroups of order `N_0` (then `P(x + y in U) = (N_0-1)/(n-1)` by transitivity of
`GL(k,2)` on nonzero vectors), the `B` cosets uniform among the `n/N_0 - 1` cosets avoiding the base point (so a given
coset is chosen with probability `<= B N_0/(n - N_0)`), and the `r` pairs uniform among disjoint pairs of the
`n' >= n/2` remaining non-base points. The expected value of `(m_0+K)/2 sum_b sum_{X_b} eps + (1/2) sum_l eps_l` is at most
`(sum_{x<y} eps_xy) [(m_0+K) B N_0(N_0-1)/(2(n - N_0)(n-1)) + r/(2 binom(n',2))] <= abar Q_1 (1 + O(N_0/n))`, using
`sum eps <= N_2 abar` and `N_2/binom(n',2) <= 4.2`. Some decomposition achieves at most the expectation; it depends only on the pair data, not on `u`.
*`MIN^pair`.* `ln(1/eps_0) = W/4 - ln(Cd/(m'+1))` and `abar = W/(2d) + 2 ln C'`; substitute. QED

**Corollary 6.6 (a restricted conditional theorem).** Assume `ZE-TC^braid` for one box `Syl^box(2^{t_0}, j_0)` at an order
`m_0 > 8 j_0` (Lemma 6.5). Then there is `n_0` such that no cluster (any `C`, any norm `N`, `W = log N`) whose split primes
realise exactly the Sylvester columns (`n = 2^k`, `r = n - 1`, `e_p = 1`) has `n >= n_0`. More generally the same holds for
every sign profile in the regime that admits a probability distribution on decompositions (`B = floor(m'/m_0)` disjoint
non-base `N_0`-sets on each of which every column trace is `∅`, everything, or a Sylvester column of `F_2^{t_0}` up to
complement, plus `r` disjoint pairs) under which every pair `{x,y}` lies in a common block with probability
`<= (1 + o(1)) B binom(N_0, 2)/N_2` and in a leftover pair with probability `O(r/N_2)`.

*Proof.* For `n >= 2900` the box is in the regime (`r = d`, t1rigour §2). By the proved bound
`n <= M <= (2/3 + eps) log R/loglog R` (two_thirds/proof.md) and `W = log N = 2 log R`, `W >= (3 - o(1)) n ln n >= W_2`
for large `n`, so Theorem 6.4 gives `MIN^pair` at the cluster's own data; Step 0 (t1rigour Thm 4.2), the pair-gap
separation (nearstrata Lemmas 3.1-3.3) and `MIN^pair` contradict rigidity (nearstrata Thm 4.1 with the referee's
`MIN^pair` form). QED

The second sentence is what the proof of Theorem 6.4 uses; for the Sylvester profile the distribution is the uniform
coset decomposition.

**What Theorem 6.4 says about the barrier.**
* For Sylvester-type profiles, `MIN^pair` at **every** large `n` (in the stated `W` range) reduces to **one** finite
  statement about one fixed LP-supercritical block. The regime is reached by products, which fake Lemma 5.1 already
  had, but embedded into the slab budget of the big box: the budget `j` is split among `~ m'/m_0` blocks.
* The mechanism needs blocks on which the profile has **few traces**. For a generic profile with many columns, every
  small set `X` has all `2^{|X|}` traces, `A_b` is the ball `{|k|_1 <= 2 j_b}`, `rho(A_b) <= 2|X|` and, by the pair-split
  count of Theorem 3.1, a block whose full cube lies in `A` has `|X| <= 4 j_b` hence `rho <= 8 j_b`: ratio `<= 8`, never
  `> 8`. So Theorem 6.4 does **not** apply to the fake pattern profiles (Paley-type `sdc_q`, `pal_q`, whose traces on small
  sets are full), nor to generic actual profiles. It isolates the property that would separate a profile from the fakes:
  small blocks with low trace complexity and an LP-supercritical tangent cone.

### 6.3 What remains

1. **`Loc(b_0)`** itself (§1.4): open; loctheory §5.3 and loccompute Prop 4.4 make it doubtful for every fixed `b_0`
   (complete XOR blocks, Paley congruences).
2. **Existence of local certificates** (Cor 6.2 item 2) in regime boxes. Proved: pairs are impossible (Theorem 3.1);
   generators must be heavy on average (Proposition 3.2); binary families cost only counts (Corollary 4.2). Open: the
   cost of heavy generators in complete blocks, which decides between "certificates possible" and "certificates
   impossible by the averaging argument" (§4.2, §5.3); the first exact case (`Y_3(2)`, degree `3.5j`) shows the zonotope
   cost. Also open: the local conditions from families of non-column integer directions (Proposition 2.6 applies, but
   `alpha_F` is not known for them).
3. **Quantitative locality** (Cor 6.2 item 3): no height bound for `c` is known even in exact models beyond zonotopes.
4. **Tight characters** (item 4): holds for Plotkin-tight (Sylvester-type) profiles, fails in general.
5. **`BLOCK`** for one LP-supercritical Sylvester block (Theorem 6.4): a finite statement on `N_0 >= 32` points,
   beyond exact reach; it would settle Sylvester-type profiles but not generic ones.
6. **Next exact steps.** `Y_3(3)` and `Y_4(2)` with the `GL(t,2)` isotypic decomposition (does the zonotope cost persist
   in `[alpha'_t, rho(Y_t(j))]`?); the support-4 certificate search for `pal11` (`L = 2`, `m = 8`) with the `Z_11` symmetry;
   `ZE-TC^braid` for the smallest candidate block `Syl^box(32, j_0)` is out of reach.

## Appendix: scripts and certificates

All in `fatpoint_checks/`; outputs `*_out.txt`. `p = 2147475331`. Mod-`p` non-memberships and mod-`p` rank
deficits are used only in the rigorous direction stated in §5.1; mod-`p` memberships are evidence unless marked exact.

| script | what | status |
| --- | --- | --- |
| `fplib.py` | models (`sdc_q`, `pal_q`, Sylvester, Paley-12; copied from round12 `fake_checks/flib.py`), boxes, reduced monomials, `T_m` mod `p`, product expansion | library (overflow bug fixed, see §5.1) |
| `t_hilb.py` | Hilbert functions of the small models (reproduces fake §3) | mod-`p` ranks (upper bounds for `rho`) |
| `zcheck.py`, `zdebug.py` | validation: explicit lattice zonotopes must lie in `T_m` | exact construction, mod-`p` test |
| `prodsearch.py` | exhaustive search of product certificates over a generator pool (balanced, support `<= 4`, entries `±1`), pruned by the width counts over all directions in `{-2..2}^d` (`{-1,0,1}^d` for `d >= 9`); membership mod `p`; survivors re-checked exactly over `Q`; zonotope test of every certificate found | negatives rigorous within the pool; positives exact |
| `run_small.sh` -> `prodsearch_small_out.txt` | `sdc7`, `pal7`, `Syl(8)` at `L = 2, 3`, orders `rho - 2 .. rho` | as above |
| `run_d11.sh` -> `prodcount_d11_out.txt` | `pal11`, `sdc11`, Paley-12 (`L = 2`, `m = 8`): width-count feasibility only | count pruning only |
| `zonofit.py` | zonotope test of printed certificates | exact |
| `uncert_budget.py` | Theorem 3.3 on the models and on Sylvester `n = 2^{12..20}` | exact integers / floating `w` |
| `locsys.py`, `run_locsys.sh` -> `locsys_out.txt` | is `J(A)_m` generated by width powers (`b_0 = 1`) or widths plus column families of size `<= 2, 3`? | exact via sandwich when equal; deficits mod `p` (evidence) |
| `yblock.py` | complete-block images `Y_t(j)`, Hilbert function, stratified sampling of products of light/heavy restricted generators (`t = 2, 3`) | negatives rigorous, sampled |
| `yblock_exh.py` | `Y_3(2)`, degree 7: **all** products of 7 heavy generators with axis counts `<= 4` and line counts `<= 5` (resp. `<= 6`) | exhaustive; negatives rigorous |
| `yblock_mix.py` | `Y_3(2)`, degree 7: mixed light/heavy products vs zonotope cost | sampled; negatives rigorous |
