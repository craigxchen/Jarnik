# Referee report on Astra's VIII.1 (the common extraction lemma)

Scope: lines 1748-1776 of `round31/astra_lemma_VIII/answer.md`, checked against the statement of (L_8) in
`expert_note/expert_note.md` §2, against `growth/ptolemy.md` §§1, 3(iv), 6 and 7, and against `round18/README.md` §1
("Balance", "Anchor"). Scripts are in `round31/ref8L_extract/`: `check_extract.py` (output `out.txt`) and
`paley_scale.py` (output `paley_out.txt`). I wrote both myself; I did not import Astra's `viii_check.py`.
Literature: none was needed for this block. I made no web requests and did not read `cg07.txt`.

## 0. Verdict

**VIII.1 is correct (Proved).** I re-derived every step and every constant. The escape from the fixed exceptional set
is sound. The datum that the lemma feeds to (L_8) is of exactly the kind that (L_8) admits, and it is also the kind
for which Theorem D derives (L_8) from V_8. There is no false step. Section 6 lists four wording corrections. None of
them changes the mathematics.

Summary of the block:

| Item | Label |
|---|---|
| The hypotheses give data admissible for the matrix form of (L_8), anchor 0, with Z unsymmetrised | Proved |
| The version needed is the matrix (L_8) with linear cw; (L_8^f), (L_8^pc) and (L_8^c)^cat do not suffice; (L_8^c)^all suffices | Proved |
| Theorem D yields that version: c = 6/69, eta_0 = 1/381, Z = Z_{1/43} | Proved in the record; proof re-read, not re-refereed |
| G_xy <= 4 lambda off the diagonal; E G >= 0; (VIII.2) | Proved |
| E m_A^2 = E_{u,v} g^{\|A\|} <= E G^2 for even \|A\|, repeated rows included | Proved; Computed (exact Fractions) |
| Expansion of the 128 pattern masses (127 summands); Markov step; 64-pattern bound on D_ij | Proved; Computed |
| Escape from Z by Schwartz-Zippel (d/M), repeated rows (28/M), union bound with 3/4 | Proved; Computed |
| Final contradiction 28(a/4+lambda_*)W < cW/256 against cW/128 - O(1); uniformity | Proved |
| Numerical margins a = 2e-5, lambda_* = 1e-16 (used in VIII.3) | Computed (exact rationals) |

## 1. The statement as I verify it (corrected wording)

**Lemma VIII.1 (corrected statement; Proved).** Assume (L_8) in the matrix form of expert note §2. That is, assume
there are c > 0, eta_0 > 0, C_0 >= 0 and a proper Zariski-closed Z in Mbar_{0,8} (all fixed, anchor = slot 0) such
that every datum (P, (n_T), (t_ij)) with Phi not in Z that is eta-balanced at scale w with eta <= eta_0 satisfies
H_Sigma >= c w - C_0. Fix a nonzero polynomial F in 8 variables, of total degree d, that vanishes on the preimage of
Z ∩ M_{0,8} under the ordered-configuration map (Section 4). Then there are lambda_* > 0, M_* and W_*, depending only on
(c, eta_0, C_0, d), with the following property.

Let X be a finite set of rationals with |X| = M >= 8. Let (Omega, mu) be a finite measure space of total mass W > 0, and
for each x in X let b_x : Omega -> {+1, -1} be measurable. Write D_xy for the mu-mass of {b_x ≠ b_y}. For an ordered
tuple (x_0, ..., x_7) and nonempty T in [7], the *cut T* is the set of atoms where b_{x_i} ≠ b_{x_0} exactly for
i in T. Let lambda <= lambda_*. Assume that every ordered tuple of 8 distinct elements of X comes with a datum of
(L_8) such that:

1. the columns are primitive, P_0 = (1,0)^T, the chain condition holds, and the point Phi of the columns in
   M_{0,8}(Q) equals the configuration [x_0 : ... : x_7]. Any Möbius normalisation applied to all eight points
   together is allowed, even one that depends on the tuple;
2. log n_T = mu(cut T) for every nonempty T in [7];
3. 0 <= log|t_ij| <= D_{x_i x_j}/2 - W/4 + lambda W for all 0 <= i < j <= 7, where t_0j = Y_j.

Then W >= W_* implies M <= M_*.

The original wording is correct. The corrected statement adds only the dependence of the constants, lambda <= lambda_*
for possibly negative lambda, and the Möbius remark in item 1.

## 2. Point (1): the datum is admissible, and which version of (L_8) it needs

**2.1 Blocks for all 127 T.** The expert note's (L_8) takes positive integers n_T for *every* nonempty T in [7],
singletons and the top block [7] included. Singletons enter no minor. The top block divides every det(P_i,P_j) with
i, j >= 1. The lemma assigns every block a mass, so its data have the shape that (L_8) prescribes.

Balance is required only on the 119 visible blocks (expert note) or on "the n_T" (ptolemy §7, round18 §1). The lemma
proves balance on all 128 patterns, including the pattern T = ∅, which is not a block. So both readings of the
balance hypothesis are met directly. The top-block inflation argument of round18 §1 ("Balance") is not needed.

**2.2 Anchor minors.** H_Sigma = Σ_j log|Y_j| + Σ_{1<=i<j<=7} log|t_ij| has 7 + 21 = 28 terms. Hypothesis 3 is
imposed for all 28 pairs, the anchor pairs (0, j) included. The cofactor sum is therefore taken over 28 pairs, which is
what the proof uses.

**2.3 Cofactors sharing primes with blocks; non-canonical blocks.** The expert note requires only t_ij ∈ Z \ {0} and
imposes no coprimality. The n_T are arbitrary positive integers subject to the chain condition. Item 532's original
data (ptolemy §1) have blocks that "may share rational primes", and residues t_ij that may be divisible by block
primes. So (L_8) as recorded allows both. Astra notes after VIII.2 that her square blocks need not be the canonical
contact integers d_T.

**2.4 Which version is needed.**
* The lemma needs the **matrix form (L_8)** (expert note §2 = ptolemy §7 = round18 §1), at anchor 0, with Z
  unsymmetrised, for arbitrary non-canonical blocks and arbitrary nonzero cofactors. If only a restricted class of
  data occurs, it needs (L_8) restricted to that class. This restriction clause is trivially valid.
* It needs the **linear** lower bound c w. The upper bound it proves is H_Sigma <= 28(a/4 + lambda) W, which is linear
  in W with a small coefficient. Therefore (L_8^f), with f(w) = o(w), does **not** suffice. Since (L_8) ⇒ (L_8^f), it
  is the stronger statement that is used.
* (L_8^pc) in lattice form and (L_8^c)^cat do not suffice either. Both concern canonical contact integers d_T, and
  the lattice form also needs pairwise-coprime cofactors of size e^{w(1±eta)}. The lemma's data satisfy neither
  condition.
* (L_8^c)^all suffices, because it implies (L_8) (expert note, with eta_0 = eta_1/2 and c = eta_1/4).

**2.5 Theorem D gives this version.**
* V_8 is applied at an arbitrary point P of X(Q) \ Z_eps, so Theorem D concerns Phi only.
* Lemma 6.2 uses (3.2). Its proof (ptolemy §3, item (iv)) uses only three things: primitivity of the columns (the
  Gromov products v_p(det) are Bruhat-Tits Gromov products from the standard vertex), the factorisation
  det = t ∏ n_T, and the two-chains condition. Primes dividing residues, anchor residues included, are absorbed by
  the bound |l_T(p) − v_p(n_T)| <= 2 max v_p(t_ij). Singletons and the top block are pendant edges and drop out.
* So the proof applies verbatim to abstract matrix data with non-canonical blocks and with cofactors that share
  primes. It gives (L_8) with Z = Z_{1/43}, eta_0 = 1/381 and H_Sigma >= ((7 − 381 eta) w − C')/69 >= (6w − C')/69.
* Lemma 6.2 needs balance only on boundary T. I re-read this proof; I did not re-referee it (it is refereed in
  `growth/referee_ptolemy.md`).

## 3. Point (2): the Gram argument, with every constant re-derived

Normalise mu to a probability P. Then G_xy = E_u b_x b_y = 1 − 2 D_xy / W.

* **G_xy <= 4 lambda for x ≠ y (Proved).**
  - M >= 8, so every pair x ≠ y occurs in slots (i, j) of some distinct tuple.
  - Hypothesis 3 with log|t| >= 0 gives D_xy >= (1/2 − 2 lambda) W, hence G_xy <= 4 lambda.
  - "Integrality" here means only that t is a nonzero integer, so |t| >= 1. The hypothesis forces the right-hand side
    to be >= 0. Identical sign functions (D = 0) are therefore excluded.
* **E_{x,y} G_xy = E_u (M^{-1} Σ_x b_x(u))^2 >= 0 (Proved).** Here x and y are independent and uniform on X, and x = y
  is allowed.
* **(VIII.2) (Proved).**
  - |G| <= 1, so E G^2 <= E|G|.
  - E|G| = 2 E G_+ − E G <= 2 E G_+.
  - E G_+ <= (1/M)·1 + (1 − 1/M)·max(4 lambda, 0).
  - Hence E G^2 <= 2/M + 8 lambda_+. The original wording requires lambda >= 0; for lambda < 0, read lambda_+.
* **E_x m_A^2 = E_{u,v} g(u,v)^{|A|} (Proved).** Here g(u,v) = M^{-1} Σ_x b_x(u) b_x(v) is in [−1, 1], and the slots
  of A are distinct coordinates of an independent uniform tuple in X^8. Coincident rows are allowed, because the
  computation is over the product measure on X^8. For even |A| >= 2 we have 0 <= g^{|A|} <= g^2, and
  E_{u,v} g^2 = E_{x,y} G_xy^2. The case |A| >= 4 holds because |g| <= 1. Exact check (Fractions, two random weighted
  systems, M = 4, |A| = 2, 4, 6, 8, all of X^|A|): equality holds, and so does the inequality against E G^2.
* **Pattern masses (Proved; Computed).**
  - The pattern T is {b_{x_i} = s_i for all i} ∪ {b_{x_i} = −s_i for all i}, with s_0 = +1 and s_i = −1 exactly for
    i in T.
  - Expanding ∏(1 + s_i b_i)/2 + ∏(1 − s_i b_i)/2 = 2^{-7} Σ_{A even} s_A b_A gives
    P(cut T) = 2^{-7}(1 + Σ_{∅≠A even} s_A m_A).
  - There are 2^7 − 1 = 127 nonempty even A. They correspond bijectively to the nonempty A' ⊆ [7] (add or drop slot 0),
    and s_A = (−1)^{|A'∩T|}.
  - So the pattern-mass vector is a Walsh transform. The script checks the identity exactly on 30 tuples. On 18,000
    sampled tuples it checks max_T |mass_T / w − 1| <= Σ|m_A| with exact integer masses.
* **Expectation and Markov step (Proved).**
  - By Cauchy-Schwarz, E Σ|m_A| <= 127 sqrt(E G^2) <= 127 sqrt(2/M + 8 lambda).
  - If this is < a/4, then P(Σ|m_A| > a) < 1/4.
  - On the complement, every one of the 128 masses is within a·w of w = W/128. With a < eta_0, the datum is
    eta-balanced with eta <= a <= eta_0 on all 127 nonempty T.
* **Separating patterns (Proved; Computed).**
  - For slots 0 and j, the separating T are those containing j: 64 of them.
  - For 1 <= i < j, the separating T contain exactly one of i, j: 2·32 = 64.
  - Hence D_{x_i x_j} <= 64(1 + a) w = (1 + a) W/2.
  - Then log|t_ij| <= (1 + a)W/4 − W/4 + lambda W = (a/4 + lambda) W.
  - The script found 0 violations of the 64-pattern bound in 18,000 tuples × 28 pairs.
* **Choice of constants (Proved).**
  - Pick 0 < a < eta_0 with 7a < c/256.
  - Then pick lambda_* > 0 with 28(a/4 + lambda_*) < c/256 and 8 lambda_* < (a/508)^2.
  - Then pick M_* with 2/M_* < (a/508)^2 − 8 lambda_* and M_* >= 4(28 + d).
  - Then W_* = 256 C_0 / c.
  - The condition "28a/4 < c/256" in the text is redundant, but harmless.

## 4. Point (3): the escape from Z (Astra asked for this first)

**Proved.**

* **Preimage of Z.** Let U ⊂ A^8 be the complement of the big diagonal, and let φ : U → M_{0,8},
  x ↦ [x_0 : ... : x_7], be the quotient by PGL_2. φ is a surjective morphism over Q.
  - Z is proper and closed in the irreducible Mbar_{0,8}, so Z ∩ M_{0,8} ≠ M_{0,8}.
  - φ^{-1}(Z ∩ M_{0,8}) is closed in U and is not all of U, because φ is surjective.
  - Its closure in the irreducible A^8 is a proper closed set, so a nonzero F vanishes on it. Its degree d depends
    on Z alone.
  - For a distinct tuple, Phi is in M_{0,8} (no boundary), so Phi ∈ Z ⇒ F(x) = 0.
* **The anchor.**
  - The anchor is the slot 0 of the *ordered* tuple, and the Z used is that of (L_8) for anchor 0. No symmetrisation
    is needed: every ordered tuple is sampled, and φ is applied to the ordered tuple.
  - In applications the columns realise the configuration after a Möbius map that depends on the tuple: y ↦ (y − y_0)/Q
    in VIII.2 and d ↦ (d/a − 1)/q with a = d_0 in VIII.3. The same map is applied to all eight points, so the class
    in M_{0,8} is still φ(x). This must hold in every application. A per-coordinate transform that is not a common
    Möbius map would break the escape.
* **Schwartz-Zippel for an arbitrary finite grid (Proved).**
  - Write F = Σ_k F_k(x_0..x_6) x_7^k with top coefficient F_e ≠ 0, where deg F_e <= d − e.
  - Then P(F = 0) <= P(F_e = 0) + e/M <= (d − e)/M + e/M by induction. Only |X| = M is used; X may be any finite
    subset of Q, or of any field.
  - Computed: for F = (harmonic cross-ratio of slots 0-3) × (cross-ratio of slots 4-7 = 2), which has d = 4, the exact
    probability on X^8 is 0.0049-0.0072 <= 0.1 = d/M for four 40-element grids: consecutive integers, random integers,
    {1/k} and squares.
* **Union bound (Proved).**
  - P(some x_i = x_j) <= C(8,2)/M = 28/M.
  - Hence P(balanced ∧ distinct ∧ F ≠ 0) >= 3/4 − (28 + d)/M > 0 for M > 4(28 + d).
  - The balance event is computed on X^8 with coincidences allowed. The hypotheses are invoked only for the tuple
    finally chosen, which is distinct. This order of operations is correct.
* **Comparison with Lemma 6.1.** Lemma 6.1 of ptolemy (c(Z) M^7 ordered tuples in Z) gives the same conclusion. Both
  constants depend on Z only. No fixed-S Subspace exceptional set enters.

## 5. Point (4): the final contradiction and uniformity

**Proved.** On the chosen tuple, H_Sigma <= 28(a/4 + lambda_*) W < cW/256, and (L_8) gives
H_Sigma >= c W/128 − C_0. So W < 256 C_0 / c = W_*.

The constants a, lambda_*, M_*, W_* depend only on (c, eta_0, C_0, d(Z)). They are fixed before X, the measure space,
the modulus, the parent integer or the grid is chosen. Nothing depends on the radius or on a template.

**Numerical margins quoted in VIII.3 (Computed, exact rationals).** These use Theorem D's c = 6/69, eta_0 = 1/381,
a = 2·10^{-5} and lambda_* = 10^{-16}.
* 28(a/4 + lambda_*) = 1.4000000028·10^{-4} < c/256 = 3.3967·10^{-4}.
* a < 1/381.
* 127 sqrt(8 lambda_*) = 3.592·10^{-6} < a/4 = 5·10^{-6}.
* The room for 2/M requires M > 2.67·10^{15}. So M_* >= 2.67·10^{15}, besides 4(28 + d).
* lambda = epsilon/2 = 5·10^{-17} <= lambda_* in VIII.2.

## 6. Corrections (wording only)

1. State the dependence: lambda_*, M_*, W_* depend only on c, eta_0, the O(1) constant C_0 of (L_8), and the degree d
   of a polynomial vanishing on the preimage of Z.
2. In (VIII.2), read lambda_+ (or assume lambda >= 0).
3. In hypothesis 1, say that the moduli point equals the class of (x_0, ..., x_7) in M_{0,8}, that is, up to one common
   Möbius map. This is what each application must verify.
4. The lemma uses the linear form c w of (L_8). It does not follow from (L_8^f) with f = o(w), nor from the
   canonical-contact or lattice restrictions (L_8^c)^cat or (L_8^pc).

## 7. Point (5): numerical tests on synthetic sign systems (Computed)

The table below lists the quantities that `check_extract.py` computes for each system:
* the exact Gram matrix, with lambda_eff = max_{x≠y} G_xy / 4;
* E G^2;
* the lemma's bound 127 sqrt(2/M + 8 lambda);
* means, over 3000 random tuples in X^8, of Σ|m_A| and of max_T |mass_T / w − 1|, using exact integer masses.

| system | M | atoms | max G (off-diag) | E G^2 | bound | mean Σ\|m_A\| | mean max dev |
|---|---|---|---|---|---|---|---|
| Sylvester Hadamard 2^11 | 2048 | 2048 | 0 | 4.9e-4 | 3.97 | 0.061 | 0.061 |
| Paley b_x(u) = χ(u+x), p = 2003 | 2003 | 2003 | −5e-4 | 5.0e-4 | 4.01 | 1.86 | 0.64 |
| uniform random | 2000 | 8192 | 0.055 | 6.2e-4 | 42.2 | 1.14 | 0.37 |
| uniform random, random masses | 2000 | 8192 | 0.058 | 6.6e-4 | 43.6 | 1.30 | 0.43 |
| biased, P(+1) = 0.65 | 2000 | 8192 | 0.146 (E G = 0.09) | 8.7e-3 | 68.9 | 3.57 | 3.13 |
| 6 clusters, 2% noise | 2000 | 8192 | 0.94 | 0.142 | 174 | 10.9 | 10.4 |

* **Hadamard.** m_A ∈ {0, 1}, and m_A = 1 exactly when the indices over A XOR to 0. A tuple is exactly balanced unless
  such a dependence occurs, which has probability about 127/2048 = 0.062. Observed: 6.0% unbalanced, 94% exactly
  balanced. This is the 1/M term of (VIII.2) at work.
* **Paley.** This is the quadratic-residue system behind VIII.2. `paley_scale.py` at p = 10007, 100003, 1000003 gives
  mean Σ|m_A| = 0.83, 0.26, 0.083. That is 0.65·127/sqrt(p) in each case, against the lemma's 127 sqrt(2/p). It also
  gives p·E m_A^2 ≈ 1.0 for |A| = 4, 6, 8, so the step E m_A^2 <= E G^2 is tight up to a factor of about 2 here.
* **Balance needs huge M.** The balance conclusion is correct, but it is reached only at very large M. Σ|m_A| <= a
  needs M of order (508/a)^2 (about 10^8 at a = 0.05). For random systems it also needs the number of atoms to be of
  order (127/a)^2, because each m_A has size about n^{-1/2} there. This is consistent with the lemma, where the atom
  count enters through lambda >= max G/4.
* **Negative controls.** These are the biased and clustered systems. The one-sided hypothesis G <= 4 lambda fails for
  small lambda, E G^2 is no longer O(1/M), and balance fails badly at the same M. For example, in the biased case
  Σ|m_A| ≈ 28 μ^2 + 70 μ^4 + ..., with μ = 0.3.
* **What forces balance.** The one-sided bound G <= 4 lambda, together with the positivity E G >= 0, is what forces
  E G^2 = O(1/M + lambda), and hence balance. M large alone does not, as the clusters show. Large negative G cannot
  occur on average, because E G >= 0.

## 8. What I did not do

* I did not re-referee the proof of (3.2) and Lemma 6.2 (ptolemy §3(iv), §6.1). I only checked that they use nothing
  beyond primitivity, the factorisation and the chain condition.
* I did not check the applications VIII.2-VIII.6. In particular, I did not verify that their columns satisfy
  hypothesis 1 beyond the common-Möbius remark in Section 4.
