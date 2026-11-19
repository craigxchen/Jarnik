# Referee report, round 31, Part III, block "model" (astra_free_III §§III.6, III.11)

Scripts: `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref3_model/`
* `five_model.py`: an independent rebuild of the III.11 five-mark model. It solves for the polynomial row itself and checks everything with exact integers. It also writes out the (L_5) anchored chart.
  `five_model.py small` uses K = 84259905319712200, so every block is fully factored. `five_model.py crt` follows Astra's CRT recipe and gives a 143-digit K. Outputs are in `small.out` and `crt.out`.
* `shear_model.py`: the uncomplemented five- and six-mark models with the shear, and the m ≤ 6 separation inequality.
* `horizontal.py`: an independent implementation of the III.6 CRT construction on the Paley double core, for m = 19, 43, 67.

All runs used prlimit --as=3.5e9. I wrote all three scripts myself. I did not run Astra's scripts.

## 0. Bottom line

* **III.11 is CORRECT.** I rebuilt it from the description and checked two explicit instances exactly. One clarification is needed (C1 below).
* **III.6 is CORRECT.** Its assertion 4 is conditional on III.2, which belongs to another block.
* **Neither section supplies the "fake with integrality" the coordinator asked for**, meaning integrality plus all divisibilities with m ≫ W/log W.
  * III.11 has fixed size. Astra proves correctly that the pair-cut method stops at m = 6.
  * III.6 grows, but it keeps only half the log divisibility, and its normalised minors are integers only approximately.
* **Novelty is modest.** Up to PGL_2(Q), the five-mark configuration is a point of the record's C1 family (`growth/ptolemy.md` Thm C1; `threshold/L8.md` §4: u_i = (s − τ_i)/a_i), specialised to s = K².
  * The record already knew that (L_5) is false (Prop C; `growth/referee_ptolemy.md` §3.3).
  * It already knew that C1 is "not a circle configuration … a point of the form-blind relaxation".
  * It already had the gap-only angle model at 10/3 (`round31/ref_square.md` §6).
  * What is new is an integer realisation of the critical-arc archimedean geometry at the same time: star norms ≈ n_i ≈ √N and a cone of width ≍ N^{-3/10}. It comes from an index-D lattice change in GL_2(Q), with split-only blocks (the s = K² trick with inert s_i).
* **Significance: peripheral.** Nothing changes for the growth rate, the uniform bound, Question A, Question B or route 5.

## 1. III.11: reconstruction and exact verification

**Set-up (checked).**
* Take s = (3,7,11,19,23), u_i = s_i², J_ij = K² + u_i u_j and U_i = (u_i, u_i² − K²).
* Then [U_i,U_j] = (u_j − u_i)(K² + u_i u_j) as a polynomial identity (asserted in the script).
* gcd(U_i) = gcd(u_i, K²) = 1 when K ≡ 1 mod s_i.

**Coprimality argument (checked).** A common prime of J_ij and J_kl divides u_iu_j − u_ku_l. The ten products are distinct (checked). If p | K, then p ∤ J. If p = s_a, then a ∈ {i,j} gives J ≡ 1 mod p. If a ∉ {i,j}, then p would need −1 to be a square mod s_a, but every s_a ≡ 3 mod 4. So the restriction to inert s_i is essential, and Astra's argument is right. J is odd, and J = K² + (s_is_j)² with gcd(K, s_is_j) = 1, so every prime factor is ≡ 1 mod 4.

**Polynomial row.** I solved the five linear conditions f(−u_iu_j) + (u_i+u_j)h(−u_iu_j) = 0 on the cycle F = (01,12,23,34,40) for six unknowns.
* The kernel has dimension exactly 1.
* Its primitive generator is **identical** to Astra's (f,h), with ratio 1.
* Res(f,h) = 2^17·3^4·5²·7²·11²·19²·23^4·47·73·139·167·191·223·227·229·257·3691, which is nonzero.

**Instance 1 (K = 84259905319712200, 17 digits).**
* Choice of K: K ≡ 0 mod 2·5·13·17·29·31·37·47 and K ≡ 1 mod each s_i, then random search.
* Blocks: all ten J are fully factored. They are squarefree, pairwise coprime and coprime to every u_i − u_j. Every prime is ≡ 1 mod 4.
* Index one: the Smith form of the anchored 4×10 valuation matrix is (1,1,1,1). Only one prime per block was used for this.
* Profile after complementing F: each row has exactly 5 marked blocks, and each pair shares exactly 2. g_ij = D J_ij/(A_iA_j) holds exactly.
* Lattice basis:
  * lattice content of (f(K²),h(K²)) is 1;
  * (w,v) is a basis of L with |det| = D;
  * the second entry of the first row equals K^7(1+o(1)).
* For four consecutive j:
  * **every column is primitive**;
  * **every g_ij divides its minor, with cofactor exactly ±(u_j − u_i) ∈ {40, 72, …, 520}**;
  * **all five norms differ from n_i**, for every j tried.
* Norms: K²(|V_i|²/n_i − 1) ∈ {−233437.37, −180536.82, −93379.82, 122263.52, 204669.86}. So |V_i|²/n_i − 1 ≍ K^{-2}, and the same holds for |V_i|⁴/N − 1.
* Width: the diameter of {2 arg V_i} is 2.906·10^{-99}. Diam·N^{3/10} = 1040.0 = 2·max|u_i−u_j|, and the minimum pair separation ·N^{3/10} = 80.0 = 2·min|u_i−u_j|. So the width is two-sided ≍ N^{-3/10}. Also diam·N^{1/4} = 1.2·10^{-14}, which is ≤ 1/2.

**Instance 2 (Astra's recipe, K with 143 digits).**
* Choice of K: K ≡ 0 mod 2·∏(primes of differences except the s_i), K ≡ 1 mod s_i, and K² ≡ −u_iu_j mod q_e but not mod q_e², for ten split primes q_e larger than every product.
* Each q_e ∥ J_e. All checks above pass.
* The constants K²(|V_i|²/n_i − 1) are the same to 10 digits as in instance 1. The model is a fixed algebraic family, and the norm defect is c_i K^{-2} with c_i ≠ 0 fixed.
* Diam·N^{3/10} = 1040.0 again.

**Isotropy.** At the distinguished primes I counted (column, prime) incidences with p | |V_i|². There are **0 of 25** marked incidences in both instances. So the model is not merely "non-isotropic at one line" (Astra's phrasing). Generically it is non-isotropic at every marked prime tested.

**Projective identity.** All cross-ratios of the V_i equal those of the U_i, as exact rationals. This must hold, since V_i = M U_i / A_i with M ∈ GL_2(Q). So the balanced star matrix is the uncomplemented one after a fixed rational projectivity, and the cycle complement is purely a choice of integral lattice and orientation.

**Other III.11 statements.**
* "Why this construction stops": the claim is 2(m−2)/C(m,2) > 1/2 iff m ∈ {5,6} (m² − 9m + 16 < 0). It is correct, and checked for m = 5..11. The averaging remark is trivially correct, since each two-against-the-rest cut separates exactly 2(m−2) pairs.
* Six-mark example (s_5 = 31, shear by K³). Checked exactly at K = 118431602337651154630:
  * all 15 minors equal J_ij(u_j − u_i);
  * the columns are primitive;
  * norm relative error is 2·10^{-34};
  * log diam/log N = −0.2612. The exact asymptotic is −4/15 = −0.2667, and diam·K^8 → 1904 = 2·max|u_i − u_j|.
* The unbalanced five-mark shear model was checked the same way (log diam/log N → −3/10).

**Correction / clarification C1 (scope of "asymptotically exact star radii").** The norm statement is true only in the weak form |V_i|² = n_i(1 + c_iK^{-2}) = n_i(1 + O(N^{-1/10})) with c_i ≠ 0.
* The absolute defect is ≍ n_i^{4/5}.
* As a relative error, N^{-1/10} is far coarser than the arc's own angular scale N^{-3/10}.
* So the model leaves open any input that uses norms to relative precision finer than N^{-1/10}. An example is |V_i|² ≡ 0 mod a positive proportion of n_i, which is isotropy at some primes.

Astra does not claim more. But "retaining … asymptotically exact star radii" should be stated with this rate.

**What a genuine circle has that the model lacks.**
* *Exact individual norms*, |V_i|² = n_i. Equivalently, at every marked prime p the reduction of V_i lies on an isotropic line, and on the same isotropic line for all columns marked at p.
* I checked Astra's proof that exact norms plus all divisibilities plus primitivity give a genuine cluster:
  * Primitivity forbids both conjugate factors, so a_ip sits on one factor π_p.
  * If V_i ≡ 0 mod π and V_j ≡ 0 mod π̄ with p | g_ij, then w = V̄_iV_j is divisible by π̄ but not π, while w̄ is divisible by π. So w − w̄ ≢ 0 mod π, which contradicts p | Im w.
  * Hence V_i = unit·∏π_p^{a_ip}, and z_i = Ḡ V_i²/n_i are genuine points.
* The model lacks this at essentially every prime (0/25 isotropic incidences). As ptolemy.md already says of C1, X_i² + Y_i² "has no reason to factor along the blocks."

**Does the model show that Question-A-type inputs cannot give more than 10/3 at five points?** Yes, in this precise sense. Consider the inputs:
* a 2×5 integer matrix with primitive columns;
* every prescribed g_ij dividing its minor with bounded cofactor;
* index one and the chain condition;
* star norms equal to n_i up to relative N^{-1/10};
* a cone of width ≤ ½N^{-1/4}.

On the ten-pair profile these inputs admit the width ≍ N^{-3/10}. So no argument that uses only these inputs can prove a five-point lower bound on the width of the form N^{-3/10+ε}, i.e. beat exponent 10/3. But 10/3 is the proved Cilleruelo–Córdoba value. The model therefore shows only that these inputs give CC and nothing more at five points. That was the record's expectation (τ_e = 10/3, round 23; ref_square §6). It does not touch the open question of whether five points lie on arcs O(R^{2/5}) infinitely often.

## 2. Comparison with the expert note §2 ((L_k), chain condition)

* **As written, Astra's star matrix is not in (L_k) normal form.**
  * (L_k) is anchored: P_0 = (1,0), and blocks n_T with T ⊆ [k−1] divide the minors with i,j ≥ 1.
  * Astra's complemented orientation has the cuts 12, 23, 34 marked on 3-sets that contain mark 0. So its divisibility pattern is the symmetric star pattern, not the anchored one.
  * Its archimedean normalisation (|V_i|² ≈ √N, narrow cone) plays no role in (L_k). (L_k) is PGL_2(Q)-invariant data plus heights.
* **After anchoring it is an (L_5) datum.** I computed it explicitly (`five_model.py`, last block).
  * Send mark 0 to ∞: x_j = 1/(t_j − t_0) with t_i = u_i − K²/u_i, cleared to integers.
  * Put P_0 = (1,0), P_j = (x_j, 1). Then Y_j = 1, the columns are primitive, and every visible block divides as prescribed. The visible blocks are the 6 pairs J_ab ⊆ [4] and the 4 triples [4]∖{j} carrying J_0j.
  * The cofactors are **the same constants for both K** (after removing a content 9 that is prime to every block): t = 23166, 31941, 33264, 8775, 10098, 1323.
  * So H_Σ = O(1), with η = O(K^{-2}/log K) → 0.
  * **The chain condition holds trivially.** The blocks are pairwise coprime, so each prime divides one n_T, which is a single chain. This is the pc case, where the record proves that chains carry no information (expert note §2, Round-18 fact 3).
* **Consistency with the record.** Everything is consistent.
  * (L_5) is proved false (Prop C; referee_ptolemy §3.3), and K·β = 2^{k−3}(k−8)+k+2 = −5 < 0 at k = 5.
  * Astra's configuration is a point of exactly the C1 shape u_i = (s − τ_i)/a_i, with s = K², τ_i = u_i², a_i = u_i. It is a class-β rational curve: each of the ten contacts is linear in z = K².
  * A single such curve does not refute (L_k), but the record already has the Zariski-dense repair at k = 5.
  * Nothing here bears on k = 8, where (L_8) is expected and K·β = +10. The model is fully compatible with the record's statement that (L_k) can hold only for k ≥ 8.

## 3. III.6 (growing m, half divisibility)

**Profile (checked independently).** I built the double core myself: columns (a,±), entries χ(a−i) off-diagonal and ±1 on the diagonal. For m = 19, 43, 67 every row has m ones, of which (m+1)/2 are in the first family. Every pair shares (m−1)/2 ones, (m+1)/4 of them in the first family, and is separated by m+1 cuts.

**Proof (checked line by line).**
* CRT period: deleting ≥ 3 ones per row gives M_i ≤ G/P².
* Candidates: X* = H_0 + G b_i, scanned over 2m+2 consecutive candidates.
  * At a selected one, the mod-p² condition forbids exactly one class, because the slope 2X_0·(2M_i/p) is a unit.
  * At a selected zero, two classes are forbidden.
  * The total forbidden is ≤ (m+1)/2 + (m−1) < 2m+2, and each prime exceeds the window.
* b_i = 2m²i + i² gives distinct |b_j − b_i|, because the d·[2m², 2m²+2m) blocks are disjoint for d < m.
* Error terms: O(m^{3−K} + m/P²) for the normalised minors, O(m³/√P + m^{-K}) for the norms, and doubled width ≪ m³/(H_0√P). All of these are right, given primes in [P, P e^{m^{-K-1}}] ⊂ [m^A, 3m^A].
* log g_{ij,S} = ½ log g_ij + O((log m)²): correct. The row-dependent O(log m) variant is also correct.

**My computation (`horizontal.py`), using consecutive split primes rather than a narrow bin.**
* All m: exact prescribed Gaussian valuations on S, with a consistent factor X ≡ r_p. g_{ij,S} divides every minor. Full g_ij divides **no** minor (0 of 171, 903 and 2211 pairs).
* m = 67, P = 10^15: max |det/g − (b_j−b_i)| = 1.1·10^{-4}, max norm error 0.038, width·N^{1/4} = 0.037.
* m = 19 and 43 are too small for the asymptotics to bite: width·N^{1/4} = 2.0 and 1.2. This is expected. The proof requires P ≥ m^{4L+60}.
* min log g_S / log g rises slowly (0.0, 0.24, 0.30) towards 1/2, because of the O(log m) deleted columns.

**Assertion 4** uses III.2 assertion 1 (random assignment, probability O(m⁴2^{-cm}) of any quartet passing with U = P/4), together with (3m³+1)² < P/4. The deduction is correct. **III.2 itself is not verified by me**, since it is outside this block.

**Full-isotropy endpoint and "exact norms ⇒ automatic isotropy" are CORRECT.**
* Exact norms: as in §1.
* Bounded mean slack: P_i | V_i, the quotient norms are positive integers, and half have norm ≤ e^{2C}. There are finitely many Gaussian integers up to units, so a positive fraction share one quotient q, and the rotation q/q̄ is common.

**Significance of III.6: peripheral, close to trivial.**
* For columns (X,1) with X ≈ √n_i, any CRT condition must have modulus below the spacing scale G ≈ √(n_i/P). So "half the log divisibility" is exactly what pure counting gives for free.
* The minors are approximately integral only after dividing by g_ij. The absolute error is g_ij·m^{-L}, which is enormous.
* Astra is candid that the full-divisibility case "is left open". III.6 is therefore not an integral fake in the sense the coordinator requested.

## 4. Novelty and significance

* **III.11.**
  * The PGL_2 content is the record's C1/(L_5)-false phenomenon, re-derived with split-only blocks via s = K².
  * The index-D lattice change that turns it into a critical-arc integer star matrix with ≈ exact radii is new to the record, and it is a clean construction. It is a routine but correct addition, and it makes the record's informal statement precise: at five points, integrality plus divisibility plus approximate norms are exactly as strong as CC.
  * It is the integer analogue of ref_square §6. It stops at m = 6 by Astra's own (correct) proof.
  * Peripheral for growth.
* **III.6.** New as a written statement, but trivial in substance (CRT ceiling). Peripheral.
* **Standing cautions.**
  * No fixed-S exceptional sets are used.
  * The Gaussian gcd is handled honestly: the model has non-isotropic lines, and it says so.
  * No drifting phase is set to zero. The common rotation in the endpoint lemma is a genuine Gaussian quotient q/q̄.
  * No template-dependent threshold is claimed to be uniform. Assertion 4 of III.6 is per-assignment, and Astra says so.
* **Status of the goal.** No change to the growth rate (2/3 + ε)·log R/log log R, to the uniform bound, to Question A (k = 8, untouched), to Question B, or to route 5.
