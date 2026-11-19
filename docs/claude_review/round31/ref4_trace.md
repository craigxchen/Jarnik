# Referee report, round 31, Part IV, block "trace"

Scope: in `astra_free_IV/answer.md`, Part IV (from line 1011): the finite-characteristic "trace" seed, "why simply assigning integer primes to this seed fails" (the symmetry exclusion), and the cone remark. The surface inequality (IV.2) is refereed elsewhere (`ref4_bmy.md`). I use it here only for calibration.

All code is mine. It is in `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref4_trace/`:

- `gf2d.py`: GF(2^d) and GF(2^d)[t] from scratch. The primitive polynomial is found by search, and arithmetic uses log tables, Euclidean gcd and exact division.
- `trace_check.py`: the full construction, recomputed from the definitions. Output is in `out_trace.txt`. It was run for q = 8, 16, 32, 64, 128, 256, 512.
- `norm_forms.py`: every binary quadratic form over F_q[t] with Q(P_a) in F_q·n_a for all marks. Output is in `out_norm.txt`.
- `odd_char.py`: the analogue for F_p-hyperplanes of F_{p^d}, for p = 3, 5, 7. Output is in `out_odd.txt`.
- `symmetry_check.py`: integer prime assignments, translations, pair separation with unequal weights, and the group theory (symbolic, plus brute force in PGL_2(F_p)). Output is in `out_sym.txt`.

I read Astra's `iv_trace_check.py` only after writing mine. It checks a single thing: that det equals its leading coefficient times the product over the *predicted* common set. It checks neither gcds, values, separation, budget nor primitivity. My checks cover all of these, and they agree with it.

## Verdicts

| # | Claim | Verdict |
|---|---|---|
| a1 | For q = 2^d, d >= 3: L_a is monic of degree q/2 with t dividing L_a. f_a(x) = x^(q/2-1) at marked x. deg(f_a - f_b) = q/4 - 1 exactly, and f_a - f_b = nonzero constant × the product of *all* common marked primes. Pair separation is q/2 > (q-1)/2. W = (q-1) log q = m log(m+1). Degree slack is 1/4, so constant cofactors fit. | **CORRECT.** Proof checked line by line. Reproduced exactly for q = 8, 16, 32, 64, 128, 256, 512 (up to 130,305 pairs), with gcds computed by Euclid. |
| a2 | This is an "integral fake at growth scale" in the F_q(t) model. | **CORRECT WITH CORRECTIONS.** See §2: the literal 1/4 normalisation needs q >= 256; the cone needs a shear of degree H >= q/4; the constant field changes with m. |
| a3 | Why W = 1·m log(m+1) does not contradict W >= 3 m log m - O(m). | Explained in §3. The integer constant 3 comes from forced residue collisions at primes of norm below M. F_q[t] has no primes of norm below q = m+1. |
| a4 | Can a norm condition be added? (referee's question) | **Partly yes; this is new and not in Part IV.** The construction is *exactly* a genuine cluster of the **split** form (hyperbola UV = N): q-1 divisors of t^(q-1)-1 in a "short interval". The only forms Q with Q(P_a) in F_q^*·n_a for these columns are inseparable (genuine points collapse) or isotropic (split). An anisotropic, circle-type norm is **not** reached. §4. |
| a5 | An analogue in odd characteristic? (referee's question) | The direct analogue (F_p-hyperplanes, p-linearised subspace polynomials) exists, with constant cofactors. It misses the critical budget by exactly q(p-2)^2/(4p^2) - 1/4, which is also half the separation deficit. So p = 2 is singled out by (p-2)^2 = 0. I did not prove that no odd-characteristic construction exists. §5. |
| b1 | For every assignment of integer primes to the cuts, d_ab depends only on a+b. | **CORRECT** (one line; checked on 12 random assignments, q = 8 to 64). |
| b2 | An integer realisation at the normalised bound III.5 would give a faithful homomorphism (F_q,+) → PGL_2(Q). This is impossible, so the trace table is excluded for large q. | **CORRECT WITH CORRECTIONS.** The hypotheses of III.7 are met, and checked in §6. "Large q" means q - 2 >= 16390, i.e. **d >= 15** with Astra's constants (d >= 14 with sharp counting). Pair separation is a genuine extra hypothesis: it fails for some prime assignments (§6). |
| b3 | An elementary abelian 2-subgroup of PGL_2(C) has order <= 4. | **CORRECT** (centraliser proof checked symbolically; brute force in PGL_2(F_p) for p <= 13 finds no order 8). |
| b4 | Remark: the rigidity survives any fixed multiplier C_1 under strict pair separation. | **CORRECT.** The margin 0.0664 W - log 2 - 2 log(4C_1) > 0 holds once W is large, and W >= log m. |

## 1. The construction (a1): proof check

Write H_a = {x in F_q : Tr(ax) = 0}, a hyperplane of the F_2-space F_q. The trace polynomial T(t) = sum_{j<d} t^(2^j) is F_2-linear, with kernel ker Tr. Hence L_a(t) = a^(-q/2) T(at) is monic of degree q/2, is 2-linearised, and has root set exactly H_a. So **n_a = L_a/t** (`trace_check.py` asserts this). The only monomials in L_a are t^(2^j), so
f_a = (t^(q/2) - L_a)/t = sum_{j<=d-2} a^(2^j - q/2) t^(2^j - 1) has degree q/4 - 1. Its leading coefficient is a^(-q/4), and a ↦ a^(-q/4) is a power of Frobenius composed with inversion, hence injective. At x in H_a \ 0 we have L_a(x) = 0, so f_a(x) = x^(q/2-1).

Two distinct hyperplanes meet in q/4 elements, of which q/4 - 1 are nonzero. f_a - f_b vanishes at those points and has degree exactly q/4 - 1, so it equals (a^(-q/4) - b^(-q/4)) · prod(t - x). Pair separation is |H_a Δ H_b| = q/2 > (q-1)/2 = L/2. The critical degree is (alpha_a + alpha_b)/2 - L/4 = q/4 - 3/4, so the slack is 1/4.

Two simplifications that Part IV does not state:

- f_a = t^(q/2-1) - n_a. So P_a = (f_a, 1) is the image of (n_a, 1) under the unimodular map [[-1, t^(q/2-1)],[0,1]], and **f_a - f_b = n_a - n_b**. Divisibility of the minor by g_ab is therefore automatic, since gcd(n_a, n_b) divides n_a - n_b. The whole content is the degree statement deg(n_a - n_b) = deg gcd(n_a, n_b) = q/4 - 1. In words: the subspace polynomials of hyperplanes agree in their top q/4 coefficients.
- When the rows have equal size and the cofactor is constant, slack = (|Δ| - L/2)/2. The degree budget and pair separation are then *the same condition*.

**Computed** (`out_trace.txt`). For q = 8, 16, 32, 64, 128, 256, 512 every assertion passes on all pairs: 21, 105, 465, 1953, 8001, 32385 and 130305 pairs. The checks are:
- the identity N = t^(q-1) + 1 for the parent;
- n_a | N;
- f_a(x) = x^(q/2-1) at marked x, and f_a(x) ≠ x^(q/2-1) at unmarked x;
- gcd(n_a, n_b) by Euclid, equal to the product over common marks;
- deg D = q/4 - 1, with a constant nonzero cofactor;
- slack exactly 1/4;
- separation exactly q/2.

## 2. In what sense this is an integral fake at growth scale (a2)

The dictionary for the F_q(t) model, with |f| = q^deg f, is:

- **Parent.** N = t^(q-1) - 1, the product of the q - 1 nonzero degree-one primes. |N| = q^(q-1), so W = (q-1) log q.
- **Profile.** Rows a in F_q^*, cuts x in F_q^*, and x is marked for a iff Tr(ax) = 0. This is the **Sylvester/Walsh Hadamard core** of order q with row 0 deleted and one physical prime per label. It is the pure Walsh profile of the record's `growth/walsh-extraction.md` §6 (model WF with F empty and b = 1).
- **Bound |t_ij| <= sqrt(d_ij)/(4N^(1/4)).** Here t_ab is a nonzero constant, so |t_ab| = 1. Also sqrt(d_ab)/N^(1/4) = q^(q/4 - (q-1)/4) = q^(1/4). In Part IV's form (IV.1) with C_1 = 1 the bound holds for every q. **With the literal Part III constant 1/4 it needs q^(1/4) >= 4, i.e. q >= 256 (d >= 8).** In degree language the factor 4 is invisible, and the slack is 1/4 degree.
- **Cone.** The unsheared columns have non-archimedean angular width |D|/(|P_a||P_b|) = q^(-(q/4-1)). This exceeds N^(-1/4) = q^(-(q-1)/4) by the factor q^(3/4). The common shear (f, 1) ↦ (f + t^H, 1) preserves all minors and gcds and gives width q^(q/4-1-2H). That is <= N^(-1/4) once H >= q/4. So Part IV's remark that "the cone is free without norms" holds here too (it is trivially correct over Z as well).
- **W = O(m log m).** Here m = q - 1 and W = m log(m+1), with constant 1.

**Caveat.** This is not a sequence of fakes over one function field. Each m = 2^d - 1 uses its own constant field F_{2^d}, and the prime norms q = m+1 grow with m. Over a fixed F_q(t), letting m grow forces higher-degree primes, and the construction does not apply.

## 3. Why constant 1 < 3 is not a contradiction (a3)

1. **Different classes.** The bound W >= 3 m log m - O(m) (`two_thirds/proof.md`, Thm 7.5) is proved for genuine clusters, i.e. lattice points with exact norms. Part IV's objects are integral fakes with no norm condition. Nothing in the record gives integral fakes a lower bound with constant 3.
2. **The mechanism is absent in F_q[t] at this scale.** That bound is driven by **forced collisions** (Cor. 3.6 and Lemma 3.3 of `two_thirds/proof.md`). Among M columns, those in a common residue class modulo a small prime ℓ, with ℓ+1 or ℓ-1 classes, give ℓ | c_ij. Summed with Mertens weights over ℓ < M, this produces the extra M log M terms.
   - Over Z there are about M/log M primes below M. The pigeonhole gains are unavoidable, and part of them also bind integral fakes, through primitive columns in P^1(F_ℓ).
   - In F_q[t] with m = q - 1, every prime has norm >= q > m, and P^1 of every residue field has >= q + 1 > m points. So **no collision is forced at any prime**, and the Mertens sum that produces the constant 3 is identically zero.
   - Also, m primes of norm about m (here q - 1 primes of norm q) are not available over Z, where π(m) ~ m/log m. The cheapest integer parent with m distinct primes already has W ~ m log m. The trace example sits at the absolute minimum of its own model.
3. **Characteristic.** The example has δ = 0, all E_ij = 0, and m = q - 1 → ∞. So the conclusion "m <= 24" of (IV.2) is **false in characteristic 2**. This is consistent, because Bogomolov–Miyaoka–Yau fails in positive characteristic and Langer's theorem is a characteristic-zero statement. The example is a sharp witness that (IV.2) is a characteristic-zero phenomenon. Part IV says as much ("the characteristics ... are essential").

## 4. Adding a norm condition over F_q(t) (a4); new

The Gaussian form X^2 + Y^2 = (X + Y)^2 is degenerate in characteristic 2, so the record's circle model has no direct characteristic-2 version over F_q(t). I asked which binary quadratic forms Q = αX^2 + βXY + γY^2 over F_q[t] give Q(P_a) = c_a n_a with c_a in F_q for **all** marks.

- **By hand.** Q(f_a, 1) ≡ Q(t^(q/2-1), 1) mod n_a. Degree bookkeeping forces Q = λ(tX^2 + Y^2) + b(XY + t^(q/2-1) Y^2) with λ, b in F_q. Then c_a = λ a^(-q/2) + b, using the identities t f_a^2 + 1 = a^(-q/2) n_a and N/n_a = L_a + a^(-q/2). Both identities are asserted in `trace_check.py` for q <= 512.
- **By machine** (`norm_forms.py`). An exact linear solve over F_q, with deg α, β, γ <= B for B up to 3q/2 (q = 8, 16), gives a solution space of **dimension exactly 2**. Every separable member (b ≠ 0) has an F_q(t)-rational isotropic vector: either (1,0) when λ = 0, or the column P_a with a^(-q/2) = b/λ.

The consequences:

- **b = 0.** The inseparable form tX^2 + Y^2 is the norm from F_q(√t). Every prime ramifies (π = π̄), so the genuine points π-part × π̄-co-part all coincide. A norm condition of this kind is vacuous.
- **λ = 0.** Q = bY(X + t^(q/2-1)Y). This is the **split form**. The genuine points are Z_a = (n_a, N/n_a) on the "hyperbola" UV = N. I checked deg(n_a - n_b) = q/4 - 1 and deg(N/n_a - N/n_b) = q/4, against size q/2, for all pairs up to q = 512.
  - So the trace construction **is a genuine cluster of the split model**: q - 1 points with D = q/2 and c = 0.
  - Equivalently, the q - 1 divisors n_a of t^(q-1) - 1 lie in a function-field short interval.
  - This is the function-field counterpart, in characteristic 2 with growing q, of the record's "split calibration": divisors in short intervals, Erdős–Rosenfeld and Erdős #887 (`lit18/direct.md`, `growth/conditional.md` Cor. 5.3, `expert_note/outreach_drafts.md` §3).
  - The divisors sit at degree q/2 - 1 against L/2 = (q-1)/2. The integer analogue would be an off-centre direction of height about q, which is outside the window of the Erdős–Rosenfeld bound 1 + C^2.
- **An anisotropic separable conic** (the true circle analogue) is **not** obtained from these columns: no such form exists. Whether some other characteristic-2 construction reaches an anisotropic norm is open. For instance one could try the constant-field extension F_{q^2}(t), but there the degree-one primes are inert, so the construction's primes cannot be split.

`growth/ptolemy.md` already says that the Plücker/star system "is blind to the quadratic form: the split form (divisors in short intervals) ... give[s] the same system". The trace seed is a concrete instance: a split-form genuine cluster that is therefore automatically an integral fake.

## 5. Odd characteristic (a5)

The natural analogue uses q = p^d, marks H_a = ker Tr(a·) up to F_p^* (that is, (q-1)/(p-1) marks), and the p-linearised subspace polynomial L_{H_a}. Then f_a = (t^(q/p) - L_{H_a})/t has degree q/p^2 - 1, f_a(x) = x^(q/p-1) on H_a \ 0, and f_a - f_b is a nonzero constant times gcd(n_a, n_b). The cofactor is nonzero because a ↦ a^(-(p-1)p^(d-2)) is injective on F_q^*/F_p^*.

`odd_char.py` confirms all of this for q = 9, 27, 81, 25, 125, 49. But the rows contain only a fraction 1/p of the primes:

- the separation is |H_a Δ H_b| = 2q(p-1)/p^2 = q/2 - q(p-2)^2/(2p^2);
- the excess over the critical degree is q(p-2)^2/(4p^2) - 1/4.

Computed excesses: 1/2 at q = 27, 2 at q = 81, 11 at q = 125, 6 at q = 49. The deficit vanishes only when (p-2)^2 = 0. Taking complements does not help, because Δ is unchanged. A Paley-type analogue (marks b, marked x when x - b is a nonzero square, so n_b = (t-b)^((q-1)/2) - 1) has minors of degree about q/2, far above the budget q/4. In short, I know of no odd-characteristic construction. I have not proved that none exists. (IV.2) does not apply there either, since it is a characteristic-zero result.

## 6. The symmetry exclusion (b)

**b1.** n_a n_b / g_ab^2 = prod{p_x : Tr(ax) ≠ Tr(bx)} = prod{p_x : Tr((a+b)x) = 1}. This depends only on a+b, for any weights. Checked with actual split primes (`symmetry_check.py`).

**Hypotheses of III.7.** The theorem is accepted in `ref3_proj.md`. I check that each of its hypotheses holds for the restriction to S_u = F_q^* \ {u} and the relabelled realisation V'_a = V_{a+u}.

1. **Same parent and the same labelled d.** d'_ab = d_{a+u,b+u} = d_ab.
2. **Orientation.** At x with Tr(ux) = 1 the cut on S_u is complemented, and otherwise it is unchanged. This was checked for every u and x. So V' realises the reoriented exponent array a'_{a,x} = a_{a+u,x}. Its g' are the g of that array, and its cofactors obey the bound III.5 with the same d. III.7 explicitly allows separately chosen orientations with the same d_ij, and `ref3_proj.md` confirms this is legitimate: A, B, C are computed from the d_ij.
3. **Pair separation** w(a_i - a_j) > W/2 on S_u, with W = log N for the whole parent.
   - It holds with equal weights (q/2 > (q-1)/2), and on average for any weights, because sum_c sep_c = (q/2) W.
   - **It is a real hypothesis.** It fails for some honest assignments: min separation/W = 0.479 to 0.488 with split primes starting at 13 (q = 16, 64). Weighting a hyperplane gives 0.26 W.
   - It holds for primes in a narrow bin. Part IV states it as a hypothesis ("Suppose the resulting array meets pair separation"), so this is not an error.
4. **Distinct directions** follow from nonzero minors.
5. **Size.** |S_u| = q - 2 must be at least the III.7 threshold: **q - 2 >= 16390, i.e. d >= 15**, with Astra's τ = 1/32 constants (per `ref3_proj.md`), and about d >= 14 with sharp counting. Part IV says only "sufficiently large q". The number should be stated. For 3 <= d <= 14 the integer question for this table remains open in the no-norm model.

**Composition.** T_u T_v and T_{u+v} (with u ≠ v) agree on dir V_a for a ∉ {0, v, u+v}. That is q - 4 >= 4 distinct points, more than the 3 that determine a projectivity. T_u^2 fixes q - 2 points, so it is the identity. T_u ≠ id, because dir V_{a+u} ≠ dir V_a. So u ↦ T_u (with T_0 = 1) is an injective homomorphism (Z/2)^d → PGL_2(Q).

**b3.** In PGL_2 over a field of characteristic not 2, every involution is conjugate to z ↦ -z. Here is why its centraliser is as stated:
- if MD = λDM then the eigenvalues force λ = ±1;
- λ = 1 gives diagonal M and λ = -1 gives antidiagonal M (sympy);
- the diagonal elements of order dividing 2 are only z and -z, and they have index 2.
So |A| <= 4. A brute-force search of PGL_2(F_p) for p = 3, 5, 7, 11, 13 finds no three pairwise commuting involutions generating order 8. The proof is correct and gives the contradiction for d >= 3.

**Generalisation (referee remark, not in Part IV).** The same argument shows the following. For a pair-separated profile above the III.7 threshold that is realised at the normalised bound, the group Aut(d) of mark permutations preserving all labelled d_ij embeds in PGL_2(Q). With distinct primes, this is the group of row permutations preserving every unoriented cut. Hence Aut(d) is cyclic of order 1, 2, 3, 4 or 6, or dihedral of order 4, 6, 8 or 12. I quote this standard list of finite subgroups of PGL_2(Q) from memory (Beauville, "Finite subgroups of PGL_2(K)", 2010) and did not re-derive it.
- The Walsh case is special because translations preserve each cut up to complement.
- Paley translations instead permute the cuts, so they do not preserve d for distinct primes, and they are not excluded this way.

**b4 (fixed C_1).** Replacing 1/4 by C_1 multiplies U by 4C_1. The III.7 margin [(1-7τ)/8 - τ] W - log 2 - 2 log(4C_1) stays positive for W large. W tends to infinity with m because the exponent rows are distinct (N >= m). This is correct.

## 7. Novelty and significance

- **Novelty.**
  - The trace seed (subspace polynomials of hyperplanes) does not occur in the record: there are no hits for "linearized" or "subspace polynomial". Mathematically it is elementary and classical: the top-coefficient sparsity of linearised polynomials.
  - Its identification as a split-form genuine cluster (§4) is the referee's.
  - The Walsh profile itself is the record's pure Sylvester core. For **genuine** clusters it is already excluded by `growth/walsh-extraction.md` Theorem B and Cor. B.2–B.3 (any C, for M >= max(28, C^(2+o(1)))).
  - The symmetry exclusion is new only for the **no-norm** integral-fake model. There it is a five-line corollary of III.7, valid from about 3.3·10^4 marks.
- **Significance for the goal.** **Negligible.** Neither item moves the growth rate log R/loglog R, the constant 2/3, or the uniform bound. Nor does either touch Questions A or B or route 5.
  - The trace seed lives in characteristic 2, with the constant field growing with m. Its mechanism (sparsity of linearised polynomials, and no primes of norm below m) has no integer counterpart.
  - The exclusion kills one highly symmetric incidence table, already dead for genuine clusters.
  - Their real value is calibration. (i) The no-norm integral-fake conditions can be met at W = m log(m+1) in positive characteristic, even by a genuine split-model cluster. (ii) Consequently, (IV.2) and any proof of (P) must use characteristic zero, or arithmetic collision gains, in an essential way.
- **Standing cautions.** None is violated: no fixed-S Subspace input, no discarded gcd or denominator cost, no phase set to zero, and the III.7 threshold is absolute (not template-dependent).
