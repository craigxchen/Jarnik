# Referee report, round 31, Part IV, block "balanced" (astra_lemma_IV/answer.md, section IV.3)

Scripts and outputs are in `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref4L_balanced/`. I wrote all of them and did not run Astra's scripts. The one exception: I compared my factorizations with the list hard-coded in `iv_shell_certificate.py`.
* `ident.py` / `ident.out`: rebuilds the cycle-complemented profile S_i from the description, then checks both identities (IV.2) symbolically over Q[P,q].
* `inst.py` / `inst.out`: exact checks of the numerical instance, done independently. It covers factorizations, primitivity, the ten divisibilities and cofactors, formula (IV.3), the inequality, and the Smith form of the valuation matrix (index one).
* `width.py` / `width.out`: the doubled-angle diameter, computed with 80-digit arithmetic.

## 0. Bottom line

| Item | Verdict |
|---|---|
| Profile: S_i = (non-cycle edges at i) ∪ (cycle edges avoiding i); \|S_i\|=5, \|S_i∩S_j\|=2, \|S_iΔS_j\|=6 | CORRECT (computed) |
| Identity 1: T_iR_j−T_jR_i=(R_j−R_i)g_ij | CORRECT (all ten vanish identically in Q[P,q]) |
| Identity 2: p_*T_i^2−n_i has no P^5, P^4q or P^3q^2 term | CORRECT. Its P^2q^3 coefficients are 15177600, 13251249, 10476400, 11379456, 10319616 |
| Integrality and primitivity (q ≡ 0 mod 2∏R_i, P prime) | CORRECT |
| V_i=T_iα+R_iβ, signed cofactors R_j−R_i, exact formula (IV.3) | CORRECT |
| O(P^3) bound, so every exponent >3/5 | CORRECT (needs q=P^{o(1)}, which the progression theorem gives) |
| DFI + Tao–Ziegler inside one subset, giving o(n_i^{3/5}) | CORRECT. Both cited theorems exist with the hypotheses used (§3) |
| Instance q=924, K=3999318414, P=15994547776559171401 | CORRECT. Every claimed property was re-checked exactly (§4) |
| Obstruction: no o(P^{5/2}) maximum defect for this pencil | CORRECT, with the scope stated by Astra. Two minor precisions in §5 |

No false statement and no gap that affects a conclusion.

## 1. Table and identities

The cycle is 12, 23, 34, 45, 51. Complementing its five cuts gives the marked sets:
* S_1={13,14,23,34,45}
* S_2={15,24,25,34,45}
* S_3={12,13,15,35,45}
* S_4={12,14,15,23,24}
* S_5={12,23,25,34,35}

Here J_e=P+qb_e with b=(145,79,−31,109,85,−47,113,19,−41,25) in lexicographic order, p_*=P+329q and T_i=P^2+d_iqP+f_iq^2.

**Identity 1.** All ten instances of T_iR_j−T_jR_i−(R_j−R_i)∏_{e∈S_i∩S_j}J_e expand to 0.

**Identity 2.** p_*T_i^2−n_i has monomials only in P^2q^3, Pq^4 and q^5. For i=2 the coefficients have denominator 4, because f_2=43223/2.

So (IV.2) holds exactly as stated, including the bound O(q^3P^2) for 1≤q≤P.

The 329 shift of p_* and the shift 0 of P are distinct from the ten b_e. So all twelve progression terms are distinct.

I did not re-derive the table from the cubic pencil on abscissas 1, 1/2, −1/9, 1/13, −1/7. The direct identity check makes that derivation irrelevant to correctness.

## 2. Integrality, primitivity, construction, defect formula

**Integrality.** If 2 | q then f_2q^2 ∈ Z, so T_i ∈ Z.

**Coprimality of T_i and R_i.** Take a prime ℓ | R_i.
* If ℓ is odd: ℓ | q and ℓ | q^2/2, so T_i ≡ P^2 mod ℓ.
* If ℓ = 2 (this occurs only for i=1,4,5): f_i is an integer and q is even, so again T_i ≡ P^2 mod ℓ.
* P is a large prime, so P^2 ≢ 0 mod ℓ, and gcd(T_i,R_i)=1.

**Primitivity.** (α,β) is a unimodular basis, so V_i=T_iα+R_iβ is primitive exactly when gcd(T_i,R_i)=1.

**Minors.** [V_i,V_j]=(T_iR_j−T_jR_i)[α,β]=(R_j−R_i)g_ij. When the J_e are pairwise coprime, g_ij=gcd(n_i,n_j). The five R_i are distinct, so the cofactors are fixed, nonzero and independent of P and q. The largest is |R_1−R_5|=1872.

**Defect formula (IV.3).** Expand |V_i|^2=T_i^2|α|^2+2T_iR_i(α·β)+R_i^2|β|^2. Then use |α|^2|β|^2−(α·β)^2=[α,β]^2=1, which gives |β|^2=(s^2+1)/p_*. Formula (IV.3) is exact.

**O(P^3).** With the centred s (|s|≤p_*/2):
* |2sT_iR_i| ≤ R_iP^3(1+o(1));
* the third term of (IV.3) is ≤ R_i^2P/4(1+o(1));
* the first term is O(q^3P^2).

With n_i=P^5(1+O(q/P)), this gives |δ_i| ≤ (R_i+o(1))n_i^{3/5}. That is exponent 3/5 with a constant up to 4680, hence every θ>3/5 literally. ✓

**Directions.** |β| ≤ |α|/2·(1+o(1)) and T_i ~ P^2. The angular estimate uses only sines of pair angles, sin = (R_j−R_i)g_ij/(|V_i||V_j|). No common phase is set to zero. ✓

## 3. The analytic input (DFI and Tao–Ziegler)

**Duke–Friedlander–Iwaniec**, Ann. of Math. 141 (1995) 423–441. The theorem is that for an irreducible quadratic of negative discriminant, the fractions ν/p over roots ν mod p, p≤x, are equidistributed in R/Z. X^2+1 has discriminant −4, so it is covered.

Fix ε<1/2. The two roots of a split prime are ν and p−ν, so at most one of them lies in (0,εp). Hence the set A_ε of primes having a root in (0,εp) has relative natural density exactly ε among all primes. Every member except 2 is ≡1 mod 4.

For any α with |α|^2=p_* and any completing β, s^2+1=p_*|β|^2≡0 mod p_*. So s mod p_* is one of the two roots, and the centred representative satisfies |s|≤εp_*. ✓

**Tao–Ziegler, "Narrow progressions in the primes"** (arXiv:1409.1327). I could not fetch arXiv directly. Search snippets and Shao's follow-up (arXiv:1509.04955) give the statement: every subset of the primes of positive relative (upper) density contains k-term progressions with common difference r ≤ (log N)^{L_k}. Shao makes the exponent explicit as L_k=(k−1)2^{k−2}.

This is stronger than the P^{o(1)} that Astra uses. The fixed modulus m=2∏R_i is handled correctly:
* take a progression a+jr in A_ε of length m·376+1;
* put q=mr and P = a+j_0r with j_0=47m;
* then every P+qb, for b∈{−47,…,329}, is a term.

Applying TZ once to the single set A_ε puts all twelve numbers (P, the J_e, p_*) in A_ε, so p_* has a small root. Every member of A_ε except 2 is split, so the J_e are split primes. **The combination is legitimate:** DFI supplies one positive-density subset of the primes, and TZ applies to any such subset.

**Final bound.** |δ_i|/P^3 ≤ 2εR_i(1+o(1)) + O(q^3/P) + O(1/P^2). The diagonal choice ε=1/j, with P_j large enough for each j, gives o(n_i^{3/5}) along one sequence. ✓

The width is two-sided ≍ N^{-3/10}. On the instance, diam{2 arg V_i}·N^{3/10} = 3744.000 = 2·1872 (`width.out`).

**Standing cautions.** The exceptional sets and thresholds are non-uniform only in the harmless sense Astra states ("progression length and prime threshold need not be uniform in j"). No uniformity in the radius is claimed. No Gaussian gcd or denominator cost is dropped: the denominators are 2 and the fixed ∏R_i, and they are absorbed into q.

## 4. The numerical instance (all exact, `inst.out`)

**p_*.** p_* = K^2+1 = P+329·924 holds exactly. p_* is prime and ≡1 mod 4. P itself is **not** prime, which the instance does not need. Coprimality of T_i and R_i is checked directly: all five gcds are 1.

**Blocks.** Factorizations from sympy agree exactly with the certificate's list. Every factor is <2^64, where sympy's isprime is deterministic.
* Prime blocks: J_12, J_13, J_24, J_25, J_35. That is five prime blocks, as claimed.
* J_14 = 1669·3595601·2665287953
* J_15 = 9510509·1681776209513
* J_23 = 41·390110921379493901
* J_34 = 709·22559305749730873
* J_45 = 10513·1521406618145077

All ten blocks are squarefree and pairwise coprime, and every prime factor is ≡1 mod 4.

**Index one.** The anchored valuation matrix, with one row per prime, has Smith form (1,1,1,1). Index one ✓.

**Stars.** The stars are (KT_i−R_i, T_i), i.e. α=(K,1), β=(−1,0), [α,β]=1, s=−K. All five are primitive. gcd(n_i,n_j)=g_ij, and the ten signed cofactors equal R_j−R_i exactly: −819, −1820, −1512, −1872, −1001, −693, −1053, 308, −52, −360.

**Defects.** (IV.3) reproduces δ_i exactly. All δ_i are positive, and 0<(1000|δ_i|)^5<n_i^3 holds at every mark. The margins are tight: |δ_i|/n_i^{3/5} ∈ [5.08·10^{-4}, 7.46·10^{-4}].

**Remark.** In this instance about 99.7% of each defect is the q^3P^2 term (q=924 is large relative to P^{1/2}/10^4). The s-term is small because K^2+1 is prime, giving a root of size √p_*. The instance therefore illustrates (IV.2)/(IV.3) and does not illustrate the DFI mechanism. It is a single parent, not evidence for the asymptotic o(·).

## 5. Obstruction: "this pencil cannot have maximum defect o(P^{5/2})"

I checked each step.

1. **Same row plane.** Proportional Plücker vectors, with ν=det(A,B)≠0, mean the 2×5 rows span the same Q-plane as the rows (T_i) and (R_i). So V'_i=AT_i+BR_i with A,B∈Q^2. ✓
2. **Bounded denominators.** The lattice L spanned by the columns (T_i,R_i) has index G = gcd of the minors (R_j−R_i)g_ij. A prime ℓ divides at most one block, because the blocks are pairwise coprime. Some pair (i,j) has that block outside S_i∩S_j, so v_ℓ(G) ≤ v_ℓ(∏|R_j−R_i|). Since M(L)⊂Z^2 and G·Z^2⊂L, the denominators of A and B divide G, which is bounded. ✓
3. **Vandermonde step.** Subtract (IV.2). This needs q^3P^2=o(P^{5/2}), i.e. q=o(P^{1/6}), which holds for q=P^{o(1)}. Rescale to x=(a−p_*)P^4, y=bP^2, z=c. The limit matrix has rows (1, 2R_i, R_i^2), a Vandermonde with doubled middle column, and it is invertible for three distinct R_i. So x, y, z = o(P^{5/2}), which gives a−p_*=o(P^{−3/2}) and b=o(P^{1/2}). ✓
4. **Conclusion.** a has denominator dividing G^2, so a=p_* for large P. Then the Gram identity ac−b^2=ν^2 gives c=(b^2+ν^2)/p_*=o(1). This contradicts c=|B|^2≥1/G^2, since B≠0 because ν≠0. ✓

**Precisions.**
* The argument only needs ν^2=o(P), not ν=P^{o(1)}.
* It is an obstruction for the fixed pencil (T_i,R_i) with fixed minor ratios. Astra says this explicitly ("different parents or different minor ratios can evade it").
* Matching sharpness at P^{5/2} would need primes p_*=K^2+1 inside the same narrow prime pattern. That is conjectural, and Astra does not claim it. Section IV.3 correctly states only the conditional improvement "root of size p_*^{1−γ} gives O(P^{3−γ}+P^2q^3)". That statement is also correct by (IV.3).

## 6. Novelty and significance

**Novelty.** I grepped the record for "norm defect", "e(theta)" and Duke–Friedlander. The only absolute-defect analysis it contains is III.11's n^{4/5} (ref3_model.md). Lowering to o(n^{3/5}) at balanced norms with fixed cofactors is new to the record. So are the DFI device for centring the basis and the half-power pencil obstruction.

The construction is a natural refinement of III.11: the same cycle complement, with a cubic pencil replacing the polynomial row and lattice change. Its ingredients (TZ, DFI) are standard.

**Significance.**
* It sharpens a barrier: no argument that tolerates absolute norm errors of size o(n_i^{3/5}) in the balanced ten-pair orientation can beat Liouville's 10/3 at five marks.
* It says nothing about e(0), genuine clusters, the closed-set version (this row is one fixed projective pencil, so it does not escape closed sets), or Question B's 15/4.
* Combined with the pencil obstruction and IV.4's 2/5 threshold, the balanced orientation leaves the window [2/5, 3/5) for defect exponents undecided.

For the goal: peripheral. It is a clean, correct negative result about a class of methods.

Sources:
* [Tao–Ziegler, arXiv:1409.1327](https://www.arxiv.org/abs/1409.1327)
* [Shao, arXiv:1509.04955](https://export.arxiv.org/abs/1509.04955)
* [DFI citation via arXiv:2112.07494](https://arxiv.org/pdf/2112.07494)
