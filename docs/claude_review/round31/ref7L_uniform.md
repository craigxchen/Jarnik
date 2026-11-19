# Referee report, round 31, block "uniform": Part VII, sections VII.1 and VII.2 (Astra)

Scope: the covering-family inequality (VII.2); the quantifier paragraph of VII.1; the passage to a fixed torus and (VII.3); the "one exceptional set" step; the theorem of VII.2 as a whole. I take the gcd estimate as given, since another referee checks it, but I re-read it for consistency. Scripts and outputs are in
`/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref7L_uniform/`:

- `cover_det.py`, `cover_det.out`: log-derivative determinant on Astra's Pell family;
- `mult_test.py`: the same determinant with multiplicities and a common place;
- `pell_sparse.py`: the Pell family pulled back by t = u - 1/u;
- `pell_terms.py`: term budgets of the Pell family at fixed n.

## Verdicts

| Claim | Verdict |
|---|---|
| (VII.1) height formula W = sum_v(max - min); W = 2D in the conjugate-norm model; h >= (1/2+delta)D iff h >= (1/4+delta/2)W | CORRECT |
| Quantifier paragraph (constant torus suffices; t -> t+s; t -> t^k) | CORRECT |
| (VII.2) covering-family inequality 4h_0+1 <= W_0 | CORRECT. The proof gives the stronger bound 4h_0+1 <= #{finite places where some ratio has a zero or pole}. |
| Passage to a fixed torus and (VII.3) W = d_F+d_G-o(d_F+d_G) | CORRECT, given the gcd estimate |
| "One exceptional set" | CORRECT WITH CORRECTIONS (wording only): the frozen family must use exact supports, not supports contained in the frozen ones |
| Theorem of VII.2 as a whole | CORRECT WITH CORRECTIONS (same wording point), conditional on the gcd estimate. The exponent 1/4 is sharp as b -> infinity (new observation, below). |

## (a) The covering-family inequality (VII.2)

Checked line by line.

**Independence of the vectors.** Restrict to a dominating component B, with B smooth on a dense open set. The evaluation map Phi: B x P^1 -> T is dominant. In characteristic 0, generic smoothness makes dPhi surjective at a general point. In the invariant frame, dPhi(d/dt) is the log-derivative vector omega_0. It is nonzero because the generic member is nonconstant; a component of constant maps has image e, since q(infinity) = e. One can therefore complete omega_0 by three coordinate derivations d_1, d_2, d_3 of K = C(B) to a basis at a general point. The determinant Delta(t) of the four vectors is then a nonzero element of K'(t). The derivations extend uniquely to the finite extension K' that splits all zeros and poles, because the characteristic is 0.

**Poles.** Write q_i = c_i prod_a (t-a)^{m_{i,a}}. Then
- d_k log q_i = d_k c_i / c_i - sum_a m_{i,a} d_k a / (t-a);
- d/dt log q_i = sum_a m_{i,a} / (t-a).

At each a, the polar parts of all four columns are therefore multiples of the single vector (m_{i,a})_i / (t-a). By multilinearity, any term of the determinant that uses two polar parts vanishes. So Delta has at most a simple pole at a, and no pole if m_{.,a} = 0. Such a place has max - min >= 1 in (VII.1). The infinite place is not charged, since all a_{i,infinity} = 0. Hence the total number of poles of Delta is at most #charged finite places <= W.

Multiplicities m_i > 1, places where several coordinates vanish to different orders, and the split model are all covered by this argument. The genuine model needs nothing extra.

**Zeros at infinity.** The identity q = e mod t^{-h_0} holds identically on the family, so the parameter derivatives of the first h_0 - 1 Laurent coefficients vanish. Hence d_k log q = O(t^{-h_0}) and d/dt log q = O(t^{-h_0-1}). This gives ord_infinity Delta >= 3h_0 + h_0 + 1. A nonzero rational function on P^1 has equally many zeros and poles, so 4h_0 + 1 <= #charged places <= W.

**Hypothesis "W <= W_0 for every member".** The generic W equals W at a general closed point, so the hypothesis passes to the generic member.

**Independent checks.**
- Pell five-mark family, n = 1 and 2, exact over Q(i), seeds (1,3,5,8), all choices of three seed derivations (`cover_det.out`). In every case Delta is nonzero, ord_infinity Delta = 4h+1 exactly (17 and 33), the denominator is squarefree and divides the product of the charged places, and the number of poles is W = 18 and 34. The numerator has degree 1. So (VII.2) is attained within one, and Theorem D is sharp on this family.
- Multiplicity test (`mult_test.py`): q_i = ((t-a_i)/(t-c))^{m_i} with m = (1,2,3,1). This family is dominant, with h = 1 and W = 10. At t = c all four coordinates have poles of different orders. Delta has exactly 5 simple poles, one at each charged place including c, and ord_infinity = 5 = 4h+1. This confirms the sharper "#places" form. W = 10 is far from tight here.

**Comparison with round 29 Theorem D.** Theorem D (round29/ff.md, line 203) gives r(L+1)+1 <= 2D. For r = 4 and h_0 = L+1 this is the same inequality, 4h_0 + 1 <= W. Theorem D's proof needed primitive cocharacters (rows in {0,±1}) for transversality. Astra's log-Wronskian proof removes that need and allows arbitrary multiplicities, the split model and any finite-type family. This is a genuine but modest generalisation of a record result, by a standard technique (the log-derivative Wronskian in the style of Brownawell–Masser and Voloch).

## (b) Quantifier paragraph of VII.1

CORRECT.
- Suppose a nonzero P(t,Y) over C(t) vanishes on every qualifying curve of large height. Real translates q(t+s) also qualify: conjugation on coefficients commutes with a real shift, and degree, primitivity and contact are unchanged. So P(t, q(t+s)) = 0. Fix t_0 with P(t_0, .) nonzero. Then P(t_0, .) vanishes on q(t_0 + R), an infinite subset of an irreducible curve, and hence on the whole curve. The constant set {P(t_0, .) = 0} is proper and does not depend on q.
- The substitution t -> t^k multiplies D, W and h by k and preserves primitivity (substitute in a Bezout identity), the norm and the image. So an offending image cannot be removed by discarding finitely many degrees.
- Excluding dominance degree by degree is indeed insufficient, because the union of countably many non-dominating families can be dense.

## (c) Passage to a fixed torus

Each step checked.

1. **Sequence by contradiction.** At step l, choose a curve with W >= l whose image is not contained in V(P_1 ... P_l). The curve is irreducible, so it is contained in none of the V(P_k), k <= l. This survives all later subsequences: the term with original index l still escapes V(P_k) for every k <= l. The final step uses a single term of large original index, so every use of "eventually" is along one fixed sequence and is legitimate.
2. **Bounded-length exponent vectors.** There are at most 8b distinct exponents. After fixing the combinatorial placement on a subsequence, every A_j, B_j is linear in the variables Z.
3. **Persistent relations.**
   - Use a diagonal subsequence over the countably many u in Z^r. The persistent set Lambda is a saturated subgroup, finitely generated, so it is orthogonal to n_l on a tail.
   - On the quotient Z^r/Lambda, every fixed nonzero vector pairs nonzero with n_l eventually.
   - **Cancellation question.** The same distinct exponent never sits in two coordinates within one A_j. So the difference of two monomials of one A_j is e_k - e_k', and n_l·(e_k - e_k') = n_k - n_k' is never 0, so the difference is not in Lambda. Hence no individual A_j or B_j loses or merges monomials in the quotient.
   - Astra clears denominators after the quotient, in the domain C[Z^{r'}]. The products are therefore nonzero, and no list can become zero.
   - A coordinate can become formally constant only if q_j = 1 identically. That curve already lies in {q_j = 1}, and the argument never uses non-constancy of individual coordinates, only W -> infinity. So this case is handled.
4. **Gcd removal in C[Z^{r'}] and bounded supports.** The Newton polytope of a factor is a Minkowski summand, so after removing a common monomial the supports are bounded. There are finitely many supports, including those of F_i - F_1, so they can be frozen on a subsequence.
5. **Smaller ring C[M].** Here M is generated by support differences within each list. Suppose f in C[M] is a unit upstairs. Then f = c Y^u with u in supp f, which lies in M, so f is a unit of C[M]. Hence "no common nonunit factor" is preserved, as Astra claims. The identity F_i G_i = F_1 G_1 survives the per-list monomial shifts.
6. **Size of n.** Each M-basis vector is a fixed integral combination of within-list support differences, and each such difference has |n·diff| <= d_F or d_G. So ||n''|| << d_F + d_G. The reverse bound is trivial. The relation is correct.
7. **Degree loss.**
   - On the tail every frozen support specialises injectively. The minimal and maximal weights of the F-list are therefore carried by surviving monomials. After the shift, the list has no common zero at 0 or infinity, and the raw degree is d_F.
   - The degree loss equals the sum over a in C^* of min_i ord_a F_i(t^n), which is the common vanishing on C^*. The same holds for G.
   - d_+ and d_- are the degrees of [F_i] and [G_i] = [1/q_i]. So W = d_F + d_G - (losses), and the gcd estimate gives (VII.3).
   - W is the intrinsic primitive height. No gcd or denominator cost is discarded: resonant common factors (example: Z_2 = Z_1^2 becomes persistent, giving the common factor Y - a) are removed formally before the count.
8. **Sanity check of the o(L).** Take n = (mk', (2m+1)k') with m -> infinity slowly. There are no persistent relations, and gcd(t^{n_1} - 1, t^{n_2} - 1) has degree k' = L/(2m+1) = o(L), not O(1). So the o(L) is the right form, and it cannot be improved to O(1) with moving data.

## (d) One exceptional set

1. **Choice of n_*.** By (VII.3) and h >= (1/4+eps)W,
   4h >= (1+4eps)(1-o(1))(d_F+d_G) > d_F+d_G
   holds on a tail. Any tail term works as n_*, and the frozen supports are injective at it.
2. **Correction (wording).** The family must be the locally closed set of coefficient lists with exactly the frozen supports: every coefficient on the support is nonzero, and the coefficients of F_i and F_1 differ exactly on S_{i1}. Otherwise contact fails. Example: if the top coefficient of F_1 is allowed to vanish, then all F_i lose that monomial too (it lies outside S_{i1}), deg F_1 drops, deg(F_i - F_1) does not, and the contact falls below h_*. Astra's phrase "scaling a monomial by a nonzero scalar preserves them" indicates that she means exact supports. The exact-support family is Q-constructible, which is all that is needed.
3. **Bounds for every exact-support member.**
   - Raw degrees are <= d_F^* and d_G^*, so W <= d_F^* + d_G^*. Common factors only lower W.
   - G represents the reciprocals because F_i G_i = F_1 G_1 and all entries are nonzero.
   - The unique top monomial of F_1 has a nonzero coefficient, and deg(F_i - F_1) <= max weight on S_{i1}. So the contact is >= h_* for every member and every lambda.
   - Constant members are allowed and harmless.
4. **No component dominates.** By (VII.2), since 4h_* + 1 > d_F^* + d_G^*. There are finitely many components. The image of a Q-constructible set under a Q-morphism is Q-constructible, so the image closure Z^* is proper, Q-closed, and contained in some V(P_k).
5. **Containment of later curves.** At t = 1 the image is the union over all coefficient lists c in the family of Phi_c((C^*)^{r''}). Later terms have the same frozen supports in the same M-coordinates, with the reference exponent chosen as a fixed element of a fixed support. So each later term is a member, and its curve Phi_{c_l}(t^{n_l}) lies in Phi_{c_l}(torus), which is in Z^*. Other coefficients and other exponent vectors do not matter. This contradicts the escape property for terms of original index >= k.

The conclusion is uniform in degree, coefficients, exponents, the roots of the parent and the size of the cancelled common factors. W_0 is ineffective, which is legitimate for an existence statement.

### Standing cautions

- No fixed-S Subspace Theorem is used.
- Gcd and denominator costs are kept, through the intrinsic W and the formal removal of persistent factors.
- No drifting phase is set to zero: all coefficients, and the scaling lambda, are free in the family.
- The threshold is not template-dependent: W_0 comes from the contradiction argument and depends only on (b, eps).

### Attempts to break the theorem

None of the known dense families is a bounded-term family with h >= (1/4+eps)W.

- **Pell five-mark family in t.** h/W = 4n/(16n+2) < 1/4. The ratios are beta-bar G_j G_k / (beta G-bar_j G-bar_k), of degree 4n+1, so the term budget is b(n) <= 4n+2 and unbounded in n (`pell_terms.py`).
- **Pell family pulled back by t = u - 1/u.** This is a bounded-term family: each block times (u+1/u) is a Laurent polynomial with at most 8 terms (as sympy counts them), and (u+i)^2 = u·beta. It is dense, with W' = 2W = 32n+4 and h' = 4n at both u = infinity and u = 0 (`pell_sparse.py`, n = 1..4, exact). So h'/W' is about 1/8: no conflict with the one-place theorem. In the two-place version VII.7, 8h' = 32n < 32n+4: no conflict there either.
- **Cube family (round 29, Theorem C).** h = W/4 exactly, which is not >= (1/4+eps)W. It is also not dominating.
- **Single templates composed with t^k.** These are bounded-term, but each has one image curve, which is a proper closed set.
- **Characteristic-2 Moore family.** Not applicable: (VII.2) fails there (VI.6, not my block).

**New sharpness observation.** Fix n and compose the Pell family with t -> t^k.
- The result is bounded-term with b <= 4n+2, for every height W = (16n+2)k.
- Every height has the same dense union of images, so for every proper closed Z some member escapes Z.
- The ratio h/W = 4n/(16n+2) tends to 1/4.

So for every eps' > 0 there is a b such that the theorem fails with 1/4 replaced by 1/4 - eps'. The exponent 1/4 is optimal as b -> infinity. Whether eps = 0 is allowed remains open.

## (e) The theorem as it should be recorded

> **Theorem (VII.2; proved, conditional on the gcd estimate refereed separately; characteristic 0).** Fix b >= 1 and eps > 0. There exist W_0 = W_0(b, eps) (ineffective) and a proper Zariski-closed Z_{b,eps}, contained in T = G_m^4 and defined over Q, with the following property. Let q = (q_2, ..., q_5) be in (C(t)^*)^4 with q(infinity) = e, and suppose each q_j = A_j/B_j for Laurent polynomials A_j, B_j in C[t, t^{-1}] with at most b monomials each. The presentations need not be reduced, and the coefficients and exponents are arbitrary. Let W = sum_v (max_i v(q_i) - min_i v(q_i)), with q_1 = 1, be the intrinsic height, and h = min_j ord_infinity(q_j - 1). If W >= W_0 and h >= (1/4 + eps)W, then q(P^1) is contained in Z_{b,eps}. The same holds after a common affine change of the source variable. For genuine five-tuples of degree D (W = 2D), take eps = delta/2.
>
> **Remark.** The exponent 1/4 cannot be lowered: the Pell five-mark family at fixed n, composed with t -> t^k, is bounded-term (b <= 4n+2) and dense with h/W = 4n/(16n+2), which tends to 1/4.
>
> **Ingredient (VII.2, proved).** A finite-type characteristic-0 family of nonconstant maps P^1 -> G_m^4 with q = e mod t^{-h_0}, whose generic member has s finite places where some coordinate has a zero or pole, satisfies 4h_0 + 1 <= s <= W if it dominates.

## Novelty and significance (blunt)

**Novelty.**
- (VII.2) is round 29 Theorem D with a cleaner proof that handles multiplicities and the split model. It also has the sharper count by places, which Astra states only implicitly.
- The VII.2 sparse theorem is new to the record: a grep for "sparse", "lacunary" and "term budget" finds nothing of this kind. It is the first exceptional set in the record that is uniform over all degrees, for families of unbounded degree, in characteristic 0.
- Methodologically it is a combination of known pieces: Zannier/Bombieri–Masser–Zannier-style passage to a fixed torus for lacunary specialisations, an Amoroso–Leroux–Sombra-type gcd bound for F(t^n) (FoCM 15 (2015); I confirmed only the abstract-level statement), and Theorem D. The combination is neat, and the "one n_*" trick is the essential idea.

**Significance for the goal.** Low. Generic primitive tuples have no bounded-term presentation (Astra says so in "what is missing"). The integer analogue of "bounded number of monomials in t^n" has no arithmetic meaning for Gaussian integers on short arcs, and Astra concedes that nothing transfers. It does not decide A_5(delta) in the function field for any delta < 1/10. Its value is that it separates characteristic 0 from characteristic 2 on the same class, and that (with the remark above) 1/4, i.e. delta = 0, is the sharp threshold inside the sparse class. This is consistent with the Vojta prediction that every delta > 0 works.
