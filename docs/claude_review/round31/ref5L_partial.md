# Referee report, round 31, Part V, block "partial": the partial-isotropy constructions

Scope: Part V, section V.3 of `astra_lemma_V/answer.md` (lines 1225-1282) and the matching summary items 5-9 in `final5.txt`. The items are: the IV.2 pencil fails isotropy; the four-block star; the flip from eight endpoint incidences to ten in a mixed orientation; the five-block curve and its change of parameter; and the specialisation lemma that all of these use. I did not referee V.1, V.2 or V.4 except where they feed into this block.

Scripts are in `round31/ref5L_partial/`. All are my own except Astra's three certificates, which I reran unchanged.

| script | what it checks |
|---|---|
| `r_pencil.py` | IV.2 pencil: the symbolic reduction mod J_ij at all 20 incidences, plus an integer instance |
| `r_star4_sym.py` | four-star: s3=s1s2 *identically* on the (v,t,k) parametrisation, [G,A] for symbolic w, and the Jacobian at three points |
| `r_star4_int.py astra` | Astra's 191-digit instance rebuilt independently: dets, incidences, prime-level isotropy, norm ratios, defects, width constant, leading direction, and the two Gaussian flips |
| `r_star4_lemma.py` | four-star at a **new shape** w=(1,-9,-5,7,8), specialised exactly as the lemma prescribes (dilation D·L, T≡1 mod L), at two sizes |
| `r_five_sym.py` | five-block curve: every polynomial identity in both charts, with d re-derived independently |
| `r_five_int.py S` | five-block curve specialised as in the lemma, at m/D' ≈ 10^4, 10^6, 10^8 |
| `glib.py` | Gaussian gcd and division, determinants and chords |

Astra's certificates were all rerun and all pass: `v_star4_certificate.py`, `v_isotropic_star.py` and `v_five_test.py`.

---

## 1. "The IV.2 pencil fails even one whole-block isotropy condition for every admissible fixed shape": CORRECT

Put P = -q b_ij in V_i = (P U_i + B x_i, L U_i). This gives V_i ≡ q^2(-b_ij(h_i - b_ij) + b_0 x_i, l(h_i - b_ij)) mod J_ij, and so

|V_i|^2 ≡ q^4·β, where β = (-b_ij(h_i-b_ij)+b_0x_i)^2 + l^2(h_i-b_ij)^2.

Here β > 0 always: if h_i ≠ b_ij the second square is positive, and otherwise the first square equals (b_0x_i)^2 > 0. For a fixed shape and q = P^{o(1)}, the integer q^4β satisfies 0 < q^4β < J_ij once P is large. So J_ij never divides |V_i|^2, and the blocks are prime, so there is no partial isotropy at a prime either.

`r_pencil.py` checks the symbolic identity at all 20 incidences for the record shape (-27,-270,135,120,40). It also finds that an integer instance (P = 10^60+7, q divisible by the primes of l·b_0·∏x_i) has 0 of 20 isotropic incidences.

The extension to "integral parameters of height P^{o(1)}" is also correct, since then β = P^{o(1)}. Scope: this covers the IV.2 pencil only, not the IV.3 balanced pencil.

## 2. The specialisation lemma (polynomial columns to integer models): CORRECT

I checked each step of the lemma.

- **No common zero of Re F_i and Im F_i.** A common zero would be a root of all four J_ij (j ≠ i), and these are pairwise coprime.
- **Second dilation.** The second dilation by L, together with T ≡ 1 mod L, makes every block and every real leading coordinate ≡ 1 at the bad primes. This gives primitivity there, coprimality to the cofactors, odd blocks, and sums of two *coprime* squares, so every prime factor is ≡ 1 mod 4.
- **Cofactors.** After dilation by D they become D^6·d_ij, which is fixed.
- **Squarefree sieve.** Each block has at most 2 roots mod p^2 off the discriminant primes. p^2 | J(T) with T ≤ X forces p ≪ X. Primes above Y cost O(X/Y + π(X)) = o(X). This is elementary, and no prime-tuple input is used (unlike IV.1, which used Tao–Ziegler).
- **Non-isotropy persists.** A nonzero linear remainder r(T) satisfies gcd(J(T), r(T)) | Res, so whole-block divisibility fails for large T.

Two cosmetic points:

- (a) The listed bad set does not include the primes of Res(J_e, remainder of |F_i|^2). A bounded prime of that kind can divide a non-isotropic block, which gives *prime-level* isotropy on a non-isotropic block. This does not violate the whole-block definition, and adding those primes to the bad set removes it. In all my instances the gcd at every non-isotropic incidence was 1.
- (b) Astra's 191-digit instance is **not** an output of the lemma. Its blocks are divisible by 13, 17 and 29, and 13 divides both J_{25} and its own cofactor t_{25} (0-based edge (1,4)), so 13^2 | det(V_2,V_5). This does no harm to anything the certificate asserts, but the instance does not illustrate the "blocks coprime to cofactors" clause.

## 3. Four complete isotropic blocks through mark 1, |V_i|^2/n_i → 1, width ≍ N^{-3/10}, escaping any fixed closed set: CORRECT, with clarifications on scale

**Algebra (reproduced).** Expanding gives G = T^4 - i s_1T^3 - s_2T^2 + i s_3T + s_4, so Re G and Im G are as stated. With s_3 = s_1s_2 we get [G,A] = s_4 + (s_3/s_1 - s_2)T^2 = s_4; `r_star4_sym.py` verifies this with symbolic real w.

The coefficient difference is c_j(u_j+T^2) - c_i(u_i+T^2) = (1/u_i - 1/u_j)(T^2+u_iu_j), with c = 1 - 1/u. This gives the ten identities, and F_1 = G with |F_1|^2 = n_1 identically.

On central edges, both F_1 and F_j vanish at T = i w_j, and the coefficient conjugates vanish at -i w_j. On satellite edges, at T = iτ with τ = ±w_jw_k: G(iτ) = ∏(τ - w_m) and Ḡ(iτ) = ∏(τ + w_m) are real and nonzero by the |w_aw_b|-distinctness hypothesis (which includes w_1 = 1), while the correction terms are purely imaginary. So exactly four blocks are isotropic.

**Instances.**

- *Astra's shape (1,-9,-7,2,4), rebuilt with my code.*
  - 27 distinct split primes and a 191-digit parent.
  - All ten determinants are exactly D^6 s_4(1/u_i-1/u_j)·J_ij, and all columns are primitive.
  - Whole-block isotropy holds at exactly the 8 central incidences. Prime-level isotropy occurs at 24 (prime, mark) pairs: the 12 central primes times 2 endpoints, and nothing on satellite blocks.
  - |V_i|^2/n_i - 1 = (0, -0.021, -0.019, -0.0014, -0.0076).
- *New shape (1,-9,-5,7,8)* (s_1 = 1, s_4 = 2520; one of 34 integer solutions with entries in [-20,20]), specialised exactly by the lemma at two sizes. The 22 bad primes run up to 107, D·L has 39 digits, and the parents have 823 and 883 digits.
  - Isotropic incidences are exactly the 8 central ones, and the gcd at the other 12 is 1.
  - All dets are as claimed.
  - |V_i|^2/n_i ∈ [0.98699, 1] at the first size and 1 - 2·10^{-6} at the second; mark 1 is exact.
  - log10(Δ·N^{3/10}) = 232.36 at both sizes, against a largest cofactor of 233 digits. So Δ ≍ N^{-3/10} holds with a constant of the size of the fixed cofactor.
  - Squarefreeness was checked only by trial division up to 10^5. The blocks are too large to factor, and existence comes from the sieve in §2.

**Size of the defects (not stated by Astra).** In both families δ_j = |F_j|^2 - n_j has degree 7 in T, while P = J ≈ T^2. So |δ_j| ≍ P^{7/2} = n_j^{7/8} at the satellites and δ_1 = 0. Measured exponents in units of P were 3.98 at the smaller size and 3.94 at the larger one, still falling, because of the factor D in the leading coefficient. "|V_i|^2/n_i → 1" is therefore true only with relative error n^{-1/8}. That is far coarser than IV.2's n^{1/4} absolute defect. Isotropy was bought by giving up norm accuracy.

**The numerical instance is not at critical width in any visible sense.** For the 191-digit parent, Δ·N^{3/10} ≈ 6.6·10^{43} (log Δ/log N ≈ -0.07). The certificate does not claim otherwise, but item 9 should not be read as a width illustration before shearing.

**Escape argument and quantifiers: correct.**

- Using `r_star4_sym.py`, s_3 - s_1s_2 vanishes **identically** on the (v,t,k) parametrisation. (Astra's certificate checks it only at one point.)
- The 4×4 determinant [D, ∂_vD, ∂_tD, ∂_kD] is nonzero at Astra's point and at two further points, (5/7,3,-4/3) and (2,1/3,5). So the direction map is dominant onto P^3, and rational shapes give a Zariski-dense set of directions.
- In the actual instance, θ_j/D_j is constant across j to about 10^{-4} relative (≈ -3.224·10^{-14}). This confirms that the leading term of (y_j - 1) is a common scalar times T^{-6}(w_j^{-2} - 1), with the drifting phase of mark 1 divided out correctly.
- The expansion φ(y(T)) = κ^m T^{-6m} φ_m(D) + O(T^{-6m-1}) is valid because y is rational in T with y(∞) = 1.
- Quantifier: the shape, and with it the fixed cofactors, depend on Z. That is the right order for A_5, whose exceptional set is fixed in advance.
- The diagonal sequence with growing (subpower) cofactors avoids every closed set. The Zariski closure of a set of Q-points of T_5 is defined over Q, so the reduction to countably many hypersurfaces is sound.
- No uniform exceptional set from a fixed-S theorem is used.

**Summary of the four-block family.**

| property | value |
|---|---|
| orientation | endpoint |
| n_i | product of 4 blocks, about N^{2/5} |
| isotropic blocks / incidences | 4 blocks (the star at mark 1), exactly 8 incidences |
| escapes closed sets | yes |
| cofactors | fixed for each shape (233 digits in my instance) |
| defects | δ_1 = 0; satellites ≍ P^{7/2} = n^{7/8} |
| width | Δ ≍ |t|·N^{-3/10} |

Mark 1 is literally a Gaussian integer of norm n_1, and the satellites agree with it at the central primes.

## 4. From eight endpoint incidences to ten in a mixed orientation: CORRECT; ten is bookkeeping

I did the flip on edges 12 and 13 (1-based) with Gaussian gcds γ_e = gcd(V_1, J_e) in Z[i]. The new columns are V_1/(γ_12 γ_13), V_2γ̄_13/γ_12, V_3γ̄_12/γ_13, and V_4, V_5 multiplied by γ̄_12 γ̄_13. I found:

- all columns primitive;
- all new g_ij = gcd(n'_i, n'_j);
- all ten signed cofactors unchanged;
- all squared chords unchanged;
- every |V_i|^2/n_i unchanged;
- 10 isotropic incidences among the 22 marked ones.

Clarifications:

1. Because the incidences counted depend on the orientation, the count itself carries no information. On a flipped block, the two new marked columns V_k γ̄_e are isotropic *because* they were multiplied by γ̄_e. The intrinsic quantity is the number of isotropic blocks, still 4, and so is the transfer bound (V.2), H ≤ N∏_{e∉E} J_e = N^{8/5+o(1)}. Ten incidences here contain exactly the information of eight.
2. The flip changes the normalisation: n'_i has 2, 4, 4, 6 and 6 blocks, so n'_i ranges from N^{1/5} to N^{3/5}. This is a legitimate mixed orientation, since each pair still has six separating blocks and the Liouville bookkeeping is the same. But it is neither the endpoint normalisation nor √N.
3. Flipping one central block gives 9 incidences, and flipping three gives 11.
4. Exact counts in the endpoint orientation are always even, because isotropy at a marked incidence propagates to the whole marked side (V.1). So "k = 1" and "k = 5" are attainable only as "at least". Astra's wording, "sets of one, two or five required incidences", is accurate. The exact counts 2, 4 and 6 off closed sets were not constructed.

## 5. Five complete endpoint blocks with fixed cofactors on one rational curve, and the change of parameter: CORRECT; no escape, as Astra says

Verified with my own code (`r_five_sym.py`):

- [G,B] = (x-1)^6 exactly.
- For edge 23 (h_23 = z_2 z_3 = 1), c_2(i) = c_3(i). The reason is u_2u_3 = 1, so (u_2-u_1)(u_1-1/u_2) = (u_3-u_1)(u_1-1/u_3), and one real d kills both L_2(i) and L_3(i). Re-deriving d = -G(i)/(B(i)c_2(i)) reproduces -418284581/1009043070 exactly, and it is real.
- [L_j,L_k] = d u_1(1/u_j - 1/u_k)(x-1)^6(x^2+h_jk^2).
- In the x-chart the isotropic incidences are exactly 12, 13, 14, 15 and 23, at both endpoints of each.

The change of parameter is sound.

- Because [G,B](1) = 0, G(1) and B(1) are R-parallel, so every L_j(1) is a real multiple of G(1). The ratios are rational: 1, -4760083/4804967, 19246317/4804967, 46962141/161851520, 3511510513/102421665.
- F_j = T^4L_j(1+1/T)/L_j(1) is monic of degree 4, and its real and imaginary parts are coprime.
- [F_j,F_k] = (nonzero rational)·J_jk, where J_jk = ((T+1)^2+h^2T^2)/(1+h^2) = (T-a)^2+b^2 with a = -1/(1+h^2) and b = h/(1+h^2). These are rational, as the lemma needs.
- The ten J are distinct because the |h| are distinct.
- In the T-chart the isotropic incidences are again exactly the 10 incidences on 12, 13, 14, 15 and 23.
- |F_1|^2 = n_1 exactly. The other defects have degree 7, so they are ≍ P^{7/2}.

Integer specialisation by the lemma (`r_five_int.py`):

- D has 52 digits, there are 46 bad primes, the parents have 3281 to 3361 digits, and the cofactors have 959 digits.
- Exactly the 10 claimed incidences are isotropic, and the gcd at the others is 1. All dets are exact, and every block is ≡ 1 mod 4.
- |V_i|^2/n_i - 1 is about 10^{-4}, 10^{-6} and 10^{-8} at the three sizes.
- log10(Δ·N^{3/10}) = 959.198 at all three sizes. So Δ ≍ N^{-3/10} with a constant of the size of the cofactor, and the exponent -3/10 is visible only once log10 N ≫ 3000.
- Squarefreeness was checked only up to 10^5. Existence comes from §2.

**Summary of the five-block family.**

| property | value |
|---|---|
| orientation | endpoint |
| n_i | 4 blocks each |
| isotropic blocks / incidences | 5 blocks (star at 1 plus edge 23, i.e. a triangle 1-2-3 plus 14, 15), 10 incidences |
| escapes closed sets | **no** — one shape, so the ratios lie on one fixed rational curve, which can be put into Z |
| cofactors | fixed and huge |
| defects | δ_1 = 0; others ≍ n^{7/8} |

The edges contain a non-star triple, so V.4's rigidity says defects ≤ P^{2-ε} would force exact norms. The defects P^{7/2} are consistent with that. For Question B this family is pointwise only: A_5 permits excluding its curve. Astra's own "EVIDENCE" paragraph correctly states that there is no density of shapes.

## 6. Novelty and significance

**Novelty.** Isotropic five-mark models are new to the record. `ref4L_quarter.md` recorded 0 of 20 isotropic incidences for IV.2; nothing earlier constructs isotropic blocks at critical width.

**Significance for Question B: small.**

- Every constructed isotropic block either passes through one exact mark (the star), or the configuration lies on a single curve.
- The H-transfer (V.2) for 4 or 5 known blocks gives H ≤ N^{8/5} or N^{3/2}, which is useless next to the 15/4 threshold. By (V.2), nine blocks would give H ≤ N^{11/10} and exponent 11/3. That contradicts A_5(δ) only for δ > 3/10 - 3/11 = 3/110, not for every δ < 1/30.
- The families also give up norm accuracy (defects n^{7/8}). So the conjunction "isotropy at a non-star set of blocks plus accurate norms" remains untested, and V.4 shows that this conjunction becomes rigid at defect ≤ P^{2-ε}.

What the block establishes is that **star-shaped local isotropy costs nothing.** Isotropy might become a lever at the first non-star configuration — a triangle, or three disjoint-ish edges — that holds on a Zariski-dense set of shapes, and no such configuration is constructed.
