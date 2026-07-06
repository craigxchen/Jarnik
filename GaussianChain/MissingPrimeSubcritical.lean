import GaussianChain.SubcriticalBound
import GaussianChain.DeterminantDivisibility
import GaussianChain.LogProductBounds

namespace GaussianChain
namespace MissingPrimeSubcritical

open SubcriticalBound
open DeterminantDivisibility
open LogProductBounds

/-- Missing-prime strengthened version of the real-diameter Ramana obstruction. -/
theorem not_injective_of_sqDist_le_of_natFloor_pow_lt_missing_prime_product
    {s N : ℕ} {K : ℝ} {z : Fin (2 * s + 1) → GaussianInt} {P : Finset ℕ}
    (hprime : ∀ p ∈ P, Nat.Prime p)
    (hcircle : ∀ i, z i * star (z i) = ((N : ℤ) : GaussianInt))
    (hN : 0 < N)
    (hpN : ∀ p ∈ P, ¬ (p : ℤ) ∣ (N : ℤ))
    (hdiam : ∀ i j, j ∈ Finset.Ioi i → gaussianSqDist (z i) (z j) ≤ K)
    (hsmall : (Nat.floor K) ^ (s * (2 * s + 1)) <
      (∏ p ∈ P, p ^ (2 * pairCollisionCount p z)) * N ^ (s * s)) :
    ¬ Function.Injective z := by
  refine not_injective_of_pair_norm_pow_lt_missing_prime_product hprime hcircle hN hpN ?_ hsmall
  intro i j hij
  exact pair_norm_natAbs_le_natFloor_of_gaussianSqDist_le (hdiam i j hij)

/-- If the strengthened missing-prime threshold is subcritical after taking natural floors, then
an injective block contains a pair whose squared chord length is larger than `K`. -/
theorem exists_pair_sqDist_gt_of_natFloor_pow_lt_missing_prime_product
    {s N : ℕ} {K : ℝ} {z : Fin (2 * s + 1) → GaussianInt} {P : Finset ℕ}
    (hprime : ∀ p ∈ P, Nat.Prime p)
    (hcircle : ∀ i, z i * star (z i) = ((N : ℤ) : GaussianInt))
    (hN : 0 < N)
    (hpN : ∀ p ∈ P, ¬ (p : ℤ) ∣ (N : ℤ))
    (hz : Function.Injective z)
    (hsmall : (Nat.floor K) ^ (s * (2 * s + 1)) <
      (∏ p ∈ P, p ^ (2 * pairCollisionCount p z)) * N ^ (s * s)) :
    ∃ i j, j ∈ Finset.Ioi i ∧ K < gaussianSqDist (z i) (z j) := by
  classical
  by_contra hnone
  have hdiam : ∀ i j, j ∈ Finset.Ioi i → gaussianSqDist (z i) (z j) ≤ K := by
    intro i j hij
    exact le_of_not_gt fun hgt => hnone ⟨i, j, hij, hgt⟩
  exact (not_injective_of_sqDist_le_of_natFloor_pow_lt_missing_prime_product
    hprime hcircle hN hpN hdiam hsmall) hz

/-- Convert the finite-field collision lower bound into a natural-number exponent lower bound. -/
theorem exponent_le_pairCollisionCount_of_real_le_lower_bound
    {n p : ℕ} [Fact p.Prime] {N : ℤ} {z : Fin n → GaussianInt} {E : ℕ}
    (hz : ∀ i, (z i).norm = N)
    (hE : (E : ℝ) ≤ ((n : ℝ) ^ 2) / (4 * (p : ℝ)) - (n : ℝ) / 2) :
    E ≤ pairCollisionCount p z := by
  have hcoll := pairCollisionCount_lower_bound_on_circle p N z hz
  exact Nat.cast_le.mp (hE.trans hcoll)

/-- If every chosen exponent is below the corresponding collision count, then the explicit
prime-power product is below the collision product. -/
theorem prod_pow_two_mul_le_of_exponent_le
    {n : ℕ} {z : Fin n → GaussianInt} {P : Finset ℕ} {E : ℕ → ℕ}
    (hone : ∀ p ∈ P, 1 ≤ p)
    (hE : ∀ p ∈ P, E p ≤ pairCollisionCount p z) :
    (∏ p ∈ P, p ^ (2 * E p)) ≤
      ∏ p ∈ P, p ^ (2 * pairCollisionCount p z) := by
  refine Finset.prod_le_prod' ?_
  intro p hp
  exact Nat.pow_le_pow_right (Nat.lt_of_lt_of_le Nat.zero_lt_one (hone p hp))
    (Nat.mul_le_mul_left 2 (hE p hp))

/-- A simplified product threshold with explicit exponents implies the collision-product
threshold required by the missing-prime obstruction. -/
theorem natFloor_pow_lt_missing_prime_product_of_exponent_le
    {s N : ℕ} {K : ℝ} {z : Fin (2 * s + 1) → GaussianInt} {P : Finset ℕ}
    {E : ℕ → ℕ}
    (hone : ∀ p ∈ P, 1 ≤ p)
    (hE : ∀ p ∈ P, E p ≤ pairCollisionCount p z)
    (hsmallE : (Nat.floor K) ^ (s * (2 * s + 1)) <
      (∏ p ∈ P, p ^ (2 * E p)) * N ^ (s * s)) :
    (Nat.floor K) ^ (s * (2 * s + 1)) <
      (∏ p ∈ P, p ^ (2 * pairCollisionCount p z)) * N ^ (s * s) := by
  have hprod := prod_pow_two_mul_le_of_exponent_le hone hE
  exact hsmallE.trans_le (Nat.mul_le_mul_right _ hprod)

/-- Windowed form of `natFloor_pow_lt_missing_prime_product_of_exponent_le`. -/
theorem window_natFloor_pow_lt_missing_prime_product_of_exponent_le
    {M s N : ℕ} {K : ℝ} {z : ℕ → GaussianInt} {P : Finset ℕ} {E : ℕ → ℕ}
    (hone : ∀ p ∈ P, 1 ≤ p)
    (hE : ∀ j, j < M - 2 * s → ∀ p ∈ P,
      E p ≤ pairCollisionCount p (fun i : Fin (2 * s + 1) => z (j + (i : ℕ))))
    (hsmallE : (Nat.floor K) ^ (s * (2 * s + 1)) <
      (∏ p ∈ P, p ^ (2 * E p)) * N ^ (s * s)) :
    ∀ j, j < M - 2 * s →
      (Nat.floor K) ^ (s * (2 * s + 1)) <
        (∏ p ∈ P,
          p ^ (2 * pairCollisionCount p
            (fun i : Fin (2 * s + 1) => z (j + (i : ℕ))))) * N ^ (s * s) := by
  intro j hj
  exact natFloor_pow_lt_missing_prime_product_of_exponent_le hone (hE j hj) hsmallE

/-- The finite-field collision lower bound supplies windowwise exponent lower bounds for every
prime in a fixed finite set. -/
theorem window_exponent_le_pairCollisionCount_of_real_le_lower_bound
    {M s N : ℕ} {z : ℕ → GaussianInt} {P : Finset ℕ} {E : ℕ → ℕ}
    (hprime : ∀ p ∈ P, Nat.Prime p)
    (hcircle : OnCircleUpTo M N z)
    (hE : ∀ p ∈ P,
      (E p : ℝ) ≤ (((2 * s + 1 : ℕ) : ℝ) ^ 2) / (4 * (p : ℝ)) -
        (((2 * s + 1 : ℕ) : ℝ) / 2)) :
    ∀ j, j < M - 2 * s → ∀ p ∈ P,
      E p ≤ pairCollisionCount p (fun i : Fin (2 * s + 1) => z (j + (i : ℕ))) := by
  intro j hj p hp
  letI : Fact p.Prime := ⟨hprime p hp⟩
  let block : Fin (2 * s + 1) → GaussianInt := fun i => z (j + (i : ℕ))
  have hz : ∀ i, (block i).norm = (N : ℤ) := by
    intro i
    exact RamanaDeterminant.norm_eq_int_of_mul_star_eq
      (hcircle (j + (i : ℕ)) (by have hi := i.isLt; omega))
  exact exponent_le_pairCollisionCount_of_real_le_lower_bound (p := p) (N := (N : ℤ))
    (z := block) hz (hE p hp)

/-- For primes in the small-prime range `4p ≤ s`, the finite-field collision lower bound for a
`2s+1` block is at least `s^2 / (2p)`. -/
theorem half_s_sq_div_le_collision_real_lower {s p : ℕ}
    (hp : 0 < p) (hle : 4 * p ≤ s) :
    ((s : ℝ) ^ 2) / (2 * (p : ℝ)) ≤
      (((2 * s + 1 : ℕ) : ℝ) ^ 2) / (4 * (p : ℝ)) -
        (((2 * s + 1 : ℕ) : ℝ) / 2) := by
  have hpR : 0 < (p : ℝ) := by exact_mod_cast hp
  have hleR : (4 : ℝ) * (p : ℝ) ≤ (s : ℝ) := by exact_mod_cast hle
  have hnonneg : 0 ≤ (2 : ℝ) * (s : ℝ) + 1 := by positivity
  have hmul : (4 : ℝ) * (p : ℝ) * ((2 : ℝ) * (s : ℝ) + 1) ≤
      (s : ℝ) * ((2 : ℝ) * (s : ℝ) + 1) := by
    exact mul_le_mul_of_nonneg_right hleR hnonneg
  rw [div_le_iff₀ (by positivity)]
  rw [sub_mul]
  field_simp [hpR.ne']
  norm_num at *
  nlinarith [hmul]

/-- Natural floor form of `half_s_sq_div_le_collision_real_lower`. -/
theorem floor_half_s_sq_div_le_collision_real_lower {s p : ℕ}
    (hp : 0 < p) (hle : 4 * p ≤ s) :
    ((Nat.floor (((s : ℝ) ^ 2) / (2 * (p : ℝ))) : ℕ) : ℝ) ≤
      (((2 * s + 1 : ℕ) : ℝ) ^ 2) / (4 * (p : ℝ)) -
        (((2 * s + 1 : ℕ) : ℝ) / 2) := by
  have hx_nonneg : 0 ≤ ((s : ℝ) ^ 2) / (2 * (p : ℝ)) := by positivity
  exact (Nat.floor_le hx_nonneg).trans (half_s_sq_div_le_collision_real_lower hp hle)

/-- Small-prime floor exponents are valid collision-count lower bounds in every window. -/
theorem window_floor_half_s_sq_div_le_pairCollisionCount
    {M s N : ℕ} {z : ℕ → GaussianInt} {P : Finset ℕ}
    (hprime : ∀ p ∈ P, Nat.Prime p)
    (hsmallPrime : ∀ p ∈ P, 4 * p ≤ s)
    (hcircle : OnCircleUpTo M N z) :
    ∀ j, j < M - 2 * s → ∀ p ∈ P,
      Nat.floor (((s : ℝ) ^ 2) / (2 * (p : ℝ))) ≤
        pairCollisionCount p (fun i : Fin (2 * s + 1) => z (j + (i : ℕ))) := by
  refine window_exponent_le_pairCollisionCount_of_real_le_lower_bound hprime hcircle ?_
  intro p hp
  exact floor_half_s_sq_div_le_collision_real_lower (hprime p hp).pos (hsmallPrime p hp)

/-- Abstract missing-prime subcritical window lower bound.

The hypothesis `hsmall` is allowed to depend on the window, because the residue collision count
is a statistic of that window. -/
theorem window_span_gt_of_missing_prime_subcritical
    {M s N : ℕ} {K B : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ} {P : Finset ℕ}
    (hprime : ∀ p ∈ P, Nat.Prime p)
    (hN : 0 < N)
    (hpN : ∀ p ∈ P, ¬ (p : ℤ) ∣ (N : ℤ))
    (hsmall : ∀ j, j < M - 2 * s →
      (Nat.floor K) ^ (s * (2 * s + 1)) <
        (∏ p ∈ P,
          p ^ (2 * pairCollisionCount p
            (fun i : Fin (2 * s + 1) => z (j + (i : ℕ))))) * N ^ (s * s))
    (hcircle : OnCircleUpTo M N z)
    (hinj : ∀ j, j < M - 2 * s →
      Function.Injective fun i : Fin (2 * s + 1) => z (j + (i : ℕ)))
    (hdiam_of_span : ∀ j, j < M - 2 * s → t (j + 2 * s) - t j ≤ B →
      ∀ i k : Fin (2 * s + 1),
        gaussianSqDist (z (j + (i : ℕ))) (z (j + (k : ℕ))) ≤ K) :
    ∀ j, j < M - 2 * s → B < t (j + 2 * s) - t j := by
  intro j hj
  by_contra hnot
  have hspan_le : t (j + 2 * s) - t j ≤ B := le_of_not_gt hnot
  let block : Fin (2 * s + 1) → GaussianInt := fun i => z (j + (i : ℕ))
  have hcircle_block : ∀ i, block i * star (block i) = ((N : ℤ) : GaussianInt) := by
    intro i
    exact hcircle (j + (i : ℕ)) (by have hi := i.isLt; omega)
  have hinj_block : Function.Injective block := hinj j hj
  obtain ⟨i, k, _hik, hlarge⟩ :=
    exists_pair_sqDist_gt_of_natFloor_pow_lt_missing_prime_product
      (z := block) hprime hcircle_block hN hpN hinj_block (hsmall j hj)
  have hsmall_pair : gaussianSqDist (block i) (block k) ≤ K :=
    hdiam_of_span j hj hspan_le i k
  exact not_lt_of_ge hsmall_pair hlarge

/-- Sliding-window cardinality bound from the missing-prime strengthened Ramana obstruction. -/
theorem card_le_of_missing_prime_subcritical_windows
    {M s N : ℕ} {a L K B : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ} {P : Finset ℕ}
    (hk : 2 * s ≤ M)
    (hB : 0 < B)
    (hprime : ∀ p ∈ P, Nat.Prime p)
    (hN : 0 < N)
    (hpN : ∀ p ∈ P, ¬ (p : ℤ) ∣ (N : ℤ))
    (hsmall : ∀ j, j < M - 2 * s →
      (Nat.floor K) ^ (s * (2 * s + 1)) <
        (∏ p ∈ P,
          p ^ (2 * pairCollisionCount p
            (fun i : Fin (2 * s + 1) => z (j + (i : ℕ))))) * N ^ (s * s))
    (hcircle : OnCircleUpTo M N z)
    (hinj : ∀ j, j < M - 2 * s →
      Function.Injective fun i : Fin (2 * s + 1) => z (j + (i : ℕ)))
    (hdiam_of_span : ∀ j, j < M - 2 * s → t (j + 2 * s) - t j ≤ B →
      ∀ i k : Fin (2 * s + 1),
        gaussianSqDist (z (j + (i : ℕ))) (z (j + (k : ℕ))) ≤ K)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    (M : ℝ) ≤ ((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) * L / B := by
  refine ArcCombinatorics.card_le_of_window_span (M := M) (k := 2 * s) (a := a)
    (L := L) (B := B) t hk hB hmem ?_
  intro j hj
  exact le_of_lt (window_span_gt_of_missing_prime_subcritical hprime hN hpN hsmall
    hcircle hinj hdiam_of_span j hj)

/-- Variant of `card_le_of_missing_prime_subcritical_windows` using finite-family distinctness. -/
theorem card_le_of_missing_prime_subcritical_windows_of_injectiveUpTo
    {M s N : ℕ} {a L K B : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ} {P : Finset ℕ}
    (hk : 2 * s ≤ M)
    (hB : 0 < B)
    (hprime : ∀ p ∈ P, Nat.Prime p)
    (hN : 0 < N)
    (hpN : ∀ p ∈ P, ¬ (p : ℤ) ∣ (N : ℤ))
    (hsmall : ∀ j, j < M - 2 * s →
      (Nat.floor K) ^ (s * (2 * s + 1)) <
        (∏ p ∈ P,
          p ^ (2 * pairCollisionCount p
            (fun i : Fin (2 * s + 1) => z (j + (i : ℕ))))) * N ^ (s * s))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hdiam_of_span : ∀ j, j < M - 2 * s → t (j + 2 * s) - t j ≤ B →
      ∀ i k : Fin (2 * s + 1),
        gaussianSqDist (z (j + (i : ℕ))) (z (j + (k : ℕ))) ≤ K)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    (M : ℝ) ≤ ((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) * L / B :=
  card_le_of_missing_prime_subcritical_windows hk hB hprime hN hpN hsmall hcircle
    (fun _j hj => nat_window_injective_of_injectiveUpTo hz hj) hdiam_of_span hmem

/-- Cardinality bound with concrete ordered-parameter geometry and missing-prime determinant
input. This is the form used for arclength-ordered arc points. -/
theorem card_le_of_param_missing_prime_subcritical_windows
    {M s N : ℕ} {a L B : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ} {P : Finset ℕ}
    (hk : 2 * s ≤ M)
    (hB : 0 < B)
    (hprime : ∀ p ∈ P, Nat.Prime p)
    (hN : 0 < N)
    (hpN : ∀ p ∈ P, ¬ (p : ℤ) ∣ (N : ℤ))
    (hsmall : ∀ j, j < M - 2 * s →
      (Nat.floor (B ^ 2)) ^ (s * (2 * s + 1)) <
        (∏ p ∈ P,
          p ^ (2 * pairCollisionCount p
            (fun i : Fin (2 * s + 1) => z (j + (i : ℕ))))) * N ^ (s * s))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hmono : ∀ p q, p ≤ q → q < M → t p ≤ t q)
    (hparam : ∀ p q, p < M → q < M → gaussianSqDist (z p) (z q) ≤ (t q - t p) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    (M : ℝ) ≤ ((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) * L / B := by
  refine card_le_of_missing_prime_subcritical_windows_of_injectiveUpTo (M := M) (s := s)
    (N := N) (a := a) (L := L) (K := B ^ 2) (B := B) (P := P)
    hk hB hprime hN hpN hsmall hcircle hz ?_ hmem
  intro j hj hspan i k
  exact sqDist_le_sq_span_of_param_sq_bound hB.le hmono hparam hj hspan i k

/-- Ordered-parameter cardinality bound using explicit lower exponents for the missing-prime
collision product. -/
theorem card_le_of_param_missing_prime_subcritical_windows_of_exponent_bound
    {M s N : ℕ} {a L B : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ} {P : Finset ℕ}
    {E : ℕ → ℕ}
    (hk : 2 * s ≤ M)
    (hB : 0 < B)
    (hprime : ∀ p ∈ P, Nat.Prime p)
    (hN : 0 < N)
    (hpN : ∀ p ∈ P, ¬ (p : ℤ) ∣ (N : ℤ))
    (hE : ∀ j, j < M - 2 * s → ∀ p ∈ P,
      E p ≤ pairCollisionCount p (fun i : Fin (2 * s + 1) => z (j + (i : ℕ))))
    (hsmallE : (Nat.floor (B ^ 2)) ^ (s * (2 * s + 1)) <
      (∏ p ∈ P, p ^ (2 * E p)) * N ^ (s * s))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hmono : ∀ p q, p ≤ q → q < M → t p ≤ t q)
    (hparam : ∀ p q, p < M → q < M → gaussianSqDist (z p) (z q) ≤ (t q - t p) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    (M : ℝ) ≤ ((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) * L / B :=
  card_le_of_param_missing_prime_subcritical_windows hk hB hprime hN hpN
    (window_natFloor_pow_lt_missing_prime_product_of_exponent_le
      (fun p hp => (hprime p hp).one_le) hE hsmallE)
    hcircle hz hmono hparam hmem

/-- Ordered-parameter cardinality bound using the finite-field collision lower bound to certify
explicit missing-prime exponents. -/
theorem card_le_of_param_missing_prime_subcritical_windows_of_real_exponent_bound
    {M s N : ℕ} {a L B : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ} {P : Finset ℕ}
    {E : ℕ → ℕ}
    (hk : 2 * s ≤ M)
    (hB : 0 < B)
    (hprime : ∀ p ∈ P, Nat.Prime p)
    (hN : 0 < N)
    (hpN : ∀ p ∈ P, ¬ (p : ℤ) ∣ (N : ℤ))
    (hE : ∀ p ∈ P,
      (E p : ℝ) ≤ (((2 * s + 1 : ℕ) : ℝ) ^ 2) / (4 * (p : ℝ)) -
        (((2 * s + 1 : ℕ) : ℝ) / 2))
    (hsmallE : (Nat.floor (B ^ 2)) ^ (s * (2 * s + 1)) <
      (∏ p ∈ P, p ^ (2 * E p)) * N ^ (s * s))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hmono : ∀ p q, p ≤ q → q < M → t p ≤ t q)
    (hparam : ∀ p q, p < M → q < M → gaussianSqDist (z p) (z q) ≤ (t q - t p) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    (M : ℝ) ≤ ((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) * L / B :=
  card_le_of_param_missing_prime_subcritical_windows_of_exponent_bound hk hB hprime hN hpN
    (window_exponent_le_pairCollisionCount_of_real_le_lower_bound hprime hcircle hE)
    hsmallE hcircle hz hmono hparam hmem

/-- Ordered-parameter cardinality bound where the missing-prime determinant product is certified
by the explicit small-prime floor exponents `⌊s^2 / (2p)⌋`. -/
theorem card_le_of_param_missing_prime_subcritical_windows_of_small_prime_floor_bound
    {M s N : ℕ} {a L B : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ} {P : Finset ℕ}
    (hk : 2 * s ≤ M)
    (hB : 0 < B)
    (hprime : ∀ p ∈ P, Nat.Prime p)
    (hsmallPrime : ∀ p ∈ P, 4 * p ≤ s)
    (hN : 0 < N)
    (hpN : ∀ p ∈ P, ¬ (p : ℤ) ∣ (N : ℤ))
    (hsmallE : (Nat.floor (B ^ 2)) ^ (s * (2 * s + 1)) <
      (∏ p ∈ P, p ^ (2 * Nat.floor (((s : ℝ) ^ 2) / (2 * (p : ℝ))))) *
        N ^ (s * s))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hmono : ∀ p q, p ≤ q → q < M → t p ≤ t q)
    (hparam : ∀ p q, p < M → q < M → gaussianSqDist (z p) (z q) ≤ (t q - t p) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    (M : ℝ) ≤ ((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) * L / B :=
  card_le_of_param_missing_prime_subcritical_windows_of_exponent_bound hk hB hprime hN hpN
    (window_floor_half_s_sq_div_le_pairCollisionCount hprime hsmallPrime hcircle)
    hsmallE hcircle hz hmono hparam hmem

/-- Logarithmic version of
`card_le_of_param_missing_prime_subcritical_windows_of_small_prime_floor_bound`. -/
theorem card_le_of_param_missing_prime_subcritical_windows_of_small_prime_floor_log_bound
    {M s N : ℕ} {a L B : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ} {P : Finset ℕ}
    (hk : 2 * s ≤ M)
    (hB : 0 < B)
    (hfloor : 0 < Nat.floor (B ^ 2))
    (hprime : ∀ p ∈ P, Nat.Prime p)
    (hsmallPrime : ∀ p ∈ P, 4 * p ≤ s)
    (hN : 0 < N)
    (hpN : ∀ p ∈ P, ¬ (p : ℤ) ∣ (N : ℤ))
    (hlog :
      ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor (B ^ 2) : ℝ) <
        (∑ p ∈ P,
          ((2 * Nat.floor (((s : ℝ) ^ 2) / (2 * (p : ℝ))) : ℕ) : ℝ) *
            Real.log (p : ℝ)) +
          ((s * s : ℕ) : ℝ) * Real.log (N : ℝ))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hmono : ∀ p q, p ≤ q → q < M → t p ≤ t q)
    (hparam : ∀ p q, p < M → q < M → gaussianSqDist (z p) (z q) ≤ (t q - t p) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    (M : ℝ) ≤ ((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) * L / B := by
  refine card_le_of_param_missing_prime_subcritical_windows_of_small_prime_floor_bound
    hk hB hprime hsmallPrime hN hpN ?_ hcircle hz hmono hparam hmem
  exact nat_pow_lt_prod_pow_mul_pow_of_log_lt hfloor hN
    (fun p hp => (hprime p hp).pos) hlog

/-- Weighted-log version of the missing-prime cardinality bound. The logarithmic determinant
gain is supplied as `s^2 * ∑ log p / p`, with the floor loss bounded by `2 * ∑ log p`. -/
theorem card_le_of_param_missing_prime_subcritical_windows_of_weighted_log_bound
    {M s N : ℕ} {a L B : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ} {P : Finset ℕ}
    (hk : 2 * s ≤ M)
    (hB : 0 < B)
    (hfloor : 0 < Nat.floor (B ^ 2))
    (hprime : ∀ p ∈ P, Nat.Prime p)
    (hsmallPrime : ∀ p ∈ P, 4 * p ≤ s)
    (hN : 0 < N)
    (hpN : ∀ p ∈ P, ¬ (p : ℤ) ∣ (N : ℤ))
    (hweightedLog :
      ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor (B ^ 2) : ℝ) <
        ((s : ℝ) ^ 2 * (∑ p ∈ P, Real.log (p : ℝ) / (p : ℝ)) -
            2 * (∑ p ∈ P, Real.log (p : ℝ))) +
          ((s * s : ℕ) : ℝ) * Real.log (N : ℝ))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hmono : ∀ p q, p ≤ q → q < M → t p ≤ t q)
    (hparam : ∀ p q, p < M → q < M → gaussianSqDist (z p) (z q) ≤ (t q - t p) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    (M : ℝ) ≤ ((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) * L / B := by
  refine card_le_of_param_missing_prime_subcritical_windows_of_small_prime_floor_log_bound
    hk hB hfloor hprime hsmallPrime hN hpN ?_ hcircle hz hmono hparam hmem
  have hfloorLower :=
    sum_two_floor_half_s_sq_div_mul_log_lower
      (P := P) (s := s)
      (fun p hp => (hprime p hp).pos)
      (fun p hp => (hprime p hp).one_le)
  have hfloorWithNorm :
      ((s : ℝ) ^ 2 * (∑ p ∈ P, Real.log (p : ℝ) / (p : ℝ)) -
            2 * (∑ p ∈ P, Real.log (p : ℝ))) +
          ((s * s : ℕ) : ℝ) * Real.log (N : ℝ) ≤
        (∑ p ∈ P,
          ((2 * Nat.floor (((s : ℝ) ^ 2) / (2 * (p : ℝ))) : ℕ) : ℝ) *
            Real.log (p : ℝ)) +
          ((s * s : ℕ) : ℝ) * Real.log (N : ℝ) := by
    simpa [add_comm, add_left_comm, add_assoc] using
      add_le_add_right hfloorLower (((s * s : ℕ) : ℝ) * Real.log (N : ℝ))
  exact hweightedLog.trans_le hfloorWithNorm

end MissingPrimeSubcritical
end GaussianChain
