import GaussianChain.ArcCombinatorics
import GaussianChain.RamanaDeterminant

namespace GaussianChain
namespace SubcriticalBound

open RamanaDeterminant

/-- The first `M` entries of a sequence lie on the Gaussian-integer circle of norm `N`. -/
def OnCircleUpTo (M N : ℕ) (z : ℕ → GaussianInt) : Prop :=
  ∀ i, i < M → z i * star (z i) = ((N : ℤ) : GaussianInt)

/-- The first `M` entries of a sequence are pairwise distinct. -/
def InjectiveUpTo {α : Type*} (M : ℕ) (f : ℕ → α) : Prop :=
  ∀ ⦃i j⦄, i < M → j < M → f i = f j → i = j

/-- The squared Euclidean chord length between two Gaussian integers, viewed in `ℂ`. -/
noncomputable def gaussianSqDist (z w : GaussianInt) : ℝ :=
  Complex.normSq ((w : ℂ) - (z : ℂ))

/-- The real squared chord length is the natural absolute Gaussian norm of the difference. -/
theorem gaussianSqDist_eq_norm_natAbs (z w : GaussianInt) :
    gaussianSqDist z w = ((w - z).norm.natAbs : ℝ) := by
  unfold gaussianSqDist
  rw [← GaussianInt.toComplex_sub]
  simp

/-- Coordinate form of the squared chord length between two Gaussian integers. -/
theorem gaussianSqDist_eq_sq_add_sq (z w : GaussianInt) :
    gaussianSqDist z w = (w.re - z.re : ℝ) ^ 2 + (w.im - z.im : ℝ) ^ 2 := by
  unfold gaussianSqDist
  rw [← GaussianInt.toComplex_sub]
  simp [Complex.normSq_apply]
  ring

/-- A real squared-chord upper bound by a natural number gives the corresponding Gaussian norm
bound over `ℕ`. -/
theorem pair_norm_natAbs_le_of_gaussianSqDist_le_nat {z w : GaussianInt} {B : ℕ}
    (h : gaussianSqDist z w ≤ (B : ℝ)) :
    (w - z).norm.natAbs ≤ B := by
  rw [gaussianSqDist_eq_norm_natAbs] at h
  exact Nat.cast_le.mp h

/-- A real squared-chord upper bound gives a Gaussian norm bound by the natural floor. -/
theorem pair_norm_natAbs_le_natFloor_of_gaussianSqDist_le {z w : GaussianInt} {K : ℝ}
    (h : gaussianSqDist z w ≤ K) :
    (w - z).norm.natAbs ≤ Nat.floor K := by
  rw [gaussianSqDist_eq_norm_natAbs] at h
  exact Nat.le_floor h

/-- Real-valued subcritical chord consequence: if the Ramana power threshold is too small for
`B`, then an injective block of `2s+1` lattice points on the circle contains a pair whose
squared chord length is larger than `B`. -/
theorem exists_pair_sqDist_gt_of_pow_lt {s N B : ℕ}
    {z : Fin (2 * s + 1) → GaussianInt}
    (hcircle : ∀ i, z i * star (z i) = ((N : ℤ) : GaussianInt))
    (hN : 0 < N) (hz : Function.Injective z)
    (hsmall : B ^ (s * (2 * s + 1)) < N ^ (s * s)) :
    ∃ i j, j ∈ Finset.Ioi i ∧ (B : ℝ) < gaussianSqDist (z i) (z j) := by
  obtain ⟨i, j, hij, hgt⟩ := exists_pair_norm_gt_of_pow_lt hcircle hN hz hsmall
  refine ⟨i, j, hij, ?_⟩
  rw [gaussianSqDist_eq_norm_natAbs]
  exact_mod_cast hgt

/-- If all strict-pair squared chord lengths are at most the natural number `B`, and `B` is below
the Ramana subcritical threshold, then the block cannot consist of distinct points. -/
theorem not_injective_of_sqDist_le_of_pow_lt {s N B : ℕ}
    {z : Fin (2 * s + 1) → GaussianInt}
    (hcircle : ∀ i, z i * star (z i) = ((N : ℤ) : GaussianInt))
    (hN : 0 < N)
    (hdiam : ∀ i j, j ∈ Finset.Ioi i → gaussianSqDist (z i) (z j) ≤ (B : ℝ))
    (hsmall : B ^ (s * (2 * s + 1)) < N ^ (s * s)) :
    ¬ Function.Injective z := by
  refine not_injective_of_pair_norm_pow_lt hcircle hN ?_ hsmall
  intro i j hij
  exact pair_norm_natAbs_le_of_gaussianSqDist_le_nat (hdiam i j hij)

/-- A pairwise-diameter version of `not_injective_of_sqDist_le_of_pow_lt`, with no ordering
condition on the two indices. -/
theorem not_injective_of_pairwise_sqDist_le_of_pow_lt {s N B : ℕ}
    {z : Fin (2 * s + 1) → GaussianInt}
    (hcircle : ∀ i, z i * star (z i) = ((N : ℤ) : GaussianInt))
    (hN : 0 < N)
    (hdiam : ∀ i j, gaussianSqDist (z i) (z j) ≤ (B : ℝ))
    (hsmall : B ^ (s * (2 * s + 1)) < N ^ (s * s)) :
    ¬ Function.Injective z :=
  not_injective_of_sqDist_le_of_pow_lt hcircle hN (fun i j _ => hdiam i j) hsmall

/-- Real-diameter version using `Nat.floor`: if all strict-pair squared chord lengths are at most
`K`, and `⌊K⌋₊` is below the Ramana subcritical threshold, the block is not injective. -/
theorem not_injective_of_sqDist_le_of_natFloor_pow_lt {s N : ℕ} {K : ℝ}
    {z : Fin (2 * s + 1) → GaussianInt}
    (hcircle : ∀ i, z i * star (z i) = ((N : ℤ) : GaussianInt))
    (hN : 0 < N)
    (hdiam : ∀ i j, j ∈ Finset.Ioi i → gaussianSqDist (z i) (z j) ≤ K)
    (hsmall : (Nat.floor K) ^ (s * (2 * s + 1)) < N ^ (s * s)) :
    ¬ Function.Injective z := by
  refine not_injective_of_pair_norm_pow_lt hcircle hN ?_ hsmall
  intro i j hij
  exact pair_norm_natAbs_le_natFloor_of_gaussianSqDist_le (hdiam i j hij)

/-- Pairwise real-diameter version of `not_injective_of_sqDist_le_of_natFloor_pow_lt`. -/
theorem not_injective_of_pairwise_sqDist_le_of_natFloor_pow_lt {s N : ℕ} {K : ℝ}
    {z : Fin (2 * s + 1) → GaussianInt}
    (hcircle : ∀ i, z i * star (z i) = ((N : ℤ) : GaussianInt))
    (hN : 0 < N)
    (hdiam : ∀ i j, gaussianSqDist (z i) (z j) ≤ K)
    (hsmall : (Nat.floor K) ^ (s * (2 * s + 1)) < N ^ (s * s)) :
    ¬ Function.Injective z :=
  not_injective_of_sqDist_le_of_natFloor_pow_lt hcircle hN (fun i j _ => hdiam i j) hsmall

/-- If the squared-diameter threshold `K` is subcritical after taking natural floors, then an
injective block contains a pair whose squared chord length is larger than `K`. -/
theorem exists_pair_sqDist_gt_of_natFloor_pow_lt {s N : ℕ} {K : ℝ}
    {z : Fin (2 * s + 1) → GaussianInt}
    (hcircle : ∀ i, z i * star (z i) = ((N : ℤ) : GaussianInt))
    (hN : 0 < N) (hz : Function.Injective z)
    (hsmall : (Nat.floor K) ^ (s * (2 * s + 1)) < N ^ (s * s)) :
    ∃ i j, j ∈ Finset.Ioi i ∧ K < gaussianSqDist (z i) (z j) := by
  classical
  by_contra hnone
  have hdiam : ∀ i j, j ∈ Finset.Ioi i → gaussianSqDist (z i) (z j) ≤ K := by
    intro i j hij
    exact le_of_not_gt fun hgt => hnone ⟨i, j, hij, hgt⟩
  exact (not_injective_of_sqDist_le_of_natFloor_pow_lt hcircle hN hdiam hsmall) hz

/-- A globally injective natural-indexed sequence is injective on every consecutive finite
window. -/
theorem nat_window_injective_of_injective {α : Type*} {n : ℕ} {f : ℕ → α}
    (hf : Function.Injective f) (j : ℕ) :
    Function.Injective fun i : Fin n => f (j + (i : ℕ)) := by
  intro i k hik
  apply Fin.ext
  have hnat : j + (i : ℕ) = j + (k : ℕ) := hf hik
  omega

/-- A finite-window version of `nat_window_injective_of_injective`. -/
theorem nat_window_injective_of_injectiveUpTo {α : Type*} {M s : ℕ} {f : ℕ → α}
    (hf : InjectiveUpTo M f) {j : ℕ} (hj : j < M - 2 * s) :
    Function.Injective fun i : Fin (2 * s + 1) => f (j + (i : ℕ)) := by
  intro i k hik
  apply Fin.ext
  have hiM : j + (i : ℕ) < M := by
    have hi : (i : ℕ) < 2 * s + 1 := i.isLt
    omega
  have hkM : j + (k : ℕ) < M := by
    have hk : (k : ℕ) < 2 * s + 1 := k.isLt
    omega
  have hnat : j + (i : ℕ) = j + (k : ℕ) := hf hiM hkM hik
  omega

/-- Abstract subcritical window lower bound.

If every window whose parameter span is at most `B` has all squared chords at most `K`, then the
Ramana subcritical threshold forces every injective `2s+1`-point window to have span larger
than `B`. -/
theorem window_span_gt_of_subcritical
    {M s N : ℕ} {K B : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hN : 0 < N)
    (hsmall : (Nat.floor K) ^ (s * (2 * s + 1)) < N ^ (s * s))
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
    have hi : (i : ℕ) < 2 * s + 1 := i.isLt
    exact hcircle (j + (i : ℕ)) (by omega)
  have hinj_block : Function.Injective block := hinj j hj
  obtain ⟨i, k, _hik, hlarge⟩ :=
    exists_pair_sqDist_gt_of_natFloor_pow_lt hcircle_block hN hinj_block hsmall
  have hsmall_pair : gaussianSqDist (block i) (block k) ≤ K :=
    hdiam_of_span j hj hspan_le i k
  exact not_lt_of_ge hsmall_pair hlarge

/-- Sliding-window cardinality bound obtained from the subcritical Ramana obstruction, with the
arc/chord geometry kept as the abstract hypothesis `hdiam_of_span`. -/
theorem card_le_of_subcritical_windows
    {M s N : ℕ} {a L K B : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hk : 2 * s ≤ M)
    (hB : 0 < B)
    (hN : 0 < N)
    (hsmall : (Nat.floor K) ^ (s * (2 * s + 1)) < N ^ (s * s))
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
  exact le_of_lt (window_span_gt_of_subcritical hN hsmall hcircle hinj hdiam_of_span j hj)

/-- Variant of `card_le_of_subcritical_windows` using distinctness on the first `M` points. -/
theorem card_le_of_subcritical_windows_of_injectiveUpTo
    {M s N : ℕ} {a L K B : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hk : 2 * s ≤ M)
    (hB : 0 < B)
    (hN : 0 < N)
    (hsmall : (Nat.floor K) ^ (s * (2 * s + 1)) < N ^ (s * s))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hdiam_of_span : ∀ j, j < M - 2 * s → t (j + 2 * s) - t j ≤ B →
      ∀ i k : Fin (2 * s + 1),
        gaussianSqDist (z (j + (i : ℕ))) (z (j + (k : ℕ))) ≤ K)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    (M : ℝ) ≤ ((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) * L / B :=
  card_le_of_subcritical_windows hk hB hN hsmall hcircle
    (fun _j hj => nat_window_injective_of_injectiveUpTo hz hj) hdiam_of_span hmem

/-- Instantiate the abstract diameter hypothesis from an ordered parameter: if squared chord
length is bounded by squared parameter separation, then a window of parameter span at most `B`
has squared diameter at most `B^2`. -/
theorem sqDist_le_sq_span_of_param_sq_bound
    {M s : ℕ} {B : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hB : 0 ≤ B)
    (hmono : ∀ p q, p ≤ q → q < M → t p ≤ t q)
    (hparam : ∀ p q, p < M → q < M → gaussianSqDist (z p) (z q) ≤ (t q - t p) ^ 2)
    {j : ℕ} (hj : j < M - 2 * s)
    (hspan : t (j + 2 * s) - t j ≤ B) :
    ∀ i k : Fin (2 * s + 1),
      gaussianSqDist (z (j + (i : ℕ))) (z (j + (k : ℕ))) ≤ B ^ 2 := by
  intro i k
  let p := j + (i : ℕ)
  let q := j + (k : ℕ)
  have hwin_lt : j + 2 * s < M := by omega
  have hp_le_end : p ≤ j + 2 * s := by
    dsimp [p]
    have hi : (i : ℕ) < 2 * s + 1 := i.isLt
    omega
  have hq_le_end : q ≤ j + 2 * s := by
    dsimp [q]
    have hk : (k : ℕ) < 2 * s + 1 := k.isLt
    omega
  have hp_lt : p < M := lt_of_le_of_lt hp_le_end hwin_lt
  have hq_lt : q < M := lt_of_le_of_lt hq_le_end hwin_lt
  have hstart_le_p : t j ≤ t p := hmono j p (by dsimp [p]; omega) hp_lt
  have hstart_le_q : t j ≤ t q := hmono j q (by dsimp [q]; omega) hq_lt
  have hp_le_finish : t p ≤ t (j + 2 * s) := hmono p (j + 2 * s) hp_le_end hwin_lt
  have hq_le_finish : t q ≤ t (j + 2 * s) := hmono q (j + 2 * s) hq_le_end hwin_lt
  have habs : |t q - t p| ≤ B := by
    rcases le_total p q with hpq | hqp
    · have hpq_t : t p ≤ t q := hmono p q hpq hq_lt
      have hnonneg : 0 ≤ t q - t p := sub_nonneg.mpr hpq_t
      rw [abs_of_nonneg hnonneg]
      have hdiff_le_span : t q - t p ≤ t (j + 2 * s) - t j := by
        nlinarith
      exact hdiff_le_span.trans hspan
    · have hqp_t : t q ≤ t p := hmono q p hqp hp_lt
      have hnonpos : t q - t p ≤ 0 := sub_nonpos.mpr hqp_t
      rw [abs_of_nonpos hnonpos]
      have hdiff_le_span : -(t q - t p) ≤ t (j + 2 * s) - t j := by
        nlinarith
      exact hdiff_le_span.trans hspan
  calc
    gaussianSqDist (z p) (z q) ≤ (t q - t p) ^ 2 := hparam p q hp_lt hq_lt
    _ ≤ B ^ 2 := by
      exact (sq_le_sq).mpr (by simpa [abs_of_nonneg hB] using habs)

/-- Cardinality bound with a concrete ordered-parameter geometry hypothesis. This is the form
used when `t` is an arclength parameter: the chord distance is no larger than parameter
separation. -/
theorem card_le_of_param_subcritical_windows
    {M s N : ℕ} {a L B : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hk : 2 * s ≤ M)
    (hB : 0 < B)
    (hN : 0 < N)
    (hsmall : (Nat.floor (B ^ 2)) ^ (s * (2 * s + 1)) < N ^ (s * s))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hmono : ∀ p q, p ≤ q → q < M → t p ≤ t q)
    (hparam : ∀ p q, p < M → q < M → gaussianSqDist (z p) (z q) ≤ (t q - t p) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    (M : ℝ) ≤ ((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) * L / B := by
  refine card_le_of_subcritical_windows_of_injectiveUpTo (M := M) (s := s) (N := N)
    (a := a) (L := L) (K := B ^ 2) (B := B)
    hk hB hN hsmall hcircle hz ?_ hmem
  intro j hj hspan i k
  exact sqDist_le_sq_span_of_param_sq_bound hB.le hmono hparam hj hspan i k

end SubcriticalBound
end GaussianChain
