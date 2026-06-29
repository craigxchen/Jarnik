import Mathlib

namespace GaussianChain
namespace ModPrimeCollision

/-- Over an integral domain, a quadratic equation of the form `y^2 = b` has at most two roots. -/
theorem sqFiber_card_le_two
    {R : Type*} [CommRing R] [NoZeroDivisors R] [Fintype R] [DecidableEq R] (b : R) :
    Fintype.card {y : R // y ^ 2 = b} ≤ 2 := by
  classical
  by_cases hnonempty : Nonempty {y : R // y ^ 2 = b}
  · let y₀ : {y : R // y ^ 2 = b} := Classical.choice hnonempty
    have hroot (y : {y : R // y ^ 2 = b}) : (y : R) = y₀ ∨ (y : R) = -y₀ := by
      exact sq_eq_sq_iff_eq_or_eq_neg.mp (by rw [y.property, y₀.property])
    let f : {y : R // y ^ 2 = b} → Bool := fun y => if (y : R) = y₀ then false else true
    have hinj : Function.Injective f := by
      intro y z hf
      by_cases hy : (y : R) = y₀
      · by_cases hz : (z : R) = y₀
        · exact Subtype.ext (hy.trans hz.symm)
        · simp [f, hy, hz] at hf
      · by_cases hz : (z : R) = y₀
        · simp [f, hy, hz] at hf
        · rcases hroot y with hyroot | hyroot
          · exact False.elim (hy hyroot)
          rcases hroot z with hzroot | hzroot
          · exact False.elim (hz hzroot)
          exact Subtype.ext (hyroot.trans hzroot.symm)
    calc
      Fintype.card {y : R // y ^ 2 = b} ≤ Fintype.card Bool :=
        Fintype.card_le_of_injective f hinj
      _ = 2 := Fintype.card_bool
  · have hempty : IsEmpty {y : R // y ^ 2 = b} := not_nonempty_iff.mp hnonempty
    rw [Fintype.card_eq_zero_iff.mpr hempty]
    norm_num

/-- Repackage solutions of `x^2 + y^2 = a` by first choosing `x`, then a square root for `y`. -/
noncomputable def zmodCircleSolutionsEquivSigma (p : ℕ) (a : ZMod p) :
    {q : ZMod p × ZMod p // q.1 ^ 2 + q.2 ^ 2 = a} ≃
      Sigma fun x : ZMod p => {y : ZMod p // y ^ 2 = a - x ^ 2} where
  toFun q := ⟨q.1.1, ⟨q.1.2, by
    calc
      q.1.2 ^ 2 = (q.1.1 ^ 2 + q.1.2 ^ 2) - q.1.1 ^ 2 := by ring
      _ = a - q.1.1 ^ 2 := by rw [q.2]⟩⟩
  invFun s := ⟨(s.1, s.2.1), by
    calc
      s.1 ^ 2 + s.2.1 ^ 2 = s.1 ^ 2 + (a - s.1 ^ 2) := by rw [s.2.2]
      _ = a := by ring⟩
  left_inv q := by
    ext <;> rfl
  right_inv s := by
    cases s with
    | mk x y =>
      ext <;> rfl

/-- For prime `p`, the finite-field circle `x^2 + y^2 = a` over `ZMod p`
has at most `2p` points. -/
theorem zmod_circle_card_le_two_mul (p : ℕ) [Fact p.Prime] (a : ZMod p) :
    Fintype.card {q : ZMod p × ZMod p // q.1 ^ 2 + q.2 ^ 2 = a} ≤ 2 * p := by
  classical
  rw [Fintype.card_congr (zmodCircleSolutionsEquivSigma p a), Fintype.card_sigma]
  calc
    (∑ x : ZMod p, Fintype.card {y : ZMod p // y ^ 2 = a - x ^ 2})
        ≤ ∑ _x : ZMod p, 2 := by
          exact Finset.sum_le_sum fun x _ => sqFiber_card_le_two (a - x ^ 2)
    _ = Fintype.card (ZMod p) * 2 := by simp
    _ = 2 * p := by
          rw [ZMod.card]
          ring

/-- Reduction of a Gaussian integer modulo `p`, as a pair of residues. -/
def gaussianResidue (p : ℕ) (z : GaussianInt) : ZMod p × ZMod p :=
  ((z.re : ZMod p), (z.im : ZMod p))

/-- Reducing the real and imaginary parts modulo `p` preserves the norm equation. -/
theorem gaussianResidue_norm (p : ℕ) (z : GaussianInt) :
    (gaussianResidue p z).1 ^ 2 + (gaussianResidue p z).2 ^ 2 =
      (z.norm : ZMod p) := by
  simp [gaussianResidue, Zsqrtd.norm]
  ring

/-- A family of Gaussian integers on `norm = N` maps into the corresponding finite-field circle. -/
noncomputable def gaussianCircleResidueMap
    {α : Type*} (p : ℕ) (N : ℤ) (z : α → GaussianInt)
    (hz : ∀ a, (z a).norm = N) :
    α → {q : ZMod p × ZMod p // q.1 ^ 2 + q.2 ^ 2 = (N : ZMod p)} :=
  fun a => ⟨gaussianResidue p (z a), by
    rw [gaussianResidue_norm, hz a]⟩

section Pigeonhole

variable {α β : Type*} [Fintype α] [Fintype β] [DecidableEq β]

/-- Ordered colliding pairs for a map `α → β`, including diagonal pairs. -/
def orderedCollisionCount (f : α → β) : ℕ :=
  Fintype.card {p : α × α // f p.1 = f p.2}

/-- Ordered colliding pairs are the disjoint union of squares of fibers. -/
noncomputable def orderedCollisionEquivSigma (f : α → β) :
    {p : α × α // f p.1 = f p.2} ≃
      Sigma fun b : β => {a : α // f a = b} × {a : α // f a = b} where
  toFun p := ⟨f p.1.1, ⟨⟨p.1.1, rfl⟩, ⟨p.1.2, p.2.symm⟩⟩⟩
  invFun s := ⟨(s.2.1.1, s.2.2.1), by rw [s.2.1.2, s.2.2.2]⟩
  left_inv p := by
    ext <;> rfl
  right_inv s := by
    cases s with
    | mk b pair =>
      cases pair with
      | mk x y =>
        cases x with
        | mk x hx =>
          cases y with
          | mk y hy =>
            subst b
            simp

/-- The ordered collision count is the sum of the squares of fiber cardinalities. -/
theorem orderedCollisionCount_eq_sum_fiber_sq (f : α → β) :
    orderedCollisionCount f =
      ∑ b : β, Fintype.card {a : α // f a = b} *
        Fintype.card {a : α // f a = b} := by
  classical
  rw [orderedCollisionCount, Fintype.card_congr (orderedCollisionEquivSigma f),
    Fintype.card_sigma]
  simp [Fintype.card_prod]

/-- The cardinality of the domain is the sum of the fiber cardinalities. -/
theorem card_eq_sum_fiber (f : α → β) :
    Fintype.card α = ∑ b : β, Fintype.card {a : α // f a = b} := by
  classical
  rw [← Fintype.card_congr (Equiv.sigmaFiberEquiv f), Fintype.card_sigma]

/-- Cauchy's inequality applied to the fiber sizes: `|α|^2 ≤ |β|` times the number
of ordered colliding pairs. -/
theorem card_sq_le_card_mul_orderedCollisionCount (f : α → β) :
    ((Fintype.card α : ℝ) ^ 2) ≤
      (Fintype.card β : ℝ) * (orderedCollisionCount f : ℝ) := by
  classical
  have hcs :=
    sq_sum_le_card_mul_sum_sq
      (s := (Finset.univ : Finset β))
      (f := fun b : β => (Fintype.card {a : α // f a = b} : ℝ))
  have hsum :
      (∑ b : β, (Fintype.card {a : α // f a = b} : ℝ)) =
        (Fintype.card α : ℝ) := by
    exact_mod_cast (card_eq_sum_fiber f).symm
  have hcoll_mul :
      (∑ b : β, (Fintype.card {a : α // f a = b} : ℝ) *
        (Fintype.card {a : α // f a = b} : ℝ)) =
        (orderedCollisionCount f : ℝ) := by
    exact_mod_cast (orderedCollisionCount_eq_sum_fiber_sq f).symm
  have hcoll :
      (∑ b : β, (Fintype.card {a : α // f a = b} : ℝ) ^ 2) =
        (orderedCollisionCount f : ℝ) := by
    simpa [pow_two] using hcoll_mul
  rw [hsum, hcoll] at hcs
  simpa using hcs

/-- Pigeonhole collision lower bound. If the map lands in at most `N` boxes, then its
ordered collision count is at least `|α|^2 / N`. -/
theorem orderedCollisionCount_lower_bound (f : α → β) {N : ℕ}
    (hβN : Fintype.card β ≤ N) (hN : 0 < N) :
    ((Fintype.card α : ℝ) ^ 2) / (N : ℝ) ≤ orderedCollisionCount f := by
  rw [div_le_iff₀ (by exact_mod_cast hN)]
  have hbase := card_sq_le_card_mul_orderedCollisionCount f
  have hboxes : (Fintype.card β : ℝ) ≤ (N : ℝ) := by exact_mod_cast hβN
  have hcoll_nonneg : 0 ≤ (orderedCollisionCount f : ℝ) := by positivity
  calc
    ((Fintype.card α : ℝ) ^ 2)
        ≤ (Fintype.card β : ℝ) * (orderedCollisionCount f : ℝ) := hbase
    _ ≤ (N : ℝ) * (orderedCollisionCount f : ℝ) := by
        exact mul_le_mul_of_nonneg_right hboxes hcoll_nonneg
    _ = (orderedCollisionCount f : ℝ) * (N : ℝ) := by ring

/-- The same pigeonhole lower bound, written in the common off-diagonal-pair normalization.
The right side is the real-valued half of the ordered off-diagonal collision count. -/
theorem offDiagonalCollisionLowerBound (f : α → β) {N : ℕ}
    (hβN : Fintype.card β ≤ N) (hN : 0 < N) :
    ((Fintype.card α : ℝ) ^ 2) / (2 * (N : ℝ)) -
        (Fintype.card α : ℝ) / 2 ≤
      ((orderedCollisionCount f : ℝ) - Fintype.card α) / 2 := by
  have h := orderedCollisionCount_lower_bound f hβN hN
  have hNpos : (0 : ℝ) < N := by exact_mod_cast hN
  have hNne : (N : ℝ) ≠ 0 := ne_of_gt hNpos
  have hhalf :
      ((Fintype.card α : ℝ) ^ 2) / (2 * (N : ℝ)) ≤
        (orderedCollisionCount f : ℝ) / 2 := by
    calc
      ((Fintype.card α : ℝ) ^ 2) / (2 * (N : ℝ))
          = (((Fintype.card α : ℝ) ^ 2) / (N : ℝ)) / 2 := by
              field_simp [hNne]
      _ ≤ (orderedCollisionCount f : ℝ) / 2 := by
          exact div_le_div_of_nonneg_right h (by norm_num)
  calc
    ((Fintype.card α : ℝ) ^ 2) / (2 * (N : ℝ)) -
        (Fintype.card α : ℝ) / 2
        ≤ (orderedCollisionCount f : ℝ) / 2 - (Fintype.card α : ℝ) / 2 := by
            exact sub_le_sub_right hhalf _
    _ = ((orderedCollisionCount f : ℝ) - Fintype.card α) / 2 := by ring

end Pigeonhole

section GaussianResiduePigeonhole

variable {α : Type*} [Fintype α]

/-- Collision lower bound for Gaussian integer points on a fixed norm circle after reduction
modulo a prime `p`. -/
theorem gaussianCircleResidue_orderedCollision_lower_bound
    (p : ℕ) [Fact p.Prime] (N : ℤ) (z : α → GaussianInt)
    (hz : ∀ a, (z a).norm = N) :
    ((Fintype.card α : ℝ) ^ 2) / (2 * (p : ℝ)) ≤
      orderedCollisionCount (gaussianCircleResidueMap p N z hz) := by
  classical
  have h := orderedCollisionCount_lower_bound
    (f := gaussianCircleResidueMap p N z hz)
    (N := 2 * p)
    (zmod_circle_card_le_two_mul p (N : ZMod p))
    (Nat.mul_pos (by norm_num) (Fact.out : Nat.Prime p).pos)
  simpa [Nat.cast_mul] using h

/-- Off-diagonal collision lower bound for Gaussian integer points on a fixed norm circle after
reduction modulo a prime `p`. -/
theorem gaussianCircleResidue_offDiagonalCollisionLowerBound
    (p : ℕ) [Fact p.Prime] (N : ℤ) (z : α → GaussianInt)
    (hz : ∀ a, (z a).norm = N) :
    ((Fintype.card α : ℝ) ^ 2) / (4 * (p : ℝ)) -
        (Fintype.card α : ℝ) / 2 ≤
      ((orderedCollisionCount (gaussianCircleResidueMap p N z hz) : ℝ) -
          Fintype.card α) / 2 := by
  classical
  have h := offDiagonalCollisionLowerBound
    (f := gaussianCircleResidueMap p N z hz)
    (N := 2 * p)
    (zmod_circle_card_le_two_mul p (N : ZMod p))
    (Nat.mul_pos (by norm_num) (Fact.out : Nat.Prime p).pos)
  convert h using 1
  · congr 1
    rw [Nat.cast_mul]
    ring_nf

end GaussianResiduePigeonhole

end ModPrimeCollision
end GaussianChain
