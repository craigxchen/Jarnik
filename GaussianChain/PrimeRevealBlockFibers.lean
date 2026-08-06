import GaussianChain.PrimeRevealProjection

namespace GaussianChain.PrimeReveal

/-- The subfamily of `A` realizing a prescribed sign pattern on `T`. -/
def blockFiber {n : ℕ} (A : Finset (SignVector n)) (T : Finset (Fin n))
    (y : T → Bool) : Finset (SignVector n) :=
  A.filter fun x ↦ restrictSigns T x = y

@[simp]
theorem mem_blockFiber {n : ℕ} (A : Finset (SignVector n))
    (T : Finset (Fin n)) (y : T → Bool) (x : SignVector n) :
    x ∈ blockFiber A T y ↔ x ∈ A ∧ restrictSigns T x = y := by
  simp [blockFiber]

/-- A realized projected pattern has a nonempty block fiber. -/
theorem blockFiber_nonempty_of_mem_projection {n : ℕ} (A : Finset (SignVector n))
    (T : Finset (Fin n)) {y : T → Bool}
    (hy : y ∈ coordinateProjection A T) :
    (blockFiber A T y).Nonempty := by
  rcases (mem_coordinateProjection A T y).mp hy with ⟨x, hxA, hxy⟩
  exact ⟨x, (mem_blockFiber A T y x).mpr ⟨hxA, hxy⟩⟩

/-- Distinct block patterns give disjoint fibers. -/
theorem blockFiber_disjoint_of_ne {n : ℕ} (A : Finset (SignVector n))
    (T : Finset (Fin n)) {y z : T → Bool} (hyz : y ≠ z) :
    Disjoint (blockFiber A T y) (blockFiber A T z) := by
  rw [Finset.disjoint_left]
  intro x hxy hxz
  have hy : restrictSigns T x = y := ((mem_blockFiber A T y x).mp hxy).2
  have hz : restrictSigns T x = z := ((mem_blockFiber A T z x).mp hxz).2
  exact hyz (hy.symm.trans hz)

/-- The realized block fibers partition `A` exactly. -/
theorem card_eq_sum_blockFiber {n : ℕ} (A : Finset (SignVector n))
    (T : Finset (Fin n)) :
    A.card = ∑ y ∈ coordinateProjection A T, (blockFiber A T y).card := by
  classical
  rw [← Finset.sum_fiberwise_eq_card_filter]
  simp [coordinateProjection, blockFiber]

/-- Every realized block pattern has at least one completion, so the number
of realized patterns is at most the number of successful full assignments. -/
theorem projection_card_le_card {n : ℕ} (A : Finset (SignVector n))
    (T : Finset (Fin n)) :
    (coordinateProjection A T).card ≤ A.card := by
  classical
  exact Finset.card_image_le

end GaussianChain.PrimeReveal
