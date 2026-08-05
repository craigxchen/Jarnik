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

/-- Every realized block pattern has at least one completion. -/
theorem projection_card_le_card {n : ℕ} (A : Finset (SignVector n))
    (T : Finset (Fin n)) :
    (coordinateProjection A T).card ≤ A.card := by
  classical
  exact Finset.card_image_le

/-- There is a realized block pattern whose fiber has at least the average size. -/
theorem exists_large_blockFiber {n : ℕ} (A : Finset (SignVector n))
    (T : Finset (Fin n)) (hA : A.Nonempty) :
    ∃ y ∈ coordinateProjection A T,
      A.card ≤ (coordinateProjection A T).card * (blockFiber A T y).card := by
  classical
  have hproj : (coordinateProjection A T).Nonempty := by
    rcases hA with ⟨x, hxA⟩
    exact ⟨restrictSigns T x, (mem_coordinateProjection A T _).mpr ⟨x, hxA, rfl⟩⟩
  by_contra h
  push_neg at h
  have hlt : ∀ y ∈ coordinateProjection A T,
      (coordinateProjection A T).card * (blockFiber A T y).card < A.card := h
  have hsum := card_eq_sum_blockFiber A T
  have hcardpos : 0 < (coordinateProjection A T).card := Finset.card_pos.mpr hproj
  have havg : ∀ y ∈ coordinateProjection A T,
      (blockFiber A T y).card < A.card := by
    intro y hy
    have := hlt y hy
    omega
  have hsumlt : ∑ y ∈ coordinateProjection A T, (blockFiber A T y).card <
      (coordinateProjection A T).card * A.card := by
    exact Finset.sum_lt_sum_of_nonempty hproj fun y hy ↦ havg y hy
  rw [← hsum] at hsumlt
  have hge : A.card ≤ (coordinateProjection A T).card * A.card := by
    nlinarith
  omega

end GaussianChain.PrimeReveal
