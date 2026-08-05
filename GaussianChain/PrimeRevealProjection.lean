import GaussianChain.PrimeRevealCube

namespace GaussianChain.PrimeReveal

/-- Restrict a full sign vector to a chosen coordinate set. -/
def restrictSigns {n : ℕ} (S : Finset (Fin n)) (x : SignVector n) : S → Bool :=
  fun i ↦ x i.1

/-- The coordinate projection of a finite sign-vector family onto `S`. -/
def coordinateProjection {n : ℕ} (A : Finset (SignVector n)) (S : Finset (Fin n)) :
    Finset (S → Bool) :=
  A.image (restrictSigns S)

@[simp]
theorem mem_coordinateProjection {n : ℕ} (A : Finset (SignVector n))
    (S : Finset (Fin n)) (y : S → Bool) :
    y ∈ coordinateProjection A S ↔ ∃ x ∈ A, restrictSigns S x = y := by
  simp [coordinateProjection]

/-- Projection cardinality is monotone under enlarging the original family. -/
theorem coordinateProjection_mono {n : ℕ} {A B : Finset (SignVector n)}
    (hAB : A ⊆ B) (S : Finset (Fin n)) :
    (coordinateProjection A S).card ≤ (coordinateProjection B S).card := by
  apply Finset.card_le_card
  intro y hy
  rcases (mem_coordinateProjection A S y).mp hy with ⟨x, hxA, rfl⟩
  exact (mem_coordinateProjection B S _).mpr ⟨x, hAB hxA, rfl⟩

/-- Projecting onto all coordinates preserves cardinality. -/
theorem coordinateProjection_univ_card {n : ℕ} (A : Finset (SignVector n)) :
    (coordinateProjection A Finset.univ).card = A.card := by
  classical
  let f : SignVector n → (Finset.univ : Finset (Fin n)) → Bool := restrictSigns Finset.univ
  have hf : Function.Injective f := by
    intro x y hxy
    funext i
    have hi := congrFun hxy ⟨i, by simp⟩
    simpa [f, restrictSigns] using hi
  simpa [coordinateProjection, f] using Finset.card_image_iff.mpr hf

/-- Forgetting one coordinate. -/
def forgetCoordinate {n : ℕ} (i : Fin n) : Finset (Fin n) :=
  Finset.univ.erase i

/-- Every projected pattern has at least one successful completion. -/
theorem exists_completion_of_mem_projection {n : ℕ} (A : Finset (SignVector n))
    (S : Finset (Fin n)) {y : S → Bool}
    (hy : y ∈ coordinateProjection A S) :
    ∃ x ∈ A, restrictSigns S x = y :=
  (mem_coordinateProjection A S y).mp hy

/-- A one-coordinate projection can have at most the same cardinality as `A`,
and at least half as many elements because every projected pattern has at most two lifts. -/
theorem card_le_two_mul_projection_forget {n : ℕ} (A : Finset (SignVector n))
    (i : Fin n) :
    A.card ≤ 2 * (coordinateProjection A (forgetCoordinate i)).card := by
  classical
  let f : SignVector n → (forgetCoordinate i → Bool) := restrictSigns (forgetCoordinate i)
  have hfiber : ∀ y, ((A.filter fun x ↦ f x = y).card) ≤ 2 := by
    intro y
    have hinj : Function.Injective fun x : {x // x ∈ A.filter fun x ↦ f x = y} ↦ x.1 i := by
      intro x z hxz
      apply Subtype.ext
      funext j
      by_cases hji : j = i
      · subst j
        exact hxz
      · have hxproj := x.2.2
        have hzproj := z.2.2
        have hjmem : j ∈ forgetCoordinate i := by
          simp [forgetCoordinate, hji]
        have hcoord := congrFun (hxproj.trans hzproj.symm) ⟨j, hjmem⟩
        simpa [f, restrictSigns] using hcoord
    calc
      (A.filter fun x ↦ f x = y).card = Fintype.card {x // x ∈ A.filter fun x ↦ f x = y} := by simp
      _ ≤ Fintype.card Bool := Fintype.card_le_of_injective _ hinj
      _ = 2 := by decide
  have hsum : A.card = ∑ y ∈ coordinateProjection A (forgetCoordinate i),
      (A.filter fun x ↦ f x = y).card := by
    rw [← Finset.sum_fiberwise_eq_card_filter]
    simp [coordinateProjection, f]
  rw [hsum]
  calc
    ∑ y ∈ coordinateProjection A (forgetCoordinate i),
        (A.filter fun x ↦ f x = y).card
      ≤ ∑ _y ∈ coordinateProjection A (forgetCoordinate i), 2 := by
          exact Finset.sum_le_sum fun y _hy ↦ hfiber y
    _ = 2 * (coordinateProjection A (forgetCoordinate i)).card := by simp [Nat.mul_comm]

/-- The corresponding lower bound on a one-coordinate projection. -/
theorem half_card_le_projection_forget {n : ℕ} (A : Finset (SignVector n))
    (i : Fin n) :
    A.card / 2 ≤ (coordinateProjection A (forgetCoordinate i)).card := by
  omega

end GaussianChain.PrimeReveal
