import GaussianChain.PrimeReveal

namespace GaussianChain.PrimeReveal

/-- An order-free partial assignment of Gaussian-prime orientations.
`none` means that the prime has not yet been revealed. -/
abbrev PartialSignVector (n : ℕ) := Fin n → Option Bool

/-- A full orientation pattern is compatible with a partial assignment when
it agrees at every revealed coordinate. -/
def Compatible {n : ℕ} (a : PartialSignVector n) (x : SignVector n) : Prop :=
  ∀ i, a i = none ∨ a i = some (x i)

/-- Reveal coordinate `i` with value `b`. -/
def reveal {n : ℕ} (a : PartialSignVector n) (i : Fin n) (b : Bool) :
    PartialSignVector n :=
  Function.update a i (some b)

/-- Successful full assignments extending the partial state. -/
def completionFiber {n : ℕ} (A : Finset (SignVector n)) (a : PartialSignVector n) :
    Finset (SignVector n) :=
  A.filter (Compatible a)

/-- The canonical completion count attached to a partial state. -/
def completionCount {n : ℕ} (A : Finset (SignVector n)) (a : PartialSignVector n) : ℕ :=
  (completionFiber A a).card

@[simp]
theorem mem_completionFiber {n : ℕ} (A : Finset (SignVector n))
    (a : PartialSignVector n) (x : SignVector n) :
    x ∈ completionFiber A a ↔ x ∈ A ∧ Compatible a x := by
  simp [completionFiber]

@[simp]
theorem reveal_apply_same {n : ℕ} (a : PartialSignVector n)
    (i : Fin n) (b : Bool) :
    reveal a i b i = some b := by
  simp [reveal]

@[simp]
theorem reveal_apply_of_ne {n : ℕ} (a : PartialSignVector n)
    (i j : Fin n) (b : Bool) (h : j ≠ i) :
    reveal a i b j = a j := by
  simp [reveal, h]

/-- Reveals at distinct primes commute. This is the formal order-invariance
of the partial Gaussian factorization. -/
theorem reveal_commute {n : ℕ} (a : PartialSignVector n)
    (i j : Fin n) (bi bj : Bool) (h : i ≠ j) :
    reveal (reveal a i bi) j bj = reveal (reveal a j bj) i bi := by
  funext k
  by_cases hki : k = i
  · subst k
    simp [reveal, h]
  · by_cases hkj : k = j
    · subst k
      simp [reveal, h, hki]
    · simp [reveal, hki, hkj]

/-- Consequently the successful completion fiber is independent of the
order in which two distinct primes are revealed. -/
theorem completionFiber_reveal_commute {n : ℕ} (A : Finset (SignVector n))
    (a : PartialSignVector n) (i j : Fin n) (bi bj : Bool) (h : i ≠ j) :
    completionFiber A (reveal (reveal a i bi) j bj) =
      completionFiber A (reveal (reveal a j bj) i bi) := by
  rw [reveal_commute a i j bi bj h]

/-- The corresponding completion counts are order-independent. -/
theorem completionCount_reveal_commute {n : ℕ} (A : Finset (SignVector n))
    (a : PartialSignVector n) (i j : Fin n) (bi bj : Bool) (h : i ≠ j) :
    completionCount A (reveal (reveal a i bi) j bj) =
      completionCount A (reveal (reveal a j bj) i bi) := by
  rw [completionCount, completionCount, reveal_commute a i j bi bj h]

/-- At an unrevealed coordinate, compatibility after revealing `b` is the
same as old compatibility together with the chosen Boolean value. -/
theorem compatible_reveal_iff {n : ℕ} (a : PartialSignVector n)
    (i : Fin n) (b : Bool) (x : SignVector n) (hi : a i = none) :
    Compatible (reveal a i b) x ↔ Compatible a x ∧ x i = b := by
  constructor
  · intro hx
    constructor
    · intro j
      by_cases hji : j = i
      · subst j
        exact Or.inl hi
      · simpa [reveal, hji] using hx j
    · have hxi := hx i
      simpa [reveal, eq_comm] using hxi
  · rintro ⟨ha, hxi⟩ j
    by_cases hji : j = i
    · subst j
      right
      simp [reveal, hxi]
    · simpa [reveal, hji] using ha j

/-- Revealing an unrevealed prime partitions the current completion fiber
into its two orientation children. -/
theorem completionFiber_eq_reveal_union {n : ℕ} (A : Finset (SignVector n))
    (a : PartialSignVector n) (i : Fin n) (hi : a i = none) :
    completionFiber A a =
      completionFiber A (reveal a i true) ∪
        completionFiber A (reveal a i false) := by
  ext x
  cases hxi : x i <;>
    simp [completionFiber, compatible_reveal_iff a i true x hi,
      compatible_reveal_iff a i false x hi, hxi]

/-- The two children of an unrevealed coordinate are disjoint. -/
theorem completionFiber_reveal_disjoint {n : ℕ} (A : Finset (SignVector n))
    (a : PartialSignVector n) (i : Fin n) (hi : a i = none) :
    Disjoint (completionFiber A (reveal a i true))
      (completionFiber A (reveal a i false)) := by
  rw [Finset.disjoint_left]
  intro x htrue hfalse
  have ht : x i = true :=
    ((compatible_reveal_iff a i true x hi).mp
      ((mem_completionFiber A (reveal a i true) x).mp htrue).2).2
  have hf : x i = false :=
    ((compatible_reveal_iff a i false x hi).mp
      ((mem_completionFiber A (reveal a i false) x).mp hfalse).2).2
  simp [ht] at hf

/-- Exact one-step dynamic-programming identity for the canonical completion
count. It holds for every unrevealed prime, with no reveal order chosen. -/
theorem completionCount_eq_children {n : ℕ} (A : Finset (SignVector n))
    (a : PartialSignVector n) (i : Fin n) (hi : a i = none) :
    completionCount A a =
      completionCount A (reveal a i true) +
        completionCount A (reveal a i false) := by
  rw [completionCount, completionCount, completionCount,
    completionFiber_eq_reveal_union A a i hi,
    Finset.card_union_of_disjoint (completionFiber_reveal_disjoint A a i hi)]

end GaussianChain.PrimeReveal
