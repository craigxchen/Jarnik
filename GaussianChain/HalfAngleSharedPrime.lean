import Mathlib

namespace GaussianChain
namespace HalfAngleSharedPrime

/--
Let `D₁,D₂` be two primitive half-angle norms and suppose their squared chord
has the integral form

`d * D₁ * D₂ = 4 * N * Δ²`.

If an odd prime `p` occurs exactly once in `N` and divides both `D₁` and `D₂`,
then `p` divides the Plücker determinant `Δ`.  The hypotheses are stated with
`N = p * N'`; `¬ p ∣ 4` isolates the odd-prime case and `¬ p ∣ N'` is the
squarefree-at-`p` condition.
-/
theorem prime_dvd_plucker_of_chord_integrality
    {p N' D₁ D₂ Δ d : ℕ}
    (hp : Nat.Prime p)
    (hp4 : ¬ p ∣ 4)
    (hN' : ¬ p ∣ N')
    (hD₁ : p ∣ D₁)
    (hD₂ : p ∣ D₂)
    (hchord : d * D₁ * D₂ = 4 * (p * N') * Δ ^ 2) :
    p ∣ Δ := by
  rcases hD₁ with ⟨a, rfl⟩
  rcases hD₂ with ⟨b, rfl⟩
  have hcancel : d * p * a * b = 4 * N' * Δ ^ 2 := by
    apply Nat.eq_of_mul_eq_mul_left hp.pos
    calc
      p * (d * p * a * b) = d * (p * a) * (p * b) := by ring
      _ = 4 * (p * N') * Δ ^ 2 := hchord
      _ = p * (4 * N' * Δ ^ 2) := by ring
  have hp_rhs₀ : p ∣ 4 * N' * Δ ^ 2 := by
    rw [← hcancel]
    exact ⟨d * a * b, by ring⟩
  have hp_rhs : p ∣ 4 * (N' * Δ ^ 2) := by
    simpa [mul_assoc] using hp_rhs₀
  rcases hp.dvd_mul.mp hp_rhs with hp_four | hp_tail
  · exact (hp4 hp_four).elim
  · rcases hp.dvd_mul.mp hp_tail with hp_N' | hp_sq
    · exact (hN' hp_N').elim
    · exact hp.dvd_of_dvd_pow hp_sq

/-- Coordinate determinant notation for the same shared-prime constraint. -/
theorem prime_dvd_det2_of_chord_integrality
    {p N' D₁ D₂ d q₁ r₁ q₂ r₂ : ℕ}
    (hp : Nat.Prime p)
    (hp4 : ¬ p ∣ 4)
    (hN' : ¬ p ∣ N')
    (hD₁ : p ∣ D₁)
    (hD₂ : p ∣ D₂)
    (hchord :
      d * D₁ * D₂ =
        4 * (p * N') * (q₁ * r₂ - q₂ * r₁) ^ 2) :
    p ∣ q₁ * r₂ - q₂ * r₁ :=
  prime_dvd_plucker_of_chord_integrality hp hp4 hN' hD₁ hD₂ hchord

end HalfAngleSharedPrime
end GaussianChain
