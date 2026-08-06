import GaussianChain.PrimeRevealBlockFibers

namespace GaussianChain.PrimeReveal

/-- The unnormalized second collision moment of the projection onto `T`.
If the realized block fibers have sizes `m_y`, this is `∑_y m_y^2`. -/
def projectionCollisionNumerator {n : ℕ} (A : Finset (SignVector n))
    (T : Finset (Fin n)) : ℕ :=
  ∑ y ∈ coordinateProjection A T, (blockFiber A T y).card ^ 2

/-- More generally, the unnormalized `k`-replica agreement moment. -/
def projectionCollisionMoment {n : ℕ} (k : ℕ) (A : Finset (SignVector n))
    (T : Finset (Fin n)) : ℕ :=
  ∑ y ∈ coordinateProjection A T, (blockFiber A T y).card ^ k

@[simp]
theorem projectionCollisionMoment_two {n : ℕ} (A : Finset (SignVector n))
    (T : Finset (Fin n)) :
    projectionCollisionMoment 2 A T = projectionCollisionNumerator A T := by
  rfl

/-- The first collision moment is exactly the size of the successful family. -/
theorem projectionCollisionMoment_one {n : ℕ} (A : Finset (SignVector n))
    (T : Finset (Fin n)) :
    projectionCollisionMoment 1 A T = A.card := by
  simp [projectionCollisionMoment, card_eq_sum_blockFiber]

/-- Every realized fiber is nonempty, so the second collision moment dominates
its first moment. Equivalently, diagonal pairs already contribute `A.card`. -/
theorem card_le_projectionCollisionNumerator {n : ℕ} (A : Finset (SignVector n))
    (T : Finset (Fin n)) :
    A.card ≤ projectionCollisionNumerator A T := by
  rw [card_eq_sum_blockFiber A T]
  apply Finset.sum_le_sum
  intro y hy
  have hpos : 0 < (blockFiber A T y).card :=
    Finset.card_pos.mpr (blockFiber_nonempty_of_mem_projection A T hy)
  dsimp [projectionCollisionNumerator]
  nlinarith

/-- The collision numerator counts ordered pairs of successful assignments
whose restrictions to `T` agree. -/
def agreeingPairCount {n : ℕ} (A : Finset (SignVector n))
    (T : Finset (Fin n)) : ℕ :=
  ((A ×ˢ A).filter fun xy ↦ restrictSigns T xy.1 = restrictSigns T xy.2).card

end GaussianChain.PrimeReveal
