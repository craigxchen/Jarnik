import GaussianChain.AsymptoticParameters
import GaussianChain.MertensMain

namespace GaussianChain
namespace ArcToParameter

open Filter
open SubcriticalBound

/--
A finite enumeration of the lattice points on a circular arc, sorted by arclength.

The fields are the elementary geometric facts supplied by an arclength parametrization:
the points lie on the Gaussian-integer circle of norm `N`, they are distinct, their
arclength coordinates are monotone in the chosen order, chord length is bounded by
arclength separation, and all coordinates lie in an interval of length `L`.
-/
structure ArclengthCircularArc (M N : ℕ) (a L : ℝ) where
  point : Fin M → GaussianInt
  coord : Fin M → ℝ
  on_circle : ∀ i : Fin M, point i * star (point i) = ((N : ℤ) : GaussianInt)
  injective : Function.Injective point
  monotone_coord : ∀ i j : Fin M, (i : ℕ) ≤ (j : ℕ) → coord i ≤ coord j
  chord_le_arclength :
    ∀ i j : Fin M, gaussianSqDist (point i) (point j) ≤ (coord j - coord i) ^ 2
  in_interval : ∀ i : Fin M, a ≤ coord i ∧ coord i ≤ a + L

namespace ArclengthCircularArc

variable {M N : ℕ} {a L : ℝ}

/-- Extend the finite arc enumeration to the sequence shape used by the descent theorem. -/
def pointSeq (A : ArclengthCircularArc M N a L) : ℕ → GaussianInt :=
  fun n => if h : n < M then A.point ⟨n, h⟩ else 0

/-- Extend the finite arclength coordinates to the sequence shape used by the descent theorem. -/
def coordSeq (A : ArclengthCircularArc M N a L) : ℕ → ℝ :=
  fun n => if h : n < M then A.coord ⟨n, h⟩ else 0

@[simp] theorem pointSeq_of_lt (A : ArclengthCircularArc M N a L) {n : ℕ} (h : n < M) :
    A.pointSeq n = A.point ⟨n, h⟩ := by
  simp [pointSeq, h]

@[simp] theorem coordSeq_of_lt (A : ArclengthCircularArc M N a L) {n : ℕ} (h : n < M) :
    A.coordSeq n = A.coord ⟨n, h⟩ := by
  simp [coordSeq, h]

/-- Circle membership transfers from the finite arc to its extended point sequence. -/
theorem onCircleUpTo_pointSeq (A : ArclengthCircularArc M N a L) :
    OnCircleUpTo M N A.pointSeq := by
  intro i hi
  rw [pointSeq_of_lt A hi]
  exact A.on_circle ⟨i, hi⟩

/--
Distinctness transfers from the finite arc to the first `M` entries of its extended
point sequence.  This is the injectivity condition needed by the parametrized theorem.
-/
theorem injectiveUpTo_pointSeq (A : ArclengthCircularArc M N a L) :
    InjectiveUpTo M A.pointSeq := by
  intro i j hi hj hzij
  have hpoint : A.point ⟨i, hi⟩ = A.point ⟨j, hj⟩ := by
    simpa [pointSeq_of_lt A hi, pointSeq_of_lt A hj] using hzij
  exact congrArg Fin.val (A.injective hpoint)

/-- Monotone finite arclength coordinates transfer to the extended coordinate sequence. -/
theorem coordSeq_monotone (A : ArclengthCircularArc M N a L) :
    ∀ i j, i ≤ j → j < M → A.coordSeq i ≤ A.coordSeq j := by
  intro i j hij hj
  have hi : i < M := lt_of_le_of_lt hij hj
  simpa [coordSeq_of_lt A hi, coordSeq_of_lt A hj] using
    A.monotone_coord ⟨i, hi⟩ ⟨j, hj⟩ hij

/--
The finite chord-versus-arclength estimate transfers to the extended sequences.
This is the precise arclength-to-parameter inequality consumed by the descent proof.
-/
theorem sqDist_le_coordSeq (A : ArclengthCircularArc M N a L) :
    ∀ i j, i < M → j < M →
      gaussianSqDist (A.pointSeq i) (A.pointSeq j) ≤ (A.coordSeq j - A.coordSeq i) ^ 2 := by
  intro i j hi hj
  simpa [pointSeq_of_lt A hi, pointSeq_of_lt A hj, coordSeq_of_lt A hi,
    coordSeq_of_lt A hj] using A.chord_le_arclength ⟨i, hi⟩ ⟨j, hj⟩

/-- Interval containment transfers from the finite arc to the extended coordinate sequence. -/
theorem coordSeq_mem_interval (A : ArclengthCircularArc M N a L) :
    ∀ i, i < M → a ≤ A.coordSeq i ∧ A.coordSeq i ≤ a + L := by
  intro i hi
  simpa [coordSeq_of_lt A hi] using A.in_interval ⟨i, hi⟩

end ArclengthCircularArc

/--
Radius-based circular-arc form of the sublogarithmic bound.

For every fixed `C > 0` and every `D > 8`, all sufficiently large radii `R`
have the following property: any arclength-ordered finite family of lattice points
on the circle of radius `R`, contained in an arc of length at most `C * sqrt R`,
has cardinality at most `D * log R / log log R`.
-/
theorem eventually_jarnik_arclength_circular_arc_sublog_of_constant_gt_eight
    (C D : ℝ) (hC : 0 < C) (hD : (8 : ℝ) < D) :
    ∀ᶠ R : ℝ in atTop,
      ∀ {M N : ℕ} {a L : ℝ},
        R ^ 2 = (N : ℝ) →
        0 < M →
        0 < N →
        L ≤ C * Real.sqrt R →
        (A : ArclengthCircularArc M N a L) →
        (M : ℝ) ≤ D * (Real.log R / Real.log (Real.log R)) := by
  have hmain :=
    MertensMain.mertens_asymptotic_parametrized_theorem C D hC hD
  filter_upwards [hmain] with R hR
  intro M N a L hR2 hM hN hL A
  exact hR (z := A.pointSeq) (t := A.coordSeq) hR2 hM hN hL
    A.onCircleUpTo_pointSeq A.injectiveUpTo_pointSeq A.coordSeq_monotone
    A.sqDist_le_coordSeq A.coordSeq_mem_interval

/-- The same radius-based circular-arc form with the rounded constant `10`. -/
theorem eventually_jarnik_arclength_circular_arc_sublog (C : ℝ) (hC : 0 < C) :
    ∀ᶠ R : ℝ in atTop,
      ∀ {M N : ℕ} {a L : ℝ},
        R ^ 2 = (N : ℝ) →
        0 < M →
        0 < N →
        L ≤ C * Real.sqrt R →
        (A : ArclengthCircularArc M N a L) →
        (M : ℝ) ≤ 10 * (Real.log R / Real.log (Real.log R)) := by
  have hmain := MertensMain.eventually_jarnik_arc_sublog_mertens C hC
  filter_upwards [hmain] with R hR
  intro M N a L hR2 hM hN hL A
  exact hR (z := A.pointSeq) (t := A.coordSeq) hR2 hM hN hL
    A.onCircleUpTo_pointSeq A.injectiveUpTo_pointSeq A.coordSeq_monotone
    A.sqDist_le_coordSeq A.coordSeq_mem_interval

end ArcToParameter
end GaussianChain
