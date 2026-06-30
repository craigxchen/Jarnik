import GaussianChain.MissingPrimeIntervalBranch
import GaussianChain.PrimeDescentDichotomy

namespace GaussianChain
namespace MainDichotomy

open SubcriticalBound
open PrimeDescentDichotomy

/-- The weighted-log condition that makes the missing-prime determinant branch fire on an
interval `m < p ≤ U`. -/
noncomputable def manyMissingWeightedLogCondition (s N m U : ℕ) (A : ℝ) : Prop :=
  ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor (A ^ 2) : ℝ) <
    ((s : ℝ) ^ 2 * MertensLower.weightedMissingPrimeInterval N m U -
        2 * MissingPrimeIntervalBranch.missingPrimeLogSum N m U) +
      ((s * s : ℕ) : ℝ) * Real.log (N : ℝ)

/-- The finite-stack "few missing primes" condition. -/
noncomputable def fewMissingStackCondition (N m₀ k : ℕ) : Prop :=
  MertensLower.weightedMissingPrimeInterval N m₀ (IntervalStack.geomLower m₀ k) <
    (k : ℝ) * (Real.log 2 / 2)

/-- If a lower bound `B` for the missing-prime mass is already enough to make the
weighted-log determinant inequality true, then the actual missing-prime mass also makes the
many-missing branch true. -/
theorem manyMissingWeightedLogCondition_of_missing_lower_bound
    {s N m U : ℕ} {A B : ℝ}
    (hmissing : B ≤ MertensLower.weightedMissingPrimeInterval N m U)
    (hineq :
      ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor (A ^ 2) : ℝ) <
        ((s : ℝ) ^ 2 * B - 2 * MissingPrimeIntervalBranch.missingPrimeLogSum N m U) +
          ((s * s : ℕ) : ℝ) * Real.log (N : ℝ)) :
    manyMissingWeightedLogCondition s N m U A := by
  have hs_nonneg : 0 ≤ (s : ℝ) ^ 2 := sq_nonneg _
  have hgain :
      ((s : ℝ) ^ 2 * B - 2 * MissingPrimeIntervalBranch.missingPrimeLogSum N m U) +
          ((s * s : ℕ) : ℝ) * Real.log (N : ℝ) ≤
        ((s : ℝ) ^ 2 * MertensLower.weightedMissingPrimeInterval N m U -
            2 * MissingPrimeIntervalBranch.missingPrimeLogSum N m U) +
          ((s * s : ℕ) : ℝ) * Real.log (N : ℝ) := by
    nlinarith [mul_le_mul_of_nonneg_left hmissing hs_nonneg]
  exact hineq.trans_le hgain

/-- Finite missing-prime dichotomy from a threshold inequality.

Either the total missing-prime mass on the geometric stack is below the stack threshold, or it is
at least that threshold; in the latter case `hthreshold` turns it into the determinant branch. -/
theorem missing_or_few_stack_dichotomy_of_threshold
    {s N m₀ k : ℕ} {A : ℝ}
    (hthreshold :
      ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor (A ^ 2) : ℝ) <
        ((s : ℝ) ^ 2 * ((k : ℝ) * (Real.log 2 / 2)) -
            2 * MissingPrimeIntervalBranch.missingPrimeLogSum N m₀
              (IntervalStack.geomLower m₀ k)) +
          ((s * s : ℕ) : ℝ) * Real.log (N : ℝ)) :
    manyMissingWeightedLogCondition s N m₀ (IntervalStack.geomLower m₀ k) A ∨
      fewMissingStackCondition N m₀ k := by
  classical
  by_cases hfew : fewMissingStackCondition N m₀ k
  · exact Or.inr hfew
  · exact Or.inl
      (manyMissingWeightedLogCondition_of_missing_lower_bound
        (s := s) (N := N) (m := m₀) (U := IntervalStack.geomLower m₀ k)
        (A := A) (B := (k : ℝ) * (Real.log 2 / 2))
        (le_of_not_gt hfew) hthreshold)

/-- Version of `missing_or_few_stack_dichotomy_of_threshold` using the crude bound
`missingPrimeLogSum ≤ (U + 1) log U`, where `U` is the geometric stack upper endpoint. -/
theorem missing_or_few_stack_dichotomy_of_crude_threshold
    {s N m₀ k : ℕ} {A : ℝ}
    (hU : 1 ≤ IntervalStack.geomLower m₀ k)
    (hthreshold :
      ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor (A ^ 2) : ℝ) <
        ((s : ℝ) ^ 2 * ((k : ℝ) * (Real.log 2 / 2)) -
            2 * (((IntervalStack.geomLower m₀ k + 1 : ℕ) : ℝ) *
              Real.log (IntervalStack.geomLower m₀ k : ℝ))) +
          ((s * s : ℕ) : ℝ) * Real.log (N : ℝ)) :
    manyMissingWeightedLogCondition s N m₀ (IntervalStack.geomLower m₀ k) A ∨
      fewMissingStackCondition N m₀ k := by
  have hlogsum :=
    MissingPrimeIntervalBranch.missingPrimeLogSum_le_succ_mul_log
      (N := N) (m := m₀) (U := IntervalStack.geomLower m₀ k) hU
  refine missing_or_few_stack_dichotomy_of_threshold
    (s := s) (N := N) (m₀ := m₀) (k := k) (A := A) ?_
  have hle :
      ((s : ℝ) ^ 2 * ((k : ℝ) * (Real.log 2 / 2)) -
            2 * (((IntervalStack.geomLower m₀ k + 1 : ℕ) : ℝ) *
              Real.log (IntervalStack.geomLower m₀ k : ℝ))) +
          ((s * s : ℕ) : ℝ) * Real.log (N : ℝ) ≤
        ((s : ℝ) ^ 2 * ((k : ℝ) * (Real.log 2 / 2)) -
            2 * MissingPrimeIntervalBranch.missingPrimeLogSum N m₀
              (IntervalStack.geomLower m₀ k)) +
          ((s * s : ℕ) : ℝ) * Real.log (N : ℝ) := by
    nlinarith [hlogsum]
  exact hthreshold.trans_le hle

/-- Split-prime logarithmic descent condition reduced to the stack upper endpoint. -/
theorem logSplit_condition_of_upper_endpoint
    {s N m₀ U : ℕ} {A : ℝ}
    (hupper :
      ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor (A ^ 2) : ℝ) <
        ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - Real.log (U : ℝ))) :
    ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ U → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) * (N' : ℤ) →
          ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor (A ^ 2) : ℝ) <
            ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - Real.log (p : ℝ)) := by
  intro p _N' hp _hmp hpU _hpN _hN' _hfactor
  have hlog_le : Real.log (p : ℝ) ≤ Real.log (U : ℝ) :=
    Real.log_le_log (by exact_mod_cast hp.pos) (by exact_mod_cast hpU)
  have hs_nonneg : 0 ≤ ((s * s : ℕ) : ℝ) := by positivity
  have hsub :
      Real.log (N : ℝ) - Real.log (U : ℝ) ≤
        Real.log (N : ℝ) - Real.log (p : ℝ) :=
    sub_le_sub_left hlog_le (Real.log (N : ℝ))
  have hmul :=
    mul_le_mul_of_nonneg_left hsub hs_nonneg
  exact hupper.trans_le hmul

/-- Inert-prime logarithmic descent condition reduced to the stack upper endpoint. -/
theorem logInert_condition_of_upper_endpoint
    {s N m₀ U : ℕ} {A : ℝ}
    (hupper :
      ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor (A ^ 2) : ℝ) <
        ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - 2 * Real.log (U : ℝ))) :
    ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ U → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) ^ 2 * (N' : ℤ) →
          ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor (A ^ 2) : ℝ) <
            ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - 2 * Real.log (p : ℝ)) := by
  intro p _N' hp _hmp hpU _hpN _hN' _hfactor
  have hlog_le : Real.log (p : ℝ) ≤ Real.log (U : ℝ) :=
    Real.log_le_log (by exact_mod_cast hp.pos) (by exact_mod_cast hpU)
  have hs_nonneg : 0 ≤ ((s * s : ℕ) : ℝ) := by positivity
  have hsub :
      Real.log (N : ℝ) - 2 * Real.log (U : ℝ) ≤
        Real.log (N : ℝ) - 2 * Real.log (p : ℝ) := by
    nlinarith [hlog_le]
  have hmul :=
    mul_le_mul_of_nonneg_left hsub hs_nonneg
  exact hupper.trans_le hmul

/-- Finite interval-stack dichotomy combining the two branches already formalized.

For a fixed geometric stack, either the missing-prime determinant estimate directly bounds the
number of points, or the stack contains a prime divisor of `N` and the appropriate split/inert
descent bound holds for that prime. -/
theorem interval_stack_missing_or_descent_dichotomy
    {M s N m₀ k : ℕ}
    {a L A : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hdichotomy :
      manyMissingWeightedLogCondition s N m₀ (IntervalStack.geomLower m₀ k) A ∨
        fewMissingStackCondition N m₀ k)
    (hk : 0 < k)
    (hm₀ : 0 < m₀)
    (h2m₀ : 2 ≤ m₀)
    (hM : 0 < M)
    (hlarge : 4 * s ≤ M)
    (hA : 0 < A)
    (hB : 0 < Nat.floor (A ^ 2))
    (hsmallPrime : 4 * IntervalStack.geomLower m₀ k ≤ s)
    (hN : 0 < N)
    (herr : ∀ i : Fin k,
      Real.log (((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) + 1) +
          2 * Real.sqrt ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) *
            Real.log ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) ≤
        (Real.log 2 / 4) * ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ))
    (hsmallSplit : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) * (N' : ℤ) →
          (Nat.floor (A ^ 2)) ^ (s * (2 * s + 1)) < N' ^ (s * s))
    (hsmallInert : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) ^ 2 * (N' : ℤ) →
          (Nat.floor (A ^ 2)) ^ (s * (2 * s + 1)) < N' ^ (s * s))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hmono : ∀ i j, i ≤ j → j < M → t i ≤ t j)
    (hparam : ∀ i j, i < M → j < M → gaussianSqDist (z i) (z j) ≤ (t j - t i) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    ((M : ℝ) ≤ ((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) * L / A) ∨
      ∃ p, Nat.Prime p ∧ m₀ < p ∧ p ≤ IntervalStack.geomLower m₀ k ∧ p ∣ N ∧
        ((p % 4 = 3 ∧
            (M : ℝ) ≤
              ((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) * (L / (p : ℝ)) / A) ∨
          (p % 4 = 1 ∧
            (M : ℝ) ≤
              2 * (((2 * s : ℕ) : ℝ) +
                ((2 * s : ℕ) : ℝ) * (L / Real.sqrt (p : ℝ)) / A))) := by
  rcases hdichotomy with hmany | hfew
  · have hk_many : 2 * s ≤ M := by omega
    exact Or.inl
      (MissingPrimeIntervalBranch.card_le_of_missing_prime_interval_weighted_log_bound
        (M := M) (s := s) (N := N) (m := m₀) (U := IntervalStack.geomLower m₀ k)
        (a := a) (L := L) (A := A) (z := z) (t := t)
        hk_many hA hB hsmallPrime hN hmany hcircle hz hmono hparam hmem)
  · exact Or.inr
      (PrimeDescentDichotomy.exists_prime_dvd_with_descent_bound_of_geometric_stack_total_missing_lt
        (M := M) (s := s) (N := N) (m₀ := m₀) (k := k)
        (a := a) (L := L) (A := A) (z := z) (t := t)
        hk hm₀ h2m₀ hM hlarge hA hN herr hfew hsmallSplit hsmallInert
        hcircle hz hmono hparam hmem)

/-- Uniform-bound version of `interval_stack_missing_or_descent_dichotomy`.

The descent alternative is compressed to a single lower-endpoint estimate using
`sqrt m₀ ≤ sqrt p` for every divisor prime found above `m₀`. -/
theorem interval_stack_missing_or_uniform_descent_bound
    {M s N m₀ k : ℕ}
    {a L A : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hdichotomy :
      manyMissingWeightedLogCondition s N m₀ (IntervalStack.geomLower m₀ k) A ∨
        fewMissingStackCondition N m₀ k)
    (hk : 0 < k)
    (hm₀ : 0 < m₀)
    (h2m₀ : 2 ≤ m₀)
    (hM : 0 < M)
    (hlarge : 4 * s ≤ M)
    (hA : 0 < A)
    (hB : 0 < Nat.floor (A ^ 2))
    (hL : 0 ≤ L)
    (hsmallPrime : 4 * IntervalStack.geomLower m₀ k ≤ s)
    (hN : 0 < N)
    (herr : ∀ i : Fin k,
      Real.log (((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) + 1) +
          2 * Real.sqrt ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) *
            Real.log ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) ≤
        (Real.log 2 / 4) * ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ))
    (hsmallSplit : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) * (N' : ℤ) →
          (Nat.floor (A ^ 2)) ^ (s * (2 * s + 1)) < N' ^ (s * s))
    (hsmallInert : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) ^ 2 * (N' : ℤ) →
          (Nat.floor (A ^ 2)) ^ (s * (2 * s + 1)) < N' ^ (s * s))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hmono : ∀ i j, i ≤ j → j < M → t i ≤ t j)
    (hparam : ∀ i j, i < M → j < M → gaussianSqDist (z i) (z j) ≤ (t j - t i) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    ((M : ℝ) ≤ ((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) * L / A) ∨
      (M : ℝ) ≤
        2 * (((2 * s : ℕ) : ℝ) +
          ((2 * s : ℕ) : ℝ) * (L / Real.sqrt (m₀ : ℝ)) / A) := by
  have hraw :=
    interval_stack_missing_or_descent_dichotomy
      (M := M) (s := s) (N := N) (m₀ := m₀) (k := k)
      (a := a) (L := L) (A := A) (z := z) (t := t)
      hdichotomy hk hm₀ h2m₀ hM hlarge hA hB hsmallPrime hN herr
      hsmallSplit hsmallInert hcircle hz hmono hparam hmem
  rcases hraw with hmany | hdesc
  · exact Or.inl hmany
  · rcases hdesc with ⟨p, _hp, hmp, _hpU, _hpN, hbranch⟩
    exact Or.inr
      (PrimeDescentDichotomy.descent_branch_le_uniform_lower_endpoint
        (M := M) (s := s) (m₀ := m₀) (p := p) (L := L) (A := A)
        (by omega) hmp hL hA hbranch)

/-- Logarithmic-threshold version of `interval_stack_missing_or_uniform_descent_bound`. -/
theorem interval_stack_missing_or_uniform_descent_bound_log
    {M s N m₀ k : ℕ}
    {a L A : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hdichotomy :
      manyMissingWeightedLogCondition s N m₀ (IntervalStack.geomLower m₀ k) A ∨
        fewMissingStackCondition N m₀ k)
    (hk : 0 < k)
    (hm₀ : 0 < m₀)
    (h2m₀ : 2 ≤ m₀)
    (hM : 0 < M)
    (hlarge : 4 * s ≤ M)
    (hA : 0 < A)
    (hB : 0 < Nat.floor (A ^ 2))
    (hL : 0 ≤ L)
    (hsmallPrime : 4 * IntervalStack.geomLower m₀ k ≤ s)
    (hN : 0 < N)
    (herr : ∀ i : Fin k,
      Real.log (((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) + 1) +
          2 * Real.sqrt ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) *
            Real.log ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) ≤
        (Real.log 2 / 4) * ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ))
    (hlogSplit : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) * (N' : ℤ) →
          ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor (A ^ 2) : ℝ) <
            ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - Real.log (p : ℝ)))
    (hlogInert : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) ^ 2 * (N' : ℤ) →
          ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor (A ^ 2) : ℝ) <
            ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - 2 * Real.log (p : ℝ)))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hmono : ∀ i j, i ≤ j → j < M → t i ≤ t j)
    (hparam : ∀ i j, i < M → j < M → gaussianSqDist (z i) (z j) ≤ (t j - t i) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    ((M : ℝ) ≤ ((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) * L / A) ∨
      (M : ℝ) ≤
        2 * (((2 * s : ℕ) : ℝ) +
          ((2 * s : ℕ) : ℝ) * (L / Real.sqrt (m₀ : ℝ)) / A) := by
  exact interval_stack_missing_or_uniform_descent_bound
    (M := M) (s := s) (N := N) (m₀ := m₀) (k := k)
    (a := a) (L := L) (A := A) (z := z) (t := t)
    hdichotomy hk hm₀ h2m₀ hM hlarge hA hB hL hsmallPrime hN herr
    (fun p N' hp hmp hpU hpN hN' hfactor =>
      SubcriticalLog.split_natFloor_sq_pow_lt_of_log_sub
        (s := s) (N := N) (p := p) (N' := N') (A := A)
        hB hp.pos hN' hfactor (hlogSplit p N' hp hmp hpU hpN hN' hfactor))
    (fun p N' hp hmp hpU hpN hN' hfactor =>
      SubcriticalLog.inert_natFloor_sq_pow_lt_of_log_sub
        (s := s) (N := N) (p := p) (N' := N') (A := A)
        hB hp.pos hN' hfactor (hlogInert p N' hp hmp hpU hpN hN' hfactor))
    hcircle hz hmono hparam hmem

/-- Finite cardinality consequence of the stack dichotomy.

If the ordinary branch has `S * L / A ≤ S` and the descent branch has
`S * (L / sqrt m₀) / A ≤ S`, where `S = 2s`, then every case gives `M ≤ 8s`.
The case `M < 4s` is handled before invoking the dichotomy. -/
theorem card_le_eight_mul_s_of_interval_stack_log
    {M s N m₀ k : ℕ}
    {a L A : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hdichotomy :
      manyMissingWeightedLogCondition s N m₀ (IntervalStack.geomLower m₀ k) A ∨
        fewMissingStackCondition N m₀ k)
    (hk : 0 < k)
    (hm₀ : 0 < m₀)
    (h2m₀ : 2 ≤ m₀)
    (hM : 0 < M)
    (hA : 0 < A)
    (hB : 0 < Nat.floor (A ^ 2))
    (hL : 0 ≤ L)
    (hsmallPrime : 4 * IntervalStack.geomLower m₀ k ≤ s)
    (hN : 0 < N)
    (htermMany :
      ((2 * s : ℕ) : ℝ) * L / A ≤ ((2 * s : ℕ) : ℝ))
    (htermDesc :
      ((2 * s : ℕ) : ℝ) * (L / Real.sqrt (m₀ : ℝ)) / A ≤
        ((2 * s : ℕ) : ℝ))
    (herr : ∀ i : Fin k,
      Real.log (((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) + 1) +
          2 * Real.sqrt ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) *
            Real.log ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) ≤
        (Real.log 2 / 4) * ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ))
    (hlogSplit : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) * (N' : ℤ) →
          ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor (A ^ 2) : ℝ) <
            ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - Real.log (p : ℝ)))
    (hlogInert : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) ^ 2 * (N' : ℤ) →
          ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor (A ^ 2) : ℝ) <
            ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - 2 * Real.log (p : ℝ)))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hmono : ∀ i j, i ≤ j → j < M → t i ≤ t j)
    (hparam : ∀ i j, i < M → j < M → gaussianSqDist (z i) (z j) ≤ (t j - t i) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    (M : ℝ) ≤ (8 * s : ℕ) := by
  by_cases hlarge : 4 * s ≤ M
  · have hbound :=
      interval_stack_missing_or_uniform_descent_bound_log
        (M := M) (s := s) (N := N) (m₀ := m₀) (k := k)
        (a := a) (L := L) (A := A) (z := z) (t := t)
        hdichotomy hk hm₀ h2m₀ hM hlarge hA hB hL hsmallPrime hN herr
        hlogSplit hlogInert hcircle hz hmono hparam hmem
    rcases hbound with hmany | hdesc
    · calc
        (M : ℝ) ≤ ((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) * L / A := hmany
        _ ≤ ((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) := by
          nlinarith [htermMany]
        _ ≤ (8 * s : ℕ) := by
          norm_num
          nlinarith
    · calc
        (M : ℝ) ≤
            2 * (((2 * s : ℕ) : ℝ) +
              ((2 * s : ℕ) : ℝ) * (L / Real.sqrt (m₀ : ℝ)) / A) := hdesc
        _ ≤ 2 * (((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ)) := by
          nlinarith [htermDesc]
        _ ≤ (8 * s : ℕ) := by
          norm_num
          nlinarith
  · have hM_le : M ≤ 8 * s := by omega
    exact_mod_cast hM_le

/-- Finite cardinality consequence with an explicit branch factor.

This is the same dichotomy as `card_le_eight_mul_s_of_interval_stack_log`, but the two
length terms are only assumed to be bounded by `q` copies of the base window size
`S = 2s`. The resulting bound is `2 * (q + 1) * S`; this is the finite form used when
`q` is a fixed factor depending on the permitted arc length. -/
theorem card_le_two_mul_succ_q_mul_two_s_of_interval_stack_log
    {M s q N m₀ k : ℕ}
    {a L A : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hdichotomy :
      manyMissingWeightedLogCondition s N m₀ (IntervalStack.geomLower m₀ k) A ∨
        fewMissingStackCondition N m₀ k)
    (hk : 0 < k)
    (hm₀ : 0 < m₀)
    (h2m₀ : 2 ≤ m₀)
    (hM : 0 < M)
    (hA : 0 < A)
    (hB : 0 < Nat.floor (A ^ 2))
    (hL : 0 ≤ L)
    (hsmallPrime : 4 * IntervalStack.geomLower m₀ k ≤ s)
    (hN : 0 < N)
    (htermMany :
      ((2 * s : ℕ) : ℝ) * L / A ≤ (q : ℝ) * ((2 * s : ℕ) : ℝ))
    (htermDesc :
      ((2 * s : ℕ) : ℝ) * (L / Real.sqrt (m₀ : ℝ)) / A ≤
        (q : ℝ) * ((2 * s : ℕ) : ℝ))
    (herr : ∀ i : Fin k,
      Real.log (((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) + 1) +
          2 * Real.sqrt ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) *
            Real.log ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) ≤
        (Real.log 2 / 4) * ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ))
    (hlogSplit : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) * (N' : ℤ) →
          ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor (A ^ 2) : ℝ) <
            ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - Real.log (p : ℝ)))
    (hlogInert : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) ^ 2 * (N' : ℤ) →
          ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor (A ^ 2) : ℝ) <
            ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - 2 * Real.log (p : ℝ)))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hmono : ∀ i j, i ≤ j → j < M → t i ≤ t j)
    (hparam : ∀ i j, i < M → j < M → gaussianSqDist (z i) (z j) ≤ (t j - t i) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    (M : ℝ) ≤ 2 * (((q + 1 : ℕ) : ℝ) * ((2 * s : ℕ) : ℝ)) := by
  by_cases hlarge : 4 * s ≤ M
  · have hbound :=
      interval_stack_missing_or_uniform_descent_bound_log
        (M := M) (s := s) (N := N) (m₀ := m₀) (k := k)
        (a := a) (L := L) (A := A) (z := z) (t := t)
        hdichotomy hk hm₀ h2m₀ hM hlarge hA hB hL hsmallPrime hN herr
        hlogSplit hlogInert hcircle hz hmono hparam hmem
    have hS_nonneg : 0 ≤ ((2 * s : ℕ) : ℝ) := by positivity
    rcases hbound with hmany | hdesc
    · calc
        (M : ℝ) ≤ ((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) * L / A := hmany
        _ ≤ ((2 * s : ℕ) : ℝ) + (q : ℝ) * ((2 * s : ℕ) : ℝ) := by
          nlinarith [htermMany]
        _ ≤ 2 * (((q + 1 : ℕ) : ℝ) * ((2 * s : ℕ) : ℝ)) := by
          norm_num
          nlinarith [hS_nonneg]
    · calc
        (M : ℝ) ≤
            2 * (((2 * s : ℕ) : ℝ) +
              ((2 * s : ℕ) : ℝ) * (L / Real.sqrt (m₀ : ℝ)) / A) := hdesc
        _ ≤ 2 * (((2 * s : ℕ) : ℝ) + (q : ℝ) * ((2 * s : ℕ) : ℝ)) := by
          nlinarith [htermDesc]
        _ ≤ 2 * (((q + 1 : ℕ) : ℝ) * ((2 * s : ℕ) : ℝ)) := by
          norm_num
          ring_nf
          exact le_rfl
  · have hM_le : M ≤ 4 * s := by omega
    have hM_real : (M : ℝ) ≤ (4 * s : ℕ) := by exact_mod_cast hM_le
    have hS_nonneg : 0 ≤ ((2 * s : ℕ) : ℝ) := by positivity
    calc
      (M : ℝ) ≤ (4 * s : ℕ) := hM_real
      _ ≤ 2 * (((q + 1 : ℕ) : ℝ) * ((2 * s : ℕ) : ℝ)) := by
        norm_num
        nlinarith [hS_nonneg]

/-- Finite cardinality consequence with separate scales for the two branches.

The many-missing determinant branch uses the fixed scale `A₀`. In the few-missing branch, once a
divisor prime `p` is found, the descent subcritical scale may be `A p`. This is the final-use
shape: one can later take `A p = C * sqrt R / sqrt p`. -/
theorem card_le_two_mul_succ_q_mul_two_s_of_interval_stack_primeScale_log
    {M s q N m₀ k : ℕ}
    {a L A₀ : ℝ} {A : ℕ → ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hdichotomy :
      manyMissingWeightedLogCondition s N m₀ (IntervalStack.geomLower m₀ k) A₀ ∨
        fewMissingStackCondition N m₀ k)
    (hk : 0 < k)
    (hm₀ : 0 < m₀)
    (h2m₀ : 2 ≤ m₀)
    (hM : 0 < M)
    (hA₀ : 0 < A₀)
    (hB₀ : 0 < Nat.floor (A₀ ^ 2))
    (hsmallPrime : 4 * IntervalStack.geomLower m₀ k ≤ s)
    (hN : 0 < N)
    (htermMany :
      ((2 * s : ℕ) : ℝ) * L / A₀ ≤ (q : ℝ) * ((2 * s : ℕ) : ℝ))
    (htermInert : ∀ p, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N →
        ((2 * s : ℕ) : ℝ) * (L / (p : ℝ)) / A p ≤
          (q : ℝ) * ((2 * s : ℕ) : ℝ))
    (htermSplit : ∀ p, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N →
        ((2 * s : ℕ) : ℝ) * (L / Real.sqrt (p : ℝ)) / A p ≤
          (q : ℝ) * ((2 * s : ℕ) : ℝ))
    (herr : ∀ i : Fin k,
      Real.log (((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) + 1) +
          2 * Real.sqrt ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) *
            Real.log ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) ≤
        (Real.log 2 / 4) * ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ))
    (hA : ∀ p, Nat.Prime p → m₀ < p → p ≤ IntervalStack.geomLower m₀ k → p ∣ N →
      0 < A p)
    (hB : ∀ p, Nat.Prime p → m₀ < p → p ≤ IntervalStack.geomLower m₀ k → p ∣ N →
      0 < Nat.floor ((A p) ^ 2))
    (hlogSplit : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) * (N' : ℤ) →
          ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor ((A p) ^ 2) : ℝ) <
            ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - Real.log (p : ℝ)))
    (hlogInert : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) ^ 2 * (N' : ℤ) →
          ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor ((A p) ^ 2) : ℝ) <
            ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - 2 * Real.log (p : ℝ)))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hmono : ∀ i j, i ≤ j → j < M → t i ≤ t j)
    (hparam : ∀ i j, i < M → j < M → gaussianSqDist (z i) (z j) ≤ (t j - t i) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    (M : ℝ) ≤ 2 * (((q + 1 : ℕ) : ℝ) * ((2 * s : ℕ) : ℝ)) := by
  by_cases hlarge : 4 * s ≤ M
  · rcases hdichotomy with hmany | hfew
    · have hk_many : 2 * s ≤ M := by omega
      have hmany_bound :
          (M : ℝ) ≤ ((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) * L / A₀ :=
        MissingPrimeIntervalBranch.card_le_of_missing_prime_interval_weighted_log_bound
          (M := M) (s := s) (N := N) (m := m₀) (U := IntervalStack.geomLower m₀ k)
          (a := a) (L := L) (A := A₀) (z := z) (t := t)
          hk_many hA₀ hB₀ hsmallPrime hN hmany hcircle hz hmono hparam hmem
      have hS_nonneg : 0 ≤ ((2 * s : ℕ) : ℝ) := by positivity
      calc
        (M : ℝ) ≤ ((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) * L / A₀ := hmany_bound
        _ ≤ ((2 * s : ℕ) : ℝ) + (q : ℝ) * ((2 * s : ℕ) : ℝ) := by
          nlinarith [htermMany]
        _ ≤ 2 * (((q + 1 : ℕ) : ℝ) * ((2 * s : ℕ) : ℝ)) := by
          norm_num
          nlinarith [hS_nonneg]
    · have hdesc :=
        exists_prime_dvd_with_descent_bound_of_geometric_stack_total_missing_lt_log_primeScale
          (M := M) (s := s) (N := N) (m₀ := m₀) (k := k)
          (a := a) (L := L) (A := A) (z := z) (t := t)
          hk hm₀ h2m₀ hM hlarge hN herr hfew hA hB hlogSplit hlogInert
          hcircle hz hmono hparam hmem
      obtain ⟨p, hp, hmp, hpU, hpN, hbranch⟩ := hdesc
      exact PrimeDescentDichotomy.descent_branch_le_of_prime_scale_terms
        (M := M) (s := s) (p := p) (q := q) (L := L) (A := A p)
        (htermInert p hp hmp hpU hpN) (htermSplit p hp hmp hpU hpN) hbranch
  · have hM_le : M ≤ 4 * s := by omega
    have hM_real : (M : ℝ) ≤ (4 * s : ℕ) := by exact_mod_cast hM_le
    have hS_nonneg : 0 ≤ ((2 * s : ℕ) : ℝ) := by positivity
    calc
      (M : ℝ) ≤ (4 * s : ℕ) := hM_real
      _ ≤ 2 * (((q + 1 : ℕ) : ℝ) * ((2 * s : ℕ) : ℝ)) := by
        norm_num
        nlinarith [hS_nonneg]

/-- Two-scale version of
`card_le_two_mul_succ_q_mul_two_s_of_interval_stack_primeScale_log`.

The split branch uses `Asplit p`, while the inert branch uses `Ainert p`. -/
theorem card_le_two_mul_succ_q_mul_two_s_of_interval_stack_twoScale_log
    {M s q N m₀ k : ℕ}
    {a L A₀ : ℝ} {Asplit Ainert : ℕ → ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hdichotomy :
      manyMissingWeightedLogCondition s N m₀ (IntervalStack.geomLower m₀ k) A₀ ∨
        fewMissingStackCondition N m₀ k)
    (hk : 0 < k)
    (hm₀ : 0 < m₀)
    (h2m₀ : 2 ≤ m₀)
    (hM : 0 < M)
    (hA₀ : 0 < A₀)
    (hB₀ : 0 < Nat.floor (A₀ ^ 2))
    (hsmallPrime : 4 * IntervalStack.geomLower m₀ k ≤ s)
    (hN : 0 < N)
    (htermMany :
      ((2 * s : ℕ) : ℝ) * L / A₀ ≤ (q : ℝ) * ((2 * s : ℕ) : ℝ))
    (htermInert : ∀ p, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N →
        ((2 * s : ℕ) : ℝ) * (L / (p : ℝ)) / Ainert p ≤
          (q : ℝ) * ((2 * s : ℕ) : ℝ))
    (htermSplit : ∀ p, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N →
        ((2 * s : ℕ) : ℝ) * (L / Real.sqrt (p : ℝ)) / Asplit p ≤
          (q : ℝ) * ((2 * s : ℕ) : ℝ))
    (herr : ∀ i : Fin k,
      Real.log (((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) + 1) +
          2 * Real.sqrt ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) *
            Real.log ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) ≤
        (Real.log 2 / 4) * ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ))
    (hAsplit : ∀ p, Nat.Prime p → m₀ < p → p ≤ IntervalStack.geomLower m₀ k → p ∣ N →
      0 < Asplit p)
    (hBs : ∀ p, Nat.Prime p → m₀ < p → p ≤ IntervalStack.geomLower m₀ k → p ∣ N →
      0 < Nat.floor ((Asplit p) ^ 2))
    (hAinert : ∀ p, Nat.Prime p → m₀ < p → p ≤ IntervalStack.geomLower m₀ k → p ∣ N →
      0 < Ainert p)
    (hBi : ∀ p, Nat.Prime p → m₀ < p → p ≤ IntervalStack.geomLower m₀ k → p ∣ N →
      0 < Nat.floor ((Ainert p) ^ 2))
    (hlogSplit : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) * (N' : ℤ) →
          ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor ((Asplit p) ^ 2) : ℝ) <
            ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - Real.log (p : ℝ)))
    (hlogInert : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) ^ 2 * (N' : ℤ) →
          ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor ((Ainert p) ^ 2) : ℝ) <
            ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - 2 * Real.log (p : ℝ)))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hmono : ∀ i j, i ≤ j → j < M → t i ≤ t j)
    (hparam : ∀ i j, i < M → j < M → gaussianSqDist (z i) (z j) ≤ (t j - t i) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    (M : ℝ) ≤ 2 * (((q + 1 : ℕ) : ℝ) * ((2 * s : ℕ) : ℝ)) := by
  by_cases hlarge : 4 * s ≤ M
  · rcases hdichotomy with hmany | hfew
    · have hk_many : 2 * s ≤ M := by omega
      have hmany_bound :
          (M : ℝ) ≤ ((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) * L / A₀ :=
        MissingPrimeIntervalBranch.card_le_of_missing_prime_interval_weighted_log_bound
          (M := M) (s := s) (N := N) (m := m₀) (U := IntervalStack.geomLower m₀ k)
          (a := a) (L := L) (A := A₀) (z := z) (t := t)
          hk_many hA₀ hB₀ hsmallPrime hN hmany hcircle hz hmono hparam hmem
      have hS_nonneg : 0 ≤ ((2 * s : ℕ) : ℝ) := by positivity
      calc
        (M : ℝ) ≤ ((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) * L / A₀ := hmany_bound
        _ ≤ ((2 * s : ℕ) : ℝ) + (q : ℝ) * ((2 * s : ℕ) : ℝ) := by
          nlinarith [htermMany]
        _ ≤ 2 * (((q + 1 : ℕ) : ℝ) * ((2 * s : ℕ) : ℝ)) := by
          norm_num
          nlinarith [hS_nonneg]
    · have hdesc :=
        exists_prime_dvd_with_descent_bound_of_geometric_stack_total_missing_lt_log_twoScale
          (M := M) (s := s) (N := N) (m₀ := m₀) (k := k)
          (a := a) (L := L) (Asplit := Asplit) (Ainert := Ainert) (z := z) (t := t)
          hk hm₀ h2m₀ hM hlarge hN herr hfew hAsplit hBs hAinert hBi
          hlogSplit hlogInert hcircle hz hmono hparam hmem
      obtain ⟨p, hp, hmp, hpU, hpN, hbranch⟩ := hdesc
      exact PrimeDescentDichotomy.descent_branch_le_of_prime_twoScale_terms
        (M := M) (s := s) (p := p) (q := q) (L := L)
        (Asplit := Asplit p) (Ainert := Ainert p)
        (htermInert p hp hmp hpU hpN) (htermSplit p hp hmp hpU hpN) hbranch
  · have hM_le : M ≤ 4 * s := by omega
    have hM_real : (M : ℝ) ≤ (4 * s : ℕ) := by exact_mod_cast hM_le
    have hS_nonneg : 0 ≤ ((2 * s : ℕ) : ℝ) := by positivity
    calc
      (M : ℝ) ≤ (4 * s : ℕ) := hM_real
      _ ≤ 2 * (((q + 1 : ℕ) : ℝ) * ((2 * s : ℕ) : ℝ)) := by
        norm_num
        nlinarith [hS_nonneg]

/-- Finite cardinality consequence with a small split-branch length term.

The many-missing and inert branches each give at most `4s` when their sliding-window
length term is bounded by the base window size `2s`.  In the split branch there is an
unavoidable factor `2` from passing to one Gaussian-prime divisibility class; if the split
length term is only `eps` times the base window size, that branch gives
`(4 + 4 * eps) s`. -/
theorem card_le_four_add_split_eps_mul_s_of_interval_stack_twoScale_log
    {M s N m₀ k : ℕ}
    {a L A₀ eps : ℝ} {Asplit Ainert : ℕ → ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hdichotomy :
      manyMissingWeightedLogCondition s N m₀ (IntervalStack.geomLower m₀ k) A₀ ∨
        fewMissingStackCondition N m₀ k)
    (hk : 0 < k)
    (hm₀ : 0 < m₀)
    (h2m₀ : 2 ≤ m₀)
    (hM : 0 < M)
    (hA₀ : 0 < A₀)
    (hB₀ : 0 < Nat.floor (A₀ ^ 2))
    (hsmallPrime : 4 * IntervalStack.geomLower m₀ k ≤ s)
    (hN : 0 < N)
    (heps : 0 ≤ eps)
    (htermMany :
      ((2 * s : ℕ) : ℝ) * L / A₀ ≤ ((2 * s : ℕ) : ℝ))
    (htermInert : ∀ p, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N →
        ((2 * s : ℕ) : ℝ) * (L / (p : ℝ)) / Ainert p ≤
          ((2 * s : ℕ) : ℝ))
    (htermSplit : ∀ p, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N →
        ((2 * s : ℕ) : ℝ) * (L / Real.sqrt (p : ℝ)) / Asplit p ≤
          eps * ((2 * s : ℕ) : ℝ))
    (herr : ∀ i : Fin k,
      Real.log (((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) + 1) +
          2 * Real.sqrt ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) *
            Real.log ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) ≤
        (Real.log 2 / 4) * ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ))
    (hAsplit : ∀ p, Nat.Prime p → m₀ < p → p ≤ IntervalStack.geomLower m₀ k → p ∣ N →
      0 < Asplit p)
    (hBs : ∀ p, Nat.Prime p → m₀ < p → p ≤ IntervalStack.geomLower m₀ k → p ∣ N →
      0 < Nat.floor ((Asplit p) ^ 2))
    (hAinert : ∀ p, Nat.Prime p → m₀ < p → p ≤ IntervalStack.geomLower m₀ k → p ∣ N →
      0 < Ainert p)
    (hBi : ∀ p, Nat.Prime p → m₀ < p → p ≤ IntervalStack.geomLower m₀ k → p ∣ N →
      0 < Nat.floor ((Ainert p) ^ 2))
    (hlogSplit : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) * (N' : ℤ) →
          ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor ((Asplit p) ^ 2) : ℝ) <
            ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - Real.log (p : ℝ)))
    (hlogInert : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) ^ 2 * (N' : ℤ) →
          ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor ((Ainert p) ^ 2) : ℝ) <
            ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - 2 * Real.log (p : ℝ)))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hmono : ∀ i j, i ≤ j → j < M → t i ≤ t j)
    (hparam : ∀ i j, i < M → j < M → gaussianSqDist (z i) (z j) ≤ (t j - t i) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    (M : ℝ) ≤ (4 + 4 * eps) * (s : ℝ) := by
  by_cases hlarge : 4 * s ≤ M
  · rcases hdichotomy with hmany | hfew
    · have hk_many : 2 * s ≤ M := by omega
      have hmany_bound :
          (M : ℝ) ≤ ((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) * L / A₀ :=
        MissingPrimeIntervalBranch.card_le_of_missing_prime_interval_weighted_log_bound
          (M := M) (s := s) (N := N) (m := m₀) (U := IntervalStack.geomLower m₀ k)
          (a := a) (L := L) (A := A₀) (z := z) (t := t)
          hk_many hA₀ hB₀ hsmallPrime hN hmany hcircle hz hmono hparam hmem
      have hs_nonneg : 0 ≤ (s : ℝ) := by positivity
      calc
        (M : ℝ) ≤ ((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) * L / A₀ := hmany_bound
        _ ≤ ((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) := by
          nlinarith [htermMany]
        _ = 4 * (s : ℝ) := by norm_num; ring
        _ ≤ (4 + 4 * eps) * (s : ℝ) := by nlinarith
    · have hdesc :=
        exists_prime_dvd_with_descent_bound_of_geometric_stack_total_missing_lt_log_twoScale
          (M := M) (s := s) (N := N) (m₀ := m₀) (k := k)
          (a := a) (L := L) (Asplit := Asplit) (Ainert := Ainert) (z := z) (t := t)
          hk hm₀ h2m₀ hM hlarge hN herr hfew hAsplit hBs hAinert hBi
          hlogSplit hlogInert hcircle hz hmono hparam hmem
      obtain ⟨p, hp, hmp, hpU, hpN, hbranch⟩ := hdesc
      rcases hbranch with hbranch | hbranch
      · rcases hbranch with ⟨_hp3, hinert⟩
        have hs_nonneg : 0 ≤ (s : ℝ) := by positivity
        calc
          (M : ℝ) ≤
              ((2 * s : ℕ) : ℝ) +
                ((2 * s : ℕ) : ℝ) * (L / (p : ℝ)) / Ainert p := hinert
          _ ≤ ((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) := by
            nlinarith [htermInert p hp hmp hpU hpN]
          _ = 4 * (s : ℝ) := by norm_num; ring
          _ ≤ (4 + 4 * eps) * (s : ℝ) := by nlinarith
      · rcases hbranch with ⟨_hp1, hsplit⟩
        calc
          (M : ℝ) ≤
              2 * (((2 * s : ℕ) : ℝ) +
                ((2 * s : ℕ) : ℝ) * (L / Real.sqrt (p : ℝ)) / Asplit p) := hsplit
          _ ≤ 2 * (((2 * s : ℕ) : ℝ) + eps * ((2 * s : ℕ) : ℝ)) := by
            nlinarith [htermSplit p hp hmp hpU hpN]
          _ = (4 + 4 * eps) * (s : ℝ) := by norm_num; ring
  · have hM_le : M ≤ 4 * s := by omega
    have hM_real : (M : ℝ) ≤ (4 * s : ℕ) := by exact_mod_cast hM_le
    have hs_nonneg : 0 ≤ (s : ℝ) := by positivity
    calc
      (M : ℝ) ≤ (4 * s : ℕ) := hM_real
      _ = 4 * (s : ℝ) := by norm_num
      _ ≤ (4 + 4 * eps) * (s : ℝ) := by nlinarith

end MainDichotomy
end GaussianChain
