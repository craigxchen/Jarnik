import Mathlib.NumberTheory.LegendreSymbol.Basic
import Mathlib.NumberTheory.SumTwoSquares
import Mathlib.NumberTheory.Zsqrtd.GaussianInt

namespace GaussianChain
namespace PrimeDescent

open scoped ComplexConjugate

/-- If a Gaussian integer has rational-prime norm, then it is irreducible in `ℤ[i]`. -/
theorem irreducible_of_norm_natAbs_eq_prime {π : GaussianInt} {p : ℕ} [Fact p.Prime]
    (hπnorm : π.norm.natAbs = p) : Irreducible π := by
  rw [irreducible_iff]
  constructor
  · intro hunit
    have hnorm1 : π.norm.natAbs = 1 := Zsqrtd.norm_eq_one_iff.mpr hunit
    have hp1 : p = 1 := by rw [← hπnorm, hnorm1]
    exact (Fact.out : Nat.Prime p).ne_one hp1
  · intro a b hmul
    have hprod : a.norm.natAbs * b.norm.natAbs = p := by
      rw [← hπnorm, hmul, Zsqrtd.norm_mul, Int.natAbs_mul]
    have hdiva : a.norm.natAbs ∣ p := ⟨b.norm.natAbs, hprod.symm⟩
    rcases (Fact.out : Nat.Prime p).eq_one_or_self_of_dvd _ hdiva with ha | ha
    · exact Or.inl (Zsqrtd.norm_eq_one_iff.mp ha)
    · right
      have hb : b.norm.natAbs = 1 := by
        rw [ha] at hprod
        exact Nat.eq_of_mul_eq_mul_left (Fact.out : Nat.Prime p).pos (by simpa using hprod)
      exact Zsqrtd.norm_eq_one_iff.mp hb

/-- Extract the integer norm from a factorization `π * conj π = p`. -/
theorem norm_eq_int_of_mul_star_eq_nat {π : GaussianInt} {p : ℕ}
    (hπ : π * star π = (p : GaussianInt)) :
    π.norm = (p : ℤ) := by
  have hcast : ((π.norm : ℤ) : GaussianInt) = (p : GaussianInt) := by
    rw [Zsqrtd.norm_eq_mul_conj]
    exact hπ
  simpa using congrArg Zsqrtd.re hcast

/-- The natural absolute norm form of `norm_eq_int_of_mul_star_eq_nat`. -/
theorem norm_natAbs_eq_of_mul_star_eq_nat {π : GaussianInt} {p : ℕ}
    (hπ : π * star π = (p : GaussianInt)) :
    π.norm.natAbs = p := by
  rw [norm_eq_int_of_mul_star_eq_nat hπ]
  simp

/-- A Gaussian factor whose product with its conjugate is a rational prime is prime in `ℤ[i]`. -/
theorem prime_of_mul_star_eq_nat_prime {π : GaussianInt} {p : ℕ} [Fact p.Prime]
    (hπ : π * star π = (p : GaussianInt)) :
    Prime π := by
  exact irreducible_iff_prime.mp
    (irreducible_of_norm_natAbs_eq_prime (norm_natAbs_eq_of_mul_star_eq_nat hπ))

/-- Fermat's two-square theorem, packaged as a Gaussian factor of a rational prime. -/
theorem exists_gaussian_factor_mul_star_eq_nat_prime_of_mod_four_ne_three
    {p : ℕ} [Fact p.Prime] (hp : p % 4 ≠ 3) :
    ∃ π : GaussianInt, π * star π = (p : GaussianInt) := by
  obtain ⟨a, b, hab⟩ := Nat.Prime.sq_add_sq (p := p) hp
  refine ⟨(⟨(a : ℤ), (b : ℤ)⟩ : GaussianInt), ?_⟩
  ext
  · simp only [Int.reduceNeg, Zsqrtd.star_mk, Zsqrtd.re_mul, neg_mul, one_mul, mul_neg,
      neg_neg, Zsqrtd.re_natCast]
    have habz : ((a : ℤ) ^ 2 + (b : ℤ) ^ 2) = (p : ℤ) := by exact_mod_cast hab
    simpa [pow_two] using habz
  · simp only [Int.reduceNeg, Zsqrtd.star_mk, Zsqrtd.im_mul, mul_neg, Zsqrtd.im_natCast]
    ring

/-- The `p ≡ 1 mod 4` split-prime version of the Gaussian factor package. -/
theorem exists_gaussian_factor_mul_star_eq_nat_prime_of_mod_four_eq_one
    {p : ℕ} [Fact p.Prime] (hp : p % 4 = 1) :
    ∃ π : GaussianInt, π * star π = (p : GaussianInt) :=
  exists_gaussian_factor_mul_star_eq_nat_prime_of_mod_four_ne_three (by omega)

/-- Divisibility can be conjugated in `ℤ[i]`. -/
theorem star_dvd_of_dvd_star {a z : GaussianInt} (h : a ∣ star z) :
    star a ∣ z := by
  rcases h with ⟨c, hc⟩
  refine ⟨star c, ?_⟩
  calc
    z = star (star z) := by simp
    _ = star (a * c) := by rw [hc]
    _ = star a * star c := by simp

/-- If a split Gaussian factor above `p` is fixed and `p` divides the norm of `z`, then either
that factor or its conjugate divides `z`. -/
theorem split_factor_dvd_or_conj_dvd_of_nat_prime_dvd_norm
    {π z : GaussianInt} {p : ℕ} [Fact p.Prime]
    (hπ : π * star π = (p : GaussianInt))
    (hpz : (p : ℤ) ∣ z.norm) :
    π ∣ z ∨ star π ∣ z := by
  have hp_dvd_normCast : (p : GaussianInt) ∣ (z.norm : GaussianInt) := by
    change ((p : ℤ) : GaussianInt) ∣ ((z.norm : ℤ) : GaussianInt)
    exact (Zsqrtd.intCast_dvd_intCast (d := -1) (p : ℤ) z.norm).2 hpz
  have hπ_dvd_normCast : π ∣ (z.norm : GaussianInt) := by
    rw [← hπ] at hp_dvd_normCast
    exact (dvd_mul_right π (star π)).trans hp_dvd_normCast
  have hπ_dvd_mul : π ∣ z * star z := by
    simpa [Zsqrtd.norm_eq_mul_conj] using hπ_dvd_normCast
  rcases (prime_of_mul_star_eq_nat_prime hπ).dvd_or_dvd hπ_dvd_mul with hz | hstarz
  · exact Or.inl hz
  · exact Or.inr (star_dvd_of_dvd_star hstarz)

/-- Split-prime pointwise descent: for `p ≡ 1 mod 4`, a point whose norm is divisible by `p`
is divisible by one of the two conjugate Gaussian factors above `p`. -/
theorem exists_split_factor_dvd_or_conj_dvd_of_nat_prime_dvd_norm
    {p : ℕ} [Fact p.Prime] (hp : p % 4 = 1) {z : GaussianInt}
    (hpz : (p : ℤ) ∣ z.norm) :
    ∃ π : GaussianInt,
      π * star π = (p : GaussianInt) ∧ (π ∣ z ∨ star π ∣ z) := by
  obtain ⟨π, hπ⟩ := exists_gaussian_factor_mul_star_eq_nat_prime_of_mod_four_eq_one hp
  exact ⟨π, hπ, split_factor_dvd_or_conj_dvd_of_nat_prime_dvd_norm hπ hpz⟩

/-- A positive natural number divisible by `p` has a positive natural quotient, in the integer
factorization form used by split-prime descent. -/
theorem exists_pos_nat_factor_of_prime_dvd_norm
    {p N : ℕ} (hN : 0 < N) (hpN : p ∣ N) :
    ∃ N' : ℕ, 0 < N' ∧ (N : ℤ) = (p : ℤ) * (N' : ℤ) := by
  rcases hpN with ⟨N', hNfac⟩
  have hN'_pos : 0 < N' := by
    by_contra hnot
    have hN'zero : N' = 0 := Nat.eq_zero_of_not_pos hnot
    rw [hN'zero, Nat.mul_zero] at hNfac
    omega
  refine ⟨N', hN'_pos, ?_⟩
  exact_mod_cast hNfac

/-- An odd rational prime is `1` or `3` modulo `4`. -/
theorem nat_prime_mod_four_eq_one_or_three_of_ne_two {p : ℕ} [Fact p.Prime]
    (hp2 : p ≠ 2) :
    p % 4 = 1 ∨ p % 4 = 3 := by
  have hp := (Fact.out : Nat.Prime p)
  have hlt : p % 4 < 4 := Nat.mod_lt p (by norm_num)
  rcases (by omega : p % 4 = 0 ∨ p % 4 = 1 ∨ p % 4 = 2 ∨ p % 4 = 3) with
    h0 | h1 | h2 | h3
  · have h2dvd : 2 ∣ p := by omega
    rcases hp.eq_one_or_self_of_dvd 2 h2dvd with htwo_one | htwo_p
    · norm_num at htwo_one
    · exact False.elim (hp2 htwo_p.symm)
  · exact Or.inl h1
  · have h2dvd : 2 ∣ p := by omega
    rcases hp.eq_one_or_self_of_dvd 2 h2dvd with htwo_one | htwo_p
    · norm_num at htwo_one
    · exact False.elim (hp2 htwo_p.symm)
  · exact Or.inr h3

/-- Number of indexed Gaussian integers divisible by a fixed Gaussian factor. -/
noncomputable def divisibleIndexCount {n : ℕ} (π : GaussianInt)
    (z : Fin n → GaussianInt) : ℕ := by
  classical
  exact (Finset.univ.filter fun i : Fin n => π ∣ z i).card

/-- The subtype of indices whose point is divisible by a fixed Gaussian factor. -/
def DivisibleIndexSubtype {n : ℕ} (π : GaussianInt) (z : Fin n → GaussianInt) : Type :=
  {i : Fin n // π ∣ z i}

noncomputable instance divisibleIndexSubtypeFintype {n : ℕ}
    (π : GaussianInt) (z : Fin n → GaussianInt) :
    Fintype (DivisibleIndexSubtype π z) := by
  classical
  dsimp [DivisibleIndexSubtype]
  infer_instance

/-- The subtype cardinal agrees with `divisibleIndexCount`. -/
theorem card_divisibleIndexSubtype_eq_divisibleIndexCount
    {n : ℕ} (π : GaussianInt) (z : Fin n → GaussianInt) :
    Fintype.card (DivisibleIndexSubtype π z) = divisibleIndexCount π z := by
  classical
  rw [divisibleIndexCount]
  simpa [DivisibleIndexSubtype] using
    (Nat.card_eq_finsetCard (Finset.univ.filter fun i : Fin n => π ∣ z i))

/-- Split-prime family descent pigeonhole: if `p ≡ 1 mod 4` divides every norm in a finite
family, then one of the two conjugate factors above `p` divides at least half the family. -/
theorem exists_split_factor_many_divisible
    {n p : ℕ} [Fact p.Prime] (hp : p % 4 = 1)
    (z : Fin n → GaussianInt) (hpz : ∀ i, (p : ℤ) ∣ (z i).norm) :
    ∃ π : GaussianInt, π * star π = (p : GaussianInt) ∧
      (n ≤ 2 * divisibleIndexCount π z ∨
        n ≤ 2 * divisibleIndexCount (star π) z) := by
  classical
  obtain ⟨π, hπ⟩ := exists_gaussian_factor_mul_star_eq_nat_prime_of_mod_four_eq_one hp
  let A : Finset (Fin n) := Finset.univ.filter fun i : Fin n => π ∣ z i
  let B : Finset (Fin n) := Finset.univ.filter fun i : Fin n => star π ∣ z i
  have hA : divisibleIndexCount π z = A.card := by simp [divisibleIndexCount, A]
  have hB : divisibleIndexCount (star π) z = B.card := by simp [divisibleIndexCount, B]
  have hunion : Finset.univ ⊆ A ∪ B := by
    intro i _hi
    have hsplit := split_factor_dvd_or_conj_dvd_of_nat_prime_dvd_norm
      (π := π) (z := z i) (p := p) hπ (hpz i)
    rcases hsplit with hdiv | hdiv
    · exact Finset.mem_union.mpr (Or.inl (by simp [A, hdiv]))
    · exact Finset.mem_union.mpr (Or.inr (by simp [B, hdiv]))
  have hcard_union : n ≤ (A ∪ B).card := by
    simpa using Finset.card_le_card hunion
  have hcard_sum : (A ∪ B).card ≤ A.card + B.card := Finset.card_union_le A B
  have hsum : n ≤ A.card + B.card := hcard_union.trans hcard_sum
  have hhalf : n ≤ 2 * A.card ∨ n ≤ 2 * B.card := by omega
  exact ⟨π, hπ, by simpa [hA, hB] using hhalf⟩

/-- Split-prime family descent pigeonhole, stated using the actual divisible-index subtypes. -/
theorem exists_split_factor_many_divisible_subtype
    {n p : ℕ} [Fact p.Prime] (hp : p % 4 = 1)
    (z : Fin n → GaussianInt) (hpz : ∀ i, (p : ℤ) ∣ (z i).norm) :
    ∃ π : GaussianInt, π * star π = (p : GaussianInt) ∧
      (n ≤ 2 * Fintype.card (DivisibleIndexSubtype π z) ∨
        n ≤ 2 * Fintype.card (DivisibleIndexSubtype (star π) z)) := by
  obtain ⟨π, hπ, hmany⟩ := exists_split_factor_many_divisible hp z hpz
  refine ⟨π, hπ, ?_⟩
  simpa [card_divisibleIndexSubtype_eq_divisibleIndexCount] using hmany

/-- Multiplying by a Gaussian factor above `p` scales the norm by `p`. -/
theorem norm_mul_of_mul_star_eq_nat
    {π w : GaussianInt} {p : ℕ} (hπ : π * star π = (p : GaussianInt)) :
    (π * w).norm = (p : ℤ) * w.norm := by
  rw [Zsqrtd.norm_mul, norm_eq_int_of_mul_star_eq_nat hπ]

/-- If `z` is divisible by a Gaussian factor above `p`, the quotient has norm scaled down by
`p`, expressed without choosing a division operation. -/
theorem exists_quotient_norm_of_dvd_split_factor
    {π z : GaussianInt} {p : ℕ} (hπ : π * star π = (p : GaussianInt))
    (hdiv : π ∣ z) :
    ∃ w : GaussianInt, z = π * w ∧ z.norm = (p : ℤ) * w.norm := by
  rcases hdiv with ⟨w, rfl⟩
  exact ⟨w, rfl, norm_mul_of_mul_star_eq_nat hπ⟩

/-- Quotient family over the subtype of indices divisible by a fixed split Gaussian factor. -/
theorem exists_split_quotient_family_on_divisible_subtype
    {n p : ℕ} (hπ : π * star π = (p : GaussianInt))
    (z : Fin n → GaussianInt) :
    ∃ w : DivisibleIndexSubtype π z → GaussianInt,
      ∀ i, z i.1 = π * w i ∧ (z i.1).norm = (p : ℤ) * (w i).norm := by
  classical
  have hq : ∀ i : DivisibleIndexSubtype π z,
      ∃ w : GaussianInt, z i.1 = π * w ∧ (z i.1).norm = (p : ℤ) * w.norm := by
    intro i
    exact exists_quotient_norm_of_dvd_split_factor hπ i.2
  choose w hw using hq
  exact ⟨w, hw⟩

/-- Split-factor quotient family with an explicitly named quotient norm. -/
theorem exists_split_quotient_family_norm_eq_on_divisible_subtype
    {n p : ℕ} [Fact p.Prime] {π : GaussianInt}
    (hπ : π * star π = (p : GaussianInt))
    {N N' : ℤ} (hNfactor : N = (p : ℤ) * N')
    (z : Fin n → GaussianInt) (hz : ∀ i, (z i).norm = N) :
    ∃ w : DivisibleIndexSubtype π z → GaussianInt,
      ∀ i, z i.1 = π * w i ∧ (w i).norm = N' := by
  obtain ⟨w, hw⟩ := exists_split_quotient_family_on_divisible_subtype hπ z
  refine ⟨w, ?_⟩
  intro i
  rcases hw i with ⟨hzw, hnorm⟩
  refine ⟨hzw, ?_⟩
  have hp_ne : (p : ℤ) ≠ 0 := by
    exact_mod_cast (Fact.out : Nat.Prime p).ne_zero
  have heq : (p : ℤ) * N' = (p : ℤ) * (w i).norm := by
    rw [← hNfactor, ← hz i.1]
    exact hnorm
  exact mul_left_cancel₀ hp_ne heq.symm

/-- Split-prime descent package: when `N = p * N'`, at least half of a finite family on norm
`N` descends through one Gaussian factor above `p` to norm `N'`. -/
theorem exists_split_descended_subfamily_norm_eq
    {n p : ℕ} [Fact p.Prime] (hp : p % 4 = 1)
    {N N' : ℤ} (hNfactor : N = (p : ℤ) * N')
    (z : Fin n → GaussianInt) (hz : ∀ i, (z i).norm = N) :
    ∃ ρ : GaussianInt,
      ρ * star ρ = (p : GaussianInt) ∧
      n ≤ 2 * Fintype.card (DivisibleIndexSubtype ρ z) ∧
      ∃ w : DivisibleIndexSubtype ρ z → GaussianInt,
        ∀ i, z i.1 = ρ * w i ∧ (w i).norm = N' := by
  have hpz : ∀ i, (p : ℤ) ∣ (z i).norm := by
    intro i
    rw [hz i, hNfactor]
    exact dvd_mul_right (p : ℤ) N'
  obtain ⟨π, hπ, hmany⟩ := exists_split_factor_many_divisible_subtype hp z hpz
  rcases hmany with hmany | hmany
  · obtain ⟨w, hw⟩ :=
      exists_split_quotient_family_norm_eq_on_divisible_subtype hπ hNfactor z hz
    exact ⟨π, hπ, hmany, w, hw⟩
  · have hπstar : star π * star (star π) = (p : GaussianInt) := by
      simpa [mul_comm] using hπ
    obtain ⟨w, hw⟩ :=
      exists_split_quotient_family_norm_eq_on_divisible_subtype hπstar hNfactor z hz
    exact ⟨star π, hπstar, hmany, w, hw⟩

/-- Ordered split-prime descent package.

The split-prime pigeonhole gives a divisible subtype of size `Md`. This theorem enumerates that
subtype in increasing order of the original `Fin M` indices, turns the quotient family into a
natural-indexed sequence, and records the facts needed by the sliding-window descent bound. -/
theorem exists_split_descended_ordered_subsequence_norm_eq
    {M p : ℕ} [Fact p.Prime] (hp : p % 4 = 1)
    {N N' : ℤ} (hNfactor : N = (p : ℤ) * N')
    (z : Fin M → GaussianInt) (hz : ∀ i, (z i).norm = N) (hM : 0 < M) :
    ∃ (Md : ℕ) (ρ : GaussianInt) (idx : ℕ → Fin M) (w : ℕ → GaussianInt),
      ρ * star ρ = (p : GaussianInt) ∧
      M ≤ 2 * Md ∧
      (∀ i, i < Md → z (idx i) = ρ * w i) ∧
      (∀ i, (w i).norm = N') ∧
      (∀ i j, i ≤ j → j < Md → (idx i : ℕ) ≤ (idx j : ℕ)) ∧
      (Function.Injective z →
        ∀ s j, j < Md - 2 * s →
          Function.Injective fun i : Fin (2 * s + 1) => w (j + (i : ℕ))) := by
  classical
  obtain ⟨ρ, hρ, hmany, wD, hwD⟩ :=
    exists_split_descended_subfamily_norm_eq hp hNfactor z hz
  let D : Finset (Fin M) := Finset.univ.filter fun i : Fin M => ρ ∣ z i
  let Md : ℕ := Fintype.card (DivisibleIndexSubtype ρ z)
  have hMd_pos : 0 < Md := by omega
  have hDcard : D.card = Md := by
    simp [D, Md, card_divisibleIndexSubtype_eq_divisibleIndexCount, divisibleIndexCount]
  let e : Fin Md ↪o Fin M := D.orderEmbOfFin hDcard
  have he_mem : ∀ i : Fin Md, e i ∈ D := fun i => D.orderEmbOfFin_mem hDcard i
  let toD : Fin Md → DivisibleIndexSubtype ρ z := fun i =>
    ⟨e i, (Finset.mem_filter.mp (he_mem i)).2⟩
  let idx : ℕ → Fin M := fun i => if hi : i < Md then e ⟨i, hi⟩ else ⟨0, hM⟩
  let w : ℕ → GaussianInt :=
    fun i => if hi : i < Md then wD (toD ⟨i, hi⟩) else wD (toD ⟨0, hMd_pos⟩)
  refine ⟨Md, ρ, idx, w, hρ, hmany, ?_, ?_, ?_, ?_⟩
  · intro i hi
    have hidx : idx i = e ⟨i, hi⟩ := by simp [idx, hi]
    have hw : w i = wD (toD ⟨i, hi⟩) := by simp [w, hi]
    rw [hidx, hw]
    exact (hwD (toD ⟨i, hi⟩)).1
  · intro i
    by_cases hi : i < Md
    · have hw : w i = wD (toD ⟨i, hi⟩) := by simp [w, hi]
      rw [hw]
      exact (hwD (toD ⟨i, hi⟩)).2
    · have hw : w i = wD (toD ⟨0, hMd_pos⟩) := by simp [w, hi]
      rw [hw]
      exact (hwD (toD ⟨0, hMd_pos⟩)).2
  · intro i j hij hj
    have hi : i < Md := lt_of_le_of_lt hij hj
    have hidx_i : idx i = e ⟨i, hi⟩ := by simp [idx, hi]
    have hidx_j : idx j = e ⟨j, hj⟩ := by simp [idx, hj]
    rw [hidx_i, hidx_j]
    exact e.monotone (by exact hij)
  · intro hz_inj s j hj
    have hsub_inj : Function.Injective fun i : DivisibleIndexSubtype ρ z => z i.1 := by
      intro i k hik
      exact Subtype.ext (hz_inj hik)
    have hwD_inj : Function.Injective wD := by
      intro a b hab
      apply hsub_inj
      change z a.1 = z b.1
      rw [(hwD a).1, (hwD b).1, hab]
    intro i k hik
    apply Fin.ext
    have hi_lt : j + (i : ℕ) < Md := by
      have hi_bound : (i : ℕ) < 2 * s + 1 := i.isLt
      omega
    have hk_lt : j + (k : ℕ) < Md := by
      have hk_bound : (k : ℕ) < 2 * s + 1 := k.isLt
      omega
    have hwi : w (j + (i : ℕ)) = wD (toD ⟨j + (i : ℕ), hi_lt⟩) := by
      simp [w, hi_lt]
    have hwk : w (j + (k : ℕ)) = wD (toD ⟨j + (k : ℕ), hk_lt⟩) := by
      simp [w, hk_lt]
    have htoD :
        toD ⟨j + (i : ℕ), hi_lt⟩ = toD ⟨j + (k : ℕ), hk_lt⟩ :=
      hwD_inj (by simpa [hwi, hwk] using hik)
    have heq :
        (⟨j + (i : ℕ), hi_lt⟩ : Fin Md) =
          ⟨j + (k : ℕ), hk_lt⟩ := by
      apply e.injective
      simpa [toD] using congrArg Subtype.val htoD
    exact Nat.add_left_cancel (Fin.ext_iff.mp heq)

/-- Multiplication by the rational Gaussian integer `p` scales the norm by `p^2`. -/
theorem norm_mul_natCast (p : ℕ) (w : GaussianInt) :
    ((p : GaussianInt) * w).norm = (p : ℤ) ^ 2 * w.norm := by
  rw [Zsqrtd.norm_mul]
  simp [pow_two]

/-- If `z` is divisible by the rational Gaussian integer `p`, the quotient has norm scaled down
by `p^2`, expressed without choosing a division operation. -/
theorem exists_quotient_norm_of_dvd_natCast
    {p : ℕ} {z : GaussianInt} (hdiv : (p : GaussianInt) ∣ z) :
    ∃ w : GaussianInt, z = (p : GaussianInt) * w ∧ z.norm = (p : ℤ) ^ 2 * w.norm := by
  rcases hdiv with ⟨w, rfl⟩
  exact ⟨w, rfl, norm_mul_natCast p w⟩

/-- Inert-prime coordinate divisibility: if `p ≡ 3 mod 4` divides `x^2 + y^2`,
then it divides both `x` and `y`. -/
theorem inert_prime_dvd_sq_add_sq_int {p : ℕ} [Fact p.Prime] (hp3 : p % 4 = 3)
    {x y : ℤ} (h : (p : ℤ) ∣ x ^ 2 + y ^ 2) :
    (p : ℤ) ∣ x ∧ (p : ℤ) ∣ y := by
  have hsum_cast : ((x ^ 2 + y ^ 2 : ℤ) : ZMod p) = 0 :=
    (ZMod.intCast_zmod_eq_zero_iff_dvd _ p).2 h
  have hsum : (x : ZMod p) ^ 2 + (y : ZMod p) ^ 2 = 0 := by
    simpa [Int.cast_add, Int.cast_pow] using hsum_cast
  have hy0 : (y : ZMod p) = 0 := by
    by_contra hy
    have hxeq : (x : ZMod p) ^ 2 = - (y : ZMod p) ^ 2 := by
      rw [eq_neg_iff_add_eq_zero]
      exact hsum
    exact (ZMod.mod_four_ne_three_of_sq_eq_neg_sq' hy hxeq) hp3
  have hx0 : (x : ZMod p) = 0 := by
    have hx_sq : (x : ZMod p) ^ 2 = 0 := by
      simpa [hy0] using hsum
    exact sq_eq_zero_iff.mp hx_sq
  exact ⟨(ZMod.intCast_zmod_eq_zero_iff_dvd x p).1 hx0,
    (ZMod.intCast_zmod_eq_zero_iff_dvd y p).1 hy0⟩

/-- If an inert rational prime divides the Gaussian norm, then it divides the Gaussian integer. -/
theorem inert_prime_dvd_gaussian_of_dvd_norm {p : ℕ} [Fact p.Prime]
    (hp3 : p % 4 = 3) {z : GaussianInt} (h : (p : ℤ) ∣ z.norm) :
    (p : GaussianInt) ∣ z := by
  have hsum : (p : ℤ) ∣ z.re ^ 2 + z.im ^ 2 := by
    convert h using 1
    simp [Zsqrtd.norm]
    ring
  have hcoords := inert_prime_dvd_sq_add_sq_int (p := p) hp3 hsum
  change ((p : ℤ) : GaussianInt) ∣ z
  exact (Zsqrtd.intCast_dvd (d := -1) (p : ℤ) z).2 hcoords

/-- A norm-equality variant of inert-prime Gaussian divisibility. -/
theorem inert_prime_dvd_gaussian_of_norm_eq {p : ℕ} [Fact p.Prime]
    (hp3 : p % 4 = 3) {z : GaussianInt} {N : ℤ}
    (hzn : z.norm = N) (hpN : (p : ℤ) ∣ N) :
    (p : GaussianInt) ∣ z := by
  exact inert_prime_dvd_gaussian_of_dvd_norm (p := p) hp3 (by simpa [hzn] using hpN)

/-- If an inert prime divides a positive natural norm represented by a Gaussian integer, then
it divides that norm to even order at least two. This supplies the natural quotient needed for
inert descent. -/
theorem exists_pos_nat_factor_of_inert_prime_dvd_norm
    {p N : ℕ} [Fact p.Prime] (hp3 : p % 4 = 3) (hN : 0 < N)
    {z : GaussianInt} (hz : z.norm = (N : ℤ)) (hpN : p ∣ N) :
    ∃ N' : ℕ, 0 < N' ∧ (N : ℤ) = (p : ℤ) ^ 2 * (N' : ℤ) := by
  have hpN_int : (p : ℤ) ∣ (N : ℤ) := by exact_mod_cast hpN
  have hdiv : (p : GaussianInt) ∣ z :=
    inert_prime_dvd_gaussian_of_norm_eq hp3 hz hpN_int
  obtain ⟨w, hzw, hnorm⟩ := exists_quotient_norm_of_dvd_natCast hdiv
  let N' : ℕ := w.norm.natAbs
  have hfactor : (N : ℤ) = (p : ℤ) ^ 2 * (N' : ℤ) := by
    rw [← hz, hnorm]
    simp [N', Int.natAbs_of_nonneg (GaussianInt.norm_nonneg w)]
  have hfactor_nat : N = p ^ 2 * N' := by
    exact_mod_cast hfactor
  have hN'_pos : 0 < N' := by
    by_contra hnot
    have hN'zero : N' = 0 := Nat.eq_zero_of_not_pos hnot
    rw [hN'zero, Nat.mul_zero] at hfactor_nat
    omega
  exact ⟨N', hN'_pos, hfactor⟩

/-- Family form of inert descent: if an inert rational prime divides the common norm, every
point in the family has a quotient by that rational Gaussian integer. -/
theorem exists_inert_quotient_family_of_norm_eq
    {ι : Type*} {p : ℕ} [Fact p.Prime] (hp3 : p % 4 = 3)
    {N : ℤ} (hpN : (p : ℤ) ∣ N)
    (z : ι → GaussianInt) (hz : ∀ i, (z i).norm = N) :
    ∃ w : ι → GaussianInt,
      ∀ i, z i = (p : GaussianInt) * w i ∧
        (z i).norm = (p : ℤ) ^ 2 * (w i).norm := by
  classical
  have hq : ∀ i, ∃ w : GaussianInt,
      z i = (p : GaussianInt) * w ∧ (z i).norm = (p : ℤ) ^ 2 * w.norm := by
    intro i
    exact exists_quotient_norm_of_dvd_natCast
      (inert_prime_dvd_gaussian_of_norm_eq hp3 (hz i) hpN)
  choose w hw using hq
  exact ⟨w, hw⟩

/-- Inert descent with an explicitly named quotient norm. -/
theorem exists_inert_quotient_family_norm_eq_of_common_norm_factor
    {ι : Type*} {p : ℕ} [Fact p.Prime] (hp3 : p % 4 = 3)
    {N N' : ℤ} (hNfactor : N = (p : ℤ) ^ 2 * N')
    (z : ι → GaussianInt) (hz : ∀ i, (z i).norm = N) :
    ∃ w : ι → GaussianInt,
      ∀ i, z i = (p : GaussianInt) * w i ∧ (w i).norm = N' := by
  have hpN : (p : ℤ) ∣ N := by
    rw [hNfactor]
    have hp_dvd_sq : (p : ℤ) ∣ (p : ℤ) ^ 2 := by
      rw [pow_two]
      exact dvd_mul_right (p : ℤ) (p : ℤ)
    exact hp_dvd_sq.trans (dvd_mul_right ((p : ℤ) ^ 2) N')
  obtain ⟨w, hw⟩ := exists_inert_quotient_family_of_norm_eq hp3 hpN z hz
  refine ⟨w, ?_⟩
  intro i
  rcases hw i with ⟨hzw, hnorm⟩
  refine ⟨hzw, ?_⟩
  have hp2_ne : (p : ℤ) ^ 2 ≠ 0 := by
    exact pow_ne_zero _ (by exact_mod_cast (Fact.out : Nat.Prime p).ne_zero)
  have heq : (p : ℤ) ^ 2 * N' = (p : ℤ) ^ 2 * (w i).norm := by
    rw [← hNfactor, ← hz i]
    exact hnorm
  exact mul_left_cancel₀ hp2_ne heq.symm

/-- If a family is obtained by multiplying a quotient family by a fixed factor and the original
family is injective, then the quotient family is injective. -/
theorem injective_quotient_of_mul_left_injective
    {ι : Type*} {a : GaussianInt} {z w : ι → GaussianInt}
    (hzfac : ∀ i, z i = a * w i)
    (hz : Function.Injective z) :
    Function.Injective w := by
  intro i j hij
  apply hz
  rw [hzfac i, hzfac j, hij]

end PrimeDescent
end GaussianChain
