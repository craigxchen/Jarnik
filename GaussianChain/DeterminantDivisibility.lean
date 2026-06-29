import GaussianChain.ModPrimeCollision
import GaussianChain.RamanaDeterminant

namespace GaussianChain
namespace DeterminantDivisibility

open ModPrimeCollision
open RamanaDeterminant
open scoped ComplexConjugate

/-- Number of unordered pairs `i < j` whose Gaussian residues modulo `p` collide. -/
def pairCollisionCount {n : ℕ} (p : ℕ) (z : Fin n → GaussianInt) : ℕ :=
  ∑ i : Fin n,
    ((Finset.Ioi i).filter fun j => gaussianResidue p (z i) = gaussianResidue p (z j)).card

section OrderedUnordered

variable {n : ℕ} {β : Type*} [DecidableEq β]

/-- Number of unordered colliding pairs `i < j` for an arbitrary finite sequence. -/
def strictPairCollisionCount (f : Fin n → β) : ℕ :=
  ∑ i : Fin n, ((Finset.Ioi i).filter fun j => f i = f j).card

/-- Colliding pairs with the first index strictly below the second. -/
def UpperCollision (f : Fin n → β) : Type _ :=
  {q : Fin n × Fin n // q.1 < q.2 ∧ f q.1 = f q.2}

instance upperCollisionFintype (f : Fin n → β) : Fintype (UpperCollision f) := by
  dsimp [UpperCollision]
  infer_instance

/-- Upper-triangular colliding pairs are the same as choosing a row and then a colliding
entry in its strict upper tail. -/
noncomputable def upperCollisionEquivSigma (f : Fin n → β) :
    UpperCollision f ≃
      Sigma fun i : Fin n => {j : Fin n // j ∈ (Finset.Ioi i).filter fun j => f i = f j} where
  toFun u := ⟨u.1.1, ⟨u.1.2, by
    rw [Finset.mem_filter]
    exact ⟨Finset.mem_Ioi.mpr u.2.1, u.2.2⟩⟩⟩
  invFun s := ⟨(s.1, s.2.1), by
    have hs := Finset.mem_filter.mp s.2.2
    exact ⟨Finset.mem_Ioi.mp hs.1, hs.2⟩⟩
  left_inv u := by
    rcases u with ⟨⟨i, j⟩, _hlt, _hij⟩
    rfl
  right_inv s := by
    rcases s with ⟨i, j, hj⟩
    rfl

/-- The subtype-cardinality form of `strictPairCollisionCount`. -/
theorem card_upperCollision_eq_strictPairCollisionCount (f : Fin n → β) :
    Fintype.card (UpperCollision f) = strictPairCollisionCount f := by
  classical
  have hrow (i : Fin n) :
      Fintype.card {j : Fin n // j ∈ (Finset.Ioi i).filter fun j => f i = f j} =
        ((Finset.Ioi i).filter fun j => f i = f j).card := by
    simpa using (Nat.card_eq_finsetCard ((Finset.Ioi i).filter fun j => f i = f j))
  rw [strictPairCollisionCount, Fintype.card_congr (upperCollisionEquivSigma f),
    Fintype.card_sigma]
  exact Finset.sum_congr rfl (fun i _ => hrow i)

/-- Ordered colliding pairs split into the diagonal plus two orientations of the strict
upper-triangular colliding pairs. -/
noncomputable def orderedCollisionEquivDiagUpper (f : Fin n → β) :
    {p : Fin n × Fin n // f p.1 = f p.2} ≃
      (Fin n) ⊕ (UpperCollision f ⊕ UpperCollision f) where
  toFun p :=
    if hdiag : p.1.1 = p.1.2 then
      Sum.inl p.1.1
    else if hlt : p.1.1 < p.1.2 then
      Sum.inr (Sum.inl ⟨p.1, ⟨hlt, p.2⟩⟩)
    else
      have hgt : p.1.2 < p.1.1 := by
        exact lt_of_le_of_ne (le_of_not_gt hlt) (fun h => hdiag h.symm)
      Sum.inr (Sum.inr ⟨(p.1.2, p.1.1), ⟨hgt, p.2.symm⟩⟩)
  invFun s :=
    match s with
    | Sum.inl i => ⟨(i, i), rfl⟩
    | Sum.inr (Sum.inl u) => ⟨u.1, u.2.2⟩
    | Sum.inr (Sum.inr u) => ⟨(u.1.2, u.1.1), u.2.2.symm⟩
  left_inv p := by
    rcases p with ⟨⟨i, j⟩, _hij⟩
    dsimp
    by_cases hdiag : i = j
    · subst j
      simp
    · by_cases hlt : i < j
      · simp [hdiag, hlt]
      · simp [hdiag, hlt]
  right_inv s := by
    rcases s with i | u | u
    · simp
    · rcases u with ⟨⟨i, j⟩, hlt, _hij⟩
      have hdiag : ¬ i = j := ne_of_lt hlt
      simp [hdiag, hlt]
    · rcases u with ⟨⟨i, j⟩, hlt, _hij⟩
      have hdiag : ¬ j = i := fun h => (ne_of_lt hlt) h.symm
      have hnlt : ¬ j < i := not_lt_of_gt hlt
      simp [hdiag, hnlt]

/-- Exact relation between ordered collisions and unordered strict-pair collisions. -/
theorem orderedCollisionCount_eq_card_add_two_mul_strictPairCollisionCount
    (f : Fin n → β) :
    orderedCollisionCount f = n + 2 * strictPairCollisionCount f := by
  classical
  rw [orderedCollisionCount, Fintype.card_congr (orderedCollisionEquivDiagUpper f),
    Fintype.card_sum, Fintype.card_sum, Fintype.card_fin,
    card_upperCollision_eq_strictPairCollisionCount]
  ring

/-- The unordered strict-pair count is half of the off-diagonal ordered collision count. -/
theorem strictPairCollisionCount_eq_orderedCollisionCount_sub_card_div_two
    (f : Fin n → β) :
    (strictPairCollisionCount f : ℝ) =
      ((orderedCollisionCount f : ℝ) - n) / 2 := by
  have hnat := orderedCollisionCount_eq_card_add_two_mul_strictPairCollisionCount f
  have hreal :
      (orderedCollisionCount f : ℝ) =
        (n : ℝ) + 2 * (strictPairCollisionCount f : ℝ) := by
    exact_mod_cast hnat
  rw [hreal]
  ring

/-- Pigeonhole lower bound for strict unordered collision pairs. -/
theorem strictPairCollisionCount_lower_bound
    [Fintype β] (f : Fin n → β) {N : ℕ}
    (hβN : Fintype.card β ≤ N) (hN : 0 < N) :
    ((n : ℝ) ^ 2) / (2 * (N : ℝ)) - (n : ℝ) / 2 ≤
      (strictPairCollisionCount f : ℝ) := by
  classical
  have h := offDiagonalCollisionLowerBound (f := f) hβN hN
  rw [Fintype.card_fin] at h
  rwa [strictPairCollisionCount_eq_orderedCollisionCount_sub_card_div_two f]

end OrderedUnordered

/-- The Gaussian residue collision count is the generic strict-pair collision count of the
Gaussian residue map. -/
theorem pairCollisionCount_eq_strictPairCollisionCount
    {n : ℕ} (p : ℕ) (z : Fin n → GaussianInt) :
    pairCollisionCount p z =
      strictPairCollisionCount (fun i => gaussianResidue p (z i)) := rfl

/-- Ordered Gaussian residue collisions split into the diagonal and twice the unordered
Gaussian residue collision count. -/
theorem orderedCollisionCount_gaussianResidue_eq_card_add_two_mul_pairCollisionCount
    {n : ℕ} (p : ℕ) (z : Fin n → GaussianInt) :
    orderedCollisionCount (fun i => gaussianResidue p (z i)) =
      n + 2 * pairCollisionCount p z := by
  rw [orderedCollisionCount_eq_card_add_two_mul_strictPairCollisionCount,
    ← pairCollisionCount_eq_strictPairCollisionCount]

/-- For points on a Gaussian norm circle, the subtype-valued residue map has the same strict
collision count as the underlying Gaussian residue map. -/
theorem strictPairCollisionCount_gaussianCircleResidueMap_eq_pairCollisionCount
    {n : ℕ} (p : ℕ) [Fact p.Prime] (N : ℤ) (z : Fin n → GaussianInt)
    (hz : ∀ i, (z i).norm = N) :
    strictPairCollisionCount (gaussianCircleResidueMap p N z hz) = pairCollisionCount p z := by
  unfold strictPairCollisionCount pairCollisionCount gaussianCircleResidueMap
  simp

/-- Pigeonhole lower bound for the unordered Gaussian residue collisions on a fixed norm circle.
This is the exponent lower bound used by the determinant divisibility step. -/
theorem pairCollisionCount_lower_bound_on_circle
    {n : ℕ} (p : ℕ) [Fact p.Prime] (N : ℤ) (z : Fin n → GaussianInt)
    (hz : ∀ i, (z i).norm = N) :
    ((n : ℝ) ^ 2) / (4 * (p : ℝ)) - (n : ℝ) / 2 ≤
      (pairCollisionCount p z : ℝ) := by
  have h := gaussianCircleResidue_offDiagonalCollisionLowerBound
    (α := Fin n) p N z hz
  rw [Fintype.card_fin] at h
  rw [← strictPairCollisionCount_eq_orderedCollisionCount_sub_card_div_two
    (gaussianCircleResidueMap p N z hz)] at h
  simpa [strictPairCollisionCount, pairCollisionCount, gaussianCircleResidueMap] using h

/-- If two Gaussian integers have the same residue pair modulo `p`, their difference is divisible
by the rational prime `p` in `ℤ[i]`. -/
theorem prime_dvd_difference_of_gaussianResidue_eq
    (p : ℕ) {z w : GaussianInt} (h : gaussianResidue p z = gaussianResidue p w) :
    (p : GaussianInt) ∣ w - z := by
  have hre : (z.re : ZMod p) = (w.re : ZMod p) := congrArg Prod.fst h
  have him : (z.im : ZMod p) = (w.im : ZMod p) := congrArg Prod.snd h
  have hdvd_re : (p : ℤ) ∣ w.re - z.re := by
    refine (ZMod.intCast_zmod_eq_zero_iff_dvd (w.re - z.re) p).1 ?_
    rw [Int.cast_sub, hre, sub_self]
  have hdvd_im : (p : ℤ) ∣ w.im - z.im := by
    refine (ZMod.intCast_zmod_eq_zero_iff_dvd (w.im - z.im) p).1 ?_
    rw [Int.cast_sub, him, sub_self]
  change ((p : ℤ) : GaussianInt) ∣ w - z
  exact (Zsqrtd.intCast_dvd (d := -1) (p : ℤ) (w - z)).2 (by
    simpa using And.intro hdvd_re hdvd_im)

/-- A rational prime power is coprime to an integer if the prime itself does not divide it. -/
theorem int_gcd_prime_pow_eq_one_of_not_dvd {p e : ℕ} [Fact p.Prime] {a : ℤ}
    (h : ¬ (p : ℤ) ∣ a) :
    Int.gcd ((p : ℤ) ^ e) a = 1 := by
  have hnotnat : ¬ p ∣ a.natAbs := by
    intro hd
    have hd' : ((p : ℤ).natAbs) ∣ a.natAbs := by
      simpa [Int.natAbs_natCast] using hd
    exact h ((Int.natAbs_dvd_natAbs).mp hd')
  have hcop : Nat.Coprime (p ^ e) a.natAbs := by
    exact ((Fact.out : Nat.Prime p).coprime_iff_not_dvd.mpr hnotnat).pow_left e
  rw [Int.gcd_def, Int.natAbs_pow, Int.natAbs_natCast]
  exact hcop.gcd_eq_one

/-- If a Gaussian integer has norm prime to `p`, multiplication by it can be cancelled from
rational `p`-power divisibility in `ℤ[i]`.

The proof is just the coordinate inverse matrix: if `u = a + bi`, then linear combinations of
the coordinates of `u * d` give `(a^2 + b^2) d.re` and `(a^2 + b^2) d.im`. -/
theorem nat_prime_pow_dvd_right_of_mul_left_norm_coprime
    {p e : ℕ} [Fact p.Prime] {u d : GaussianInt}
    (hu : ¬ (p : ℤ) ∣ u.norm)
    (hdiv : ((p : GaussianInt) ^ e) ∣ u * d) :
    ((p : GaussianInt) ^ e) ∣ d := by
  have hdivInt : (((p : ℤ) : GaussianInt) ^ e) ∣ u * d := by
    simpa using hdiv
  change (((p : ℤ) : GaussianInt) ^ e) ∣ d
  rw [← Int.cast_pow] at hdivInt ⊢
  rw [Zsqrtd.intCast_dvd] at hdivInt ⊢
  rcases hdivInt with ⟨hre, him⟩
  constructor
  · have hlin : ((p : ℤ) ^ e) ∣ u.re * (u * d).re - (-u.im) * (u * d).im :=
      dvd_sub (dvd_mul_of_dvd_right hre _) (dvd_mul_of_dvd_right him _)
    have hcoord : u.re * (u * d).re - (-u.im) * (u * d).im = u.norm * d.re := by
      simp [Zsqrtd.norm, Zsqrtd.re_mul, Zsqrtd.im_mul]
      ring
    rw [hcoord] at hlin
    exact Int.dvd_of_dvd_mul_right_of_gcd_one hlin
      (int_gcd_prime_pow_eq_one_of_not_dvd (p := p) (e := e) hu)
  · have hlin : ((p : ℤ) ^ e) ∣ u.re * (u * d).im + (-u.im) * (u * d).re :=
      dvd_add (dvd_mul_of_dvd_right him _) (dvd_mul_of_dvd_right hre _)
    have hcoord : u.re * (u * d).im + (-u.im) * (u * d).re = u.norm * d.im := by
      simp [Zsqrtd.norm, Zsqrtd.re_mul, Zsqrtd.im_mul]
      ring
    rw [hcoord] at hlin
    exact Int.dvd_of_dvd_mul_right_of_gcd_one hlin
      (int_gcd_prime_pow_eq_one_of_not_dvd (p := p) (e := e) hu)

/-- Integer prime nondivisibility is preserved by powers. -/
theorem int_prime_not_dvd_pow_of_not_dvd {p k : ℕ} [Fact p.Prime] {a : ℤ}
    (h : ¬ (p : ℤ) ∣ a) :
    ¬ (p : ℤ) ∣ a ^ k := by
  intro hpow
  have hpowNat : p ∣ a.natAbs ^ k := by
    have hpowNat' : ((p : ℤ).natAbs) ∣ (a ^ k).natAbs :=
      (Int.natAbs_dvd_natAbs).mpr hpow
    simpa [Int.natAbs_natCast, Int.natAbs_pow] using hpowNat'
  have hbaseNat : p ∣ a.natAbs :=
    (Fact.out : Nat.Prime p).dvd_of_dvd_pow hpowNat
  have hbaseInt : (p : ℤ) ∣ a := by
    have hbaseNat' : ((p : ℤ).natAbs) ∣ a.natAbs := by
      simpa [Int.natAbs_natCast] using hbaseNat
    exact (Int.natAbs_dvd_natAbs).mp hbaseNat'
  exact h hbaseInt

/-- Norm divisibility induced by Gaussian divisibility. -/
theorem norm_natAbs_dvd_of_dvd {a b : GaussianInt} (h : a ∣ b) :
    a.norm.natAbs ∣ b.norm.natAbs := by
  rcases h with ⟨c, rfl⟩
  rw [Zsqrtd.norm_mul, Int.natAbs_mul]
  exact dvd_mul_right _ _

/-- Norm of a rational prime power inside `ℤ[i]`. -/
theorem norm_natAbs_natCast_pow (p e : ℕ) :
    (((p : GaussianInt) ^ e).norm).natAbs = p ^ (2 * e) := by
  change (Zsqrtd.normMonoidHom (d := -1) ((p : GaussianInt) ^ e)).natAbs = p ^ (2 * e)
  rw [map_pow]
  simp [Zsqrtd.normMonoidHom, Int.natAbs_pow, Int.natAbs_mul, Int.natAbs_natCast, pow_mul]
  ring_nf

/-- A finite-product divisibility helper: if every selected factor is divisible by `a`,
then `a` to the number of selected factors divides the whole product. -/
theorem pow_card_filter_dvd_prod
    {ι M : Type*} [CommMonoid M] (a : M) (s : Finset ι) (P : ι → Prop)
    [DecidablePred P] (f : ι → M)
    (h : ∀ x ∈ s, P x → a ∣ f x) :
    a ^ (s.filter P).card ∣ ∏ x ∈ s, f x := by
  have hfilter :
      a ^ (s.filter P).card ∣ ∏ x ∈ s.filter P, f x := by
    rw [← Finset.prod_const]
    refine Finset.prod_dvd_prod_of_dvd (fun _ => a) f ?_
    intro x hx
    exact h x (Finset.mem_of_mem_filter x hx) (Finset.mem_filter.mp hx).2
  have hfilterProd :
      (∏ x ∈ s.filter P, f x) ∣ ∏ x ∈ s, f x := by
    rw [← Finset.prod_filter_mul_prod_filter_not s P f]
    exact dvd_mul_right _ _
  exact hfilter.trans hfilterProd

/-- Every colliding unordered pair contributes one factor of `p` to the Vandermonde pair product. -/
theorem pow_pairCollisionCount_dvd_vandermondeProduct
    {n : ℕ} (p : ℕ) (z : Fin n → GaussianInt) :
    (p : GaussianInt) ^ pairCollisionCount p z ∣
      ∏ i : Fin n, ∏ j ∈ Finset.Ioi i, (z j - z i) := by
  classical
  unfold pairCollisionCount
  rw [← Finset.prod_pow_eq_pow_sum]
  refine Finset.prod_dvd_prod_of_dvd
    (fun i : Fin n =>
      (p : GaussianInt) ^
        ((Finset.Ioi i).filter fun j => gaussianResidue p (z i) = gaussianResidue p (z j)).card)
    (fun i : Fin n => ∏ j ∈ Finset.Ioi i, (z j - z i)) ?_
  intro i _hi
  exact pow_card_filter_dvd_prod (p : GaussianInt) (Finset.Ioi i)
    (fun j => gaussianResidue p (z i) = gaussianResidue p (z j))
    (fun j => z j - z i)
    (by
      intro j _hj hcollision
      exact prime_dvd_difference_of_gaussianResidue_eq p hcollision)

/-- Collision divisibility transferred through Ramana's determinant identity. This is the local
divisibility statement before cancellation of the row-scaling and norm factors. -/
theorem pow_pairCollisionCount_dvd_prod_pow_mul_det_ramana
    {s : ℕ} {N : ℤ} {z : Fin (2 * s + 1) → GaussianInt} (p : ℕ)
    (hcircle : ∀ i, z i * star (z i) = (N : GaussianInt)) :
    (p : GaussianInt) ^ pairCollisionCount p z ∣
      (∏ i, z i ^ s) * (ramanaMatrix s z).det := by
  rw [prod_pow_mul_det_ramana_eq_normPow_mul_vandermondeProduct hcircle]
  exact dvd_mul_of_dvd_right (pow_pairCollisionCount_dvd_vandermondeProduct p z) _

/-- Collision divisibility after cancelling the row-scaling product. The hypothesis `p ∤ N`
is exactly what makes every row-scaling factor invertible modulo all powers of the rational
prime `p`. -/
theorem pow_pairCollisionCount_dvd_det_ramana_of_prime_not_dvd_common_norm
    {s : ℕ} {N : ℤ} {z : Fin (2 * s + 1) → GaussianInt} (p : ℕ) [Fact p.Prime]
    (hcircle : ∀ i, z i * star (z i) = (N : GaussianInt))
    (hpN : ¬ (p : ℤ) ∣ N) :
    (p : GaussianInt) ^ pairCollisionCount p z ∣ (ramanaMatrix s z).det := by
  let u : GaussianInt := ∏ i, z i ^ s
  have hu_norm : u.norm = N ^ (s * (2 * s + 1)) := by
    dsimp [u]
    rw [norm_prod_pow, prod_norm_pow_eq_common hcircle]
  have hu_coprime : ¬ (p : ℤ) ∣ u.norm := by
    rw [hu_norm]
    exact int_prime_not_dvd_pow_of_not_dvd (p := p) hpN
  exact nat_prime_pow_dvd_right_of_mul_left_norm_coprime hu_coprime
    (pow_pairCollisionCount_dvd_prod_pow_mul_det_ramana p hcircle)

/-- Norm-valued form of the cancelled local determinant divisibility. -/
theorem pow_two_mul_pairCollisionCount_dvd_det_ramana_norm_natAbs
    {s : ℕ} {N : ℤ} {z : Fin (2 * s + 1) → GaussianInt} (p : ℕ) [Fact p.Prime]
    (hcircle : ∀ i, z i * star (z i) = (N : GaussianInt))
    (hpN : ¬ (p : ℤ) ∣ N) :
    p ^ (2 * pairCollisionCount p z) ∣ (ramanaMatrix s z).det.norm.natAbs := by
  have hdiv :=
    pow_pairCollisionCount_dvd_det_ramana_of_prime_not_dvd_common_norm p hcircle hpN
  have hnorm := norm_natAbs_dvd_of_dvd hdiv
  simpa [norm_natAbs_natCast_pow] using hnorm

/-- Lower-bound hook for the absolute norm of the Ramana determinant at one missing prime. -/
theorem pow_two_mul_pairCollisionCount_le_det_ramana_norm_natAbs
    {s : ℕ} {N : ℤ} {z : Fin (2 * s + 1) → GaussianInt} (p : ℕ) [Fact p.Prime]
    (hcircle : ∀ i, z i * star (z i) = (N : GaussianInt))
    (hpN : ¬ (p : ℤ) ∣ N)
    (hdet : (ramanaMatrix s z).det ≠ 0) :
    p ^ (2 * pairCollisionCount p z) ≤ (ramanaMatrix s z).det.norm.natAbs := by
  have hdiv := pow_two_mul_pairCollisionCount_dvd_det_ramana_norm_natAbs p hcircle hpN
  have hpos : 0 < (ramanaMatrix s z).det.norm.natAbs := by
    exact Int.natAbs_pos.mpr (fun hnorm => hdet (GaussianInt.norm_eq_zero.mp hnorm))
  exact Nat.le_of_dvd hpos hdiv

/-- Pairwise-coprime prime-power divisibility, packaged for finite sets of rational primes. -/
theorem prod_prime_powers_dvd_of_each {P : Finset ℕ} {e : ℕ → ℕ} {N : ℕ}
    (hprime : ∀ p ∈ P, Nat.Prime p)
    (hdiv : ∀ p ∈ P, p ^ e p ∣ N) :
    (∏ p ∈ P, p ^ e p) ∣ N := by
  classical
  revert hprime hdiv
  refine Finset.induction_on P ?base ?step
  · intro _ _
    simp
  · intro a s has ih hprime hdiv
    rw [Finset.prod_insert has]
    have hsprime : ∀ p ∈ s, Nat.Prime p := fun p hp => hprime p (Finset.mem_insert_of_mem hp)
    have hsdiv : ∀ p ∈ s, p ^ e p ∣ N := fun p hp => hdiv p (Finset.mem_insert_of_mem hp)
    have hprod_dvd : (∏ p ∈ s, p ^ e p) ∣ N := ih hsprime hsdiv
    have ha_dvd : a ^ e a ∣ N := hdiv a (Finset.mem_insert_self a s)
    have hcop : (a ^ e a).Coprime (∏ p ∈ s, p ^ e p) := by
      rw [Nat.coprime_prod_right_iff]
      intro q hq
      have hqa : q ≠ a := by
        intro h
        exact has (by simpa [h] using hq)
      have hnotdvd : ¬ a ∣ q := by
        intro hadvd
        rcases (hprime q (Finset.mem_insert_of_mem hq)).eq_one_or_self_of_dvd a hadvd with ha1 | haq
        · exact (hprime a (Finset.mem_insert_self a s)).ne_one ha1
        · exact hqa haq.symm
      have hbase : a.Coprime q :=
        (hprime a (Finset.mem_insert_self a s)).coprime_iff_not_dvd.mpr hnotdvd
      exact (hbase.pow_left (e a)).pow_right (e q)
    simpa using hcop.mul_dvd_of_dvd_of_dvd ha_dvd hprod_dvd

/-- Product form of the cancelled determinant divisibility over any finite set of missing
rational primes. -/
theorem prod_pow_two_mul_pairCollisionCount_dvd_det_ramana_norm_natAbs
    {s : ℕ} {N : ℤ} {z : Fin (2 * s + 1) → GaussianInt}
    {P : Finset ℕ}
    (hprime : ∀ p ∈ P, Nat.Prime p)
    (hcircle : ∀ i, z i * star (z i) = (N : GaussianInt))
    (hpN : ∀ p ∈ P, ¬ (p : ℤ) ∣ N) :
    (∏ p ∈ P, p ^ (2 * pairCollisionCount p z)) ∣
      (ramanaMatrix s z).det.norm.natAbs := by
  refine prod_prime_powers_dvd_of_each hprime ?_
  intro p hp
  letI : Fact p.Prime := ⟨hprime p hp⟩
  exact pow_two_mul_pairCollisionCount_dvd_det_ramana_norm_natAbs p hcircle (hpN p hp)

/-- Product lower-bound hook for the absolute norm of the Ramana determinant over finitely many
missing rational primes. -/
theorem prod_pow_two_mul_pairCollisionCount_le_det_ramana_norm_natAbs
    {s : ℕ} {N : ℤ} {z : Fin (2 * s + 1) → GaussianInt}
    {P : Finset ℕ}
    (hprime : ∀ p ∈ P, Nat.Prime p)
    (hcircle : ∀ i, z i * star (z i) = (N : GaussianInt))
    (hpN : ∀ p ∈ P, ¬ (p : ℤ) ∣ N)
    (hdet : (ramanaMatrix s z).det ≠ 0) :
    (∏ p ∈ P, p ^ (2 * pairCollisionCount p z)) ≤
      (ramanaMatrix s z).det.norm.natAbs := by
  have hdiv :=
    prod_pow_two_mul_pairCollisionCount_dvd_det_ramana_norm_natAbs
      hprime hcircle hpN
  have hpos : 0 < (ramanaMatrix s z).det.norm.natAbs := by
    exact Int.natAbs_pos.mpr (fun hnorm => hdet (GaussianInt.norm_eq_zero.mp hnorm))
  exact Nat.le_of_dvd hpos hdiv

/-- Enhanced Ramana obstruction from finitely many missing rational primes.

If the pair-norm bound is smaller than the circle contribution times the determinant divisibility
contribution from the missing primes, then the block cannot consist of distinct points. -/
theorem not_injective_of_pair_norm_pow_lt_missing_prime_product
    {s N : ℕ} {z : Fin (2 * s + 1) → GaussianInt} {B : ℕ} {P : Finset ℕ}
    (hprime : ∀ p ∈ P, Nat.Prime p)
    (hcircle : ∀ i, z i * star (z i) = ((N : ℤ) : GaussianInt))
    (hN : 0 < N)
    (hpN : ∀ p ∈ P, ¬ (p : ℤ) ∣ (N : ℤ))
    (hB : ∀ i j, j ∈ Finset.Ioi i → (z j - z i).norm.natAbs ≤ B)
    (hsmall : B ^ (s * (2 * s + 1)) <
      (∏ p ∈ P, p ^ (2 * pairCollisionCount p z)) * N ^ (s * s)) :
    ¬ Function.Injective z := by
  intro hz
  have hNgi : (((N : ℤ) : GaussianInt) ≠ 0) := by
    intro hzero
    have hNzero_int : (N : ℤ) = 0 := by
      simpa using congrArg Zsqrtd.re hzero
    have hNzero : N = 0 := by
      exact_mod_cast hNzero_int
    omega
  have hdet : (ramanaMatrix s z).det ≠ 0 :=
    det_ramana_ne_zero_of_injective hcircle hNgi hz
  have hlower :=
    prod_pow_two_mul_pairCollisionCount_le_det_ramana_norm_natAbs
      (s := s) (N := (N : ℤ)) (z := z) (P := P) hprime hcircle hpN hdet
  have hupper :=
    det_norm_common_nat_norm_pow_le_of_pair_norm_le
      (s := s) (N := N) (z := z) (B := B) hcircle hN hB
  have hprod_upper :
      (∏ p ∈ P, p ^ (2 * pairCollisionCount p z)) * N ^ (s * s) ≤
        B ^ (s * (2 * s + 1)) := by
    calc
      (∏ p ∈ P, p ^ (2 * pairCollisionCount p z)) * N ^ (s * s)
          ≤ (ramanaMatrix s z).det.norm.natAbs * N ^ (s * s) :=
            Nat.mul_le_mul_right _ hlower
      _ = N ^ (s * s) * (ramanaMatrix s z).det.norm.natAbs := by ring
      _ ≤ B ^ (s * (2 * s + 1)) := hupper
  exact not_lt_of_ge hprod_upper hsmall

end DeterminantDivisibility
end GaussianChain
