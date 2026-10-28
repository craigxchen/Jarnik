# Referee report on `round24/search.md` (twisted Pell ansatz: B.C = 30, two conjugate branches of contact 4, alpha = 15/4)

Referee: adversarial. All scripts: `round24/referee_search_checks/` (independent code; the author's scripts were not
reused, except that `qsearch_m.py` was read to see which `m` it skips). Notes: `round24/referee_search_notes.md` (R1–R17).
Labels: **proved** (proof checked or re-derived here), **exact** (exact computer algebra), **numerical** (stated
precision), **heuristic**.

## 0. Verdict

**The mathematics survives. The headline interpretation does not.**

* **Normal form (Lemmas 1–3): correct.** Re-derived with all quantifiers. One scope remark: the class is the
  *ray-generic* minimal balanced class (§2).
* **Theorem A (real family `M`): correct.** I rebuilt `M` from the formulas of §5.1 with my own code.
  - All exact steps (a)–(d) re-derived symbolically, including the full factorisation of the double-root discriminant.
    The only relevant factor is `F`; the other factors are degenerations.
  - Class, rays and contact re-verified at 7 values of `m` (60 digits) and on a 399-point grid of `(-3,-1)`.
  - The Jacobian of the general system (3.1) has rank 28 (80 digits). Its kernel is exactly the span of the family
    tangent and the gauge.
  - I also supply the short argument, missing from the write-up, that makes "for every `m` in `(-3,-1)`" a proof rather
    than an 8-point check (§3.4).
* **Theorem B (finitely many `Q`-members of `M`): correct as stated for the component `M`.**
  - The discriminant, the 8 odd-order points, Riemann–Hurwitz and Faltings all re-derived.
  - The birationality of `m -> (u,v)` over `Q`, which the proof uses but does not check, was verified here.
  - Step 2 needs a slightly more careful argument, supplied in §4.2; the conclusion is unchanged.
  - The geometric monodromy is in fact the full `S_4`.
  - The rational search was extended to height 300: no nondegenerate rational `m`.
* **Two major problems (interpretation and scope; neither breaks a proof).**
  1. **Scope.** The task's Result line, "only finitely many of its [the class's] members are defined over `Q`
     (proved)", overclaims. Finiteness is proved for the component `M` (and its `S_5`-images) only. Whether the class
     has other components, and whether any of them carries infinitely many `Q`-curves, is open: it is the author's own
     Q24.1. `search.md` §0 item 4 and §5.3 state this correctly; the Result line and the last paragraph of the summary
     do not.
  2. **Attribution to the anisotropy is wrong.** Theorem B is *form-independent*, so it is not an anisotropy effect.
     - The proof uses two inputs. The first is that the 30 places are rational over the splitting field of the torus.
       That holds in **every** `Q`-form: in the split form the same Galois argument makes every ray-generic place
       `Q`-rational, as in split Thm 2(i).
     - The second is the `S_4` monodromy of the places over the invariant curve. That belongs to the **complex**
       component `M_C`, which the split and twisted forms share.
     - Running the proof of Theorem B on the split `Q`-form of `M_C` gives the same conclusion (§5.2): finitely many
       split `Q`-members.
     - So the comparison "split: `W` forced; twisted: finitely many" compares *different classes* (`B.C = 15` against
       `B.C = 30`). The only twist-specific input is still the round-23 parity lemma, which pushes the twisted search
       from `B.C = 15` up to `B.C = 30`.
     - The status label **"anisotropy gain"** is therefore not supported. The correct label is **map value**: the
       literal Pell transplant at `B.C = 30` has a real one-parameter component with only finitely many `Q`-members, in
       either form.
* **Minor points** (§6): the rational-search sentence silently excludes `m = +-1, +-3` (`m = -3` *does* give a
  rational, degenerate `b = 0`); the "for every `m`" claim lacked its argument; Lemma 1's "`b_r` rational" holds only
  up to a common factor; the two-branch covering bound is cited from a source that did not re-derive it; and there is
  an unexplained phrase, "the other branch is `m -> -m`".
* **Bottom line (unchanged).** No growth improvement and no uniform bound: `M <= (2/3+eps) log R/loglog R` stands. No
  owner constraint is violated.

## 1. What is claimed

1. **Lemma 1 (normal form).** A `Q`-curve of the ray-generic minimal balanced class (`E_A.C = 1` for all 30 rays) has
   all 30 places in `P^1(Q(i))`, and `u_j = gamma_j P_j / conj(gamma_j P_j)`.
2. **Lemma 2 (branch condition).** Contact `>= 4` at both branches `+-sqrt m` holds iff `Im(gamma_j P_j) = kappa_j (t^2-m)^4`.
   Then `mu_0 = 4` exactly, and the curve is V-null.
3. **Lemma 3 (s-form).** The condition is (1.3): `e_k(R_j)` is real for `k != 4`.
4. **Count (heuristic).** The expected dimension of the class is 1.
5. **Theorem A.** An explicit real one-parameter family `M`, `m in (-3,-1)`, lies in the class (over `R`). It is a full
   component (Jacobian rank 28) and sweeps a surface.
6. **Theorem B.** `M` has only finitely many members defined over `Q`. The proof goes through Lemma 1, rational `(u,v)`,
   the quartic `G_m(Y)`, its 8 odd-order discriminant points, Riemann–Hurwitz (genus `>= 3`) and Faltings.
7. **Interpretation.** The anisotropy enters through the `Q(i)`-rationality of the places and the `S_4` monodromy; split:
   `Z ⊇ W` is forced, twisted: only finitely many curves are forced.

## 2. Normal form (Lemmas 1–3, §§1–2)

**2.1 Lemma 1. Correct (proved; re-derived).**
* *Valuation vectors.* For `sigma in Gal(Qbar/Q)` and a boundary place `p` of `P^1_Q`, `v(sigma p) = chi(sigma) v(p)` in
  `Z^5/Z1`, where `chi` is the quadratic character of `Q(i)`. This is because Galois acts on the cocharacters of `T_5` by
  `chi`.
* *Places are `Q(i)`-rational.* If `sigma` fixes `i`, then `sigma p` has the same vector as `p`. The 30 vectors are
  distinct rays, so `sigma p = p`. Hence `p in P^1(Q(i))`.
* *Conjugate pairing.* Complex conjugation sends `1_A` to `-1_A = 1_{A^c}` mod `Z1`.
* *Zeros and poles of `u_j = x_j/x_5`.* `ord_p u_j = v_j - v_5`. This is `+1` exactly on the rays `1_A` with `j in A`,
  `5 notin A`; it is `-1` on the conjugates; it is `0` otherwise. Hence `u_j = c_j P_j/Pbar_j`, `N(c_j) = 1`, and
  Hilbert 90 gives `c_j = gamma_j/gammabar_j`.
* *Moving places off infinity.* The stabiliser of `{+-sqrt m}` in `PGL_2(Q)` is infinite. The poles `-a/b` of its
  elements form an infinite set, so the 30 places can be moved off `infinity`.

*Scope remark (minor).* The minimal balanced class also contains curves with *chain* places, whose vector is
`1_A + 1_B` with `A ⊂ B`. Their count is the same but they have fewer places. Lemma 1 and everything after it concern the
**ray-generic** subclass only. §1.2 states this, but the Result line does not.

*Not twist-specific.* The same argument with `chi = 1` shows that in the **split** form every ray-generic place of a
`Q`-curve is `Q`-rational; split Thm 2(i) has 16 rational places. This matters in §5.

**2.2 Lemma 2. Correct.**
* `Im(gamma_j P_j)` has coefficients in `Q` and degree `<= 8`. If `(t - sqrt m)^4` divides it, then so does
  `(t + sqrt m)^4`, by Galois conjugation over `Q`.
* `u_a/u_b` has degree 8: its zeros are at the `2^3 = 8` sets `A` with `a in A`, `b notin A`. So `2 mu_0 <= 8`, and
  equality forces V-null.
* `kappa_j != 0`. Otherwise `P_j` would be real, so its root set would be closed under conjugation, which is impossible
  because `r_T != rbar_S`.

**2.3 Lemma 3 and §2. Correct.**
* `s = (sqrt m + t)/(sqrt m - t)` has determinant `2 sqrt m > 0`, so it preserves the upper half plane.
* `c'_j = u_j(s = infinity) = 1`, because the branch maps to `e`. No phase is set to zero by assumption.
* 2.2 (four roots in each half plane): the argument principle with `arg Q_j` confined to `[0, pi]`, and `arg(+-infinity) = 0`.
* 2.3: the continuous lift from `-infinity` gives `sum Arg rho = arg Q_j(0)`. Re-derived: the `+-pi` contributions of
  the 4+4 roots cancel.
* The count "278 patterns in 15 classes" was recomputed by `c8_signpat.py`: **exact**, confirmed.

**2.4 Count (§3).** It is correctly labelled heuristic. Lemma 2 shows that the drop from 2 (split, `B.C = 15`) to 1
(twisted, `B.C = 30`) comes from **parity** (doubling). The split class with `B.C = 30` and contact `(4,4)` has the same
count: `30 - 28 - 1 = 1`. So §3.3's "anisotropy in the count" is the parity lemma again, nothing new. §3.3's own wording
("parity doubles the places") is accurate. The Result line "(split: 2)" compares different classes.

## 3. Theorem A (the real family `M`)

**3.1 Exact steps (a)–(b). Confirmed** (`c1_sym.py`, sympy, with `b_3 = -(b_0+b_1+b_2)`, `t^2 = prod b`, `a_r = t/b_r`).
* `Im P_j` is a nonzero constant for `j = 1,2,3` (contact 4 for the three pulled-back coordinates).
* `Im n` contains only `X^3` and `X`.
* `L` has a numerator with only even powers `X^2, ..., X^8`, so `L(0) = 0` and `L` is even. `t` cancels everywhere.

**3.2 Step (c): the double-root condition. Confirmed, and completed** (`c1b_disc.py`, flint factorisation).
The discriminant of the cubic `L/Y` factors exactly as

    disc = - b_0^4 b_1^4 b_2^4 b_3^20 · (four cubics) · F_24 · (three quartics)^2 .

* The four cubics are `|p_r|^2 = 0`, i.e. `b_r(b_r^3 + prod_{s != r} b_s) = 0`, where `sigma = 0` and the leading
  coefficient of the cubic vanishes. Degenerate.
* The three quartics divide `A_0 = (L/Y)(0)` and occur squared. On them `L = A_3 Y^2 (Y + e)^2`, i.e. `D = 0`. The conic
  `z^2 = -X^2/sigma` is then a **line pair**, and each line has a single rational branch, giving `alpha >= 4`. Degenerate
  for this class.
* `F_24` is the only relevant factor. **Exact:** `F_24 = F(e_2, e_3, e_4)` with exactly the `F` printed in §5.1(c).
So within the reflection-symmetric `x`-structured ansatz, `M` is the only nondegenerate family. The ansatz without the
reflection symmetry was not analysed (by either of us).

**3.3 Step (d). Confirmed** (`c2_phi.py`).
* `Phi(u,v) = F/e_2^12` is the quintic `-256u^4 - 256u^3v^2 + 1024u^3v + 32u^3 - 992u^2v^2 - 96u^2v + 288uv^3 + 37uv^2 - 27v^4 - 4v^3`.
  It is irreducible over `Q`.
* `Phi(u(m), v(m)) = 0` identically.
* `Phi(u, ku) = -u^3 · (quadratic in u)`, whose discriminant is `(k-4)^5 (9k-4)^3`, as claimed.
* **New check, needed for Theorem B step 3.** `gcd_{m'}(num(u(m) - u(m')), num(v(m) - v(m'))) = m - m'`. So
  `m -> (u,v)` has degree 1: it is birational over `Q`, and `Phi` is geometrically irreducible of genus 0.
* The sentence "The other branch is `m -> -m`" in §5.1(d) has no counterpart in this computation: there is only one
  branch. It is harmless but should be removed or explained.

**3.4 The class, contact and nondegeneracy. Confirmed; "for every `m`" now proved.**
* *Pointwise check* (`c4_construct.py`, own implementation, 60 digits) at `m = -1.01, -1.1, -1.3333, -1.5, -2, -2.5, -2.95`:
  - 30 places, and all 30 proper subsets of `[5]` occur exactly once.
  - Contact slopes `4.0000` at `w -> infinity` and at `w -> 0` for all four `u_j`, in the real conic coordinate
    `w = z - sqrt(kappa) X`.
  - The minimal separations agree with the author's (3.18 at `-1.1`, `8.8e-5` at `-2.95`).
* *Grid* (`c7_scan.py`): 399 equally spaced `m` in `(-3,-1)`, 40 digits. All nondegenerate; `e > 0`, `D > 0`, `e < D`,
  `sigma < 0` throughout.
* *Why "for every `m` in `(-3,-1)`" holds* (**proved here**; the write-up only checks 8 points):
  - The 16 `X`-values of places are `+-t/b_r +- i b_r`.
  - They are pairwise distinct iff the `b_r^2` are distinct and nonzero. On `(-3,-1)` this holds because `disc_Y G_m`
    and `u` have no zero there, and `Y_r = -b_r^2/v`.
  - Each value has `Re X_0 = +-t/b_r != 0`. So it is neither a branch point of `X` (those have `X_0^2 = -D < 0`) nor a
    zero of `beta = X(X^2+e)`. Hence `H` vanishes at exactly one of the two conic points over each root of `n`.
  - `kappa_j != 0` is automatic (§2.2). The four `b_r` stay real on the interval, since no discriminant zero lies in it.
  - So the configuration never degenerates on `(-3,-1)`, and the ray labels are locally constant there.
* *Jacobian* (`c5_jac.py`, 80 digits, at `m = -1.5` and `-2.2`):
  - The 15 places in the `w`-coordinate satisfy (3.1) to `5e-47`. This is limited by the floating double root.
  - The `28 x 30` Jacobian has smallest singular values `0.089` (resp. `0.011`) against largest values `1e12` (resp.
    `3e9`). So the **rank is 28**.
  - The family tangent `d rho/dm` and the gauge vector are both annihilated, to `1e-15`, and are independent.
  - So the real solution set is locally exactly `M` × gauge: **`M` is a full component** (numerical, with a wide margin).
* *Sweep* (`c9_sweep.py`): the map `(m, w) -> (arg u_j)` has rank 2 at 4 sampled points. That the image is not a
  hypersurface is automatic: one parameter plus the curve parameter gives dimension at most 2.

**Verdict on Theorem A.** Correct, with the labels the author gives. With 3.4 the nondegeneracy on the whole interval
becomes a proof (exact discriminant information plus elementary inequalities), apart from the fixed ray labelling,
which is read off at one point.

## 4. Theorem B (finitely many `Q`-members)

**4.1 Steps 4–9. Correct** (`c3_genus.py`, exact).
* *The quartic.* From `W^4 + vW^2 - v^2W + uv^2 = 0` and `Y = -W^2/v`: `Y^2 - Y + u = W`, hence `(Y^2-Y+u)^2 + vY = 0`.
  Checked by hand and by resultant: `Res_W(b-quartic, Y + W^2/v)` is proportional to `G_m`.
* *`b` and `Y` determine each other.* `b = Y^2 - Y + u` is a polynomial in `Y`, so `b` is rational iff `Y` is. The two
  quartics have the same splitting field.
* *Discriminant.* In my normalisation `G = Y^4 - 2Y^3 + (1+2u)Y^2 + (v-2u)Y + u^2`, cleared to leading coefficient
  `256(m-3)^4(m-1)^4`:

      disc_Y = -2^47 (m-3)^12 (m-1)^16 (m+1)^2 (m+3) c(m)^3 (m^4+70m^2-48m+9).

  The author's exponents `42, 34` at `m = 3, 1` differ from mine by even amounts, which reflects only the normalisation.
  The **parities agree**, and parity does not depend on rescaling: rescaling multiplies the discriminant by a 6th power.
* *The 8 odd-order points.* They are the roots of `(m+3) c (m^4+70m^2-48m+9)`, and they are distinct:
  `disc(c) = -26624`, `disc(quartic) = 311427072` and `Res(c, quartic) = 2^23` are all nonzero. The leading
  coefficient does not vanish at them. At `m = infinity` the order is 4 (even).
* *Odd monodromy.* Odd order of the discriminant implies odd local monodromy. The reason: `prod_{i<j}(Y_i - Y_j)`
  changes sign under continuation around the point.
* *Riemann–Hurwitz.* `2g - 2 >= -2|N| + 8·|N|/2`, so `g >= |N| + 1 >= 3`.
* *Faltings.* A component of the ordered-root curve that is not defined over `Q` has only finitely many rational points
  (they lie on its intersection with a conjugate). The other components are curves over `Q` of genus `>= 3`.
* *Strengthening.* `Gal(G_m/Q) = S_4` at `m = -2, -3/2, 5, 7/3` (sympy). So the group over `Q(m)` is `S_4`. The
  geometric monodromy is a normal subgroup containing an odd element, hence equals `S_4`. "S_4-type monodromy" is
  therefore literally true.

**4.2 Steps 1–3. Correct after a small repair.**
* Step 3 needs `m -> (u,v)` to be birational over `Q`. Verified in 3.3 above.
* Step 2 as written ("the centring `X = x + q_1/2` and the scaling are defined over `Q`, hence `p_r in Q(i)`, so
  `b_r in Q`") silently assumes that the member is a *real* point of the explicit family. Over `Q` one must also handle
  the other points of the complex component `M_C` (complex `b_r`). Repair, valid for every point of `M_C`:
  - In the algebraic normal form the 16 place `X`-values are `+-t/b_r + i b_r` (labels `S`) and `+-t/b_r - i b_r`
    (their partners `S^c`).
  - Let `X'` be a `Q`-coordinate with the branch image at `infinity`, so `X = lambda (X' - X'_0)` with `lambda`,
    `X'_0` in `C`.
  - For a twisted `Q`-member the `X'`-value of `S^c` is the complex conjugate of that of `S` (Lemma 1).
  - The 8 base values sum to `2i sum b_r = 0` in `X`, and so do their partners. So `X'_0` equals its own conjugate and
    lies in `Q(i)`; hence `X'_0 in Q`.
  - Then `(t/b_r +- i b_r)/lambda` is a conjugate pair in `Q(i)`. So `b_r/lambda in Q` for all `r`, with one common
    `lambda`.
  - `Y_r` and `(u,v)` are invariant under common scaling, so they are rational.
* With this, **Theorem B is proved** for the Zariski closure of the explicit family. Its finitely many exceptional
  parameter points, where the parametrisation degenerates, contribute finitely many curves.

**4.3 Search.** My independent search (`c6_qsearch.py`, flint factorisation of `G_m` over `Q`) covered all `m = p/q`
with `|p|, q <= 60` (4405 values) and then `|p|, q <= 300` (109589 values). The **only** `m` with a rational root is
`m = -3`. There `u = 0` and `G = Y (3Y-4)(3Y-1)^2/27`, so `b = 0`: the degenerate endpoint where `t = 0`. The author's
script skips `m in {+-1, +-3}`, but §5.3 says "no `m` ... gives even one rational root". That should read: no `m` other
than the degenerate values `+-1, +-3`. Minor.

## 5. Bookkeeping, and where the anisotropy actually enters

**5.1 Balanced, through `e`, over `Q`, exponent 15/4. Correct for the class as defined.**
* Every cut has weight `E_A.C + E_{A^c}.C = 2`, so the curve is balanced (each toric height is `2h(t) + O(1)`).
* `-K` is base-point free (nef toric), so `h_{-K} = 30 h(t) + O(1)`.
* `dist ~ |t - t_0|^4`, and Liouville at a real quadratic `t_0` gives `dist ≍ H_{-K}^{-4/15}` infinitely often, never
  better. At `gamma = 15/4`, i.e. `delta = 1/30`, such points are A-points for every `delta < 1/30` ✓.
* A real member of `M` is a curve over `R` only. For a `Q`-member one needs in addition: `-1/sigma` not a square (else
  the branches are rational and `alpha = 15/2`), the conic to have a `Q`-point, and `b_r in Q·lambda`. The write-up
  says this.
* §1.2's claim that two-branch curves "cannot cover `T_5`" cites `split.md` 5.2. That source says it did *not* re-derive
  the two-point form of the covering inequality. The claim is not used in Theorems A and B. Minor: it should carry the
  same label.

**5.2 The anisotropy claim. Not supported (major).** §5.5 and §6.1 say that the anisotropy enters "through Lemma 1 (all
30 places `Q(i)`-rational) plus `S_4`-type monodromy". Both ingredients are form-independent.
* *Place rationality.* In the split form `chi = 1`, so every ray-generic place of a `Q`-curve is `Q`-rational (§2.1).
  This is a *stronger* rationality condition than `Q(i)`-rationality.
* *Monodromy.* `(u, v, G_m)` are functions on the complex component `M_C`. Over `C` the split and twisted 4-folds are
  the same toric variety, so `M_C` is a family of curves on both.
* *Theorem B for the split form* (**proved**, same steps as 4.1–4.2):
  - For a split `Q`-member, every place `X'`-value is in `Q`.
  - The same centring argument gives `X'_0 in Q`. Then `(t/b_r + i b_r)/lambda` and `(t/b_r - i b_r)/lambda` are both
    in `Q`. Hence `i b_r/lambda in Q`, i.e. `b_r in (lambda/i)Q`.
  - So `Y_r = b_r^2/(-e_2) in Q` and `(u,v) in Q^2` (both scale-invariant), so `m in Q` and `G_m` splits completely.
  - The `S_4`/genus/Faltings argument applies verbatim. So `M_C` has only finitely many split `Q`-members as well.
  - Under `b -> i b` (or any common scaling) `u = e_4/e_2^2` and `v = e_3^2/e_2^3` are unchanged; this is the
    algebraic content.
* *Consequence.* Theorem B cannot distinguish split from twisted.
  - The contrast drawn in §5.5 ("split: `W`, rational in its places; twisted: finitely many") sets a `B.C = 15` family
    (split Thm 2) against a `B.C = 30` family (`M`).
  - The only twist-specific step separating them is the round-23 parity lemma, which forbids `B.C = 15` on the twisted
    form.
  - Within `B.C = 30`, nothing in this round separates the forms.
  - Whether the split form has *other* `B.C = 30` components with infinitely many `Q`-members is not known either, and
    it would not matter for the twisted question.
* *Correct reading.* The literal Pell transplant, forced up to `B.C = 30` by parity, lands on a component whose
  place-monodromy over its rational invariant curve is `S_4`. That kills it over `Q` in every form. This is **map
  value** (one component eliminated), not an anisotropy gain. The status `anisotropy_gain` should be downgraded.

**5.3 Owner constraints.** Audited:
* No fixed-`S` Subspace exceptional sets are used.
* No gcd or denominator cost is dropped: `gamma_j` is kept, and the heights are full toric heights.
* No common phase is set to zero: `c'_j = 1` is forced by `u_j(branch) = 1`, a property of the curve, not a drifting
  phase.
* No threshold is claimed uniform: Theorem B is ineffective and no radius bound is asserted.

No violation.

## 6. Evidence versus proof; minor issues

**6.1 Labels.** The author labels evidence honestly.
* Theorem A's nondegeneracy is stated as checked at 8 points. §3.4 above upgrades it to a proof on the whole interval.
* The component claim (rank 28) is numerical, with a wide margin; it is reproduced here.
* The expected-dimension table is heuristic and is labelled so.
* The blind searches (§5.4) are explicitly **not** used as evidence of emptiness, and the earlier "evidence of
  emptiness" reading is withdrawn. Correct, and an important methodological correction. It also applies to
  `toric.md` §7, as the author notes.

**6.2 Minor issues (none affects a proof).**
1. *Ray-generic scope.* The Result line's "the class" should read "the ray-generic minimal balanced class". Chain places
   (`1_A + 1_B`) are not covered.
2. *"`b_r` rational" holds only up to a common factor.* The Result's "a Q-member needs all four `b_r` rational (Lemma 1)"
   is true only after normalising the scale. Theorem B uses only the scale-invariant `Y_r`, so no harm is done.
3. *Search sentence.* §5.3 says no `m` up to height 60 gives a rational `b`; the script excludes `m = +-1, +-3`, and at
   `m = -3` there is a (degenerate) rational root `b = 0` (§4.3).
4. *Discriminant exponents.* `(m-3)^42 (m-1)^34` depend on the normalisation. Only the parities matter, and these were
   confirmed.
5. *"The other branch is `m -> -m`"* (§5.1(d)) is unexplained. `Phi` is irreducible and the parametrisation is
   birational (§3.3).
6. *Covering bound.* §1.2 cites the two-branch covering bound from `split.md` 5.2, which did not re-derive it. Label it.
7. *Step 2 of Theorem B.* Repaired in §4.2. The written version covers only real points of the family.

**6.3 What I could not refute.**
* I tried to break Lemma 1: the Galois twist bookkeeping, the conjugate pairing `1_A <-> 1_{A^c}`, and the zero/pole
  rule.
* I tried to break Lemma 2: the degree count `8`, `kappa != 0` and V-null.
* I tried to break Theorem A: an independent construction, all discriminant factors, rank, kernel and the whole-interval
  nondegeneracy.
* I tried to break Theorem B: the resultant, the discriminant, the odd points, the leading coefficient, birationality,
  RH and the descent step.
* Each of them survived, after the repair to step 2.

## 7. Corrected statement

> **Round 24 search (corrected).** Let `X = Lbar_5^tw`. Consider the *ray-generic* minimal balanced class of rational
> curves through `e`: `E_A.C = 1` for all 30 rays, 30 distinct places, and two conjugate real-quadratic branches of
> contact `mu_0 = 4`, so `alpha_C = 15/4`.
> 1. **(Normal form; proved.)** A `Q`-curve of this class has all 30 places in `P^1(Q(i))`. It has the form
>    `u_j = gamma_j P_j/conj(gamma_j P_j)` with `Im(gamma_j P_j) = kappa_j (t^2-m)^4`, `kappa_j in Q^*`. Equivalently,
>    in the `s`-form, `e_k(R_j)` is real for `k != 4`. The class is V-null with `mu_0 = 4` exactly. *(The analogous
>    statement in the split form, with `Q` in place of `Q(i)`, holds by the same proof.)*
> 2. **(Theorem A; proved, with the component claim numerical.)** The explicit family `M`, `m in (-3,-1)`, of §5.1
>    consists of real curves of this class for **every** `m` in the interval.
>    - It is the only nondegenerate family inside the reflection-symmetric `x`-structured ansatz. The other
>      double-root factors are `sigma = 0` or a line-pair conic.
>    - The general system has Jacobian rank 28 along it (numerical, 80 digits). So `M` is an irreducible component of
>      the real solution set, of dimension 1 modulo gauge, and it sweeps a real surface in `T_5(R)`.
> 3. **(Theorem B; proved.)** Only finitely many curves of the complex component `M_C` (hence of `M` and of its
>    `S_5`-images) are defined over `Q`. The proof goes through the scale-invariants `Y_r = b_r^2/(-e_2)`. These are
>    the roots of `G_m(Y) = (Y^2-Y+u(m))^2 + v(m) Y`, whose discriminant has exactly 8 odd-order points. The geometric
>    monodromy is `S_4`, the Galois closure has genus `>= 3`, and Faltings applies.
>    - The bound is ineffective.
>    - None has `|p|, q <= 300`, apart from the degenerate `m = -3`.
>    - **The same finiteness holds for the split `Q`-form of `M_C`.** So Theorem B is not an effect of the anisotropy.
> 4. **(Not proved.)** Finiteness for the whole class is not proved: other components may exist (Q24.1). The
>    expected dimension 1 is heuristic. The only twist-specific input remains the round-23 parity lemma, which forces
>    `B.C >= 30`.
>
> **Status:** map value (one component of the twisted Pell-type class eliminated over `Q`). **Not** an anisotropy gain.
> No growth improvement and no uniform bound: `M <= (2/3+eps) log R/loglog R` stands.

## 8. Scripts (`round24/referee_search_checks/`)

| script | content | label |
|---|---|---|
| `c1_sym.py` (`c1_L.pkl`) | symbolic (a),(b): `Im P_j` constant, `Im n` odd cubic, `L` even, `L(0)=0` | exact |
| `c1b_disc.py` (`c1b_factors.pkl`) | discriminant of `L/Y`, full factorisation (flint); `F_24` isolated | exact |
| `c2_phi.py` | `F_24 = F(e_2,e_3,e_4)`; `Phi(u,v)`; `Phi(u(m),v(m)) = 0`; pencil discriminant; birationality of `m -> (u,v)` | exact |
| `c3_genus.py` | `G_m` from the resultant; `disc_Y`; odd-order points; leading coefficient; coprimality; `m = infinity` | exact |
| `c4_construct.py` | independent construction from `m`; places, rays, contact slopes (60 digits) | numerical |
| `c5_jac.py` | residual of (3.1); `28 x 30` Jacobian SVD; tangent and gauge in kernel (80 digits) | numerical |
| `c6_qsearch.py` (`c6_q300.out`) | rational roots of `G_m`, height `<= 60` and `<= 300`; Galois groups at 4 values | exact |
| `c7_scan.py` | 399-point scan of `(-3,-1)`: nondegeneracy and signs of `e, D, e-D, sigma` | numerical |
| `c8_signpat.py` | 278 half-plane patterns, 15 classes | exact |
| `c9_sweep.py` | sweep rank 2 | numerical |
