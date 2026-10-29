## 4. Part (b): the syl8 `e = 6` obstruction, identified

**Theorem 4.1 (exact).** Let `A` = M8 (syl8, `L = 2`, `d = 7`, the 7 patterns are the complements of the 7 lines of the
Fano plane on the non-base points `0..6`), and `sigma = k_0 + ... + k_6 = -k_8` the **base-point coordinate**. Then
1. `sigma(A) = [-3, 3]` and `sigma^7 in J(A)_7` (exact);
2. `sigma^6 notin J(A)_6` (exact over `Q`, `t15`);
3. `sigma^6 in J(A_l)_6` for every `l in L_pair` (proved: Theorem 2.3, with the ranges below; also checked exactly for
   four shifts, `t15`);
4. the defect `cap_{l in L_pair} J(A_l)_6 / J(A)_6` is one-dimensional (mod p, `t1`; the invariant part is computed
   exactly over `Q`, `t2`), hence it is spanned by the class of `sigma^6`. The "one-dimensional defect that is not a
   product of column forms" of round 27 (`t6_obstruction`) is the **sixth power of the base coordinate**.

*Proof of 3.* `sigma(e_x) = 1`, so `sigma(A_{+-e_x})` has at most 6 values. For `l = e_x - e_y`, `sigma(l) = 0`; the top
level `{k in A : sigma(k) = 3}` is exactly the set of 7 indicator vectors `1_l` of Fano lines (exhaustive enumeration),
and `1_{l'} - 1_l` has weight `|l Delta l'| = 4` for distinct lines, so no top point survives a weight-2 move: the top
level is rigid and `sigma(A_{e_x - e_y}) subset [-2, 2]` (in fact equal). By central symmetry the bottom level behaves
the same way for `-l`. Theorem 2.3 applies with `R = 7`. QED

**Why `|sigma| <= 3`.** Each non-base point lies on 3 lines, hence in `7 - 3 = 4` patterns, so `sum_p s_p = 4 sigma`
with `|s_p| <= 2`: `|sigma| <= 14/4`. This is Proposition 2.5 for the Sylvester profile (`f = 4/7` with uniform
weights; the averaging bound is attained). The top level is the "design": the Fano lines are the 3-sets meeting
every line complement in at most 2 points.

**Why the weight-4 shifts remove it.** For two lines `l, l'` (meeting in one point) the shift `1_{l'} - 1_l` has weight
4 and `sigma = 0`, and carries the top point `1_l` to the top point `1_{l'}`; on `A_{1_{l'} - 1_l}` the base coordinate
again takes 7 values and the range argument no longer applies. This matches the killer's observation that
`(SHIFT_6)` holds with all `{-1,0,1}^7` shifts and that "the rescue needs shifts of weight `>= 4`" (round 27 §5.2(b)).

**Symmetry.** The coordinate symmetry group of M8 is `GL(3,2)` (order 168; the base point is fixed, `t2`), it fixes
`sigma`, and the defect is its trivial representation. Over `Q`: `dim J(A)_6^G = 2`, `dim J(A_l)_6^G = 3` for both
orbit representatives, and the invariant obstruction found by linear algebra is `F == c . sigma^6 mod J(A)_6^G`
(`t5`). It already lies in the first-order closures `J1(A, l)` (`t4`).

**Predictions from the same mechanism (not computed: `d >= 15`).** Sylvester boxes with `L = 2^{m-1}` have the rigid
subspace level (§2.3; IP-verified for `(n, L) = (16, 2), (16, 4), (32, 8), (64, 8), (64, 16), (128, 16)`), hence
`reach(A, L_pair) <= 2h = 4L - 2` whenever `(N)` holds there. For `(n, L) = (16, 2)` that is `reach <= 6`.

**The next failure of M8 (`e = 7`).** `delta(7) = 34` on `A`; after deleting the 14 points `|sigma| = 3` the subset
`A'` has `reach = 7` and `delta(7) = 28`, and that whole defect is spanned by 70 powers `phi^7` of ternary forms
(`t9`). Removing the extreme levels of these forms one at a time does not reduce it (28, 27, 26, 29 over four steps).
