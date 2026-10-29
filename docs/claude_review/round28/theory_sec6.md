## 6. Consequences for the pre-registered decision rule

The rule is unchanged (it was fixed before the data). What the theory contributes to reading the `d = 7 -> 8 -> 9` data:

1. **The theory does not settle the rule.** It proves neither `reach <= zeta + c` on regime boxes (false in general,
   Cor 2.2; conditional and design-dependent in the form of Cor 2.4) nor any mechanism forcing `reach > zeta`.
2. **Use `zeta_pair`.** The criterion's hypothesis is `Y^l != 1` for `l in L_pair`; the fair comparison is with
   zonotopes whose generators lie in `L_pair`. On the four `L >= 2` models `zeta_pair = zeta_all` (§1), but at `d = 8, 9`
   the two can differ, and only `zeta_pair` is the same-hypothesis competitor.
3. **Design effects.** A failure at `e = 2h(sigma)` with a rigid top level of the base coordinate (Sylvester boxes with
   `L` a power of 2; M8), or more generally a defect spanned by powers of implied slabs with rigid level sets (§3.4),
   is predicted by Theorem 2.3. Such failures are arithmetic coincidences of small symmetric profiles; they are evidence
   about designs, not about generic regime profiles. Recording, for each failure, `h(sigma)`, the rigidity of its top
   level (`t7`/`t10` style) and the dimension of the span of power obstructions (`t6`) makes this visible.
4. **Vacuous favourable cases.** Product (block) profiles satisfy the favourable pattern trivially (Cor 2.2), and
   vertex cuts (`L = |P| - 1`, like M8') are cube-bound rather than slab-bound. Neither resembles the regime
   (`j/|P| -> 0`, character gap). Favourable evidence should come from gap-satisfying, slab-bound (`L <= |P|/2`) profiles.
5. **Constant chain versus chains.** The computed "reach" is the constant-chain reach. Corollary 4.1 works verbatim
   with any subset or chain (Theorem 3.1, §5.1), and `reach_chain >= zeta_L` by construction. So "reach <= zeta" on
   every box would show that the *constant chain* never beats zonotopes, not that the interpolation route is closed.
   I do **not** propose to change the rule; I flag that the cheap extra statistic `reach(A cap {|sigma| <= h-1})`
   (one more rank computation per box) measures whether the BCO is the only obstruction (on M8: 6 -> 7 = `zeta_pair`).
6. **Probabilities (judgement).** The round-27 estimate (favourable trend ~15%; proof of `(SHIFT)` at regime size within
   a few rounds 5–10% conditional on that) is not moved much. The theory removes one reason for pessimism (the M8 failure is a design
   effect, and part of the sdc7 failure is a rigid-level effect of the same kind) and one reason for optimism (M8' succeeds because its cuts are weak and
   design-free). If the `d = 8, 9` data come out unfavourable on slab-bound, non-design, gap-satisfying boxes, the
   referral to experts is the right call; nothing in this note gives an internal route to a proof of `(SHIFT)` at
   regime size.
