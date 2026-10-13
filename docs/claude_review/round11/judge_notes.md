# Judge notes, round 11

## Context
- round10 README §8: recommended ZE-real test. §8.4 model conditions (i) slabs cut, (ii) 2w < m, (iii) |A| > binom.
- determinant.md Thm 5.7: ZE(n') for all n'>=n_1 => n < 4 max(n_1, 300 ln 4r +1)+4 on squarefree norms.
- referee Prop R1: uncut box => ZE holds.

## Reading log
### compute.md (read)
- Prop 2.1 zonotope cert; Prop 2.2 zeta(A)<m in regime (uses char ineq, |v|_1<=2D); Prop 2.3 section principle;
  Prop 2.4 n-1=2; Prop 3.1 Hadamard zeta <= ((n-1)/n) max h_A.
- M8: H_8, T=1, L=2, m=8, |A|=1037, dimV=42; M8': L=3, m=11, |A|=1877, dimV=7. Thm 5.2 exact via Lemma 5.1 (factor + Bayer saturation).
- Prop 6.2: skew double core q=7, T=1,L=2,m=8, dimV=1 -> real hypersurface through e w/ nondeg points.
- ZE-braid proposed (§5.3). Sym M8 variant unfinished.
- Key structural: in M8' each P_c = w^a H_c(w-1), H_c form deg m. Bs subset e + V_C(Sigma_A).
### referee_compute.md (read)
- All exact confirmed with own code (flint). Prop 2.4 repaired: every n=3 slab-box is zonotope, zeta=rho.
- Major: ZE-braid false at Sylvester n=2^k>=2^16 (k-col subtorus in Bs, Prop 5.2 kappa=0). (I check numerics: w~1.676e8 at n=2^16; 16*w ~2.68e9 < m=3.16e9.)
- Major: M8, M8' fail (ii) h_A<m (8=8, 12>11); regime h_A<=m/2. For H_8 (ii) unreachable at T=1 and T=2,L<=2.
- So Thm 5.2 = correct theorem on LP-subcritical toy models; weak evidence.
- JUDGE THOUGHT: multi-column subtori in regime Bs are DEGENERATE (they are the Prop 5.2 strata); so ZE must allow
  rational subtori of dim up to ~ m/w. Algebraic replacement: "ZE-tor": Bs near e = union of rational subtori through e.
### theory.md (read)
- Lemma 1.1: codim(V_H in V(S,m)) <= (U ln2 + E_S ln(UB))/ln(2H+1). Lemma 1.2: dichotomy on any S in B_T, cost<=3mW/8, H<=e^{mW/40}.
- Prop 2.5: tangent cone criterion (2.1) inconclusive near low-rank strata (eta < |u|/(10 n^1.5)); need X_v<W for some v perp L.
- Thm 3.1: widen slab of light p (log p<=W/8) to |<k,beta_p>|<=m/e_p; no UM_p vector in V_0(A_p). SL(P,p) fails for actual.
- Thm 3.3: SL(n) (8 columns) => n=O(loglog R) squarefree.
- Prop 3.5 obstructions: levels (h<m), a'_p in Q-span of cheap family.
### referee_theory.md (read)
- Fatal: SL(n) false (Prop R1 twins/duplicates: l_p capped by narrow twin slab, w+4T<m; barrier profiles admissible, fail at all but 1 (Paley fake) / 2 (sdc) columns). Doubled Hadamard: 0 columns.
- Major: parity. |l_p|<=m/2 -> odd m gives <=m levels -> SL fails for all profiles. Lemma R2(b): even m, UM_p only at bottom level.
  m odd for 4076/8137 n in [64,8200].
- Repair: Thm R3 grouped widening (Q,lambda) light: sum lambda_q log q <= W_P/8. SL*(n) open. Thm R4: SL* => n=O(loglog R).
- SL* first test = barrier profiles (Q = label/twin class). q=7 illustration: grouped widening gives verified UM_p.
- JUDGE: tropical reading: SL-full fails at p iff (l_p(k))_k in trop(Poly_{<m}|_S over C_p) [realisability of tropical linear spaces].
### subcluster.md (read §0-5)
- Lemma 1.0: delta in L, nondeg rel. L. Lemma 1.1 universal admissibility (ord_L).
- Props 2.1-2.2: dim Z_comb <= n-2; fibres circle or <=2T pts. Uninformative.
- Thm 3.1 push-forward; (c) nondeg J => injective. Prop 3.2 product push-forwards in O(log) window of characters.
- Thm 4.1 (a) rho_hat>=300 ln4r => visible to Z_comb; (b) M>=(rho+1)(rho+2) => visible to Z_univ via linear forms in points (D=2).
- Cor 4.2 blind to univ => M<(300 ln4r+2)^2.
- Thm 4.4 reg(S) <= E_eps h_S(<eps,.>) (new; strengthens threshold Thm B).
- Prop 5.1: ZE-tor => blind => rank target equiv B(theta). Thm 5.4: ZE-col(univ) => M=O((loglog R)^2) all norms.
- Remark 5.5: z_x = Q(hat a_x - hat a_b) Q quadratic, for low-rank clusters with W >> E ln(MB).
### referee_subcluster.md (read)
- Prop R: thin family (one measure on Z^3, ord 2) has nondeg zeros near 0 -> ZE-col false as stated; ZE-col* (V1)/(V2) at actual delta.
- ZE-tor must be local (lifted ball, rational subspaces). "nondeg subclusters add nothing" unsupported. "empty below 2137" wrong: univ family LP-supercritical from n=440, admissible from 468 (r=n).
- §6.5 TC = {0} trivially because m<=2T.

## JUDGE FINDINGS (to verify)
J1. SL* FALSE for Sylvester-Hadamard n=2^k>=2^14 (2^13 with exact m): identity
    1_{a.x=0} = 1 - 1_{b.x=0} + 1_{(a+b).x=0} - 2*1_{a.x=1,b.x=1}; so a'_p = -a'_b + a'_{a+b} - 2 sum_{x in V''} e_x,
    |V''| = n/4; cheap cost <= 2*2floor(w/2) + (n/4)*2T < m. Choose pair {b,a+b} disjoint from Q (light Q small).
    Prop 3.5(ii) (refereed) => no UM_p for any light (Q,lambda), p in Q. Resolves referee's open question negatively.
J2. Chebyshev: |Phi_c(u)| <= |c|_1 * 2 (R/2)^m, R = D/2 => kappa <= 2^{1-m}. (2.1) needs kappa > (D/2)|u|/(m+1), |u|>= ~e^{-W/2}
    => needs W >~ 2 ln2 (m-1) ~ 1.03 n^2 at T=n. Near-extremal W ~ 3n ln n: (2.1) NEVER holds. Lemma 2.4 "very weak requirement" wrong.
J3. T=1 reparametrisation: Thm 2.1 (a),(c),(d) and Thm 4.1 work at T=1 (m ~ c n, regime n >~ 3500?), Liouville W >= ~57+8lnC.
    Then m << W, actual delta at e^{-W/4} deep inside germ (room e^{(W-2.77m)/4}); tangent cone becomes the right object.
    Relative expansion at nearest stratum point: only absolute distance needed -> Prop 2.5 obstacle disappears (heuristic).
## Judge check results
- sl_star_hadamard_out.txt: identity exact k=3..7 (0 failures); cost<m: n=2^13 (exact m 49388999 vs 47.08e6), all k>=14 via proved m>=0.7n(n-1). m values match referee.
- scale_T_Tn_out.txt: T=n: W_TC ~ 1.03 n^2: 242 x W_min at n=2137. (2.1) vacuous for near-extremal.
- scale_T_T1_out.txt: T=1: m ~ 0.642 n; m>=4w from n=2900 (r=n), 3623 (n ln n), 5809 (n^2); W_1 (Thm4.1, C=1) = 51.75 at n=3000;
  log room for (2.1) at W_min: 6178 (n=3000).
- Near strata (judge heuristic): Phi_c vanishes to order >= m - d_F transversally along cheap span L_F (d_F = cost). (2.1) near L needs
  eta^{m-d_F} >~ e^{-W/4} 2^{O(m)} -> at T=1, d_F<=m/2: eta >~ n^{-2.3}; at T=n: eta ~ 1 impossible. Pair ineq gives eta >= e^{-(d_xy-W/2)/2}/(C sqrt n).
  Quantitative Lojasiewicz step is the risk.
- README written: §1,2,3,4,5 (J1-J4),6. J4 = corollary of Thm 5.2 factor data (Phi_c = m! prod tangent forms * prod H_f). GB run in background (tc_m8.py) is an independent check only.
- TODO: §7 meaning, §8 next step (T=1 ZE-TC_eps programme), §0 bottom line + scores, §9 files.
- README complete (all sections). m'=ceil(4w) refinement added (scale_m4w_out.txt): W_1(m') 53..631, cheap long dirs <=4, rho0 beats pair eta from n~1e4.
- Corrections applied: removed wrong "new" remark on Remark 5.5 (it's character-level); softened "at most four cheap columns" to long directions (sketch).
