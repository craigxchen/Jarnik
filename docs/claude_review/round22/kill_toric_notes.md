# kill_toric notes (scratch)
## Read so far
- toric.md read in full; gauss.md, kill_gauss.md, expert_note read. construct.md Thm B/B1 read (ord<=2^{1-M} sum_P w_P, homogeneous slice).
## Findings (candidate issues)
F1 (major?): Thm3.4(b)/Remark(1)/§0 "toric Liouville with polynomials of degree<=7 cannot improve CC 1/10": data t_max(m<=7) only bounds
   single f with max cut width <=7 (via f f^*) ; f in H^0(7B) can have width up to 14; and finite families: min over convex hull of
   width/ord vectors corresponds to PRODUCTS (large widths) -> not covered by m<=7 data. Overclaim.
F2 (minor): alpha_e(B)=2 lower bound over all T_5(Q) not proved in text; easy proof via f_ab=(x_a-x_b)(x_a^-1-x_b^-1) in H^0(B), ord 2.
F3 (minor): "no toric argument" = meta-claim; Prop 6.1 tautological (invariants of transitive action constant).
F4 (minor, lit): McKinnon-Satriano Trans AMS 374 (2021) 3557-3577 (arXiv 2004.05212) not cited: split toric, "P better approximated by
   a rational curve than any Zariski dense sequence" (Q-factorial terminal); erratum 2023: Vojta => canonical boundedness when -K big.
   Monahan-Satriano 2308.11847 horospherical, conditional on Vojta. Huang ANT journal ref NOT confirmed by search (arXiv 1809.02001 v5 Apr 2021).
F5 (minor): Lemma1.2 "smooth Q-conic": closure of 1-PS orbit can be singular at boundary pts; only smooth at e is used.
F6 (minor): Thm 3.3 / §6.1 "every argument with aux Laurent poly is Thm 3.2": exact toric height -psi_Q(trop) <= width form; on
   non-squarefree profiles can be smaller; cap still holds since adversary profiles can be squarefree. wording.
## Verified OK (by hand)
- Lemma1.3 dictionary, 15/4 = 4/(1+2/30); Prop2.3 80kappa; alpha_bal(k)=(2^k-2)/2^{k-2}; Prop 4.3 = -(K_Y+E).C>=0 ; Prop5.1.
- Thm 3.3 minimax and inf_F = 2/tau - 1/2 (ff^* in wR_5, w = max width; Fekete sup).
- §7 count 30-2 vs 28; 2-2d; k=4 1+2d; argument principle 4/4.
## Lit (search snippets, record queries)
- Q "McKinnon Roth Seshadri constants ... Inventiones": Invent. Math. 200 (2015) 513-583, arXiv 1306.2976. OK.
- Q "An analogue of Liouville's theorem and an application to cubic surfaces": Eur. J. Math. 2 (2016) 929-959, arXiv 1306.2977;
  snippet: "when the locus being approximated is defined over the ground field, Liouville's Theorem is strictly better"; extension with
  asymptotic base locus (Ballay quotes: gamma with pi*L - gamma E effective, B = image of stable base locus). OK.
- Q Huang ANT: arXiv 1809.02001 v5 30 Apr 2021, universal torsor, "general fixed point", split toric, McKinnon conj under hypotheses,
  minimal centred primitive collections. Journal ref ANT 15 (2021) 461-512 NOT confirmed by any snippet (unverified, not refuted).
- Q McKinnon-Satriano: Trans AMS 374(5) (2021) 3557-3577, arXiv 2004.05212. Define CANONICALLY BOUNDED at P:
  alpha_{P,{x_i}}(-K_X) >= dim X for all Zariski dense sequences. Erratum 2023: if -K big, Vojta Main Conj => canonically bounded
  at P (original Prop 1.3 overclaimed). => toric Prop 5.1 is (snippet-level) a known statement; (McK^ess)/"not in literature" is
  wrong-ish; the expert question is "is Lbar_5^tw canonically bounded at e (norm-one points)?" Also Lehmann-McKinnon 2403.02480
  "essentially bounded", Monahan-Satriano 2308.11847 (horospherical, Vojta-conditional).
- canonical boundedness at e => A_5(delta) for all delta>0 (A-points not dense => closure is Z; dense => Zariski-dense sequence with
  alpha <= 4/(1+2delta) < 4).
## Computations
- basics_own.out: Ehrhart k=3,4,5; B^{k-1}=6,20,70; vertices = roots; widths; alpha_bal(k)=4(1-2^{1-k}) k=3..8; V in 3R_5; Vandermonde
  products give floor(10m/3) m<=9. ALL CONFIRMED.
- reg_small.out (C, p=2^31-1, own): k=5 m<=5 full at floor(10m/3), deficient at -1; k=4 m<=9 3m; k=3 m<=3 3m. CONFIRMED.
- k4curve_own.out: 5 parameter points: ImP const, 14 rays once, BC 14, contact 4, 7/2, D_V.C=24. CONFIRMED.
- pte_own k=4 calibration: 100 starts, 100 converged, 24 nondegenerate.
- reg_m6/m7/m8: FULL (5321, 9381, 15421). m=8 new: t_max(8)=26.
- pte_k5b.out: 44/80 converged, 0 nondegenerate (all coincident roots). Long run killed.
- Report written: kill_toric.md. Verdict: survives as restricted_theorem with corrections F1-F6.
