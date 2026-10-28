# Judge notes round 22

Read: expert_note, round20 README/gauss/kill_gauss, round10 rubric (a)(b)(c), overall weighs b,c over a.

Key prior facts:
- A-point: dist(y,e) <= H_B^{-(1/4+delta/2)}, H_B = R''^2. gamma = 4/(1+2delta).
- balanced Liouville delta=1/30 <-> gamma 15/4; CC delta=1/10 <-> 10/3; counting delta=0 <-> 4.
- Lemma K1 (killer, unrefereed by referee): gcd invisible.
- construct.md Thm B (mean width): ord_e f <= 2^{1-k} sum_{A subset [k]} width_Q(1_A). (sum over all A incl. complement pairs -> 2 per cut).
- gauss 4.4: balanced exclusion iff 2^{k-2} ord_e f > sum_cuts width(1_A).

Round 20 scores: gauss 5, composite 4, wild 3, arakelov 2. Decision stop.

To do:
- read lit.md, kill_lit.md, toric.md, kill_toric.md
- judge checks: (1) -K not ample at k=4,5 (support function); (2) (-K)^4 = 70 Ehrhart; (3) rank t_max(m) small m reproduce e.g. m<=4; (4) maybe check the minimax/width claim F1 of kill_toric; (5) McKinnon-Satriano canonical boundedness (TAMS 2021) relation.

## After reading all four files
- lit: verdict no_new_result survives; F1 major (-K not ample; X' Gorenstein Fano of root polytope). Prop 1.3 = row 20 restated; tau<=4-2^{3-k} new; ranks 3,6,10,13,16.
- toric: survives restricted_theorem; F1 overclaim (m<=7 data cap only single f with max cut width <=7, killer m=8 -> width<=8); F4 McKinnon-Satriano canonical boundedness.
- Killers disagree slightly on MS erratum wording: kill_lit "Vojta does not imply CB in general; does for toric"; kill_toric "if -K big Vojta => CB". Consistent if toric => -K big. Snippet-level.
- My WebSearch (2026-10-28): MS TAMS 2021 use canonical boundedness as a HYPOTHESIS (split toric surfaces, conditional on Vojta); no unconditional CB theorem found. Pasten 2205.07841: unconditional Vojta progress thin.
- Judge observation: tau_e, eps_e, balanced Liouville 15/4, covering threshold 4 are GEOMETRIC -> same for the SPLIT LM_5 at e over Q. There H_{-K}(u) = lcm(u) * (arch factor ~1 near e). So the regime tau < dim X at a rational point already occurs split; integer analogue A_5^split. Calibration only (neither implies the circle problem).
- Open toric point: weighted mean-width ord_e f <= (1/6) sum_{|A|=2} width(1_A) -> tau = 10/3. Equality on products of root differences. Test on random supports (evidence).
