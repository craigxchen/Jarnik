judge notes
- contradiction S4: coordinator/prove (r-1)m<=|S|-2+eps h vs bridge r m<=|S|-1+eps h. Reconcile: bridge = Vojta for SNC pair pi*bd+E; implies the other (subtract m). Bridge stronger, correct.
- contradiction balance faked (prove 2.5) vs kill_prove: check fake_lattice independently.
- M(7,0) upper: kill_bridge says 8 (Plotkin); search balanced clique 7 => <=7 proved.
- plan own check: saturated-lattice filter DFS at D=8 (and 7,9): is every rank-5 6-pt code killed? max filtered clique?
- D5 filtered max 3 (=M), D6 5 (=M), D7 5 => M(7,0)<=5 proved (judge, unrefereed); D8 running
- satdfs D8 (no symmetry) found 6-pt rank-5 code surviving filter B=2:
  rows [1,1,1,1,1,-1,-1,-1],[1,1,1,-1,-1,1,1,-1],[1,1,-1,1,-1,1,-1,1],[1,-1,1,1,-1,-1,1,1],[1,-1,-1,-1,1,1,1,1],[-1,1,1,-1,1,1,-1,1]
  saturated basis (index 1) of weight-4 balanced vectors; LLL coeffs<=4 no kill; Gordan feasible. LM 400 starts: all converged sols degenerate; best nondeg residual 1.8e-4 (evidence non-realisable). exp dim = 16-2-15 = -1.
- vertical-line killer sol 2 rechecked numerically: 5 zero-sum sets incl empty, rank 4.
- running sym8.py D=8, D=9 (symmetry reduced filtered DFS).
- D8 sym DFS FINAL: filtered max 8 => M(8,0)=8 [J]; 21 rank-5 6pt classes survive filter (classes_D8_c0.pkl)
- D9 FINAL filtered max 5 => M(9,0)<=5 [J]. realcls D8 partial: nondeg residual >=1e-4 or none for classes 0-8 (evidence)
