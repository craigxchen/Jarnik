# finishing session notes (p3data.md)
Missing in draft: §0 Verdict, §4.2 (M=12,20,24).
Existing: q_lean M12 (8 inst), M20 (2 inst; t0 prime failed); ils M8 t0, M12 t0 X24, M20 t0,t1 X40 (crash fixed by guard).
ILS ratios W_ils/W_cut: M8 ~1.1-1.2, M12 1.33, M20 1.75 (r~90=4.5M). rwl1 worse (2.3-3.0 at M20).
Concern: growing ratio could be heuristic weakness. Plan: equal-weight SA code search (C):
  for tau >= w0, find S (M x r) with d_xy >= ceil(r/2 + D_xy/tau); W = r tau; LP polish weights (lp_S).
  exact equal-weight LBs: Plotkin-with-ceilings, Cor 6.4.  Perron 2 lam_max always valid LB.
  Then actual primes: realise_fixed_S with robust prime generation (no list limit).

## results so far
* eq-SA at M=8 (test8.py): within 0-15% of proved optima; 6-10 s/instance.
* eq-SA at M=20: no code found up to r_P+30 for any tau (test20.py) -> SA weak.
* SA calibration (test_unif.py): uniform targets d >= r/2+1 at M=20, r=38 (two Paley copies feasible):
  SA best V = 36-45; r=40,44,48: V 23-45; excess 1.5 at r=57..66: V 46-100.  => local search cannot find
  near-tight (Hadamard-type) codes at M=20.
* Hadamard-pool MILP (hpool.py) M=20 t0 X40: K=10 pool 190 -> W 439 (73 cols, 300s); K=30 pool 570 -> 579.6 (400s). Worse than ILS 363.
* Exact equal-weight Plotkin-with-ceilings LB (eqlb.py -> eqlb_out.txt/json): LB_eq/2lam medians
  M8 1.17/1.12, M12 1.19/1.15, M20 1.17/1.20, M24 1.18/1.20, M32 1.22/1.21 (X=2M/4M), flat; tau* ~ w0.
* heuristic count correction for §6.2: S-side local-CLT count: W=O(M log M) forces avg slack O(1) columns
  (Plotkin identity), so >= half pairs in O(1) windows; random-like count 2^{Mr}(c/sqrt r)^{M^2/4} needs
  r >~ (M/8) log2 r; same count predicts no Hadamard matrices (fails by 2^700 at n=32) -> decides nothing.
* Fixed unsupported §3 claim (ILS M8 calibration): ils2 M8 on six proved optima: gaps 1.1,1.2,8.5,4.4,4.2,0.0 % (median 2.7%).
  Realisation failed on 2 of 6 M8 ILS solutions (t4,t5 X16) -- ties/delta; not needed for table (M8 rows use §4.1 q_run data).
* SA calibration on Paley-stack targets (test_unif_out.txt): M8 solved, M12 1/6, M20/M24 never.
* Structured-start SA (test_struct.py) M20 t0 X40: polished 323.7 (1.61 x 2lam) < ILS 362.9.
* Remark 4.2: equal-weight ties (complementary D=0 pairs) unrealisable by distinct primes; best_ub refines on D+0.05[D=0].
* Example 6.5 (one collision layer cheap: W <= 2 lam + 2 M w0 via Kronecker of Hadamards).
* char MILP (|c|<=12) on M12 profile: min 0.12 at a pair, 1.3 s.

## final (finishing session)
* best_ub (structured SA + polish + ILS on D+0.05[D=0]) beats ILS at M=20 by 4-13%, at M=24 by 18-25%.
* Table 4.2b medians best/2lam: M8 1.32/1.22, M12 1.45/1.37, M20 1.59/1.64, M24 1.82/1.72; unif/best 2.4,2.3-2.7,2.6-2.7,2.9/3.7.
* 40 realised actual-prime profiles (M>=12): exact pairs + exact supp-4 OK; min s_c (|c|<=12) in [0.02,0.33], always at a pair.
* best_ub M24 ran under load ~20 (other session's mono.py), ~1200 s/job.
* Pitfall: `while pgrep -f "best_ub.py 12"` matched its own shell; queue never started.
