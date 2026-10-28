# Judge notes round 23
- Read: expert_note §3, round22 README, split.md + referee, width.md (+ referee summary).
- Hand re-derivation (judge): Lemma 2.4 n=5 step: X<=2(a-p)+q, Y<=(a+b-q)/2+p => X/3+2Y/3 <= a+b/3. OK.
- Plan checks: J1 recursion/peel exact (own code) + referee 7m/2 correction; J2 Thm 2 at new params (sympy) ;
  J3 PR #7218 C=7/2 claim (gh).
- J1 done (j1_recursion.out): rho_n==F_n on [0,70]^2 n=3..5; rho_5(m,m)=floor(10m/3); peel with C_bal(4) slice bound = 7m/2 (referee right), ThmB(4) = 15m/4, true = floor(10m/3).
- J2 done (j2_pell.out): Thm 2 norm identity + profile (1 const place, 15 distinct cuts, B.C=15) at (1,2,4) and 4 NEW points (2,5,9),(1,4,6),(3,4,10),(2,7,8); m-czech offsets -> 0, 2sqrt56=14.9666, 15, 18, 30; C -> 30; gcd alternates 1/9.
- J3: 887.lean fetched raw (dl/887.lean): parts.i/ii 'research open'; rosenfeld_4 'research solved' via z^2=56x^2+65. C=7/2 claim not in file (PR body only; gh REST blocked).
- J4 (judge derivation, unrefereed): Gaussian analog of split Thm 1: A_5(delta<1/30) => exists K=K(delta,Z) with
  sup_C M_infty(C) <= K (absolute CG). Proof: random 5-subset of the M points on arc C sqrt R; reduced radius
  log R''^2 = sum_p (max_S a - min_S a) log p (inert & 2 go to gcd) ; E <= (15/16+O(1/M)) log R^2 ;
  dist invariant = C R^{-1/2}; need (15/16+eta)(1+2delta)<1 ; Northcott on distinct tuples with R''<R_0 ;
  Schwartz-Zippel on {zeta_1..zeta_M}^5 vs F (complex coords zeta_j). No unit class needed.
  CONTRAPOSITIVE: if M_infty(C) -> infinity as C -> infinity, A_5(delta) false for every delta<1/30.
  => need lit check: is sup_C M_infty(C) known finite/infinite? Record: M_infty(C) >= 8 (C>5.9).
- J4 exact (j4_gauss_avg.out): 10-point cluster n=1176852625, span/sqrtR=8.2386 (chord); identity N''=prod p^(max a-min a) on all 252 5-subsets: True; mean ratio 0.9765 (M=10 too small to be informative for the 15/16 bound; check is of the identity only).
- Lit: CG Conjecture 14 (R^{1-eps}) ; Oganesyan 2107.09991 WITHDRAWN (record lit18 #26). sup_C M_infty(C) status: not found settled (search-negative).
- Scores plan: split (a)7 (b)0 (c)2 overall 4 ; width (a)9 (b)0 (c)2 overall 4. Decision: stop.
