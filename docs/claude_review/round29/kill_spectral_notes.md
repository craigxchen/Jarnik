# kill_spectral notes
## Plan: check (1) Li-Liu account vs PDF; (2) proofs 2.0-2.7; (3) numerics: ll_kernel, kappa, resonance; (4) lit; (5) F1/F2.
## Initial issues spotted while reading
- Prop2.0 consequence overstated: Q^eps-lossy bound on k-th moment, k growing, not vacuous (M^k <= Q^eps => M <= Q^{eps/k}).
- Lemma 2.3: unit u is redundant; every norm-1 elt of Q(i) is P/Pbar, P primitive unique up to sign. "only u=1 for non-trivial I" is a normalisation, not a consequence.
- heuristic kappa: my estimate E ~ (2/pi^2)(W/R)C^2 (1/2 log R + log C) ~ 0.13-0.15 (W/R)C^2 logR vs their 0.25-0.33. CHECK.
- Cor 2.5(b) heuristic m~(R/W)^2 from random hits: consistent.
## PDF check (pdftotext to stdout)
- (2.28),(2.30), Thm 6.1, Thm 6.2 ||T||<=||T||_HS<<(logN)^2, p.12 "they do not measure cancellation", Lemma 3.1 0<=P<=D<<log N, Prop A.4, App B 905 and ||C||/delta=0.826864 delta=0.4323, Sec 8 "central advance" all CONFIRMED verbatim.
- (7.16) Q_N = e^2/(2pi) int |V|^2 du/u + O(1/logN); V from Vtilde=(1+s)Phitilde => V = Phi - u Phi'. Lemma 2.7 consistent.
- script defs of D, K, P_N, X1 match (2.11),(2.12).
- spectral.md table row N=2018 not in ll_kernel_out.txt (only notes). 
## kernel recomputation (kernel_indep.py, mpmath Hurwitz exact deriv, Moebius prime zeta)
- P, K agree with lane to 6 digits; D agrees after including prime powers (my first bug). E(2018,1/9)=0.2373 confirms unlisted row.
- CAVEAT: sigma=0.05 lies OUTSIDE Hardy strip (alpha_N=1/logN > 0.05 for all N<=2e5); in-strip points: E=0.234(20014),-0.059(60022),-0.015,-0.0115(200006). "cancel to 1e-2" only at the largest N; no decay law visible. Weak evidence.
## kappa (kappa_indep.py, kappa_out.txt), R=7e3..1e6
- A/((W/R)C^2) = 2.6..4.35; slope in logR ~0.19-0.20 (W=R) matches my corrected heuristic (2/pi^2)(log R + 2 log(C/2)) + O(1) intercept. kappa_eff 0.24-0.44 (W=R/10, R=7e3 gives 0.44 outside lane's [0.25,0.33]). Asymptotic constant heuristically 2/pi^2=0.203. Lane's "kappa in [0.25,0.33]" = finite-size effective value; order log R confirmed.
## resonance (resonance_indep.py): W=5R: E(pred arc)=484,980,3904 >= pred 178,402,1788 CONFIRMED. W=R: sup_nontriv=140,238,740 (sqrtR=100,225,1000), mean~46-57; sup/sqrtR decreasing 1.4,1.06,0.74.
## small checks: lambda_k formula exact (also with 2: (-1)^{k e2}, q=3mod4 even exps factor 1). Prop2.6 0 collisions.
## structural: zbar P' y: E^W = sum_{P'} #{y: |P'|^2|y|^2 in window, arg y in shrunken sector} -> self-similar (circles of radius R/|P'|), lattices similar to Z[i]; "SL2(Z)\SL2(R) equidistribution" framing unnecessary.
## lit: Li-Liu "conditional proof" not found by WebSearch (same authors: arXiv 2606.05224 (1+1.9)); GGM Crelle 537 (2001) 105-126 confirmed via secondary (Chirre AMS ref list).
