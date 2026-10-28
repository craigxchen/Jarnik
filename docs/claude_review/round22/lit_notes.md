# lit notes (round 22)
## Normalisation (derived)
- L = -K_{Lbar_5^tw} = sum_rho E_rho ample (permutohedral variety Fano: Voskresenskii-Klyachko). h_L([zeta]) = W' = log N' = log R''^2 (gauss Lemma 1.1).
- dist(y,e) ~ max_j |zeta_j/zeta_1 - 1| ~ diam/R''.
- A-point at delta: dist <= H^{-(1+2delta)/4}. McKinnon-Roth alpha (dist^gamma H bounded): gamma = 4/(1+2 delta).
- Liouville balanced: diam >= 2 R''^{7/15} -> dist >= H^{-4/15} -> gamma_L = 15/4. General distinct (CC): gamma >= 10/3.
- Counting / Vojta: gamma = 4 = dim. Vojta (D=E) => alpha_ess(-K) >= n at any rational point (generic computation).
- Needed: essential (off Z), balanced alpha > 15/4 strictly.
- Best approx: diagonal conics C_{1_A} (two distinct values), -K.C = 2, theta=1/2, alpha=2. So eps_e(-K) <= 2, alpha_e(-K) <= 2.
- 1-PS conic C_v, v distinct coords: -K.C_v = 2(max v - min v) >= 8, alpha along it = 2(max-min) >= 8 (theta <= 1/8).
- Split analogues: P^1xP^1 at (1,1): dense theta = 1/4 < count 1/2 (Liouville sharp, repulsion); split dP6 at e theta=1/3 vs count 1/2.
- Twisted k=5: Liouville theta 4/15 > count 1/4 -> Liouville NOT sharp predicted. k=4: Liouville 2/7 < count 1/3 (and attained Thm 4.1).
## Candidate papers
- McKinnon 2007 JAG conjecture; McKinnon-Roth Invent 2015; McKinnon-Roth EJM 2016 (Liouville analogue, cubic surfaces)
- McKinnon-Satriano toric (Trans AMS?); Lehmann-McKinnon-Satriano 2403.02480 (split surfaces via Vojta); Manon-McKinnon-Satriano 2504.05213 (Coba conj, wonderful varieties)
- Huang: dP6 (IJNT?), Acta Arith toric surface, ANT 15 (2021) toric; Huang-Wilsch 2509.03998 integral; Huang-de Saxce 2510.09954 flag; Huang 2111.01509 equidistribution toric geometric sieve
## Key derivation: tau_e(-K) = mean width bound
- Liouville-with-base-locus (McKinnon-Roth EJM extension): y outside B(pi^*L - gamma E) => dist^gamma H_L >= c. So alpha_ess >= gamma for gamma < tau_e(L) (effective threshold).
- Polytope of -K on Lbar_k: width in direction 1_A = 2 for all proper A. construct.md Thm B: ord_e f <= 2^{1-k} sum_A width_Q(1_A) => tau_e(-K) <= 2^{1-k}(2^k-2)*2 = 4 - 2^{3-k}.
- k=5: 15/4 = balanced Liouville gamma. k=4: 7/2 (= dist H^{-2/7}, Thm 4.1 exponent 3/7) consistent. k=3: 3.
- So every Seshadri/tau-type theorem caps at 15/4. Need > 15/4. Vojta: >= 4 = n. Regime tau_x(-K) < n is unusual (P^n: tau = n+1 > n).
- Inf alpha_e(-K) <= 2 (diagonal conics). Any lower bound on alpha_e over all points <= 2.
- Split analogues: tau sharp (dP6 split: dense gamma=3=tau; natural 2).
## Access: arxiv, ar5iv, queensu, smf, uwspace blocked. Mirror abstracts [A], WebSearch snippets [S].
- Monahan-Satriano 2308.11847: split horospherical, assuming Vojta, constructs curve through x whose points approximate better than any Zariski dense sequence (weak McKinnon). [S]
- Snippet: MR Invent Thm 6.3-ish: alpha >= eps/2 for rational point (cited in 2308.11847) [S, exact form unverified].
## tau computation (lit_checks/)
- |mP cap M| k=5: 21,131,471,1251 -> Ehrhart leading 70/24 -> (-K)^4 = 70, vol^{1/4}=2.89. k=3 hexagon 7,19,37,61.
- max ord_e of section of m(-K) = reg(mP cap M) (Hilbert regularity; ord>=r iff c perp Poly_{<r}).
- mod-p reg (upper bound on true max ord): k=3: 3m (m<=6) = bound 3. k=4: 3m (m<=4) < 3.5. k=5: 3,6,10,13 (m=1..4): max ratio 10/3 at m=3.
- Vandermonde prod_{i<j}(w_i-w_j): widths |P||P^c| <= 6 = 2*3 -> section of 3(-K), ord 10 -> ratio 10/3 attained over Q. Exactly CC exponent delta=1/10.
- So Liouville-with-base-locus (MR EJM) at e w.r.t. -K: alpha_ess >= 10/3 (off diagonals) = CC; upper limit tau<=15/4 (Thm B). Need >15/4 on balanced points: impossible by any tau-type theorem.
- k=4: tau-level 3 < 7/2 (balanced Liouville), consistent w/ Thm 4.1 since forced-div/pair on balanced points uses restricted heights.
- background: m=5 r=16 spanning test (reg_k5_m5.out).
